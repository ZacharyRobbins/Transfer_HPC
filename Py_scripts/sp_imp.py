# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:16:55 2022

@author: zacha
"""
import os
os.chdir('C:/Users/345578/Documents/GitHub/Transfer_HPC/Py_scripts/')
import SP_func as sp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy
#scipy.optimize.curve_fit(f, xdata, ydata, kwargs)
#def rate_of_spread(wind, ### Site level wind 'm/s'
#                   sum_fuel, ### the sum fuel from f() fuelchar
#                   Fuelbulkdensity, ### the FBD from f() fuelchar ( 0.40- 6.4)t-ac 
#                   fuel_sav, ### The surface area volume from f() fuelchar 1672- 2054 (this is ft-1)
#                   fuel_eff_moist, ### Thefuel effective moisture from f() fuelchar
#                   fuel_mef, ### The moisture extinction factor from f() fuelchar ## Commonly 15%dry or 40% wet (maybe 0.15-0.25)
#                   effect_wspeed, ### The effective windspeed from f() effective winspeed
#                   mineral_total,## UD portion of fuel that cant be burned.  Thonicke 0.01
#                   SF_val_part_dens, ## UD I believe this is roughly the 1/ Packing ratio of fuel   0.11-0.00143
#                   SF_val_fuel_energy, ## UD Roughly energy per fuel  fuel (18 000) kJ kg−1 other places (20 000 kJ kg−1),
#                   SF_val_miner_damp):  ## UD Another fuel char: dampening. Thonicke 0.41739
        
ROS_front,ROS_back,ir=sp.rate_of_spread(10,
                                        2,
                                        0.1,
                                        50.16,
                                        .10,
                                        0.15,
                                        4,0.,
                                        15,1/0.01,18000,0.41739)

#ROS in #m/min

ROS_front,ROS_back,ir=sp.rate_of_spread(60, ### Site level wind #Meters per min 
                   .2, ### the sum fuel from f() fuelchar  
                   0.1, ### the FBD from f() fuelchar ( 0.40- 6.4)t-ac ##FATES os 0.4 * 0.1
                   50, ### The surface area volume from f() fuelchar 1672- 2054 (this is ft-1)## FATES is cm-1
                   0.155, ### The fuel effective moisture from f() fuelchar
                   0.16, ### The moisture extinction factor from f() fuelchar ## Commonly 15%dry or 40% wet (maybe 0.15-0.25)
                   300, ### The effective windspeed from f() effective winspeed
                   0.01 ,## UD portion of fuel that cant be burned.  Thonicke 0.01
                   1/0.01, ## UD I believe this is roughly the Packing ratio of fuel   0.11-0.00143
                   18000, ## UD Roughly energy per fuel  fuel (18 000) kJ kg−1 other places (20 000 kJ kg−1),
                   0.41739)  ## UD Another fuel char: dampening. Thonicke 0.41739

0.40/1/0.005
fuel_sav=1672*0.03
Fuelbulkdensity=0.4
SF_val_part_dens=1/0.0011
beta=Fuelbulkdensity/SF_val_part_dens ### Ithink this is 0.00143 in scott
beta_op=0.200395*((fuel_sav)**(-0.8189))
beta_ratio = beta/beta_op 


xi = (np.exp((0.792 + 3.7597 * (fuel_sav**0.5)) * (beta+0.1))) / (192+7.9095 * fuel_sav)

  grassfraction=np.min(grass_fraction,1-tree_fraction)
    bare_fraction=1-(tree_fraction+grassfraction)
    total_tree area
    effective_wspeed= wind *(Treefraction*0.04+(grass_fraction+bare_fraction)*0.6)


#def fire_danger_index(Site_Ni,Daily_Temp_C,Daily_Rainfall,Daily_rh,
#                      NI_param_a,NI_param_b):  
    
Sample=pd.DataFrame({
    'WS':np.random.randint(0,20,10000), #m.s
    'sum_fuel': np.random.randint(0,20,10000), # 'kg C / m2'
    'FBD': np.random.randint(0,20,10000)/100, #'kg biomass/m3'
    'fuel_sav':np.random.randint(10,500,10000),  #'per m'* 0.092903
    'fuel_eff_moist':np.random.randint(1,40,10000)/100,
    'treefraction':np.random.randint(1,100,10000)/100,
    'SF_val_part_dens':np.random.randint(10,1000,10000),
    'SF_val_fuel_energy': np.random.randint(10000,20000,10000),
    'SF_val_miner_damp' : np.random.randint(0,100,10000)/100
    } 
)
Sample['fuel_MEF']=Sample['fuel_eff_moist']+np.random.randint(1,20,10000)/100
Sample = Sample[Sample['fuel_MEF'] >= 0]

Sample['grassfraction']=1-Sample['treefraction']
Sample['effws']=Sample['WS']*((Sample['grassfraction']*0.6)+(Sample['treefraction']*0.4))

for i in list(range(0,len(Sample))):
    print(i)
    ROS_front,ROS_back,ir=sp.rate_of_spread(Sample.iloc[i]['WS'],
                                            Sample.iloc[i]['sum_fuel'],
                                            Sample.iloc[i]['FBD'],
                                            Sample.iloc[i]['fuel_sav'],
                                            Sample.iloc[i]['fuel_eff_moist'],
                                            Sample.iloc[i]['fuel_MEF'],
                                            Sample.iloc[i]['effws'],
                                            0.01,
                                            Sample.iloc[i]['SF_val_part_dens'],
                                            Sample.iloc[i]['SF_val_fuel_energy'],
                                            Sample.iloc[i]['SF_val_miner_damp'],
                                            )
    if i ==0:
        ROS_f_vec=ROS_front 
        ROS_b_vec=ROS_back         
        ir_vec=ir
    else:
        ROS_f_vec=np.append(ROS_f_vec,ROS_front)         
        ROS_b_vec=np.append(ROS_b_vec,ROS_back)      
        ir_vec=np.append(ir_vec,ir)
#Cal= fuelmef and windspeed based on tree fraction 


Sample['ROS_f']=ROS_f_vec

Sample['ROS_b']=ROS_f_vec
Sample['IR']=ir_vec
plt.plot(Sample['WS'],Sample['ROS_f'], 'o', color='black')
plt.plot(Sample['sum_fuel'],Sample['ROS_f'], 'o', color='black')
plt.plot(Sample['FBD'],Sample['ROS_f'], 'o', color='black')
plt.plot(Sample['fuel_sav'],Sample['ROS_f'], 'o', color='black')
plt.plot(Sample['fuel_MEF']-Sample['fuel_eff_moist'],Sample['ROS_f'], 'o', color='black')
'SF_val_miner_damp'

