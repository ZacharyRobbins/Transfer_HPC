# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:16:55 2022

@author: zacha
"""
import os
os.chdir('C:/Users/zacha/Documents/GitHub/Transfer_HPC')
import SP_func as sp
#def rate_of_spread(wind, ### Site level wind 
#                   sum_fuel, ### the sum fuel from f() fuelchar
#                   Fuelbulkdensity, ### the FBD from f() fuelchar ( 0.40- 6.4)t-ac 
#                   fuel_sav, ### The surface area volume from f() fuelchar 1672- 2054 (this is ft-1)
#                   fuel_eff_moist, ### Thefuel effective moisture from f() fuelchar
#                   fuel_mef, ### The moisture extinction factor from f() fuelchar ## Commonly 15%dry or 40% wet (maybe 0.15-0.25)
#                   effect_wspeed, ### The effective windspeed from f() effective winspeed
#                   mineral_total,## UD portion of fuel that cant be burned.  Thonicke 0.01
#                   SF_val_part_dens, ## UD I believe this is roughly the Packing ratio of fuel   0.11-0.00143
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

ROS_front,ROS_back,ir=sp.rate_of_spread(600, ### Site level wind #Meters per min 
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
SF_val_part_dens=1/0.005
beta=Fuelbulkdensity/SF_val_part_dens ### Ithink this is 0.00143 in scott
beta_op=0.200395*((fuel_sav)**(-0.8189))
beta_ratio = beta/beta_op 


 xi = (np.exp((0.792 + 3.7597 * (fuel_sav**0.5)) * (beta+0.1))) / (192+7.9095 * fuel_sav)
