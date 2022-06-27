# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:16:55 2022

@author: zacha
"""
import os
os.chdir('C:/Users/zacha/Documents/GitHub/Transfer_HPC/Py_scripts/')
import SP_func as sp
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import scipy


#######Climate inputs############



######## The parameters for the fuel classes#########
######### These will not change within run ##########
NFSC=pd.DataFrame({'TissueType':['tw_sf','sb_sf','lb_sf','tr_sf','deadleaves','livegrass'],
      'FBD':[15.4,16.8,19.6,999.0,4.0,4.0],
      'SAV':[13.00,3.58,0.98,0.20,66.00,66.00],
      'low_moisture_Coeff':[1.12, 1.09,0.98,0.80,1.15,1.15],
      'low_moisture_Slope':[0.62,0.72 ,0.85, 0.80, 0.62, 0.62],
      'mid_moisture':[0.72, 0.51, 0.38, 1.00, 0.80, 0.80],
      'mid_moisture_Coeff':[2.35, 1.47, 1.06, 0.80, 3.20, 3.20],
      'mid_moisture_Slope':[2.35, 1.47, 1.06, 0.80, 3.20, 3.20],
      'min_moisture':[0.18, 0.12, 0.00, 0.00, 0.24, 0.24]
      })

##################Blank stand dataframe################################
SV_TW=0.022 # Starting biomass twigs
SV_SB=0.0464 #Starting biomass small branch
SV_LB=0.12998#Starting biomass large branch
SV_TR=0.6108 # Starting biomass trunk
SV_DL=0.3318 #Starting biomass dead leaves
SV_LG=0.02035 #Starting biomass live grass

Vegframe=pd.DataFrame({'TissueType':['tw_sf','sb_sf','lb_sf','tr_sf','deadleaves','livegrass'],
          'Biomass':[SV_TW,SV_SB,SV_LB,SV_TR,SV_DL,SV_LG],
          'MEF':[0,0,0,0,0,0],
          'alphaFMC':[0.0,0.0,0.0,0.0,0.0,0.0],
          'fuelmoisture':[0,0,0,0,0,0],
          'Fuelfrac':[0,0,0,0,0,0],
          'litter_moisture':[0,0,0,0,0,0],
          'BF':[0,0,0,0,0,0],
          'tau_b':[0,0,0,0,0,0]          
          })
#### Other user defined parameters

NI_param_a=12 ## NI_parameter_a
NI_param_b=76 ## NI_parameter_b
 ### Percent of the stand that is tree fraction 
cg_strikes=1
ED_val_nignitions=1
SFval_max_durat=24
SF_val_fuel_energy=18000

#SF_val_miner_damp=.30
#SF_val_drying_ratio=100000
#mineral_total=.02

#SF_val_fdi_alpha=2.000 ### NI to FDI 
#SF_val_part_dens=1/0.0011
#SF_val_durat_slope=1.1
#tree_fraction=1.0


######### Begin##############
#############################
#Site_Ni=10000
#wind=20
##Weather 
#fire_danger_index(Site_Ni,Daily_Temp_C,Daily_Rainfall,Daily_rh,
#                       NI_param_a,NI_param_b)


Outy=pd.DataFrame(columns=["Run","tree_fraction","SF_val_drying_ratio","mineral_total","SF_val_miner_damp","SF_val_part_dens","SF_val_durat_slope",'wind',"Site_Ni","Fractionburned","Fire_intensity","AB",'ROS_front',"fuel_moisture"])


for r in list(range(1,100000)):
    print(r)
    ###Generate Parameters
    wind=np.random.uniform(2,10)
    Site_Ni=np.random.uniform(3000,12000)
    SF_val_drying_ratio=np.random.uniform(100,5000)
    mineral_total=np.random.uniform(.01,.1)                                   
    SF_val_miner_damp =np.random.uniform(0.3,0.5) 
    SF_val_fdi_alpha=np.random.uniform(1.0,10.0) 
    SF_val_part_dens=np.random.uniform(1/0.011,1/0.0011) 
    SF_val_durat_slope=np.random.uniform(0.0,.2) 
    tree_fraction=np.random.uniform(0.0,1.0)                          
    
    
    
    
    effect_wspeed=sp.wind_effect(tree_fraction,0, wind)
    FBD,SAV,MEF,fuel_moisture,sum_fuel,Vegframe=sp.fuel_char(Vegframe,NFSC,SF_val_drying_ratio,Site_Ni)
    
    Vegframe,TFC_ROS,tau_si1=sp.ground_fuel_consumption(Vegframe,NFSC, sum_fuel)
    
    ROS_front,ROS_back,ir=sp.rate_of_spread(wind, ### Site level wind 
                       sum_fuel, ### the sum fuel from f() fuelchar
                       FBD, ### the FBD from f() fuelchar ( 0.40- 6.4)t-ac 
                       SAV, ### The surface area volume from f() fuelchar 1672- 2054 (this is ft-1)
                       fuel_moisture, ### Thefuel effective moisture from f() fuelchar
                       MEF, ### The moisture extinction factor from f() fuelchar ## Commonly 15%dry or 40% wet (maybe 0.15-0.25)
                       effect_wspeed, ### The effective windspeed from f() effective winspeed
                       mineral_total,## UD portion of fuel that cant be burned.  Thonicke 0.01
                       SF_val_part_dens, ## UD I believe this is roughly the Packing ratio of fuel   0.11-0.00143
                       SF_val_fuel_energy, ## UD Roughly energy per fuel  fuel (18 000) kJ kg−1 other places (20 000 kJ kg−1),
                       SF_val_miner_damp)  ## UD Another fuel char: dampening. Thonicke 0.41739
    
    Fractionburned,Fire_intensity,AB=sp.area_burnt_intensity(SF_val_fdi_alpha, ##relative importance of the nestrov index 
                             Site_Ni, ## Current Nesterov index
                             cg_strikes, ##Ration of cloud to ground strikes
                             ED_val_nignitions, ###Valid ignitions
                             SFval_max_durat, ### Max duration of fires
                             SF_val_durat_slope, ### relationship between fire weather and duration of fire
                             effect_wspeed, ### effective windspeed. 
                             tree_fraction, ## Fraction of a patch that is tree
                             ROS_front, ## Rate of spread forward m/min
                             ROS_back, ## Rate of spread forward m/min
                             TFC_ROS, ### total fuel consumed by rate of spread. KgC/m2
                             SF_val_fuel_energy ### the energy value of fuel kJ/kg
                             ) 
    
    Inny=pd.DataFrame({"Run":r,
                       "tree_fraction":tree_fraction,
                       "SF_val_drying_ratio":SF_val_drying_ratio,
                       "mineral_total":mineral_total,
                       "SF_val_miner_damp":SF_val_miner_damp,
                       "SF_val_part_dens":SF_val_part_dens,
                       "SF_val_durat_slope":SF_val_durat_slope,
                       'wind':wind,
                       "Site_Ni":Site_Ni,
                       "Fractionburned":Fractionburned,
                       "Fire_intensity":Fire_intensity,
                       "AB":AB,
                       'ROS_front':ROS_front,
                       "fuel_moisture":fuel_moisture
                       })
    Outy=Outy.append(Inny)
Outy.to_csv("FireSpread_100k_2.csv")
