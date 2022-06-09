# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:19:23 2022

@author: zacha
"""

import numpy as np

def rate_of_spread(wind, ### Site level wind 
                   sum_fuel, ### the sum fuel from f() fuelchar
                   Fuelbulkdensity, ### the FBD from f() fuelchar ( 0.40- 6.4)t-ac 
                   fuel_sav, ### The surface area volume from f() fuelchar 1672- 2054 (this is ft-1)
                   fuel_eff_moist, ### Thefuel effective moisture from f() fuelchar
                   fuel_mef, ### The moisture extinction factor from f() fuelchar ## Commonly 15%dry or 40% wet (maybe 0.15-0.25)
                   effect_wspeed, ### The effective windspeed from f() effective winspeed
                   mineral_total,## UD portion of fuel that cant be burned.  Thonicke 0.01
                   SF_val_part_dens, ## UD I believe this is roughly the Packing ratio of fuel   0.11-0.00143
                   SF_val_fuel_energy, ## UD Roughly energy per fuel  fuel (18 000) kJ kg−1 other places (20 000 kJ kg−1),
                   SF_val_miner_damp):  ## UD Another fuel char: dampening. Thonicke 0.41739
        
    q_dry = 581.0 ### heat of pre ignition of dry fuels. 
    ### remove mineral content from net fulel load (Thonicke 2010)
    sum_fuel=sum_fuel*(1-mineral_total)
    beta=Fuelbulkdensity/SF_val_part_dens ### Ithink this is 0.00143 in scott
    beta_op=0.200395*(fuel_sav)**(-0.8189)
    beta_ratio = beta/beta_op  
    ### heat of pre ignition 
    # ---heat of pre-ignition---
    #   !  Equation A4 in Thonicke et al. 2010
    #   !  Rothermal EQ12= 250 Btu/lb + 1116 Btu/lb * fuel_eff_moist
    #   !  conversion of Rothermal (1972) EQ12 in BTU/lb to current kJ/kg 
    #   !  q_ig in kJ/kg 
    q_ig = q_dry +2594.0 * fuel_eff_moist
    #! ---effective heating number---
    #! Equation A3 in Thonicke et al. 2010.  
    eps = np.exp(-4.528/ fuel_sav)     
    #! Equation A7 in Thonicke et al. 2010
    b = 0.15988 * (fuel_sav**0.54)
    #! Equation A8 in Thonicke et al. 2010
    c = 7.47 * (np.exp(-0.8711 * (fuel_sav**0.55)))
    #! Equation A9 in Thonicke et al. 2010. 
    e = 0.715 * (np.exp(-0.01094 * fuel_sav))
    phi_wind = c * ((3.281*effect_wspeed)**b)*(beta_ratio**(-e)) 
    #! ---propagating flux----
    # Equation A2 in Thonicke et al.2010 and Eq. 42 Rothermal 1972
    #! xi (unitless)       
    xi = (np.exp((0.792 + 3.7597 * (fuel_sav**0.5)) * (beta+0.1))) / (192+7.9095 * fuel_sav)      
    ## ! ---reaction intensity----
    ##   ! Equation in table A1 Thonicke et al. 2010. 
    a = 8.9033 * (fuel_sav**(-0.7913))
    a_beta = np.exp(a*(1.0-beta_ratio))  #!dummy variable for reaction_v_opt equation
    ## ! Equation in table A1 Thonicke et al. 2010.
    ## ! reaction_v_max and reaction_v_opt = reaction velocity in units of per min
    ## ! reaction_v_max = Equation 36 in Rothermal 1972 and Fig 12 
    reaction_v_max  = 1.0 / (0.0591 + 2.926* (fuel_sav**(-1.5)))
    ## ! reaction_v_opt =  Equation 38 in Rothermal 1972 and Fig 11
    reaction_v_opt = reaction_v_max*(beta_ratio**a)*a_beta
    # ! mw_weight = relative fuel moisture/fuel moisture of extinction
    # ! average values for litter pools (dead leaves, twigs, small and large branches) plus grass
    mw_weight = fuel_eff_moist/fuel_mef
    # ! Equation in table A1 Thonicke et al. 2010. 
    # ! moist_damp is unitless
    moist_damp = max(0.0,
                     (1.0 - (2.59 * mw_weight) + (5.11 * (mw_weight**2.0)) - (3.52*(mw_weight**3.0))))   
    # ! ir = reaction intenisty in kJ/m2/min
    #  currentPatch%sum_fuel converted from kgC/m2 to kgBiomass/m2 for ir calculation
    ir = reaction_v_opt*(sum_fuel/0.4)*SF_val_fuel_energy*moist_damp*SF_val_miner_damp   
    ##if fuel or effective heat or heat of pre-igintion is zero 
    if (Fuelbulkdensity<=0.0 or eps<= 0 or q_ig <=0.0) :
        ROS_front =0.0
    else: 
        ROS_front=(ir*xi*(1.0+phi_wind)) / (Fuelbulkdensity*eps*q_ig)
    ROS_back=ROS_front*np.exp(-0.012*wind)    ## check wind
    return(ROS_front,ROS_back,ir)


def fire_danger_index(Site_Ni,Daily_Temp_C,Daily_Rainfall,Daily_rh,
                      NI_param_a,NI_param_b):    
    if (Daily_Rainfall >3.0): ### if rain > 3.0 reset NI
        Site_Ni=0
    else:
        yipsolon=(NI_param_a*Daily_Temp_C)/(NI_param_b+Daily_Temp_C)+np.log(Daily_rh/100)
        dewpoint=(NI_param_b*yipsolon)/(NI_param_a-yipsolon)
        d_NI=(Daily_Temp_C-dewpoint)*Daily_Temp_C
        if(d_NI<0.0): 
            d_NI=0.0
        Site_Ni+d_NI
    return(Site_Ni)

def fuelcalc(SAV,Site_Ni,SF_val_drying_ratio):
    fuel_mef= 0.524 - 0.066 * np.log(fuel_sav) 
    fuel_eff_moist= 2.71828**(-1.0 * (fuel_sav/SF_val_drying_ratio)*Site_Ni) 
    return(fuel_mef,fuel_eff_moist)


MEF,FuelMoisture=fuelcalc(126,1000,66000)
print(MEF,FuelMoisture)
