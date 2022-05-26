# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd

### All User defined parameter list 
   
    #NI_param_a
    #NI_param_b
    
### All drivers 
    # Daily_Temp_C
    # Daily_Rainfall
    # Daily_rh,
    
    
### NSFC is for each tissue category     
     #There are SIX fuel classes
     # 1:4) four CWD_AG pools (twig, s branch, l branch, trunk), 5) dead leaves and 6) live grass
     # NCWD =4  NFSC = 6
     # tw_sf = 1, lb_sf = 3, tr_sf = 4, 


NFSC={TissueType:['tw_sf','sb_sf','lb_sf','tr_sf','deadleaves','livegrass'],
     
     
     SF_val_min_moisture:
     SF_val_low_moisture_Coeff:
     SF_val_low_moisture_Slope:
     SF_val_mid_mid moisture:
     SF_val_mid_moisture_Coeff:
     SF_val_mid_moisture_Slope:
     FC_ground:
     tau_b}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
### For the subroutuine fire model 

### Patch from oldest to youngest. 
    #intialize burned area and fire to zero
### For each site
    # calculate4 fire data index
    # calculate wind effect
    # chars of fuel
    # rate of spread 
    # ground fuel consumption
    # area burnt intensity
    # crown scorching
    # crown damage
    #cambila damage kill
    # post fire mortlaity 
###~~~~~~~~~~~~~~~~~~~~~~####
###~~~~~Functions~~~~~~~~####
###~~~~~~~~~~~~~~~~~~~~~~####


##################    
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
#####################





### Perhaps we remove this and just have it be 
#fuel bulk desisty
#fuel surface area volume
#fuel moisture extinstion factor and effective mosiiture. 



def fuel_char(drying_ratio, SAV, FBD):
    litt_c= currentPatch%litter(element_pos(carbon12_element))
    ##using each cohort
    ## caculate the cohorts that are grass
    livegrass=livegrass+leaforgan_carbon
    sumfuel=leaf_fines+agcwd+livegrass

    #There are SIX fuel classes
     # 1:4) four CWD_AG pools (twig, s branch, l branch, trunk), 5) dead leaves and 6) live grass
     # NCWD =4  NFSC = 6
     # tw_sf = 1, sb_sf lb_sf = 3, tr_sf = 4, dl_sf = 5, lg_sf = 6,
    ### calculate fraction of fuels in each class
    
    fuel_tw_sd=agcwd/sumfeul
    #fuel_dl_sf=agcwd/sumfeul
    fuel_dl_sf=leaf_fines/sumfuel
    fuel_lg=livegrass/sumfuel
    
    ### This could be a subfunction
    ### Calculate the drying ratio of no live grass fuels
    alpha_FMC(tw_sf:dl_sf)      = SF_val_SAV(tw_sf:dl_sf)/SF_val_drying_ratio
    ### calcualte fuel moisture for each class
    tw_sf:dl_sf = exp(-1*((tw_sf:dl)sf)*Site_Ni)
    
    ## Fuelbulkdensity=currentPatch%fuel_bulkd +sum(fuelfraction(f1-f5)*Weight of each pool)?
    ## Fuel_Savy=currentPatch%fuel_sav+sum(fuelfraction(f1-f5)*Weight of each pool)### surface area volume
    ## Fuel_mef=currentPatch%fuel_mef  
    ## Fuel_eff_moist 
    ## currentPatch%fuel_eff_moist
    
    ##Basically we need to sum fuel bulk density,surface area volume andmoisture of extinction factor.7
    
    
    
    
    currentPatch%fuel_sav = sum(SF_val_SAV(1:nfsc))/(nfsc) ! make average sav to avoid crashing code. 


def wind_effect():
    total grass area
    tree fraction
    grass_fraction
    bare_fraction
    ## Convert wind into m/min
    ## Calculate total tree area
    ## Calculate total grass area
    ## calculate tree fraction and grass fraction 
    ##check continuity
    grassfraction=np.min(grass_fraction,1-tree_fraction)
    bare_fraction=1-(tree_fraction+grassfraction)
    total_tree area
    effective_wspeed= wind *(Treefraction*0.04+(grass_fraction+bare_fraction)*0.6)



def rate_of_spread(sum_fuel, ### the sum fuel from f() fuelchar
                   mineral_total,## portion of fuel that cant be burned. 
                   
        ):
    q_dry = 581.0 ### heat of pre ignition of dry fuels. 
    ### remove mineral content from net fulel load (Thonicke 2010)
    sum_fuel=sum_fuel(1-mineral_total)
    beta=Fuelbulkdensity/SF_val_part_dens
    bet_op=0.200395*(fuel_sav)**(-0.8189)
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
    eps = exp(-4.528/ fuel_sav)     
    #! Equation A7 in Thonicke et al. 2010
    b = 0.15988_r8 * (fuel_sav**0.54_r8)
    #! Equation A8 in Thonicke et al. 2010
    c = 7.47 * (exp(-0.8711 * (fuel_sav**0.55_r8)))
    #! Equation A9 in Thonicke et al. 2010. 
    e = 0.715 * (exp(-0.01094 * fuel_sav))
    phi_wind = c * ((3.281*effect_wspeed)**b)*(beta_ratio**(-e)) 
    #! ---propagating flux----
    # Equation A2 in Thonicke et al.2010 and Eq. 42 Rothermal 1972
    #! xi (unitless)       
    xi = (exp((0.792_r8 + 3.7597_r8 * (fuel_sav**0.5_r8)) * (beta+0.1_r8))) / &
            (192_r8+7.9095_r8 * fuel_sav)      
    ## ! ---reaction intensity----
    ##   ! Equation in table A1 Thonicke et al. 2010. 
    a = 8.9033_r8 * (currentPatch%fuel_sav**(-0.7913_r8))
    a_beta = exp(a*(1.0_r8-beta_ratio))  !dummy variable for reaction_v_opt equation
    ## ! Equation in table A1 Thonicke et al. 2010.
    ## ! reaction_v_max and reaction_v_opt = reaction velocity in units of per min
    ## ! reaction_v_max = Equation 36 in Rothermal 1972 and Fig 12 
    reaction_v_max  = 1.0_r8 / (0.0591_r8 + 2.926_r8* (currentPatch%fuel_sav**(-1.5_r8)))
    ## ! reaction_v_opt =  Equation 38 in Rothermal 1972 and Fig 11
    reaction_v_opt = reaction_v_max*(beta_ratio**a)*a_beta
    # ! mw_weight = relative fuel moisture/fuel moisture of extinction
    # ! average values for litter pools (dead leaves, twigs, small and large branches) plus grass
    mw_weight = currentPatch%fuel_eff_moist/currentPatch%fuel_mef
    # ! Equation in table A1 Thonicke et al. 2010. 
    # ! moist_damp is unitless
    moist_damp = max(0.0_r8,(1.0_r8 - (2.59_r8 * mw_weight) + (5.11_r8 * (mw_weight**2.0_r8)) - &
            (3.52_r8*(mw_weight**3.0_r8))))   
    # ! ir = reaction intenisty in kJ/m2/min
    #  currentPatch%sum_fuel converted from kgC/m2 to kgBiomass/m2 for ir calculation
    ir = reaction_v_opt*(sum_fuel/0.45_r8)*SF_val_fuel_energy*moist_damp*SF_val_miner_damp   
    ##if fuel or effective heat or heat of pre-igintion is zero 
    if fuel density <=0.0 .or. eps.<=0 .or. q_ig<=0.0:
        ROS_front =0.0
    else: 
        ROS_front=(ir*xi*(1.0+phi_wind)) / (fuel_bulkd*eps*q_ig)
    ROS_back=ROS_front*exp(-0.012*wind)    


    
def ground_fuel_consumption(NFSC,###DataFrame of tissues
                            moist ### Calculated litter moisture
                            agcwd ### pool of above ground course woody debris
                            leaf_fines ###pool of leaves
                            livegrass #### pool of live grass.
                            ):
    ## moist ! effective fuel moisture
    ## ta_b(nfsc) lethal heating rate for each fuel class
    ## fc_ground(nfsc) ! Propriont of fuel consumbed
    for 1:nfsc
        if moist<=SF_val_min_moisture():
            pool_burntfraction_litter=1.0 ### All fuel in pool is consumed
        if moist >SF_val_mid_moisture .and. moist <=SF_val_mid_mid moisture(c):
            pool_burntfraction_litter=min(1.0_r8,SF_val_low_moisture_Coeff(c)- 
                  SF_val_low_moisture_Slope(c)*moist)) 
        if moist > SF_val_mid_moisture and moist <=1.0:
            pool_burntfraction_litter= currentPatch%burnt_frac_litter(c) = 
            max(0.0_r8,min(1.0_r8,SF_val_mid_moisture_Coeff- SF_val_mid_moisture_Slope(c)*moist))
        if moist > 1.0:
            pool_burntfraction_litter=0.0
    ###Sum pools based on the ag_csd and leaf_fines ## and live grass
    FC_ground(nfsc)= pool_burntfraction_litter*agcwd or leaf_fines .or.livegrass
    ###calculate tau_b for each pool
    for pool:
        tau_b[i]=39.4*(fuel_frac*sumfuel/0.45/19)*(1-((1-burnt_frac_litter)**0.5))
    tau_b(tr_sf)=0.0
    ### Calculate residence time (cap at 8 mins)
    tau1=min(8.0,sum(tau_b))
    ### Remove 
    TFC_ROS=sum(FC_ground)-FC_ground(tr_sf) ### Total fuel consumed per ROS
    
    
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~!           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#    
def area_burnt_intensity(SF_val_fdi_alpha, ##relative importance of the nestrov index 
                         Site_Ni, ## Current Nesterov index
                         cg_strikes, ##Ration of cloud to ground strikes
                         ED_val_nignitions, ###Valid ignitions
                         SFval_max_durat ### Max duration of fires
                         SF_val_durat_slopeFDI ### relationship between fire weather and duration of fire
                         effect_wspeed ### effective windspeed. 
                         tree_fraction_patch ## Fraction of a patch that is tree
                         ROS_front ## Rate of spread forward m/min
                         ROS_back ## Rate of spread forward m/min
                         TFS_ROS ### total fuel consumed by rate of spread. KgC/m2
                         SF_val_fuel_energy ### the energy value of fuel kJ/kg
                         ): 
    forest_grassland_lengthtobreadth_threshold = 0.55
    FDI=1.0-exp(-SF_val_fdi_alpha*Site_Ni)    ### Fuel~weather relaitionship
    ### strikes
    NF = ED_val_nignitions * years_per_day * cg_strikes ##Gross number of ignitions* 1/365 *how many cloud to ground strikes
    ### Fire duration calculation
    Fire_duration=(SFval_max_durat+1)/(1+SFval_max_durat*exp(SF_val_durat_slope*FDI))
    #### The effective of fuel type of fire movement: 
     if (effect_wspeed*16.67<.1): #16.67m/min = 1km/hr
        lb=1
    ### Forest fuels    
    else:
        if(tree_fraction_patch > forest_grassland_lengthtobreadth_threshold):
            lb=((1.0+8.729*((1.0-(exp(-0.03*16.67*effect_wspeed)))**2.155)))
    ### Grass fuels 
        else:
            lb=1.1*((16.67*effect_wspeed)**0.464)
    ###Lenght of the major axis of fire spread
    db=ROS_back*Fire_duration #m
    df=ROS_front*Fire_duration #m 
    if lb>0.0:
        size_of_fire = ((3.1416/(4.0*lb))*((df+db)**2.0))
        AB=size of fire*FDI*NF ### Area burned is equal to size of fire by fuel dryness index and number of fires in a day. 
        Fractionburned=min((0.99,AB/1000000.0))## in KM ### return to line 757
    else: 
        Fractionburned=0.0
    ### Calculating intensity 
    ROS=ROS_Front/60 #m/min to m/sec
    W=TFS_ROS/0.45 #kgC/m2 to kgbiomass/m2 of burned area. 
    Fire_intensity=SF_val_fuel_energy*W*ROS*Fractionburned
    return(Fractionburned,Fire_intensity)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
area_burnt_intensity()


