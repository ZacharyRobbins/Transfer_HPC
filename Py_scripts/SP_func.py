# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:19:23 2022

@author: zacha
"""

import numpy as np
grass_fraction=0
def wind_effect(tree_fraction,grass_fraction, wind ):
    grassfraction=min(grass_fraction,1-tree_fraction)
    bare_fraction=1-(tree_fraction+grass_fraction)
    effective_wspeed= wind *(tree_fraction*0.4+(grass_fraction+bare_fraction)*0.6)
    return(effective_wspeed)

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
        Site_Ni=Site_Ni+d_NI
    return(Site_Ni)

#def fuelcalc(SAV,Site_Ni,SF_val_drying_ratio):
#    fuel_mef= 0.524 - 0.066 * np.log(fuel_sav) 
#    fuel_eff_moist= 2.71828**(-1.0 * (fuel_sav/SF_val_drying_ratio)*Site_Ni) 
#    return(fuel_mef,fuel_eff_moist)


def fuel_char(Vegframe,NFSC,SF_val_drying_ratio,Site_Ni):
    #litt_c= currentPatch%litter(element_pos(carbon12_element))
    ##using each cohort
    ## caculate the cohorts that are grass
    #livegrass=livegrass+leaforgan_carbon
    sumfuel=sum(Vegframe['Biomass'])

    #There are SIX fuel classes
     # 1:4) four CWD_AG pools (twig, s branch, l branch, trunk), 5) dead leaves and 6) live grass
     # NCWD =4  NFSC = 6
     # tw_sf = 1, sb_sf lb_sf = 3, tr_sf = 4, dl_sf = 5, lg_sf = 6,
    ### calculate fraction of fuels in each class
   # sumfuel=4.0
    
    #'tw_sf','sb_sf','lb_sf','tr_sf'
    #coarsefuel=Vegframe[Vegframe.TissueType=='tw_sf'].Biomass.values+Vegframe[Vegframe.TissueType=='sb_sf'].Biomass.values+Vegframe[Vegframe.TissueType=='lb_sf'].Biomass.values+Vegframe[Vegframe.TissueType=='tr_sf'].Biomass.values
    #fuel_dl_sf=agcwd/sumfeul
    #fuel_dl_sf=Vegframe[Vegframe.TissueType=='deadleaves'].Biomass.values/sumfuel
    #fuel_lg= Vegframe[Vegframe.TissueType=='livegrass'].Biomass.values/sumfuel
    FBD=0
    MEF=0
    SAV=0
    fuel_moisture=0
    for tt in Vegframe.TissueType:
      #  print(tt)
        Vegframe.loc[(Vegframe.TissueType==tt),'Fuelfrac']=Vegframe[Vegframe.TissueType==tt].Biomass.values/sumfuel
        Vegframe.loc[(Vegframe.TissueType==tt),'MEF']= 0.524 - 0.066 * np.log(NFSC[NFSC.TissueType==tt].SAV.values)
        Vegframe.loc[(Vegframe.TissueType==tt),'alphaFMC']=NFSC[NFSC.TissueType==tt].SAV.values/SF_val_drying_ratio
        Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture']=np.exp(-1*Vegframe.loc[(Vegframe.TissueType==tt),'alphaFMC']*Site_Ni)
        if(tt=='livegrass'):
                 Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture']=np.exp(-1*Vegframe.loc[(Vegframe.TissueType=='deadleaves'),'alphaFMC']*Site_Ni).values
        if(tt!='tr_sf'):
            FBD=FBD+(Vegframe.loc[(Vegframe.TissueType==tt),'Fuelfrac'].values*NFSC[NFSC.TissueType==tt].FBD.values)
            SAV=SAV+(Vegframe.loc[(Vegframe.TissueType==tt),'Fuelfrac'].values*NFSC[NFSC.TissueType==tt].SAV.values)
            MEF=MEF+(Vegframe.loc[(Vegframe.TissueType==tt),'Fuelfrac'].values*Vegframe.loc[(Vegframe.TissueType==tt),'MEF'].values)
            fuel_moisture=fuel_moisture+(Vegframe.loc[(Vegframe.TissueType==tt),'Fuelfrac'].values*Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture'].values)
        Vegframe.loc[(Vegframe.TissueType==tt),'litter_moisture']=Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture']/Vegframe.loc[(Vegframe.TissueType==tt),'MEF'].values
    return(FBD,SAV,MEF,fuel_moisture,sumfuel,Vegframe)  




def ground_fuel_consumption(Vegframe,###DataFrame of tissues
                            NFSC, ### Params
                            sumfuel):      
    for tt in Vegframe.TissueType:  
       # print(tt)                   
        if Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture'].values <=NFSC[NFSC.TissueType==tt].min_moisture.values:
            Vegframe.loc[(Vegframe.TissueType==tt),'BF']=1.0 ### All fuel in pool is consumed
        if (Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture'].values>=NFSC[NFSC.TissueType==tt].min_moisture.values) and (Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture'].values <= NFSC[NFSC.TissueType==tt].mid_moisture.values):
            Vegframe.loc[(Vegframe.TissueType==tt),'BF']=min(1.0,NFSC[NFSC.TissueType==tt].low_moisture_Coeff.values- 
                  NFSC[NFSC.TissueType==tt].low_moisture_Slope.values*Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture'].values) 
        if Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture'].values>=NFSC[NFSC.TissueType==tt].mid_moisture.values and Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture'].values <=1.0 :
             Vegframe.loc[(Vegframe.TissueType==tt),'BF']=min(1.0,NFSC[NFSC.TissueType==tt].mid_moisture_Coeff.values- 
                  NFSC[NFSC.TissueType==tt].mid_moisture_Slope.values*Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture'].values)
        if Vegframe.loc[(Vegframe.TissueType==tt),'fuelmoisture'].values > 1.0:
            Vegframe.loc[(Vegframe.TissueType==tt),'BF']=0.0
        if(tt!='tr_sf'):
            Vegframe.loc[(Vegframe.TissueType==tt),'tau_b']=39.4*(Vegframe.loc[(Vegframe.TissueType==tt),'Fuelfrac']*sumfuel/0.45/19)*(1-((1-Vegframe.loc[(Vegframe.TissueType==tt),'BF'])**0.5))
    
    TFC_ROS=sum(Vegframe["BF"]*Vegframe['Biomass'])-(Vegframe.loc[(Vegframe.TissueType=='tr_sf'),'BF']*Vegframe.loc[(Vegframe.TissueType=='tr_sf'),'Biomass'])
    ### Calculate residence time (cap at 8 mins)
    tau_si1=min(8.0,sum(Vegframe['tau_b']))
    return(Vegframe,TFC_ROS,tau_si1)



def area_burnt_intensity(SF_val_fdi_alpha, ##relative importance of the nestrov index 
                         Site_Ni, ## Current Nesterov index
                         cg_strikes, ##Ration of cloud to ground strikes
                         ED_val_nignitions, ###Valid ignitions
                         SFval_max_durat, ### Max duration of fires
                         SF_val_durat_slope, ### relationship between fire weather and duration of fire
                         effect_wspeed, ### effective windspeed. 
                         tree_fraction_patch, ## Fraction of a patch that is tree
                         ROS_front, ## Rate of spread forward m/min
                         ROS_back, ## Rate of spread forward m/min
                         TFC_ROS, ### total fuel consumed by rate of spread. KgC/m2
                         SF_val_fuel_energy ### the energy value of fuel kJ/kg
                         ): 
    forest_grassland_lengthtobreadth_threshold = 0.55
    FDI=1.0-np.exp(-SF_val_fdi_alpha*Site_Ni)    ### Fuel~weather relaitionship
    ### strikes
    NF = ED_val_nignitions * 365* cg_strikes ##Gross number of ignitions* 1/365 *how many cloud to ground strikes
    ### Fire duration calculation
    Fire_duration=(SFval_max_durat+1)/(1+SFval_max_durat*np.exp(SF_val_durat_slope*FDI))
    #### The effective of fuel type of fire movement: 
    if (effect_wspeed*16.67<.1): #16.67m/min = 1km/hr
        lb=1
    ### Forest fuels    
    else:
        if(tree_fraction_patch > forest_grassland_lengthtobreadth_threshold):
            lb=((1.0+8.729*((1.0-(np.exp(-0.03*16.67*effect_wspeed)))**2.155)))
    ### Grass fuels 
        else:
            lb=1.1*((16.67*effect_wspeed)**0.464)
    ###Lenght of the major axis of fire spread
    db=ROS_back*Fire_duration #m
    df=ROS_front*Fire_duration #m 
    if lb>0.0:
        size_of_fire = ((3.1416/(4.0*lb))*((df+db)**2.0))
        AB =size_of_fire*FDI*NF ### Area burned is equal to size of fire by fuel dryness index and number of fires in a day. 
        Fractionburned=min((0.99,AB/1000000.0))## in KM ### return to line 757
    else: 
        Fractionburned=0.0
    ### Calculating intensity 
    ROS=ROS_front/60 #m/min to m/sec
    W=TFC_ROS/0.45 #kgC/m2 to kgbiomass/m2 of burned area. 
    Fire_intensity=SF_val_fuel_energy*W*ROS*Fractionburned
   # print(Fractionburned,Fire_intensity,AB)
    return(Fractionburned,Fire_intensity,AB)










