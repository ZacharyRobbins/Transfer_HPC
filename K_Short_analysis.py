# -*- coding: utf-8 -*-
"""
Created on Thu May 26 08:57:48 2022

@author: 345578
"""


import pandas as pd 
import numpy as np
import geopandas as gpd

Ign_Dir='C:/Users/345578/Desktop/California_chaparral_project/Validation_Data/Ignitions/'
StudyArea=gpd.read_file(
    "C:/Users/345578/Desktop/California_chaparral_project/Rough_Study_Area.gpkg"
    ).to_crs({'init': 'epsg:3857'})
#print(StudyArea.crs)
Area=StudyArea['geometry'].area
Ign =gpd.read_file(Ign_Dir+"Ign_SA.gpkg")

### Ignitions per year. 
nIgnbyYear=(Ign
    .query('FIRE_SIZE >1')
    .groupby('FIRE_YEAR')
    ['OBJECTID'].count()
    )
nIgnbyYear=pd.DataFrame(nIgnbyYear).rename(columns={'OBJECTID':'peryear'})
nIgnbyYear["perm2"]=nIgnbyYear['peryear']/Area.values
nIgnbyYear["perkm2"]=nIgnbyYear["perm2"]*10**6
list(Ign.columns)

#### Fire days per year
Ign['FireDays']=1+(Ign['CONT_DOY']-Ign['DISCOVERY_DOY'])
Ign['FireDays'].values[Ign['FireDays'].values < 0] = Ign['FireDays'].values[Ign['FireDays'].values < 0]+360

nFiredaysperYear=(Ign
    .query('FIRE_SIZE >1')
    .groupby('FIRE_YEAR')
    ['FireDays'].sum()
    )

nFiredaysperYear=pd.DataFrame(nFiredaysperYear).rename(columns={'FireDays':'peryear'})
nFiredaysperYear["perm2"]=nFiredaysperYear['peryear']/Area.values
nFiredaysperYear["perkm2"]=nFiredaysperYear["perm2"]*10**6


### Burned area per year per unit (Fire rotation period)

Area_Burned_Year=(Ign
    .query('FIRE_SIZE >1')
    .groupby('FIRE_YEAR')
    ['FIRE_SIZE'].sum()
    )
Area_Burned_Year=pd.DataFrame(Area_Burned_Year).rename(columns={'FIRE_SIZE':'Acperyear'})
Area_Burned_Year["haperyear"]=Area_Burned_Year['Acperyear']/2.54
Area_Burned_Year['FRP']=Area.values/(Area_Burned_Year["haperyear"]*10000)
meanFRP=(Area.values*27)/(sum(Area_Burned_Year["haperyear"])*10000)

### Fire size distribution ###
Ign['FIRE_SIZE_ha']=Ign['FIRE_SIZE']/2.54
bins = [0, 10, 50,100,500,1000,5000, 10000,50000,100000]
labels = [1,2,3,4,5,6,7,8,9]
Ign['FireClass'] = pd.cut(Ign['FIRE_SIZE_ha'], bins=bins, labels=labels)
Fire_size_Class=pd.DataFrame(Ign
    .query('FIRE_SIZE >1')
    .groupby('FireClass')
    ['FireClass'].count()
    )
Fire_size_Class['Percentage']=100*Fire_size_Class['FireClass']/sum(Fire_size_Class['FireClass'])
