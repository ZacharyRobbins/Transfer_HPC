# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 09:51:52 2022

@author: 345578
"""



import pandas as pd 
import numpy as np
import geopandas as gpd
import os 
import matplotlib.pyplot as plt

#filter( lambda x:x ==  os.listdir('C:/Users/345578/Documents/GitHub/Transfer_HPC/SF_Wrap/SF_runs614/')
#any('run' in x for x in os.listdir('C:/Users/345578/Documents/GitHub/Transfer_HPC/SF_Wrap/SF_runs614/'))

os.chdir('C:/Users/345578/Documents/GitHub/Transfer_HPC/SF_Wrap/SF_runs614/')


Runs_list = [nm  for nm in os.listdir() if "run" in nm]
Params_list= [nm  for nm in os.listdir() if "params" in nm]
Out=pd.DataFrame(columns=["FB"])

for r in list(range(0,1118)):
    first=pd.read_csv(Runs_list[r])
    second=pd.read_csv(Params_list[r])
    In=pd.DataFrame({r:first['Fractionburned'].values},index=pd.DatetimeIndex(first['Date'])).groupby(pd.Grouper(freq="y")).sum()
    if r ==0:
        Out=In
        Out2=second
    else:
        Out=pd.concat([In,Out],1,ignore_index=True)
        Out2=pd.concat([second,Out2],0,ignore_index=True)


RunTotals=1/pd.DataFrame(Out.sum(axis=0)/10)
list(Out.columns)
#ut2.rename(columns={"Unnamed: 0"})
hah=Out2
reverse(RunTotals)
1/0.003

check=pd.concat([Out2,RunTotals],1,ignore_index=True)
check=check.rename(columns={0:"Run",1:'SF_val_drying_ratio',2:'mineral_total',3:'SF_val_fdi_alpha',
                      4:'SF_val_part_dens',5:'SF_val_durat_slope',6:"FRI"})



plt.xlim(0,1000)
plt.scatter(check['FRI'],check['SF_val_drying_ratio']) #(1000, 30000)
plt.axvline(82)
plt.xlim(0,10000)
plt.scatter(check['FRI'],check['mineral_total'])
plt.axvline(82)
plt.xlim(0,10000)
plt.scatter(check['FRI'],check['SF_val_fdi_alpha']) #(Same?)
plt.axvline(82)
plt.xlim(0,10000)
plt.scatter(check['FRI'],check['SF_val_part_dens']) #(go 10 1000)
plt.axvline(82)
plt.xlim(0,10000)
plt.scatter(check['FRI'],check['SF_val_durat_slope']) #(go, 0.5)
plt.axvline(82)
    


plt.hist(check['SF_val_drying_ratio'])


Ign_Dir='C:/Users/345578/Desktop/California_chaparral_project/Validation_Data/Ignitions/'
StudyArea=gpd.read_file(
    "C:/Users/345578/Desktop/California_chaparral_project/StudyArea_climate.shp"
   ).to_crs({'init': 'epsg:3857'})

#print(StudyArea.crs)
Area=StudyArea['geometry'].area
Ign =gpd.read_file(Ign_Dir+"KS_subbed.gpkg")