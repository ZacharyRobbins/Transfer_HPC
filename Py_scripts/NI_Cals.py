# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 09:10:13 2022

@author: 345578
"""


import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from bisect import bisect_left
import geopandas as gpd
from functools import partial 
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
        else: 
          Site_Ni=Site_Ni+d_NI
    return(Site_Ni)
###
def NI_proc(DailyTmean,Precip,DailyRhmean,a,b):
   i=0
   Site_Ni=0
   NI_vec=0
   for i in list(range(0,len(Precip['value']))):
       Site_Ni=fire_danger_index(Site_Ni,
                              np.array(DailyTmean)[i],
                              np.array(Precip['value'])[i]*10,
                              np.array(DailyRhmean)[i],
                              a,
                              b)
       NI_vec=np.append(NI_vec,Site_Ni)
   return(NI_vec)

def:

#####################
### Datasets and cleaning ###################################
######


Dir='C:/Users/345578/Desktop/California_chaparral_project/Practice_Climate/'
Precip=pd.read_csv(Dir+"Eco1_CA_Chap_pr_Table.csv")
rhmax=pd.read_csv(Dir+"Eco1_CA_Chap_rmax_Table.csv")
rhmin=pd.read_csv(Dir+"Eco1_CA_Chap_rmin_Table.csv")
windspeed=pd.read_csv(Dir+"Eco1_CA_Chap_th_Table.csv")
winddirection=pd.read_csv(Dir+"Eco1_CA_Chap_vs_Table.csv")
tmmn=pd.read_csv(Dir+"Eco1_CA_Chap_tmmn_Table.csv")
tmmx=pd.read_csv(Dir+"Eco1_CA_Chap_tmmx_Table.csv")
DailyTmean=((tmmx['value']-273.15)+(tmmn['value']-273.15))/2
DailyRhmean=((rhmax['value']+rhmin['value'])/2)
Ign_Dir='C:/Users/345578/Desktop/California_chaparral_project/Validation_Data/Ignitions/'
Ign =gpd.read_file(Ign_Dir+"KS_subbed.gpkg")
nIgnperday=pd.DataFrame((Ign[['DISCOVERY_DATE', 'FIRE_SIZE']]
    .query('FIRE_SIZE >1')
     .groupby('DISCOVERY_DATE')
     ['FIRE_SIZE'].count()))
#nIgnperday=pd.DataFrame(nIgnperday).columns=['nIgn']
#nIgnperday['Date']=pd.to_datetime(nIgnperday['Date'])
nIgnperday['Dates']=pd.to_datetime(nIgnperday.index.values)
##########
### Setup a partial function.
### Shorten up Climate ins. 
NI_part=partial(NI_proc, DailyTmean=DailyTmean, Precip=Precip, DailyRhmean=DailyRhmean)
######
Outy=pd.DataFrame(columns=['a','b','No',"Low","Med","Max"])

runs=1000
for i in list(range(1,runs)):
    print(i)
    a=np.round(np.random.uniform(50,200))# 20:200
    b=np.round(np.random.uniform(500,1000))
    NI_out=NI_part(a=200,b=200)
    #plt.plot(NI_out)
    NI_df=pd.DataFrame({'Dates':pd.to_datetime(Precip.date),"NI":NI_out[1:]}).query('Dates> (19920101)')
    Joined=pd.merge(NI_df, nIgnperday, how='left', on='Dates')
    Joined['FIRE_SIZE']=Joined['FIRE_SIZE'].fillna(0)
    #plt.plot(Joined['FIRE_SIZE'],Joined['NI'])
    Joined=Joined.sort_values(by=['NI'])
    Joined['Cumsum']=Joined['FIRE_SIZE'].cumsum()
    
    No=  Joined.loc[bisect_left(Joined['Cumsum'],Joined['Cumsum'].quantile(0.05)),'NI'] ### Aiming for 300
    Low=  Joined.loc[bisect_left(Joined['Cumsum'],Joined['Cumsum'].quantile(0.20)),'NI'] ### Aiming for 1000
    Med=  Joined.loc[bisect_left(Joined['Cumsum'],Joined['Cumsum'].quantile(0.65)),'NI'] ### Aiming for 4000     
    Max=max(Joined['NI']) ### Aiming in the range of 10,000
    indy=pd.Series(1)
    Inny=pd.DataFrame({'a':a,'b':b,'No':No,"Low":Low,"Med":Med,"Max":Max},index=indy)
    Outy=Outy.append(Inny,ignore_index=True)


Outy['Ratio']=Outy['Max']+Outy['Med']
