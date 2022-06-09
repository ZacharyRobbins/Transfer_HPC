# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:12:27 2022

@author: 345578
"""


# The threshold values of these danger classes were determined by comparing fire
# occurrence to the index values over a period of 10 years. For instance, 
# the "no fire danger" class upper limit (300) corresponds to the index values 
# during the reference period under which there was no fire occurrence.
# The "low fire danger" class upper limit (1000) corresponds to the index values
# under which 25 % of all fires occurred, and the "medium fire danger" class upper 
# limit (4000) to the index values under which 65 % of all fires occurred. 
# According to this method, the Nesterov index can be adjusted for every region 
# under consideration (Chandler et al. 1983).

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
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
#####################

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


i=0
Site_Ni=0
NI_vec=0
for i in list(range(0,len(Precip['value']))):
    Site_Ni=fire_danger_index(Site_Ni,
                              np.array(DailyTmean)[i],
                              np.array(Precip['value'])[i]*10,
                              np.array(DailyRhmean)[i],
                              80.5,
                              500)
    NI_vec=np.append(NI_vec,Site_Ni )
   
#plt.plot(np.array(Precip['value']))
together=pd.DataFrame({'NI':NI_vec[1:]},index=pd.DatetimeIndex(Precip.date))
g=together.groupby(pd.Grouper(freq="y"))
Annual=g.mean()
np.array(Precip.date)
plt.plot(together['NI'], color='black')   
plt.axhline(y=1000, color='r', linestyle='-',label="isture of Exhaustion")
plt.axhline(y=1500, color='r', linestyle='-',label="Moisture of Exhaustion")
plt.axhline(y=2000, color='r', linestyle='-',label="Moisture of Exhaustion")


np.quantile(NI_vec,.50)
np.quantile(NI_vec,.90)
np.quantile(NI_vec,.997)
#plt.plot(np.array(Precip['value']))
