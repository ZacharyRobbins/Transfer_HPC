# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 09:34:36 2022

@author: 345578
"""
import numpy as  np
import matplotlib.pyplot as plt
SAV=np.array(list(range(0,55)))
MEF= 0.524 - 0.066 * np.log(SAV) 
plt.plot(SAV,MEF)

SF_val_drying_ratio=70000

##Livegrass moisture


out1= 2.71828**(-1.0 * (SAV/SF_val_drying_ratio) * 3000)   
out2= 2.71828**(-1.0 * (SAV/SF_val_drying_ratio) * 5000)   
out3= 2.71828**(-1.0 * (SAV/SF_val_drying_ratio) * 10000)   


plt.plot(SAV,MEF)
plt.plot(SAV,out1)
plt.plot(SAV,out2)
plt.plot(SAV,out3)

#2.71828**(-1.0 * (SAV/SF_val_drying_ratio) * 3000)   

SAV=40.7
NI=np.array(list(range(0,10000)))
MEF= 0.524 - 0.066 * np.log(SAV) 
print(MEF)
plt.plot(SAV,MEF)
out1= 2.71828**(-1.0 * (SAV/66000)*NI)   
out2= 2.71828**(-1.0 * (SAV/250000) * NI)   
out3= 2.71828**(-1.0 * (SAV/500000) * NI)   

plt.ylim(0,1)
plt.plot(NI,out1,label="low SF_val_drying_ratio")
plt.plot(NI,out2,label="med SF_val_drying_ratio")
plt.plot(NI,out3,label="high SF_val_drying_ratio")
plt.ylabel("Moisture %")
plt.xlabel("Nesterov Index")
plt.axhline(y=MEF, color='r', linestyle='-',label="Moisture of Exhaustion")
plt.legend()


SAV=126
NI=np.array(list(range(0,10000)))
MEF= 0.524 - 0.066 * np.log(SAV) 
print(MEF)
out1= 2.71828**(-1.0 * (SAV/50000)*NI)   
out2= 2.71828**(-1.0 * (SAV/250000) * NI)   
out3= 2.71828**(-1.0 * (SAV/500000) * NI)   
plt.ylim(0,1)
plt.ylabel("Moisture %")
plt.xlabel("Nesterov Index")
plt.plot(NI,out1,label="low SF_val_drying_ratio")
plt.plot(NI,out2,label="med SF_val_drying_ratio")
plt.plot(NI,out3,label="high SF_val_drying_ratio")
plt.axhline(y=MEF, color='r', linestyle='-',label="Moisture of Exhaustion")
plt.legend()

