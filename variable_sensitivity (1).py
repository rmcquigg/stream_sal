# -*- coding: utf-8 -*-
"""
Created on Sun May 23 16:44:57 2021

@author: apaldor, rmcquiggan
"""
#This bit below is Anner's code
import numpy as np
import pandas as pd
#Define coefficients
rho_s = 1025
d_rho = 25
g = 9.8
v0 = 1.2 
T = 44280
H0=3.35
B0 = 795
h0=8.5
A0 = 6241
Qf=224
a=33333
b=a

def flhws(v0,H0,B0,h0,A0,Qf,a,rho_s=1025,d_rho=25,T=44280,b=a):
    g=9.8
    K = 1.6e-7*(h0**(0.69)*g**(1.12)*T**(2.24)/(H0**(0.59)*b**(1.1)*B0**(0.13)))
    E0 = v0 * T / np.pi
    NR = d_rho/rho_s * g*h0*Qf*T / (A0*E0*v0**2)
    D0 = 1400*(E0*v0*h0)/a*(np.sqrt(NR))
    beta = (K*a*Qf) / (D0*A0)
    LHWS = (a*np.log((1/beta) + 1))/1000
    return round(LHWS,2)

#LHWS=flhws(v0=v0,H0=H0,B0=B0,h0=h0,A0=A0,Qf=Qf,a=a,T=T)
#print(LHWS,'km')

#_________________________________________
#RWM
#Calibrate by variable
#You can change the values in the % list
#Variables go by whatver you defined above
caldf=pd.DataFrame(columns=['%'])
caldf['%']=[10,25,50,75,100,125,150,175,190]
vardict={'v0':v0,'H0':H0,'B0':B0,'h0':h0,'A0':A0,'Qf':Qf,'a':a}
vardict_new={'v0':v0,'H0':H0,'B0':B0,'h0':h0,'A0':A0,'Qf':Qf,'a':a}
fit=[]

for l in list(vardict)[0:]:
    caldf[str(l)+'_new']=vardict[l]*(caldf['%']/100)
    var=[]
    for v in caldf[str(l)+'_new']:
        vardict_new[l]=v
        res=flhws(**vardict_new)
        var.append(res)
        vardict_new[l]=vardict[l]
    caldf[str(l)]=var
    trend=np.polyfit(np.log(caldf[str(l)+'_new']),caldf[str(l)],1)
    fit.append("".join([str(round(trend[0],2)),'/x']))
    
rates=pd.DataFrame(columns=['vars','rate'])
rates['rate']=fit
rates['vars']=list(vardict)[0:]
#list(vardict.values())[0:]
    
#Copy df to a csv on shared drive
caldf.to_csv('G:/Shared drives/Pumping-Induced-Streamflow-Depletion/Writing/solution_cal_by_parameter_REV1.csv',index=False)
rates.to_csv('G:/Shared drives/Pumping-Induced-Streamflow-Depletion/Writing/variable_rate_of_change.csv',index=False)
