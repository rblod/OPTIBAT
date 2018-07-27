#!/usr/bin/env python
import numpy as np
import random
import netCDF4
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import griddata


def ncf_get_var(fname,var):
        nc= netCDF4.Dataset(fname)
        myvar=nc.variables[var][0,0,:]
        return myvar

def ncf_get_var_t(fname,var,indx):
        nc= netCDF4.Dataset(fname)
        myvar=nc.variables[var][indx,0,:]
        return myvar

niter=int(sys.argv[1])


fig=plt.figure(1)
plt.subplot(211)
for i in range(niter):
	fname='for_avg_'+str(i+1)+'.nc'
	h=   ncf_get_var(fname,'h')
	plt.plot(-h,label=str(i+1))

leg = plt.legend(loc='best', ncol=1, mode=None, shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)	


plt.subplot(212)
for i in range(niter):
	fname='for_std_'+str(i+1)+'.nc'
	h=   ncf_get_var(fname,'h')
	plt.plot(-h,label=str(i+1))

leg = plt.legend(loc='best', ncol=1, mode=None, shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)	

plt.show()     

num=10

fig=plt.figure(2)
ax = fig.add_subplot(111)
for i in range(100):
	h=ncf_get_var_t('for_cat_'+str(num)+'.nc','h',i)
	Cg=ncf_get_var_t('for_cat_'+str(num)+'.nc','Cg',i)
	plt.plot(-h,Cg,'k+')
ax.set_xlabel('h')	
ax.set_ylabel('Cg')	
ax.set_title('at step '+str(num))

fig=plt.figure(3)
ax = fig.add_subplot(111)
for i in range(100):
	h=ncf_get_var_t('for_cat_'+str(num)+'.nc','h',i)
	epb=ncf_get_var_t('for_cat_'+str(num)+'.nc','epb',i)
	plt.plot(-h,epb,'k+')
ax.set_xlabel('h')	
ax.set_ylabel('epb')	
ax.set_title('Forcast at step '+str(num))


plt.show()     