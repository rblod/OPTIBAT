#!/opt/local/bin/python

import os,sys
from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib import rcParams
from random import randrange
import scipy.io

#- offshore side (west preferred ...)
side='west'

#- input data :
#  depth (positive)
#  tide   (anomaly relative to reference level)
#  hs 
#  period, dir (degree, relative to the local mark, trigonometric convetion) 
datafile='Bathy_Evol.mat'


#- load datafile
mat = scipy.io.loadmat(datafile)

#- duplicate data to have 5 points in lonshore direction
X=np.repeat(mat["X"],5,0)

sizeX=X.shape[1]
sizeY=X.shape[0]

#- load hs on 5 points
hs=np.repeat(np.transpose(mat["hs"][:,:]),sizeY,1)
#- load period on 5 points 
Tp=np.repeat(mat["Tp"][:,:],sizeY,1)
Tp=np.repeat(Tp[:,:],31,0)

#- load depth and revert direction to have 
dep=np.repeat(mat["dep"][:,np.newaxis,:],sizeY,1)
dep=dep[:,:,::-1]

#- load tide, here 0
Tide=np.zeros((31,sizeY))*1.

#- load direction, here 0
Dire=Tide

#plt.plot(hs[:,1])
#plt.show()

#- grid data with depth (h) and coordinates (x_rho and eta_rho)
for ll in range(31):
	fname='depth_file_'+str(ll).zfill(2)+'.nc'	
	f = netcdf(fname,'w', format='NETCDF4', clobber=True) #'w' stands for write
	f.createDimension('xi_rho', X.shape[1])
	f.createDimension('eta_rho', X.shape[0])
	f.createDimension('time_counter', None)
	
	longitude = f.createVariable('xi_rho', 'f4', 'xi_rho')
	latitude = f.createVariable('eta_rho', 'f4', 'eta_rho')  
	time_counter = f.createVariable('time_counter', 'i4', 'time_counter')

	h = f.createVariable('h', 'f4', ('eta_rho', 'xi_rho'))

	f.variables['xi_rho'][:]=X[1,:]
	f.variables['eta_rho'][:]=X[1,0:5]
	f.variables['h'][:,:]=dep[ll,:,:]
	f.close()

#- boundary forcing data with tide, hs, period, dir (with suffix correponding to side
#- time counter is number of records (1,2,3 etc), not time vector
for ll in range(31):
	fname='bryfile_'+str(ll).zfill(2)+'.nc'	
	f = netcdf(fname,'w', format='NETCDF4', clobber=True) #'w' stands for write
	f.createDimension('xi_rho', X.shape[1])
	f.createDimension('eta_rho', X.shape[0])
	f.createDimension('time_counter', None)
	
	longitude = f.createVariable('xi_rho', 'f4', 'xi_rho')
	latitude = f.createVariable('eta_rho', 'f4', 'eta_rho')  
	time_counter = f.createVariable('time_counter', 'i4', 'time_counter')

	tide = f.createVariable('tide_'+side, 'f4', ('time_counter', 'eta_rho'))
	hsid = f.createVariable('hs_'+side, 'f4', ('time_counter', 'eta_rho'))
	period = f.createVariable('period_'+side, 'f4', ('time_counter', 'eta_rho'))
	cdir = f.createVariable('dir_'+side, 'f4', ('time_counter', 'eta_rho'))

	f.variables['eta_rho'][:]=X[1,0:5]
	f.variables['tide_'+side][0,:]=Tide[ll,:]
	f.variables['hs_'+side][0,:]=hs[ll,:]
	f.variables['period_'+side][0,:]=Tp[ll,:]
	f.variables['dir_'+side][0,:]=Dire[ll,:]
	f.close()

