#!/opt/local/bin/python
# a sample script to create different forcing and depth file for a shoreface like

import os,sys
from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib import rcParams
from random import randrange
import scipy.io
from math import *

# domaine
xsize=61     #  number of points in x-direction
ysize=5      # number of points in y direction, at least 3
dx=20.       # resolution (m)
side='west'  # boundary side (left)

# constant wave data
wkb_hs = 0.7*2/0.7         
wkb_ang = -10.0    # degree relative to the shore normal, trigonometric convention  
wkb_prd = 8.3      # second

# variable tide data
time=6*3600         # 6 hours in seconds
nstp=time/(5*60)+1    # 1 sampl;e each 5mn 
wkb_tide=np.zeros(nstp)
for i in range(nstp):
	wkb_tide[i]=cos(i*pi/nstp)


# bathymetry      # axis down, ie positive when wet 
xs=350      # inner surf zone
db=80       # distance from xs to sand bar
eps=2.
nbathy=4    # here 4 sets of bathy as an example
h=np.zeros([nbathy,ysize,xsize])

for i in range(xsize):
	h0=12.0-(0.0125*dx*i)        # original fonction
	h[0,:,i]=h0
	h[1,:,i]=12.0-(0.0108*dx*i)  #original fonction down
	h[2,:,i]=12.0-(0.014*dx*i)   #original fonction up
	xx=dx*(xsize-1)-i*dx
	per=eps*exp(-(((xx-xs)/db)**2))  #add a bar
	h[3,:,i]=h0-per


# fill bathy files
for nfile in range(4):
	fname='depth_file_'+str(nfile)+'.nc'
	f = netcdf(fname,'w', format='NETCDF4', clobber=True) #'w' stands for write
	f.createDimension('xi_rho', xsize)
	f.createDimension('eta_rho', ysize)
	f.createDimension('time_counter', None)
	longitude = f.createVariable('xi_rho', 'f4', 'xi_rho')
	latitude = f.createVariable('eta_rho', 'f4', 'eta_rho')
	time_counter = f.createVariable('time_counter', 'i4', 'time_counter')
	hh = f.createVariable('h', 'f4', ('eta_rho', 'xi_rho'))

	xi_rho=np.zeros(xsize)
	for i in range(xsize):
		xi_rho[i]=(i*dx)

	eta_rho=np.zeros(ysize)
	for i in range(ysize):
		eta_rho[i]=(i*dx)

	f.variables['xi_rho'][:]=xi_rho    #np.tile(xi_rho,(ysize,1))
	f.variables['eta_rho'][:]=eta_rho  #np.tile(eta_rho,(1,xsize))
	f.variables['h'][:,:]=h[nfile,:,:]
	f.close()

# boundary file, we have to fill al the boundary even for duplicate values
fname='shoreface_bry.nc'
f = netcdf(fname,'w', format='NETCDF4', clobber=True) #'w' stands for write
f.createDimension('xi_rho', xsize)
f.createDimension('eta_rho', ysize)
f.createDimension('time_counter', None)
longitude = f.createVariable('xi_rho', 'f4', 'xi_rho')
latitude = f.createVariable('eta_rho', 'f4', 'eta_rho')
time_counter = f.createVariable('time_counter', 'i4', 'time_counter')

xi_rho=np.zeros(xsize)
for i in range(xsize):
	xi_rho[i]=(i*dx)

eta_rho=np.zeros(ysize)
for i in range(ysize):
	eta_rho[i]=(i*dx)

tide = f.createVariable('tide_'+side, 'f4', ('time_counter', 'eta_rho'))
hsid = f.createVariable('hs_'+side, 'f4', ('time_counter', 'eta_rho'))
period = f.createVariable('period_'+side, 'f4', ('time_counter', 'eta_rho'))
cdir = f.createVariable('dir_'+side, 'f4', ('time_counter', 'eta_rho'))

f.variables['xi_rho'][:]=xi_rho    #np.tile(xi_rho,(ysize,1))
f.variables['eta_rho'][:]=eta_rho  #np.tile(eta_rho,(1,xsize))
for i in range(nstp):
	f.variables['time_counter'][i]=nstp
	f.variables['tide_'+side][i,:]=wkb_tide[i]
	f.variables['hs_'+side][i,:]=wkb_hs
	f.variables['period_'+side][i,:]=wkb_prd
	f.variables['dir_'+side][i,:]=wkb_ang
f.close()









