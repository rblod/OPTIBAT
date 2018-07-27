
#!/usr/bin/env python
import numpy as np
import random
import netCDF4
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import griddata


def ncf_get_var(fname,var):
	nc= netCDF4.Dataset(fname)
	myvar=nc.variables[var][:]
	return myvar

def ncf_copy_file(fname1,fname2,z2):
	xi_rho=ncf_get_var(fname1,'xi_rho')
	eta_rho=ncf_get_var(fname1,'eta_rho')
	#
	nc2=netCDF4.Dataset(fname2,'w', format='NETCDF4', clobber=True)
	nc2.createDimension('xi_rho', np.size(xi_rho))
	nc2.createDimension('eta_rho',np.size(eta_rho))
	nc2.createDimension('time_counter', None)
	nc2.createVariable('xi_rho', 'f4', 'xi_rho')
	nc2.createVariable('eta_rho', 'f4', 'eta_rho')
	nc2.createVariable('time_counter', 'i4', 'time_counter')
	nc2.createVariable('h', 'f4', ('eta_rho', 'xi_rho'))
	nc2.variables['xi_rho'][:]=xi_rho
	nc2.variables['eta_rho'][:]=eta_rho
	nc2.variables['h'][:,:]=z2
	nc2.close()

def genbathy(X):

#	X=np.arange(60)
	Y=5

	# interval number bars
	ib=np.array([0,5])
	b=random.random()
	bnb=int(round((ib[1]-ib[0])*b + ib[0]))   #number of bars


	# slope range
	ib=np.array([0.01, 0.05])
	b=random.random()
	sl=((ib[1]-ib[0])*b + ib[0])       #linear beach slope

	# bar amplitude range
	ib=np.array([0.25, 6])
	b=np.random.rand(bnb,1)
	bamp=sorted(((ib[1]-ib[0])*b+ ib[0]))

	# bar position range
	b=np.random.rand(bnb,1)
	ib=np.array([1*np.size(X)/5., 3*np.size(X)/5.])
	bx=sorted(np.round((ib[1]-ib[0])*b+ ib[0])) 
	bx = [int(i) for i in bx]

	# bar width range
	ib=np.array([50, 600])
	b=np.random.rand(bnb,1)
	bw=sorted(np.round((ib[1]-ib[0])*b + ib[1])) 
	bw = [int(i) for i in bw]

	z=-X*sl;    
	for ib in range (0,bnb) :
		bw[ib]=3.*bx[ib];
		while z[bx[ib]]+bamp[ib]*np.exp(  -((X[bx[ib]]-bx[ib])**2)/(bw[ib]) )>z[bx[ib]]/3 :
			bamp[ib]=bamp[ib]/2
		z=z+bamp[ib]*np.exp(  -((X-bx[ib])**2)/(bw[ib]) )
	return z

niter=int(sys.argv[1])
filein=sys.argv[2]
doplot=int(sys.argv[3])

X1=ncf_get_var(filein,'xi_rho')
X2=np.arange(X1.max())
yindx=np.size(ncf_get_var(filein,'eta_rho'))

for i in range(niter):
	z2=genbathy(X2)
	z1=griddata(X2,z2,X1)
	
	if  doplot == 1:
		if i%10 == 0 :
			plt.plot(X1,z1)
			plt.plot(X2,z2)
			plt.show()

	z1[0]=z1[1]
	z1[-1]=z1[-2]
	z1=-z1[::-1]
	z3=np.tile(z1,(yindx,1))
	fileout='bathy_'+str(i+1).zfill(3)+'.nc'
	nc=ncf_copy_file(filein,fileout,z3)


