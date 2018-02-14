import numpy as np; import scipy; from netCDF4 import Dataset; from bathy_smoother import *; from shutil import copyfile

copyfile('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Grid_Nesting/Nested.nc', 'Nested_grid.nc')
rr = Dataset('Nested_grid.nc', 'a', format='NETCDF4')
h = rr.variables['h'][:]; mask_rho=rr.variables['mask_rho'][:]

nc = Dataset('Nested_without_refinement.nc', 'r', format='NETCDF4');mask_rho2=nc.variables['mask_rho'][:]
for i in range(len(h)):
	for j in range(len(h.T)):
		if i<10 or abs(i-len(h))<11 or j<10 or abs(j-len(h.T))<11:
			mask_rho[i,j]=mask_rho2[i,j]

for i in range(len(h)):
	for j in range(len(h.T)):
		if i<15 or abs(i-len(h))<16 or j<15 or abs(j-len(h.T))<16:
			mint = min(i,j,abs(i-len(h))-1, abs(j-len(h.T))-1)
			if mint<5:
				h[i,j]=nc.variables['h'][i,j]
			else:
				mint=mint-5
				h[i,j]=(0.1*mint)*h[i,j]+(1-0.1*mint)*nc.variables['h'][i,j]
nc.close()

rr.variables['h'][:]=h[:]
rr.variables['mask_rho'][:]= mask_rho[:]
rr.variables['mask_u'][:] = mask_rho[:,1:]*mask_rho[:,:-1]
rr.variables['mask_v'][:] = mask_rho[1:,:]*mask_rho[:-1,:]
rr.variables['mask_psi'][:] = mask_rho[1:,1:]*mask_rho[:-1,1:]*mask_rho[1:,:-1]*mask_rho[:-1,:-1]
#rr.close()

ncdata3 = Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Coarsest_improved.nc', 'r', format='NETCDF3')
hp = ncdata3.variables['h'][:]; mask_rhop = ncdata3.variables['mask_rho'][:]

#plot(32+arange(len(h[0,::5])),h[0,::5]); plot(arange(len(hp[0,:])),hp[0,:])

#plot(35+arange(len(h[0,::5])),h[0,::5]); plot(arange(len(hp[35,:])),hp[35,:])

#plot(32+arange(len(mask_rho[0,::5])),mask_rho[0,::5]); plot(arange(len(mask_rhop[0,:])),mask_rhop[0,:])

RoughMat = bathy_tools.RoughnessMatrix(h, mask_rho)
print 'Max Roughness value is: ', RoughMat.max()

rx0_max = 0.35
depth = bathy_smoothing.smoothing_Positive_rx0(mask_rho, h, rx0_max)

RoughMat = bathy_tools.RoughnessMatrix(depth, mask_rho)
print 'Max Roughness value is: ', RoughMat.max()

rr.variables['h'][:]=depth
rr.close()

