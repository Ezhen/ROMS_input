import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from shutil import copyfile
from bathy_smoother import *

#parent_Imin = 16 ;
#parent_Imax = 43 ;
#parent_Jmin = 35 ;
#parent_Jmax = 66
#:parent_Imin = 13 ;
#:parent_Imax = 41 ;
#:parent_Jmin = 21 ;
#:parent_Jmax = 49 ;


rr = Dataset('/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Parent_river_boundary_6m.nc', 'r', format='NETCDF4')
h = rr.variables['h'][:]
msk = rr.variables['mask_rho'][:]
rr.close()

copyfile('/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Child_river_boundary_update_bath.nc', '/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Child_river_boundary_update_bath_smoothed_boundary.nc')
nc = Dataset('/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Child_river_boundary_update_bath_smoothed_boundary.nc', 'a', format='NETCDF4')
hf = nc.variables['h'][:]
mask_rho = nc.variables['mask_rho'][:]

RoughMat = bathy_tools.RoughnessMatrix(hf, mask_rho)
print 'Max Roughness value is: ', RoughMat.max()

print 'smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)'
rx0_max = 0.35
hff = bathy_smoothing.smoothing_Positive_rx0(mask_rho, hf, rx0_max)

RoughMat = bathy_tools.RoughnessMatrix(hff, mask_rho)
print 'Max Roughness value is: ', RoughMat.max()

mask_rho[hff[5:-6,5:-6] < 4] = 0
hff[hff < 4] = 4

nc.variables['h'][:] = hff
nc.variables['mask_rho'][5:-6,5:-6] = mask_rho[5:-6,5:-6]
"""
for i in range(len(hf)):
	for j in range(len(hf.T)):
		if i<5 or len(hf)-i<5 or j<5 or len(hf.T)-j<5:
			if nc.variables['mask_rho'][i,j]==1:
				nc.variables['h'][i,j] = h[21+i/5,13+j/5]
				#nc.variables['mask_rho'][i,j] = msk[21+i/5,13+j/5]
"""
mask_rho[:] = nc.variables['mask_rho'][:]
nc.variables['mask_u'][:] = mask_rho[:,1:]*mask_rho[:,:-1]
nc.variables['mask_v'][:] = mask_rho[1:,:]*mask_rho[:-1,:]
nc.variables['mask_psi'][:] = mask_rho[1:,1:]*mask_rho[:-1,1:]*mask_rho[1:,:-1]*mask_rho[:-1,:-1]

nc.close()


