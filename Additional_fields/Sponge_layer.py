import netCDF4; from shutil import copyfile; import numpy as np

copyfile('/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Parent_river_boundary_6m.nc', '/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Parent_river_boundary_6m_sponge.nc')
rr = netCDF4.Dataset('/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Parent_river_boundary_6m_sponge.nc', 'a', format='NETCDF4')
msk = rr.variables['mask_rho'][:]; h = rr.variables['h'][:]

rr.createVariable('visc_factor','f4',('eta_rho','xi_rho',))
rr.variables['visc_factor'].long_name = "horizontal viscosity sponge factor"
rr.variables['visc_factor'].valid_min = 0.
rr.variables['visc_factor'].coordinates = "lon_rho lat_rho"

rr.createVariable('diff_factor','f4',('eta_rho','xi_rho',))
rr.variables['diff_factor'].long_name = "horizontal diffusivity sponge factor"
rr.variables['diff_factor'].valid_min = 0.
rr.variables['diff_factor'].coordinates = "lon_rho lat_rho"

m=[]
for i in range(len(msk)):
	for j in range(len(msk.T)):
		if i<5 or j<5 or len(msk)-i<6 or len(msk.T)-j<6:
			a = 5 - min(list([i,j,len(msk)-i,len(msk.T)-j]))
			rr.variables['visc_factor'][i,j] = a
			rr.variables['diff_factor'][i,j] = a
		else:
			rr.variables['visc_factor'][i,j] = 1
			rr.variables['diff_factor'][i,j] = 1
rr.close()
