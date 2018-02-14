import netCDF4; from shutil import copyfile; import numpy as np

copyfile('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Coarsest_improved.nc', 'Drag_map.nc')
rr = netCDF4.Dataset('Drag_map.nc', 'a', format='NETCDF4')
msk = rr.variables['mask_rho'][:]; h = rr.variables['h'][:]

#rr.createVariable('bttm_length','f4',('eta_rho','xi_rho',))
rr.createVariable('ZoBot','f4',('eta_rho','xi_rho',))


rr1 = netCDF4.Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Validation_tides/Tide_local.nc', 'r', format='NETCDF4')
rr2 = netCDF4.Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Validation_tides/UV/Major_Axis_diff.nc', 'r', format='NETCDF4')

rr.variables['ZoBot'][:]=0.0024

for i in range(1,len(msk)-1):
	for j in range(1,len(msk.T)-1):
		if msk[i,j]==1:
			if msk[i-1,j-1]==0 or msk[i-1,j]==0 or msk[i-1,j+1]==0 or msk[i,j-1]==0 or msk[i,j+1]==0 or msk[i+1,j-1]==0 or msk[i+1,j]==0 or msk[i+1,j+1]==0:
				rr.variables['ZoBot'][i,j]=0.75


			"""
			a=rr2.variables['tide_Cmax'][0,i,j]/rr1.variables['tide_Cmax'][0,i,j]+1
			print i,j,a
			if a<1.0001:
				rr.variables['ZoBot'][i,j]=0.0024
			else:
				rr.variables['ZoBot'][i,j]=0.0024*a**5.2+0.002
			#rr.variables['bttm_length'][i,j]=(0.05*(1.0+np.tanh(h[i,j]/50.0)))*1000.
			"""
rr.close(); rr1.close(); rr2.close()
