from netCDF4 import Dataset; import numpy as np; from shutil import copyfile

copyfile('DIRECT_FLUXES.nc', 'DIRECT_FLUXES_MY_MOMENTUM.nc')

nc = Dataset('DIRECT_FLUXES_MY_MOMENTUM.nc', 'a', format='NETCDF4')
rr = Dataset('sustr_svstr_from_ROMS_simulation.nc', 'r', format='NETCDF4')
for i in range(len(rr.variables['ocean_time'][:])/3):
	nc.variables['sustr'][i,:,:-1]=rr.variables['sustr'][i*3]
	nc.variables['svstr'][i,:-1,:]=rr.variables['svstr'][i*3]
	nc.variables['sustr'][i,:,-1]=nc.variables['sustr'][i,:,-2]
	nc.variables['svstr'][i,-1,:]=nc.variables['svstr'][i,-2,:]

nc.close()
rr.close()
