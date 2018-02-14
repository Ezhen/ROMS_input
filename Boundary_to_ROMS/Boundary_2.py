from netCDF4 import Dataset
import numpy as np

#nc = Dataset('Forcing.nc', 'a', format='NETCDF4')
nc = Dataset('/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/Boundary.nc', 'a', format='NETCDF4')

a=nc.variables['zeta_time'][:]

nc.variables['v3d_time'][:] = nc.variables['v3d_time'][:]/86400
nc.variables['v2d_time'][:] = nc.variables['v2d_time'][:]/86400
nc.variables['salt_time'][:] = nc.variables['salt_time'][:]/86400
nc.variables['temp_time'][:] = nc.variables['temp_time'][:]/86400
nc.variables['zeta_time'][:] = nc.variables['zeta_time'][:]/86400

nc.close()

