from netCDF4 import Dataset
#import numpy as np

nc = Dataset('Forcing.nc', 'a', format='NETCDF4'); rr = Dataset('Nested_grid.nc', 'r', format='NETCDF4')
nc.variables['mask_rho'][:] =  rr.variables['mask_rho'][:]
nc.close(); rr.close()
