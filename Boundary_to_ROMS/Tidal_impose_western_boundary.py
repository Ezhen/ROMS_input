from netCDF4 import Dataset
import numpy as np
from shutil import copyfile
import pyroms

copyfile('/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Tides_parent_untouched.nc', '/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Tides_parent_touched.nc')

nc = Dataset('/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Tides_parent_touched.nc', 'a', format='NETCDF4')
dstgrd = pyroms.grid.get_ROMS_grid('ParentRiver')

nc.variables['mask_rho'][:] = dstgrd.hgrid.mask_rho
nc.variables['tide_Eamp'][:,28,0:2] = nc.variables['tide_Eamp'][:,37,13:15]
nc.variables['tide_Ephase'][:,28,0:2] = nc.variables['tide_Ephase'][:,37,13:15]
nc.variables['tide_Cphase'][:,28,0:2] = nc.variables['tide_Cphase'][:,37,13:15]
nc.variables['tide_Cangle'][:,28,0:2] = nc.variables['tide_Cangle'][:,37,13:15]
nc.variables['tide_Cmin'][:,28,0:2] = nc.variables['tide_Cmin'][:,37,13:15]
nc.variables['tide_Cmax'][:,28,0:2] = nc.variables['tide_Cmax'][:,37,13:15]

nc.close()
