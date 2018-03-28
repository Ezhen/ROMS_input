from netCDF4 import Dataset
import numpy as np
from shutil import copyfile
import pyroms

copyfile('/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/ParentRiver_Initial.nc', '/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/ParentRiver_Initial_western_boundary.nc')

nc = Dataset('/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/ParentRiver_Initial_western_boundary.nc', 'a', format='NETCDF4')
rr = Dataset('/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Boundary_ParentRiver_100_days_western_boundary.nc', 'r', format='NETCDF4')

nc.variables['ubar'][0,:,0] = rr.variables['ubar_west'][0]
nc.variables['u'][0,:,:,0] = rr.variables['u_west'][0]

rr.close()
nc.close()
