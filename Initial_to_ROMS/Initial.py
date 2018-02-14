from netCDF4 import Dataset
import numpy as np
from numpy import shape
import pyroms
import pyroms_toolbox
from datetime import datetime
import subprocess
import os
import commands
import _iso
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
from nco import Nco
import _remapping
from get_nc_Grid_HYCOM import get_nc_Grid_HYCOM
from flood import flood
from make_remap_grid_file import make_remap_grid_file
from remap import remap
from remap_uv import remap_uv
from shutil import copyfile

src_grd_file = '/media/sf_Swap-between-windows-linux/New_Grid/TEMP_SALT_CURR_2004_2014.nc' #'/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/TSC.nc'
file = '/media/sf_Swap-between-windows-linux/New_Grid/TEMP_SALT_CURR_2004_2014.nc'
dst_dir='./'

# load the grid
srcgrd = get_nc_Grid_HYCOM(src_grd_file)
#dstgrd = pyroms.grid.get_ROMS_grid('FINER')
dstgrd = pyroms.grid.get_ROMS_grid('FINER')

# make remap grid file for scrip
make_remap_grid_file(srcgrd)
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_PUSSY_t.nc'
grid2_file = 'remap_grid_FINER_rho.nc'
interp_file1 = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_FINER_to_PUSSY_bilinear_rho_to_t.nc'
map1_name = 'PUSSY to FINER Bilinear Mapping'
map2_name = 'FINER to PUSSY Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

grid1_file = 'remap_grid_PUSSY_t.nc'
grid2_file = 'remap_grid_FINER_u.nc'
interp_file1 = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_u.nc'
interp_file2 = 'remap_weights_FINER_to_PUSSY_bilinear_u_to_t.nc'
map1_name = 'PUSSY to FINER Bilinear Mapping'
map2_name = 'FINER to PUSSY Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

grid1_file = 'remap_grid_PUSSY_t.nc'
grid2_file = 'remap_grid_FINER_v.nc'
interp_file1 = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_v.nc'
interp_file2 = 'remap_weights_FINER_to_PUSSY_bilinear_v_to_t.nc'
map1_name = 'PUSSY to FINER Bilinear Mapping'
map2_name = 'FINER to PUSSY Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

# load the grid
src_grd = get_nc_Grid_HYCOM(src_grd_file)
#dst_grd = pyroms.grid.get_ROMS_grid('FINER')
dst_grd = pyroms.grid.get_ROMS_grid('FINER')

# Triggering of variables remapping

#zero=3288;l_time=3289
zero=731;l_time=732			#2006-01-01 - 731; 2009-01-01 - 1827

#dst_grd = pyroms.grid.get_ROMS_grid('FINER')
dst_grd = pyroms.grid.get_ROMS_grid('FINER')
remap(zero,l_time, file, 'votemper', src_grd, dst_grd, dst_dir=dst_dir)
remap(zero,l_time, file, 'vosaline', src_grd, dst_grd, dst_dir=dst_dir)
remap_uv(zero,l_time, file,  src_grd, dst_grd, dst_dir=dst_dir)


# For merging of files you have to install nco library: sudo aptitude install nco



zeta = remap(zero,l_time, file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)

# merge file
ic_file = 'Initial_FINER.nc'

out_file = 'sshFINER.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)

out_file = 'votemperFINER.nc' 
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)

out_file = 'vosalineFINER.nc' 
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)

out_file = 'uFINER.nc' 
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)

out_file = 'vFINER.nc' 
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)

#os.remove(out_file)
#print 'file is being copied'
#copyfile('Initial.nc', '/home/eivanov/COAWST/Data/ROMS/Initial/Initial.nc')
