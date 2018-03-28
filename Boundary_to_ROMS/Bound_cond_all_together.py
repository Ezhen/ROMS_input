from netCDF4 import Dataset
import numpy as np
from numpy import shape
import pyroms
import pyroms_toolbox
from datetime import datetime
import subprocess
from shutil import copyfile
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
from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv
from remap_bdry_uv_zeros import remap_bdry_uv_zeros

src_grd_file = '/media/sf_Swap-between-windows-linux/DATA_INPUT_ROMS/Mercator/TEMP_SALT_CURR_2004_2014.nc'
file = '/media/sf_Swap-between-windows-linux/DATA_INPUT_ROMS/Mercator/TEMP_SALT_CURR_2004_2014.nc'

grid_name = 'ParentRiver'

# load the grid
srcgrd = get_nc_Grid_HYCOM(src_grd_file)
dstgrd = pyroms.grid.get_ROMS_grid(grid_name)

# make remap grid file for scrip
make_remap_grid_file(srcgrd)
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_PUSSY_t.nc'
grid2_file = 'remap_grid_%s_rho.nc' %(grid_name)
interp_file1 = 'remap_weights_PUSSY_to_%s_bilinear_t_to_rho.nc' %(grid_name)
interp_file2 = 'remap_weights_%s_to_PUSSY_bilinear_rho_to_t.nc' %(grid_name)
map1_name = 'PUSSY to %s Bilinear Mapping' %(grid_name)
map2_name = '%s to PUSSY Bilinear Mapping' %(grid_name)
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

grid1_file = 'remap_grid_PUSSY_t.nc'
grid2_file = 'remap_grid_%s_u.nc' %(grid_name)
interp_file1 = 'remap_weights_PUSSY_to_%s_bilinear_t_to_u.nc' %(grid_name)
interp_file2 = 'remap_weights_%s_to_PUSSY_bilinear_u_to_t.nc' %(grid_name)
map1_name = 'PUSSY to %s Bilinear Mapping' %(grid_name)
map2_name = '%s to PUSSY Bilinear Mapping' %(grid_name)
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

grid1_file = 'remap_grid_PUSSY_t.nc'
grid2_file = 'remap_grid_%s_v.nc' %(grid_name)
interp_file1 = 'remap_weights_PUSSY_to_%s_bilinear_t_to_v.nc' %(grid_name)
interp_file2 = 'remap_weights_%s_to_PUSSY_bilinear_v_to_t.nc' %(grid_name)
map1_name = 'PUSSY to %s Bilinear Mapping' %(grid_name)
map2_name = '%s to PUSSY Bilinear Mapping' %(grid_name)
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

# load the grid
src_grd = get_nc_Grid_HYCOM(src_grd_file)
dst_grd = pyroms.grid.get_ROMS_grid(grid_name)


# Triggering of variables remapping

#zero=731;l_time=735
zero=731; l_time = 831 #l_time=1827 


remap_bdry(zero,l_time, file, 'votemper', src_grd, dst_grd, grid_name)
remap_bdry(zero,l_time, file, 'vosaline', src_grd, dst_grd, grid_name)
remap_bdry_uv(zero,l_time, file,  src_grd, dst_grd, grid_name)
#####remap_bdry_uv_zeros(zero,l_time, file,  src_grd, dst_grd, grid_name)
remap_bdry(zero,l_time, file, 'ssh', src_grd, dst_grd, grid_name)

#For merging of files you have to install nco library: sudo aptitude install nco

# merge file
ic_file = '/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Boundary_ParentRiver_100_days.nc'

out_file = 'votemper_bdry.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 'vosaline_bdry.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 'u_bdry.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 'v_bdry.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 'ssh_bdry.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)


os.remove('remap_weights_PUSSY_to_%s_bilinear_t_to_rho.nc' %(grid_name))
os.remove('remap_weights_PUSSY_to_%s_bilinear_t_to_u.nc' %(grid_name))
os.remove('remap_weights_PUSSY_to_%s_bilinear_t_to_v.nc' %(grid_name))
os.remove('remap_grid_%s_rho.nc' %(grid_name))
os.remove('remap_grid_%s_u.nc' %(grid_name))
os.remove('remap_grid_%s_v.nc' %(grid_name))
os.remove('remap_grid_PUSSY_t.nc')
os.remove('flood.pyc')
os.remove('get_nc_Grid_HYCOM.pyc')
os.remove('Grid_HYCOM.pyc')
os.remove('make_remap_grid_file.pyc')
os.remove('remap_bdry.pyc')
os.remove('remap_bdry_uv.pyc')
os.remove('remap_bdry_uv_zeros.pyc')
