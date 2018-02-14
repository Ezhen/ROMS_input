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

src_grd_file = '/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/TSC.nc'
file = '/media/sf_Swap-between-windows-linux/New_Grid/TEMP_SALT_CURR_2004_2014.nc'
#file = '/media/sf_Swap-between-windows-linux/New_Grid/Climatology_Climatology/Climatology_Climatology.nc'
dst_dir='./'

# load the grid
srcgrd = get_nc_Grid_HYCOM(src_grd_file)
dstgrd = pyroms.grid.get_ROMS_grid('COARSEST')

# make remap grid file for scrip
make_remap_grid_file(srcgrd)
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_PUSSY_t.nc'
grid2_file = 'remap_grid_COARSEST_rho.nc'
interp_file1 = 'remap_weights_PUSSY_to_COARSEST_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_COARSEST_to_PUSSY_bilinear_rho_to_t.nc'
map1_name = 'PUSSY to COARSEST Bilinear Mapping'
map2_name = 'COARSEST to PUSSY Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

grid1_file = 'remap_grid_PUSSY_t.nc'
grid2_file = 'remap_grid_COARSEST_u.nc'
interp_file1 = 'remap_weights_PUSSY_to_COARSEST_bilinear_t_to_u.nc'
interp_file2 = 'remap_weights_COARSEST_to_PUSSY_bilinear_u_to_t.nc'
map1_name = 'PUSSY to COARSEST Bilinear Mapping'
map2_name = 'COARSEST to PUSSY Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

grid1_file = 'remap_grid_PUSSY_t.nc'
grid2_file = 'remap_grid_COARSEST_v.nc'
interp_file1 = 'remap_weights_PUSSY_to_COARSEST_bilinear_t_to_v.nc'
interp_file2 = 'remap_weights_COARSEST_to_PUSSY_bilinear_v_to_t.nc'
map1_name = 'PUSSY to COARSEST Bilinear Mapping'
map2_name = 'COARSEST to PUSSY Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

# load the grid
src_grd = get_nc_Grid_HYCOM(src_grd_file)
dst_grd = pyroms.grid.get_ROMS_grid('COARSEST')


# Triggering of variables remapping

#zero=3288;l_time=3439
zero=731; l_time=1827 #3653 #l_time=30


#remap_bdry(zero,l_time, file, 'votemper', src_grd, dst_grd, dst_dir=dst_dir)
#remap_bdry(zero,l_time, file, 'vosaline', src_grd, dst_grd, dst_dir=dst_dir)
#remap_bdry_uv(zero,l_time, file,  src_grd, dst_grd, dst_dir=dst_dir)
#remap_bdry_uv_zeros(zero,l_time, file,  src_grd, dst_grd, dst_dir=dst_dir)
#remap_bdry(zero,l_time, file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)

#For merging of files you have to install nco library: sudo aptitude install nco
"""
# merge file
ic_file = 'Boundary_COARSEST.nc'

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
"""
#print 'file is being copied'
#copyfile('Boundary.nc', '/home/eivanov/COAWST/Data/ROMS/Boundary/Boundary.nc')
#os.remove(ic_file)

