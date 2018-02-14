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
from get_nc_Grid_Nest import get_nc_Grid_Nest
from remap_nested import remap_nested
from remap_nested_uv import remap_nested_uv
from make_remap_grid_file import make_remap_grid_file
#from make_remap_grid_file_CHILD import make_remap_grid_file_CHILD
from shutil import copyfile



src_grd_file = 'Meteo_rude.nc'
file = '/media/sf_Swap-between-windows-linux/New_Grid/Meteo_coastline_shifted/Climatology_Meteo.nc'
#file ='/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Meteo_Climatology/t2m.nc'
dst_dir='./'

# load the grid
srcgrd = get_nc_Grid_HYCOM(src_grd_file)
dstgrd = pyroms.grid.get_ROMS_grid('TRANSCHILD')

# make remap grid file for scrip
make_remap_grid_file(srcgrd)
pyroms.remapping.make_remap_grid_file_2(dstgrd, Cpos='rho')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_PUSSY_t.nc'
grid2_file = 'remap_grid_TRANSCHILD_rho.nc'
interp_file1 = 'remap_weights_PUSSY_to_TRANSCHILD_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_TRANSCHILD_to_PUSSY_bilinear_rho_to_t.nc'
map1_name = 'PUSSY to TRANSCHILD Bilinear Mapping'
map2_name = 'TRANSCHILD to PUSSY Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

# load the grid
src_grd = get_nc_Grid_HYCOM(src_grd_file)
dst_grd = pyroms.grid.get_ROMS_grid('TRANSCHILD')

# Triggering of variables remapping 733
tart=0 ; tend=3652*4 #tart=3288*4; tend=3300*4#; tend=3439*4
remap_nested(tart, tend, file, 'msl', src_grd, dst_grd, dst_dir=dst_dir)
remap_nested(tart, tend, file, 'tcc', src_grd, dst_grd, dst_dir=dst_dir)
remap_nested_uv(1,tart, tend, file, src_grd, dst_grd, dst_dir=dst_dir)
remap_nested(tart, tend, file, 'd2m', src_grd, dst_grd, dst_dir=dst_dir)
remap_nested(tart, tend, file, 'ssr', src_grd, dst_grd, dst_dir=dst_dir)
remap_nested(tart, tend, file, 'str', src_grd, dst_grd, dst_dir=dst_dir)
remap_nested(tart, tend, file, 'tp', src_grd, dst_grd, dst_dir=dst_dir)
remap_nested(tart, tend, file, 't2m', src_grd, dst_grd, dst_dir=dst_dir)



