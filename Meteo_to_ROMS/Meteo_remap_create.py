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
from remap_main import remap_main
from remap_main_uv import remap_main_uv
from make_remap_grid_file import make_remap_grid_file
#from make_remap_grid_file_FINER import make_remap_grid_file_FINER
from shutil import copyfile



src_grd_file = '/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Meteo_Climatology/Meteo_2006_2009_bulk_cleaned.nc'
file ='/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Meteo_Climatology/Meteo_2006_2009_bulk_cleaned.nc'
dst_dir='./'

# load the grid
srcgrd = get_nc_Grid_HYCOM(src_grd_file)
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
#dst_grd = get_nc_Grid_Nest('FINER_grid.nc')
dst_grd = pyroms.grid.get_ROMS_grid('FINER')

# Triggering of variables remapping 733
tart=0; tend=1096*8
#remap_main(tart, tend, file, 'slhf', src_grd, dst_grd, dst_dir=dst_dir)
#remap_main(tart, tend, file, 'sshf', src_grd, dst_grd, dst_dir=dst_dir)
#remap_main(tart, tend, file, 'str', src_grd, dst_grd, dst_dir=dst_dir)
#remap_main(tart, tend, file, 'ssr', src_grd, dst_grd, dst_dir=dst_dir)
#remap_main(tart, tend, file, 'e', src_grd, dst_grd, dst_dir=dst_dir)
#remap_main(tart, tend, file, 'tp', src_grd, dst_grd, dst_dir=dst_dir)

#remap_main_uv(2,tart, tend, file, src_grd, dst_grd, dst_dir=dst_dir)
#remap_main(tart, tend, file, 'shflux', src_grd, dst_grd, dst_dir=dst_dir)
#remap_main(tart, tend, file, 'swflux', src_grd, dst_grd, dst_dir=dst_dir)

remap_main_uv(1,tart, tend, file, src_grd, dst_grd, dst_dir=dst_dir)
remap_main(tart, tend, file, 'tcc', src_grd, dst_grd, dst_dir=dst_dir)
remap_main(tart, tend, file, 'd2m', src_grd, dst_grd, dst_dir=dst_dir)
remap_main(tart, tend, file, 'msl', src_grd, dst_grd, dst_dir=dst_dir)
remap_main(tart, tend, file, 'ssr', src_grd, dst_grd, dst_dir=dst_dir)
remap_main(tart, tend, file, 'str', src_grd, dst_grd, dst_dir=dst_dir)
remap_main(tart, tend, file, 'tp', src_grd, dst_grd, dst_dir=dst_dir)
remap_main(tart, tend, file, 't2m', src_grd, dst_grd, dst_dir=dst_dir)



