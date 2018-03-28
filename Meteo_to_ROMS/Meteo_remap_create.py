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



src_grd_file = '/media/sf_Swap-between-windows-linux/DATA_INPUT_ROMS/Meteo/Meteo_2006_2009_bulk_cleaned.nc'
file ='/media/sf_Swap-between-windows-linux/DATA_INPUT_ROMS/Meteo/Meteo_2006_2009_bulk_cleaned.nc'

#grid_name = 'ParentRiver'; grid_flag = 'Parent'
grid_name = 'ChildRiver'; grid_flag = 'Child'

# load the grid
srcgrd = get_nc_Grid_HYCOM(src_grd_file)
dstgrd = pyroms.grid.get_ROMS_grid(grid_name)

# make remap grid file for scrip
make_remap_grid_file(srcgrd)

if grid_flag == 'Parent':
	pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
	pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
	pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')
elif grid_flag == 'Child':
	pyroms.remapping.make_remap_grid_file_2(dstgrd, Cpos='rho')
	pyroms.remapping.make_remap_grid_file_2(dstgrd, Cpos='u')
	pyroms.remapping.make_remap_grid_file_2(dstgrd, Cpos='v')

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

# Triggering of variables remapping 733
tart=0; tend=100*8 #1096*8
#remap_main(tart, tend, file, 'slhf', src_grd, dst_grd)
#remap_main(tart, tend, file, 'sshf', src_grd, dst_grd)
#remap_main(tart, tend, file, 'str', src_grd, dst_grd)
#remap_main(tart, tend, file, 'ssr', src_grd, dst_grd)
#remap_main(tart, tend, file, 'e', src_grd, dst_grd)
#remap_main(tart, tend, file, 'tp', src_grd, dst_grd)

#remap_main_uv(2,tart, tend, file, src_grd, dst_grd)
#remap_main(tart, tend, file, 'shflux', src_grd, dst_grd)
#remap_main(tart, tend, file, 'swflux', src_grd, dst_grd)

remap_main_uv(1,tart, tend, file, src_grd, dst_grd, grid_name)
remap_main(tart, tend, file, 'tcc', src_grd, dst_grd, grid_name)
remap_main(tart, tend, file, 'd2m', src_grd, dst_grd, grid_name)
remap_main(tart, tend, file, 'msl', src_grd, dst_grd, grid_name)
remap_main(tart, tend, file, 'ssr', src_grd, dst_grd, grid_name)
remap_main(tart, tend, file, 'str', src_grd, dst_grd, grid_name)
remap_main(tart, tend, file, 'tp', src_grd, dst_grd, grid_name)
remap_main(tart, tend, file, 't2m', src_grd, dst_grd, grid_name)


ic_file = '/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Meteo_ChildRiver_100_days.nc'
try:
	os.remove(ic_file)
except:
	pass

out_file = 'wind%s.nc' %(grid_name)
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 'tcc%s.nc' %(grid_name)
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 'd2m%s.nc' %(grid_name)
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 'msl%s.nc' %(grid_name)
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 'ssr%s.nc' %(grid_name)
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 'str%s.nc' %(grid_name)
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 'tp%s.nc' %(grid_name)
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)

out_file = 't2m%s.nc' %(grid_name)
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
os.remove('get_nc_Grid_HYCOM.pyc')
os.remove('get_nc_Grid_Nest.pyc')
os.remove('Grid_HYCOM.pyc')
os.remove('Grid_Nest.pyc')
os.remove('make_remap_grid_file.pyc')
os.remove('nc_create_roms_bdry_file_mine.pyc')
os.remove('remap_main.pyc')
os.remove('remap_main_uv.pyc')
