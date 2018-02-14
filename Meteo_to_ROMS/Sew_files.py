from nco import Nco
import subprocess
import os
import commands
import _iso

ic_file = 'Forcing.nc'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/windFINER.nc'
out_file = 'windFINER.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)
print 'wind'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/mslFINER.nc'
out_file = 'mslFINER.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)
print 'msl'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/tccFINER.nc'
out_file = 'tccFINER.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)
print 'tcc'


#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/d2mFINER.nc'
out_file = 'd2mFINER.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)
print 'd2m'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/t2mFINER.nc'
out_file = 't2mFINER.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)
print 't2m'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/ssrFINER.nc'
out_file = 'ssrFINER.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)
print 'ssr'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/tpFINER.nc'
out_file = 'tpFINER.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)
print 'tp'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/strFINER.nc'
out_file = 'strFINER.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)
print 'str'
"""
#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/sstFINER.nc'
out_file = 'sstFINER.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
#os.remove(out_file)
print 'sst'
"""
