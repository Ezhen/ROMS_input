from nco import Nco
import subprocess
import os
import commands
import _iso

ic_file = 'Forcing_nested.nc'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/windTRANSCHILD.nc'
out_file = 'windTRANSCHILD.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
print 'wind'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/mslTRANSCHILD.nc'
out_file = 'mslTRANSCHILD.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
print 'msl'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/tccTRANSCHILD.nc'
out_file = 'tccTRANSCHILD.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
print 'tcc'


#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/d2mTRANSCHILD.nc'
out_file = 'd2mTRANSCHILD.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
print 'd2m'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/t2mTRANSCHILD.nc'
out_file = 't2mTRANSCHILD.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
print 't2m'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/ssrTRANSCHILD.nc'
out_file = 'ssrTRANSCHILD.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
print 'ssr'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/tpTRANSCHILD.nc'
out_file = 'tpTRANSCHILD.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
print 'tp'

#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/strTRANSCHILD.nc'
out_file = 'strTRANSCHILD.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
print 'str'
"""
#out_file = '/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/sstTRANSCHILD.nc'
out_file = 'sstTRANSCHILD.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
print 'sst'
"""
