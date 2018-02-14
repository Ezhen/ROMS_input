from netCDF4 import Dataset
import numpy as np
from shutil import copyfile


nc = Dataset('Forcing_nested.nc', 'a', format='NETCDF4')
lw = nc.variables['lwrad'][:]/43200. 
nc.variables['lwrad'][:] = lw
print 'lwrad is written into the output file'
nc.close()
'''
#copyfile('/media/sf_Swap-between-windows-linux/Meteo_2000_2014/Forcing.nc', 'Forcing.nc')
#['tcc','u10','v10','t2m','d2m';'par','sshf','slhf','ssrd','ssr','str','ewss','nsss','e','ro','ssrc','tp', 'strd']
'''
nc = Dataset('Forcing_nested.nc', 'a', format='NETCDF4')

d2m = nc.variables['Qair'][:]-273.15
rh=d2m

t2m = nc.variables['Tair'][:]-273.15
nc.variables['Tair'][:] = t2m
print 'Tair is written into the output file'
nc.close()

nc = Dataset('Forcing_nested.nc', 'a', format='NETCDF4')
for i in range(len(d2m[:,0,0])):
	rh[i,:,:]=100-5*(t2m[i,:,:]-d2m[i,:,:])
	print i, 'Qair'
nc.variables['Qair'][:] = rh
print 'Qair is written into the output file'
nc.close()



nc = Dataset('Forcing_nested.nc', 'a', format='NETCDF4')
sw = nc.variables['swrad'][:]/43200. 
nc.variables['swrad'][:] = sw
print 'swr is written into the output file'
nc.close()




nc = Dataset('Forcing_nested.nc', 'a', format='NETCDF4')
rain = nc.variables['rain'][:]*1000/43200. 
nc.variables['rain'][:] = rain
print 'rain is written into the output file'
nc.close()





nc = Dataset('Forcing_nested.nc', 'a', format='NETCDF4')
Pair = nc.variables['Pair'][:]/100.
nc.variables['Pair'][:] = Pair
print 'Pair is written into the output file'
nc.close()

print 'Everything is done!!!'

#print 'file is being copied'
#copyfile('Forcing.nc', '/home/eivanov/COAWST/Data/ROMS/Forcing/Forcing_2000_2014.nc')
#os.remove(ic_file)
