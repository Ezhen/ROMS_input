from netCDF4 import Dataset; import numpy as np; from shutil import copyfile

nc = Dataset('Forcing.nc', 'a', format='NETCDF4')

d2m = nc.variables['Qair'][:]
rh=np.copy(d2m)

t2m = nc.variables['Tair'][:]
nc.variables['Tair'][:] = t2m-273.16
print 'Tair is written into the output file'
nc.close()

nc = Dataset('Forcing.nc', 'a', format='NETCDF4')
for i in range(len(d2m[:,0,0])):
	E=6.11*10**(7.5*(d2m[i]/(237.7+t2m[i])))
	Es=6.11*10**(7.5*(t2m[i]/(237.7+t2m[i])))
	nc.variables['Qair'][i]=100*(E/Es)
	print i, 'Qair'
print 'Qair is written into the output file'
nc.close()

nc = Dataset('Forcing.nc', 'a', format='NETCDF4')
Pair = nc.variables['Pair'][:]/100.
nc.variables['Pair'][:] = Pair
print 'Pair is written into the output file'
nc.close()

print 'Everything is done!!!'

#print 'file is being copied'
#copyfile('Forcing.nc', '/home/eivanov/COAWST/Data/ROMS/Forcing/Forcing_2000_2014.nc')
#os.remove(ic_file)
