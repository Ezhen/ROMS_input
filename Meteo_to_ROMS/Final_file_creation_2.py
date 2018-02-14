from netCDF4 import Dataset; import numpy as np; from shutil import copyfile

nc = Dataset('Forcing.nc', 'a', format='NETCDF4')
a=np.zeros((len(nc.variables['swrad'][:,0,0]),len(nc.variables['swrad'][0,:,0]),len(nc.variables['swrad'][0,0,:])))
b=np.zeros((len(nc.variables['swrad'][:,0,0]),len(nc.variables['swrad'][0,:,0]),len(nc.variables['swrad'][0,0,:])))
c=np.zeros((len(nc.variables['swrad'][:,0,0]),len(nc.variables['swrad'][0,:,0]),len(nc.variables['swrad'][0,0,:])))
for i in range(len(nc.variables['swrad'][:,0,0])):
	if i%4==0:
		a[i] = nc.variables['swrad'][i]/10800.  
		b[i] = nc.variables['lwrad'][i]/10800. 
		c[i] = nc.variables['rain'][i]*1000/10800. 
	else:
		a[i] = (nc.variables['swrad'][i]-nc.variables['swrad'][i-1])/10800.
		b[i] = (nc.variables['lwrad'][i]-nc.variables['lwrad'][i-1])/10800. 
		c[i] = (nc.variables['rain'][i]-nc.variables['rain'][i-1])*1000/10800.
	print i, 'lwrad before:', nc.variables['lwrad'][i,109,45]/10800., 'lwrad after:', b[i,109,45]
nc.variables['swrad'][:]=a; nc.variables['lwrad'][:]=b; nc.variables['rain'][:]=c
print 'swr is written into the output file','lwrad is written into the output file','rain is written into the output file'
nc.close()


