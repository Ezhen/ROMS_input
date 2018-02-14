from netCDF4 import Dataset; import numpy as np; from shutil import copyfile; spval=-32767

nc = Dataset('/media/sf_Swap-between-windows-linux/FORCING_2006_2009_ATTEMPT_2.nc', 'r', format='NETCDF4')

rr = Dataset('Meteo_2006_2008_all_improved.nc', 'w', format='NETCDF4')
rr.createDimension('time', len(nc.variables['sshf'][:,0,0]))
rr.createDimension('latitude', len(nc.variables['sshf'][0,:,0]))
rr.createDimension('longitude', len(nc.variables['sshf'][0,0,:]))

rr.createVariable('time', 'f4', ('time')); rr.createVariable('latitude', 'f4', ('latitude')); rr.createVariable('longitude', 'f4', ('longitude'))
rr.variables['time'][:] = nc.variables['time'][:]; rr.variables['latitude'][:] = nc.variables['latitude'][:]; rr.variables['longitude'][:] = nc.variables['longitude'][:]

def func_zhopa(name):
	rr.createVariable(name, 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
	for i in range(len(nc.variables[name][:,0,0])):
		if i%4==0:
			rr.variables[name][i] = nc.variables[name][i]/10800. 
		else:
			rr.variables[name][i] = (nc.variables[name][i]-nc.variables[name][i-1])/10800.
		print i, '%s before:' %(name), rr.variables[name][i,34,34]/10800., '%s after:' %(name), rr.variables[name][i,34,34]
	print '%s is written into the output file' %(name)	

func_zhopa('sshf'); func_zhopa('slhf'); func_zhopa('ssr'); func_zhopa('str'); func_zhopa('ewss'); func_zhopa('nsss'); func_zhopa('e'); func_zhopa('tp')
rr.close()

rr = Dataset('Meteo_2006_2008_all_improved.nc', 'a', format='NETCDF4')
################### TOTAL HEAT FLUX #########################################
rr.createVariable('shflux', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['tp'][:,0,0])):
	rr.variables['shflux'][i]=rr.variables['ssr'][i]+rr.variables['str'][i]+rr.variables['sshf'][i]+rr.variables['slhf'][i]

################### TOTAL FRESHWATER FLUX #########################################
rr.createVariable('swflux', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['tp'][:,0,0])):
	rr.variables['swflux'][i]=1000*(-rr.variables['e'][i]-rr.variables['tp'][i])*800
rr.close()
nc.close()
