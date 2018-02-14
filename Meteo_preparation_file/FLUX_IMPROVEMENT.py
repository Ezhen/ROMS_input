from netCDF4 import Dataset; import numpy as np; from shutil import copyfile; var=['sshf','slhf','ssr','str','ewss','nsss','e','tp','tcc','d2m','t2m','msl','u10','v10']; spval=-32767

nc = Dataset('/home/eivanov/coawst_data_prrocessing/ECMWF_Download/Bulk_2006_2009.nc', 'r', format='NETCDF4')


rr = Dataset('Meteo_2006_2009_bulk.nc', 'w', format='NETCDF4')
#rr = Dataset('Qair_correct.nc', 'w', format='NETCDF4')
rr.createDimension('time', len(nc.variables['ssr'][:,0,0]))
rr.createDimension('latitude', len(nc.variables['ssr'][0,:,0]))
rr.createDimension('longitude', len(nc.variables['ssr'][0,0,:]))

rr.createVariable('time', 'f4', ('time')); rr.createVariable('latitude', 'f4', ('latitude')); rr.createVariable('longitude', 'f4', ('longitude'))
rr.variables['time'][:] = nc.variables['time'][:]; rr.variables['latitude'][:] = nc.variables['latitude'][:]; rr.variables['longitude'][:] = nc.variables['longitude'][:]

rr.createVariable('u10', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval); rr.createVariable('v10', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
rr.variables['u10'][:] = nc.variables['u10'][:]; rr.variables['v10'][:] = nc.variables['v10'][:]
"""
################### SENSIBLE HEAT FLUX #########################################
rr.createVariable('sshf', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['sshf'][:,0,0])):
	if i%4==0:
		rr.variables['sshf'][i] = nc.variables['sshf'][i]/10800.
	else:
		rr.variables['sshf'][i] = (nc.variables['sshf'][i]-nc.variables['sshf'][i-1])/10800.
	print i, 'sshf before:', rr.variables['sshf'][i,34,34]/10800., 'sshf after:', rr.variables['sshf'][i,34,34]
print 'sshf is written into the output file'


################### LATENT HEAT FLUX #########################################
rr.createVariable('slhf', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['slhf'][:,0,0])):
	if i%4==0:
		rr.variables['slhf'][i] = nc.variables['slhf'][i]/10800.
	else:
		rr.variables['slhf'][i] = (nc.variables['slhf'][i]-nc.variables['slhf'][i-1])/10800.
	print i, 'slhf before:', nc.variables['slhf'][i,34,34]/10800., 'slhf after:', rr.variables['slhf'][i,34,34]
print 'slhf is written into the output file'
"""

################### SHORTWAVE RADIATION #########################################
rr.createVariable('ssr', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['ssr'][:,0,0])):
	if i%4==0:
		rr.variables['ssr'][i] = nc.variables['ssr'][i]/10800.
	else:
		rr.variables['ssr'][i] = (nc.variables['ssr'][i]-nc.variables['ssr'][i-1])/10800.
	print i, 'ssr before:', nc.variables['ssr'][i,34,34]/10800., 'ssr after:', rr.variables['ssr'][i,34,34]
print 'ssr is written into the output file'


################### LONGWAVE RADIATION #########################################
rr.createVariable('str', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['str'][:,0,0])):
	if i%4==0:
		rr.variables['str'][i] = nc.variables['str'][i]/10800.
	else:
		rr.variables['str'][i] = (nc.variables['str'][i]-nc.variables['str'][i-1])/10800.
	print i, 'str before:', nc.variables['str'][i,34,34]/10800., 'str after:', rr.variables['str'][i,34,34]
print 'str is written into the output file'

"""
################### EAST-WEST SURFACE STRESS #########################################
rr.createVariable('ewss', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['ewss'][:,0,0])):
	if i%4==0:
		rr.variables['ewss'][i] = nc.variables['ewss'][i]/10800.
	else:
		rr.variables['ewss'][i] = (nc.variables['ewss'][i]-nc.variables['ewss'][i-1])/10800.
	print i, 'ewss before:', nc.variables['ewss'][i,34,34]/10800., 'ewss after:', rr.variables['ewss'][i,34,34]
print 'ewss is written into the output file'


################### NORTH-SOUTH SURFACE STRESS #########################################
rr.createVariable('nsss', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['nsss'][:,0,0])):
	if i%4==0:
		rr.variables['nsss'][i] = nc.variables['nsss'][i]/10800.
	else:
		rr.variables['nsss'][i] = (nc.variables['nsss'][i]-nc.variables['nsss'][i-1])/10800.
	print i, 'nsss before:', nc.variables['nsss'][i,34,34]/10800., 'nsss after:', rr.variables['nsss'][i,34,34]
print 'nsss is written into the output file'


################### EVAPORATION #########################################
rr.createVariable('e', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['e'][:,0,0])):
	if i%4==0:
		rr.variables['e'][i] = nc.variables['e'][i]*1000/10800. 
	else:
		rr.variables['e'][i] = (nc.variables['e'][i]-nc.variables['e'][i-1])*1000/10800.
	print i, 'e before:', nc.variables['e'][i,34,34]*1000/10800., 'e after:', rr.variables['e'][i,34,34]
print 'e is written into the output file'

"""
################### PRECIPITATION #########################################
rr.createVariable('tp', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['tp'][:,0,0])):
	if i%4==0:
		rr.variables['tp'][i] = nc.variables['tp'][i]*1000/10800. 
	else:
		rr.variables['tp'][i] = (nc.variables['tp'][i]-nc.variables['tp'][i-1])*1000/10800.
	print i, 'tp before:', nc.variables['tp'][i,34,34]*1000/10800., 'tp after:', rr.variables['tp'][i,34,34]
print 'tp is written into the output file'

################### TOTAL CLOUD COVER #########################################
rr.createVariable('tcc', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['tcc'][:,0,0])):
	#if i%4==0:
	#	rr.variables['tcc'][i] = nc.variables['tcc'][i]
	#else:
	#	rr.variables['tcc'][i] = nc.variables['tcc'][i]-nc.variables['tcc'][i-1]
	rr.variables['tcc'][i] = nc.variables['tcc'][i]
	#print i, 'tcc before:', nc.variables['tcc'][i,34,34], 'tcc after:', rr.variables['tcc'][i,34,34]
print 'tcc is written into the output file'

################### DEW POINT #########################################
rr.createVariable('d2m', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['d2m'][:,0,0])):
	E=6.11*10**(7.5*((nc.variables['d2m'][i]-273.15)/(237.7+(nc.variables['d2m'][i]-273.15))))		####################################
	Es=6.11*10**(7.5*((nc.variables['t2m'][i]-273.15)/(237.7+(nc.variables['t2m'][i]-273.15))))
	rr.variables['d2m'][i] = 100*(E/Es)
	print i, 'd2m before:', nc.variables['d2m'][i,34,34], 'd2m after:', rr.variables['d2m'][i,34,34]
print 'd2m is written into the output file'

################### AIR TEMPERATURE #########################################
rr.createVariable('t2m', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['t2m'][:,0,0])):
	rr.variables['t2m'][i] = nc.variables['t2m'][i]-273.15
	print i, 't2m before:', nc.variables['t2m'][i,34,34], 't2m after:', rr.variables['t2m'][i,34,34]
print 't2m is written into the output file'

################### AIR PRESSURE #########################################
rr.createVariable('msl', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['msl'][:,0,0])):
	rr.variables['msl'][i] = nc.variables['msl'][i]/100.
	print i, 'msl before:', nc.variables['msl'][i,34,34], 'msl after:', rr.variables['msl'][i,34,34]
print 'msl is written into the output file'
rr.close()
"""
rr = Dataset('Meteo_2006_2008_all_improved.nc', 'a', format='NETCDF4')
################### TOTAL HEAT FLUX #########################################
rr.createVariable('shflux', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['msl'][:,0,0])):
	rr.variables['shflux'][i]=rr.variables['ssr'][i]+rr.variables['str'][i]+rr.variables['sshf'][i]+rr.variables['slhf'][i]

################### TOTAL FRESHWATER FLUX #########################################
rr.createVariable('swflux', 'f4', ('time', 'latitude', 'longitude'), fill_value=spval)
for i in range(len(nc.variables['msl'][:,0,0])):
	rr.variables['swflux'][i]=(-rr.variables['e'][i]-rr.variables['tp'][i])*8*10800/10.
rr.close()
"""
nc.close()
