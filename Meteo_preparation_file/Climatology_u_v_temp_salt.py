import numpy as np; import matplotlib as mpl; from netCDF4 import Dataset; import datetime

d2 = lambda x: datetime.datetime(1985,1,1,0,0,0) + datetime.timedelta(days=x)

rr = Dataset('/media/sf_Swap-between-windows-linux/New_Grid/TEMP_SALT_CURR_2004_2014.nc', 'r', format='NETCDF3')
varm=['votemper','vosaline','vomecrty','vozocrtx']
aa=9; length=3655

def lol(v,j,p):
	m=rr.variables[v][j-7:j+7,p,:,:]
	n=m[0,:,:].copy()
	for r in range(len(m)-1):
		n=np.matrix(m[r])+np.matrix(n)
	n=n/14.
	return n
	
for pp in range(len(varm)):
	var=varm[pp]
	varC=np.zeros((length, aa, 100, 102))
	a=rr.variables[var][0:365].mask.copy()
	print np.shape(a)
	for j in range(7,length-7-1):
		for p in range(0,aa):
			varC[j,p]=lol(var,j,p)
		if j%1000==0:
			print var,j
	for i in range(0,8):
		for p in range(0,aa):
			varC[i,p]=rr.variables[var][i,p]
			varC[length-1-i,p]=rr.variables[var][length-1-i,p]
	hf=varC.copy()							# 0, 14616, 29233
	#for i in range(len(hf)):
	#	for p in range(0,aa):
	#		hf[i,p]=np.matrix(rr.variables[var][i,p])-np.matrix(varC[i,p])
	lfc=np.zeros((365, aa, 100, 102)); clim=np.zeros((365, aa, 100, 102)); kk=0
	for p in range(len(lfc[0,:,0,0])):
		for i in range(365*1+1,365*9+3):#		for i in range(length-2):
			print rr.variables['time'][i]
			if d2(int(rr.variables['time'][i])).month==2 and d2(int(rr.variables['time'][i])).day==29:
				pass
			else:
				lfc[kk,p]=np.matrix(lfc[kk,p])+np.matrix(varC[i,p]) #lfc[kk,p]=lfc[kk,p]+hf[i,p]; kk=kk+1
				if kk==364:
					kk=0
				else:
					kk=kk+1
				print i,p,kk
	lfc=lfc/8.	#lfc=lfc/10.						# year 2007-01-01 - 8035 - 1096; year 2008-01-01 - 8400 - 1461
	for m in range(len(lfc)):
		for p in range(0,aa):
			clim[m,p]=np.matrix(lfc[m,p])#+np.matrix(varC[1096+m,p,:,:])
	nc = Dataset('/media/sf_Swap-between-windows-linux/New_Grid/'+var+'.nc', 'w', format='NETCDF3_64BIT')
	nc.createDimension('time', 365); nc.createDimension('depth',aa); nc.createDimension('latitude', 100); nc.createDimension('longitude', 102)
	nc.createVariable('time', 'f4', ('time')); nc.createVariable('depth', 'f4', ('depth')); nc.createVariable('latitude', 'f4', ('latitude')); nc.createVariable('longitude', 'f4', ('longitude'))
	nc.createVariable(var, 'f4', ('time', 'depth', 'latitude', 'longitude'))
	nc.variables['time'][:]=np.arange(0,365); nc.variables['depth'][:]=rr.variables['depth'][0:aa]; nc.variables['latitude'][:]=rr.variables['lat'][:]; nc.variables['longitude'][:]=rr.variables['lon'][:]
	clim=np.ma.masked_array(clim,a)
	print type(clim)
	nc.variables[var][:]=clim[:]
	nc.close()
