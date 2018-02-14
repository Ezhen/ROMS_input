from netCDF4 import Dataset; import numpy as np; from shutil import copyfile; import sys,os,shutil

rr = Dataset('/media/sf_Swap-between-windows-linux/METEO_2006_2008.nc', 'r', format='NETCDF3'); msk = rr.variables['sst'][1,:,:]; rr.close()

copyfile('Meteo_2006_2009_bulk.nc', 'Meteo_2006_2009_bulk_cleaned.nc')
nc = Dataset('Meteo_2006_2009_bulk_cleaned.nc', 'a', format='NETCDF4')
#varm=['sshf','slhf','ssr','str','ewss','nsss','e','tp','tcc','d2m','t2m','msl','u10','v10','shflux','swflux']; lmt=276.17
varm=['ssr','str','tp','tcc','d2m','t2m','msl','u10','v10']; lmt=276.17
nc.variables['time'] = nc.variables['time'][:]
for m in range(len(varm)):
	for k in range(len(nc.variables['time'][:])):
		var = nc.variables[varm[m]][k]
		varC=var.copy()
		for i in range(1,len(msk[:,0])-2):
			for j in range(1,len(msk[0,:])-2):
				if msk[i,j]<lmt:
					a=[]
					if msk[i-1,j-1]>lmt:
						a.append(var[i-1,j-1])
					if msk[i-1,j]>lmt:
						a.append(var[i-1,j])
					if msk[i-1,j+1]>lmt:
						a.append(var[i-1,j+1])
					if msk[i,j-1]>lmt:
						a.append(var[i,j-1])
					if msk[i,j+1]>lmt:
						a.append(var[i,j+1])
					if msk[i+1,j-1]>lmt:
						a.append(var[i+1,j-1])
					if msk[i+1,j]>lmt:
						a.append(var[i+1,j])
					if msk[i+1,j+1]>lmt:
						a.append(var[i+1,j+1])
					if msk[i-2,j-2]>lmt:
						a.append(var[i-2,j-2])
					if msk[i-1,j-2]>lmt:
						a.append(var[i-1,j-2])
					if msk[i,j-2]>lmt:
						a.append(var[i,j-2])
					if msk[i+1,j-2]>lmt:
						a.append(var[i+1,j-2])
					if msk[i+2,j-2]>lmt:
						a.append(var[i+2,j-2])
					if msk[i-2,j-1]>lmt:
						a.append(var[i-2,j-1])
					if msk[i+2,j-1]>lmt:
						a.append(var[i+2,j-1])
					if msk[i-2,j]>lmt:
						a.append(var[i-2,j])
					if msk[i+2,j]>lmt:
						a.append(var[i+2,j])
					if msk[i-2,j+1]>lmt:
						a.append(var[i-2,j+1])
					if msk[i+2,j+1]>lmt:
						a.append(var[i+2,j+1])
					if msk[i-2,j+2]>lmt:
						a.append(var[i-2,j+2])
					if msk[i-1,j+2]>lmt:
						a.append(var[i-1,j+2])
					if msk[i,j+2]>lmt:
						a.append(var[i,j+2])
					if msk[i+1,j+2]>lmt:
						a.append(var[i+1,j+2])
					if msk[i+2,j+2]>lmt:
						a.append(var[i+2,j+2])
					if len(a)>0:
						varC[i,j]=sum(a)/len(a)
					else:
						pass
		sys.stdout.write("%s is being treated at the time step %s\r" % (varm[m],k) ); sys.stdout.flush()
		nc.variables[varm[m]][k] = varC[:,:]
nc.close()
