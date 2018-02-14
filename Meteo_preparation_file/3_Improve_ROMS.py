from pylab import *; from netCDF4 import Dataset; from shutil import copyfile; from scipy import signal; import sys

nc = Dataset('/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/Forcing.nc', 'a', format='NETCDF3')
crap = nc.variables['mask_rho'][:,:]; lent=len(nc.variables['tair_time'][:])
#varm=['Uwind','Vwind','Pair','cloud','Qair','Tair','swrad','rain','lwrad']
varm=['tide_Cangle','tide_Cmax','tide_Cmin','tide_Cphase','tide_Eamp','tide_Ephase']


ii=[];jj=[]
for i in range(1,len(crap)-2):
	for j in range(1,len(crap.T)-2):
		if crap[i,j]==0 and sum(crap[i-1:i+2,j-1:j+2])>0:
			ii.append(i);jj.append(j)

for i in range(13):
	mask = 1-crap[ii[i]-1:ii[i]+2,jj[i]-1:jj[i]+2];new_mask=np.tile(mask, (13, 1,1))
	#nc.variables[varm[m]][:,ii[i],jj[i]]=sum(np.tensordot(nc.variables[varm[m]][:,ii[i]-1:ii[i]+2,jj[i]-1:jj[i]+2],kernel,axes=([1,2],[0,1])))/sum(kernel)
	#nc.variables[varm[m]][:,ii[i],jj[i]]=sum(np.einsum('ijk,jk',nc.variables[varm[m]][:,ii[i]-1:ii[i]+2,jj[i]-1:jj[i]+2],kernel))/sum(kernel)
	for m in range(len(varm)):
		#for k in range(lent):
		a=np.nanmean(np.ma.array(np.around(nc.variables[varm[m]][:,ii[i]-1:ii[i]+2,jj[i]-1:jj[i]+2],decimals=1),mask=new_mask),axis=(1,2))
		nc.variables[varm[m]][:,ii[i],jj[i]]=a.filled()
		print i,m
		#	nc.variables[varm[m]][k,ii[i],jj[i]] = np.nanmean(nc.variables[varm[m]][k,ii[i]-1:ii[i]+2,jj[i]-1:jj[i]+2])
		#	sys.stdout.write("Cell %s Var %s Time %s \r" % (str(i),str(m),str(k)) ); sys.stdout.flush()
		#nc.variables[varm[m]][k,ii[i],jj[i]]=sum(np.matrix(nc.variables[varm[m]][k,ii[i]-1:ii[i]+2,jj[i]-1:jj[i]+2])*np.matrix(kernel))/sum(kernel)
	#if i%2==0:
	#	nc.close()
	#	nc = Dataset('/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/Forcing.nc', 'a', format='NETCDF3')
nc.close()
		

