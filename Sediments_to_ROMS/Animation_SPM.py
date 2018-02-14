import matplotlib.pyplot as plt; import numpy as np; from netCDF4 import Dataset; from datetime import datetime, timedelta; from mpl_toolkits.basemap import Basemap
import matplotlib as mpl; import os, sys; import calendar;  from ezhen.plotbcz import *

lo = lambda x: datetime(2006,1,1,0,0,0) + timedelta(seconds=x)
cmap=mpl.cm.Greys_r
mintemp=0; maxtemp=320; division=10
clevs = np.arange(mintemp,maxtemp,division)
#loc=np.array([0,5,10,17.8,31.6,56.2,100,178,316,340])
loc=np.array([40,50,60,70,80,90,100,110,120,130,140,150,160])
i=np.arange(889,891)


rr = Dataset('SPM_total.nc', 'r', format='NETCDF4')

ncdata1 = Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Coarsest.nc', 'r', format='NETCDF4')
mask_p=ncdata1.variables['mask_rho'][:]; lons_p, lats_p = ncdata1.variables['lon_rho'][:], ncdata1.variables['lat_rho'][:]
ncdata1.close()


def PRINT_PNG(w_p,n):
	t=rr.variables['ocean_time'][n]
	CS1 = m1.contourf(x1,y1,w_p,loc,cmap=cmap)
	file_name = os.path.abspath("SPM_%s" %(str(lo(t).hour))+".png"); fig.savefig(file_name, dpi=200)

fig, ax, cax, m1 = grid_instance(llcrnrlon=1.5, urcrnrlon=4, llcrnrlat=50.5, urcrnrlat=52.5, lat_ts=51.5, r='i', discr=0.5, caxx='vertical')		# build the grid
#fig, ax, cax, m1 = grid_instance(llcrnrlon=-3, urcrnrlon=6.0, llcrnrlat=49, urcrnrlat=55.0, lat_ts=51.5, r='i', discr=0.5, caxx='vertical')
xz,yz = bcz_bound()
x, y = m1(xz, yz)
m1.plot(x,y,color='black',linewidth=1.0)
x1, y1 = m1(lons_p, lats_p)
m1.plot(x1[0],y1[0],color='k',linewidth=1); m1.plot(x1[-1],y1[-1],color='k',linewidth=1); m1.plot(x1[:,0],y1[:,0],color='k',linewidth=1); m1.plot(x1[:,-1],y1[:,-1],color='k',linewidth=1)

#norm = mpl.colors.Normalize(vmin=mintemp, vmax=maxtemp)
bounds=np.arange(mintemp,maxtemp+division,division)
#mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, spacing=loc, boundaries=[-10] + bounds + [10], orientation='vertical')
mpl.colorbar.ColorbarBase(cax, cmap=cmap, boundaries=loc, orientation='vertical')

for k in range(len(i)):
	ww_p = (rr.variables['total'][i[k]])*1000
	w_p = np.ma.masked_where(abs(mask_p-1),ww_p)
	PRINT_PNG(w_p,i[k])

