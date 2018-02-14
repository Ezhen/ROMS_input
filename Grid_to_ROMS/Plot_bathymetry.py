import matplotlib.gridspec as gridspec; import matplotlib.pyplot as plt; import matplotlib as mpl; import cmocean
import numpy as np; from netCDF4 import Dataset; from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable

y_bcz=np.array([51.37361, 51.37361, 51.37268, 51.33611, 51.32416, 51.31485, 51.27638, 51.24972, 51.21334, 51.09403, 51.09111, 51.09111, 51.09111, 51.09361, 51.09433, 51.26917, 51.55472, 51.55777, 51.55777, 51.61306, 51.61306, 51.80500, 51.87000, 51.87000, 51.55167, 51.48472, 51.45000, 51.37944, 51.37361, 51.37361]); x_bcz=np.array([3.36472, 3.36472, 3.36491, 3.17972, 3.13166, 3.10403, 3.02000, 2.95528, 2.86305, 2.55555, 2.54166, 2.54166, 2.54166, 2.54361, 2.54298, 2.39028, 2.23973, 2.23812, 2.23812, 2.25333, 2.25333, 2.48167, 2.53944, 2.53944, 3.08139, 3.21222, 3.29639, 3.35389, 3.36472, 3.36472])

vl_name = ['MP 0','MP 3','MP 4']
vl_lat = [51.39445555555555, 51.389602777777775, 51.41835555555555]
vl_lon = [3.045783333333333,3.1987277777777776,3.2985777777777776]

cp = cmocean.cm.deep #mpl.cm.rainbow

ncdata1 = Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Coarsest.nc', 'r', format='NETCDF4')
lats_p = ncdata1.variables['lat_rho'][:]; lons_p = ncdata1.variables['lon_rho'][:]
ncdata2 = Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Finer.nc', 'r', format='NETCDF4')
lats_n = ncdata2.variables['lat_rho'][:]; lons_n = ncdata2.variables['lon_rho'][:]
ncdata3 = Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Second_Level_Nesting/Windfarm_bath_4.nc', 'r', format='NETCDF4')
lats_v = ncdata3.variables['lat_rho'][:]; lons_v = ncdata3.variables['lon_rho'][:]
ncdata4 = Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Second_Level_Nesting/MeetNet_4.nc', 'r', format='NETCDF4')
lats_x = ncdata4.variables['lat_rho'][:]; lons_x = ncdata4.variables['lon_rho'][:]

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
gs1 = gridspec.GridSpec(20, 30)
gs1.update(left=0.05, right=0.95, bottom=0.05, top = 0.95, wspace=0.02)
ax1 = plt.subplot(gs1[0:18,0:18])
ax2 = plt.subplot(gs1[0:10,18:])
ax3 = plt.subplot(gs1[10:,18:])
m1 = Basemap(projection='merc',llcrnrlat=49,urcrnrlat=55.0,llcrnrlon=-3,urcrnrlon=6,lat_ts=51.5,resolution='h', ax=ax1)
m1.drawparallels(np.arange(49,55.0,1),labels=[1,0,0,1],fontsize=10); m1.drawmeridians(np.arange(-3.0,6.0,1),labels=[1,0,0,1],fontsize=10); m1.drawcoastlines(); m1.drawmapboundary(fill_color='aqua')
vmin1 = 4; vmax1 = 88
cax = make_axes_locatable(ax1).append_axes("bottom", size=0.2, pad=0.3); norm = mpl.colors.Normalize(vmin=vmin1, vmax=vmax1); 
mpl.colorbar.ColorbarBase(cax, cmap=cp, norm=norm, orientation='horizontal')
w1=ncdata1.variables['h'][:,:]; w2=ncdata2.variables['h'][:,:]
x1, y1 = m1(lons_p, lats_p); x2, y2 = m1(lons_n, lats_n); x22, y22 = m1(x_bcz, y_bcz); m1.plot(x22,y22,color='black',linewidth=0.5)
clevs1 = np.arange(vmin1,vmax1,0.2); m1.contourf(x1,y1,w1,clevs1,cmap=cp); m1.contourf(x2,y2,w2,clevs1,cmap=cp)
m1.drawcountries(); m1.fillcontinents(color='#ddaa66',lake_color='#9999FF')
m1.plot(x1[0],y1[0],color='k',linewidth=0.5); m1.plot(x1[-1],y1[-1],color='k',linewidth=0.5); m1.plot(x1[:,0],y1[:,0],color='k',linewidth=0.5); m1.plot(x1[:,-1],y1[:,-1],color='k',linewidth=0.5)
m1.plot(x2[0],y2[0],color='k',linewidth=0.5); m1.plot(x2[-1],y2[-1],color='k',linewidth=0.5); m1.plot(x2[:,0],y2[:,0],color='k',linewidth=0.5); m1.plot(x2[:,-1],y2[:,-1],color='k',linewidth=0.5)


m2 = Basemap(projection='merc',llcrnrlat=51,urcrnrlat=52,llcrnrlon=1.8,urcrnrlon=3.7,lat_ts=51.4,resolution='h', ax=ax2)
m2.drawparallels(np.arange(51,52,0.2),labels=[0,0,0,0],fontsize=10); m2.drawmeridians(np.arange(1.8,3.7,0.25),labels=[0,0,0,0],fontsize=10); m2.drawcoastlines(); m2.drawmapboundary(fill_color='aqua')
vmin2 = 4; vmax2 = 52
cax2 = make_axes_locatable(ax2).append_axes("right", size=0.18, pad=0.18); norm2 = mpl.colors.Normalize(vmin=vmin2, vmax=vmax2)
mpl.colorbar.ColorbarBase(cax2, cmap=cp, norm=norm2, orientation='vertical')
x25, y25 = m2(lons_p, lats_p); x3, y3 = m2(lons_n, lats_n); x4, y4 = m2(lons_v, lats_v); x41, y41 = m2(lons_x, lats_x); x42, y42 = m2(x_bcz, y_bcz); m2.plot(x42,y42,color='black',linewidth=0.5)
w3=ncdata3.variables['h'][:,:]
clevs2 = np.arange(vmin2,vmax2,0.2); m2.contourf(x25,y25,w1,clevs2,cmap=cp); m2.contourf(x3,y3,w2,clevs2,cmap=cp); m2.contourf(x4,y4,w3,clevs2,cmap=cp)
m2.drawcountries(); m2.fillcontinents(color='#ddaa66',lake_color='#9999FF')
m2.plot(x4[0],y4[0],color='k',linewidth=0.5); m2.plot(x4[-1],y4[-1],color='k',linewidth=0.5); m2.plot(x4[:,0],y4[:,0],color='k',linewidth=0.5); m2.plot(x4[:,-1],y4[:,-1],color='k',linewidth=0.5)

w4=ncdata4.variables['h'][:,:]
m2.contourf(x41,y41,w4,clevs2,cmap=cp)
m2.drawcountries(); m2.fillcontinents(color='#ddaa66',lake_color='#9999FF')
m2.plot(x41[0],y41[0],color='k',linewidth=0.5); m2.plot(x41[-1],y41[-1],color='k',linewidth=0.5); m2.plot(x41[:,0],y41[:,0],color='k',linewidth=0.5); m2.plot(x41[:,-1],y41[:,-1],color='k',linewidth=0.5)

x43,y43=m2(vl_lon,vl_lat)
m2.scatter(x43[0],y43[0],marker='o',c='r',s=40,label='MP0',edgecolor='w')
m2.scatter(x43[1],y43[1],marker='o',c='y',s=40,label='MP3',edgecolor='w')
m2.scatter(x43[2],y43[2],marker='o',c='g',s=40,label='MP4',edgecolor='w')
aaaa = ax2.legend(ncol=1,title='MeetNet',prop={'family':'Times New Roman','size':10},scatterpoints =1,labelspacing=1, loc=2,fontsize = 'large')
plt.setp(aaaa.get_title(), fontsize=12, family='Times New Roman')



m3 = Basemap(projection='merc',llcrnrlat=51.45,urcrnrlat=51.8,llcrnrlon=2.6,urcrnrlon=3.2,lat_ts=51.65,resolution='i', ax=ax3)
m3.drawparallels(np.arange(51.5,51.85,0.1),labels=[0,0,0,0],fontsize=10); m3.drawmeridians(np.arange(2.6,3.2,0.1),labels=[0,0,0,0],fontsize=10); m2.drawcoastlines(); m3.drawmapboundary(fill_color='aqua')
vmin3 = 4; vmax3 = 40
x5, y5 = m3(lons_n, lats_n); x6, y6 = m3(lons_v, lats_v)
cax3 = make_axes_locatable(ax3).append_axes("right", size=0.2, pad=0.2); norm3 = mpl.colors.Normalize(vmin=vmin3, vmax=vmax3); 
mpl.colorbar.ColorbarBase(cax3, cmap=cp, norm=norm3, orientation='vertical')
clevs3 = np.arange(vmin3,vmax3,0.2); m3.contourf(x5,y5,w2,clevs3,cmap=cp); m3.contourf(x6,y6,w3,clevs3,cmap=cp)
m3.drawcountries(); m3.fillcontinents(color='#ddaa66',lake_color='#9999FF')
m3.plot(x6[0],y6[0],color='k',linewidth=0.5); m3.plot(x6[-1],y6[-1],color='k',linewidth=0.5); m3.plot(x6[:,0],y6[:,0],color='k',linewidth=0.5); m3.plot(x6[:,-1],y6[:,-1],color='k',linewidth=0.5)


y = []; x = []; f = []
for line in open('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Second_Level_Nesting/All_operational_wind_farms_together.txt','r').readlines():
	x.append(float(line.split( )[1])); y.append(float(line.split( )[2])); f.append(float(line.split( )[3]))
x = np.array(x); y = np.array(y); f = np.array(f)


for i in range(len(f)):
	xx,yy=m3(x[i],y[i])
	if f[i]==1:
		m3.plot(xx,yy,marker='o',color='y',markersize=5)#,label='1')
	if f[i]==2:
		m3.plot(xx,yy,marker='o',color='k',markersize=5)#,label='2')
	if f[i]==3:
		m3.plot(xx,yy,marker='o',color='m',markersize=5)#,label='3')
	if f[i]==4:
		m3.plot(xx,yy,marker='o',color='w',markersize=5)#,label='4')

m3.scatter(0,0,marker='o',c='m',s=40,label='Nobelwind',edgecolor='w')
m3.scatter(0,0,marker='o',c='k',s=40,label='C-Power',edgecolor='w')
m3.scatter(0,0,marker='o',c='y',s=40,label='Belwind',edgecolor='w')
m3.scatter(0,0,marker='o',c='w',s=40,label='Northwind',edgecolor='k')
aaa = ax3.legend(title='Wind farms',prop={'family':'Times New Roman','size':10}, scatterpoints =1, labelspacing=1, loc=1,fontsize = 'large')
plt.setp(aaa.get_title(), fontsize=12, family='Times New Roman')

fig.savefig('Grids_bathymetry_200.png', dpi=200)
fig.savefig('Grids_bathymetry_300.png', dpi=300)
plt.show()
