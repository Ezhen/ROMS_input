from scipy.interpolate import griddata
from netCDF4 import Dataset
import numpy as np
from shutil import copyfile

rr = Dataset('/media/sf_Swap-between-windows-linux/DATA_INPUT_ROMS/Bathymetry/GEBCO.nc', 'r', format='NETCDF4')
lats = rr.variables['lat'][:]
lons = rr.variables['lon'][:]
topo = rr.variables['elevation'][:,:]*(-1)

copyfile('Child_river_boundary.nc', 'Child_river_boundary_update_bath.nc')
nc = Dataset('Child_river_boundary_update_bath.nc', 'a', format='NETCDF4')
lt = nc.variables['lat_rho'][:]
ln = nc.variables['lon_rho'][:]


a1 = np.where(lats<lt.min())[0][-5]
a2 = np.where(lats>lt.max())[0][5]

b1 = np.where(lons<ln.min())[0][-5]
b2 = np.where(lons>ln.max())[0][5]


lats = rr.variables['lat'][a1:a2]
lons = rr.variables['lon'][b1:b2]
topo = rr.variables['elevation'][a1:a2,b1:b2]*(-1)



for k in range(len(lt)):
	print k
	for m in range(len(lt.T)):
		lon_int, lat_int, index_i, index_j, topo_int = [], [], [], [], []
		for i in range(np.where(lats<lt[k,m]-0.025)[0][-1],np.where(lats>lt[k,m]+0.025)[0][0]):		
			for j in range(np.where(lons<ln[k,m]-0.025)[0][-1],np.where(lons>ln[k,m]+0.025)[0][0]):
				lat_int.append(lats[i])
				lon_int.append(lons[j])
				index_i.append(i)
				index_j.append(j)
		topo_int=np.zeros((len(lon_int),len(lat_int)))
		for s in range(len(index_i)):
			for t in range(len(index_j)):
				topo_int[s,t]=topo[index_i[s],index_j[t]]
		lon, lat = np.meshgrid(lon_int, lat_int)
		arr=np.zeros((len(lon_int)*len(lat_int),2))
		arr[:,0]=np.ravel(lon)
		arr[:,1]=np.ravel(lat)
		h = np.ravel(topo_int)
		nc.variables['h'][k,m] = griddata(arr,h, (ln[k,m], lt[k,m]), method='cubic')

nc.close()
rr.close()


