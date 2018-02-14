import numpy as np; import matplotlib.pyplot as plt; from netCDF4 import Dataset; from scipy import spatial

def coord(xx,yy):
	aa=np.array((list(lon.flatten()), list(lat.flatten()))).T
	idxx=spatial.KDTree(aa).query([xx,yy])[1]
	a=int(idxx/width); b=int(idxx-(idxx/width)*width)
	return a,b


y = []; x = []; f = []

for line in open('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Second_Level_Nesting/All_operational_wind_farms_together.txt','r').readlines():
	x.append(float(line.split( )[1]))
	y.append(float(line.split( )[2]))
	f.append(float(line.split( )[3]))

x = np.array(x); y = np.array(y); f = np.array(f)

ncdata1 = Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Finer.nc', 'r', format='NETCDF4')
lat = ncdata1.variables['lat_rho'][:]; lon = ncdata1.variables['lon_rho'][:]; hs = ncdata1.variables['h'][:];  width=len(hs.T)

a,b = [],[]
for i in range(len(x)):
	print i
	aa,bb = coord(x[i],y[i])
	a.append(aa); b.append(bb)

# min(a) = 28; max(a) = 37
# min(b) = 34; max(b) = 55
# I move the boundary lon to 25, 40
# I move the boundary lat to 31, 58
