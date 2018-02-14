import matplotlib.pyplot as plt; import numpy as np; from netCDF4 import Dataset; from scipy import spatial; from shutil import copyfile; import matplotlib.path as mplPath

xx=[]; yy=[]; zz=[]

def coord(xx,yy):
	idxx=spatial.KDTree(cs).query([xx,yy])[1]
	return hz[idxx]

for line in open('Utm_BCZ.txt','r').readlines():
	xx.append(float(line.split()[0]))
	yy.append(float(line.split()[1]))
	zz.append(float(line.split()[2]))

#xx=np.array(xx); yy=np.array(yy); zz=np.array(zz)
xx=np.array(xx[::5]); yy=np.array(yy[::5]); zz=np.array(zz[::5])

y_bcz=np.array([51.37361, 51.37361, 51.37268, 51.33611, 51.32416, 51.31485, 51.27638, 51.24972, 51.21334, 51.09403, 51.09111, 51.09111, 51.09111, 51.09361, 51.09433, 51.26917, 51.55472, 51.55777, 51.55777, 51.61306, 51.61306, 51.80500, 51.87000, 51.87000, 51.55167, 51.48472, 51.45000, 51.37944, 51.37361, 51.37361])
x_bcz=np.array([3.36472, 3.36472, 3.36491, 3.17972, 3.13166, 3.10403, 3.02000, 2.95528, 2.86305, 2.55555, 2.54166, 2.54166, 2.54166, 2.54361, 2.54298, 2.39028, 2.23973, 2.23812, 2.23812, 2.25333, 2.25333, 2.48167, 2.53944, 2.53944, 3.08139, 3.21222, 3.29639, 3.35389, 3.36472, 3.36472])
coord_bcz=np.zeros((len(y_bcz),2))
coord_bcz[:,0]=y_bcz
coord_bcz[:,1]=x_bcz
bbPath = mplPath.Path(coord_bcz)

#copyfile('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Second_Level_Nesting/Windfarm.nc', 'Windfarm_bath_3.nc')
#ncdata1 = Dataset('Windfarm_bath_3.nc', 'a', format='NETCDF4')
copyfile('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Second_Level_Nesting/MeetNet.nc', 'MeetNet_3.nc')
ncdata1 = Dataset('MeetNet_3.nc', 'a', format='NETCDF4')
x = ncdata1.variables['lat_rho'][:]; y = ncdata1.variables['lon_rho'][:]; hs = ncdata1.variables['h'][:];  width=len(hs.T)
aa=np.array((list(x.flatten()), list(y.flatten()))).T
h_old = np.array((list(hs.flatten())))


#aaa = 51.73; bbb = 3.1 
aaa = x.max(); bbb = y.max()
cs=np.zeros((len(xx),2)); ee,dd,hz = [],[],[]
for i in range(len(xx)):
	if xx[i] < (x.min()-0.01) or xx[i] > (aaa+0.01) or yy[i] < (y.min()-0.01) or yy[i] > (bbb + 0.01):
		pass
	else:
		ee.append(xx[i]); dd.append(yy[i]); hz.append(zz[i])

cs = np.array((ee, dd)).T

for i in range(len(aa)):
	if bbPath.contains_point((aa[i,0],aa[i,1]))==1:
		h_new = coord(aa[i,0],aa[i,1])
		print i, h_old[i], h_new
		o,p = np.where(hs==h_old[i])
		ncdata1.variables['h'][o,p] = h_new*(-1)
	else:
		o,p = np.where(hs==h_old[i])
		ncdata1.variables['h'][o,p] = h_old[i]

ncdata1.close()
