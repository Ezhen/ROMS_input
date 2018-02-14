import matplotlib.pyplot as plt; import numpy as np; from netCDF4 import Dataset; from scipy import spatial; from shutil import copyfile; import utm; import matplotlib.path as mplPath
import pyroms
from bathy_smoother import *

def coord(xx,yy):
	idxx=spatial.KDTree(cs).query([xx,yy])[1]
	return hz[idxx]

y_bcz=np.array([51.37361, 51.37361, 51.37268, 51.33611, 51.32416, 51.31485, 51.27638, 51.24972, 51.21334, 51.09403, 51.09111, 51.09111, 51.09111, 51.09361, 51.09433, 51.26917, 51.55472, 51.55777, 51.55777, 51.61306, 51.61306, 51.80500, 51.87000, 51.87000, 51.55167, 51.48472, 51.45000, 51.37944, 51.37361, 51.37361])
x_bcz=np.array([3.36472, 3.36472, 3.36491, 3.17972, 3.13166, 3.10403, 3.02000, 2.95528, 2.86305, 2.55555, 2.54166, 2.54166, 2.54166, 2.54361, 2.54298, 2.39028, 2.23973, 2.23812, 2.23812, 2.25333, 2.25333, 2.48167, 2.53944, 2.53944, 3.08139, 3.21222, 3.29639, 3.35389, 3.36472, 3.36472])
coord_bcz=np.zeros((len(y_bcz),2))
coord_bcz[:,0]=y_bcz
coord_bcz[:,1]=x_bcz
bbPath = mplPath.Path(coord_bcz)


xx=[]; yy=[]; zz=[]

for line in open('/home/eivanov/coawst_data_prrocessing/New_grid/Real_bathymetry/BCP GM 20X20 VH_COPCO.xyz','r').readlines():
	xx.append(float(line.split()[0]))
	yy.append(float(line.split()[1]))
	zz.append(float(line.split()[2]))

xx=np.array(xx); yy=np.array(yy); zz=np.array(zz)

cs=np.zeros((len(xx),2)); ee,dd,hz = [],[],[]
for i in range(len(xx)):
	e, d = utm.to_latlon(xx[i], yy[i], 31, northern='True')
	ee.append(e)
	dd.append(d)

f = open('Utm_BCZ.txt', 'w')
for i in range(len(xx)): 
	f.write("%s	%s	%s\n" %(ee[i],dd[i],-1*zz[i]))
f.close()

