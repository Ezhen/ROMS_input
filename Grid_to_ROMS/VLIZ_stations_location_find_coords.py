import numpy as np; import matplotlib.pyplot as plt; from netCDF4 import Dataset; from scipy import spatial

def coord(xx,yy):
	aa=np.array((list(lon.flatten()), list(lat.flatten()))).T
	idxx=spatial.KDTree(aa).query([xx,yy])[1]
	a=int(idxx/width); b=int(idxx-(idxx/width)*width)
	return a,b


y1=51+23/60.+40.04/3600.; x1=3+2./60.+44.82/3600. # MP0 (Measuring Pile 0, Oostdyckbank)
y2=51+23/60.+22.57/3600.; x2=3+11/60.+55.42/3600. # MP3 (Measuring Pile 3, Westhinder)  
y3=51+25/60.+6.08/3600. ; x3=3+17/60.+54.88/3600. # MP4 (Measuring Pile 4, Westhinder)


ncdata1 = Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Coarsest.nc', 'r', format='NETCDF4')
#ncdata1 = Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Finer.nc', 'r', format='NETCDF4')
#ncdata1 = Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Second_Level_Nesting/MeetNet_4.nc', 'r', format='NETCDF4')
lat = ncdata1.variables['lat_rho'][:]; lon = ncdata1.variables['lon_rho'][:]; hs = ncdata1.variables['h'][:];  width=len(hs.T); ncdata1.close()


print coord(x1,y1)
print coord(x2,y2)
print coord(x3,y3)

# Parent grid
#(42, 19) - 8.9
#(40, 18) - 9.1
#(39, 18) - 6.6

# First level of nesting
#(38, 19) - 8.3
#(30, 13) - 8.8
#(23, 12) - 5.8

# I move the boundary lon to 20, 41
# I move the boundary lat to 9, 22


# Second level of nesting
# With new bathymetry
#(92, 55) - 9.15
#(55, 24) - 9.49
#(21, 16) - 10.49

