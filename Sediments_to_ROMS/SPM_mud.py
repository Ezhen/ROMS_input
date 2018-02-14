import netCDF4; import numpy as np

rr = netCDF4.Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Coarsest.nc', 'r', format='NETCDF4')
s = rr.variables['Cs_w'][:]
ss=np.zeros((15))
for i in range(15):
	ss[i]=s[i+1]-s[i]
print ss
