# Reverse negative depths, minimum depth establishement

import matplotlib.pyplot as plt; import numpy as np; from netCDF4 import Dataset; from shutil import copyfile

#copyfile('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Second_Level_Nesting/MeetNet_3.nc', 'MeetNet_4.nc')
#ncdata1 = Dataset('MeetNet_4.nc', 'a', format='NETCDF4')

copyfile('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Second_Level_Nesting/Backup_Windfarm_bath_3.nc', 'Windfarm_bath_4.nc')
ncdata1 = Dataset('Windfarm_bath_4.nc', 'a', format='NETCDF4')

h = ncdata1.variables['h'][:]

for i in range(len(h)):
	for j in range(len(h.T)):
		if h[i,j]<0:
			h[i,j] = -1 * h[i,j]

for i in range(len(h)):
	for j in range(len(h.T)):
		if h[i,j] < 4.0:
			h[i,j] = 4.0

for i in range(1,len(h)-1):
	for j in range(1,len(h.T)-1):
		A = h[i-1:i+2,j-1:j+2]
		if len(np.where(A>h[i,j])[0]) == 0 or len(np.where(A>h[i,j])[0]) == 8:
			h[i,j] = (np.sum(A)-h[i,j])/8

for i in range(len(h)):
	for j in range(len(h.T)):
		ncdata1.variables['h'][i,j] = h[i,j]

ncdata1.close()
