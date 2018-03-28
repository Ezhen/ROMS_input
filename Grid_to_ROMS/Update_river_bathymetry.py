from netCDF4 import Dataset
from shutil import copyfile
"""
copyfile('Parent_river_boundary.nc', 'Parent_river_boundary_2.nc')
n1 = Dataset('Parent_river_boundary_2.nc', 'a', format='NETCDF4')

def deep_river(a1,a2,b1,b2):
	for i in range(a1,a2):
		for j in range(b1,b2):
			if n1.variables['mask_rho'][i,j]==1:
				n1.variables['h'][i,j] = 6.0

#deep_river(62,70,30,41)
#deep_river(22,31,1,13)
#deep_river(18,19,12,20)
#deep_river(55,60,35,43)
deep_river(27,39,0,16)

n1.close()
"""
copyfile('/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Parent_river_boundary.nc', '/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Parent_river_boundary_6m.nc')

n1 = Dataset('/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Parent_river_boundary_6m.nc', 'a', format='NETCDF4')
msk = n1.variables['mask_rho'][:]

for i in range(len(msk)):
	for j in range(len(msk.T)):
		if msk[i,j]==1 and n1.variables['h'][i,j]<6:
			n1.variables['h'][i,j] = 6.0

n1.close()
