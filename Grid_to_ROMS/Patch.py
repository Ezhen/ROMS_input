from netCDF4 import Dataset
from shutil import copyfile

copyfile('Child15.nc', 'Child15_2.nc')
n1 = Dataset('Child15_2.nc', 'a', format='NETCDF4')


n1.variables['h'][53:60,0:7] = 9.5


n1.close()
