from netCDF4 import Dataset
from pylab import *

rr=Dataset('u_bdry.nc', 'r', format='netCDF4')


ue=rr.variables['u_east'][:]
rr.close()
print shape(ue)

for i in range(365):
    for j in range(15):
        for k in range(112):
            if ue[i,j,k]<-0.1 or ue[i,j,k]>0.1:
                print i,j,k



