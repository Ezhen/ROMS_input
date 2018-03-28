import matplotlib.pyplot as plt; import numpy as np; from netCDF4 import Dataset; from bathy_smoother import *; import pyroms


dstgrd = pyroms.grid.get_ROMS_grid('Parent')

h = dstgrd.vgrid.h.copy()

hraw = h
print 'check bathymetry roughness'
RoughMat = bathy_tools.RoughnessMatrix(dstgrd.vgrid.h, dstgrd.hgrid.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()

print 'smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)'
rx0_max = 0.35
dstgrd.vgrid.h = bathy_smoothing.smoothing_Positive_rx0(dstgrd.hgrid.mask_rho, h, rx0_max)

RoughMat = bathy_tools.RoughnessMatrix(dstgrd.vgrid.h, dstgrd.hgrid.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()

print 'vertical coordinate'
theta_b = 3
theta_s = 7
Tcline=10
N = 20
vgrd = pyroms.vgrid.s_coordinate_4(dstgrd.vgrid.h, theta_b, theta_s, Tcline, N, hraw=hraw)

grd = pyroms.grid.ROMS_Grid('Parent', dstgrd.hgrid, vgrd)

pyroms.grid.write_ROMS_grid(grd, filename='Parent_smoothed.nc')



