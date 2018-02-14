import numpy as np
import pyroms
from netCDF4 import Dataset
from mpl_toolkits.basemap import pyproj
from Grid_Nest import Grid_Nest


def get_nc_Grid_Nest(grdfile, name='NESTED'):

    nc = Dataset(grdfile)
    lon = nc.variables['lon_rho'][:]
    lat = nc.variables['lat_rho'][:]
    #depth = nc.variables['depth'][:]
    h = nc.variables['h'][:]
    print "3"
    nc.close()	

    #lon,lat=np.meshgrid(lon,lat)

    lon_t = lon[:,:]
    lat_t = lat[:,:]

    lon_vert = 0.5 * (lon[:,1:] + lon[:,:-1])
    lon_vert = 0.5 * (lon_vert[1:,:] + lon_vert[:-1,:])

    lat_vert = 0.5 * (lat[1:,:] + lat[:-1,:])
    lat_vert = 0.5 * (lat_vert[:,1:] + lat_vert[:,:-1])

    rr = Dataset('Nested_grid.nc', 'a', format='NETCDF4')
    mask_t = rr.variables['mask_rho'][:]

    z_t = 1
    #h=1


    geod = pyproj.Geod(ellps='WGS84')
    az_forward, az_back, dx = geod.inv(lon_vert[:,:-1], lat_vert[:,:-1], lon_vert[:,1:], lat_vert[:,1:])
    angle=rr.variables['angle'][:]
    rr.close()
    #angle = 0.5 * (az_forward[1:,:] + az_forward[:-1,:])
    #angle = (90 - angle) * np.pi/180.


    return Grid_Nest(lon_t, lat_t, lon_vert, lat_vert, mask_t, z_t, h, angle, name)


