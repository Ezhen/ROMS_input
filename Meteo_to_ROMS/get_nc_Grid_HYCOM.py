import numpy as np
import pyroms
from netCDF4 import Dataset
from mpl_toolkits.basemap import pyproj
from Grid_HYCOM import Grid_HYCOM


def get_nc_Grid_HYCOM(grdfile, name='PUSSY'):

    nc = Dataset(grdfile)
    lon = nc.variables['longitude'][:]
    lat = nc.variables['latitude'][:]
    #depth = nc.variables['depth'][:]
    print "2"
    var = nc.variables['t2m'][0,:,:]
    nc.close()	

    lon,lat=np.meshgrid(lon,lat)

    lon_t = lon[:,:]
    lat_t = lat[:,:]

    lon_vert = 0.5 * (lon[:,1:] + lon[:,:-1])
    lon_vert = 0.5 * (lon_vert[1:,:] + lon_vert[:-1,:])

    lat_vert = 0.5 * (lat[1:,:] + lat[:-1,:])
    lat_vert = 0.5 * (lat_vert[:,1:] + lat_vert[:,:-1])

    mask_t = np.ones((len(var[:,0]), len(var[0,:])))

    z_t = 1
    h=1


    geod = pyproj.Geod(ellps='WGS84')
    az_forward, az_back, dx = geod.inv(lon_vert[:,:-1], lat_vert[:,:-1], lon_vert[:,1:], lat_vert[:,1:])
    angle = 0.5 * (az_forward[1:,:] + az_forward[:-1,:])
    angle = (90 - angle) * np.pi/180.


    return Grid_HYCOM(lon_t, lat_t, lon_vert, lat_vert, mask_t, z_t, h, angle, name)


