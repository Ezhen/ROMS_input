import numpy as np
import pyroms
from netCDF4 import Dataset
from mpl_toolkits.basemap import pyproj
from Grid_HYCOM import Grid_HYCOM


def get_nc_Grid_HYCOM(grdfile, name='PUSSY'):

    nc = Dataset(grdfile)
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    depth = nc.variables['depth'][:]
    var = nc.variables['votemper'][0,:,:,:]
    nc.close()	

    lon,lat=np.meshgrid(lon,lat)

    lon_t = lon[:,:]
    lat_t = lat[:,:]

    lon_vert = 0.5 * (lon[:,1:] + lon[:,:-1])
    lon_vert = 0.5 * (lon_vert[1:,:] + lon_vert[:-1,:])

    lat_vert = 0.5 * (lat[1:,:] + lat[:-1,:])
    lat_vert = 0.5 * (lat_vert[:,1:] + lat_vert[:,:-1])

    mask_t = np.array(~var[:].mask, dtype='int')

    z_t = np.tile(depth,(mask_t.shape[2],mask_t.shape[1],1)).T
    print z_t
    depth_bnds = np.zeros(len(depth)+1)
    for i in range(1,len(depth)):
        depth_bnds[i] = 0.5 * (depth[i-1] + depth[i])
    depth_bnds[-1] = 1000

    bottom = pyroms.utility.get_bottom(var[::-1,:,:], mask_t[0], spval=var.fill_value)
    nlev = len(depth)
    bottom = (nlev-1) - bottom
    h = np.zeros(mask_t[0,:].shape)
    for i in range(mask_t[0,:].shape[1]):
        for j in range(mask_t[0,:].shape[0]):
            if mask_t[0,j,i] == 1:
                h[j,i] = depth_bnds[bottom[j,i]+1]


    geod = pyproj.Geod(ellps='WGS84')
    az_forward, az_back, dx = geod.inv(lon_vert[:,:-1], lat_vert[:,:-1], lon_vert[:,1:], lat_vert[:,1:])
    angle = 0.5 * (az_forward[1:,:] + az_forward[:-1,:])
    angle = (90 - angle) * np.pi/180.


    return Grid_HYCOM(lon_t, lat_t, lon_vert, lat_vert, mask_t, z_t, h, angle, name)


