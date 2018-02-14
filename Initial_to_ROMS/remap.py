from netCDF4 import Dataset
import numpy as np
from numpy import shape
import pyroms
import pyroms_toolbox
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import _remapping
from get_nc_Grid_HYCOM import get_nc_Grid_HYCOM
from flood import flood

class nctime(object):
    pass

def remap(zero,l_time, src_file, src_varname, src_grd, dst_grd, dmax=20, cdepth=0, kk=0,dst_dir='./'):
    nctime.long_name = 'time'
    nctime.units = 'hours since 2004-01-01 00:00:00'
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    cw=int(len(dst_grd.vgrid.Cs_r))
    cdf = Dataset(src_file)
    if src_varname=='ssh':
	src_var = np.zeros((1,Mp,Lp))
    else:
	src_var = cdf.variables[src_varname][0:1]
    ndim = np.ndim(src_var)
    spval= -32767
    dst_file = src_file.rsplit('/')[-1]
    dst_file =  src_varname + dst_grd.name + '.nc'
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)
    nc = Dataset(dst_file, 'a', format='NETCDF3_64BIT')
    print ndim
    if ndim == 4:
	time = cdf.variables['time'][zero:l_time]
        src_var = src_var[0,:,:,:]
    elif ndim == 3:
	time = cdf.variables['time'][zero:l_time]
	src_var = np.zeros((len(time),Mp,Lp))
    time=time-time[0]
    if src_varname == 'ssh':
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'zeta'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'free-surface'
        units = 'meter'
        field = 'free-surface, scalar, series'
    elif src_varname == 'votemper':
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'temp'
        dimensions = ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        long_name = 'potential temperature'
        units = 'Celsius'
        field = 'temperature, scalar, series'
    elif src_varname == 'vosaline':
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'salt'
        dimensions = ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        long_name = 'salinity'
        units = 'PSU'
        field = 'salinity, scalar, series'

    if ndim == 4:
    	zlevel = -z[::-1,0,0]
	nzlevel = len(zlevel)
	dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
	dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)

    nc.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
    nc.variables[dst_varname].long_name = long_name
    nc.variables[dst_varname].units = units
    nc.variables[dst_varname].field = field
    dxy=5; cdepth=0; kk=0

    
    if ndim == 4:
        dst_var=np.zeros((l_time-zero, cw, Mp, Lp))
	for i in range(len(time)):
		if src_varname == 'votemper':
        		src_varz = flood(cdf.variables[src_varname][zero+i,:,:,:]-273.15, src_grd,  spval=spval)
		else:
        		src_varz = flood(cdf.variables[src_varname][zero+i,:,:,:], src_grd,  spval=spval)
		dst_varz = pyroms.remapping.remap(src_varz, wts_file, spval=spval)
		dst_cont = pyroms.remapping.z2roms(dst_varz[::-1,:,:], dst_grdz, dst_grd, Cpos=Cpos, spval=spval, flood=True)
		dst_var[i,:,:,:]=dst_cont

    else:
        dst_var=np.zeros((l_time-zero, Mp, Lp))
	for i in range(len(time)):
		#src_varz = src_var[i,:,:]
		xx=np.arange(-7,7,1);yy=np.arange(-7,7,1)
		R=np.zeros((len(xx),len(yy)))
		for j in range(len(xx)):
			for k in range(len(yy)):
				dst_var[i,j+17,k+31]=2*np.exp(-0.05*(xx[j]**2+yy[k]**2))

    nc.variables['ocean_time'][:] = time/24.
    print time
    print np.shape(dst_var)
    nc.variables[dst_varname][:] = dst_var
    nc.close()
