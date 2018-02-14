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
from get_nc_Grid_Nest import get_nc_Grid_Nest
from make_remap_grid_file import make_remap_grid_file
from nc_create_roms_bdry_file_mine import nc_create_roms_bdry_file_mine

#142 137 112 82


def remap_main(tart, tend, src_file, src_varname, src_grd, dst_grd, dmax=20, cdepth=0, kk=0,dst_dir='./'):
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    cdf = Dataset(src_file)
    src_var_dim = cdf.variables[src_varname]
    src_var = cdf.variables[src_varname][0]
    spval= -32767
    dst_file = src_varname + dst_grd.name + '.nc'
    nc_create_roms_bdry_file_mine(dst_file, dst_grd)
    nc = Dataset(dst_file, 'a', format='NETCDF3_64BIT')
    l_time=tend-tart
    time=cdf.variables['time'][tart:tend]
    print time[0], time[-1], shape(time)
    time=time-time[0]
    time = time
    time[0]=0
    #time = (cdf.variables['time'][tart:tend]-920424)/24.
    #time = (cdf.variables['time'][tart:tend]-876579)/24.
    #src_var = src_var[tart:tend,:,:]
    if src_varname == 'Pair' or src_varname == 'sp' or src_varname == 'msl':
        dimensions_time = ('pair_time')
	dst_varname_time = 'pair_time'
	dst_varname_time_long_name = 'pair time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_calendar = 'LEAP'
	dst_varname_time_field = 'pair_time field'
        Bpos = 't'
        Cpos = 'rho'
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'Pair'
        dimensions = ('pair_time', 'er', 'xr')
        long_name = 'surface air pressure'
        units = 'millibar'
        field = 'Pair, scalar, series'
	print np.shape(dst_grd.hgrid.lat_rho)
        nc.createVariable('lon', 'f8', ('er', 'xr'))
        nc.variables['lon'].long_name = 'longitude of RHO-points'
        nc.variables['lon'].units = 'degree_east'
        nc.variables['lon'][:] = dst_grd.hgrid.lon_rho
        nc.createVariable('lat', 'f8', ('er', 'xr'))
        nc.variables['lat'].long_name = 'latitude of RHO-points'
        nc.variables['lat'].units = 'degree_north'
        nc.variables['lat'][:] = dst_grd.hgrid.lat_rho
    elif src_varname == 'cloud' or src_varname == 'tcc':
        dimensions_time = ('cloud_time')
	dst_varname_time = 'cloud_time'
	dst_varname_time_long_name = 'cloud time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_calendar = 'LEAP'
	dst_varname_time_field = 'cloud_time field'
        Bpos = 't'
        Cpos = 'rho'
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'cloud'
        dimensions = ('cloud_time', 'er', 'xr')
        long_name = 'cloud fraction'
        units = 'nondimensional'
        field = 'cloud, scalar, series'
    elif src_varname == 'Uwind' or src_varname == 'u10':
        dimensions_time = ('wind_time')
	dst_varname_time = 'wind_time'
	dst_varname_time_long_name = 'wind time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_calendar = 'LEAP'
	dst_varname_time_field = 'wind_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'Uwind'
        dimensions = ('wind_time', 'er', 'xr')
        long_name = 'surface u-wind component'
        units = 'meter second-1'
        field = 'Uwind, scalar, series'
    elif src_varname == 'Vwind' or src_varname == 'v10':
        dimensions_time = ('wind_time')
	dst_varname_time = 'wind_time'
	dst_varname_time_long_name = 'wind time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_calendar = 'LEAP'
	dst_varname_time_field = 'wind_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'Vwind'
        dimensions = ('wind_time', 'er', 'xr')
        long_name = 'surface v-wind component'
        units = 'meter second-1'
        field = 'Vwind, scalar, series'
    elif src_varname == 'Tair' or src_varname == 't2m':
        dimensions_time = ('tair_time')
	dst_varname_time = 'tair_time'
	dst_varname_time_long_name = 'tair time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_calendar = 'LEAP'
	dst_varname_time_field = 'tair_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'Tair'
        dimensions = ('tair_time', 'er', 'xr')
        long_name = 'surface air temperature'
        units = 'Celcius'
        field = 'Tair, scalar, series'
    elif src_varname == 'Qair' or src_varname == 'd2m':
        dimensions_time = ('qair_time')
	dst_varname_time = 'qair_time'
	dst_varname_time_long_name = 'qair time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_calendar = 'LEAP'
	dst_varname_time_field = 'qair_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'Qair'
        dimensions = ('qair_time', 'er', 'xr')
        long_name = 'surface air relative humidity'
        units = 'percentage'
        field = 'Qair, scalar, series'
    elif src_varname == 'swrad' or src_varname == 'ssr' or src_varname == 'ssrd' or src_varname == 'ssrc':
        dimensions_time = ('srf_time')
	dst_varname_time = 'srf_time'
	dst_varname_time_long_name = 'srf time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_calendar = 'LEAP'
	dst_varname_time_field = 'srf_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'swrad'
        dimensions = ('srf_time', 'er', 'xr')
        long_name = 'solar shortwave radiation flux'
        units = 'Watt/m2'
        field = 'swrad, scalar, series'
    elif src_varname == 'tp':
        dimensions_time = ('rain_time')
	dst_varname_time = 'rain_time'
	dst_varname_time_long_name = 'rain time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_calendar = 'LEAP'
	dst_varname_time_field = 'rain_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'rain'
        dimensions = ('rain_time', 'er', 'xr')
        long_name = 'rain fall rate'
        units = 'kg/ m2 /s'
        field = 'rain, scalar, series'
    elif src_varname == 'par':
	print 'AAAAAAAAAAAAA'
        dimensions_time = ('ocean_time')
	dst_varname_time = 'ocean_time'
	dst_varname_time_long_name = 'ocean time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_field = 'ocean_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'PARout'
        dimensions = ('ocean_time', 'er', 'xr')
        long_name = 'Photosynthetically Available Radiation (PAR)'
        units = 'Watts meter-2'
        field = 'PARout, scalar, series'
    elif src_varname == 'sshf':
        dimensions_time = ('shf_time')
	dst_varname_time = 'shf_time'
	dst_varname_time_long_name = 'shf time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_field = 'shf_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'sensible'
        dimensions = ('shf_time', 'er', 'xr')
        long_name = 'net sensible heat flux'
        units = 'watt meter-2'
        field = 'sensible heat flux, scalar, series'
    elif src_varname == 'slhf':
        dimensions_time = ('lhf_time')
	dst_varname_time = 'lhf_time'
	dst_varname_time_long_name = 'lhf time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_field = 'lhf_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'latent'
        dimensions = ('lhf_time', 'er', 'xr')
        long_name = 'net latent heat flux'
        units = 'watt meter-2'
        field = 'latent heat flux, scalar, series'
    elif src_varname == 'strd':
        dimensions_time = ('lrf_time')
	dst_varname_time = 'lrf_time'
	dst_varname_time_long_name = 'lrf time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_field = 'lrf_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'lwrad_down'
        dimensions = ('lrf_time', 'er', 'xr')
        long_name = 'downwelling longwave radiation flux'
        units = 'watt meter-2'
        field = 'downwelling longwave radiation, scalar, series'
    elif src_varname == 'str':
        dimensions_time = ('lrf_time')
	dst_varname_time = 'lrf_time'
	dst_varname_time_long_name = 'lrf time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_field = 'lrf_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'lwrad'
        dimensions = ('lrf_time', 'er', 'xr')
        long_name = 'net longwave radiation flux'
        units = 'watt meter-2'
        field = 'longwave radiation, scalar, series'
    elif src_varname == 'shflux':
        dimensions_time = ('shf_time')
	dst_varname_time = 'shf_time'
	dst_varname_time_long_name = 'shf time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_field = 'shf_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'shflux'
        dimensions = ('shf_time', 'er', 'xr')
        long_name = 'surface net heat flux'
        units = 'watt meter-2'
        field = 'surface heat flux, scalar, series'
    elif src_varname == 'swflux':
        dimensions_time = ('swf_time')
	dst_varname_time = 'swf_time'
	dst_varname_time_long_name = 'swf time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_field = 'swf time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'swflux'
        dimensions = ('swf_time', 'er', 'xr')
        long_name = 'surface net freswater flux, (E-P)'
        units = 'centimeter day-1'
        field = 'surface net salt flux, scalar, series'
        nc.createVariable('lon', 'f8', ('er', 'xr'))
        nc.variables['lon'].long_name = 'longitude of RHO-points'
        nc.variables['lon'].units = 'degree_east'
        nc.variables['lon'][:] = dst_grd.hgrid.lon_rho
        nc.createVariable('lat', 'f8', ('er', 'xr'))
        nc.variables['lat'].long_name = 'latitude of RHO-points'
        nc.variables['lat'].units = 'degree_north'
        nc.variables['lat'][:] = dst_grd.hgrid.lat_rho
    elif src_varname == 'e':
        dimensions_time = ('evap_time')
	dst_varname_time = 'evap_time'
	dst_varname_time_long_name = 'evap time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_field = 'evap_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'evaporation'
        dimensions = ('evap_time', 'er', 'xr')
        long_name = 'evaporation rate'
        units = 'kilogram meter-2 second-1'
        field = 'evaporation, scalar, series'
    elif src_varname == 'ro':
        dimensions_time = ('runoff_time')
	dst_varname_time = 'runoff_time'
	dst_varname_time_long_name = 'runoff time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_field = 'runoff_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'Runoff'
        dimensions = ('runoff_time', 'er', 'xr')
        long_name = 'Fresh water runoff from land'
        units = 'kg/s/m^2'
        field = 'runoff, scalar, series'
    elif src_varname == 'sst' or 'analysed_sst':
        dimensions_time = ('sst_time')
	dst_varname_time = 'sst_time'
	dst_varname_time_long_name = 'sst time'
	dst_varname_time_units = 'days since 2006-01-01 00:00:00'
	dst_varname_time_calendar = 'LEAP'
	dst_varname_time_field = 'sst_time field'
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'
        dst_varname = 'sst'
        dimensions = ('sst_time', 'er', 'xr')
        long_name = 'surface surface temperature'
        units = 'Kelvin'
        field = 'sst, scalar, series'


    nc.createDimension(dimensions_time, l_time)		###

    nc.createVariable(dst_varname_time, 'f4', dimensions_time)
    nc.variables[dst_varname_time].long_name = dst_varname_time_long_name
    nc.variables[dst_varname_time].units = dst_varname_time_units
    nc.variables[dst_varname_time].field = dst_varname_time_field

    nc.createVariable(dst_varname, 'f4', dimensions, fill_value=spval)
    nc.variables[dst_varname].long_name = long_name
    nc.variables[dst_varname].units = units
    nc.variables[dst_varname].field = field
    nc.variables[dst_varname].coordinates = 'lon lat'
    dxy=5; cdepth=0; kk=0
    
    dst_var=np.zeros((l_time, Mp,Lp))
    for i in range(len(time)):
	src_varz = cdf.variables[src_varname][i,:,:]
	dst_cont = pyroms.remapping.remap(src_varz, wts_file, spval=spval)
	#print np.shape(dst_var), np.shape(dst_cont)
	dst_var[i,:,:]=dst_cont
	print src_varname, i

    #time_2=np.concatenate((time,time+365,time+365+365,time+365+365+365,time+365+365+365+365),axis=0)		###
    time_2=time/24.		

    nc.variables[dst_varname_time][:] = time_2

    #dst_var_2=np.concatenate((dst_var,dst_var,dst_var,dst_var,dst_var),axis=0)		###
    dst_var_2=dst_var
    nc.variables[dst_varname][:] = dst_var_2

    print src_varname, dst_varname, 'is written into the separated netcdf file'
    nc.close()
