import numpy as np
import os
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import matplotlib.pyplot as plt
import time
from datetime import datetime
from matplotlib.dates import date2num, num2date
from flood import flood

import pyroms
import pyroms_toolbox
import _remapping



def remap_bdry(zero,l_time,src_file, src_varname, src_grd, dst_grd, grid_name):

    print src_file
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    cw=int(len(dst_grd.vgrid.Cs_r))
    # create boundary file
    dst_file =  src_varname + '_bdry' + '.nc'
    print '\nCreating boundary file', dst_file
    if os.path.exists(dst_file) is True:
        os.remove(dst_file)
    pyroms_toolbox.nc_create_roms_bdry_file(dst_file, dst_grd)

    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF4')
    cdf = netCDF.Dataset(src_file) 
    spval= -32767
    time = cdf.variables['time'][zero:l_time]
    time=time-time[0]
    if src_varname == 'ssh':
	src_var=np.zeros((len(time),Mp, Lp)); ndim=3
    else:
	ndim=4

    if src_varname == 'ssh':
        pos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_%s_bilinear_t_to_rho.nc' %(grid_name)
        dst_varname = 'zeta'
        dimensions_time = ('zeta_time')
	dst_varname_time = 'zeta_time'
	long_name = "free-surface time"
	unitstime = "days since 2006-01-01 00:00:00 GMT" 
	calendar = "gregorian" 
        long_name = 'free-surface'
        dst_varname_north = 'zeta_north'
        dimensions_north = ('zeta_time', 'xi_rho')
        long_name_north = 'free-surface north boundary condition'
        field_north = 'zeta_north, scalar, series'
        dst_varname_south = 'zeta_south'
        dimensions_south = ('zeta_time', 'xi_rho')
        long_name_south = 'free-surface south boundary condition'
        field_south = 'zeta_south, scalar, series'
        dst_varname_east = 'zeta_east'
        dimensions_east = ('zeta_time', 'eta_rho')
        long_name_east = 'free-surface east boundary condition'
        field_east = 'zeta_east, scalar, series'
        dst_varname_west = 'zeta_west'
        dimensions_west = ('zeta_time', 'eta_rho')
        long_name_west = 'free-surface west boundary condition'
        field_west = 'zeta_west, scalar, series'
        units = 'meter'
    elif src_varname == 'temperature' or src_varname == 'votemper':
        pos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_%s_bilinear_t_to_rho.nc' %(grid_name)
        dst_varname = 'temperature'
	dimensions_time = ('temp_time')
	dst_varname_time = 'temp_time'
	long_name = "potential temperature time"
	unitstime = "days since 2004-01-01 00:00:00 GMT"
	calendar = "gregorian"
        dst_varname_north = 'temp_north'
        dimensions_north = ('temp_time', 's_rho', 'xi_rho')
        long_name_north = 'potential temperature north boundary condition'
        field_north = 'temp_north, scalar, series'
        dst_varname_south = 'temp_south'
        dimensions_south = ('temp_time', 's_rho', 'xi_rho')
        long_name_south = 'potential temperature south boundary condition'
        field_south = 'temp_south, scalar, series'
        dst_varname_east = 'temp_east'
        dimensions_east = ('temp_time', 's_rho', 'eta_rho')
        long_name_east = 'potential temperature east boundary condition'
        field_east = 'temp_east, scalar, series'
        dst_varname_west = 'temp_west'
        dimensions_west = ('temp_time', 's_rho', 'eta_rho')
        long_name_west = 'potential temperature west boundary condition'
        field_west = 'temp_west, scalar, series'
        units = 'Celsius'
    elif src_varname == 'salinity' or src_varname == 'vosaline':
        pos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_PUSSY_to_%s_bilinear_t_to_rho.nc' %(grid_name)
        dst_varname = 'salinity'
        dimensions_time = ('salt_time')
	dst_varname_time = 'salt_time'
	long_name = "surface net heat flux time" 
	unitstime = "days since 2004-01-01 00:00:00 GMT" 
	calendar = "gregorian" 
        dst_varname_north = 'salt_north'
        dimensions_north = ('salt_time', 's_rho', 'xi_rho')
        long_name_north = 'salinity north boundary condition'
        field_north = 'salt_north, scalar, series'
        dst_varname_south = 'salt_south'
        dimensions_south = ('salt_time', 's_rho', 'xi_rho')
        long_name_south = 'salinity south boundary condition'
        field_south = 'salt_south, scalar, series'
        dst_varname_east = 'salt_east'
        dimensions_east = ('salt_time', 's_rho', 'eta_rho')
        long_name_east = 'salinity east boundary condition'
        field_east = 'salt_east, scalar, series'
        dst_varname_west = 'salt_west'
        dimensions_west = ('salt_time', 's_rho', 'eta_rho')
        long_name_west = 'salinity west boundary condition'
        field_west = 'salt_west, scalar, series'
        units = 'PSU'
    else:
        raise ValueError, 'Undefined src_varname'


    if ndim == 4:
        # build intermediate zgrid
        zlevel = -z[::-1,0,0]
        nzlevel = len(zlevel)
        dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
        dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)


    # create variable in boudary file
    print 'Creating dimension'
    nc.createDimension(dst_varname_time, l_time-zero)			###
    print l_time-zero, len(time)

    print 'Creating variable', dst_varname_north
    nc.createVariable(dst_varname_north, 'f8', dimensions_north, fill_value=spval)
    nc.variables[dst_varname_north].long_name = long_name_north
    nc.variables[dst_varname_north].units = units
    nc.variables[dst_varname_north].field = field_north
    #nc.variables[dst_varname_north]._FillValue = spval

    print 'Creating variable', dst_varname_south
    nc.createVariable(dst_varname_south, 'f8', dimensions_south, fill_value=spval)
    nc.variables[dst_varname_south].long_name = long_name_south
    nc.variables[dst_varname_south].units = units
    nc.variables[dst_varname_south].field = field_south
    #nc.variables[dst_varname_south]._FillValue = spval

    print 'Creating variable', dst_varname_east
    nc.createVariable(dst_varname_east, 'f8', dimensions_east, fill_value=spval)
    nc.variables[dst_varname_east].long_name = long_name_east
    nc.variables[dst_varname_east].units = units
    nc.variables[dst_varname_east].field = field_east
    #nc.variables[dst_varname_east]._FillValue = spval

    print 'Creating variable', dst_varname_west
    nc.createVariable(dst_varname_west, 'f8', dimensions_west, fill_value=spval)
    nc.variables[dst_varname_west].long_name = long_name_west
    nc.variables[dst_varname_west].units = units
    nc.variables[dst_varname_west].field = field_west
    #nc.variables[dst_varname_west]._FillValue = spval

    print 'Creating variable', dst_varname_time
    nc.createVariable(dst_varname_time, 'f4', dimensions_time)
    nc.variables[dst_varname_time].long_name = long_name
    nc.variables[dst_varname_time].units = unitstime
    nc.variables[dst_varname_time].calendar = calendar


    
    if ndim == 4:
        dst_var_north=np.zeros((l_time-zero, cw, 1, Lp))
        dst_var_south=np.zeros((l_time-zero, cw, 1, Lp))
        dst_var_east=np.zeros((l_time-zero, cw, Mp, 1))
        dst_var_west=np.zeros((l_time-zero, cw, Mp, 1))
	for i in range(len(time)):
       		if src_varname == 'votemper':
        		src_varz = flood(cdf.variables[src_varname][zero+i,:,:,:]-273.15, src_grd,  spval=spval)
       		else:
        		src_varz = flood(cdf.variables[src_varname][i,:,:,:], src_grd,  spval=spval)
		dst_varz = pyroms.remapping.remap(src_varz, wts_file, spval=spval)
        	dst_var_north_cont = pyroms.remapping.z2roms(dst_varz[::-1, Mp-1:Mp, :], dst_grdz, dst_grd, Cpos=Cpos, spval=spval, flood=True, irange=(0,Lp), jrange=(Mp-1,Mp))
        	dst_var_south_cont = pyroms.remapping.z2roms(dst_varz[::-1, 0:1, :], dst_grdz, dst_grd, Cpos=Cpos, spval=spval, flood=True, irange=(0,Lp), jrange=(0,1))
        	dst_var_east_cont = pyroms.remapping.z2roms(dst_varz[::-1, :, Lp-1:Lp], dst_grdz, dst_grd, Cpos=Cpos, spval=spval, flood=True, irange=(Lp-1,Lp), jrange=(0,Mp))
        	dst_var_west_cont = pyroms.remapping.z2roms(dst_varz[::-1, :, 0:1], dst_grdz, dst_grd, Cpos=Cpos, spval=spval, flood=True, irange=(0,1), jrange=(0,Mp))
		dst_var_north[i,:,:,:]=dst_var_north_cont
		dst_var_south[i,:,:,:]=dst_var_south_cont
		dst_var_east[i,:,:,:]=dst_var_east_cont
		dst_var_west[i,:,:,:]=dst_var_west_cont
		print i
    else:
	dst_var_north=np.zeros((l_time-zero, Lp))
        dst_var_south=np.zeros((l_time-zero, Lp))
        dst_var_east=np.zeros((l_time-zero, Mp))
        dst_var_west=np.zeros((l_time-zero, Mp))
	for i in range(len(time)):
		#src_varz = src_var[i,:,:]
		src_varz = np.zeros((len(src_var[0,:,0]),len(src_var[0,0,:])))
		dst_varz = pyroms.remapping.remap(src_varz, wts_file, spval=spval)
		dst_var_north[i,:] = np.zeros((len(dst_varz[-1, :])))
        	dst_var_south[i,:] = np.zeros((len(dst_varz[0, :])))
        	dst_var_east[i,:] = np.zeros((len(dst_varz[:, -1])))
        	dst_var_west[i,:] = np.zeros((len(dst_varz[:, 0])))
		print dst_var_south[i,:]

    # write data in destination file
    print 'write data in destination file'
    #time_2=np.concatenate((time,time+365,time+365+365,time+365+365+365,time+365+365+365+365),axis=0)		####
    time_2=time					####
    nc.variables[dst_varname_time][:] = time_2
    print  time


    #dst_var_north_2=np.concatenate((dst_var_north,dst_var_north,dst_var_north,dst_var_north,dst_var_north),axis=0) ####	
    #dst_var_south_2=np.concatenate((dst_var_south,dst_var_south,dst_var_south,dst_var_south,dst_var_south),axis=0) ####	
    #dst_var_east_2=np.concatenate((dst_var_east,dst_var_east,dst_var_east,dst_var_east,dst_var_east),axis=0)	 ####
    #dst_var_west_2=np.concatenate((dst_var_west,dst_var_west,dst_var_west,dst_var_west,dst_var_west),axis=0)	 ####
    #nc.variables[dst_varname_north][:] = np.squeeze(dst_var_north_2)	 ####
    #nc.variables[dst_varname_south][:] = np.squeeze(dst_var_south_2)	 ####
    #nc.variables[dst_varname_east][:] = np.squeeze(dst_var_east_2)	 ####
    #nc.variables[dst_varname_west][:] = np.squeeze(dst_var_west_2)	 ####
    nc.variables[dst_varname_north][:] = dst_var_north
    nc.variables[dst_varname_south][:] = dst_var_south
    nc.variables[dst_varname_east][:] = dst_var_east
    nc.variables[dst_varname_west][:] = dst_var_west

    # close file
    nc.close()
    cdf.close()

    if src_varname == 'ssh':
        return dst_varz
