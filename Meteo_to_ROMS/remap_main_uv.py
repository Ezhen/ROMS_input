from netCDF4 import Dataset
import numpy as np
from numpy import shape
import pyroms
import pyroms_toolbox
from datetime import datetime
import subprocess
import os
import commands
import _iso
from nco import Nco
import _remapping
from get_nc_Grid_Nest import get_nc_Grid_Nest
from get_nc_Grid_HYCOM import get_nc_Grid_HYCOM
from nc_create_roms_bdry_file_mine import nc_create_roms_bdry_file_mine


def remap_main_uv(number,tart, tend, src_file, src_grd, dst_grd, grid_name):
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    dxy=5; cdepth=0; kk=0; Bpos = 't'; Cpos = 'rho'; spval = -32767
    cdf = Dataset(src_file)
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    if number==1:
    	dst_file = 'wind%s.nc' %(grid_name)
    	nc_create_roms_bdry_file_mine(dst_file, dst_grd)
	print 'lol'
    	cdf = Dataset(src_file, 'r', format='NETCDF4')
    	l_time=tend-tart
    	time = cdf.variables['time'][tart:tend]
	print time[0], time[-1], shape(time)
 	time=time-time[0]
	time = time
    	time[0]=0
    	src_varu = cdf.variables['u10'][0,:,:]
    	src_varv = cdf.variables['v10'][0,:,:]
  	wts_file = 'remap_weights_PUSSY_to_%s_bilinear_t_to_rho.nc' %(grid_name)
	nc = Dataset(dst_file, 'a', format='NETCDF4')
   	nc.createDimension('wind_time', l_time)				###
  	nc.createVariable('wind_time', 'f4', 'wind_time')
 	nc.variables['wind_time'].long_name = 'wind time'
  	nc.variables['wind_time'].units = 'days since 2006-01-01 00:00:00'
  	nc.variables['wind_time'].calendar = 'LEAP'
  	nc.variables['wind_time'].field = 'wind_time field'
        #time_2=np.concatenate((time,time+365,time+365+365,time+365+365+365,time+365+365+365+365),axis=0)				###
	time_2=time
  	nc.variables['wind_time'][:] = time_2/24.
  	nc.createVariable('Uwind', 'f4', ('wind_time', 'er', 'xr'), fill_value=spval)
  	nc.variables['Uwind'].long_name = 'surface u-wind component'
  	nc.variables['Uwind'].units = 'meter second-1'
 	nc.variables['Uwind'].field = 'Uwind, scalar, series'
   	nc.variables['Uwind'].coordinates = 'lon lat'
	nc.createVariable('Vwind', 'f4', ('wind_time', 'er', 'xr'), fill_value=spval)
 	nc.variables['Vwind'].long_name = 'surface v-wind component'
 	nc.variables['Vwind'].units = 'meter second-1'
 	nc.variables['Vwind'].field = 'Vwind, scalar, series'
 	nc.variables['Vwind'].coordinates = 'lon lat'
	src_angle = pyroms.remapping.remap(src_grd.angle,  'remap_weights_PUSSY_to_%s_bilinear_t_to_rho.nc' %(grid_name), spval=spval)
  	dst_angle = dst_grd.hgrid.angle_rho
  	angle = dst_angle - src_angle
   	angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
   	eitheta = np.exp(-1j*angle[:,:,:])
   	dst_u=np.zeros((l_time, Mp,Lp))
   	dst_v=np.zeros((l_time, Mp,Lp))
   	for m in range(len(time)):
        	dst_uz = pyroms.remapping.remap(cdf.variables['u10'][m,:,:], wts_file, spval=spval)
        	dst_vz = pyroms.remapping.remap(cdf.variables['v10'][m,:,:], wts_file, spval=spval)
        	U = dst_uz + dst_vz*1j      
		U = U * eitheta
        	dst_uz = np.real(U)
        	dst_vz = np.imag(U)
        	dst_u_cont = 0.5 * (dst_uz[:,:] + dst_uz[:,:])
       		dst_v_cont = 0.5 * (dst_vz[:,:] + dst_vz[:,:])
		dst_u[m,:,:]=dst_u_cont[0,:,:]
		dst_v[m,:,:]=dst_v_cont[0,:,:]
		print m
        #dst_u_2=np.concatenate((dst_u,dst_u,dst_u,dst_u,dst_u),axis=0)				###
	dst_u_2=dst_u
        #dst_v_2=np.concatenate((dst_v,dst_v,dst_v,dst_v,dst_v),axis=0)				###
	dst_v_2=dst_v
    	#nc.variables['Uwind'][:] = np.around(dst_u,1)
    	nc.variables['Uwind'][:] = dst_u_2 					###
   	#nc.variables['Vwind'][:] = np.around(dst_v,1)
   	nc.variables['Vwind'][:] = dst_v_2 					###
	print 'surface wind is written into the separated netcdf file'
   	nc.close()
    if number==2:
    	dst_file = 'momentum%s.nc' %(grid_name)
    	#pyroms_toolbox.nc_create_roms_bdry_file(dst_file, dst_grd)
	nc_create_roms_bdry_file_mine(dst_file, dst_grd)
    	cdf = Dataset(src_file, 'r', format='NETCDF4')
    	l_time=tend-tart
    	time = (cdf.variables['time'][tart:tend])
 	time=time-time[0]
    	time[0]=0
    	src_varu = cdf.variables['nsss'][tart:tend,:,:]
    	src_varv = cdf.variables['ewss'][tart:tend,:,:]
  	wts_file = 'remap_weights_PUSSY_to_%s_bilinear_t_to_rho.nc'
	nc = Dataset(dst_file, 'a', format='NETCDF4')
   	nc.createDimension('sms_time', l_time)
  	nc.createVariable('sms_time', 'f4', 'sms_time')
 	nc.variables['sms_time'].long_name = 'wind time'
  	nc.variables['sms_time'].units = 'days since 2006-01-01 00:00:00'
  	nc.variables['sms_time'].field = 'sms_time field'
	print len(time)
  	nc.variables['sms_time'][:] = time/24.
  	nc.createVariable('sustr', 'f4', ('sms_time', 'eu', 'xu'), fill_value=spval)
  	nc.variables['sustr'].long_name = 'surface u-momentum stress'
  	nc.variables['sustr'].units = 'newton meter-2'
 	nc.variables['sustr'].field = 'surface u-momentum stress, scalar, series'
   	nc.variables['sustr'].coordinates = 'lon lat'
	#nc.variables['sustr'].coordinates = 'lon lat'
	nc.createVariable('svstr', 'f4', ('sms_time', 'ev', 'xv'), fill_value=spval)
 	nc.variables['svstr'].long_name = 'surface v-momentum stress'
 	nc.variables['svstr'].units = 'newton meter-2'
 	nc.variables['svstr'].field = 'surface v-momentum stress, scalar, series'
   	nc.variables['svstr'].coordinates = 'lon lat'
	#nc.variables['svstr'].coordinates = 'lon lat'
	src_angle = pyroms.remapping.remap(src_grd.angle,  'remap_weights_PUSSY_to_%s_bilinear_t_to_rho.nc' %(grid_name), spval=spval)
  	dst_angle = dst_grd.hgrid.angle_rho
  	angle = dst_angle - src_angle
   	angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
   	eitheta = np.exp(-1j*angle[0,:,:])#[:,:,:])
   	dst_u=np.zeros((l_time, Mp,Lp-1))
   	dst_v=np.zeros((l_time, Mp-1,Lp))
   	for m in range(len(time)):
        	dst_uz = pyroms.remapping.remap(src_varu[m,:,:], wts_file, spval=spval)
        	dst_vz = pyroms.remapping.remap(src_varv[m,:,:], wts_file, spval=spval)
        	U = dst_uz + dst_vz*1j      
		U = U * eitheta
        	dst_uz = np.real(U)
        	dst_vz = np.imag(U)
        	dst_u_cont = 0.5 * (dst_uz[:,:-1] + dst_uz[:,1:])
       		dst_v_cont = 0.5 * (dst_vz[:-1,:] + dst_vz[1:,:])
		dst_u[m,:,:]=dst_u_cont[:,:]
		dst_v[m,:,:]=dst_v_cont[:,:]
		print 'sustr, svstr', m	
    	#nc.variables['sustr'][:] = np.around(dst_u,1)
    	nc.variables['sustr'][:] = dst_u
   	#nc.variables['svstr'][:] = np.around(dst_v,1)
   	nc.variables['svstr'][:] = dst_v
	print 'surface momentum stress is written into the separated netcdf file'
   	nc.close()
