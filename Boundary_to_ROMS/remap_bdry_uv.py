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
from netCDF4 import Dataset

import pyroms
import pyroms_toolbox
import _remapping

from flood import flood

class nctime(object):
    pass

def remap_bdry_uv(zero,l_time, src_file, src_grd, dst_grd, dxy=20, cdepth=0, kk=2, dst_dir='./'):

    # get time
    # get dimensions
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    cw=int(len(dst_grd.vgrid.Cs_r))

    # create destination file
    dst_fileu = 'u_bdry' + '.nc'
    print '\nCreating destination file', dst_fileu
    if os.path.exists(dst_fileu) is True:
        os.remove(dst_fileu)
    pyroms_toolbox.nc_create_roms_bdry_file(dst_fileu, dst_grd)
    dst_filev = 'v_bdry'  '.nc'
    print 'Creating destination file', dst_filev
    if os.path.exists(dst_filev) is True:
        os.remove(dst_filev)
    pyroms_toolbox.nc_create_roms_bdry_file(dst_filev, dst_grd)
    cdf = Dataset(src_file)
    # open destination file
    ncu = Dataset(dst_fileu, 'a', format='NETCDF4')
    ncv = Dataset(dst_filev, 'a', format='NETCDF4')

    time = cdf.variables['time'][zero:l_time]
    time=time-time[0]
    print time
    src_varu = cdf.variables['vozocrtx'][0,:,:,:]
    src_varv = cdf.variables['vomecrty'][0,:,:,:]

    #get missing value
    spval = -32767 #src_varu._FillValue

    # get weights file
    wts_file = 'remap_weights_PUSSY_to_COARSEST_bilinear_t_to_rho.nc'

    # build intermediate zgrid
    zlevel = -src_grd.z_t[::-1,0,0]
    nzlevel = len(zlevel)
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
    dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)

    # create variable in destination file

    ncu.createDimension('v3d_time', l_time-zero)			####
    ncu.createDimension('v2d_time', l_time-zero)			####

    ncu.createVariable('v3d_time', 'f8', ('v3d_time'))
    ncu.variables['v3d_time'].long_name = '3D momentum time'
    ncu.variables['v3d_time'].units = 'days since 2004-01-01 00:00:00 GMT'
    ncu.variables['v3d_time'].calendar = 'gregorian'
    #ncu.variables['v3d_time']._FillValue = spval

    ncu.createVariable('v2d_time', 'f8', ('v2d_time'))
    ncu.variables['v2d_time'].long_name = '2D momentum time'
    ncu.variables['v2d_time'].units = 'days since 2004-01-01 00:00:00 GMT'
    ncu.variables['v2d_time'].calendar = 'gregorian'
    #ncu.variables['v2d_time']._FillValue = spval

    print 'Creating variable u_north'
    ncu.createVariable('u_north', 'f8', ('v3d_time', 's_rho', 'xi_u'), fill_value=spval)
    ncu.variables['u_north'].long_name = '3D u-momentum north boundary condition'
    ncu.variables['u_north'].units = 'meter second-1'
    ncu.variables['u_north'].field = 'u_north, scalar, series'
    #ncu.variables['u_north']._FillValue = spval

    print 'Creating variable u_south'
    ncu.createVariable('u_south', 'f8', ('v3d_time', 's_rho', 'xi_u'), fill_value=spval)
    ncu.variables['u_south'].long_name = '3D u-momentum south boundary condition'
    ncu.variables['u_south'].units = 'meter second-1'
    ncu.variables['u_south'].field = 'u_south, scalar, series'
    #ncu.variables['u_south']._FillValue = spval

    print 'Creating variable u_east'
    ncu.createVariable('u_east', 'f8', ('v3d_time', 's_rho', 'eta_u'), fill_value=spval)
    ncu.variables['u_east'].long_name = '3D u-momentum east boundary condition'
    ncu.variables['u_east'].units = 'meter second-1'
    ncu.variables['u_east'].field = 'u_east, scalar, series'
    #ncu.variables['u_east']._FillValue = spval
    print 'Creating variable u_west'
    ncu.createVariable('u_west', 'f8', ('v3d_time', 's_rho', 'eta_u'), fill_value=spval)
    ncu.variables['u_west'].long_name = '3D u-momentum west boundary condition'
    ncu.variables['u_west'].units = 'meter second-1'
    ncu.variables['u_west'].field = 'u_east, scalar, series'
    #ncu.variables['u_west']._FillValue = spval

    # create variable in destination file
    print 'Creating variable ubar_north'
    ncu.createVariable('ubar_north', 'f8', ('v2d_time', 'xi_u'), fill_value=spval)
    ncu.variables['ubar_north'].long_name = '2D u-momentum north boundary condition'
    ncu.variables['ubar_north'].units = 'meter second-1'
    ncu.variables['ubar_north'].field = 'ubar_north, scalar, series'
    #ncu.variables['ubar_north']._FillValue = spval

    print 'Creating variable ubar_south'
    ncu.createVariable('ubar_south', 'f8', ('v2d_time', 'xi_u'), fill_value=spval)
    ncu.variables['ubar_south'].long_name = '2D u-momentum south boundary condition'
    ncu.variables['ubar_south'].units = 'meter second-1'
    ncu.variables['ubar_south'].field = 'ubar_south, scalar, series'
    #ncu.variables['ubar_south']._FillValue = spval

    print 'Creating variable ubar_east'
    ncu.createVariable('ubar_east', 'f8', ('v2d_time', 'eta_u'), fill_value=spval)
    ncu.variables['ubar_east'].long_name = '2D u-momentum east boundary condition'
    ncu.variables['ubar_east'].units = 'meter second-1'
    ncu.variables['ubar_east'].field = 'ubar_east, scalar, series'
    #ncu.variables['ubar_east']._FillValue = spval
    print 'Creating variable ubar_west'
    ncu.createVariable('ubar_west', 'f8', ('v2d_time', 'eta_u'), fill_value=spval)
    ncu.variables['ubar_west'].long_name = '2D u-momentum west boundary condition'
    ncu.variables['ubar_west'].units = 'meter second-1'
    ncu.variables['ubar_west'].field = 'ubar_east, scalar, series'
    #ncu.variables['ubar_west']._FillValue = spval

    ncv.createDimension('v3d_time', l_time-zero)			###
    ncv.createDimension('v2d_time', l_time-zero)			###

    ncv.createVariable('v3d_time', 'f8', ('v3d_time'))
    ncv.variables['v3d_time'].long_name = '3D momentum time'
    ncv.variables['v3d_time'].units = 'days since 2004-01-01 00:00:00 GMT'
    ncv.variables['v3d_time'].calendar = 'gregorian'
    #ncu.variables['v3d_time']._FillValue = spval

    ncv.createVariable('v2d_time', 'f8', ('v2d_time'))
    ncv.variables['v2d_time'].long_name = '2D momentum time'
    ncv.variables['v2d_time'].units = 'days since 2004-01-01 00:00:00 GMT'
    ncv.variables['v2d_time'].calendar = 'gregorian'
    #ncu.variables['v2d_time']._FillValue = spval


    print 'Creating variable v_north'
    ncv.createVariable('v_north', 'f8', ('v3d_time', 's_rho', 'xi_v'), fill_value=spval)
    ncv.variables['v_north'].long_name = '3D v-momentum north boundary condition'
    ncv.variables['v_north'].units = 'meter second-1'
    ncv.variables['v_north'].field = 'v_north, scalar, series'
    #ncv.variables['v_north']._FillValue = spval

    print 'Creating variable v_south'
    ncv.createVariable('v_south', 'f8', ('v3d_time', 's_rho', 'xi_v'), fill_value=spval)
    ncv.variables['v_south'].long_name = '3D v-momentum south boundary condition'
    ncv.variables['v_south'].units = 'meter second-1'
    ncv.variables['v_south'].field = 'v_south, scalar, series'
    #ncv.variables['v_south']._FillValue = spval

    print 'Creating variable v_east'
    ncv.createVariable('v_east', 'f8', ('v3d_time', 's_rho', 'eta_v'), fill_value=spval)
    ncv.variables['v_east'].long_name = '3D v-momentum east boundary condition'
    ncv.variables['v_east'].units = 'meter second-1'
    ncv.variables['v_east'].field = 'v_east, scalar, series'
    #ncv.variables['v_east']._FillValue = spval
    print 'Creating variable v_west'
    ncv.createVariable('v_west', 'f8', ('v3d_time', 's_rho', 'eta_v'), fill_value=spval)
    ncv.variables['v_west'].long_name = '3D v-momentum west boundary condition'
    ncv.variables['v_west'].units = 'meter second-1'
    ncv.variables['v_west'].field = 'v_east, scalar, series'
    #ncv.variables['v_west']._FillValue = spval

    print 'Creating variable vbar_north'
    ncv.createVariable('vbar_north', 'f8', ('v2d_time', 'xi_v'), fill_value=spval)
    ncv.variables['vbar_north'].long_name = '2D v-momentum north boundary condition'
    ncv.variables['vbar_north'].units = 'meter second-1'
    ncv.variables['vbar_north'].field = 'vbar_north, scalar, series'
    #ncv.variables['vbar_north']._FillValue = spval

    print 'Creating variable vbar_south'
    ncv.createVariable('vbar_south', 'f8', ('v2d_time', 'xi_v'), fill_value=spval)
    ncv.variables['vbar_south'].long_name = '2D v-momentum south boundary condition'
    ncv.variables['vbar_south'].units = 'meter second-1'
    ncv.variables['vbar_south'].field = 'vbar_south, scalar, series'
    #ncv.variables['vbar_south']._FillValue = spval

    print 'Creating variable vbar_east'
    ncv.createVariable('vbar_east', 'f8', ('v2d_time', 'eta_v'), fill_value=spval)
    ncv.variables['vbar_east'].long_name = '2D v-momentum east boundary condition'
    ncv.variables['vbar_east'].units = 'meter second-1'
    ncv.variables['vbar_east'].field = 'vbar_east, scalar, series'
    #ncv.variables['vbar_east']._FillValue = spval
    print 'Creating variable vbar_west'
    ncv.createVariable('vbar_west', 'f8', ('v2d_time', 'eta_v'), fill_value=spval)
    ncv.variables['vbar_west'].long_name = '2D v-momentum west boundary condition'
    ncv.variables['vbar_west'].units = 'meter second-1'
    ncv.variables['vbar_west'].field = 'vbar_east, scalar, series'
    #ncv.variables['vbar_west']._FillValue = spval





    dst_u_north=np.ma.zeros((l_time-zero, cw, Lp-1))
    dst_u_south=np.ma.zeros((l_time-zero, cw, Lp-1))
    dst_u_east=np.ma.zeros((l_time-zero, cw, Mp))
    dst_u_west=np.ma.zeros((l_time-zero, cw, Mp))

    dst_v_north=np.ma.zeros((l_time-zero, cw,  Lp))
    dst_v_south=np.ma.zeros((l_time-zero, cw,  Lp))
    dst_v_east=np.ma.zeros((l_time-zero, cw, Mp-1))
    dst_v_west=np.ma.zeros((l_time-zero, cw, Mp-1))

    dst_ubar_north=np.ma.zeros((l_time-zero, Lp-1))
    dst_ubar_south=np.ma.zeros((l_time-zero, Lp-1))
    dst_ubar_east=np.ma.zeros((l_time-zero, Mp))
    dst_ubar_west=np.ma.zeros((l_time-zero, Mp))


    dst_vbar_north=np.ma.zeros((l_time-zero, Lp))
    dst_vbar_south=np.ma.zeros((l_time-zero, Lp))
    dst_vbar_east=np.ma.zeros((l_time-zero, Mp-1))
    dst_vbar_west=np.ma.zeros((l_time-zero, Mp-1))

    for m in range(len(time)):
        src_uz = flood(cdf.variables['vozocrtx'][m,:,:,:], src_grd,  spval=spval)
        src_vz = flood(cdf.variables['vomecrty'][m,:,:,:], src_grd,  spval=spval)
        dst_uz = pyroms.remapping.remap(src_uz, wts_file, spval=spval)
        dst_vz = pyroms.remapping.remap(src_vz, wts_file, spval=spval)

	dst_u_north_cont = pyroms.remapping.z2roms(dst_uz[::-1, Mp-2:Mp, 0:Lp], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=True, irange=(0,Lp), jrange=(Mp-2,Mp))
	#print np.shape(dst_u_north)
	dst_u_south_cont = pyroms.remapping.z2roms(dst_uz[::-1, 0:2, 0:Lp], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=True, irange=(0,Lp), jrange=(0,2))
	dst_u_east_cont = pyroms.remapping.z2roms(dst_uz[::-1, 0:Mp, Lp-2:Lp], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=True, irange=(Lp-2,Lp), jrange=(0,Mp))
	dst_u_west_cont = pyroms.remapping.z2roms(dst_uz[::-1, 0:Mp, 0:2], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=True, irange=(0,2), jrange=(0,Mp))

	dst_v_north_cont = pyroms.remapping.z2roms(dst_vz[::-1, Mp-2:Mp, 0:Lp], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=True, irange=(0,Lp), jrange=(Mp-2,Mp))
	dst_v_south_cont = pyroms.remapping.z2roms(dst_vz[::-1, 0:2, 0:Lp], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=True, irange=(0,Lp), jrange=(0,2))
	dst_v_east_cont = pyroms.remapping.z2roms(dst_vz[::-1, 0:Mp, Lp-2:Lp], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=True, irange=(Lp-2,Lp), jrange=(0,Mp))
	dst_v_west_cont = pyroms.remapping.z2roms(dst_vz[::-1, 0:Mp, 0:2], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=True, irange=(0,2), jrange=(0,Mp))
	#print dst_u_south_cont.min(), dst_u_south_cont.max()

	src_angle = pyroms.remapping.remap(src_grd.angle, 'remap_weights_PUSSY_to_COARSEST_bilinear_t_to_rho.nc', spval=spval)
	dst_angle = dst_grd.hgrid.angle_rho
	angle = dst_angle - src_angle
	angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))

	U_north = dst_u_north_cont + dst_v_north_cont*1j
	eitheta_north = np.exp(-1j*angle[:,Mp-2:Mp, 0:Lp])
	U_north = U_north * eitheta_north
	dst_u_north_cont = np.real(U_north)
	dst_v_north_cont = np.imag(U_north)

	U_south = dst_u_south_cont + dst_v_south_cont*1j
	eitheta_south = np.exp(-1j*angle[:,0:2, 0:Lp])
	U_south = U_south * eitheta_south
	dst_u_south_cont = np.real(U_south)
	dst_v_south_cont = np.imag(U_south)

	U_east = dst_u_east_cont + dst_v_east_cont*1j
	eitheta_east = np.exp(-1j*angle[:,0:Mp, Lp-2:Lp])
	U_east = U_east * eitheta_east
	dst_u_east_cont = np.real(U_east)
	dst_v_east_cont = np.imag(U_east)

	U_west = dst_u_west_cont + dst_v_west_cont*1j
	eitheta_west = np.exp(-1j*angle[:,0:Mp, 0:2])
	U_west = U_west * eitheta_west
	dst_u_west_cont = np.real(U_west)
	dst_v_west_cont = np.imag(U_west)

	dst_u_north_cont = 0.5 * np.squeeze(dst_u_north_cont[:,-1,:-1] + dst_u_north_cont[:,-1,1:])
	dst_v_north_cont = 0.5 * np.squeeze(dst_v_north_cont[:,:-1,:] + dst_v_north_cont[:,1:,:])
	dst_u_south_cont = 0.5 * np.squeeze(dst_u_south_cont[:,0,:-1] + dst_u_south_cont[:,0,1:])
	dst_v_south_cont = 0.5 * np.squeeze(dst_v_south_cont[:,:-1,:] + dst_v_south_cont[:,1:,:])
	dst_u_east_cont = 0.5 * np.squeeze(dst_u_east_cont[:,:,:-1] + dst_u_east_cont[:,:,1:])
	dst_v_east_cont = 0.5 * np.squeeze(dst_v_east_cont[:,:-1,-1] + dst_v_east_cont[:,1:,-1])
	dst_u_west_cont = 0.5 * np.squeeze(dst_u_west_cont[:,:,:-1] + dst_u_west_cont[:,:,1:])
	dst_v_west_cont = 0.5 * np.squeeze(dst_v_west_cont[:,:-1,0] + dst_v_west_cont[:,1:,0])
	#print dst_u_south_cont.min(), dst_u_south_cont.max()

	idxu_north = np.where(dst_grd.hgrid.mask_u[-1,:] == 0)
	idxv_north = np.where(dst_grd.hgrid.mask_v[-1,:] == 0)
	idxu_south = np.where(dst_grd.hgrid.mask_u[0,:] == 0)
	idxv_south = np.where(dst_grd.hgrid.mask_v[0,:] == 0)
	idxu_east = np.where(dst_grd.hgrid.mask_u[:,-1] == 0)
	idxv_east = np.where(dst_grd.hgrid.mask_v[:,-1] == 0)
	idxu_west = np.where(dst_grd.hgrid.mask_u[:,0] == 0)
	idxv_west = np.where(dst_grd.hgrid.mask_v[:,0] == 0)
	for n in range(dst_grd.vgrid.N): 
		dst_u_north_cont[n, idxu_north[0]] = spval
        	dst_v_north_cont[n, idxv_north[0]] = spval
        	dst_u_south_cont[n, idxu_south[0]] = spval
        	dst_v_south_cont[n, idxv_south[0]] = spval
        	dst_u_east_cont[n, idxu_east[0]] = spval
        	dst_v_east_cont[n, idxv_east[0]] = spval
        	dst_u_west_cont[n, idxu_west[0]] = spval
        	dst_v_west_cont[n, idxv_west[0]] = spval

	z_u_north = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:-1] + dst_grd.vgrid.z_w[0,:,-1,1:])
	z_v_north = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:] + dst_grd.vgrid.z_w[0,:,-2,:])
	z_u_south = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:-1] + dst_grd.vgrid.z_w[0,:,0,1:])
	z_v_south = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:] + dst_grd.vgrid.z_w[0,:,1,:])
	z_u_east = 0.5 * (dst_grd.vgrid.z_w[0,:,:,-1] + dst_grd.vgrid.z_w[0,:,:,-2])
	z_v_east = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,-1] + dst_grd.vgrid.z_w[0,:,1:,-1])
	z_u_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:,0] + dst_grd.vgrid.z_w[0,:,:,1])
	z_v_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,0] + dst_grd.vgrid.z_w[0,:,1:,0])


	dst_ubar_north_cont = np.ma.zeros(dst_u_north_cont.shape[1])
	dst_ubar_south_cont = np.ma.zeros(dst_u_south_cont.shape[1])
	dst_ubar_east_cont = np.ma.zeros(dst_u_east_cont.shape[1])
	dst_ubar_west_cont = np.ma.zeros(dst_u_west_cont.shape[1])
	dst_vbar_north_cont = np.ma.zeros(dst_v_north_cont.shape[1])
	dst_vbar_south_cont = np.ma.zeros(dst_v_south_cont.shape[1])
	dst_vbar_east_cont = np.ma.zeros(dst_v_east_cont.shape[1])
	dst_vbar_west_cont = np.ma.zeros(dst_v_west_cont.shape[1])


	for i in range(dst_u_north_cont.shape[1]):
        	dst_ubar_north_cont[i] = (dst_u_north_cont[:,i] * np.diff(z_u_north[:,i])).sum() / -z_u_north[0,i]
	for i in range(dst_v_north_cont.shape[1]):
        	dst_vbar_north_cont[i] = (dst_v_north_cont[:,i] * np.diff(z_v_north[:,i])).sum() / -z_v_north[0,i]
	for i in range(dst_u_south_cont.shape[1]):
        	dst_ubar_south_cont[i] = (dst_u_south_cont[:,i] * np.diff(z_u_south[:,i])).sum() / -z_u_south[0,i]
	for i in range(dst_v_south_cont.shape[1]):
        	dst_vbar_south_cont[i] = (dst_v_south_cont[:,i] * np.diff(z_v_south[:,i])).sum() / -z_v_south[0,i]
	for j in range(dst_u_east_cont.shape[1]):
        	dst_ubar_east_cont[j] = (dst_u_east_cont[:,j] * np.diff(z_u_east[:,j])).sum() / -z_u_east[0,j]
        	dst_ubar_west_cont[j] = (dst_u_west_cont[:,j] * np.diff(z_u_west[:,j])).sum() / -z_u_west[0,j]
	for j in range(dst_v_east_cont.shape[1]):
        	dst_vbar_east_cont[j] = (dst_v_east_cont[:,j] * np.diff(z_v_east[:,j])).sum() / -z_v_east[0,j]
        	dst_vbar_west_cont[j] = (dst_v_west_cont[:,j] * np.diff(z_v_west[:,j])).sum() / -z_v_west[0,j]

	#print dst_ubar_south_cont.min()

	dst_ubar_north[m,:] = np.ma.masked_where(dst_grd.hgrid.mask_u[-1,:] == 0, dst_ubar_north_cont)   
	dst_ubar_south[m,:] = np.ma.masked_where(dst_grd.hgrid.mask_u[0,:] == 0, dst_ubar_south_cont)   
	dst_ubar_east[m,:] = np.ma.masked_where(dst_grd.hgrid.mask_u[:,-1] == 0, dst_ubar_east_cont)   
	dst_ubar_west[m,:] = np.ma.masked_where(dst_grd.hgrid.mask_u[:,0] == 0, dst_ubar_west_cont)   
	dst_vbar_north[m,:] = np.ma.masked_where(dst_grd.hgrid.mask_v[-1,:] == 0, dst_vbar_north_cont)   
	dst_vbar_south[m,:] = np.ma.masked_where(dst_grd.hgrid.mask_v[0,:] == 0, dst_vbar_south_cont)   
	dst_vbar_east[m,:] = np.ma.masked_where(dst_grd.hgrid.mask_v[:,-1] == 0, dst_vbar_east_cont)   
	dst_vbar_west[m,:] = np.ma.masked_where(dst_grd.hgrid.mask_v[:,0] == 0, dst_vbar_west_cont)  


	#dst_ubar_north[m,:]=dst_ubar_north_cont
	#dst_ubar_south[m,:]=dst_ubar_south_cont
	#dst_ubar_east[m,:]=dst_ubar_east_cont
	#dst_ubar_west[m,:]=dst_ubar_west_cont

	#dst_vbar_north[m,:]=dst_vbar_north_cont
	#dst_vbar_south[m,:]=dst_vbar_south_cont
	#dst_vbar_east[m,:]=dst_vbar_east_cont
	#dst_vbar_west[m,:]=dst_vbar_west_cont

        #print m
	dst_u_north[m,:,:]=np.ma.masked_where(np.tile(dst_grd.hgrid.mask_u[-1,:],(cw,1)) == 0, dst_u_north_cont)
	dst_u_south[m,:,:]=np.ma.masked_where(np.tile(dst_grd.hgrid.mask_u[0,:],(cw,1)) == 0,dst_u_south_cont)
	dst_u_east[m,:,:]=np.ma.masked_where(np.tile(dst_grd.hgrid.mask_u[:,-1],(cw,1)) == 0,dst_u_east_cont)
	dst_u_west[m,:,:]=np.ma.masked_where(np.tile(dst_grd.hgrid.mask_u[:,0],(cw,1)) == 0,dst_u_west_cont)

	dst_v_north[m,:,:]=np.ma.masked_where(np.tile(dst_grd.hgrid.mask_v[-1,:],(cw,1)) == 0,dst_v_north_cont)
	dst_v_south[m,:,:]=np.ma.masked_where(np.tile(dst_grd.hgrid.mask_v[0,:],(cw,1)) == 0,dst_v_south_cont)
	dst_v_east[m,:,:]=np.ma.masked_where(np.tile(dst_grd.hgrid.mask_v[:,-1],(cw,1)) == 0,dst_v_east_cont)
	dst_v_west[m,:,:]=np.ma.masked_where(np.tile(dst_grd.hgrid.mask_v[:,0],(cw,1)) == 0,dst_v_west_cont)
	print m, dst_u_south[m].min(), dst_u_south[m].max(), dst_ubar_south[m,:].min(), dst_ubar_south[m,:].max()#, dst_ubar_south[-1], sum(dst_u_south[:,-1])
	#print np.shape(dst_u_south_cont)
	#if not int(dst_ubar_south[m,:].min()) == spval:
	#	break


    # write data in destination file
    print 'write data in destination file'
    ncu.variables['v2d_time'][:] = time####np.concatenate((time,time+365,time+365+365,time+365+365+365,time+365+365+365+365),axis=0)*60*60*24####
    ncu.variables['v3d_time'][:] = time####np.concatenate((time,time+365,time+365+365,time+365+365+365,time+365+365+365+365),axis=0)*60*60*24####
    ncu.variables['u_north'][:] = dst_u_north####np.concatenate((dst_u_north,dst_u_north,dst_u_north,dst_u_north,dst_u_north),axis=0) ####dst_u_north####
    ncu.variables['u_south'][:] = dst_u_south####np.concatenate((dst_u_south,dst_u_south,dst_u_south,dst_u_south,dst_u_south),axis=0)####dst_u_south####
    ncu.variables['u_east'][:] = dst_u_east####np.concatenate((dst_u_east,dst_u_east,dst_u_east,dst_u_east,dst_u_east),axis=0)####dst_u_east####
    ncu.variables['u_west'][:] = dst_u_west####np.concatenate((dst_u_west,dst_u_west,dst_u_west,dst_u_west,dst_u_west),axis=0)####dst_u_west####
    ncu.variables['ubar_north'][:] = dst_ubar_north####np.concatenate((dst_ubar_north,dst_ubar_north,dst_ubar_north,dst_ubar_north,dst_ubar_north),axis=0)####dst_ubar_north####
    ncu.variables['ubar_south'][:] = dst_ubar_south####np.concatenate((dst_ubar_south,dst_ubar_south,dst_ubar_south,dst_ubar_south,dst_ubar_south),axis=0)####dst_ubar_south####
    ncu.variables['ubar_east'][:] = dst_ubar_east####np.concatenate((dst_ubar_east,dst_ubar_east,dst_ubar_east,dst_ubar_east,dst_ubar_east),axis=0)####dst_ubar_east####
    ncu.variables['ubar_west'][:] = dst_ubar_west####np.concatenate((dst_ubar_west,dst_ubar_west,dst_ubar_west,dst_ubar_west,dst_ubar_west),axis=0)####dst_ubar_west####

    ncv.variables['v2d_time'][:] = time###np.concatenate((time,time+365,time+365+365,time+365+365+365,time+365+365+365+365),axis=0)*60*60*24####
    ncv.variables['v3d_time'][:] = time###np.concatenate((time,time+365,time+365+365,time+365+365+365,time+365+365+365+365),axis=0)*60*60*24####
    ncv.variables['v_north'][:] = dst_v_north####np.concatenate((dst_v_north,dst_v_north,dst_v_north,dst_v_north,dst_v_north),axis=0)####dst_v_north####
    ncv.variables['v_south'][:] = dst_v_south####np.concatenate((dst_v_south,dst_v_south,dst_v_south,dst_v_south,dst_v_south),axis=0)####dst_v_south####
    ncv.variables['v_east'][:] = dst_v_east####np.concatenate((dst_v_east,dst_v_east,dst_v_east,dst_v_east,dst_v_east),axis=0)####dst_v_east####
    ncv.variables['v_west'][:] = dst_v_west####np.concatenate((dst_v_west,dst_v_west,dst_v_west,dst_v_west,dst_v_west),axis=0)####dst_v_west####
    ncv.variables['vbar_north'][:] = dst_vbar_north####np.concatenate((dst_vbar_north,dst_vbar_north,dst_vbar_north,dst_vbar_north,dst_vbar_north),axis=0)####dst_vbar_north####
    ncv.variables['vbar_south'][:] = dst_vbar_south####np.concatenate((dst_vbar_south,dst_vbar_south,dst_vbar_south,dst_vbar_south,dst_vbar_south),axis=0)####dst_vbar_south####
    ncv.variables['vbar_east'][:] = dst_vbar_east####np.concatenate((dst_vbar_east,dst_vbar_east,dst_vbar_east,dst_vbar_east,dst_vbar_east),axis=0)####dst_vbar_east####
    ncv.variables['vbar_west'][:] = dst_vbar_west####np.concatenate((dst_vbar_west,dst_vbar_west,dst_vbar_west,dst_vbar_west,dst_vbar_west),axis=0)####dst_vbar_west####

    # close file
    ncu.close()
    ncv.close()
    cdf.close()
