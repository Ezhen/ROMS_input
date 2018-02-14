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
from get_nc_Grid_HYCOM import get_nc_Grid_HYCOM
from flood import flood

class nctime(object):
    pass

def remap_uv(zero,l_time, src_file, src_grd, dst_grd, dmax=20, cdepth=0, kk=0, dst_dir='./'):
    Mp,Np = dst_grd.vgrid.h.shape
    dxy=5
    cdepth=0
    kk=0
    nctime.long_name = 'time'
    nctime.units = 'seconds since 2006-01-01 00:00:00'
    cdf = Dataset(src_file)
 
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    dst_file = src_file.rsplit('/')[-1]
    dst_fileu = 'u' + dst_grd.name + '.nc'
    if os.path.exists(dst_fileu) is True:
        os.remove(dst_fileu)
    pyroms_toolbox.nc_create_roms_file_3(dst_fileu, dst_grd, nctime)
    dst_filev =  'v' + dst_grd.name + '.nc'
    if os.path.exists(dst_filev) is True:
        os.remove(dst_filev)
    pyroms_toolbox.nc_create_roms_file_3(dst_filev, dst_grd, nctime)

    # open destination file
    ncu = Dataset(dst_fileu, 'a', format='NETCDF3_64BIT')
    ncv = Dataset(dst_filev, 'a', format='NETCDF3_64BIT')

    time = cdf.variables['time'][zero:l_time]
    time=time-time[0]
    print time

    #get missing value
    spval = -32767 #src_varu._FillValue

    # get weights file
    wts_file = 'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc'

    # build intermediate zgrid
    zlevel = -src_grd.z_t[::-1,0,0]
    nzlevel = len(zlevel)
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
    dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)

    # create variable in destination file
    #print 'Creating variable u'
    ncu.createVariable('u', 'f8', ('ocean_time', 'sc_r', 'eu', 'xu'), fill_value=spval)
    ncu.variables['u'].long_name = '3D u-momentum component'
    ncu.variables['u'].units = 'meter second-1'
    ncu.variables['u'].field = 'u-velocity, scalar, series'
    #ncu.variables['u_north']._FillValue = spval
    # create variable in destination file
    #print 'Creating variable ubar'
    ncu.createVariable('ubar', 'f8', ('ocean_time', 'eu', 'xu'), fill_value=spval)
    ncu.variables['ubar'].long_name = '2D u-momentum component'
    ncu.variables['ubar'].units = 'meter second-1'
    ncu.variables['ubar'].field = 'ubar-velocity,, scalar, series'
    #ncu.variables['ubar_north']._FillValue = spval

    #print 'Creating variable v'
    ncv.createVariable('v', 'f8', ('ocean_time', 'sc_r', 'ev', 'xv'), fill_value=spval)
    ncv.variables['v'].long_name = '3D v-momentum component'
    ncv.variables['v'].units = 'meter second-1'
    ncv.variables['v'].field = 'v-velocity, scalar, series'
    #ncv.variables['v_north']._FillValue = spval
    #print 'Creating variable vbar'
    ncv.createVariable('vbar', 'f8', ('ocean_time', 'ev', 'xv'), fill_value=spval)
    ncv.variables['vbar'].long_name = '2D v-momentum component'
    ncv.variables['vbar'].units = 'meter second-1'
    ncv.variables['vbar'].field = 'vbar-velocity,, scalar, series'
    #ncv.variables['vbar_north']._FillValue = spval

    # remaping
    #print 'remapping and rotating u and v from', src_grd.name, 'to', dst_grd.name
    #print 'time =', time

    src_angle = pyroms.remapping.remap(src_grd.angle,  'remap_weights_PUSSY_to_FINER_bilinear_t_to_rho.nc', spval=spval)
    dst_angle = dst_grd.hgrid.angle_rho
    angle = dst_angle - src_angle
    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
    eitheta = np.exp(-1j*angle[:,:,:])


    dst_u=np.zeros((l_time-zero, 15, Mp, Np-1))
    dst_v=np.zeros((l_time-zero, 15, Mp-1, Np))
    dst_ubar=np.zeros((l_time-zero, Mp, Np-1))
    dst_vbar=np.zeros((l_time-zero, Mp-1, Np))
    for m in range(len(time)):
        src_uz = flood(cdf.variables['vomecrty'][zero+m,:,:,:], src_grd,  spval=spval)
        src_vz = flood(cdf.variables['vozocrtx'][zero+m,:,:,:], src_grd,  spval=spval)
        dst_uz = pyroms.remapping.remap(src_uz, wts_file, spval=spval)
        dst_vz = pyroms.remapping.remap(src_vz, wts_file, spval=spval)
        dst_u_cont = pyroms.remapping.z2roms(dst_uz[::-1,:,:], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=True)
        dst_v_cont = pyroms.remapping.z2roms(dst_vz[::-1,:,:], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=True)
        U = dst_u_cont + dst_v_cont*1j
        U = U * eitheta
        dst_u_cont = np.real(U)
        dst_v_cont = np.imag(U)
        dst_u_cont = 0.5 * (dst_u_cont[:,:,:-1] + dst_u_cont[:,:,1:])
        dst_v_cont = 0.5 * (dst_v_cont[:,:-1,:] + dst_v_cont[:,1:,:])
        idxu = np.where(dst_grd.hgrid.mask_u == 0)
        idxv = np.where(dst_grd.hgrid.mask_v == 0)
        for n in range(dst_grd.vgrid.N):
        	dst_u_cont[n,idxu[0], idxu[1]] = spval
        	dst_v_cont[n,idxv[0], idxv[1]] = spval
        z_u = 0.5 * (dst_grd.vgrid.z_w[0,:,:,:-1] + dst_grd.vgrid.z_w[0,:,:,1:])
        z_v = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,:] + dst_grd.vgrid.z_w[0,:,1:,:])
        dst_ubar_cont = np.zeros((dst_u_cont.shape[1], dst_u_cont.shape[2]))
        dst_vbar_cont = np.zeros((dst_v_cont.shape[1], dst_v_cont.shape[2]))
        for i in range(dst_ubar_cont.shape[1]):
      		  for j in range(dst_ubar_cont.shape[0]):
			if z_u[0,j,i] !=0:
        			dst_ubar_cont[j,i] = (dst_u_cont[:,j,i] * np.diff(z_u[:,j,i])).sum() / -z_u[0,j,i]
			else:	
				dst_ubar_cont[j,i] = dst_ubar_cont[j,i-1]
        for i in range(dst_vbar_cont.shape[1]):
        	for j in range(dst_vbar_cont.shape[0]):
			if z_v[0,j,i] !=0:
        			dst_vbar_cont[j,i] = (dst_v_cont[:,j,i] * np.diff(z_v[:,j,i])).sum() / -z_v[0,j,i]
			else:
				dst_vbar_cont[j,i]= dst_vbar_cont[j,i-1]
        for i in range(dst_ubar_cont.shape[1]):
        	for j in range(dst_ubar_cont.shape[0]):
			if z_u[0,j,i] !=0:
        			dst_ubar_cont[j,i] = (dst_u_cont[:,j,i] * np.diff(z_u[:,j,i])).sum() / -z_u[0,j,i]
			else:
				dst_ubar_cont[j,i]=dst_ubar_cont[j,i-1]
        for i in range(dst_vbar_cont.shape[1]):
        	for j in range(dst_vbar_cont.shape[0]):
			if z_v[0,j,i] !=0:
        			dst_vbar_cont[j,i] = (dst_v_cont[:,j,i] * np.diff(z_v[:,j,i])).sum() / -z_v[0,j,i]
			else:
				dst_vbar_cont[j,i]= dst_vbar_cont[j,i-1]
        dst_ubar_cont[idxu[0], idxu[1]] = spval
        dst_vbar_cont[idxv[0], idxv[1]] = spval

	dst_ubar[m,:,:]=dst_ubar_cont
	dst_vbar[m,:,:]=dst_vbar_cont
	dst_u[m,:,:,:]=dst_u_cont
	dst_v[m,:,:,:]=dst_v_cont


    # write data in destination file
    #print 'write data in destination file'
    ncu.variables['ocean_time'][:] = time/24.
    ncu.variables['u'][:] = dst_u
    ncu.variables['ubar'][:] = dst_ubar

    ncv.variables['ocean_time'][:] = time/24.
    ncv.variables['v'][:] = dst_v
    ncv.variables['vbar'][:] = dst_vbar

    print time
    print np.shape(dst_u)
    print np.shape(dst_v)
    # close destination file
    ncu.close()
    ncv.close()
