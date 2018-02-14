import numpy as np
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF


def nc_create_roms_bdry_file_mine(filename, grd):

    # create file
    nc = netCDF.Dataset(filename, 'w', format='NETCDF3_64BIT')
    nc.Description = 'ROMS file'
    nc.Author = 'pyroms_toolbox.nc_create_roms_file'
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'ROMS file'

    nc.createDimension('xr', np.size(grd.hgrid.mask_rho,1))
    nc.createDimension('xu', np.size(grd.hgrid.mask_u,1))
    nc.createDimension('xv', np.size(grd.hgrid.mask_v,1))
    nc.createDimension('xpsi', np.size(grd.hgrid.mask_psi,1))
    nc.createDimension('er', np.size(grd.hgrid.mask_rho,0))
    nc.createDimension('eu', np.size(grd.hgrid.mask_u,0))
    nc.createDimension('ev', np.size(grd.hgrid.mask_v,0))
    nc.createDimension('epsi', np.size(grd.hgrid.mask_psi,0))
    nc.createDimension('sc_r', grd.vgrid.N)
    nc.createDimension('sc_w', grd.vgrid.Np)

    # write time and grid information
    nc.createVariable('theta_s', 'f8', ())
    nc.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
    nc.variables['theta_s'][:] = grd.vgrid.theta_s

    nc.createVariable('theta_b', 'f8', ())
    nc.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
    nc.variables['theta_b'][:] = grd.vgrid.theta_b

    nc.createVariable('Tcline', 'f8', ())
    nc.variables['Tcline'].long_name = 'S-cordinate surface/bottom layer width'
    nc.variables['Tcline'].units = 'meter'
    nc.variables['Tcline'][:] = grd.vgrid.Tcline

    nc.createVariable('hc', 'f8', ())
    nc.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
    nc.variables['hc'].units = 'meter'
    nc.variables['hc'][:] = grd.vgrid.hc

    nc.createVariable('sc_r', 'f8', ('sc_r'))
    nc.variables['sc_r'].long_name = 'S-coordinate at RHO-points'
    nc.variables['sc_r'].valid_min = '-1'
    nc.variables['sc_r'].valid_max = '0'
    nc.variables['sc_r'].field = 'sc_r,scalar'
    nc.variables['sc_r'][:] = grd.vgrid.s_rho

    nc.createVariable('sc_w', 'f8', ('sc_w'))
    nc.variables['sc_w'].long_name = 'S-coordinate at W-points'
    nc.variables['sc_w'].valid_min = '-1'
    nc.variables['sc_w'].valid_max = '0'
    nc.variables['sc_w'].field = 'sc_w,scalar'
    nc.variables['sc_w'][:] = grd.vgrid.s_w

    nc.createVariable('Cs_r', 'f8', ('sc_r'))
    nc.variables['Cs_r'].long_name = 'S-coordinate stretching curves at RHO-points'
    nc.variables['Cs_r'].valid_min = '-1'
    nc.variables['Cs_r'].valid_max = '0'
    nc.variables['Cs_r'].field = 'Cs_r,scalar'
    nc.variables['Cs_r'][:] = grd.vgrid.Cs_r

    nc.createVariable('Cs_w', 'f8', ('sc_w'))
    nc.variables['Cs_w'].long_name = 'S-coordinate stretching curves at W-points'
    nc.variables['Cs_w'].valid_min = '-1'
    nc.variables['Cs_w'].valid_max = '0'
    nc.variables['Cs_w'].field = 'Cs_w,scalar'
    nc.variables['Cs_w'][:] = grd.vgrid.Cs_w

    nc.createVariable('h', 'f8', ('er', 'xr'))
    nc.variables['h'].long_name = 'bathymetry at RHO-points'
    nc.variables['h'].units ='meter'
    nc.variables['h'].coordinates = 'lon_rho lat_rho'
    nc.variables['h'].field = 'bath, scalar'
    nc.variables['h'][:] = grd.vgrid.h


    nc.close()
