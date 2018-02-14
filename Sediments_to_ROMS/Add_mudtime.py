import netCDF4

rr = netCDF4.Dataset('Initial_sediment.nc', 'r', format='NETCDF4')
s_rho = rr.variables['s_rho'][:]
s_w = rr.variables['s_w'][:]
Cs_r = rr.variables['Cs_r'][:]
Cs_w = rr.variables['Cs_w'][:]
theta_s = rr.variables['theta_s'][:]
theta_b = rr.variables['theta_b'][:]
rr.close()

rr = netCDF4.Dataset('Boundary_sediment.nc', 'a', format='NETCDF4')
rr.variables['s_rho'][:] = s_rho
rr.variables['s_w'][:] = s_w
rr.variables['Cs_r'][:] = Cs_r
rr.variables['Cs_w'][:] = Cs_w
rr.variables['theta_s'][:] = theta_s
rr.variables['theta_b'][:] = theta_b
rr.close()
