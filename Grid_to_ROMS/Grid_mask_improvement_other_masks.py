import netCDF4
rr = netCDF4.Dataset('Parent15_Thames_deep.nc', 'a', format='NETCDF4')
mask_rho = rr.variables['mask_rho'][:]
h = rr.variables['h'][:]
for i in range(len(h)):
	for j in range(len(h.T)):
		if h[i,j]<4.0 and mask_rho[i,j]==1:
			h[i,j]=4.0

rr.variables['mask_u'][:] = mask_rho[:,1:]*mask_rho[:,:-1]
rr.variables['mask_v'][:] = mask_rho[1:,:]*mask_rho[:-1,:]
rr.variables['mask_psi'][:] = mask_rho[1:,1:]*mask_rho[:-1,1:]*mask_rho[1:,:-1]*mask_rho[:-1,:-1]

rr.close() 
