import netCDF4; from shutil import copyfile; import numpy as np; spval= -32767

copyfile('/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/Boundary.nc', '/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Sediments/Boundary_sediment.nc')
rr = netCDF4.Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Sediments/Boundary_sediment.nc', 'a', format='NETCDF4')

rr.createDimension('mud_time',len(rr.variables['temp_north'][:,0,0]))

rr.createVariable('mud_north_01', 'f8', ('mud_time','s_rho','xi_rho',), fill_value=spval)
rr.variables['mud_north_01'].long_name = 'suspended cohesive sediment northern boundary condition, size class 01'
rr.variables['mud_north_01'].units = 'kilogram meter-3'
rr.variables['mud_north_01'].time = 'mud_time'
rr.variables['mud_north_01'].field = 'mud_north_01, scalar, series'

rr.createVariable('mud_north_02', 'f8', ('mud_time','s_rho','xi_rho',), fill_value=spval)
rr.variables['mud_north_02'].long_name = 'suspended cohesive sediment northern boundary condition, size class 02'
rr.variables['mud_north_02'].units = 'kilogram meter-3'
rr.variables['mud_north_02'].time = 'mud_time'
rr.variables['mud_north_02'].field = 'mud_north_02, scalar, series'

rr.createVariable('mud_north_03', 'f8', ('mud_time','s_rho','xi_rho',), fill_value=spval)
rr.variables['mud_north_03'].long_name = 'suspended cohesive sediment northern boundary condition, size class 03'
rr.variables['mud_north_03'].units = 'kilogram meter-3'
rr.variables['mud_north_03'].time = 'mud_time'
rr.variables['mud_north_03'].field = 'mud_north_03, scalar, series'

rr.createVariable('mud_north_04', 'f8', ('mud_time','s_rho','xi_rho',), fill_value=spval)
rr.variables['mud_north_04'].long_name = 'suspended cohesive sediment northern boundary condition, size class 04'
rr.variables['mud_north_04'].units = 'kilogram meter-3'
rr.variables['mud_north_04'].time = 'mud_time'
rr.variables['mud_north_04'].field = 'mud_north_04, scalar, series'

rr.createVariable('mud_time', 'f4', ('mud_time'))
rr.variables['mud_time'].long_name = 'mud_time'
rr.variables['mud_time'].units = "days"
rr.variables['mud_time'].calendar = "gregorian"
rr.variables['mud_time'][:]=rr.variables['temp_time'][:]

rr.createVariable('mud_south_01', 'f8', ('mud_time','s_rho','xi_rho',), fill_value=spval)
rr.variables['mud_south_01'].long_name = 'suspended cohesive sediment southern boundary condition, size class 01'
rr.variables['mud_south_01'].units = 'kilogram meter-3'
rr.variables['mud_south_01'].time = 'mud_time'
rr.variables['mud_south_01'].field = 'mud_south_01, scalar, series'

rr.createVariable('mud_south_02', 'f8', ('mud_time','s_rho','xi_rho',), fill_value=spval)
rr.variables['mud_south_02'].long_name = 'suspended cohesive sediment southern boundary condition, size class 02'
rr.variables['mud_south_02'].units = 'kilogram meter-3'
rr.variables['mud_south_02'].time = 'mud_time'
rr.variables['mud_south_02'].field = 'mud_south_02, scalar, series'

rr.createVariable('mud_south_03', 'f8', ('mud_time','s_rho','xi_rho',), fill_value=spval)
rr.variables['mud_south_03'].long_name = 'suspended cohesive sediment southern boundary condition, size class 03'
rr.variables['mud_south_03'].units = 'kilogram meter-3'
rr.variables['mud_south_03'].time = 'mud_time'
rr.variables['mud_south_03'].field = 'mud_south_03, scalar, series'

rr.createVariable('mud_south_04', 'f8', ('mud_time','s_rho','xi_rho',), fill_value=spval)
rr.variables['mud_south_04'].long_name = 'suspended cohesive sediment southern boundary condition, size class 04'
rr.variables['mud_south_04'].units = 'kilogram meter-3'
rr.variables['mud_south_04'].time = 'mud_time'
rr.variables['mud_south_04'].field = 'mud_south_04, scalar, series'

rr.createVariable('mud_west_01', 'f8', ('mud_time','s_rho','eta_rho',), fill_value=spval)
rr.variables['mud_west_01'].long_name = 'suspended cohesive sediment western boundary condition, size class 01'
rr.variables['mud_west_01'].units = 'kilogram meter-3'
rr.variables['mud_west_01'].time = 'mud_time'
rr.variables['mud_west_01'].field = 'mud_west_01, scalar, series'

rr.createVariable('mud_west_02', 'f8', ('mud_time','s_rho','eta_rho',), fill_value=spval)
rr.variables['mud_west_02'].long_name = 'suspended cohesive sediment western boundary condition, size class 02'
rr.variables['mud_west_02'].units = 'kilogram meter-3'
rr.variables['mud_west_02'].time = 'mud_time'
rr.variables['mud_west_02'].field = 'mud_west_02, scalar, series'

rr.createVariable('mud_west_03', 'f8', ('mud_time','s_rho','eta_rho',), fill_value=spval)
rr.variables['mud_west_03'].long_name = 'suspended cohesive sediment western boundary condition, size class 03'
rr.variables['mud_west_03'].units = 'kilogram meter-3'
rr.variables['mud_west_03'].time = 'mud_time'
rr.variables['mud_west_03'].field = 'mud_west_03, scalar, series'

rr.createVariable('mud_west_04', 'f8', ('mud_time','s_rho','eta_rho',), fill_value=spval)
rr.variables['mud_west_04'].long_name = 'suspended cohesive sediment western boundary condition, size class 04'
rr.variables['mud_west_04'].units = 'kilogram meter-3'
rr.variables['mud_west_04'].time = 'mud_time'
rr.variables['mud_west_04'].field = 'mud_west_04, scalar, series'

rr.createVariable('mud_east_01', 'f8', ('mud_time','s_rho','eta_rho',), fill_value=spval)
rr.variables['mud_east_01'].long_name = 'suspended cohesive sediment eastern boundary condition, size class 01'
rr.variables['mud_east_01'].units = 'kilogram meter-3'
rr.variables['mud_east_01'].time = 'mud_time'
rr.variables['mud_east_01'].field = 'mud_east_01, scalar, series'

rr.createVariable('mud_east_02', 'f8', ('mud_time','s_rho','eta_rho',), fill_value=spval)
rr.variables['mud_east_02'].long_name = 'suspended cohesive sediment eastern boundary condition, size class 02'
rr.variables['mud_east_02'].units = 'kilogram meter-3'
rr.variables['mud_east_02'].time = 'mud_time'
rr.variables['mud_east_02'].field = 'mud_east_02, scalar, series'

rr.createVariable('mud_east_03', 'f8', ('mud_time','s_rho','eta_rho',), fill_value=spval)
rr.variables['mud_east_03'].long_name = 'suspended cohesive sediment eastern boundary condition, size class 03'
rr.variables['mud_east_03'].units = 'kilogram meter-3'
rr.variables['mud_east_03'].time = 'mud_time'
rr.variables['mud_east_03'].field = 'mud_east_03, scalar, series'

rr.createVariable('mud_east_04', 'f8', ('mud_time','s_rho','eta_rho',), fill_value=spval)
rr.variables['mud_east_04'].long_name = 'suspended cohesive sediment eastern boundary condition, size class 04'
rr.variables['mud_east_04'].units = 'kilogram meter-3'
rr.variables['mud_east_04'].time = 'mud_time'
rr.variables['mud_east_04'].field = 'mud_east_04, scalar, series'

rr.close()

print 'Draft is prepared'
rr = netCDF4.Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Sediments/Boundary_sediment.nc', 'a', format='NETCDF4')
var = ['mud_north_01','mud_north_02','mud_north_03','mud_north_04','mud_south_01','mud_south_02','mud_south_03','mud_south_04','mud_west_01','mud_west_02','mud_west_03','mud_west_04','mud_east_01','mud_east_02','mud_east_03','mud_east_04']
val = [0.00125,0.00125,0.00125,0.00125,0,0,0,0,0,0,0,0,0,0,0,0]

def func(varr,vall):
	t = rr.variables['temp_'+varr.split('_')[1]][0]
	aa = np.tile(np.full(t.shape,vall),(len(rr.variables['temp_north'][:,0,0]),1,1))
	msk_3d = np.zeros(aa.shape, dtype=bool)
	msk=t.mask; msk_3d[:,:] = t.mask==True
	cc = np.expand_dims(aa, axis=0)
	bb = np.ma.array(cc, mask=msk_3d)
	print varr,np.shape(rr.variables[varr][:]), np.shape(bb)
	rr.variables[varr][:] = bb
for i in range(16):
	func(var[i],val[i])
rr.close()

