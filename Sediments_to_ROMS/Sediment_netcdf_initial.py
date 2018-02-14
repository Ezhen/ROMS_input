import netCDF4; from shutil import copyfile; import numpy as np; spval= -32767

copyfile('/media/sf_Swap-between-windows-linux/New_Grid/Interannual_simulation/Initial.nc', '/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Sediments/Initial_sediment.nc')
rr = netCDF4.Dataset('/home/eivanov/coawst_data_prrocessing/GRID_RECTANGULAR/Sediments/Initial_sediment.nc', 'a', format='NETCDF4')
msk = rr.variables['mask_rho'][:]; h = rr.variables['h'][:]

rr.createDimension('Nbed', 2)

rr.createVariable('mud_01', 'f8', ('ocean_time','s_rho','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mud_01'].long_name = 'suspended cohesive sediment, size class 01'
rr.variables['mud_01'].units = 'kilogram meter-3'
rr.variables['mud_01'].time = 'ocean_time'
rr.variables['mud_01'].field = 'mud_01, scalar, series'

rr.createVariable('mudfrac_01','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mudfrac_01'].long_name = 'cohesive sediment fraction, size class 01'
rr.variables['mudfrac_01'].units = 'nondimensional'
rr.variables['mudfrac_01'].time = 'ocean_time'
rr.variables['mudfrac_01'].field = 'mudfrac_01, scalar, series'

rr.createVariable('mudmass_01','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mudmass_01'].long_name = 'cohesive sediment mass, size class 01'
rr.variables['mudmass_01'].units = 'kilogram meter-2'
rr.variables['mudmass_01'].time = 'ocean_time'
rr.variables['mudmass_01'].field = 'mudfrac_01, scalar, series'

rr.createVariable('mud_02', 'f8', ('ocean_time','s_rho','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mud_02'].long_name = 'suspended cohesive sediment, size class 02'
rr.variables['mud_02'].units = 'kilogram meter-3'
rr.variables['mud_02'].time = 'ocean_time'
rr.variables['mud_02'].field = 'mud_02, scalar, series'

rr.createVariable('mudfrac_02','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mudfrac_02'].long_name = 'cohesive sediment fraction, size class 02'
rr.variables['mudfrac_02'].units = 'nondimensional'
rr.variables['mudfrac_02'].time = 'ocean_time'
rr.variables['mudfrac_02'].field = 'mudfrac_02, scalar, series'

rr.createVariable('mudmass_02','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mudmass_02'].long_name = 'cohesive sediment mass, size class 02'
rr.variables['mudmass_02'].units = 'kilogram meter-2'
rr.variables['mudmass_02'].time = 'ocean_time'
rr.variables['mudmass_02'].field = 'mudfrac_02, scalar, series'

rr.createVariable('mud_03', 'f8', ('ocean_time','s_rho','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mud_03'].long_name = 'suspended cohesive sediment, size class 03'
rr.variables['mud_03'].units = 'kilogram meter-3'
rr.variables['mud_03'].time = 'ocean_time'
rr.variables['mud_03'].field = 'mud_03, scalar, series'

rr.createVariable('mudfrac_03','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mudfrac_03'].long_name = 'cohesive sediment fraction, size class 03'
rr.variables['mudfrac_03'].units = 'nondimensional'
rr.variables['mudfrac_03'].time = 'ocean_time'
rr.variables['mudfrac_03'].field = 'mudfrac_03, scalar, series'

rr.createVariable('mudmass_03','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mudmass_03'].long_name = 'cohesive sediment mass, size class 03'
rr.variables['mudmass_03'].units = 'kilogram meter-2'
rr.variables['mudmass_03'].time = 'ocean_time'
rr.variables['mudmass_03'].field = 'mudfrac_03, scalar, series'

rr.createVariable('mud_04', 'f8', ('ocean_time','s_rho','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mud_04'].long_name = 'suspended cohesive sediment, size class 04'
rr.variables['mud_04'].units = 'kilogram meter-3'
rr.variables['mud_04'].time = 'ocean_time'
rr.variables['mud_04'].field = 'mud_04, scalar, series'

rr.createVariable('mudfrac_04','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mudfrac_04'].long_name = 'cohesive sediment fraction, size class 04'
rr.variables['mudfrac_04'].units = 'nondimensional'
rr.variables['mudfrac_04'].time = 'ocean_time'
rr.variables['mudfrac_04'].field = 'mudfrac_04, scalar, series'

rr.createVariable('mudmass_04','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['mudmass_04'].long_name = 'cohesive sediment mass, size class 04'
rr.variables['mudmass_04'].units = 'kilogram meter-2'
rr.variables['mudmass_04'].time = 'ocean_time'
rr.variables['mudmass_04'].field = 'mudfrac_04, scalar, series'

rr.createVariable('bed_thickness','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['bed_thickness'].long_name = 'sediment layer thickness'
rr.variables['bed_thickness'].units = 'meter'
rr.variables['bed_thickness'].time = 'ocean_time'
rr.variables['bed_thickness'].field = 'bed thickness, scalar, series'

rr.createVariable('bed_age','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['bed_age'].long_name = 'sediment layer age'
rr.variables['bed_age'].units = 'day'
rr.variables['bed_age'].time = 'ocean_time'
rr.variables['bed_age'].field = 'bed age, scalar, series'

rr.createVariable('bed_porosity','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['bed_porosity'].long_name = 'sediment layer porosity'
rr.variables['bed_porosity'].units = 'nondimensional'
rr.variables['bed_porosity'].time = 'ocean_time'
rr.variables['bed_porosity'].field = 'porosity, scalar, series'

rr.createVariable('bed_tau_crit','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['bed_tau_crit'].long_name = 'tau_critical in each layer'
rr.variables['bed_tau_crit'].units = 'Newtons meter-2'
rr.variables['bed_tau_crit'].time = 'ocean_time'
rr.variables['bed_tau_crit'].field = 'bed_tau_crit, scalar, series'

rr.createVariable('bed_biodiff','f8',('ocean_time', 'Nbed','eta_rho','xi_rho',), fill_value=spval)
rr.variables['bed_biodiff'].long_name = 'biodiffusivity at bottom of each layer'
rr.variables['bed_biodiff'].units = 'meter2 second-1'
rr.variables['bed_biodiff'].time = 'ocean_time'
rr.variables['bed_biodiff'].field = 'bed biodiffusivity, scalar, series'

rr.createVariable('grain_diameter','f8',('ocean_time', 'eta_rho','xi_rho',), fill_value=spval)
rr.variables['grain_diameter'].long_name = 'sediment median grain diameter size'
rr.variables['grain_diameter'].units = 'meter'
rr.variables['grain_diameter'].time = 'ocean_time'
rr.variables['grain_diameter'].field = 'grain diameter, scalar, series'

rr.createVariable('grain_density','f8',('ocean_time', 'eta_rho','xi_rho',), fill_value=spval)
rr.variables['grain_density'].long_name = 'sediment median grain density'
rr.variables['grain_density'].units = 'kilogram meter-3'
rr.variables['grain_density'].time = 'ocean_time'
rr.variables['grain_density'].field = 'grain density, scalar, series'

rr.createVariable('settling_vel','f8',('ocean_time', 'eta_rho','xi_rho',), fill_value=spval)
rr.variables['settling_vel'].long_name = 'sediment median grain settling velocity'
rr.variables['settling_vel'].units = 'meter second-1'
rr.variables['settling_vel'].time = 'ocean_time'
rr.variables['settling_vel'].field = 'settling vel, scalar, series'

rr.createVariable('erosion_stress','f8',('ocean_time', 'eta_rho','xi_rho',), fill_value=spval)
rr.variables['erosion_stress'].long_name = 'sediment median critical erosion stress'
rr.variables['erosion_stress'].units = 'meter2 second-2'
rr.variables['erosion_stress'].time = 'ocean_time'
rr.variables['erosion_stress'].field = 'erosion stress, scalar, series'

var = ['mud_01','mudfrac_01','mudmass_01','mud_02','mudfrac_02','mudmass_02','mud_03','mudfrac_03','mudmass_03','mud_04','mudfrac_04','mudmass_04','bed_thickness','bed_age','bed_porosity','bed_tau_crit','bed_biodiff','grain_diameter','grain_density','settling_vel','erosion_stress']
#val = [0.00125,0.25,0,0.00125,0.25,0,0.00125,0.25,0,0.00125,0.25,0,2,0,0,2,0,0.00075,0.75,0.75,0.75]
val = [0,0.25,0,0,0.25,0,0,0.25,0,0.00125,0.25,0,2,0,0,2,0,0.00075,0.75,0.75,0.75]
num = [15,2,2,15,2,2,15,2,2,15,2,2,2,2,2,2,2,0,0,0,0]

def func(varr,vall,numm):
	if numm>0:
		aa = np.tile(np.full(msk.shape,vall),(numm,1,1))
	elif numm==0:
		aa = np.full(msk.shape,vall)
	if varr=='bed_tau_crit':
		aa[0] = 0.5; aa[1] = 2
	if varr=='bed_thickness':
		aa[0] = 0.01; aa[1] = 0.01
	msk_3d = np.zeros(aa.shape, dtype=bool)
	msk_3d[:,:] = msk[:,:] == 0
	cc = np.expand_dims(aa, axis=0)
	bb = np.ma.array(cc, mask=msk_3d)
	print varr,np.shape(rr.variables[varr][:]), np.shape(bb), vall
	rr.variables[varr][:] = bb
for i in range(21):
	func(var[i],val[i],num[i])
rr.close()
