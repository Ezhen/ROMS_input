import numpy as np
from datetime import datetime
from netCDF4 import Dataset
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import scipy
import pyroms
import pyroms_toolbox
from shutil import copyfile

grd = pyroms.grid.get_ROMS_grid("GRIDRECT")

nud=np.zeros((len(grd.hgrid.mask_rho),len(grd.hgrid.mask_rho.T)))
nud_rho=np.zeros((15,len(grd.hgrid.mask_rho),len(grd.hgrid.mask_rho.T)))

for i in range(len(nud)):
	for j in range(len(nud.T)):
		if grd.hgrid.mask_rho[i,j]==1:
			a=min(len(nud)-i,len(nud.T)-j,i,j)
			if a<10:
					if i<50:
						nud[i,j]=1-a*0.1
					else:
						nud[i,j]=10-a

nud_rho=np.tile(nud,(15,1))
def nc_create_roms_file111(filename):
	pyroms.grid.write_ROMS_grid(grd, filename)			# NetCDF file creation
	nc = netCDF.Dataset(filename, 'a', format='NETCDF3_64BIT')
	nc.Description = 'ROMS file'
	nc.Author = 'Evgeny Ivanov'
	nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	nc.title = 'ROMS Nudging File'

	nc.createVariable('M2_NudgeCoef', 'f', ('eta_rho', 'xi_rho'))
	nc.variables['M2_NudgeCoef'].long_name = "2D momentum inverse nudging coefficients"
	nc.variables['M2_NudgeCoef'].units = "day-1"
	nc.variables['M2_NudgeCoef'].coordinates = "xi_rho eta_rho"
	nc.variables['M2_NudgeCoef'][:] = nud

	nc.createVariable('M3_NudgeCoef', 'f', ('s_rho', 'eta_rho', 'xi_rho'))
	nc.variables['M3_NudgeCoef'].long_name = "3D momentum inverse nudging coefficients"
	nc.variables['M3_NudgeCoef'].units = "day-1"
	nc.variables['M3_NudgeCoef'].coordinates = "s_rho xi_rho eta_rho"
	nc.variables['M3_NudgeCoef'][:] = nud_rho
	
	nc.createVariable('temp_NudgeCoef', 'f', ('s_rho', 'eta_rho', 'xi_rho'))
	nc.variables['temp_NudgeCoef'].long_name = "temp inverse nudging coefficients"
	nc.variables['temp_NudgeCoef'].units = "day-1"
	nc.variables['temp_NudgeCoef'].coordinates = "s_rho xi_rho eta_rho"
	nc.variables['temp_NudgeCoef'][:] = nud_rho		

	nc.createVariable('salt_NudgeCoef', 'f', ('s_rho', 'eta_rho', 'xi_rho'))
	nc.variables['salt_NudgeCoef'].long_name = "salt inverse nudging coefficients"
	nc.variables['salt_NudgeCoef'].units = "day-1"
	nc.variables['salt_NudgeCoef'].coordinates = "s_rho xi_rho eta_rho"
	nc.variables['salt_NudgeCoef'][:] = nud_rho							

	nc.close()

nc_create_roms_file111('Nudging.nc')
#print 'file is being copied'
#copyfile('Nudging.nc', '/home/eivanov/COAWST/Data/ROMS/Nudging/Nudging.nc')
