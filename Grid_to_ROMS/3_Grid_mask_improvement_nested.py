import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from netCDF4 import Dataset
from shutil import copyfile

copyfile('/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Child_river_boundary_update_bath_smoothed_boundary.nc', '/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Child_river_boundary_update_bath_smoothed_boundary_mask_updated.nc')
rr = Dataset('/home/eivanov/coawst_data_prrocessing/Temporal/Grid_Creation/Child_river_boundary_update_bath_smoothed_boundary_mask_updated.nc', 'a', format='NETCDF4')

mask_rho = rr.variables['mask_rho'][:]

figure, ax = plt.subplots(figsize=(20,10))
ax.imshow(mask_rho, origin='lower', interpolation='none')

def onclick(event):
	print('button=%d, x=%d, y=%d, xdata=%d, ydata=%d' %
		(event.button, event.x, event.y, int(event.xdata), int(event.ydata)))
	if event.button==1:
		mask_rho[int(event.ydata),int(event.xdata)]=0
		mask_rho[int(event.ydata)+1,int(event.xdata)]=0
		mask_rho[int(event.ydata)-1,int(event.xdata)]=0
		mask_rho[int(event.ydata),int(event.xdata)+1]=0
		mask_rho[int(event.ydata),int(event.xdata)-1]=0
		mask_rho[int(event.ydata)+1,int(event.xdata)-1]=0
		mask_rho[int(event.ydata)+1,int(event.xdata)+1]=0
		mask_rho[int(event.ydata)-1,int(event.xdata)-1]=0
		mask_rho[int(event.ydata)-1,int(event.xdata)+1]=0
		print event.ydata, event.xdata, mask_rho[int(event.ydata),int(event.xdata)]
	elif event.button==3:
		mask_rho[int(event.ydata),int(event.xdata)]=1
		mask_rho[int(event.ydata)+1,int(event.xdata)]=1
		mask_rho[int(event.ydata)-1,int(event.xdata)]=1
		mask_rho[int(event.ydata),int(event.xdata)+1]=1
		mask_rho[int(event.ydata),int(event.xdata)-1]=1
		mask_rho[int(event.ydata)+1,int(event.xdata)-1]=1
		mask_rho[int(event.ydata)+1,int(event.xdata)+1]=1
		mask_rho[int(event.ydata)-1,int(event.xdata)-1]=1
		mask_rho[int(event.ydata)-1,int(event.xdata)+1]=1
	else: 
		print('Finishing the job')
		rr.variables['mask_rho'][:] = mask_rho
		rr.variables['mask_u'][:] = mask_rho[:,1:]*mask_rho[:,:-1]
		rr.variables['mask_v'][:] = mask_rho[1:,:]*mask_rho[:-1,:]
		rr.variables['mask_psi'][:] = mask_rho[1:,1:]*mask_rho[:-1,1:]*mask_rho[1:,:-1]*mask_rho[:-1,:-1]
		rr.close()
	ax.imshow(mask_rho, origin='lower', interpolation='none')
	figure.canvas.draw()

cid = figure.canvas.mpl_connect('button_press_event', onclick)
plt.show()
