import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.mlab import griddata
import pyroms
import netCDF4
from shapely.geometry import Polygon


dst_grd = pyroms.grid.get_ROMS_grid('BCZ')

figure, ax = plt.subplots(figsize=(20,10))
ax.imshow(dst_grd.hgrid.mask_rho, origin='lower', interpolation='none')

def onclick(event):
	print('button=%d, x=%d, y=%d, xdata=%d, ydata=%d' %
		(event.button, event.x, event.y, int(event.xdata), int(event.ydata)))
	if event.button==1:
		dst_grd.hgrid.mask_rho[int(event.ydata),int(event.xdata)]=0
		#dst_grd.hgrid.mask_rho[int(event.ydata)+1,int(event.xdata)]=0
		#dst_grd.hgrid.mask_rho[int(event.ydata)-1,int(event.xdata)]=0
		#dst_grd.hgrid.mask_rho[int(event.ydata),int(event.xdata)+1]=0
		#dst_grd.hgrid.mask_rho[int(event.ydata),int(event.xdata)-1]=0
		#dst_grd.hgrid.mask_rho[int(event.ydata)+1,int(event.xdata)-1]=0
		#dst_grd.hgrid.mask_rho[int(event.ydata)+1,int(event.xdata)+1]=0
		#dst_grd.hgrid.mask_rho[int(event.ydata)-1,int(event.xdata)-1]=0
		#dst_grd.hgrid.mask_rho[int(event.ydata)-1,int(event.xdata)+1]=0
		print event.ydata, event.xdata, dst_grd.hgrid.mask_rho[int(event.ydata),int(event.xdata)]
	elif event.button==3:
		dst_grd.hgrid.mask_rho[int(event.ydata),int(event.xdata)]=1
		#dst_grd.hgrid.mask_rho[int(event.ydata)+1,int(event.xdata)]=1
		#dst_grd.hgrid.mask_rho[int(event.ydata)-1,int(event.xdata)]=1
		#dst_grd.hgrid.mask_rho[int(event.ydata),int(event.xdata)+1]=1
		#dst_grd.hgrid.mask_rho[int(event.ydata),int(event.xdata)-1]=1
		#dst_grd.hgrid.mask_rho[int(event.ydata)+1,int(event.xdata)-1]=1
		#dst_grd.hgrid.mask_rho[int(event.ydata)+1,int(event.xdata)+1]=1
		#dst_grd.hgrid.mask_rho[int(event.ydata)-1,int(event.xdata)-1]=1
		#dst_grd.hgrid.mask_rho[int(event.ydata)-1,int(event.xdata)+1]=1
	else: 
		pyroms.grid.write_ROMS_grid(dst_grd, filename='Grid_improved_mask.nc')
	ax.imshow(dst_grd.hgrid.mask_rho, origin='lower', interpolation='none')
	figure.canvas.draw()

cid = figure.canvas.mpl_connect('button_press_event', onclick)
plt.show()
