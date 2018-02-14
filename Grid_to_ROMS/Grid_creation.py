import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.mlab import griddata
import scipy
import pyroms
import pyroms_toolbox
from bathy_smoother import *
import netCDF4
from shapely.geometry import Polygon


Lp=250
Mp=250


rr = netCDF4.Dataset('GEBCO_2.nc', 'r', format='NETCDF4')
lats = rr.variables['lat'][:]
lons = rr.variables['lon'][:]
topo = rr.variables['elevation'][:,:]*(-1)
'''
Lp=70
Mp=52
lon1=1.000 ; lat1=49.740
lon2=5.900 ; lat2=51.300
lon3=4.000 ; lat3=53.500
lon0=-0.750 ; lat0=52.050
-4.99033200764 -16.8064107162 9.92032801679 -0.057273881369 46.6325836395 55.1210682254 52.3870158488 62.3190665216
'''
'''
lon0=-16.8064 ; lat0=55.121
lon2=9.9203 ; lat2=52.387
lon1=-0.0573 ; lat1=62.319
lon3=-4.9903 ; lat3=46.6325
'''
Lp=250
Mp=250
lon0=-16.5 ; lat0=55.9
lon1=-5.05 ; lat1=46.2
lon2=10 ; lat2=53.65
lon3=-1.35 ; lat3=61.7
lonp = np.array([lon0, lon1, lon2, lon3])
latp = np.array([lat0, lat1, lat2, lat3])
beta = np.array([1, 1, 1, 1])

#fig, ax = plt.subplots(figsize=(20,10))
map = Basemap(projection='merc', llcrnrlon=-17, llcrnrlat=46, urcrnrlon=11, urcrnrlat=62, lat_ts=52, resolution='h')
map.drawcoastlines()
map.drawcountries()


hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mp+3,Lp+3), proj=map)

#hgrd.x_vert= hgrd.x_vert[210:560,0:270]
#hgrd.y_vert= hgrd.y_vert[210:560,0:270]


lonv, latv = map(hgrd.x_vert, hgrd.y_vert, inverse=True)
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)

map.plot(hgrd.x, hgrd.y, '-k')
map.plot(hgrd.x.T, hgrd.y.T, '-k')
plt.show()

def point_inside_polygon(x,y,poly):	# Ex.: a=[(0,0),(0,2),(1,1),(1,0)] if x=0.4 y=1.2 function gives "True" if x=0.8 y=1.8 function gives "False"
	n = len(poly)
	inside =False
	p1x,p1y = poly[0]
	for i in range(n+1):
		p2x,p2y = poly[i % n]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xinters:
						inside = not inside
		p1x,p1y = p2x,p2y
	return inside

def PolygonArea(corners):		# Shoelace formula; Ex.: a=[(0,0),(0,2),(1,1),(1,0)] Function gives 1.5
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

def Area(vertices):
    n = len(vertices) # of corners
    a = 0.0
    for i in range(n):
        j = (i + 1) % n
        a += abs(vertices[i][0] * vertices[j][1]-vertices[j][0] * vertices[i][1])
    result = a / 2.0
    return result

def Angle(corners):			
	lb=(((corners[1][0]-corners[0][0])**2)+((corners[1][1]-corners[0][1])**2))**0.5
	ub=(((corners[2][0]-corners[1][0])**2)+((corners[2][1]-corners[1][1])**2))**0.5
	cb=(((corners[2][0]-corners[0][0])**2)+((corners[2][1]-corners[0][1])**2))**0.5
	frac_up=cb*cb-lb*lb-ub*ub
	frac_down=2*lb*ub
	alpha1=np.math.acos(-frac_up/frac_down)
	alpha1=np.math.fabs(90-np.math.degrees(alpha1))
	lb=(((corners[2][0]-corners[1][0])**2)+((corners[2][1]-corners[1][1])**2))**0.5
	ub=(((corners[3][0]-corners[2][0])**2)+((corners[3][1]-corners[2][1])**2))**0.5
	cb=(((corners[3][0]-corners[1][0])**2)+((corners[3][1]-corners[1][1])**2))**0.5
	frac_up=cb*cb-lb*lb-ub*ub
	frac_down=2*lb*ub
	alpha2=np.math.acos(-frac_up/frac_down)
	alpha2=np.math.fabs(90-np.math.degrees(alpha2))
	lb=(((corners[3][0]-corners[2][0])**2)+((corners[3][1]-corners[2][1])**2))**0.5
	ub=(((corners[0][0]-corners[3][0])**2)+((corners[0][1]-corners[3][1])**2))**0.5
	cb=(((corners[0][0]-corners[2][0])**2)+((corners[0][1]-corners[2][1])**2))**0.5
	frac_up=cb*cb-lb*lb-ub*ub
	frac_down=2*lb*ub
	alpha3=np.math.acos(-frac_up/frac_down)
	alpha3=np.math.fabs(90-np.math.degrees(alpha3))
	lb=(((corners[0][0]-corners[3][0])**2)+((corners[0][1]-corners[3][1])**2))**0.5
	ub=(((corners[1][0]-corners[0][0])**2)+((corners[1][1]-corners[0][1])**2))**0.5
	cb=(((corners[1][0]-corners[3][0])**2)+((corners[1][1]-corners[3][1])**2))**0.5
	frac_up=cb*cb-lb*lb-ub*ub
	frac_down=2*lb*ub
	alpha4=np.math.acos(-frac_up/frac_down)
	alpha4=np.math.fabs(90-np.math.degrees(alpha4))
	alpha=(alpha1+alpha2+alpha3+alpha4)/4.
	return alpha

def Width(corners):			
	ls=[((corners[0][0]+corners[1][0])/2.,(corners[0][1]+corners[1][1])/2.)]
	rs=[((corners[2][0]+corners[3][0])/2.,(corners[2][1]+corners[3][1])/2.)]
	width=(((rs[0][0]-ls[0][0])**2)+((rs[0][1]-ls[0][1])**2))**0.5
	return width

def Height(corners):
	us=[((corners[1][0]+corners[2][0])/2.,(corners[1][1]+corners[2][1])/2.)]
	ds=[((corners[0][0]+corners[3][0])/2.,(corners[0][1]+corners[3][1])/2.)]
	height=(((us[0][0]-ds[0][0])**2)+((us[0][1]-ds[0][1])**2))**0.5
	return height

def Aspect(width,height):			# Calculation of aspect; Ex.: a=[(0,0),(0,2),(1,1),(1,0)] Function gives 0.7453559924999299
	aspect=width/height
	return aspect
	
	
count_p=np.zeros((len(hgrd.lat_vert)-1,len(hgrd.lon_vert.T)-1))
depth=np.zeros((len(hgrd.lat_vert)-1,len(hgrd.lon_vert.T)-1))
aspect=np.zeros((len(hgrd.lat_vert)-1,len(hgrd.lon_vert.T)-1))
angle=np.zeros((len(hgrd.lat_vert)-1,len(hgrd.lon_vert.T)-1))
area=np.zeros((len(hgrd.lat_vert)-1,len(hgrd.lon_vert.T)-1))
area_sq=np.zeros((len(hgrd.lat_vert)-1,len(hgrd.lon_vert.T)-1))
width1=np.zeros((len(hgrd.lat_vert)-1,len(hgrd.lon_vert.T)-1))
height1=np.zeros((len(hgrd.lat_vert)-1,len(hgrd.lon_vert.T)-1))

for i in range(len(hgrd.lat_vert)-1):
	for j in range(len(hgrd.lon_vert.T)-1):
		poly=[(hgrd.lat_vert[i,j],hgrd.lon_vert[i,j]),(hgrd.lat_vert[i,j+1],hgrd.lon_vert[i,j+1]),(hgrd.lat_vert[i+1,j+1],hgrd.lon_vert[i+1,j+1]),(hgrd.lat_vert[i+1,j],hgrd.lon_vert[i+1,j])]
		poly_cartesian=[(hgrd.y_vert[i,j],hgrd.x_vert[i,j]),(hgrd.y_vert[i,j+1],hgrd.x_vert[i,j+1]),(hgrd.y_vert[i+1,j+1],hgrd.x_vert[i+1,j+1]),(hgrd.y_vert[i+1,j],hgrd.x_vert[i+1,j])]
		for k in range(np.where(lats>min(np.array(poly)[:,0]))[0][0],np.where(lats<max(np.array(poly)[:,0]))[0][-1]):
			for l in range(np.where(lons>min(np.array(poly)[:,1]))[0][0],np.where(lons<max(np.array(poly)[:,1]))[0][-1]):
				if point_inside_polygon(lats[k],lons[l],poly)==True:
					count_p[i,j]=count_p[i,j]+1
					depth[i,j]=topo[k,l]+depth[i,j]
				else:
					pass
		width1[i,j]=Width(poly_cartesian)/1000
		height1[i,j]=Height(poly_cartesian)/1000
		aspect[i,j]=round((width1[i,j]/height1[i,j]),3)
		#angle[i,j]=round((Angle(poly_cartesian)),2)
		#area[i,j]=round(Area(poly_cartesian)/1000000,0)

def in_circle(center_x, center_y, radius, x, y):		#	verification if the topo point is inside a sphere
    square_dist = (center_x - x) ** 2 + (center_y - y) ** 2
    return square_dist <= radius ** 2

for i in range(len(depth)):
	for j in range(len(depth.T)):
		if count_p[i,j]>0:
			depth[i,j]=depth[i,j]/count_p[i,j]
		else:
			pass


for k in range(len(hgrd.lat_rho)):
	for l in range(len(hgrd.lon_rho.T)):
		if count_p[k,l]==0:
			lon_int=[]
			lat_int=[]
			index_i=[]
			index_j=[]
			for i in range(np.where(lats<hgrd.lat_rho[k,l]-0.01)[0][-1],np.where(lats>hgrd.lat_rho[k,l]+0.01)[0][0]):		
				for j in range(np.where(lons<hgrd.lon_rho[k,l]-0.01)[0][-1],np.where(lons>hgrd.lon_rho[k,l]+0.01)[0][0]):
					if in_circle(hgrd.lon_rho[k,l], hgrd.lat_rho[k,l], 1, lons[j],lats[i])==True:
						lat_int.append(lats[i])
						lon_int.append(lons[j])
						index_i.append(i)
						index_j.append(j)
			topo_int=np.zeros((len(lon_int),len(lat_int)))
			for s in range(len(index_i)):
				for t in range(len(index_j)):
					topo_int[s,t]=topo[index_i[s],index_j[t]]
			b=np.ravel(topo_int)
			lon, lat = np.meshgrid(lon_int, lat_int)
			lon_lat=np.zeros((len(lon_int)*len(lat_int),2))
			lon_lat[:,0]=np.ravel(lon)
			lon_lat[:,1]=np.ravel(lat)
			depth[k,l] = round(scipy.interpolate.griddata(lon_lat,b,(hgrd.lon_rho[k,l],hgrd.lat_rho[k,l]),method='linear'),1)
		else:
			#print 'pass: it has depth points'
			pass
		if depth[k,l]<1.0:
			hgrd.mask_rho[k,l]=0				# Masking of land points and depths between 0 and 1.5 meters
		#print k,l



figure, ax = plt.subplots(figsize=(20,10))
ax.imshow(hgrd.mask_rho, origin='lower', interpolation='none')


print 'fix minimum depth'
hmin = 2
#topo = pyroms_toolbox.change(topo, '<', hmin, hmin)

print 'insure that depth is always deeper than hmin'
depth = pyroms_toolbox.change(depth, '<', hmin, hmin)

print 'set depth to hmin where masked'
idx = np.where(hgrd.mask_rho == 0)
depth[idx] = hmin

print 'save raw bathymetry'
hraw = depth.copy()

print 'check bathymetry roughness'
RoughMat = bathy_tools.RoughnessMatrix(depth, hgrd.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()

print 'smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)'
rx0_max = 0.35
depth = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, depth, rx0_max)

print 'check bathymetry roughness again'
RoughMat = bathy_tools.RoughnessMatrix(depth, hgrd.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()

print 'vertical coordinate'
theta_b = 3
theta_s = 2
Tcline=10
N = 15
vgrd = pyroms.vgrid.s_coordinate_4(depth, theta_b, theta_s, Tcline, N, hraw=hraw)

print "improving angles"
for i in range(len(hgrd.angle_rho)):
	for j in range(len(hgrd.angle_rho.T)):
		if hgrd.angle_rho[i,j]>1.55 and hgrd.angle_rho[i,j]<1.59:
			hgrd.angle_rho[i,j]=-1.55
			print "lalala"


print 'ROMS grid'
grd_name = 'LOL'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

print 'write grid to netcdf file'
pyroms.grid.write_ROMS_grid(grd, filename='Coarser_grid.nc')
