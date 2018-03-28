from netCDF4 import Dataset
import numpy as np
from shutil import copyfile
import pyroms

YY=[]; MM=[]; DD=[]; QQ=[]
for line in open('/home/eivanov/coawst_data_prrocessing/ROMS_input/Rivers_to_ROMS/Scheldt_Anwers.txt','r').readlines():
	YY.append(float(line[0:4]))
	MM.append(float(line[5:7]))
	DD.append(float(line[8:10]))
	QQ.append(float(line[11:16]))

M=MM[365:730]; D=DD[365:730]; Q=np.zeros((365)); n=0

for i in range(len(YY)):
	if MM[i]==2 and DD[i]==29:
		pass
	else:
		if MM[i]==12 and DD[i]==31:
			Q[n]=Q[n]+QQ[i]
			n=0
		else:
			Q[n]=Q[n]+QQ[i]
			n=n+1
		print YY[i],"-",MM[i],"-",DD[i], M[n-1], D[n-1], "\n"
for k in range(len(M)):
	Q[k]=np.round((Q[k]/31.),1)
Scheldt=Q

Scheldt_temp=[]
for line in open('/home/eivanov/coawst_data_prrocessing/ROMS_input/Rivers_to_ROMS/Temperature_all.txt','r').readlines():
	Scheldt_temp.append(float(line.split("\t")[2]))


Scheldt_3 = np.concatenate((Scheldt,Scheldt,Scheldt),axis=0)
Scheldt_temp_3 = np.concatenate((Scheldt_temp,Scheldt_temp,Scheldt_temp),axis=0)

dis, temp, salt = np.zeros((365*3,15)),np.zeros((365*3,15)),np.zeros((365*3,15))
for m in range(365*3):
	for n in range(15): 
		dis[m,n]=Scheldt_3[m]
		temp[m,n]=Scheldt_temp_3[m]
		salt[m,n]=1.0


d = pyroms.grid.get_ROMS_grid('ParentRiver')
a = []
for i in range(15):
	a.append((d.vgrid.Cs_w[i+1]-d.vgrid.Cs_w[i])*6*5000)

b=[]; b0=1
for i in range(15):
	b.append(b0*np.exp(1+i/7.))

b=b/np.sum(b)

u = np.zeros((365*3,15))
for m in range(365*3):
	for n in range(15): 
		u[m,n] = np.sum(Scheldt_3[m]*b[n])/a[n]


#for m in range(365*3):
#	print np.sum(u[m]), Scheldt_3[m]


copyfile('/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Boundary_ParentRiver_100_days.nc', '/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Boundary_ParentRiver_100_days_western_boundary.nc')


nc = Dataset('/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Boundary_ParentRiver_100_days_western_boundary.nc', 'a', format='NETCDF4')

for i in range(100):
	nc.variables['u_west'][i,:,28] = u[i]
nc.close()


nc = Dataset('/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Boundary_ParentRiver_100_days_western_boundary.nc', 'a', format='NETCDF4')

dst_grd = pyroms.grid.get_ROMS_grid('ParentRiver')
z_u_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:,0] + dst_grd.vgrid.z_w[0,:,:,1])

for i in range(100):
	nc.variables['ubar_west'][i,28] = (nc.variables['u_west'][i,:,28]  * np.diff(z_u_west[:,28])).sum() / -z_u_west[0,28]

nc.close()
