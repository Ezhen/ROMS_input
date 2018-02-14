import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from netCDF4 import Dataset
from shutil import copyfile

N=15; Nr=2

nc = Dataset('Rivers_nested.nc', 'w', format='NETCDF3_64BIT')
nc.Description = 'Discharges of major rivers'
nc.Author = 'Evgeny Ivanov'
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'Discharges of major rivers'
nc.createDimension('river_time', 365*10)					###
####nc.createDimension('river_time', 365)
nc.createDimension('xi_rho', 137)
nc.createDimension('xi_u', 136)
nc.createDimension('xi_v', 137)
nc.createDimension('eta_rho', 157)
nc.createDimension('eta_u', 157)
nc.createDimension('eta_v', 156)
nc.createDimension('s_rho', N)
nc.createDimension('river', Nr)


YY=[]; MM=[]; DD=[]; QQ=[]

nc.createVariable('river_time', 'f', ('river_time'))
nc.variables['river_time'].long_name = "river discharge time"
nc.variables['river_time'].units = "days since 2004-01-01 00:00:00"
nc.variables['river_time'][:]=np.linspace(0,3649,3650)	###

nc.createVariable('river', 'i4', ('river'))
nc.variables['river'].long_name = "river runoff identification number: 1)Thames, 2)Scheldt"
nc.variables['river'][:]=[1,2]

nc.createVariable('river_Xposition', 'f', ('river'))
nc.variables['river_Xposition'].long_name = "river XI-position at RHO-points"
nc.variables['river_Xposition'].valid_min = 1
nc.variables['river_Xposition'].valid_max = 82
nc.variables['river_Xposition'][:]=[121,6]#[120,2]

nc.createVariable('river_Eposition', 'f', ('river'))
nc.variables['river_Eposition'].long_name = "river ETA-position at RHO-points"
nc.variables['river_Eposition'].valid_min = 1
nc.variables['river_Eposition'].valid_max = 112
nc.variables['river_Eposition'][:]=[146,11]#[154,8]

nc.createVariable('river_direction', 'f', ('river'))
nc.variables['river_direction'].long_name = "river runoff direction"
nc.variables['river_direction'][:]=[0,0]#[1,0]

nc.createVariable('river_temp', 'f', ('river_time', 's_rho', 'river'))
nc.variables['river_temp'].long_name = "river runoff potential temperature"
nc.variables['river_temp'].units = "Celsius"

nc.createVariable('river_salt', 'f', ('river_time', 's_rho', 'river'))
nc.variables['river_salt'].long_name = "river runoff salinity"
nc.variables['river_salt'].units = "PSU"

nc.createVariable('river_transport', 'f', ('river_time', 'river'))
nc.variables['river_transport'].long_name = "river runoff vertically integrated mass transport"
nc.variables['river_transport'].units = "meter3 second-1"
nc.variables['river_transport'].time = "river_time"


nc.createVariable('river_Vshape', 'f', ('s_rho', 'river'))
nc.variables['river_Vshape'].long_name = "river runoff mass transport vertical profile"
vshape = np.zeros([N, Nr])
for k in range(N):
    vshape[k,:] = k

area = sum(vshape[:,0])
vshape = (1.0/area)*vshape

nc.variables['river_Vshape'][:] = vshape


YY=[]; MM=[]; DD=[]; QQ=[]
for line in open('Scheldt_Anwers.txt','r').readlines():
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

YY=[]; MM=[]; DD=[]; QQ=[]
for line in open('Thames_London.txt','r').readlines():
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
Thames=Q



Thames_temp=[]; Scheldt_temp=[]

for line in open('Temperature_all.txt','r').readlines():
	Thames_temp.append(float(line.split("\t")[1]))
	Scheldt_temp.append(float(line.split("\t")[2]))

Thames_t=np.zeros((365,15))
Scheldt_t=np.zeros((365,15))

nc.variables['river_transport'][:,0]=np.concatenate((-Thames,-Thames,-Thames,-Thames,-Thames,-Thames,-Thames,-Thames,-Thames,-Thames),axis=0)
nc.variables['river_transport'][:,1]=np.concatenate((Scheldt,Scheldt,Scheldt,Scheldt,Scheldt,Scheldt,Scheldt,Scheldt,Scheldt,Scheldt),axis=0)

for m in range(365):
	for n in range(len(vshape)): 
		Thames_t[m,n]=Thames_temp[m]
		Scheldt_t[m,n]=Scheldt_temp[m]


nc.variables['river_temp'][:,:,0]=np.concatenate((Thames_t,Thames_t,Thames_t,Thames_t,Thames_t,Thames_t,Thames_t,Thames_t,Thames_t,Thames_t),axis=0)
nc.variables['river_temp'][:,:,1]=np.concatenate((Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t),axis=0)

for i in range(0,2):
	nc.variables['river_salt'][:,:,i]=np.zeros((3650,15))+1.0

nc.close()

copyfile('Rivers_nested.nc', '/home/eivanov/COAWST/Data/ROMS/Rivers/Rivers_10_years_salt_nested.nc')
