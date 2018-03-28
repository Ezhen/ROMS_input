import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from netCDF4 import Dataset
from shutil import copyfile

N=15
Nr=3

nc = Dataset('Rivers.nc', 'w', format='NETCDF3_64BIT')
nc.Description = 'Discharges of major rivers'
nc.Author = 'Evgeny Ivanov'
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'Discharges of major rivers'
nc.createDimension('river_time', 365*10)					###
####nc.createDimension('river_time', 365)
nc.createDimension('xi_rho', 78)
nc.createDimension('xi_u', 77)
nc.createDimension('xi_v', 78)
nc.createDimension('eta_rho', 112)
nc.createDimension('eta_u', 112)
nc.createDimension('eta_v', 111)
nc.createDimension('s_rho', N)
nc.createDimension('river', Nr)


YY=[]; MM=[]; DD=[]; QQ=[]

for line in open('Maas_Dordrecht.txt','r').readlines():
	YY.append(float(line[0:4]))
	MM.append(float(line[5:7]))
	DD.append(float(line[8:10]))
	QQ.append(float(line[11:16]))

M=MM[365:730]; D=DD[365:730]; Q=np.zeros((365)); n=0; time=np.zeros((365*10))

Maas=np.zeros((365));Scheldt=np.zeros((365));Thames=np.zeros((365));Rhine=np.zeros((365));Seine=np.zeros((365))

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
	time[k]=k+1
Maas=Q

nc.createVariable('river_time', 'f', ('river_time'))
nc.variables['river_time'].long_name = "river discharge time"
nc.variables['river_time'].units = "days since 2004-01-01 00:00:00"
nc.variables['river_time'][:]=np.linspace(0,3649,3650)	###

#nc.createVariable('Maas', 'f', ('river_time'))
#nc.variables['Maas'].long_name = "Discharge of Maas at Dordrecht"
#nc.variables['Maas'].units = "m3 s-1"
#nc.variables['Maas'].coordinates = "time"
#nc.variables['Maas'][:] = Q

nc.createVariable('river', 'i4', ('river'))
nc.variables['river'].long_name = "river runoff identification number: 1)Thames, 3)Rhine+Maas 4)Seine"
nc.variables['river'][:]=[1,3,4]

nc.createVariable('river_Xposition', 'f', ('river'))
nc.variables['river_Xposition'].long_name = "river XI-position at RHO-points"
nc.variables['river_Xposition'].valid_min = 1
nc.variables['river_Xposition'].valid_max = 82
nc.variables['river_Xposition'][:]=[33,17,7]

nc.createVariable('river_Eposition', 'f', ('river'))
nc.variables['river_Eposition'].long_name = "river ETA-position at RHO-points"
nc.variables['river_Eposition'].valid_min = 1
nc.variables['river_Eposition'].valid_max = 112
nc.variables['river_Eposition'][:]=[63,26,95]

nc.createVariable('river_direction', 'f', ('river'))
nc.variables['river_direction'].long_name = "river runoff direction"
nc.variables['river_direction'][:]=[1,0,0]

nc.createVariable('river_temp', 'f', ('river_time', 's_rho', 'river'))
nc.variables['river_temp'].long_name = "river runoff potential temperature"
nc.variables['river_temp'].units = "Celsius"

nc.createVariable('river_salt', 'f', ('river_time', 's_rho', 'river'))
nc.variables['river_salt'].long_name = "river runoff salinity"
nc.variables['river_salt'].units = "PSU"


nc.createVariable('river_Vshape', 'f', ('s_rho', 'river'))
nc.variables['river_Vshape'].long_name = "river runoff mass transport vertical profile"
vshape = np.zeros([N, Nr])
for k in range(N):
    vshape[k,:] = k

area = sum(vshape[:,0])
vshape = (1.0/area)*vshape

nc.variables['river_Vshape'][:] = vshape


YY=[]; MM=[]; DD=[]; QQ=[]
for line in open('Rhine_Delft.txt','r').readlines():
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
Rhine=Q

nc.createVariable('river_transport', 'f', ('river_time', 'river'))
nc.variables['river_transport'].long_name = "river runoff vertically integrated mass transport"
nc.variables['river_transport'].units = "meter3 second-1"
nc.variables['river_transport'].time = "river_time"
#nc.variables['river_transport'][:] = Q

#nc.createVariable('Rhine', 'f', ('river_time'))
#nc.variables['Rhine'].long_name = "Discharge of Rhine at Delft"
#nc.variables['Rhine'].units = "m3 s-1"
#nc.variables['Rhine'].coordinates = "time"
#nc.variables['Rhine'][:] = Q



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


#nc.createVariable('Scheldt', 'f', ('river_time'))
#nc.variables['Scheldt'].long_name = "Discharge of Scheldt at Anwers"
#nc.variables['Scheldt'].units = "m3 s-1"
#nc.variables['Scheldt'].coordinates = "time"
#nc.variables['Scheldt'][:] = Q



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

#nc.createVariable('Thames', 'f', ('river_time'))
#nc.variables['Thames'].long_name = "Discharge of Thames at London"
#nc.variables['Thames'].units = "m3 s-1"
#nc.variables['Thames'].coordinates = "time"
#nc.variables['Thames'][:] = Q



YY=[]; MM=[]; DD=[]; QQ=[]
for line in open('Seine_Le_Havre.txt','r').readlines():
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
Seine=Q



Thames_temp=[]; Scheldt_temp=[]; Maas_temp=[]; Rhine_temp=[]; Seine_temp=[]

for line in open('Temperature_all.txt','r').readlines():
	Thames_temp.append(float(line.split("\t")[1]))
	Scheldt_temp.append(float(line.split("\t")[2]))
	Maas_temp.append(float(line.split("\t")[3]))
	Rhine_temp.append(float(line.split("\t")[4]))
	Seine_temp.append(float(line.split("\t")[5]))


#nc.createVariable('Maas_temp', 'f', ('river_time'))
#nc.createVariable('Rhine_temp', 'f', ('river_time'))
#nc.createVariable('Scheldt_temp', 'f', ('river_time'))
#nc.createVariable('Thames_temp', 'f', ('river_time'))
#nc.variables['Thames_temp'].units = "Celcius"; nc.variables['Scheldt_temp'].units = "Celcius"; nc.variables['Maas_temp'].units = "Celcius"; nc.variables['Rhine_temp'].units = "Celcius"
#nc.variables['Thames_temp'][:] = Thames_temp
#nc.variables['Scheldt_temp'][:] = Scheldt_temp
#nc.variables['Maas_temp'][:] = Maas_temp
#nc.variables['Rhine_temp'][:] = Rhine_temp
'''
Thames_t=np.zeros((365,15))
Scheldt_t=np.zeros((365,15))
Maas_t=np.zeros((365,15))
Rhine_t=np.zeros((365,15))
Seine_t=np.zeros((365,15))
'''
Thames_t=np.zeros((365,15))
Scheldt_t=np.zeros((365,15))
Maas_t=np.zeros((365,15))
Rhine_t=np.zeros((365,15))
Seine_t=np.zeros((365,15))
'''
nc.variables['river_transport'][:,0]=-Thames
nc.variables['river_transport'][:,1]=Scheldt
nc.variables['river_transport'][:,2]=Rhine+Maas
#nc.variables['river_transport'][:,2]=Maas
#nc.variables['river_transport'][:,3]=-Rhine
nc.variables['river_transport'][:,3]=Seine
'''

nc.variables['river_transport'][:,0]=np.concatenate((-Thames,-Thames,-Thames,-Thames,-Thames,-Thames,-Thames,-Thames,-Thames,-Thames),axis=0)
#nc.variables['river_transport'][:,1]=np.concatenate((Scheldt,Scheldt,Scheldt,Scheldt,Scheldt,Scheldt,Scheldt,Scheldt,Scheldt,Scheldt),axis=0)
nc.variables['river_transport'][:,1]=np.concatenate((Rhine+Maas,Rhine+Maas,Rhine+Maas,Rhine+Maas,Rhine+Maas,Rhine+Maas,Rhine+Maas,Rhine+Maas,Rhine+Maas,Rhine+Maas),axis=0)
#nc.variables['river_transport'][:,2]=Maas
#nc.variables['river_transport'][:,3]=-Rhine
nc.variables['river_transport'][:,2]=np.concatenate((Seine,Seine,Seine,Seine,Seine,Seine,Seine,Seine,Seine,Seine),axis=0)


	
'''
nc.variables['river_transport'][:,0]=10000
nc.variables['river_transport'][:,1]=10000
nc.variables['river_transport'][:,2]=-10000
'''
for m in range(365):
	for n in range(len(vshape)): 
		Thames_t[m,n]=Thames_temp[m]
		Scheldt_t[m,n]=Scheldt_temp[m]
		Maas_t[m,n]=Maas_temp[m]
		Rhine_t[m,n]=Rhine_temp[m]
		Seine_t[m,n]=Seine_temp[m]

'''
nc.variables['river_temp'][:,:,0]=Thames_t
nc.variables['river_temp'][:,:,1]=Scheldt_t
nc.variables['river_temp'][:,:,2]=0.15*Maas_t+0.85*Rhine_t
#nc.variables['river_temp'][:,:,2]=Maas_t
#nc.variables['river_temp'][:,:,3]=Rhine_t
nc.variables['river_temp'][:,:,3]=Seine_t
'''

nc.variables['river_temp'][:,:,0]=np.concatenate((Thames_t,Thames_t,Thames_t,Thames_t,Thames_t,Thames_t,Thames_t,Thames_t,Thames_t,Thames_t),axis=0)
#nc.variables['river_temp'][:,:,1]=np.concatenate((Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t,Scheldt_t),axis=0)
nc.variables['river_temp'][:,:,1]=np.concatenate((0.15*Maas_t+0.85*Rhine_t,0.15*Maas_t+0.85*Rhine_t,0.15*Maas_t+0.85*Rhine_t,0.15*Maas_t+0.85*Rhine_t,0.15*Maas_t+0.85*Rhine_t,0.15*Maas_t+0.85*Rhine_t,0.15*Maas_t+0.85*Rhine_t,0.15*Maas_t+0.85*Rhine_t,0.15*Maas_t+0.85*Rhine_t,0.15*Maas_t+0.85*Rhine_t),axis=0)
#nc.variables['river_temp'][:,:,2]=Maas_t
#nc.variables['river_temp'][:,:,3]=Rhine_t
nc.variables['river_temp'][:,:,2]=np.concatenate((Seine_t,Seine_t,Seine_t,Seine_t,Seine_t,Seine_t,Seine_t,Seine_t,Seine_t,Seine_t),axis=0)

for i in range(0,3):
	nc.variables['river_salt'][:,:,i]=np.zeros((3650,15))+1.
'''
nc.variables['river_temp'][:,:,0]=30
nc.variables['river_temp'][:,:,1]=30
nc.variables['river_temp'][:,:,2]=30
'''



nc.close()

copyfile('Rivers.nc', '/home/eivanov/coawst_data_prrocessing/Temporal/Input_files_to_ROMS/Rivers__western_boundary.nc')
