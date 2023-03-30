##delft 2 oil
import numpy as np
from netCDF4 import Dataset as dat
from scipy.interpolate import griddata, interp1d
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import matplotlib.dates as dates
from scipy.interpolate import interp2d


limn=1

path = '/home/fernando/fernando_noronha/2020/summer/mercator/download/fma/'
path = '/mnt/c/Users/fernando.barreto/Downloads/mercatord/'

#file=dat(path+'MYOCEAN_AZUL_FORECAST_20200413.nc')
file=dat(path+'MYOCEAN_AZUL_FORECAST_20200108.nc')

input = path+'MYOCEAN_AZUL_FORECAST_'
endfname='.nc'

dateinioil='2020-01-01 12:00:00'

dateini='2019-11-01'
dateend='2020-04-29'
dateini='2020-01-01'
dateend='2020-04-01'

a=[dateini, dateend]

date=dates.datestr2num(a)

ntimes=len(np.arange(date[0], date[1]))

date = np.arange(date[0], date[1])

time_oil=np.arange(ntimes) + dates.datestr2num(dateini) 
time_oil=np.arange(ntimes) + dates.datestr2num(dateinioil) 

time_oil=time_oil*24*60 

lat= np.ma.filled(file['latitude'][::],fill_value=0)

lon= np.ma.filled(file['longitude'][::],fill_value=0)

[lon,lat]=np.meshgrid(lon,lat)

u = np.ma.filled(file['uo'][::],fill_value=0)
v = np.ma.filled(file['vo'][::],fill_value=0)

utim=np.zeros([ntimes,limn, u.shape[2], u.shape[3]])
vtim=np.zeros([ntimes,limn, v.shape[2], v.shape[3]])
ttim=np.zeros([ntimes,limn, v.shape[2], v.shape[3]])
stim=np.zeros([ntimes,limn, v.shape[2], v.shape[3]])


latm= np.ma.filled(file['latitude'][::],fill_value=0)

lonm= np.ma.filled(file['longitude'][::],fill_value=0)
# np.arange(latm.min(), latm.max(), np.diff(latm)[0]/3)
# np.arange(lonm.min(), lonm.max(), np.diff(lonm)[0]/3)

[lonm,latm]=np.meshgrid(np.arange(lonm.min(), lonm.max(), np.diff(lonm)[0]/3),np.arange(latm.min(), latm.max(), np.diff(latm)[0]/3))


depth = np.ma.filled(file['uo'][::],fill_value=19999999)
depth = depth[0,0,:,:]
depth[depth<19999999]=1;depth[depth>=19999999]=-1000
nlat,nlon=depth.shape

hgeb=interp2d(lon[0,:],lat[:,0],depth,kind='linear')
depth=hgeb(lonm[0,:],latm[:,0])


layer =  -np.ma.filled(file['depth'][::],fill_value=0)
layer=layer[0:limn]

GEBCO=True
if GEBCO:
  print('bathymetry from GEBCO')
  gebco=dat('/mnt/c/Users/fernando.barreto/Desktop/ROMS/make_grid/gebco_222.nc')
  latgebco=gebco['lat'][:]
  longebco=gebco['lon'][:]
  elevation=-gebco['elevation'][:]
  from scipy.interpolate import interp2d
  hgeb=interp2d(longebco,latgebco,elevation,kind='linear')
  #lonbat,latbat=np.meshgrid(lonbat,latbat)
  h=hgeb(lonm[0,:],latm[:,0])
  #h=griddata((lonbat.ravel(),latbat.ravel()),gebco.ravel(),(lon.ravel(),lat.ravel())).reshape(lat.shape)

plt.pcolor(depth)
plt.savefig('/mnt/c/Users/fernando.barreto/Downloads/dep.png', dpi=200, transparent=False, bbox_inches="tight") 

if GEBCO:
 depth=h.copy()
 depth[np.where(depth>-5)]=100


coastvalue=depth.copy()
coastvalue[coastvalue>=100]=0;
#coastvalue[coastvalue<0]=1000
coastvalue[coastvalue<0]=0


for i in range(ntimes):
  fileu=input +dates.num2date(date[i]+1).strftime("%Y%m%d")+endfname
  filev=input +dates.num2date(date[i]+1).strftime("%Y%m%d")+endfname
  filet=input +dates.num2date(date[i]+1).strftime("%Y%m%d")+endfname
  files=input +dates.num2date(date[i]+1).strftime("%Y%m%d")+endfname
  fileu = dat(fileu)
  filev = dat(filev)
  filet = dat(filet)
  files = dat(files)
  u = np.ma.filled(fileu['uo'][:,0:limn,::],fill_value=0)
  v = np.ma.filled(filev['vo'][:,0:limn,::],fill_value=0)
#  t = np.ma.filled(filet['thetao'][:,0:limn,::],fill_value=20)
#  s = np.ma.filled(files['so'][:,0:limn,::],fill_value=35)
  t = np.zeros(u.shape)+20
  s = np.zeros(u.shape)+35
  utim[i,::] = u
  vtim[i,::] = v
  ttim[i,::] = t
  stim[i,::] = s

#counttimeh=np.zeros([1])
#counttimeh=np.array([10.5*60])
counttimeh=np.array([dates.datestr2num(dateinioil)*24*60])

u,v=utim.copy(), vtim.copy()

with open('time_delft.txt', 'w') as f:
   np.savetxt(f, np.array(len(time_oil)).reshape(1,), fmt = '%i')
   np.savetxt(f, np.array(int(lon.shape[0])).reshape(1,), fmt = '%i')
   np.savetxt(f, np.array(int(lon.shape[1])).reshape(1,), fmt = '%i')   
   np.savetxt(f, np.array(int(u.shape[1])).reshape(1,), fmt = '%i')      
   np.savetxt(f, time_oil,fmt='%10.4f')
   np.savetxt(f, layer,fmt='%10.4f')
   np.savetxt(f, counttimeh,fmt='%10.4f')
   np.savetxt(f, np.array(int(lonm.shape[0])).reshape(1,), fmt = '%i')
   np.savetxt(f, np.array(int(lonm.shape[1])).reshape(1,), fmt = '%i')    
#   np.savetxt(f, time_oil,fmt="%s")


for i in range(len(time_oil)):
  with open(str(i+1) + 'v.txt', 'w') as f:
    for j in range(v.shape[1]):
      np.savetxt(f, np.squeeze(vtim[i,j, :,:]),fmt='%10.8f')

for i in range(len(time_oil)):
  with open(str(i+1) + 'u.txt', 'w') as f:
    for j in range(u.shape[1]):
      np.savetxt(f, np.squeeze(utim[i,j, :,:]),fmt='%10.8f')

for i in range(len(time_oil)):
  with open(str(i+1) + 'temp.txt', 'w') as f:
    for j in range(u.shape[1]):
      np.savetxt(f, np.squeeze(ttim[i,j, :,:]),fmt='%10.8f')

for i in range(len(time_oil)):
  with open(str(i+1) + 'salt.txt', 'w') as f:
    for j in range(u.shape[1]):
      np.savetxt(f, np.squeeze(stim[i,j, :,:]),fmt='%10.8f')

with open('lat_delft.txt', 'w') as f:
     np.savetxt(f, lat,fmt='%10.8f')


with open('lon_delft.txt', 'w') as f:
     np.savetxt(f, lon,fmt='%10.8f')

with open('lat_dmask.txt', 'w') as f:
     np.savetxt(f, latm,fmt='%10.8f')


with open('lon_dmask.txt', 'w') as f:
     np.savetxt(f, lonm,fmt='%10.8f')     


with open('depth.txt', 'w') as f:
  np.savetxt(f, depth, fmt='%10.1f')


with open('coastvalue.txt', 'w') as f:
  np.savetxt(f, coastvalue, fmt='%10.1f')



