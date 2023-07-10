##delft 2 oil
import numpy as np
from netCDF4 import Dataset as dat
from scipy.interpolate import griddata, interp1d, interp2d
import matplotlib.pyplot as plt
import datetime
import datetime as dt
import pandas as pd
import matplotlib.dates as dates
import xarray as xr
limn=1

path = '/mnt/c/Users/fernando.barreto/Desktop/boletins/current_mercator/'

path = '/mnt/c/Users/fernando.barreto/Downloads/'
#path='/home/oceanpact/Desktop/datadrive/dados/corrente/'

file=dat(path+'cmems_mod_glo_phy_anfc_merged-uv_PT1H-i_1684162784618.nc')
#filex=xr.open_dataset(path+'global-analysis-forecast-phy-001-024-hourly-t-u-v-ssh_1664979631192.nc')

yesterday = dt.datetime.today()
strys=yesterday.strftime("%Y%m%d")

#file=dat(path+'previsao_corrente_'+strys+'.nc') # tome


time=file['time'][::] + dates.datestr2num(file['time'].units[12:])*24   #minutes since 2021-10-02
time=time*60

lat= np.ma.filled(file['latitude'][::],fill_value=0)

lon= np.ma.filled(file['longitude'][::],fill_value=0)

[lon,lat]=np.meshgrid(lon,lat)


# utim=np.ma.filled(file['uo'][:,0:limn,::],fill_value=0)
# vtim=np.ma.filled(file['vo'][:,0:limn,::],fill_value=0)
# ttim=np.ma.filled(file['uo'][:,0:limn,::],fill_value=20)
# stim=np.ma.filled(file['uo'][:,0:limn,::],fill_value=35)
utim=np.ma.filled(file['utotal'][:,0:limn,::],fill_value=0)
vtim=np.ma.filled(file['vtotal'][:,0:limn,::],fill_value=0)
ttim=np.ma.filled(file['utotal'][:,0:limn,::],fill_value=20)
stim=np.ma.filled(file['utotal'][:,0:limn,::],fill_value=35)
ttim=np.ones(ttim.shape)*20
stim=np.ones(ttim.shape)*30

#depth = np.ma.filled(file['uo'][::],fill_value=19999999)
depth = np.ma.filled(file['utotal'][::],fill_value=19999999)
depth = depth[0,0,:,:]
depth[depth<19999999]=1000;depth[depth>=19999999]=-1000

layer =  -np.ma.filled(file['depth'][::],fill_value=0)
layer=layer[0:limn]


time_oil=time

u,v=utim.copy(), vtim.copy()



latm= np.ma.filled(file['latitude'][::],fill_value=0)

lonm= np.ma.filled(file['longitude'][::],fill_value=0)

[lonm,latm]=np.meshgrid(np.arange(lonm.min(), lonm.max(), np.diff(lonm)[0]/3),np.arange(latm.min(), latm.max(), np.diff(latm)[0]/3))


#depth = np.ma.filled(file['uo'][::],fill_value=19999999)
depth = np.ma.filled(file['utotal'][::],fill_value=19999999)
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

counttimeh=np.zeros([1])

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


