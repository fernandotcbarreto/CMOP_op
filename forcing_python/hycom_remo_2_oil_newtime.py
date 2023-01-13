##delft 2 oil
import numpy as np
from netCDF4 import Dataset as dat
from scipy.interpolate import griddata, interp1d
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import matplotlib.dates as dates
import xarray as xr
from scipy.interpolate import interp2d

limn=20   #usar quando for arquivo 3D, arquivos 2D não precisa ignorar

#inifile='/mnt/c/Users/fernando.barreto/Documents/CMOP_OIL_SPILL-master/hycom_remo/input.dat'

lista=['/mnt/c/Users/fernando.barreto/Downloads/hycom_uv_20181104_analise.nc']

avgfile = xr.open_mfdataset(lista,  concat_dim='time', combine='nested')
avgfile=avgfile.sel(time=~avgfile.get_index("time").duplicated())
file=dat(lista[0])
time=file['time'][::]   #minutes since 2021-10-02

dateini=file['time'].units[-10:]
time_oil=time + dates.datestr2num(dateini) 
time_oil=time_oil*24*60 

lat= np.ma.filled(file['Y'][::],fill_value=0)

lon= np.ma.filled(file['X'][::],fill_value=0)

[lon,lat]=np.meshgrid(lon,lat)


utim=np.ma.filled(file['u'][:,::],fill_value=0)
a=np.isnan(utim)
utim[a]=0
vtim=np.ma.filled(file['v'][:,::],fill_value=0)
vtim[a]=0
ttim=vtim*0+25
stim=vtim*0+35


depth = np.ma.filled(file['u'][::],fill_value=19999999)
depth = depth[0,0,:,:]
depth[depth<19999999]=1000;depth[depth>=19999999]=-1000



#interp to another domain
latN= -21    #latnorth
latS= -25      #latSouth
lonW= -48.5    #lonWest
lonE= -38  #lonEast


X=lon.copy()
Y=lat.copy()
minlon = X[0,:] - lonW 
iml = np.where(np.absolute(minlon)==np.absolute(minlon).min())[0][0]
maxlon = X[0,:] - lonE
imxl = np.where(np.absolute(maxlon)==np.absolute(maxlon).min())[0][0]
 
minlat = Y[:,0] - latS
imla = np.where(np.absolute(minlat)==np.absolute(minlat).min())[0][0]
maxlat = Y[:,0] - latN
imxla = np.where(np.absolute(maxlat)==np.absolute(maxlat).min())[0][0]


lat=lat[imla:imxla,iml:imxl]
lon=lon[imla:imxla,iml:imxl]
depth=depth[imla:imxla,iml:imxl]
utim=utim[:,:,imla:imxla,iml:imxl]
vtim=vtim[:,:,imla:imxla,iml:imxl]
utim=np.squeeze(utim)
vtim=np.squeeze(vtim)

#plt.pcolor(lon,lat,depth);plt.colorbar();plt.show()

layer=file['depth'][:]

layer=0



counttimeh=np.array([0])

u,v=utim.copy(), vtim.copy()


[lonm,latm]=np.meshgrid(np.arange(lon.min(), lon.max(), np.diff(lon[0,:])[0]/3),np.arange(lat.min(), lat.max(), np.diff(lat[:,0])[0]/3))


depth = np.ma.filled(file['u'][::],fill_value=19999999)
depth = depth[0,0,:,:]
depth[depth<19999999]=1000;depth[depth>=19999999]=-1000
nlat,nlon=depth.shape
depth=depth[imla:imxla,iml:imxl]

hgeb=interp2d(lon[0,:],lat[:,0],depth,kind='linear')

depth=hgeb(lonm[0,:],latm[:,0])


#layer =  -np.ma.filled(file['depth'][::],fill_value=0)
#layer=layer[0:limn]

GEBCO=False
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

if GEBCO:
 depth=h.copy()
 depth[np.where(depth>-5)]=100


coastvalue=depth.copy()
coastvalue[coastvalue>=100]=0;
#coastvalue[coastvalue<0]=1000
coastvalue[coastvalue<0]=1000

with open('time_delft.txt', 'w') as f:
   np.savetxt(f, np.array(len(time_oil)).reshape(1,), fmt = '%i')
   np.savetxt(f, np.array(int(lon.shape[0])).reshape(1,), fmt = '%i')
   np.savetxt(f, np.array(int(lon.shape[1])).reshape(1,), fmt = '%i')   
   np.savetxt(f, np.array(1).reshape(1,), fmt = '%i')      
   np.savetxt(f, time_oil,fmt='%10.4f')
   np.savetxt(f, np.array(layer).reshape(1,),fmt='%10.4f')
   np.savetxt(f, counttimeh,fmt='%10.4f')
   np.savetxt(f, np.array(int(lonm.shape[0])).reshape(1,), fmt = '%i')
   np.savetxt(f, np.array(int(lonm.shape[1])).reshape(1,), fmt = '%i')
#   np.savetxt(f, time_oil,fmt="%s")


for i in range(len(time_oil)):
  with open(str(i+1) + 'v.txt', 'w') as f:
      np.savetxt(f, np.squeeze(vtim[i, :,:]),fmt='%10.8f')

for i in range(len(time_oil)):
  with open(str(i+1) + 'u.txt', 'w') as f:
      np.savetxt(f, np.squeeze(utim[i, :,:]),fmt='%10.8f')

for i in range(len(time_oil)):
  with open(str(i+1) + 'temp.txt', 'w') as f:
      np.savetxt(f, np.squeeze(ttim[i, :,:]),fmt='%10.8f')

for i in range(len(time_oil)):
  with open(str(i+1) + 'salt.txt', 'w') as f:
      np.savetxt(f, np.squeeze(stim[i, :,:]),fmt='%10.8f')

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


