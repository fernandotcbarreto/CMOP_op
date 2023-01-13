##delft 2 oil
import numpy as np
from netCDF4 import Dataset as dat
from scipy.interpolate import griddata, interp1d
import matplotlib.pyplot as plt
import datetime
import datetime as dt
import pandas as pd
import matplotlib.dates as dates

#base='/home/fernando/emergencia/gfs/'
#fileu=base + 'emergencia_gfs_u2.nc'
#filev=base + 'emergencia_gfs_v2.nc'
base='/mnt/c/Users/fernando.barreto/Desktop/boletins/wind_gfs/download_gfs_python/'
#fileu=base + 'GFS_F_UWIND_20220916.nc'
#filev=base + 'GFS_F_VWIND_20220916.nc'


yesterday = dt.datetime.today() - datetime.timedelta(days=2)
strys=yesterday.strftime("%Y%m%d")
base='/home/oceanpact/Desktop/datadrive/dados/vento/'
fileu=base+'previsao_vento_u_'+strys+'.nc' # tome
filev=base+'previsao_vento_v_'+strys+'.nc' # tome


uwind=dat(fileu)

vwind=dat(filev)

time=uwind['time'][::]  + dates.datestr2num(uwind['time'].units[11:]) #minutes since 2021-10-02 PRESTAR ATENÇÃO SE É MINUTES MSMM
time=time*24*60

lon=uwind['lon'][::] -360

lat=uwind['lat'][::]

long, latg = np.meshgrid(lon,lat)


u10a=uwind['ugrd10m'][::]

v10a=vwind['vgrd10m'][::]

#2021-11-04 13:20:00

u10=u10a[:]

v10=v10a[:]

time_wind = time

#counttimeh = counttimeh*60

with open('time_era_5.txt', 'w') as f:
   np.savetxt(f, np.array(len(time_wind)).reshape(1,), fmt = '%i')                 #number of time series
   np.savetxt(f, np.array(int(long.shape[0])).reshape(1,), fmt = '%i')             #number of latitudes
   np.savetxt(f, np.array(int(long.shape[1])).reshape(1,), fmt = '%i')             #number of longitudes
   np.savetxt(f, time_wind,fmt='%10.4f')                                           #time gauging points

for i in range(len(time_wind)):
  with open(str(i+1) + 'u_wind.txt', 'w') as f:
      np.savetxt(f, np.squeeze(u10[i,:,:]),fmt='%12.8f')

for i in range(len(time_wind)):
  with open(str(i+1) + 'v_wind.txt', 'w') as f:
      np.savetxt(f, np.squeeze(v10[i,:,:]),fmt='%12.8f')



with open('lat_era_5.txt', 'w') as f:
     np.savetxt(f, latg,fmt='%12.8f')



with open('lon_era_5.txt', 'w') as f:
     np.savetxt(f, long,fmt='%12.8f')

