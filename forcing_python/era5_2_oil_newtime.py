##delft 2 oil
import numpy as np
from netCDF4 import Dataset as dat
from scipy.interpolate import griddata, interp1d
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import matplotlib.dates as dates
import xarray as xr
base='/home/fernando/fernando_noronha/2020/summer/era5/download/'

fileu=[base+'download_u_10_2019.nc', base+'download_u_10_2020.nc']

filev=[base+'download_v_10_2019.nc', base+'download_v_10_2020.nc']



uwind=xr.open_mfdataset(fileu)

vwind=xr.open_mfdataset(filev)


time=np.concatenate((dat(base+'download_u_10_2019.nc')['time'][:],dat(base+'download_u_10_2020.nc')['time'][:]), axis=0)

time = time + dates.datestr2num('1900-01-01')*24

dateinioil='2019-11-03 12:00:00'   ##deprecated ignore it

dates.datestr2num(dateinioil)*24

dates.num2date(time[0]/24)

#time=uwind['time'][::]

lon=uwind['longitude'][::]

lat=uwind['latitude'][::]

long, latg = np.meshgrid(lon,lat)


u10a=uwind['u10'][::]

v10a=vwind['v10'][::]

timeinioil = '2019-11-01 12:00:00'


u10=u10a[:]

v10=v10a[:]

timemod = time 

time_wind = timemod*60

counttimeh=np.array([dates.datestr2num(dateinioil)*24*60])


with open('time_era_5.txt', 'w') as f:
   np.savetxt(f, np.array(len(time_wind)).reshape(1,), fmt = '%i')                 #number of time series
   np.savetxt(f, np.array(int(long.shape[0])).reshape(1,), fmt = '%i')             #number of latitudes
   np.savetxt(f, np.array(int(long.shape[1])).reshape(1,), fmt = '%i')             #number of longitudes
   np.savetxt(f, time_wind,fmt='%10.4f')                                           #time gauging points
   np.savetxt(f, counttimeh,fmt='%10.4f')                                           #time gauging points


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

