import numpy as np
from netCDF4 import Dataset as dat
from scipy.interpolate import griddata, interp1d
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import matplotlib.dates as dates

from netCDF4 import Dataset as dat
import xarray as xr
import numpy as np

file=dat('/mnt/c/Users/fernando.barreto/Downloads/ADCIRC_BG_uv_20211104_graderegular/ADCIRC_BG_uv_20211103_graderegular.nc', 'r+')


time=file['time2'][::]/24   #hours since 2021-10-02


dateini = '2021-11-03'


time_oil=time + dates.datestr2num(dateini) 
time_oil=time_oil*24*60 

 

lat= np.ma.filled(file['latitude'][::],fill_value=0)

lon= np.ma.filled(file['longitude'][::],fill_value=0)


utim=np.ma.filled(file['u'][:,::],fill_value=0)
a=np.isnan(utim)
utim[a]=0
vtim=np.ma.filled(file['v'][:,::],fill_value=0)
vtim[a]=0
ttim=np.zeros(vtim.shape)+27
stim=np.zeros(vtim.shape)+35

 

latm= np.ma.filled(file['latitude'][::],fill_value=0)

lonm= np.ma.filled(file['longitude'][::],fill_value=0)


depth = np.ma.filled(file['u'][::],fill_value=19999999)
depth = depth[0,:,:]
depth[~np.isnan(depth)]=1;depth[np.isnan(depth)]=-1000


#plt.pcolor(lon, lat, depth);plt.show()
coastvalue=depth.copy()
coastvalue[coastvalue>=100]=0;
coastvalue[coastvalue<0]=1000
#coastvalue[coastvalue<0]=0


layer=0



counttimeh=np.array([0])

u,v=utim.copy(), vtim.copy()

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


with open('depth.txt', 'w') as f:
  np.savetxt(f, depth, fmt='%10.1f')


with open('lat_dmask.txt', 'w') as f:
     np.savetxt(f, latm,fmt='%10.8f')


with open('lon_dmask.txt', 'w') as f:
     np.savetxt(f, lonm,fmt='%10.8f')     


with open('depth.txt', 'w') as f:
  np.savetxt(f, depth, fmt='%10.1f')


with open('coastvalue.txt', 'w') as f:
  np.savetxt(f, coastvalue, fmt='%10.1f')


