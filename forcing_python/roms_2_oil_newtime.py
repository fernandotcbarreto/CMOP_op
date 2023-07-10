from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from netCDF4 import Dataset 
from scipy.interpolate import griddata, interp1d
import numpy as np
import matplotlib.pyplot as plt
import conda
import os
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
from mpl_toolkits.basemap import Basemap
from matplotlib.pylab import *
from matplotlib import dates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import cmocean
import xarray as xr 
import pandas as pd
import glob 
import datetime as dt


#set limits to trim roms input to speed interpolation

latN= -20.5    #latnorth
latS= -24.5      #latSouth
lonW= -43.8    #lonWest
lonE= -38.5  #lonEast
#



path='/home/oceanpact/Desktop/datadrive/Previsao_MODELAGEM/'

yesterday = dt.datetime.today() - datetime.timedelta(days=2)
#yesterday = dt.datetime.today()

strys=yesterday.strftime("%Y-%m-%d")

file=path+'previsao_'+strys+'_ROMS.nc'

lista = sorted(glob.glob(file)) #SUL 1/36


avgfile = xr.open_mfdataset(lista,  concat_dim='ocean_time', combine='nested')
avgfile=avgfile.sel(ocean_time=~avgfile.get_index("ocean_time").duplicated())

grd = avgfile
x_roms = np.array(grd.variables['lon_rho'][:])
y_roms = np.array(grd.variables['lat_rho'][:])
msk_roms = np.array(grd.variables['mask_rho'][:])
msk_romsv = np.array(grd.variables['mask_v'][:])
h_roms = np.array(grd.variables['h'][:])[0]
h_roms1=np.array(h_roms.copy())
sh2 = np.array(y_roms.shape)
etamax, ximax = sh2


theta_b = np.array(avgfile.variables['theta_b'])[0]
theta_s = np.array(avgfile.variables['theta_s'])[0]
tcline = np.array(avgfile.variables['Tcline'])[0]
klevels = len(np.array(avgfile.variables['v_northward'][0,:,0,0]))
Vtransform = np.array(avgfile.variables['Vtransform'])[0]
Vstretching = np.array(avgfile.variables['Vstretching'])[0]
Spherical = True

if Vstretching==4:
  scoord = s_coordinate_4(h_roms, theta_b, theta_s, tcline, klevels)   #zeta is not used in the computation 
elif Vstretching==2:
  scoord = s_coordinate_2(h_roms, theta_b, theta_s, tcline, klevels)
elif Vstretching==1:
  scoord = s_coordinate(h_roms, theta_b, theta_s, tcline, klevels)

zr = -scoord.z_r[:]
 
zc=np.array([0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150])
#zc=np.array([0, 30])

#zc=zc[::-1]

uavg=np.array(avgfile.variables['u_eastward'][:])
#uavg=np.array(avgfile.resample(ocean_time='1D').mean('ocean_time').variables['Uwind'][:])

vavg=np.array(avgfile.variables['v_northward'][:])
#vavg=np.array(avgfile.resample(ocean_time='1D').mean('ocean_time').variables['v_northward'][:])

uwind=np.array(avgfile.variables['Uwind'][:])
#uavg=np.array(avgfile.resample(ocean_time='1D').mean('ocean_time').variables['Uwind'][:])

vwind=np.array(avgfile.variables['Vwind'][:])
#vavg=np.array(avgfile.resample(ocean_time='1D').mean('ocean_time').variables['v_northward'][:])


time=np.array(avgfile.variables['ocean_time'][:])

tempavg=np.array(avgfile.variables['temp'][:])


saltavg=np.array(avgfile.variables['salt'][:])

# minlon = x_roms[0,:] - lonW 
# iml = np.where(np.absolute(minlon)==np.absolute(minlon).min())[0][0]
# maxlon = x_roms[0,:] - lonE
# imxl = np.where(np.absolute(maxlon)==np.absolute(maxlon).min())[0][0]
 
# minlat = y_roms[:,0] - latS
# imla = np.where(np.absolute(minlat)==np.absolute(minlat).min())[0][0]
# maxlat = y_roms[:,0] - latN
# imxla = np.where(np.absolute(maxlat)==np.absolute(maxlat).min())[0][0]

# x_roms=x_roms[imla:imxla,iml:imxl]
# y_roms=y_roms[imla:imxla,iml:imxl]
# uavg=uavg[:,:,imla:imxla,iml:imxl]
# vavg=vavg[:,:,imla:imxla,iml:imxl]
# tempavg=tempavg[:,:,imla:imxla,iml:imxl]
# saltavg=saltavg[:,:,imla:imxla,iml:imxl]
# uwind=uwind[:,imla:imxla,iml:imxl]
# vwind=vwind[:,imla:imxla,iml:imxl]

##time=time[:]/(24*60*60)
intu=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
intv=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
itemp=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
isalt=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
zeta=avgfile['zeta'][:]


for j in range(intu.shape[2]):
  for k in range(intu.shape[3]):
    if (zr[-1,j,k] > zc.min()):
      zr[-1,j,k] = zc.min()

UNDEF=np.nan

for i in range(intu.shape[0]):
  for j in range(intu.shape[2]):
    for k in range(intu.shape[3]):
      intu[i,:,j,k] = np.interp(-zc, -zr[:,j,k], uavg[i,:,j,k], right=UNDEF, left=UNDEF)
      intv[i,:,j,k] = np.interp(-zc, -zr[:,j,k], vavg[i,:,j,k], right=UNDEF, left=UNDEF)
      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)
      isalt[i,:,j,k] = np.interp(-zc, -zr[:,j,k], saltavg[i,:,j,k], right=UNDEF, left=UNDEF)

layer=-zc

lon=x_roms.copy()
lat=y_roms.copy()

depth = intu[0,0,:].copy()
a=np.isnan(depth)
depth[~a]=1000;depth[a]=-1000

#plt.pcolor(depth);plt.colorbar();plt.show()

utim=intu.copy()
vtim=intv.copy()
ttim=itemp.copy()
stim=isalt.copy()

utim[np.isnan(utim)]=0
vtim[np.isnan(vtim)]=0
ttim[np.isnan(ttim)]=np.max(np.ma.masked_invalid(ttim))
stim[np.isnan(stim)]=np.max(np.ma.masked_invalid(stim))
#uwind1[np.isnan(uwind1)]=0
#vwind1[np.isnan(vwind1)]=0
utim=utim.round(3)
vtim=vtim.round(3)
ttim=ttim.round(3)
stim=stim.round(3)

u,v=utim.copy(), vtim.copy()


times=list()
for i in range(len(time)):
  times.append(dates.datestr2num(pd.to_datetime(str(time[i])).strftime("%Y%m%d-%H%M%S")))

timer=np.array(times)*24*60



time_oil=timer

u,v=utim.copy(), vtim.copy()

with open('time_delft.txt', 'w') as f:
   np.savetxt(f, np.array(len(time_oil)).reshape(1,), fmt = '%i')
   np.savetxt(f, np.array(int(lon.shape[0])).reshape(1,), fmt = '%i')
   np.savetxt(f, np.array(int(lon.shape[1])).reshape(1,), fmt = '%i')   
   np.savetxt(f, np.array(int(u.shape[1])).reshape(1,), fmt = '%i')      
   np.savetxt(f, time_oil,fmt='%10.4f')
   np.savetxt(f, layer,fmt='%10.4f')
#   np.savetxt(f, counttimeh,fmt='%10.4f')
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


with open('depth.txt', 'w') as f:
  np.savetxt(f, depth, fmt='%10.1f')



