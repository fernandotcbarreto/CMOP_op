
##delft 2 oil
import numpy as np
from netCDF4 import Dataset as dat
from scipy.interpolate import griddata, interp1d
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import matplotlib.dates as dates




from cmop_shape import seed_from_shapefile 
lon, lat=seed_from_shapefile(number=500, shapefile="/mnt/c/Users/fernando.barreto/Downloads/oleo_shape/oleo_virtual.shp")
lat=lat.round(4).astype('str')
lon=lon.round(4).astype('str')


lat=['-5','-3.51858386', '-7']
lon=['-35.','-35.89492945', '-33']
#lat=['-7','-7', '-7', '-7', '-7']
#lon=['-33','-33', '-33', '-33', '-33']
#lat=['-22.4240','-22.639', '-21.9067','-22.4686', '-20.0421', '-24.6487','-19.5305','-22.0883', '-22.7018']   #p53, 51,52, 26,vitoria,77 , peroa, p50, trident
#lon=['-39.9577','-40.09', '-39.7360', '-40.0288', '-39.5248', '-42.5145','-39.2570', '-39.8278', '-40.6771']
#spillt=np.array([48, 48, 48, 48, 48, 48, 48, 48, 48])  #hours of spill hours
#spilli=np.array([0, 24, 48])  #hours of spill hours
spilli=np.array([0, 0,0])  #hours of spill hours
spillt=np.array([24, 24, 24])  #hours of spill hours

try:
 spilli.shape
except:
  spilli=np.zeros(len(spillt))

#spillt=np.array([0,0,0,0, 0])  #hours of spill hours
#spillt=np.zeros(len(lat))
#depth=np.array([-48, -48, -48, -48, -48, -48, -48, -48, -48])  #hours of spill hours
#depth=np.array([0,-20,-40, -60, -100])  #hours of spill hours
depth=np.array([0,-20,-40]) 
#depth=np.zeros(len(lat))
vol=50


spillt=spillt*60
spilli=spilli*60

dt=900
#dt=3600


with open('testsurf.txt', 'w') as f:
 for i in range(len(spillt)):
#  a=np.arange(0,spillt[i]+(dt/60), dt/60)
  a=np.arange(spilli[i],spilli[i]+spillt[i]+(dt/60), dt/60)
  for ii in range(len(a)):
   lines='      1      '+str(a[ii])+'        '+lon[i]+'         '+lat[i]+'          '+str(depth[i])+'          '+str(vol/len(a))
 #  print(lines)
   f.write(lines)
   f.write('\n')

f.close()


#dateinioil='2019-11-03 12:00:00'

#counttimeh=np.array([dates.datestr2num(dateinioil)*24*60])

#chemis='    568          0.0270          0.0230         0.0295         0.0230           0.0573           0.0189          0.0075          0.0157          0.0202          0.0149          0.0014          0.0305          0.0537          0.0841          0.0612         0.0490          0.0567          0.0685          0.347          0.0023          0.0034          0.0026         0.0026         0.00000000          0.00000000'
chemis='568  0.011370 0.018143 0.011684 0.001947 0.018841 0.005193 0.020469 0.005592 0.022788 0.013767 0.035908 0.001238 0.038241 0.000061 0.002063 0.041303 0.040400 0.026897 0.001165 0.026539 0.020380 0.000910 0.075818 0.003286 0.555997'   
#np.sum(np.fromstring(chemis, dtype=float, sep=' '))
lines=[]
for juij in range(len(lat)):
 lines.append(chemis)
# lines = [  '    568          0.0270          0.0230         0.0295         0.0230           0.0573           0.0189          0.0075          0.0157          0.0202          0.0149          0.0014          0.0305          0.0537          0.0841          0.0612         0.0490          0.0567          0.0685          0.347          0.0023          0.0034          0.0026         0.0026         0.00000000          0.00000000'
# ,  '    568          0.0270          0.0230         0.0295         0.0230           0.0573           0.0189          0.0075          0.0157          0.0202          0.0149          0.0014          0.0305          0.0537          0.0841          0.0612         0.0490          0.0567          0.0685          0.347          0.0023          0.0034          0.0026         0.0026         0.00000000          0.00000000'
# ,  '    568          0.0270          0.0230         0.0295         0.0230           0.0573           0.0189          0.0075          0.0157          0.0202          0.0149          0.0014          0.0305          0.0537          0.0841          0.0612         0.0490          0.0567          0.0685          0.347          0.0023          0.0034          0.0026         0.0026         0.00000000          0.00000000'
# ,  '    568          0.0270          0.0230         0.0295         0.0230           0.0573           0.0189          0.0075          0.0157          0.0202          0.0149          0.0014          0.0305          0.0537          0.0841          0.0612         0.0490          0.0567          0.0685          0.347          0.0023          0.0034          0.0026         0.0026         0.00000000          0.00000000'
# ,  '    568          0.0270          0.0230         0.0295         0.0230           0.0573           0.0189          0.0075          0.0157          0.0202          0.0149          0.0014          0.0305          0.0537          0.0841          0.0612         0.0490          0.0567          0.0685          0.347          0.0023          0.0034          0.0026         0.0026         0.00000000          0.00000000'
# ,  '    568          0.0270          0.0230         0.0295         0.0230           0.0573           0.0189          0.0075          0.0157          0.0202          0.0149          0.0014          0.0305          0.0537          0.0841          0.0612         0.0490          0.0567          0.0685          0.347          0.0023          0.0034          0.0026         0.0026         0.00000000          0.00000000'
# ,  '    568          0.0270          0.0230         0.0295         0.0230           0.0573           0.0189          0.0075          0.0157          0.0202          0.0149          0.0014          0.0305          0.0537          0.0841          0.0612         0.0490          0.0567          0.0685          0.347          0.0023          0.0034          0.0026         0.0026         0.00000000          0.00000000'
# ,  '    568          0.0270          0.0230         0.0295         0.0230           0.0573           0.0189          0.0075          0.0157          0.0202          0.0149          0.0014          0.0305          0.0537          0.0841          0.0612         0.0490          0.0567          0.0685          0.347          0.0023          0.0034          0.0026         0.0026         0.00000000          0.00000000'
# ,  '    568          0.0270          0.0230         0.0295         0.0230           0.0573           0.0189          0.0075          0.0157          0.0202          0.0149          0.0014          0.0305          0.0537          0.0841          0.0612         0.0490          0.0567          0.0685          0.347          0.0023          0.0034          0.0026         0.0026         0.00000000          0.00000000']

with open('testfrac.txt', 'w') as f:
 for i in range(len(spillt)):
  a=np.arange(0,spillt[i]+(dt/60), dt/60)
  for ii in range(len(a)):
   f.write(lines[i])
   f.write('\n')

f.close()



ff=pd.read_csv('testsurf.txt', header=None, delim_whitespace=True)

bb=ff.sort_values(by=1).index
ff=ff.iloc[bb]

np.savetxt('testsurf.txt', ff.values, fmt='%5.5f')



ff=pd.read_csv('testfrac.txt', header=None, delim_whitespace=True)
ff=ff.iloc[bb]

np.savetxt('testfrac.txt', ff.values, fmt='%5.5f')