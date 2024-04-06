import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors, colorbar
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import colorsys
import dask
import os.path
import math
from plotting_routines import plot_scalarfield2 

graphdir='../graphs/'
datadir='../SHiELD_SRC/test/CI/BATCH-CI/'
figformat='png'
rad2deg=180/np.pi

#--------------------------------------------------------------------------------------------------------
# test case
tc = -4

#N
N=192

# 1d advection scheme
hord = 0

#grid type
gtype = 2  # 0-equiedge; 2-equiangular

# duogrid parameter
dg = 2

# 2d advection scheme
adv = 2

#mf
mf = 0

#----------------------------------------------------------------------------------

basename='geobalance'
Tf = 12
alpha = 45  # rotation angle


error_div = np.zeros((N,N,6))
lat = np.zeros((N+1,N+1,6))
lon = np.zeros((N+1,N+1,6))

#--------------------------------------------------------------------------------------------------------
# grid name
if gtype==0:
   gname = 'equiedge'
elif gtype==2:
   gname = 'equiangular'

# duogrid name
if dg==1:
    dg = 'dg1'
elif dg==2:
    dg = 'dg2'
else:
    dg = 'kinked'

# Directory where the netcdf files are
filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
    +".g"+str(gtype)+"."+dg+".adv"+str(adv)+".hord"+str(hord)+'.mf'+str(mf)+".tf"+str(Tf)+"/"

# store the data
for p in range(1,7):
   z = open(filepath+'error_div_'+str(p)+'.dat', 'rb')
   z = np.fromfile(z, dtype=np.float64)
   z = np.reshape(z, (N,N))
   error_div[:,:,p-1] = z

   z = open(filepath+'lon_'+str(p)+'.dat', 'rb')
   z = np.fromfile(z, dtype=np.float64)
   z = np.reshape(z, (N+1,N+1))
   lon[:,:,p-1] = z*rad2deg
   

   z = open(filepath+'lat_'+str(p)+'.dat', 'rb')
   z = np.fromfile(z, dtype=np.float64)
   z = np.reshape(z, (N+1,N+1))
   lat[:,:,p-1] = z*rad2deg

title="C"+str(N)+"-g"+str(gtype)+"."+dg+'.adv'+str(adv)+'.hord'+str(hord)+'.mf'+str(mf)
filename='div_'+"g"+str(gtype)+"_"+str(N)+"."+dg+'.adv'+str(adv)+'.hord'+str(hord)+'.mf'+str(mf)
colormap='seismic'
map_projection='mercator'
dpi=100

print(np.amin(error_div), np.amax(error_div))
#error_div [error_div==0] = np.amin(error_div*0.00001)
#error_div = np.log10(error_div)
#qmin = -8
#qmax = -2.5
qmax_abs = np.amax(abs(error_div))
qmin, qmax = -qmax_abs, qmax_abs
plot_scalarfield2(error_div, lon, lat, map_projection, title, filename, colormap, dpi, figformat, qmin, qmax)
print(np.amin(error_div), np.amax(error_div), qmax_abs)
#------------------------------------------------------------------------------------------------
