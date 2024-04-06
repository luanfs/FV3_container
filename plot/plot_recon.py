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
tc = -10

#N
N=192

# 1d advection scheme
hord = 8

#grid type
gtype = 2  # 0-equiedge; 2-equiangular

# duogrid parameter
dg = 2

#midpoint
mp = 'c'
#mp = 's'

#----------------------------------------------------------------------------------

basename='recon-test'
alpha = 45  # rotation angle


error_recon = np.zeros((N,N,6))
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
+".g"+str(gtype)+"."+mp+"."+dg+'.hord'+str(hord)+"/"

# store the data
for p in range(1,7):
   z = open(filepath+'error_ppm_hord'+str(hord)+'_'+str(p)+'.dat', 'rb')
   z = np.fromfile(z, dtype=np.float64)
   z = np.reshape(z, (N,N))
   error_recon[:,:,p-1] = z

   z = open(filepath+'lon_'+str(p)+'.dat', 'rb')
   z = np.fromfile(z, dtype=np.float64)
   z = np.reshape(z, (N+1,N+1))
   lon[:,:,p-1] = z*rad2deg
   

   z = open(filepath+'lat_'+str(p)+'.dat', 'rb')
   z = np.fromfile(z, dtype=np.float64)
   z = np.reshape(z, (N+1,N+1))
   lat[:,:,p-1] = z*rad2deg

title="C"+str(N)+"-g"+str(gtype)+"."+mp+"."+dg+'.hord'+str(hord)
filename='recon_'+"g"+str(gtype)+"_"+str(N)+"."+mp+"."+dg+'.hord'+str(hord)
colormap='jet'
map_projection='mercator'
dpi=100

print(np.amin(error_recon), np.amax(error_recon))
error_recon [error_recon==0] = np.amin(error_recon*0.00001)
error_recon = np.log10(error_recon)
qmin = -8
qmax = -2.5
plot_scalarfield2(error_recon, lon, lat, map_projection, title, filename, colormap, dpi, figformat, qmin, qmax)
print(np.amin(error_recon), np.amax(error_recon))
#------------------------------------------------------------------------------------------------



