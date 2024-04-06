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

# 1d advection scheme
#hords = (0,8)
hords = (0, )

#grid type 0-equiedge; 2-equiangular
gtypes = (0,0,0,0)
#gtypes = (2,2,2,2)
#gtypes = (2,2)

#mass fixer
mfs = (0,1,0,1)

#adv scheme
advs = (1,1,2,2)
#advs = (2,2)

#duogrid
dgs = (2,2,2,2)

# values of N
#Ns = (48, 96, 192, 384, 768 )
Ns = (48, 96, 192, 384)
#Ns = (48, 96, 192, )
#Ns = (48, 96, )
ngrids = len(Ns)
#--------------------------------------------------------------------------------------------------------

basename='geobalance'
Tf=12
alpha = 45  # rotation angle
M = len(gtypes)
K = len(hords)
errors = np.zeros((len(Ns),K,M))

for n in range(0,len(Ns)):
    N = Ns[n]
    for k in range(0,K):
        hord = hords[k]
        for m in range(0,M):
            gtype = gtypes[m]
            dg    = dgs[m]
            mf    = mfs[m]
            adv   = advs[m]

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

            errors[n,k,m] = np.amax(abs(z))

#colors = ('green', 'red', 'blue', 'purple', 'orange', 'black', 'brown', 'gray')
#colors = ('green', 'red', 'blue', 'yellow', 'green', 'red', 'blue', 'yellow')
colors = ('green', 'red', 'blue', 'red', 'blue','green', 'red', 'blue', 'red', 'blue')
markers = ('*','o','x','*', '+', 'x', '*', '+','*','+','x','*', '+', 'x', '*', '+',)
lines_style = ('-','--','.-','--','-','--')

for k in range(0,K):
    hord = str(hords[k])
    for m in range(0,M):
        gtype = str(gtypes[m])
        dg    = dgs[m]
        mf    = str(mfs[m])
        adv   = advs[m]

        if dg==1:
          dg = 'dg1'
        elif dg==2:
          dg = 'dg2'
        else:
          dg = 'kinked'

        # convergence rate
        error = errors[:,k,m]
        n = len(Ns)-1
        CR = (np.log(error[n-1])-np.log(error[n]))/np.log(2.0)
        CR = str("{:2.1f}".format(CR))
        plt.loglog(Ns, error, lines_style[m], color=colors[m], marker=markers[m], \
            label = 'g'+str(gtype)+'.'+dg+'.adv'+str(adv)+'.hord'+hord+'.mf'+str(mf)+" - order "+CR)

# plot reference lines
eref = 50*np.amin(errors)
order1 = [eref, eref/2.0]
order2 = [eref, eref/4.0]
plt.loglog(Ns[ngrids - 2:ngrids], order1, '-',  color='black', label='1st order')
plt.loglog(Ns[ngrids - 2:ngrids], order2, '-.', color='black', label='2nd order')

# Label
plt.xlabel('$N$')
plt.ylabel('Error')
plt.xlim(0, 1000)
#plt.ylim(10**(-6), 10**(-2))
title = 'Divergence error'
plt.title(title)
plt.legend()
plt.grid(True, which="both")
filename = graphdir+'divergence_g'+str(gtype)+'.mf'+mf
plt.savefig(filename+'.'+figformat, format=figformat)

plt.close()
#plt.show()
