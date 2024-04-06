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

graphdir='../graphs/'
datadir='../SHiELD_SRC/test/CI/BATCH-CI/'
figformat='png'

#--------------------------------------------------------------------------------------------------------
# test case
tc  = -2

#grid type 0-equiedge; 2-equiangular
gtypes = (2, 2, 2, 2)
gtypes = (0, 0, 0, 0)
#gtypes = (2, 2, 2, 0, 0, 0)
#gtypes = (0, 0, 2, 2)
#gtypes = (2, 2)

# midpoints
mps=('s','c','s','c')

#duogrid scheme
#dgs = (1, 1, 2)
#dgs = (1, 1, 2, 1, 1, 2)
#dgs = (2, 2, 2, 2)
#dgs = (2, 2)
dgs = (1, 1, 2, 2)
# values of N
Ns = (48, 96, 192, 384, 768)
#Ns = (48, 96, 192, 384, )
#Ns = (48, 96, 192 )
ngrids = len(Ns)
#--------------------------------------------------------------------------------------------------------

basename='duogrid-test'
alpha = 45  # rotation angle
Tf = 12

# initialize error arrays and counter
M = len(gtypes)
error_linf = np.zeros((len(Ns),M,5))

#--------------------------------------------------------------------------------------------------------
# Loop over all schemes - to get errors
for m in range(0, M):
    gtype = gtypes[m]
    dg = dgs[m]
    mp = mps[m]

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

    n = 0
    # Now let us loop over all N values and compute the errors
    #------------------------------------------------------------------------------------------------
    for N in Ns:
        # Directory where the netcdf files are
        filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
        +".g"+str(gtype)+"."+mp+"."+dg+"/"
        
        file = filepath+"error_duogrid.txt"
        file = np.loadtxt(file)
        error_linf[n,m,:] = file
        # print errors
        print(m, N, error_linf[n,m,:])

        # update counter
        n = n + 1
#------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------
fields = ('ha','hb','u', 'uc', 'ud')
for f in range(0,len(fields)):
    # get max errors - for using same plotting scale
    emax = np.amax(error_linf[:,:,f])
    emax = 1.5*emax

    # get min errors
    emin = np.amin(error_linf[:,:,f])
    emin =  0.7*emin
    #print(emin, emax)

    #------------------------------------------------------------------------------------------------
    # Now let us plot the error graph
    colors = ('green', 'green', 'red', 'red', 'red', 'blue')
    markers = ('*','+','x','*','+','x')
    lines_style = ('-','--','-','--','--','--',)


    #-------------------------------------------------------------------------------------------------
    # Plot ONLY error graph
    dpi=100

    for m in range(0,M):
       gtype = gtypes[m]
       dg = dgs[m]
       mp = mps[m]

       # subtitle
       name = 'g'+str(gtype)+str(mp)+'.dg'+str(dg)

       # Get the maximum error and plot
       plt.title('$L_{\infty}$ error for '+fields[f] + ' - duogrid interpolation')
       plt.ylim(emin, emax)
       plt.xlim(40, 1000)
       error = error_linf[:,m,f]
       # convergence rate
       n = len(Ns)-1
       CR = (np.log(error[n-1])-np.log(error[n]))/np.log(2.0)
       CR = str("{:2.1f}".format(CR))
       plt.loglog(Ns, error, lines_style[m], color=colors[m], marker=markers[m], \
       label = name +" - order "+CR)

    # plot reference lines
    eref = emin*100
    order2 = [eref, eref/4.0]
    order4 = [eref, eref/16.0]
    plt.loglog(Ns[ngrids - 2:ngrids], order2, '-',  color='black', label='2nd order')
    plt.loglog(Ns[ngrids - 2:ngrids], order4, '-.', color='black', label='4th order')

    plt.xlabel('$N$')
    plt.ylabel('Error')
    plt.ylim(10**(-12),10**(-3.5))
    #plt.title(title)
    plt.legend()
    plt.grid(True, which="both")
    title = "Duogrid interpolation - degree = 3"
    filename = graphdir+"error_tc"+str(tc)+"_"+fields[f]
    plt.savefig(filename+'.'+figformat, format=figformat)
    plt.close()
#------------------------------------------------------------------------------------------------
