import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors, colorbar


#----------------------------------------------------------------------------------------------
# directories
datadir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/solo_sw/'
graphdir = '/gpfs/f5/scratch/Luan.Santos/gfdl_w/graphs_solo_sw/'
#----------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------
#simulation parameters
N=96
Tf=100
gtype = 0
hords = (8,)
advs  = (1,2)
testname='modon'
alpha=0
dg='dg1'
basename= "C"+str(N)+".sw."+testname+".alpha"+str(alpha)+".g"+str(gtype)+"."+str(dg)
#----------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------
# Plot the scalar field q given in a A grid.
#----------------------------------------------------------------------------------------------
def plot_scalarfield(q, title, filename, filepath, colormap, qmin, qmax):
    # Figure quality
    dpi = 100

    # Figure format
    figformat='png'

    # Map projection
    map_projection = "mercator"
    #map_projection = "sphere"

    # map projection
    if map_projection == "mercator":
        plt.figure(figsize=(1600/dpi, 1000/dpi), dpi=dpi)
        plateCr = ccrs.PlateCarree()
    elif map_projection == "sphere":
        plt.figure(figsize=(1000/dpi,1000/dpi), dpi=dpi) 
        plateCr = ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0)

    projection=ccrs.PlateCarree(central_longitude=0)
    plateCr._threshold = plateCr._threshold/10.

    ax = plt.axes(projection=plateCr)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 19, 'color': 'black'}
    gl.ylabel_style = {'size': 19, 'color': 'black'}

 
    # Color of each cubed panel
    colors = ('black','black','black','black','black','black')
    N = np.shape(q)[0]
    # plot for each tile
    for tile in range(0,6):
        # Grid file to be opened
        grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

        # Load the file
        grid = xr.open_dataset(grid_file , decode_times=False)

        # Get grid
        lon = grid['grid_lon']
        lat = grid['grid_lat']

        # Plot cube edges
        A_lon, A_lat = lon[0,0], lat[0,0]
        B_lon, B_lat = lon[N, 0], lat[N, 0]
        C_lon, C_lat = lon[N, N], lat[N, N]
        D_lon, D_lat = lon[0, N], lat[0, N]
        lw = 0.2
        plt.rcParams["axes.axisbelow"] = True

        ax.plot([A_lon, B_lon], [A_lat, B_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([B_lon, C_lon], [B_lat, C_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([C_lon, D_lon], [C_lat, D_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([D_lon, A_lon], [D_lat, A_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)

        # Plot scalar field
        im = ax.pcolormesh(lon, lat, q[:,:,tile], alpha=1, transform=ccrs.PlateCarree(), \
        zorder=10, vmin = qmin, vmax=qmax,  cmap=colormap)

    plt.title(title,  fontsize=19)
    #print(qmin,qmax)
    # Plot colorbar
    cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.046, pad=0.04, format='%.1e')
    #cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.09, pad=0.04, shrink=0.9, format='%.0e')

    cb=plt.colorbar(im, cax=cax, extend='both',**kw)

    ticks = np.linspace(qmin, qmax, num=5)
    cb.set_ticks(ticks)
    cb.ax.tick_params(labelsize=22)
    plt.savefig(filename+'.'+figformat, format=figformat)
    plt.close()
    #exit()
#-----------------------------------------------------------------------------------------


# grid name
if gtype==0:
    gname = 'equiedge'
elif gtype==2:
    gname = 'equiangular'


# Create filepaths list
filepaths = []
for hord in hords: 
    datas = []
    for adv in advs:
        filepath = datadir+basename+'.adv'+str(adv)+'.hord'+str(hord)+"/rundir/"
        filepaths.append(filepath)


# Get the number of plots
atmos_file = filepaths[0]+"atmos_daily.tile1.nc"
data = xr.open_dataset(atmos_file, decode_times=False)
times = data.time.values
nplots = len(times)
dtplot = times[1]-times[0]

#-----------------------------------------------------------------------------------------
# Arrays to store the data that will be plotted
# FV3 data
h = np.zeros((N,N,6,nplots+1,len(advs)))
u = np.zeros((N,N,6,nplots+1,len(advs)))
v = np.zeros((N,N,6,nplots+1,len(advs)))
vort = np.zeros((N,N,6,nplots+1,len(advs)))



# This loop over tiles computes the maximum errors
time = 0  
tgap=5
#nplots = 100
for t in range(0,nplots+1,tgap):
    for k, filepath in enumerate(filepaths):
        #print(time, k, filepath)
        for tile in range(0,6):
            # Files to be opened
            atmos_file = filepath+"atmos_daily.tile"+str(tile+1)+".nc"
            grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

            # Load the files
            data = xr.open_dataset(atmos_file, decode_times=False)
            grid = xr.open_dataset(grid_file , decode_times=False)

            # Variable to be plotted (ps = fluid depth in the SW model)
            if t>=1:
               h[:,:,tile,t,k] = data['ps'][t-1,:,:].values
               u[:,:,tile,t,k] = data['ucomp'][t-1,:,:].values
               v[:,:,tile,t,k] = data['vcomp'][t-1,:,:].values
               vort[:,:,tile,t,k] = data['vort'][t-1,:,:].values
            else: #get ic
               h[:,:,tile,t,k] = data['ps_ic'][:,:].values
               u[:,:,tile,t,k] = data['ua_ic'][:,:].values
               v[:,:,tile,t,k] = data['va_ic'][:,:].values


    # time update
    time = time + dtplot 

# plot
time = 0  
for t in range(0,nplots+1,tgap):
    ##############################################################################################
    fref = 5000
    field_adv1 = h[:,:,:,t,0]-fref
    field_adv2 = h[:,:,:,t,1]-fref

    fmin = min(np.amin(field_adv1), np.amin(field_adv2))
    fmax = max(np.amax(field_adv1), np.amax(field_adv2))
    fabs = max(abs(fmin),abs(fmax))
    fmin, fmax = -fabs, fabs

    adv=1
    field='h'
    title = field+'_'+testname+'_t'+str(t)+'_'+field+'_'+basename+'.adv'+str(adv)+'.hord'+str(hord)
    filename = graphdir+title
    plot_scalarfield(field_adv1, title, filename, filepaths[0], 'seismic', fmin, fmax)

    adv=2
    field='h'
    title = field+'_'+testname+'_t'+str(t)+'_'+field+'_'+basename+'.adv'+str(adv)+'.hord'+str(hord)
    filename = graphdir+title
    plot_scalarfield(field_adv2, title, filename, filepaths[0], 'seismic', fmin, fmax)

    diff = (field_adv1-field_adv2)/(field_adv1+fref)
    diffabs = max(abs(np.amin(diff)), np.amax(diff))
    #plot_scalarfield(diff, title, filename, filepaths[0], 'seismic', -diffabs, diffabs)
    ##############################################################################################

    ##############################################################################################
    fref = 0
    field_adv1 = vort[:,:,:,t,0]-fref
    field_adv2 = vort[:,:,:,t,1]-fref

    fmin = min(np.amin(field_adv1), np.amin(field_adv2))
    fmax = max(np.amax(field_adv1), np.amax(field_adv2))
    fabs = max(abs(fmin),abs(fmax))
    fmin, fmax = -fabs, fabs


    adv=1
    field='vort'
    title = field+'_'+testname+'_t'+str(t)+'_'+field+'_'+basename+'.adv'+str(adv)+'.hord'+str(hord)
    filename = graphdir+title
    plot_scalarfield(field_adv1, title, filename, filepaths[0], 'seismic', fmin, fmax)

    adv=2
    field='vort'
    title = field+'_'+testname+'_t'+str(t)+'_'+field+'_'+basename+'.adv'+str(adv)+'.hord'+str(hord)
    filename = graphdir+title
    plot_scalarfield(field_adv2, title, filename, filepaths[0], 'seismic', fmin, fmax)
    ##############################################################################################


    # time update
    time = time + dtplot 

    print(t, diffabs)

