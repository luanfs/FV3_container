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
from reference_solution import href_Agrid

graphdir='../graphs/'
datadir='../SHiELD_SRC/test/CI/BATCH-CI/'
figformat='png'

#--------------------------------------------------------------------------------------------------------
# test case
tc  = 6

# advection scheme
hord = 8

#grid type 0-equiedge; 2-equiangular
#gtypes = (0, 0, 2, 2)
#gtypes = (0, 0, 2, 2)
#gtypes = (2, 2, 2, 2)
gtypes = (2, 2, 2)
gtypes = (2, 2)

#duogrid
#dgs = (0, 1, 1, 2)
#dgs = (1, 1, 2)
dgs = (1, 1, 2)
dgs = (1, 2, 2)

#midpoints flag
#mps = (False, False, True, True)
mps = (False, False, False)
#mps = (False, True, True)
#mps = (True, True)

#metric tensor
#mts = (2, 2, 2, 2)
mts = (1, 1, 1, 1)
mts = (1, 1, 2)

#departure point
dps = (1, 1, 1, 1)
#dps = (2, 2, 2, 2)
dps = (1, 1, 2)

#divergence damping (only for sw)
if tc == 2:
   dds = (0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12)
else:
   dds = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# contravariant wind init? (only for tc<=1)
if tc <=1:
   #cws = ('.false.','.true.','.false.','.true.')
   cws = ('.false.','.false.','.false.','.false.')
   #cws = ('.true.','.true.','.true.','.true.')

# values of N
#Ns = (48, 96, 192, 384, 768)
#Ns = (48, 96, 192, 384)
#Ns = (192, )
#Ns = (384, )
#Ns = (48, 96, )
Ns = (48, )
#Ns = (96, )
#Ns = (768,)
#--------------------------------------------------------------------------------------------------------


# basename used in outputs
if tc==1 :
    basename='cosine-zonal'
    hmin, hmax = 0.5, 1.6
    alpha = 45 # rotation angle
    Tf = 12

elif tc==2 :
    basename='geobalance'
    hmin, hmax = 1000.0, 3000.0
    alpha = 45 # rotation angle
    #alpha = 0 # rotation angle
    Tf = 5

elif tc==-3:
    basename='constant'
    hmin, hmax = 999.0, 1000.1
    alpha = 45 # rotation angle
    Tf = 12

elif tc==-4:
    basename='gaussian-zonal'
    hmin, hmax = 980, 1120 # for Gaussian hill
    alpha = 45 # rotation angle
    Tf = 12

elif tc==-5:
    basename='geobalance'
    hmin, hmax = 1000.0, 3000.0
    alpha = 45 # rotation angle
    Tf = 12

elif tc==-6:
    basename='gaussian-ndiv'
    hmin, hmax = 0.05, 1.05 # for Gaussian hills
    alpha = 0 # rotation angle
    Tf = 12

elif tc==-7:
    basename='gaussian-div'
    hmin, hmax = 0.05, 1.05 # for Gaussian hills
    alpha = 0 # rotation angle
    Tf = 12

else:
   print('ERROR: invalid initial condition')
   exit()


#-----------------------------------------------------------------------------------------
nplots = int(Tf)
dtplot = Tf/nplots
if tc==-3 or tc==1 or tc==2 or tc==-4 or tc==-5: # reference solution is given for all timesteps
    ts, te = 0, nplots
    time = 0
else: # reference solution only at final time
    ts, te = nplots, nplots
    time = Tf


# initialize error arrays and counter
M = len(gtypes)


error_linf_max= np.zeros((M))
#--------------------------------------------------------------------------------------
# Now let us loop over all N values and compute the errors
#------------------------------------------------------------------------------------------------
for N in Ns:
    # Loop over all schemes - to get errors
    error_linf_max[:]=0

    # Arrays to store the data that will be plotted
    h_fv3   = np.zeros((N,N,6,M,nplots+1))
    h_ref   = np.zeros((N,N,6,M,nplots+1))
    h_error = np.zeros((N,N,6,M,nplots+1))
 
    # grid plotting
    lats = []
    lons = []

    for m in range(0, M):
        gtype = gtypes[m]
        dg = dgs[m]
        mp = mps[m]
        mt = str(mts[m])
        dp = str(dps[m])
        dd = str(dds[m])

        # grid name
        if gtype==0:
            gname = 'equiedge'
        elif gtype==2:
            gname = 'equiangular'

        # duogrid name
        if dg==1:
           dg = 'duogrid1'
        elif dg==2:
           dg = 'duogrid2'
        else:
           dg = 'kinked'

        # midpoints name
        if mp:
            mp = 'mp1'
        else:
            mp = 'mp0'

        if tc<=1:
            cw = cws[m]

        if tc>1:
       	  filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
          +".g"+str(gtype)+"."+dg+"."+mp+".mt"+mt+".dp"+dp+".hord"+str(hord)+'.dd'+str(dd)+".tf"+str(Tf)+"/"
        else:
       	  filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
          +".g"+str(gtype)+"."+dg+"."+mp+".mt"+mt+".dp"+dp+".hord"+str(hord)+'.cw'+cw+".tf"+str(Tf)+"/"
        print(filepath)
        # Loop over time
        for t in range(ts,te+1):
            # Loop over tiles
            for tile in range(0,6):
                # Files to be opened
                atmos_file = filepath+"atmos_daily.tile"+str(tile+1)+".nc"
                grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

                # Check if they exist
                files_exist = os.path.exists(grid_file) and os.path.exists(atmos_file)
                if files_exist:
                    # Load the files
                    data = xr.open_dataset(atmos_file, decode_times=False, chunks={'time': 5}).squeeze()
                    grid = xr.open_dataset(grid_file , decode_times=False, chunks={'time': 5}).squeeze()
                    areas = grid.area

                    # Variable to be plotted (ps = fluid depth in the SW model)
                    if t>=1:
                        h_fv3[:,:,tile,m,t] = data['ps'][t-1,:,:].values
                    else: #get ic
                        h_fv3[:,:,tile,m,t] = data['ps_ic'][:,:].values

                    # Get reference solution
                    if tc==-4 or tc==1:
                        if t==0 or t==nplots:
                            h_ref[:,:,tile,m,t] = data['ps_ic'][:,:].values
                        else:
                            h_ref[:,:,tile,m,t] = href_Agrid(grid, t, dtplot, tc, alpha)

                    if tc==-6 or tc==-7:
                        if t==0 or t==nplots:
                            h_ref[:,:,tile,m,t] = data['ps_ic'][:,:].values
         
                    elif tc==-3 or tc==2 or tc==-5:
                        h_ref[:,:,tile,m,t] = data['ps_ic'][:,:].values

                    # Compute all errors
                    h_error[:,:,tile,m,t]= h_ref[:,:,tile,m,t] - h_fv3[:,:,tile,m,t]

                    if len(lats)<6:
                        lats.append(grid['grid_lat'])
                        lons.append(grid['grid_lon'])
                else:
                    h_error[:,:,tile,m,t]= 0.0
                    h_fv3[:,:,tile,m,t]= 0.0
                    #print("error ",m,t)

        #print(m,N,t,np.amax(h_error[:,:,:,m,t]))
        # print errors
        #print(m, N, gname, dg)
   

    # Plot range
    emax = max(abs(np.amin(h_error)), abs(np.amax(h_error)))
    emin = -emax
    #print(emin,emax)
    # Loop over time
    #time = 0
    time = 0
    for t in range(ts,te+1):
        Time = str("{:.2e}".format(time))
        print(Time)
        #-------------------------------------------------------------------------------------------------
        # Plot ONLY graph of error
        # panel dimensions
        ncols = 1
        nlines = math.ceil((M) / ncols)
        dpi=50
        fig = plt.figure(figsize=(nlines*500/dpi, ncols*1000/dpi),dpi=dpi)
        gs = fig.add_gridspec(nlines, ncols, width_ratios=[1] * ncols)
        m = -1
        # Loop over all schemes
        for i in range(0, nlines):
            for j in range(0, ncols):
                m = m+1
                if(m < M):
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
                       dg = 'duogrid1'
                    elif dg==2:
                       dg = 'duogrid2'
                    else:
                       dg = 'kinked'

                    # midpoints name
                    if mp:
                       mp = 'mp1'
                    else:
                       mp = 'mp0'

                    dmax = str("{:.2e}".format(np.amax(abs(h_error[:,:,:,m,t]))))

                    #filepath = "../BATCH-CI/C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
                    #+".g"+str(gtype)+"."+dg+".tf"+str(Tf)+"/"

                    ax = fig.add_subplot(gs[i, j], projection=ccrs.PlateCarree())  # Use PlateCarree projectioni
                    #ax.set_global()
                    plateCr = ccrs.PlateCarree()        
                    projection=ccrs.PlateCarree(central_longitude=0)
                    plateCr._threshold = plateCr._threshold/10.
                    #ax.stock_img()
                    ax.gridlines(draw_labels=True)
                    # plot for each tile
                    for tile in range(0,6):
                        # Get grid
                        lon = lons[tile]
                        lat = lats[tile]

                        # Plot cube edges
                        A_lon, A_lat = lon[0,0], lat[0,0]
                        B_lon, B_lat = lon[N, 0], lat[N, 0]
                        C_lon, C_lat = lon[N, N], lat[N, N]
                        D_lon, D_lat = lon[0, N], lat[0, N]
                        lw = 0.2
                        plt.rcParams["axes.axisbelow"] = True

                        ax.plot([A_lon, B_lon], [A_lat, B_lat], '-', linewidth=lw, color='black', transform=ccrs.Geodetic(), zorder=11)
                        ax.plot([B_lon, C_lon], [B_lat, C_lat], '-', linewidth=lw, color='black', transform=ccrs.Geodetic(), zorder=11)
                        ax.plot([C_lon, D_lon], [C_lat, D_lat], '-', linewidth=lw, color='black', transform=ccrs.Geodetic(), zorder=11)
                        ax.plot([D_lon, A_lon], [D_lat, A_lat], '-', linewidth=lw, color='black', transform=ccrs.Geodetic(), zorder=11)

                        # Plot scalar field
                        im = ax.pcolormesh(lon, lat, h_error[:,:,tile,m,t],alpha=1,transform=ccrs.PlateCarree(),\
                        zorder=10, vmin=emin, vmax=emax,cmap='seismic')

                    # subtitle
                    subtitle = gname+' - '+dg+' - '+mp+' - max='+dmax
                    ax.set_title(subtitle, fontsize=20)


        # Set a common title
        if tc==-3 or tc==1 or tc==2  or tc==-6:
            title = "TC"+str(tc)+" - $\\alpha$ = "+str(alpha)+", C"+str(N) +', '+Time+' days'
        else:
            title = "TC"+str(tc)+", C"+str(N)+', '+Time+' days'


        # Plot colorbar
        # Add a single colorbar
        cbar_ax = fig.add_axes([0.92, 0.2, 0.02, 0.7])  # Define the position and size of the vertical colorbar
        cbar = fig.colorbar(im, cax=cbar_ax, fraction=0.046, pad=0.04, shrink=0.8, format='%.1e')
        cbar.ax.tick_params(labelsize=12)
        fig.suptitle(title+'\n\n', fontsize=30)
        #plt.tight_layout()  # Ensure proper spacing between subplotsp
        filename = graphdir+"cs_"+str(N)+"error_tc"+str(tc)+'_alpha'+str(alpha)+"_t"+str(t)
        plt.savefig(filename+'.'+figformat, format=figformat)
        plt.close()

        # time update
        time = time + dtplot


        # Loop over all schemes - to print errors
#------------------------------------------------------------------------------------------------
