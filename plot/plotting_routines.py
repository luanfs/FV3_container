#----------------------------------------------------------------------------------------------
# Module for plotting routines
#----------------------------------------------------------------------------------------------
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors, colorbar
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import math
graphdir = '../graphs/'
#----------------------------------------------------------------------------------------------
# Plot the scalar field q given in a A grid.
#----------------------------------------------------------------------------------------------
def plot_scalarfield(q, map_projection, title, filename, filepath, colormap, qmin, qmax,\
    dpi, figformat):
    # map projection
    if map_projection == "mercator":
        plt.figure(figsize=(1832/dpi, 977/dpi), dpi=dpi)
        #plt.figure(figsize=(1000/dpi, 1000/dpi), dpi=dpi)
        plateCr = ccrs.PlateCarree()        
    elif map_projection == "sphere":
        plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi) 
        plateCr = ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0)
        plateCr = ccrs.Orthographic(central_longitude=0.25*180, central_latitude=180/6.0)

    projection=ccrs.PlateCarree(central_longitude=0)
    plateCr._threshold = plateCr._threshold/10.

    ax = plt.axes(projection=plateCr)
    ax.set_global()
    ax.stock_img()
    ax.gridlines(draw_labels=True)

    # Color of each cubed panel
    colors = ('black','black','black','black','black','black')
    N = np.shape(q)[0]
    # plot for each tile
    for tile in range(0,6):
        # Grid file to be opened
        grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

        # Load the file
        grid = xr.open_dataset(grid_file , decode_times=False, chunks={'time': 5}).squeeze()

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

        #ax.plot([A_lon, B_lon], [A_lat, B_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        #ax.plot([B_lon, C_lon], [B_lat, C_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        #ax.plot([C_lon, D_lon], [C_lat, D_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        #ax.plot([D_lon, A_lon], [D_lat, A_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)

        # Plot scalar field
        im = ax.pcolormesh(lon, lat, q[:,:,tile], alpha=1, transform=ccrs.PlateCarree(), \
        zorder=10, vmin = qmin, vmax=qmax,  cmap=colormap)
    #plt.xlim(45-25,45+25)
    #plt.ylim(35-25,45+25)
    plt.title(title)
    # Plot colorbar
    cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.046, pad=0.04, shrink=0.8, format='%.1e')
    cb=plt.colorbar(im, cax=cax, extend='both',**kw)
    cb.ax.tick_params(labelsize=8)
    plt.savefig(filename+'.'+figformat, format=figformat)
    plt.close()
#-----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
# Plot the scalar field q given in a A grid in a panel for all time steps
#----------------------------------------------------------------------------------------------
def plot_scalarfield_panel(q, map_projection, filename, filepath, title, colormap, qmin, qmax,\
    nplots, dt, dpi, figformat):

    # panel dimensions
    ncols = 2
    nlines = math.ceil((nplots+1) / ncols)
    fig = plt.figure(figsize=(nlines*500/dpi, ncols*4000/dpi),dpi=dpi)
    gs = fig.add_gridspec(nlines, ncols, width_ratios=[1] * ncols)
    k = 0
    time = 0
    for i in range(0, nlines):
        for j in range(0, ncols):
            if k<=nplots:
                ax = fig.add_subplot(gs[i, j], projection=ccrs.PlateCarree())  # Use PlateCarree projectioni
                #ax.set_global()
                plateCr = ccrs.PlateCarree()        
                projection=ccrs.PlateCarree(central_longitude=0)
                plateCr._threshold = plateCr._threshold/10.
                #ax.stock_img()
                ax.gridlines(draw_labels=True)
                # plot for each tile
                for tile in range(0,6):
                    # Grid file to be opened
                    grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

                    # Load the file
                    grid = xr.open_dataset(grid_file , decode_times=False, chunks={'time': 5}).squeeze()

                    # Get grid
                    lon = grid['grid_lon']
                    lat = grid['grid_lat']

                    # Plot scalar field
                    im = ax.pcolormesh(lon, lat, q[:,:,tile,k],alpha=1,transform=ccrs.PlateCarree(),\
                    zorder=10, vmin = qmin, vmax=qmax,  cmap=colormap)

                dmin = str("{:.1e}".format(np.amin(abs(q[:,:,:,k]))))
                dmax = str("{:.1e}".format(np.amax(abs(q[:,:,:,k]))))
                Time = str("{:.1e}".format(time))
                if np.amin(abs(q[:,:,:,k]))>0.0000000000001:
                    subtitle =str(Time)+" days, min = "+ dmin+", max = "+dmax
                else:
                    subtitle =str(Time)+" days, max = "+dmax

                # Add coastlines to the subfigure
                #ax.add_feature(ccrs.COASTLINE, edgecolor='k', linewidth=0.8)
                #ax.add_feature(cfeature.COASTLINE, edgecolor='k', linewidth=0.8)
     
                # Add a subtitle to the subfigure
                ax.set_title(subtitle)

                # time update
                time = time + dt

            #else:
            #    print('bye!')
            k = k+1

    #plt.tight_layout()  # Ensure proper spacing between subplotsp

    # Add a single colorbar
    #cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.5])  # Define the position and size of the vertical colorbar
    #cbar = fig.colorbar(im, cax=cbar_ax, fraction=0.046, pad=0.04, shrink=0.8, format='%.1e')
    #cbar.ax.tick_params(labelsize=8)

    #fig.suptitle(title+'\n\n', fontsize=30)

    plt.savefig(filename+'.'+figformat, format=figformat)
    #plt.show()
    plt.close()
    #exit()

    #lt.title(title)
#-----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------
# Plot the scalar field q given in a A grid.
#----------------------------------------------------------------------------------------------
def plot_scalarfield2(q, lons, lats, map_projection, title, filename, colormap, dpi, figformat, qmin, qmax):
    # map projection
    if map_projection == "mercator":
        plt.figure(figsize=(1832/dpi, 977/dpi), dpi=dpi)
        plateCr = ccrs.PlateCarree()        
    elif map_projection == "sphere":
        plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi) 
        plateCr = ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0)
       # plateCr = ccrs.Orthographic(central_longitude=0.25*180, central_latitude=180/6.0)

    projection=ccrs.PlateCarree(central_longitude=0)
    plateCr._threshold = plateCr._threshold/10.

    ax = plt.axes(projection=plateCr)
    #ax.coastlines()
    ax.set_global()
    ax.stock_img()
    #ax.gridlines(draw_labels=True)
    # Color of each cubed panel
    colors = ('black','black','black','black','black','black')

    N = np.shape(q)[0]
    # plot for each tile
    for tile in range(0,6):
        # Get grid
        lon = lons[:,:,tile]
        lat = lats[:,:,tile]

        # Plot cube edges
        A_lon, A_lat = lon[0,0], lat[0,0]
        B_lon, B_lat = lon[N, 0], lat[N, 0]
        C_lon, C_lat = lon[N, N], lat[N, N]
        D_lon, D_lat = lon[0, N], lat[0, N]
        lw = 0.2
        plt.rcParams["axes.axisbelow"] = True

        #ax.plot([A_lon, B_lon], [A_lat, B_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        #ax.plot([B_lon, C_lon], [B_lat, C_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        #ax.plot([C_lon, D_lon], [C_lat, D_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        #ax.plot([D_lon, A_lon], [D_lat, A_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)

        # Plot scalar field
        im = ax.pcolormesh(lon, lat, q[:,:,tile], alpha=1, transform=ccrs.PlateCarree(), \
        zorder=10, vmin = qmin, vmax=qmax,  cmap=colormap)
 

    plt.title(title)
    # Plot colorbar
    cax,kw = colorbar.make_axes(ax,orientation='vertical' , fraction=0.046, pad=0.04, shrink=0.8, format='%.1e')
    cb=plt.colorbar(im, cax=cax, extend='both',**kw)
    cb.ax.tick_params(labelsize=8)
    plt.savefig(graphdir+filename+'.'+figformat, format=figformat)
    plt.close()
#-----------------------------------------------------------------------------------------
