import xarray as xr
import numpy as np
from plotting_routines  import plot_scalarfield, plot_scalarfield_panel
from reference_solution import href_Agrid
#-------------------------------------------------------------------
# Plotting parameters
# Figure quality
dpi = 100

# Map projection
map_projection = "mercator"
#map_projection = "sphere"

graphdir='../graphs/'
datadir='../SHiELD_SRC/test/CI/BATCH-CI/'
figformat='png'
#figformat='pdf'
#-------------------------------------------------------------------


#-------------------------------------------------------------------
# test parameters
tc = -7

# 1dadvection scheme
hord = 8

#grid type
gtype = 2  # 0-equiedge; 2-equiangular

# duogrid parameter
dg = 2

# 2d advection scheme
adv = 2

# mass fixer
mf = 1

# values of N
#Ns = (48,)
#Ns = (96,)
Ns=(192,)
#Ns = (384,)


# basename used in outputs
if tc==1 :
    basename='cosine-zonal'
    alpha = 45 # rotation angle
    Tf = 12
    plot_error = True
    hmin, hmax = 0.5, 1.5

elif tc==-3 :
    basename='gaussian-zonal'
    alpha = 45 # rotation angle
    Tf = 12
    plot_error = True
    hmin, hmax = 0, 1.0

elif tc==-4 :
    basename='geobalance'
    alpha = 45 # rotation angle
    Tf = 12
    plot_error = True
    hmin, hmax = 1000.0, 3000.0
    
elif tc==-5 :
    basename='twogaussians-ndiv'
    alpha = 0 # rotation angle
    Tf = 12
    plot_error = True
    hmin, hmax = 0.0, 1.0
    
elif tc==-6 :
    basename='twogaussians-div'
    alpha = 0 # rotation angle
    Tf = 12
    plot_error = True
    hmin, hmax = 0.0, 5.5

elif tc==-7 :
    basename='cylinder-zonal'
    alpha = 45 # rotation angle
    Tf = 12
    plot_error = True
    hmin, hmax = 0.0, 1.0

elif tc==-8 :
    basename='cylinder-div'
    alpha = 0 # rotation angle
    Tf = 12
    plot_error = True
    hmin, hmax = 0.0, 1.0

else:
    print('ERROR: invalid initial condition')
    exit()

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

# adv name
if adv==1:
   advname = 'PL'
elif adv==2:
   advname = 'LT'


#-------------------------------------------------------------------
# Loop over grids
for N in Ns:
    # Directory where the netcdf files are
    filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
    +".g"+str(gtype)+"."+dg+".adv"+str(adv)+".hord"+str(hord)+".mf"+str(mf)+".tf"+str(Tf)+"/"

    # Get the number of plots
    atmos_file = filepath+"atmos_daily.tile1.nc"
    data = xr.open_dataset(atmos_file, decode_times=False, chunks={'time': 1}).squeeze()
    nplots = int(len(data.time.values))
    dtplot = Tf/nplots

    # Arrays to store the data that will be plotted
    h_fv3   = np.zeros((N,N,6,nplots+1))
    h_ref   = np.zeros((N,N,6,nplots+1))
    h_error = np.zeros((N,N,6,nplots+1))

    #-----------------------------------------------------------------------------------------
    # This loop over tiles computes the maximum errors
    for t in range(0,nplots+1):        
        for tile in range(0,6):
            # Files to be opened
            atmos_file = filepath+"atmos_daily.tile"+str(tile+1)+".nc"
            grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

            # Load the files
            data = xr.open_dataset(atmos_file, decode_times=False, chunks={'time': 5}).squeeze()
            grid = xr.open_dataset(grid_file , decode_times=False, chunks={'time': 5}).squeeze()

            # Variable to be plotted
            if t==0:
               h_fv3[:,:,tile,t] = data['ps_ic'][:,:].values
            else:
               h_fv3[:,:,tile,t] = data['ps'][t-1,:,:].values
            #print(np.amin(h_fv3[:,:,tile,t]),np.amax(h_fv3[:,:,tile,t]))
            # Get reference solution
            if tc==1 or tc==-3:
                if t==0 or t==nplots:
                    h_ref[:,:,tile,t] = data['ps_ic'][:,:].values
                else:
                    h_ref[:,:,tile,t] = href_Agrid(grid, t, dtplot, tc, alpha)

            if tc==-5 or tc==-6 or tc==-7:
                if t==0 or t==nplots:
                    h_ref[:,:,tile,t] = data['ps_ic'][:,:].values
 
            elif  tc==-4:
                h_ref[:,:,tile,t] = data['ps_ic'][:,:].values

    h_error = h_ref - h_fv3
    h_error_min, h_error_max = np.amin(h_error[:,:,:,nplots]), np.amax(h_error[:,:,:,nplots])
    h_error_abs = max(abs(h_error_min), abs(h_error_max))
    h_error_abs = 3*10**(-3)
    h_error_abs = 1.1*10**(-2)
    hmin_fv3, hmax_fv3 = hmin, hmax
    hmin_fv3, hmax_fv3 =  np.amin(h_fv3), np.amax(h_fv3)
    hmin_fv3, hmax_fv3 =  0.0, 1.0
    #-----------------------------------------------------------------------------------------

    #-----------------------------------------------------------------------------------------
    if plot_error:
       # Plot the h error
       if tc==1 or tc==-3 or tc==-4:
          ts, te = 0, nplots
          time = 0.0
       else: # reference solution only at final time
          ts, te = nplots, nplots
          time = Tf

       for t in range(ts,te+1):
          dmax = str("{:.2e}".format(np.amax(abs(h_error[:,:,:,t]))))
          Time = str("{:.2e}".format(time))

          if tc==1 or tc==-3 or tc==-4:
             title = "Fluid depth $h$ error - Test case "+str(tc)+", $\\alpha$ = "+str(alpha)\
             +", time = "+str(Time)+" days, max = "+dmax+"\n"\
            +"Grid = "+gname+", "+advname+".hord"+str(hord)+".mf"+str(mf)+", N = "+str(N)\
             +"\n"
          else:
             title = "Fluid depth $h$ error - Test case "+str(tc)\
             +", time = "+str(Time)+" days, max = "+dmax+"\n"\
            +"Grid = "+gname+", "+advname+".hord"+str(hord)+".mf"+str(mf)+", N = "+str(N)\
             +"\n"

          # Filename
          filename = graphdir+"h_error_tc"+str(tc)+'_t'+str(t)+'_alpha'+str(alpha)+'_C'+str(N)\
          +"_g"+str(gtype)+"_"+dg+"_adv"+str(adv)+"_hord"+str(hord)+"_mf"+str(mf)+"_tf"+str(Tf)

          # plot
          if t==te:
            plot_scalarfield(h_error[:,:,:,t], map_projection, title, filename, filepath, \
             'seismic', -h_error_abs, h_error_abs, dpi, figformat)
          print('plotted the file '+filename)

          # time update
          time = time + dtplot

    #-----------------------------------------------------------------------------------------

    #exit()

    #-----------------------------------------------------------------------------------------
    # Plot the scalar field
    #-----------------------------------------------------------------------------------------
    dtplot = Tf/nplots
    time = 0
    gap = 1
    for t in range(0,nplots+1,gap):
        # finish the plotting
        # get min/max
        hminfv3, hmaxfv3 = np.amin(h_fv3[:,:,:,t]), np.amax(h_fv3[:,:,:,t])
        dmin = str("{:.2e}".format(hminfv3))
        dmax = str("{:.2e}".format(hmaxfv3))
        Time = str("{:.2e}".format(time))

        # title
        if tc==1 or tc==-3 or tc==-4:
            title = "Fluid depth $h$ - Test case "+str(tc)+", $\\alpha$ = "+str(alpha)\
            +", time = "+Time+" days"+", min = "+dmin+", max = "+dmax+"\n"\
            +"Grid = "+gname+", "+advname+".hord"+str(hord)+".mf"+str(mf)+", N = "+str(N)\
            +"\n"
        else:
            title = "Fluid depth $h$ - Test case "+str(tc)\
            +", time = "+Time+" days"+", min = "+dmin+", max = "+dmax+"\n"\
            +"Grid = "+gname+", "+advname+".hord"+str(hord)+".mf"+str(mf)+", N = "+str(N)\
            +"\n"


        # file
        filename = graphdir+"h_tc"+str(tc)+'_t'+str(t)+'_alpha'+str(alpha)+'_C'+str(N)\
          +"_g"+str(gtype)+"_"+dg+"_adv"+str(adv)+"_hord"+str(hord)+"_mf"+str(mf)+"_tf"+str(Tf)

        # plot
        plot_scalarfield(h_fv3[:,:,:,t], map_projection, title, filename, filepath, \
        'jet', hmin_fv3, hmax_fv3, dpi, figformat)
        print('plotted the file '+filename)

        # time update
        time = time + gap*dtplot

    #-----------------------------------------------------------------------------------------