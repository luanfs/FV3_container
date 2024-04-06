import xarray as xr
import numpy as np
from plotting_routines  import plot_scalarfield, plot_scalarfield_panel
from reference_solution import topography
#-------------------------------------------------------------------
# Plotting parameters
# Figure quality
dpi = 100

# Map projection
map_projection = "mercator"
#map_projection = "sphere"

graphdir='../graphs/'
datadir='../SHiELD_SRC/test/CI/BATCH-CI/'
endgamedir='../SHiELD_SRC/test/endgame/dump/'
figformat='png'
sec2day = 86400
#figformat='pdf'
#-------------------------------------------------------------------


#-------------------------------------------------------------------
# test parameters
tc = 2

# advection scheme
hord = 8

gtype = 2 # 0-equiedge; 2-equiangular
#gtype = 0

# duogrid parameter
dg = 1
#dg = 2

# advection scheme
adv=1
#adv=2

#mass fixer
mf = 1

#divergence damping (only for sw)
if tc==2:
  dd = 0.12
else:
  dd = 0.12
  #dd = 0

# values of N
Ns = (48,)
#Ns = (96,)
Ns=(192,)
#Ns = (384,)

# basename used in outputs
if tc==2:
    basename='geobalance'
    alpha = 45 # rotation angle
    Tf = 5
    hmin, hmax = 1000.0, 3000.0
    umin, umax = -10, 40
    vmin, vmax = -25, 25
    pvmin, pvmax = -3.1e-08, 3.1e-08
    vortmin, vortmax = -3.6e-05, 4.7e-05

elif tc==5:
    basename='mountain'
    alpha = 0 # rotation angle
    Tf = 14
    hmin, hmax = 5000.0, 6000.0
    umin, umax = -10, 40
    vmin, vmax = -25, 25
    pvmin, pvmax = -3.1e-08, 3.1e-08
    vortmin, vortmax = -3.6e-05, 4.7e-05

elif tc==6:
    basename='rhwave'
    alpha = 0 # rotation angle
    Tf = 105
    hmin, hmax = 8000.0, 10500.0
    umin, umax = -70, 100
    vmin, vmax = -70, 100
    pvmin, pvmax =-2.2e-08, 2.2e-08
    vortmin, vortmax = -9.2e-05, 9.2e-05

elif tc==7:
    basename='galewski'
    alpha = 0 # rotation angle
    Tf = 8
    hmin, hmax = 8400, 10500
    umin, umax = -20, 85
    vmin, vmax = -45, 45
    pvmin, pvmax = -1.5e-08, 2.4e-08
    vortmin, vortmax = -8.9e-05, 9.8e-05


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



#-------------------------------------------------------------------
# Loop over grids
for N in Ns:
    # Directory where the netcdf files are
    filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
    +".g"+str(gtype)+"."+dg+".adv"+str(adv)+".hord"+str(hord)+'.dd'+str(dd)+".mf"+str(mf)+".tf"+str(Tf)+"/"

    # Get the number of plots
    atmos_file = filepath+"atmos_daily.tile1.nc"
    data = xr.open_dataset(atmos_file, decode_times=False, chunks={'time': 1}).squeeze()
    nplots = int(len(data.time.values))
    dtplot = Tf/nplots

    #-----------------------------------------------------------------------------------------
    # Arrays to store the data that will be plotted
    # FV3 data
    h_fv3    = np.zeros((N,N,6,nplots+1))
    u_fv3    = np.zeros((N,N,6,nplots+1))
    v_fv3    = np.zeros((N,N,6,nplots+1))
    pv_fv3   = np.zeros((N,N,6,nplots+1))
    vort_fv3 = np.zeros((N,N,6,nplots+1))

    # Exact or endgame data
    h_ref = np.zeros((N,N,6,nplots+1)) # h at agrid
    u_ref = np.zeros((N,N,6,nplots+1)) # u at agrid
    v_ref = np.zeros((N,N,6,nplots+1)) # v at agrid

    h_error = np.zeros((N,N,6,nplots+1))
    u_error = np.zeros((N,N,6,nplots+1))
    v_error = np.zeros((N,N,6,nplots+1))

    #-----------------------------------------------------------------------------------------
    # This loop over tiles computes the maximum errors
    time = 0  
    for t in range(0,nplots+1):
        for tile in range(0,6):
            # Files to be opened
            atmos_file = filepath+"atmos_daily.tile"+str(tile+1)+".nc"
            grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

            # Load the files
            data = xr.open_dataset(atmos_file, decode_times=False, chunks={'time': 5}).squeeze()
            grid = xr.open_dataset(grid_file , decode_times=False, chunks={'time': 5}).squeeze()

            # Variable to be plotted (ps = fluid depth in the SW model)
            if t>=1:
               h_fv3[:,:,tile,t] = data['ps'][t-1,:,:].values
               u_fv3[:,:,tile,t] = data['ucomp'][t-1,:,:].values
               v_fv3[:,:,tile,t] = data['vcomp'][t-1,:,:].values
            else: #get ic
               h_fv3[:,:,tile,t] = data['ps_ic'][:,:].values
               u_fv3[:,:,tile,t] = data['ua_ic'][:,:].values
               v_fv3[:,:,tile,t] = data['va_ic'][:,:].values
            if tc==5:
               b = topography(grid, tc)
               h_fv3[:,:,tile,t] =  h_fv3[:,:,tile,t] + b

            # get vorticity and pv
            if t>=1:
               vort_fv3[:,:,tile,t] = data['vort'][t-1,:,:].values
               pv_fv3[:,:,tile,t] = data['pv'][t-1,:,:].values

            # Get reference solution
            if  tc == 2:
               h_ref[:,:,tile,t] = data['ps_ic'][:,:].values
               u_ref[:,:,tile,t] = data['ua_ic'][:,:].values
               v_ref[:,:,tile,t] = data['va_ic'][:,:].values
            elif tc>=2 and t>0:
               # endgame data
               # depth
               filename = "tc"+str(tc)+"_g"+str(gtype)+".s_N"+str(N)+"_h_t"+str(int(time*sec2day))\
               +"_face"+str(tile+1)+'.dat'
               print(filename)
               z = open(endgamedir+filename, 'rb')
               z = np.fromfile(z, dtype=np.float64)
               h_ref[:,:,tile,t] = np.reshape(z, (N,N))

               filename = "tc"+str(tc)+"_g"+str(gtype)+".s_N"+str(N)+"_u_t"+str(int(time*sec2day))\
               +"_face"+str(tile+1)+'.dat'
               z = open(endgamedir+filename, 'rb')
               z = np.fromfile(z, dtype=np.float64)
               u_ref[:,:,tile,t] = np.reshape(z, (N,N))

               #v
               filename = "tc"+str(tc)+"_g"+str(gtype)+".s_N"+str(N)+"_v_t"+str(int(time*sec2day))\
               +"_face"+str(tile+1)+'.dat'
               z = open(endgamedir+filename, 'rb')
               z = np.fromfile(z, dtype=np.float64)
               v_ref[:,:,tile,t] = np.reshape(z, (N,N))

        # time update
        time = time + dtplot 
    #-----------------------------------------------------------------------------------------
    h_error = h_ref - h_fv3
    u_error = u_ref - u_fv3
    v_error = v_ref - v_fv3

    h_error_min, h_error_max = np.amin(h_error[:,:,:,nplots]), np.amax(h_error[:,:,:,nplots])
    u_error_min, u_error_max = np.amin(u_error[:,:,:,nplots]), np.amax(u_error[:,:,:,nplots])
    v_error_min, v_error_max = np.amin(v_error[:,:,:,nplots]), np.amax(v_error[:,:,:,nplots])
 
    h_error_abs = max(abs(h_error_min), abs(h_error_max))
    u_error_abs = max(abs(u_error_min), abs(u_error_max))
    v_error_abs = max(abs(v_error_min), abs(v_error_max))

    #print(h_error_abs/np.amax(h_ref),u_error_abs/np.amax(u_ref),v_error_abs/np.amax(v_ref))
    #exit()
    #hmin_fv3, hmax_fv3 = hmin, hmax
    #hmin_fv3, hmax_fv3 = umin, umax
    #hmin_fv3, hmax_fv3 = vmin, vmax
    #-----------------------------------------------------------------------------------------


    fields = [h_fv3,u_fv3,v_fv3,vort_fv3,pv_fv3]
    field_names = ('h','u','v','vort','pv')
    fmins = [hmin, umin, vmin, vortmin, pvmin]
    fmaxs = [hmax, umax, vmax, vortmax, pvmax]

    efields = [h_error, u_error, v_error]
    ferror_maxs = [h_error_abs, u_error_abs, v_error_abs]
 
    #-----------------------------------------------------------------------------------------
    for efield, name, ferror_max in zip(efields, field_names, ferror_maxs):
       # Plot the error
       ts, te = 0, nplots
       time = 0.0

       for t in range(ts,te+1):
          dmax = str("{:.2e}".format(np.amax(abs(h_error[:,:,:,t]))))
          Time = str("{:.2e}".format(time))

          if tc==2:
             title = name+" error - Test case "+str(tc)+", $\\alpha$ = "+str(alpha)\
             +", time = "+str(Time)+" days, max = "+dmax+"\n"\
             +"Grid = "+gname+", "+dg+", adv="+str(adv)+", hord="+str(hord)+", divdamp="+str(dd)+", mf"+str(mf)+", N="+str(N)\
             +"\n"
          else:
             title = name+" error - Test case "+str(tc)\
             +", time = "+str(Time)+" days, max = "+dmax+"\n"\
             +"Grid = "+gname+", "+dg+", adv="+str(adv)+", hord="+str(hord)+", divdamp="+str(dd)+", mf"+str(mf)+", N="+str(N)\
             +"\n"

          # Filename
          filename = graphdir+name+"_error_tc"+str(tc)+'_t'+str(t)+'_alpha'+str(alpha)+'_C'+str(N)\
          +"_g"+str(gtype)+"_"+dg+"_adv"+str(adv)+"_hord"+str(hord)+"_dd"+str(dd)+"_mf"+str(mf)+"_tf"+str(Tf)

          # plot
          ferror_max2 = 0.071
          plot_scalarfield(efield[:,:,:,t], map_projection, title, filename, filepath, \
             'seismic', -ferror_max2, ferror_max2, dpi, figformat)
          #plot_scalarfield(efield[:,:,:,t], map_projection, title, filename, filepath, \
          #   'seismic', -ferror_max, ferror_max, dpi, figformat)
          print('plotted the file '+filename)

          # time update
          time = time + dtplot

    #-----------------------------------------------------------------------------------------
    exit()

    #-----------------------------------------------------------------------------------------
    # Plot the scalar field
    #-----------------------------------------------------------------------------------------
    dtplot = Tf/nplots
    time = 0
    if tc==6:
       gap=10
    else:
       gap = 1

    for field, name, fmin, fmax in zip(fields, field_names, fmins, fmaxs):
      time=0
      for t in range(0,nplots+1,gap):
        # finish the plotting
        # get min/max
        dmin = str("{:.2e}".format(fmin))
        dmax = str("{:.2e}".format(fmax))
        Time = str("{:.2e}".format(time))

        # title
        title = name+" - Test case "+str(tc)+", $\\alpha$ = "+str(alpha)\
            +", time = "+Time+" days"+", min = "+dmin+", max = "+dmax+"\n"\
            +"Grid = "+gname+", "+dg+", adv="+str(adv)+", hord="+str(hord)+", divdamp="+str(dd)+", mf"+str(mf)+", N="+str(N)\
            +"\n"

        # file
        filename = graphdir+name+"_tc"+str(tc)+'_t'+str(t)+'_alpha'+str(alpha)+'_C'+str(N)\
          +"_g"+str(gtype)+"_"+dg+"_adv"+str(adv)+"_hord"+str(hord)+"_dd"+str(dd)+"_mf"+str(mf)+"_tf"+str(Tf)

        # plot
        plot_scalarfield(field[:,:,:,t], map_projection, title, filename, filepath, \
        'jet', fmin, fmax, dpi, figformat)
        print('plotted the file '+filename)

        # time update
        time = time + gap*dtplot

    #-----------------------------------------------------------------------------------------
