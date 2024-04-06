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
endgamedir='../SHiELD_SRC/test/endgame/dump/'
figformat='png'

#--------------------------------------------------------------------------------------------------------
# test case
tc = 2

# 1d advection scheme
#hords = (5, 8)
hords = (8, )

# 2d advection scheme
advs = (1,)
advs = (1,2)
#advs = (1,1,1,1)
#advs = (2,2,2,2)

#grid type 0-equiedge; 2-equiangular
#gtypes = (2,)
#gtypes = (0, )
#gtypes = (2,2,2,2)
#gtypes = (0,0,0,0)
gtypes = (0,2)
#gtypes = (0,2,0,2)
#gtypes = (0,2)

#duogrid
#dgs = (2, 2, 2, 2)
dgs = (1, 1, 1, 1)
#dgs = (0, 0, 0, 0)

#mass fixers
mfs = (1, 1, 1, 1)

#divergence damping (only for sw)
if tc == 2:
   dds = (0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12)
else:
   #dds = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   dds = (0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12)
#dds = (0, 0.12, 0, 0.12)

# values of N
#Ns = (48, 96, 192, 384, 768)
#Ns = (48, 96, 192, 384)
Ns = (48, 96, 192, )
#Ns = (48, 96, )
ngrids = len(Ns)
#--------------------------------------------------------------------------------------------------------


# basename used in outputs
if tc==2 :
    basename='geobalance'
    alpha = 0
    Tf = 1
 
else:
   print('ERROR: invalid initial condition')

# initialize error arrays and counter
M = len(advs)
nplots = Tf+1
error_linf_h, error_l1_h, error_l2_h = np.zeros((len(Ns),M,len(hords),len(gtypes),nplots)), np.zeros((len(Ns),M,len(hords),len(gtypes),nplots)), np.zeros((len(Ns),M,len(hords),len(gtypes),nplots))
error_linf_u, error_l1_u, error_l2_u = np.zeros((len(Ns),M,len(hords),len(gtypes),nplots)), np.zeros((len(Ns),M,len(hords),len(gtypes),nplots)), np.zeros((len(Ns),M,len(hords),len(gtypes),nplots)) 
error_linf_v, error_l1_v, error_l2_v = np.zeros((len(Ns),M,len(hords),len(gtypes),nplots)), np.zeros((len(Ns),M,len(hords),len(gtypes),nplots)), np.zeros((len(Ns),M,len(hords),len(gtypes),nplots)) 

den_linf_h, den_l1_h, den_l2_h = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))),np.zeros((len(Ns),M,len(hords),len(gtypes)))
den_linf_u, den_l1_u, den_l2_u = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))),np.zeros((len(Ns),M,len(hords),len(gtypes)))
den_linf_v, den_l1_v, den_l2_v = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))),np.zeros((len(Ns),M,len(hords),len(gtypes)))

#--------------------------------------------------------------------------------------------------------
# Loop over all schemes - to get errors
for g in range(0, len(gtypes)):
  gtype = gtypes[g]
  for k in range(0, len(hords)):
    hord = hords[k]
    for m in range(0, M):
        dg = dgs[m]
        dd = str(dds[m])
        adv = advs[m]
        mf = str(mfs[m])

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

        # Now let us loop over all N values and compute the errors
        #------------------------------------------------------------------------------------------------
        for n in range(0, len(Ns)):
            N = Ns[n]
            # Directory where the netcdf files are
            filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
                 +".g"+str(gtype)+"."+dg+".adv"+str(adv)+".hord"+str(hord)+'.dd'+str(dd)+".mf"+mf+".tf"+str(Tf)+"/"

            print(filepath)
            #-----------------------------------------------------------------------------------------
            # Arrays to store the data that will be plotted
            # FV3 data
            h_fv3    = np.zeros((N,N,6,nplots))
            u_fv3    = np.zeros((N,N,6,nplots))
            v_fv3    = np.zeros((N,N,6,nplots))

            # Endgame data
            h_ref = np.zeros((N,N,6,nplots)) # h at agrid
            u_ref = np.zeros((N,N,6,nplots)) # u at agrid
            v_ref = np.zeros((N,N,6,nplots)) # v at agrid

            h_error = np.zeros((N,N,6,nplots))
            u_error = np.zeros((N,N,6,nplots))
            v_error = np.zeros((N,N,6,nplots))

            # areas (for l1 and l2 norms)
            areas = np.zeros((N,N,6)) 
            #-----------------------------------------------------------------------------------------

            time = 0
            dtplot = 1
            sec2day = 86400
            for t in range(0,nplots):
              for tile in range(0,6):
                # Files to be opened
                atmos_file = filepath+"atmos_daily.tile"+str(tile+1)+".nc"
                grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

                # Load the files
                data = xr.open_dataset(atmos_file, decode_times=False, chunks={'time': 5}).squeeze()
                grid = xr.open_dataset(grid_file , decode_times=False, chunks={'time': 5}).squeeze()

                areas[:,:,tile] = grid['area'][:,:].values

                # Variable to be plotted (ps = fluid depth in the SW model)
                if t>=1:
                   #h_fv3[:,:,tile,t] = data['ps'][t-1,:,:].values
                   #u_fv3[:,:,tile,t] = data['ucomp'][t-1,:,:].values
                   #v_fv3[:,:,tile,t] = data['vcomp'][t-1,:,:].values
                   h_fv3[:,:,tile,t] = data['ps'][:,:].values
                   u_fv3[:,:,tile,t] = data['ucomp'][:,:].values
                   v_fv3[:,:,tile,t] = data['vcomp'][:,:].values
                else: #get ic
                   h_fv3[:,:,tile,t] = data['ps_ic'][:,:].values
                   u_fv3[:,:,tile,t] = data['ua_ic'][:,:].values
                   v_fv3[:,:,tile,t] = data['va_ic'][:,:].values

                if t>0:
                # endgame data
                   # depth
                   filename = "tc"+str(tc)+"_g"+str(gtype)+".s_N"+str(N)+"_h_t"+str(int(time*sec2day))\
                   +"_face"+str(tile+1)+'.dat'
                   #print(filename)
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

              if t>0:
                 error_h = abs(h_fv3[:,:,:,t] - h_ref[:,:,:,t])
                 error_u = abs(u_fv3[:,:,:,t] - u_ref[:,:,:,t])
                 error_v = abs(v_fv3[:,:,:,t] - v_ref[:,:,:,t])



                 error_linf_h[n,m,k,g,t] = np.amax(error_h)/np.amax(h_ref[:,:,:,t])
                 error_linf_u[n,m,k,g,t] = np.amax(error_u)/np.amax(u_ref[:,:,:,t])
                 error_linf_v[n,m,k,g,t] = np.amax(error_v)#/np.amax(v_ref[:,:,:,t])

                 print(N,n,error_linf_h[n,m,k,g,t] )
                 #print( np.amax(abs(h_fv3[:,:,:,t] - h_ref[:,:,:,t]))/np.amax(abs(h_ref[:,:,:,t])))
                 #print( np.amax(abs(u_fv3[:,:,:,t] - u_ref[:,:,:,t]))/np.amax(abs(u_ref[:,:,:,t])))
                 #print( np.amax(abs(v_fv3[:,:,:,t] - v_ref[:,:,:,t])))#/np.amax(abs(v_ref[:,:,:,t])))
              
              # time update
              time = time + dtplot 
              
times = np.linspace(0,nplots,nplots+1)

#--------------------------------------------------------------------------------------------------------
# Loop over all schemes - to get errors
for g in range(0, len(gtypes)):
  gtype = gtypes[g]
  for k in range(0, len(hords)):
    hord = hords[k]
    for m in range(0, M):
        dg = dgs[m]
        dd = str(dds[m])
        adv = advs[m]
        mf = str(mfs[m])

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
        print(advname+dg+gname+dd)
        print(error_linf_h[2,m,k,g,:])
        print(error_linf_u[2,m,k,g,:])
        print(error_linf_v[2,m,k,g,:])
        print()
