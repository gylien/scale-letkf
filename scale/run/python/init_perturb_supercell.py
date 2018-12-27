import numpy as np
from netCDF4 import Dataset

from datetime import datetime

import sys
import os

OVERW = True

#exp = "8km_sc"
exp = "2km_CZ2003"

#top = "/data6/honda/SCALE-LETKF/scale_lt_devel_20181002/OUTPUT"
#top = "/home/honda/work/OUTPUT"

top = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT"

### pertub only below Z < ZMAX
##ZMAX = 30

# ensemble size
MEMBER = 20

#SCALE_NP = 1
NPX = 16
NPY = 8
SCALE_NP = NPX * NPY

# variable list
# MOMX/Y/Z needs a consideration for a staggered grid system (2018/10/13)
VAR_LIST = ["RHOT","MOMX","MOMY","MOMZ"]

# Horizontal buffer size (grids) where no perturbations are added
HBUF = 2
# Vertical buffer size (grids) where no perturbations are added
VTBUF = 1 # top
VBBUF = 1 # bottom

def nc_name(top,exp,time,typ,m,p):
   if m == 0:
      mem = "mean"
   else:
      mem = str(m).zfill(4)

   if typ == "fcst":
      fh = "history.pe"
   else:
      fh = "init.pe"
   
   return(os.path.join(top,exp,time.strftime('%Y%m%d%H%M%S'),typ,mem,fh + str(p).zfill(6) + ".nc"))
 

def main(time):

  typ = "anal"
  m = 1
  p = 0


  # process loop
  for p in range(SCALE_NP):
  #for p in range(1): DEBUG
     print(p)
     fn_nat = nc_name(top,exp,time,typ,0,p) # nature run (m=0)

     try:
       nc_nat = Dataset(fn_nat, "r", format="NETCDF4")
     except:
       print("Failed to read: ",fn_nat)

     rank_idxs = nc_nat.rankidx
     rank_i = rank_idxs[0]
     rank_j = rank_idxs[1]

     dens3d = nc_nat.variables["DENS"][:,:,:]


     if rank_i == 0:
       imin = HBUF 
     else:
       imin = 1
 
     if rank_i == (NPX - 1):
       imax = dens3d.shape[0] - HBUF
     else:
       imax = dens3d.shape[0]

     if rank_j == 0:
       jmin = HBUF 
     else:
       jmin = 1
 
     if rank_j == (NPY - 1):
       jmax = dens3d.shape[1] - HBUF
     else:
       jmax = dens3d.shape[1]


     kmin = VBBUF 
     kmax = dens3d.shape[2] - VTBUF
    
     ISIZE = imax - imin 
     JSIZE = jmax - jmin 
     KSIZE = kmax - kmin

#     #rsize = (MEMBER,dens3d.shape[0],dens3d.shape[1],dens3d.shape[2]) # random numer array size
#     rsize = (MEMBER,dens3d.shape[0]-2*HBUF,dens3d.shape[1]-2*HBUF,dens3d.shape[2]-VTBUF-VBBUF) # random numer array size
     rsize = (MEMBER,ISIZE,JSIZE,KSIZE) # random numer array size

     var3d_nat = nc_nat.variables[VAR_LIST[0]][imin:imax,jmin:jmax,kmin:kmax] # reference data from nature run
     #nc_nat.close()

     #
     fn_mem = nc_name(top,exp,time,typ,1,p)
     nc_mem = Dataset(fn_mem, "r", format="NETCDF4")
     #
     var3d = nc_mem.variables[VAR_LIST[0]][imin:imax,jmin:jmax,kmin:kmax]
     nc_mem.close()

     if not OVERW:
       if np.abs(np.min(var3d_nat) - np.min(var3d)) > 0.1 or np.abs(np.max(var3d_nat) - np.max(var3d)) > 0.1:
          print("Ensemble data seem to be perturbed!")
          print("CHECK")
          sys.exit()


     # var loop
     for vname in VAR_LIST:
        print(vname)
 
     
        #if varname == "RHOT":
        sigma = 3.0 # (K) or (m/s)
        #sigma = 0.5 # (K) or (m/s)
        #sigma = 0.0 # (K) or (m/s)
        rand3d = np.random.normal(loc=0.0,scale=sigma,size=rsize)


        # ensemble member loop
        for m in range(1,MEMBER+1):
        #for m in range(1,2): # DEBUG
           fn_mem = nc_name(top,exp,time,typ,m,p)
           print(fn_mem) 
           nc_mem = Dataset(fn_mem, "r+", format="NETCDF4")

           var3d = nc_nat.variables[vname][imin:imax,jmin:jmax,kmin:kmax]


           if vname == "RHOT":
             nc_mem.variables[vname][imin:imax,jmin:jmax,kmin:kmax] =(var3d[:,:,:] / dens3d[imin:imax,jmin:jmax,kmin:kmax] + rand3d[m-1,:,:,:]) * dens3d[imin:imax,jmin:jmax,kmin:kmax]
           if vname == "MOMX":
             dens3d_tmp = (dens3d[imin:imax,jmin:jmax,kmin:kmax] + dens3d[imin-1:imax-1,jmin:jmax,kmin:kmax] ) * 0.5
             nc_mem.variables[vname][imin:imax,jmin:jmax,kmin:kmax] =(var3d[:,:,:] / dens3d_tmp + rand3d[m-1,:,:,:]) * dens3d_tmp
           if vname == "MOMY":
             dens3d_tmp = (dens3d[imin:imax,jmin:jmax,kmin:kmax] + dens3d[imin:imax,jmin-1:jmax-1,kmin:kmax] ) * 0.5
             nc_mem.variables[vname][imin:imax,jmin:jmax,kmin:kmax] =(var3d[:,:,:] / dens3d_tmp + rand3d[m-1,:,:,:]) * dens3d_tmp
           if vname == "MOMZ":
             dens3d_tmp = (dens3d[imin:imax,jmin:jmax,kmin:kmax] + dens3d[imin:imax,jmin:jmax,kmin-1:kmax-1] ) * 0.5
             nc_mem.variables[vname][imin:imax,jmin:jmax,kmin:kmax] =(var3d[:,:,:] / dens3d_tmp + rand3d[m-1,:,:,:]) * dens3d_tmp

           nc_mem.close()


#########

time = datetime(2000,1,1,0,35,0)
main(time)

