import numpy as np
from netCDF4 import Dataset

from datetime import datetime

import sys
import os

OVERW = True

#exp = "2km_CZ2003"
#exp = "2000m_InSnd_LT_SN14_Mac_0523"
exp = "2000m_InSnd_LT_SN14_Mac_0605"


top = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT"


# ensemble size
MEMBER = 80

#SCALE_NP = 1
NPX = 8
NPY = 8
SCALE_NP = NPX * NPY

# No perturbation processes
NOPT_PRC=2 

# variable list
# MOMX/Y/Z needs a consideration for a staggered grid system (2018/10/13)
VAR_LIST = ["RHOT","MOMX","MOMY","MOMZ"]
VAR_LIST = ["RHOT"]

# Vertical buffer size (grids) where no perturbations are added
VTBUF = 30 # top (below 5 km)

IHALO = 2
JHALO = 2
KHALO = 2

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
  p = 0


  # process loop
  for p in range(SCALE_NP):
     print(p)
     fn_nat = nc_name(top,exp,time,typ,0,p) # nature run (m=0)

     try:
       nc_nat = Dataset(fn_nat, "r", format="NETCDF4")
     except:
       print("Failed to read: ",fn_nat)

     rank_idxs = nc_nat.rankidx
     rank_i = rank_idxs[0]
     rank_j = rank_idxs[1]

     # No perturbation for the processes near the edge
     if rank_i < NOPT_PRC or rank_j < NOPT_PRC:
        continue

     if rank_i > (NPX - NOPT_PRC - 1) or rank_j > (NPY - NOPT_PRC - 1):
        continue


     dens3d = nc_nat.variables["DENS"][:,:,:]

     imin = 0 # IHALO
     jmin = 0 # JHALO
     kmin = 0 # KHALO

     imax = dens3d.shape[0] #- IHALO
     jmax = dens3d.shape[1] #- JHALO
     kmax = dens3d.shape[2] #- KHALO - VTBUF


     ISIZE = imax - imin 
     JSIZE = jmax - jmin 
     KSIZE = kmax - kmin


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
 
     
        if vname == "RHOT":
          sigma = 0.1 # (K)
        elif vname == "MOMX" or vname == "MOMY" or vname == "MOMZ":
          sigma = 3.0 # (m/s)

        rand3d = np.random.normal(loc=0.0,scale=sigma,size=rsize)
        rand3d -= np.mean(rand3d,axis=0) # ensure mean is zero
 

        # ensemble member loop
        for m in range(1,MEMBER+1):
           fn_mem = nc_name(top,exp,time,typ,m,p)
           print(fn_mem) 
           nc_mem = Dataset(fn_mem, "r+", format="NETCDF4")

           #var3d = nc_nat.variables[vname][imin:imax,jmin:jmax,kmin:kmax]
           var3d = nc_mem.variables[vname][:,:,:]


           if vname == "RHOT":
             nc_mem.variables[vname][:,:,:-VTBUF] =(var3d[:,:,:-VTBUF] / dens3d[:,:,:-VTBUF] + rand3d[m-1,:,:,:-VTBUF]) * dens3d[:,:,:-VTBUF]
           elif vname == "MOMX":
             # DENS is constant in horizontal
             dens3d_tmp = dens3d[:,:,:-VTBUF] 
             nc_mem.variables[vname][:,:,:-VTBUF] = (var3d[:,:,:-VTBUF] / dens3d_tmp + rand3d[m-1,:,:,:-VTBUF]) * dens3d_tmp
           elif vname == "MOMY":
             dens3d_tmp = dens3d[:,:,:-VTBUF] 
             nc_mem.variables[vname][:,:,:-VTBUF] = (var3d[:,:,:-VTBUF] / dens3d_tmp + rand3d[m-1,:,:,:-VTBUF]) * dens3d_tmp
           elif vname == "MOMZ":
             dens3d_tmp = (dens3d[:,:,:-VTBUF] + dens3d[:,:,1:-VTBUF+1] ) * 0.5 # zh-level dens3d
             nc_mem.variables[vname][:,:,:-VTBUF] = (var3d[:,:,:-VTBUF] / dens3d_tmp + rand3d[m-1,:,:,:-VTBUF]) * dens3d_tmp

           nc_mem.close()


#########

time = datetime(2000, 1, 1, 0, 0, 0)
main(time)

