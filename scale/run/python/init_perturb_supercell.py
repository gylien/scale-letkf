import numpy as np
from netCDF4 import Dataset

from datetime import datetime

import sys
import os

# pertub only below Z < ZMAX
ZMAX = 30

# ensemble size
MEMBER = 4

SCALE_NP = 4

# variable list
# MOMX/Y/Z needs a consideration for a staggered grid system (2018/10/13)
VAR_LIST = ["RHOT"]

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
 

def main():

  top = "/data6/honda/SCALE-LETKF/scale_lt_devel_20180731/OUTPUT"
  exp = "8km_sc"
  time = datetime(2000,1,1,0,0)
  typ = "anal"
  m = 1
  p = 0


  # process loop
  for p in range(SCALE_NP):
  #for p in range(1): DEBUG
     print(p)
     fn_nat = nc_name(top,exp,time,typ,0,p) # nature run (m=0)
     nc_nat = Dataset(fn_nat, "r", format="NETCDF4")

     dens3d = nc_nat.variables["DENS"][:,:,:]
     rsize = (MEMBER,dens3d.shape[0],dens3d.shape[1],dens3d.shape[2]) # random numer array size

     var3d_nat = nc_nat.variables[VAR_LIST[0]][:,:,:] # reference data from nature run
     nc_nat.close()
     #
     fn_mem = nc_name(top,exp,time,typ,1,p)
     nc_mem = Dataset(fn_mem, "r", format="NETCDF4")
     #
     var3d = nc_mem.variables[VAR_LIST[0]][:,:,:]
     nc_mem.close()

     if np.abs(np.min(var3d_nat) - np.min(var3d)) > 0.1 or np.abs(np.max(var3d_nat) - np.max(var3d)) > 0.1:
        print("Ensemble data seem to be perturbed!")
        print("CHECK")
        sys.exit()


     # var loop
     for vname in VAR_LIST:
        print(vname)
 
     
        #if varname == "RHOT":
        #sigma = 3.0 # (K) or (m/s)
        sigma = 0.1 # (K) or (m/s)
        rand3d = np.random.normal(loc=0.0,scale=sigma,size=rsize)

        # ensemble member loop
        for m in range(1,MEMBER):
        #for m in range(1,2): # DEBUG
           fn_mem = nc_name(top,exp,time,typ,m,p)
           print(fn_mem) 
           nc_mem = Dataset(fn_mem, "r+", format="NETCDF4")

           var3d = nc_mem.variables[vname][:,:,:]


           if vname == "RHOT":
             nc_mem.variables[vname][:,:,:] = (var3d[:,:,:] / dens3d[:,:,:] + rand3d[m-1,:,:,:]) / dens3d[:,:,:]
           nc_mem.close()


#########

main()

