import netCDF4 
from datetime import datetime
import numpy as np
import os
import sys


#TOP = "/Users/takumihonda/Documents/python/OFP_TEST"
TOP = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT"
EXP = "1000m_InSnd_LT"
time = datetime(2000,1,1,0,0,0)
ftyp = "fcst"
# Number of processes in X direction
NPX = 16
# Number of processes in Y direction
NPY = 12
# Number of processes (X) that will NOT be combined # Total: NPX - "2" * NPBX
NPBX = 3 
# Number of processes (Y) that will NOT be combined # Total: NPY - "2" * NPBY
NPBY = 3 

ZMAX = 30 # Largest index number in Z direction

mem = "mean"
HALO = 2

INFO = {"TOP":TOP, "EXP":EXP, "time": time, "ftyp":ftyp, "NPX":NPX, "NPY":NPY, "mem":mem, "HALO":HALO, "ZMAX":ZMAX, "NPBX":NPBX, "NPBY":NPBY}

var4d_names = ["U", "V", "W", "T", "PRES",\
#               "QV","QHYD",
               "QV", "QC", "QR", "QI", "QS", "QG", \
               "QCRG_C", "QCRG_R", "QCRG_I", "QCRG_S", "QCRG_G",\
#               "QCRG_TOT", \
               "Ex", "Ey", "Ez", "Epot", "PosFLASH", "NegFLASH", "LTpath"]


def combine(INFO):

   path = os.path.join(INFO["TOP"], INFO["EXP"], INFO["time"].strftime('%Y%m%d%H%M%S'), 
                       INFO["ftyp"], INFO["mem"])

   if INFO["ftyp"] == "fcst":
      fh = os.path.join(path, "history.pe") # file header

  
   # get information from the head file & define root (combined) NetCDF
   fn = fh + str(0).zfill(6) + ".nc"
   nc_r, DIMS = new_netcdf(fn, INFO)
   

   # define variables
   VARS = var_def(nc_r)


   # main process loop
   for p in range(INFO["NPX"]*INFO["NPY"]):
      fn = fh + str(p).zfill(6) + ".nc"
      nc = netCDF4.Dataset(fn, 'r', format='NETCDF4')

      xrank = nc.rankidx[0]
      yrank = nc.rankidx[1]

      if p % 10 == 0:
         print("mem:",INFO["mem"],"process:",p)

      if xrank < INFO["NPBX"] or xrank > (INFO["NPX"] - 1 - INFO["NPBX"]) or \
         yrank < INFO["NPBY"] or yrank > (INFO["NPY"] - 1 - INFO["NPBY"]):
         continue

      imin = DIMS["xdim"] * (xrank - INFO["NPBX"])
      jmin = DIMS["ydim"] * (yrank - INFO["NPBY"])

      # var loop
      for nvar in var4d_names:
        var_local = nc.variables[nvar][:,:,:,:] # t,z,y,x 

        VARS[nvar][:,:,jmin:jmin+DIMS["ydim"],imin:imin+DIMS["xdim"]] = \
        var_local[:,:INFO["ZMAX"],:,:]
 
        # attribute
        if xrank == INFO["NPBX"] and yrank == INFO["NPBY"]:
          nc_r.variables[nvar].long_name = nc.variables[nvar].long_name
          nc_r.variables[nvar].units = nc.variables[nvar].units

      nc.close()

   nc_r.close()

#####################
def var_def(nc):
 
    VARS = {}
    for n, nvar in enumerate(var4d_names):
       print("define:", nvar)

       VARS.setdefault(nvar, nc.createVariable(nvar, np.dtype('double').char, 
                       ('time','z','y','x')))

    return(VARS)

def new_netcdf(fn, INFO): 

   nc = netCDF4.Dataset(fn, 'r', format='NETCDF4')

   cxg = nc.variables["CXG"][:] 
   cyg = nc.variables["CYG"][:] 
   cz = nc.variables["CZ"][:] 
   zlevs = nc.variables["z"][:] 
   NX = len(nc.variables["x"][:])
   NY = len(nc.variables["y"][:])
   NZ = len(nc.variables["z"][:])
   times = nc.variables["time"][:] 

   nc.close()

   # Determine zmax
   if len(zlevs) < INFO["ZMAX"]:
      INFO["ZMAX"] = len(zlevs)

   # Define new NetCDF file
   if INFO["ftyp"] == "fcst":
      fn_r =  os.path.join(INFO["TOP"], INFO["EXP"], INFO["time"].strftime('%Y%m%d%H%M%S'),
                           INFO["ftyp"], "history_")

   print("New file:",fn_r + INFO["mem"] + '.nc')

   nc_r = netCDF4.Dataset(fn_r + INFO["mem"] + '.nc', 'w', format='NETCDF4')
   nc_r.createDimension('x', len(cxg) - 2 * INFO["HALO"] - 2 * NX * INFO["NPBX"])  
   nc_r.createDimension('y', len(cyg) - 2 * INFO["HALO"] - 2 * NY * INFO["NPBY"])
   nc_r.createDimension('z', len(zlevs[:INFO["ZMAX"]]))  
   nc_r.createDimension('time', len(times))

   tdim = nc_r.createVariable('time', np.dtype('double').char, ('time',))
   tdim.long_name = 'time of test variable'
   tdim.units = 'seconds since ' + INFO["time"].strftime("%Y-%m-%d %H:%M:%S")
   tdim[:] = times[:] 

   imin_g = INFO["HALO"] + INFO["NPBX"] * NX
   imax_g = len(cxg) - INFO["HALO"] - INFO["NPBX"] * NX
   xdim = nc_r.createVariable('x', np.dtype('double').char, ('x',))
   xdim.long_name = 'X dimension'
   xdim.units = 'm'
   xdim[:] = cxg[imin_g:imax_g]

   jmin_g = INFO["HALO"] + INFO["NPBY"] * NY
   jmax_g = len(cyg) - INFO["HALO"] - INFO["NPBY"] * NY
   ydim = nc_r.createVariable('y', np.dtype('double').char, ('y',))
   ydim.long_name = 'Y dimension'
   ydim.units = 'm'
   ydim[:] = cyg[jmin_g:jmax_g]

   zdim = nc_r.createVariable('z', np.dtype('double').char, ('z',))
   zdim.long_name = 'Z dimension'
   zdim.units = 'm'
   zdim[:] = zlevs[:INFO["ZMAX"]]

   DIMS = {"tdim":len(tdim), "zdim":len(zdim), 
           "gydim":len(ydim), "gxdim":len(xdim),
           "xdim":NX, "ydim":NY}

   return(nc_r, DIMS)


MEMBER = 80
MEMBER = 1
mem_list = [str(x).zfill(4) for x in range(1,MEMBER+1)] 
#mem_list.append("mean")
#mem_list = ["0002"]
print(mem_list)
for mem in mem_list:
  INFO["mem"] = mem
  combine(INFO)


