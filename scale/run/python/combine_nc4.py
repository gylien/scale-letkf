from netCDF4 import Dataset
from datetime import datetime
from datetime import timedelta
import numpy as np
import os
import sys


def combine(INFO, time):

   print("Hello from combine")

   path = os.path.join(INFO["TOP"], INFO["EXP"], INFO["time"].strftime('%Y%m%d%H%M%S'), 
                       INFO["ftyp"], INFO["mem"])

   if INFO["ftyp"] == "fcst":
      fh = os.path.join(path, "history.pe") # file header
   elif INFO["ftyp"] == "anal" or INFO["ftyp"] == "gues":
      fh = os.path.join(path, "init.pe") # file header

  
   # get information from the head file & define root (combined) NetCDF
   fn = fh + str(0).zfill(6) + ".nc"
   nc_r, DIMS = new_netcdf(fn, INFO)
   
   # define variables
   VARS = var_def(nc_r,INFO)

   # main process loop
   for p in range(INFO["NPX"]*INFO["NPY"]):
      fn = fh + str(p).zfill(6) + ".nc"
      nc = Dataset(fn, 'r', format='NETCDF4')

      xrank = nc.rankidx[0]
      yrank = nc.rankidx[1]

      if p % 10 == 0:
         print("mem:",INFO["mem"],"process:",p)

      if xrank < INFO["NPBX"] or xrank > (INFO["NPX"] - 1 - INFO["NPBX"]) or \
         yrank < INFO["NPBY"] or yrank > (INFO["NPY"] - 1 - INFO["NPBY"]):
         continue

      imin = DIMS["xdim"] * (xrank - INFO["NPBX"])
      jmin = DIMS["ydim"] * (yrank - INFO["NPBY"])

      imax = imin + DIMS["xdim"]
      jmax = jmin + DIMS["ydim"]

      # var loop
      for nvar in var4d_names:
        var_local = nc.variables[nvar][:] # t,z,y,x 


        if INFO["ftyp"] == "fcst":
          VARS[nvar][:,:,jmin:jmax,imin:imax] = var_local[:,:INFO["ZMAX"],:,:]
 
        elif INFO["ftyp"] == "anal" or INFO["ftyp"] == "gues":
          VARS[nvar][jmin:jmax,imin:imax,:] = var_local[:,:,:INFO["ZMAX"]]

        # attribute
        if xrank == INFO["NPBX"] and yrank == INFO["NPBY"]:
          nc_r.variables[nvar].long_name = nc.variables[nvar].long_name
          nc_r.variables[nvar].units = nc.variables[nvar].units

      nc.close()

   nc_r.close()

#####################
def var_def(nc,INFO):
 
    VARS = {}
    for n, nvar in enumerate(var4d_names):
       print("define:", nvar)

       if INFO["ftyp"] == "fcst":
          VARS.setdefault(nvar, nc.createVariable(nvar, np.dtype('double').char, 
                         ('time','z','y','x')))
       else:
          if nvar == "MOMX":
            VARS.setdefault(nvar, nc.createVariable(nvar, np.dtype('double').char, ('y','xh','z')))
          elif nvar == "MOMY":
            VARS.setdefault(nvar, nc.createVariable(nvar, np.dtype('double').char, ('yh','x','z')))
          elif nvar == "MOMZ":
            VARS.setdefault(nvar, nc.createVariable(nvar, np.dtype('double').char, ('y','x','zh')))
          else:
            VARS.setdefault(nvar, nc.createVariable(nvar, np.dtype('double').char, ('y','x','z')))

    return(VARS)

def new_netcdf(fn, INFO): 

   nc = Dataset(fn, 'r', format='NETCDF4')

   cxg = nc.variables["CXG"][:] 
   cyg = nc.variables["CYG"][:] 
   cz = nc.variables["CZ"][:] 
   zlevs = nc.variables["z"][:] 
   NX = len(nc.variables["x"][:])
   NY = len(nc.variables["y"][:])
   NZ = len(nc.variables["z"][:])
   if INFO["ftyp"] == "fcst":
     times = nc.variables["time"][:] 

   nc.close()

   if INFO["ftyp"] == "fcst":
      fhead = "history_"
   elif INFO["ftyp"] == "anal" or INFO["ftyp"] == "gues":
      fhead = "init_"

   # Determine zmax
   if len(zlevs) < INFO["ZMAX"]:
      INFO["ZMAX"] = len(zlevs)

   # Define new NetCDF file
   fn_r =  os.path.join(INFO["TOP"], INFO["EXP"], INFO["time"].strftime('%Y%m%d%H%M%S'),
                        INFO["ftyp"], fhead)

   print("New file:",fn_r + INFO["mem"] + '.nc')


   nc_r = Dataset(fn_r + INFO["mem"] + '.nc', 'w', format='NETCDF4')
   nc_r.createDimension('x', len(cxg) - 2 * INFO["HALO"] - 2 * NX * INFO["NPBX"])  
   nc_r.createDimension('y', len(cyg) - 2 * INFO["HALO"] - 2 * NY * INFO["NPBY"])
   nc_r.createDimension('z', len(zlevs[:INFO["ZMAX"]]))  

   nc_r.createDimension('xh', len(cxg) - 2 * INFO["HALO"] - 2 * NX * INFO["NPBX"])  
   nc_r.createDimension('yh', len(cyg) - 2 * INFO["HALO"] - 2 * NY * INFO["NPBY"])
   nc_r.createDimension('zh', len(zlevs[:INFO["ZMAX"]]))  

   if INFO["ftyp"] == "fcst":
     nc_r.createDimension('time', len(times))
   else:
     nc_r.createDimension('time', 1)

   tdim = nc_r.createVariable('time', np.dtype('double').char, ('time',))
   tdim.long_name = 'time of test variable'
   tdim.units = 'seconds since ' + INFO["time"].strftime("%Y-%m-%d %H:%M:%S")
   if INFO["ftyp"] == "fcst":
     tdim[:] = times[:] 
   else:
     tdim[:] = 0.0

   imin_g = INFO["HALO"] + INFO["NPBX"] * NX
   imax_g = len(cxg) - INFO["HALO"] - INFO["NPBX"] * NX
   xdim = nc_r.createVariable('x', np.dtype('double').char, ('x',))
   xdim.long_name = 'X dimension'
   xdim.units = 'm'
   xdim[:] = cxg[imin_g:imax_g]

   xdimh = nc_r.createVariable('xh', np.dtype('double').char, ('xh',))
   xdimh.long_name = 'X half dimension'
   xdimh.units = 'm'
   xdimh[:] = cxg[imin_g:imax_g] + (cxg[1] - cxg[0]) * 0.5


   jmin_g = INFO["HALO"] + INFO["NPBY"] * NY
   jmax_g = len(cyg) - INFO["HALO"] - INFO["NPBY"] * NY
   ydim = nc_r.createVariable('y', np.dtype('double').char, ('y',))
   ydim.long_name = 'Y dimension'
   ydim.units = 'm'
   ydim[:] = cyg[jmin_g:jmax_g]

   ydimh = nc_r.createVariable('yh', np.dtype('double').char, ('yh',))
   ydimh.long_name = 'Y half dimension'
   ydimh.units = 'm'
   ydimh[:] = cyg[jmin_g:jmax_g] + (cyg[1] - cyg[0]) * 0.5


   zdim = nc_r.createVariable('z', np.dtype('double').char, ('z',))
   zdim.long_name = 'Z dimension'
   zdim.units = 'm'
   zdim[:] = zlevs[:INFO["ZMAX"]]

   zdimh = nc_r.createVariable('zh', np.dtype('double').char, ('zh',))
   zdimh.long_name = 'Z half dimension'
   zdimh.units = 'm'
   zdimh[:] = (zlevs[:INFO["ZMAX"]] + zlevs[1:INFO["ZMAX"]+1]) * 0.5

   DIMS = {"tdim":len(tdim), "zdim":len(zdim), 
           "gydim":len(ydim), "gxdim":len(xdim),
           "xdim":NX, "ydim":NY}

   return(nc_r, DIMS)

########################

TOP = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT"
EXP = "2000m_InSnd_LT_SN14_Mac_0605_PAWR"
time = datetime(2000,1,1,0,0,0)

ftyp = "fcst"
#ftyp = "anal"
#ftyp = "gues"

# Number of processes in X direction
NPX = 8
# Number of processes in Y direction
NPY = 8

# Number of processes (X) that will NOT be combined # Total: NPX - "2" * NPBX
NPBX = 0
# Number of processes (Y) that will NOT be combined # Total: NPY - "2" * NPBY
NPBY = 0

ZMAX = 39 # Largest index number in Z direction

mem = "mean"
HALO = 2

INFO = {"TOP":TOP, "EXP":EXP, "time": time, "ftyp":ftyp, "NPX":NPX, "NPY":NPY, "mem":mem, "HALO":HALO, "ZMAX":ZMAX, "NPBX":NPBX, "NPBY":NPBY}

var4d_names = ["U", "V", "W", "T", "PRES",\
               "QV", "QC", "QR", "QI", "QS", "QG", \
               "QCRG_C", "QCRG_R", "QCRG_I", "QCRG_S", "QCRG_G",\
               "Ex", "Ey", "Ez", "Epot", "PosFLASH", "NegFLASH", "LTpath"]

#var4d_names = ["MOMX", "MOMY", "MOMZ", "RHOT", "DENS",\
#               "QV", "QC", "QR", "QI", "QS", "QG", \
#               "CDNS_QC", "CDNS_QR", "CDNS_QI", "CDNS_QS", "CDNS_QG"]

#var4d_names = ["QG"]


etime = datetime(2000,1,1,1,10,0)
stime = etime # DEBUG

dt = timedelta(seconds=30)

# MPI
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


print("myrank",size,rank)


MEMBER = 80
MEMBER = 2
mem_list = [str(x).zfill(4) for x in range(1,MEMBER+1)] 
mem_list.append("mean")
#mem_list = ["mean", "0001"]

if size > len(mem_list):
   print("Make sure MPI size & mem_list")
else:
   dmem = len(mem_list) // size
   sidx = dmem * rank
   eidx = dmem * (rank + 1)

   if (rank + 1) == size:
      eidx = len(mem_list)

   print("Myrank:", rank, "Mem:",mem_list[sidx:eidx])


time = stime
while time <= etime:

   for mem in mem_list[sidx:eidx]:
      INFO["mem"] = mem
      INFO["time"] = time
      combine(INFO, time)
   
   time += dt
   
print("myrank",rank,"finished")
sys.exit()
