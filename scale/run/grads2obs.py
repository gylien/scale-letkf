import numpy as np
from netCDF4 import Dataset

from datetime import datetime
from datetime import timedelta

import sys
import os

import struct

OVERW = True


exp = "2000m_InSnd_LT_SN14_Mac_0523"
top = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT" 

time = datetime(2000, 1, 1, 0, 40, 0)

def get_info(top,exp,time):

   fn_ref = os.path.join(top,exp,time.strftime('%Y%m%d%H%M%S'),"fcst/mean/history.pe000000.nc")

   HALO = 2

   nc = Dataset(fn_ref, 'r', format='NETCDF4')

   CZ = nc.variables["CZ"][HALO:-HALO]
   CXG = nc.variables["CXG"][HALO:-HALO]
   CYG = nc.variables["CYG"][HALO:-HALO]

   x = nc.variables["x"][:]
   y = nc.variables["y"][:]

   NPX = len(CXG) // len(x)
   NPY = len(CYG) // len(y)

   times = nc.variables["time"][:]


   INFO = {"CXG":CXG, "CYG":CYG, "CZ":CZ, "NPX":NPX, "NPY":NPY, "TIMES":times, 
           "DX":CXG[1]-CXG[0], "DY": CYG[1]-CYG[0]}

   return(INFO)



# file specification
XMAX = 144
YMAX = 144
ZMAX = 40
TMAX = 12
OTMAX = 12 

DX = 2000.0
DY = 2000.0
x0 = 1000.0
y0 = 1000.0

CZ_l = np.array([250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 
       5250, 5750, 6250, 6750, 7250, 7750, 8250, 8750, 9250, 9750, 10250, 10750, 
    11250, 11750, 12250, 12750, 13250, 13750, 14250, 14750, 15250, 15750, 
    16250, 16750, 17250, 17750, 18250, 18750, 19250, 19750])


DT = 300 # [s] history interval







def main(itime):

  SIGMA_VR = 5.0 # dBZ
  SIGMA_Z = 3.0  # m/s

  typ = "fcst"
  mem = "mean"

  NVAR = 2 # Z & VR
  RADAR_ZMAX = 10000.0

  INFO = get_info(top,exp,itime)
  INFO["NVAR"] =  NVAR

  ZMAX_SCALE = len(INFO["CZ"][INFO["CZ"] < RADAR_ZMAX])
  INFO["ZMAX_SCALE"] =  ZMAX_SCALE

  MIN_RADAR_REF = 15.0
  INFO["MIN_RADAR_REF"] =  MIN_RADAR_REF

  # Make sure the followings are consistent with obssim
  RADAR_X = 100.0e3 # radar location X [m]
  RADAR_Y = 100.0e3 # radar location Y [m]
  RADAR_Z = 0.0 # radar location Z [m]

  INFO["RADAR_X"] = RADAR_X
  INFO["RADAR_Y"] = RADAR_Y
  INFO["RADAR_Z"] = RADAR_Z

  INFO["SIGMA_VR"] = SIGMA_VR
  INFO["SIGMA_Z"] = SIGMA_Z

  otop = os.path.join(top,exp,"obs")
  os.makedirs(otop, exist_ok=True)

  fhd = "radar_"

  citime = itime.strftime('%Y%m%d%H%M%S')

  ofn = os.path.join(top,exp,citime,typ,mem,fhd + citime + "_" + mem + ".dat")

  
  # t level in grads
  for it in range(0,len(INFO["TIMES"])):

    z3d, vr3d =  read_grads_ZVR(ofn,it,INFO)

 
    idxs = np.where(z3d > INFO["MIN_RADAR_REF"])
 
    ctime = (itime + timedelta(seconds=INFO["TIMES"][it])).strftime('%Y%m%d%H%M%S')

    rfn = os.path.join(top,exp,citime,typ,mem,"rand_" + ctime + "_" + mem + ".npz")
    rand3d = prep_rand(rfn=rfn,mu=0.0,sigma=1.0,INFO=INFO)


    if OVERW:
       z3d += rand3d[0,:,:,:] * INFO["SIGMA_Z"]
       vr3d += rand3d[1,:,:,:] * INFO["SIGMA_VR"]

    obs_f = os.path.join(otop,"radar_" + ctime + '.dat')

    write_VR_obs(obs_f,z3d,vr3d,idxs,INFO)


def prep_rand(rfn,mu,sigma,INFO):

   XMAX = len(INFO["CXG"])
   YMAX = len(INFO["CYG"])

   if OVERW:   
     tmp1d = np.random.normal(mu, sigma, INFO["NVAR"] * INFO["ZMAX_SCALE"] * YMAX * XMAX)
     rand3d = np.reshape(tmp1d, (INFO["NVAR"],INFO["ZMAX_SCALE"],YMAX,XMAX))
     np.savez(rfn,rand3d)
   else:
     rand3d = np.load(rfn)['arr_0']

   return(rand3d)


def read_grads_ZVR(ofn,it,INFO):

   XMAX = len(INFO["CXG"])
   YMAX = len(INFO["CYG"])
   ZMAX = len(INFO["CZ"])

   count = XMAX * YMAX * ZMAX

   try:
      infile = open(ofn)
   except:
      print("Failed to open")
      print(ofn)
      sys.exit()

   rec = count * INFO["NVAR"] * it
   infile.seek(rec*4)
   tmp3d = np.fromfile(infile, dtype=np.dtype('<f4'), count=count)  #little endian
   z3d = np.reshape(tmp3d, (ZMAX,YMAX,XMAX))

   rec += count
   infile.seek(rec*4)
   tmp3d = np.fromfile(infile, dtype=np.dtype('<f4'), count=count)  #little endian
   vr3d = np.reshape(tmp3d, (ZMAX,YMAX,XMAX))

   return(z3d[0:INFO["ZMAX_SCALE"],:,:],vr3d[0:INFO["ZMAX_SCALE"],:,:])

def write_VR_obs(fn,z3d,vr3d,idxs,INFO):
   # idxs is a tuple that contains QC grid information based on nature run (without noise) data

   ID_Z = 4001
   ID_VR = 4002
  
   TYP = 22

   ofile = open(fn,'wb')

   ofile.seek(0)

   RADAR_LOC = [INFO["RADAR_X"], INFO["RADAR_Y"], INFO["RADAR_Z"]]

   for location in RADAR_LOC:
     obs_seq_ndarray = np.array(tuple([location]),dtype='>f4')

     ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
     obs_seq_ndarray.tofile(ofile)
     ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))

#   obs_seq_ndarray = np.array(tuple([RADAR_X,RADAR_Y,RADAR_Z]),dtype='>f4')
#
#   ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
#   obs_seq_ndarray.tofile(ofile)
#   ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))

   kidx, jidx, iidx = idxs


   print("Num Obs: ", len(kidx), fn)

#   for (k,j,i), val in np.ndenumerate(z3d):
   for p in range(len(kidx)):

#      print(k,j,i,val)
      obs_x = INFO["CXG"][0] + INFO["DX"] * iidx[p]
      obs_y = INFO["CYG"][0] + INFO["DY"] * jidx[p]
      obs_z = INFO["CZ"][kidx[p]]


#   No reflectivity obs
#
      obs_seq_ndarray = np.array(tuple([ID_Z,obs_x,obs_y,obs_z,
                                 z3d[kidx[p],jidx[p],iidx[p]],INFO["SIGMA_Z"],TYP]),
                                 dtype='>f4')
      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
      obs_seq_ndarray.tofile(ofile)
      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))

      obs_seq_ndarray = np.array(tuple([ID_VR,obs_x,obs_y,obs_z,
                                 vr3d[kidx[p],jidx[p],iidx[p]],INFO["SIGMA_VR"],TYP]),
                                 dtype='>f4')

      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
      obs_seq_ndarray.tofile(ofile)
      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))

   ofile.close()

#########

time = datetime(2000, 1, 1, 0, 40, 0)
main(time)


