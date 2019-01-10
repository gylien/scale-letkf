import numpy as np
from netCDF4 import Dataset

from datetime import datetime
from datetime import timedelta

import sys
import os

import struct

#top = "/data6/honda/SCALE-LETKF/scale_lt_devel_20181002/OUTPUT"
#exp = "8km_sc"

exp = "2km_CZ2003"
#top = "/home/honda/work/OUTPUT"
top = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT" 

MIN_RADAR_REF = 15.0
RADAR_ZMAX = 10000.0

OVERW = True
#OVERW = False

# file specification
XMAX = 144
YMAX = 144
ZMAX = 40
TMAX = 12
OTMAX = 12 
NVAR = 2 # Z & VR

DX = 2000.0
DY = 2000.0
x0 = 1000.0
y0 = 1000.0

CZ_l = np.array([250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 
       5250, 5750, 6250, 6750, 7250, 7750, 8250, 8750, 9250, 9750, 10250, 10750, 
    11250, 11750, 12250, 12750, 13250, 13750, 14250, 14750, 15250, 15750, 
    16250, 16750, 17250, 17750, 18250, 18750, 19250, 19750])

ZMAX_SCALE = len(CZ_l[CZ_l < RADAR_ZMAX])


DT = 300 # [s] history interval

SIGMA_VR = 5.0 # dBZ
SIGMA_Z = 3.0  # m/s



# Make sure the followings are consistent with obssim
RADAR_X = 100.0e3 # radar location X [m]
RADAR_Y = 100.0e3 # radar location Y [m]
RADAR_Z = 0.0 # radar location Z [m]



def main(itime):

  typ = "fcst"
  mem = "mean"

  otop = os.path.join(top,exp,"obs")
  os.makedirs(otop, exist_ok=True)

  fhd = "radar_"

  citime = itime.strftime('%Y%m%d%H%M%S')

  ofn = os.path.join(top,exp,citime,typ,mem,fhd + citime + "_" + mem + ".dat")

  
  # t level in grads
  for it in range(1,OTMAX+1):

    z3d, vr3d =  read_grads_ZVR(ofn,it)
 
    idxs = np.where(z3d > MIN_RADAR_REF)
 
    ctime = (itime + timedelta(seconds=DT*(it-1))).strftime('%Y%m%d%H%M%S')

    rfn = os.path.join(top,exp,citime,typ,mem,"rand_" + ctime + "_" + mem + ".npz")
    rand3d = prep_rand(rfn=rfn,mu=0.0,sigma=1.0)


    if OVERW:
       z3d += rand3d[0,:,:,:] * SIGMA_Z
       vr3d += rand3d[1,:,:,:] * SIGMA_VR

    obs_f = os.path.join(otop,"radar_" + ctime + '.dat')

    write_VR_obs(obs_f,z3d,vr3d,idxs)


def prep_rand(rfn,mu,sigma):


   if OVERW:   
     tmp1d = np.random.normal(mu, sigma, NVAR * ZMAX_SCALE * YMAX * XMAX)
     rand3d = np.reshape(tmp1d, (NVAR,ZMAX_SCALE,YMAX,XMAX))
     np.savez(rfn,rand3d)
   else:
     rand3d = np.load(rfn)['arr_0']

   return(rand3d)


def read_grads_ZVR(ofn,it):

  count = XMAX * YMAX * ZMAX

  try:
     infile = open(ofn)
  except:
     print("Failed to open")
     print(ofn)
     sys.exit()

  rec = count * NVAR * (it - 1)
  infile.seek(rec*4)
  tmp3d = np.fromfile(infile, dtype=np.dtype('<f4'), count=count)  #little endian
  z3d = np.reshape(tmp3d, (ZMAX,YMAX,XMAX))

  rec += count
  infile.seek(rec*4)
  tmp3d = np.fromfile(infile, dtype=np.dtype('<f4'), count=count)  #little endian
  vr3d = np.reshape(tmp3d, (ZMAX,YMAX,XMAX))


  return(z3d[0:ZMAX_SCALE,:,:],vr3d[0:ZMAX_SCALE,:,:])

def write_VR_obs(fn,z3d,vr3d,idxs):
   # idxs is a tuple that contains QC grid information based on nature run (without noise) data

   ID_Z = 4001
   ID_VR = 4002
  
   TYP = 22

   ofile = open(fn,'wb')

   ofile.seek(0)

   RADAR_LOC = [RADAR_X, RADAR_Y, RADAR_Z]

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
      obs_x = x0 + DX * iidx[p]
      obs_y = y0 + DY * jidx[p]
      obs_z = CZ_l[kidx[p]]

#   No reflectivity obs
#
      obs_seq_ndarray = np.array(tuple([ID_Z,obs_x,obs_y,obs_z,
                                 z3d[kidx[p],jidx[p],iidx[p]],SIGMA_Z,TYP]),
                                 dtype='>f4')
      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
      obs_seq_ndarray.tofile(ofile)
      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))

      obs_seq_ndarray = np.array(tuple([ID_VR,obs_x,obs_y,obs_z,
                                 vr3d[kidx[p],jidx[p],iidx[p]],SIGMA_VR,TYP]),
                                 dtype='>f4')

      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
      obs_seq_ndarray.tofile(ofile)
      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))

   ofile.close()

#########

time = datetime(2000, 1, 1, 0, 35, 0)
main(time)


