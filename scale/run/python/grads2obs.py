import numpy as np
from netCDF4 import Dataset

from datetime import datetime
from datetime import timedelta

import sys
import os

import struct

#top = "/data6/honda/SCALE-LETKF/scale_lt_devel_20181002/OUTPUT"
#exp = "8km_sc"

exp = "2000m_InSnd_LT_SN14_Mac_0605"
top = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT" 

HALO = 2

VR_MIN_RADAR_REF = 15.0
CLR_RADAR_REF = 5.0
RADAR_ZMAX = 10000.0

THIN = 2

OVERW = True
#OVERW = False

# file specification
XMAX = 192
YMAX = 192
ZMAX = 40
NVAR = 2 # Z & VR

DX = 2000.0
DY = 2000.0
x0 = 1000.0
y0 = 1000.0

#CZ_l = np.array([250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 
#       5250, 5750, 6250, 6750, 7250, 7750, 8250, 8750, 9250, 9750, 10250, 10750, 
#    11250, 11750, 12250, 12750, 13250, 13750, 14250, 14750, 15250, 15750, 
#    16250, 16750, 17250, 17750, 18250, 18750, 19250, 19750])
#
#ZMAX_SCALE = len(CZ_l[CZ_l < RADAR_ZMAX])


#DT = 300 # [s] history interval
DT = 300 # [s] history interval

#SIGMA_VR = 5.0 # dBZ
#SIGMA_Z = 3.0  # m/s

# Aksoy et al. (2009MWR)
SIGMA_VR = 2.0 # dBZ
SIGMA_Z = 2.0  # m/s


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


  # Read original NetCDF
  fn = os.path.join(top,exp,citime,"fcst/mean/history.pe000000.nc")
  nc = Dataset(fn, 'r', format='NETCDF4')

  CZ = nc.variables["CZ"][HALO:-HALO]
  CX = nc.variables["CX"][HALO:-HALO]
  CY = nc.variables["CY"][HALO:-HALO]
  times = nc.variables["time"][:]
  nc.close()
 
 
  ZMAX_SCALE = len(CZ[CZ < RADAR_ZMAX])
  GRID = {"ZMAX_SCALE": ZMAX_SCALE, "DX":CX[1]-CX[0], "DY":CY[1]-CY[0], 
          "x0":CX[0], "y0":CY[0], "CZ":CZ}


  # t level in grads
  for it in range(1,len(times)):

    tmp_z3d, tmp_vr3d =  read_grads_ZVR(ofn, it)
 
    z3d = tmp_z3d[:ZMAX_SCALE,:,:]
    vr3d = tmp_vr3d[:ZMAX_SCALE,:,:]

    #idxs = np.where(z3d > MIN_RADAR_REF)
    idxs = np.where(z3d < 100000000.0 ) # dummy
 
    z3d_org = np.copy(z3d)

    ctime = (itime + timedelta(seconds=times[it])).strftime('%Y%m%d%H%M%S')


    rfn = os.path.join(top,exp,citime,typ,mem,"rand_" + ctime + "_" + mem + ".npz")
    rand3d = prep_rand(rfn=rfn, mu=0.0, sigma=1.0, ZMAX_SCALE=GRID["ZMAX_SCALE"])


    if OVERW:
       z3d += rand3d[0,:,:,:] * SIGMA_Z
       vr3d += rand3d[1,:,:,:] * SIGMA_VR

    obs_f = os.path.join(otop,"radar_" + ctime + '.dat')

    write_VR_obs(obs_f, z3d, vr3d, idxs, z3d_org, GRID)


def prep_rand(rfn, mu, sigma, ZMAX_SCALE = -1):


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


  return(z3d, vr3d)

def write_VR_obs(fn,z3d,vr3d,idxs, z3d_org, GRID):
   # idxs is a tuple that contains QC grid information based on nature run (without noise) data
   # idxs is not used (06/10/2019)

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

   #print("Num Obs: ", int(len(kidx) / (THIN * THIN)), fn)

   nobs_z = 0
   nobs_zclr = 0
   nobs_vr = 0


#   for (k,j,i), val in np.ndenumerate(z3d):
   for p in range(len(kidx)):

      # Thinning for clear sky (Aksoy et al. 2009MWR)
      if z3d_org[kidx[p],jidx[p],iidx[p]] == CLR_RADAR_REF: 
         nobs_zclr += 1
         if iidx[p] % THIN != 0 or jidx[p] % THIN != 0:
            continue

      obs_x = GRID["x0"] + GRID["DX"] * iidx[p]
      obs_y = GRID["y0"] + GRID["DY"] * jidx[p]
      obs_z = GRID["CZ"][kidx[p]]

#     Reflectivity obs
#
      obs_seq_ndarray = np.array(tuple([ID_Z,obs_x,obs_y,obs_z,
                                 z3d[kidx[p],jidx[p],iidx[p]],SIGMA_Z,TYP]),
                                 dtype='>f4')
      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
      obs_seq_ndarray.tofile(ofile)
      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
      nobs_z += 1


      if z3d[kidx[p],jidx[p],iidx[p]] > VR_MIN_RADAR_REF:
         obs_seq_ndarray = np.array(tuple([ID_VR,obs_x,obs_y,obs_z,
                                    vr3d[kidx[p],jidx[p],iidx[p]],SIGMA_VR,TYP]),
                                    dtype='>f4')

         ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
         obs_seq_ndarray.tofile(ofile)
         ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))

         nobs_vr += 1

   ofile.close()

   print("No obs (z, vr):", nobs_z, nobs_vr, fn)

#########

time = datetime(2000, 1, 1, 0, 40, 0)
main(time)


