import numpy as np
from netCDF4 import Dataset

from datetime import datetime
from datetime import timedelta

import sys
import os

import struct

OVERW = True


SCALE_NP = 1
KHALO = 2

# variable list
# MOMX/Y/Z needs a consideration for a staggered grid system (2018/10/13)
VAR_LIST = ["RHOT","MOMX","MOMY","MOMZ"]

# Horizontal buffer size (grids) where no perturbations are added
HBUF = 2
## Horizontal buffer size (grids) where no perturbations are added
VTBUF = 1 # top
VBBUF = 1 # bottom


RADAR_X = 10000.0 # radar location X [m]
RADAR_Y = 0.0 # radar location Y [m]
RADAR_Z = 0.0 # radar location Z [m]

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

def write_VR_obs(fn,data,xs,ys,zs):

   ID_VR = 4002
  
   ERR = 1.0 # CHECK 
   TYP = 22

   ofile = open(fn,'wb')

   ofile.seek(0)

   obs_seq_ndarray = np.array(tuple([RADAR_X,RADAR_Y,RADAR_Z]),dtype='>f4')
   ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
   obs_seq_ndarray.tofile(ofile)
   ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))

   for idx in range(len(data)):
      obs_seq_ndarray = np.array(tuple([ID_VR,xs[idx],ys[idx],zs[idx],data[idx],ERR,TYP]),dtype='>f4')

      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))
      obs_seq_ndarray.tofile(ofile)
      ofile.write(struct.pack('>i',obs_seq_ndarray.nbytes))

   ofile.close()

def main(itime):

  top = "/data6/honda/SCALE-LETKF/scale_lt_devel_20181002/OUTPUT"
  exp = "8km_sc"
  typ = "fcst"
  m = 1

  obsint = 600 # [s]

  RADAR = True


  QMIN = 0.13 * 1.e-3 #  kg/kg

  # process loop
  for p in range(SCALE_NP):
  #for p in range(1): DEBUG
     print(p)
     fn_nat = nc_name(top,exp,itime,typ,0,p) # nature run (m=0)
     nc_nat = Dataset(fn_nat, "r", format="NETCDF4")


     if p == 0:
       times = nc_nat.variables['time'][:]
       cxg = nc_nat.variables['CXG'][:]
       cyg = nc_nat.variables['CYG'][:]
       czg = nc_nat.variables['CZ'][KHALO:-KHALO]

#     if RADAR:
# 
#
#        dist3d = np.zeros((len(czg),len(cyg),len(cxg))) 
#        for zidx in range(len(czg)):
#           print(zidx)

        #dist3d[0:len(czg),:,:] = czg[:]
#        sys.exit()

     for tidx, time in enumerate(times):
        if int(np.mod(time,obsint)) != 0:
           continue


        # Radar observation
        if RADAR:
           DENS = nc_nat.variables["DENS"][tidx,:,:,:]

           QC = nc_nat.variables["QC"][tidx,:,:,:]
           QR = nc_nat.variables["QR"][tidx,:,:,:]
           QS = nc_nat.variables["QS"][tidx,:,:,:]
           QG = nc_nat.variables["QG"][tidx,:,:,:]

           QHYD = QC + QR + QS + QG

           idxs = np.where(QHYD > QMIN)
           if len(idxs[0]) < 1:
              continue

           ctime = itime + timedelta(seconds=obsint*tidx)
           radar_vrf = "vr_" + ctime.strftime('%Y%m%d%H%M%S') + '.dat'

           zidxs = idxs[0] 
           yidxs = idxs[1] 
           xidxs = idxs[2] 

           U = nc_nat.variables["U"][tidx,:,:,:]
           V = nc_nat.variables["V"][tidx,:,:,:]
           W = nc_nat.variables["W"][tidx,:,:,:]

           dist = np.sqrt(np.square(cxg[xidxs] - RADAR_X) + 
                          np.square(cyg[yidxs] - RADAR_Y) + 
                          np.square(czg[zidxs] - RADAR_Z))

           VR = (cxg[xidxs] / dist * U[zidxs,yidxs,xidxs] +
                 cyg[yidxs] / dist * V[zidxs,yidxs,xidxs] +
                 czg[zidxs] / dist * W[zidxs,yidxs,xidxs])

           print(VR.shape)

           write_VR_obs(radar_vrf,VR,cxg[xidxs],cyg[yidxs],czg[zidxs])


     sys.exit()
     print(times)


#########

time = datetime(2000,1,1,0,0)
main(time)

