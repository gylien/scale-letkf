import numpy as np

from netCDF4 import Dataset
from datetime import datetime
from datetime import timedelta

import os

def read_var(fn, nvar):
   print(fn)

   nc = Dataset(fn, "r", format="NETCDF4")

   if nvar == "QH":
     var3d = nc.variables["QC"][:]
     var3d += nc.variables["QR"][:]
     var3d += nc.variables["QI"][:]
     var3d += nc.variables["QS"][:]
     var3d += nc.variables["QG"][:]
   else:
     var3d = nc.variables[nvar][:]

   x = nc.variables["x"][:] * 0.001
   y = nc.variables["y"][:] * 0.001

   import matplotlib.pyplot as plt
   fig, ax1= plt.subplots(1, 1, figsize=(8.5,7.5))

   print(var3d.shape)

   imin = 30
   imax = -70

   jmin = 30
   jmax = -70

   #var3d[np.abs(var3d) < 1.5] = np.nan
   SHADE = ax1.contourf(x[imin:imax],y[jmin:jmax],var3d[imin:imax,jmin:jmax,10], cmap="jet")
   fig.colorbar(SHADE)

   plt.show()

   return(var3d)


top = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT"
exp = "2000m_InSnd_LT_SN14_Mac_0523_NODA"
exp = "2000m_InSnd_LT_SN14_Mac_0523_DA"
#exp = "2000m_InSnd_LT_SN14_Mac_0523"

time = datetime(2000, 1, 1, 0, 0, 0)
time = datetime(2000, 1, 1, 0, 40, 30)

typ = "anal"

#fn = os.path.join(top, exp, time.strftime('%Y%m%d%H%M%S'), typ, "init_mean.nc")
##nvar = "MOMZ"
#MOMY_a = read_var(fn, nvar)

nvar = "QH"
fn = os.path.join(top, exp, time.strftime('%Y%m%d%H%M%S'), "gues", "init_mean.nc")
MOMY_g = read_var(fn, nvar)

#print(np.max(MOMY_a - MOMY_g))


