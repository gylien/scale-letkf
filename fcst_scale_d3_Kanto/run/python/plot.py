import numpy as np
#from datetime import datetime
#from datetime import timedelta
import os
import sys

from plot_common import get_dims, prep_proj_multi, read_ens2d

import matplotlib.pyplot as plt




dir = "/work/hp150019/f22013/SCALE-LETKF/scale-5.3.2/OUTPUT/TEST_exp_d2/20190130000000/fcst"

fn = os.path.join(dir, "history_0001.pe000000.nc")

#nvar = "U"
nvar = "RAIN"
fac = 3600.0
nens = 2
zlev = 10
tlev = 2

#mem_l = []
#for m in range(1,nens+1):
#    mem_l.append('{:04d}'.format(m))

mem_l = ["0001", "0002"]

lon_ll = 120.0
lon_ur = 150.0
lat_ur = 45.0
lat_ll = 25.0

blon = 135.0
blat = 35.0
lat2 = 40.0

PLOT_DIMS = {"nvar":nvar, "mem_l":mem_l, "zlev":zlev, "tlev":tlev, 
             "lon_ll":lon_ll, "lon_ur":lon_ur, 
             "lat_ll":lat_ll, "lat_ur":lat_ur, 
             "blon":blon, "blat":blat, "lat2":lat2, 
             }


DIMS = get_dims(fn, dir)
evar = read_ens2d(DIMS, dir, PLOT_DIMS)


#fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2, 3, figsize=(13.0, 8.0))
fig, (ax1) = plt.subplots(1, 1, figsize=(10.0, 8.0))
fig.subplots_adjust(left=0.06, bottom=0.03, right=0.91, top=0.89, wspace=0.2, hspace=0.1)

ax_l = [ax1]
m_l = prep_proj_multi(ax_l, PLOT_DIMS)

var = np.mean(evar[:,:,:], axis=0)

x1, y1 = m_l[0](DIMS["lons"], DIMS["lats"])

#cmap = plt.cmap.jet
cmap = "jet"

SHADE1 = ax1.contourf(x1, y1, var*fac, cmap=cmap)


# color bar
pos = ax1.get_position()
cb_h = pos.height
ax_cb = fig.add_axes([pos.x1+0.01, pos.y0, 0.015, cb_h])

CB1 = plt.colorbar(SHADE1, cax=ax_cb, orientation = 'vertical')

plt.show()

