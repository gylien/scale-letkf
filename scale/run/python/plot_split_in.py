import numpy as np
from datetime import datetime
#from datetime import timedelta
import os
import sys

from plot_common import get_gdims, prep_proj_multi, read_ens2d_split, def_cmap

import matplotlib.pyplot as plt


def main(nvar, INFO, itime, levs, quick):
  
    if INFO["SCALE_NP"] < INFO["SCALE_NP_ORG"]:
      dir = os.path.join(INFO["top"], INFO["exp"], INFO["ctime"], 
                        "fcst_sno_np" + str(INFO["SCALE_NP"]).zfill(5))
    else:
      dir = os.path.join(INFO["top"], INFO["exp"], INFO["ctime"], 
                        "fcst")
  
    
    lon_ll = 125.0
    lon_ur = 145.0
    lat_ur = 44.0
    lat_ll = 26.0
    
    blon = 135.0
    blat = 35.0
    lat2 = 40.0
    
    PLOT_DIMS = {"mem_l":INFO["mem_l"], "zlev":zlev, "tlev":tlev, 
                 "lon_ll":lon_ll, "lon_ur":lon_ur, 
                 "lat_ll":lat_ll, "lat_ur":lat_ur, 
                 "blon":blon, "blat":blat, "lat2":lat2, 
                 "dir":dir,
                 }
    
    
    GDIMS = get_gdims(scale_np, PLOT_DIMS)
    
  #  imin_g = 0
  #  imax_g = len(GDIMS["CXG"]) - GDIMS["IHALO"] * 2
  #  
  #  jmin_g = 0
  #  jmax_g = len(GDIMS["CYG"]) - GDIMS["JHALO"] * 2
  #  
  #  PLOT_DIMS["imin_g"] = imin_g
  #  PLOT_DIMS["jmin_g"] = jmin_g
  #  PLOT_DIMS["imax_g"] = imax_g
  #  PLOT_DIMS["jmax_g"] = jmax_g
    
    PLOT_DIMS["tlev"] = levs[0]
    PLOT_DIMS["zlev"] = levs[1]
    #U = read2d_split(dir, GDIMS, PLOT_DIMS)
    
    evar = read_ens2d_split(nvar, GDIMS, PLOT_DIMS)
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12.0, 5.0))
    fig.subplots_adjust(left=0.06, bottom=0.03, right=0.98, top=0.93, 
                        wspace=0.2, hspace=0.1)
    
    ax_l = [ax1, ax2, ax3]
    m_l = prep_proj_multi(ax_l, PLOT_DIMS, quick)
    
    x1, y1 = m_l[0](GDIMS["lon2d"], GDIMS["lat2d"])
    
    #cmap = plt.cmap.jet
    cmap = "jet"
  
    fac = 1.0

    if nvar == "CAPE":
      tvar = nvar 
      unit = r'(J kg$^{-1}$)'
    elif nvar == "U" or nvar == "V":
      tvar = nvar
    elif nvar == "RAIN":
      tvar = nvar 
      unit = "(mm)"
      fac = 3600
    else:
      tvar = nvar
  
    cmap_f, cmap_s, cnorm_f, cnorm_s, levs, levs_s = def_cmap(nvar, np.max(np.abs(evar)))

    tit_l = ["Ensemble mean\n" + tvar, 
             "Ensemble spread\n" + tvar, 
             "Maximum\n" + tvar]
    num_l = ["(a)","(b)","(c)"]
    tit = "D2 ensemble forecast"
  
    for idx, ax in enumerate(ax_l):
    
      levels = levs

      if idx == 0:
         var = np.mean(evar[:,:,:]*fac, axis=0)
         cmap = cmap_f
         cnorm = cnorm_f
      elif idx == 1:
         var = np.std(evar[:,:,:]*fac, axis=0, ddof=1)
         levels = levs_s
         cmap = cmap_s
         cnorm = cnorm_s
      elif idx == 2:
         var = evar.max(axis=0)*fac
         cmap = cmap_f
         cnorm = cnorm_f
  
      SHADE1 = ax.contourf(x1, y1, var, cmap=cmap, 
                           levels=levels, norm=cnorm, 
                           extend='both')
    
    
    # color bar
      pos = ax.get_position()
      cb_h = pos.height
      #ax_cb = fig.add_axes([pos.x1+0.01, pos.y0, 0.015, cb_h])
      ax_cb = fig.add_axes([pos.x0+0.01, pos.y0-0.05, pos.width-0.01, 0.02])
    
      CB1 = plt.colorbar(SHADE1, cax=ax_cb, orientation = 'horizontal')
      CB1.ax.tick_params(labelsize=8)    

      ax.set_title(tit_l[idx], size=14, loc = 'center')
      ax.text(0.01, 1.03, num_l[idx],
              verticalalignment='bottom', horizontalalignment='left',
              transform=ax.transAxes, color='k', fontsize=12)
  
      ax.text(1.00, 1.01, unit,
              verticalalignment='bottom', horizontalalignment='right',
              transform=ax.transAxes, color='k', fontsize=12)

    fig.suptitle(tit, fontsize=18)
  
  
    time_info = "Init: " + itime.strftime('%m/%d/%Y %H:%M:%S UTC') + \
                "\nFT=" + "{0:.1f}".format(GDIMS["time"][tlev]/3600) + "h"
     
    if nvar != "CAPE" and nvar != "RAIN":
       time_info += "\nZ=" + "{0:.1f}".format(GDIMS["z"][zlev]/1000) + "km" 
 

    fig.text(0.99, 0.97, time_info,
              verticalalignment='top', horizontalalignment='right',
              color='k', fontsize=10)
  
    ofig = "000TEST.png"

    if not quick:
#      os.makedirs(p, exist_ok=True)
      print(ofig)
      plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
      plt.clf()
      plt.close('all')
    else:
      print(ofig)
      plt.show()

  

###################

top = "/work/hp150019/f22013/SCALE-LETKF/scale-5.3.2/OUTPUT"
exp = "TEST_exp_d2"
ctime = "20190130000000"
scale_np_org = 256
scale_np = 4

itime = datetime(2019, 1, 30, 0, 0)

nvar = "CAPE"
#nvar = "U"
#nvar = "RAIN"

tlev = 1
zlev = 20
levs = [tlev, zlev]

quick = True
#quick = False

mem_l = []

mems = 8
for m in range(mems):
    mem = str(m+1).zfill(4)
    mem_l.append(mem)

INFO = {"top":top, "exp":exp, "ctime":itime.strftime('%Y%m%d%H%M%S'), 
        "SCALE_NP":scale_np, "SCALE_NP_ORG":scale_np_org, "mem_l":mem_l}
main(nvar, INFO, itime, levs, quick)
