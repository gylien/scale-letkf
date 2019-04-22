import numpy as np
from datetime import datetime
#from datetime import timedelta
import os
import sys

from plot_common import get_gdims, prep_proj_multi, read_ens2d_split, def_cmap

import matplotlib.pyplot as plt

#
# Plot 3 panel figure by reading split SCALE NetCDF files
#

def main(nvar_l, INFO, itime, tlev_l, zlev, quick):
  
    if INFO["SCALE_NP"] < INFO["SCALE_NP_ORG"]:
      dir = os.path.join(INFO["top"], INFO["exp"], INFO["ctime"], 
                        "fcst_sno_np" + str(INFO["SCALE_NP"]).zfill(5))
    else:
      dir = os.path.join(INFO["top"], INFO["exp"], INFO["ctime"], 
                        "fcst")
  
    # Display domain setting   
    lon_ll = 125.0
    lon_ur = 145.0
    lat_ur = 44.0
    lat_ll = 26.0
    
    # Display domain setting (Map projection)
    blon = 135.0
    blat = 35.0
    lat2 = 40.0
    
    PLOT_DIMS = {"mem_l":INFO["mem_l"], "zlev":zlev, "tlev":tlev_l[0], 
                 "lon_ll":lon_ll, "lon_ur":lon_ur, 
                 "lat_ll":lat_ll, "lat_ur":lat_ur, 
                 "blon":blon, "blat":blat, "lat2":lat2, 
                 "dir":dir, "nvar":nvar_l[0],
                 }
    
    
    GDIMS, PLOT_DIMS = get_gdims(scale_np, PLOT_DIMS)
    
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
    

    for nidx, nvar in enumerate(nvar_l):
      for tidx, tlev in enumerate(tlev_l):
 
          # Initiate figure
          fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12.0, 5.0))
          fig.subplots_adjust(left=0.06, bottom=0.03, right=0.98, top=0.93, 
                              wspace=0.2, hspace=0.1)
                  
          ax_l = [ax1, ax2, ax3]

          # Set up map projection
          m_l = prep_proj_multi(ax_l, PLOT_DIMS, quick)
          m_l_ = np.copy(m_l)

          # Get figure coordinate   
          x1, y1 = m_l[0](GDIMS["lon2d"], GDIMS["lat2d"])
 
          # Read ensemble forecast data (t=tlev) and store in a 4-d array evar(e,z,y,x)
          PLOT_DIMS["tlev"] = tlev
          evar = read_ens2d_split(nvar, GDIMS, PLOT_DIMS)
          
          fac = 1.0
      
          # Set figure name and unit
          unit = "(" + PLOT_DIMS["unit"] + ")"
          if nvar == "CAPE":
            tvar = nvar 
            unit = r'(J kg$^{-1}$)'
          elif nvar == "U" or nvar == "V":
            tvar = nvar
          elif nvar == "RAIN":
            tvar = nvar 
            if tlev > 0:
              fac = GDIMS["time"][tlev] - GDIMS["time"][tlev-1] 
            else: 
              print("!Unit for RAIN is not accurate! ",tlev)
              fac = 3600
            unit = "(mm/" + str(fac/3600) + "h)"
          else:
            tvar = nvar
        
          # Get color map and shading invervals
          # cmap_f, cnorm_f, & levs: For "nvar" itself
          # cmap_s, cnorm_s, & levs_s: For ensemble spread of "nvar" 
          cmap_f, cmap_s, cnorm_f, cnorm_s, levs, levs_s = def_cmap(nvar, np.max(np.abs(evar))*fac)
      
          # Set panel title
          tit_l = ["Ensemble mean\n" + tvar, 
                   "Ensemble spread\n" + tvar, 
                   "Maximum\n" + tvar]
          num_l = ["(a)","(b)","(c)"]
          # Set figure title
          tit = "D2 ensemble forecast"
        
          for idx, ax in enumerate(ax_l):
          
            levels = levs
            if idx == 0:
               # Get ensemble mean
               var = np.mean(evar[:,:,:]*fac, axis=0)
               cmap = cmap_f
               cnorm = cnorm_f
            elif idx == 1:
               # Get ensemble spread
               var = np.std(evar[:,:,:]*fac, axis=0, ddof=1)
               levels = levs_s
               cmap = cmap_s
               cnorm = cnorm_s
            elif idx == 2:
               # Get maximum among the ensemble
               var = evar.max(axis=0)*fac
               cmap = cmap_f
               cnorm = cnorm_f
        
            # Draw color shading
            SHADE1 = m_l[idx].contourf(x1, y1, var, cmap=cmap, 
                                 levels=levels, norm=cnorm, 
                                 extend='both')
          
          
            # color bar
            pos = ax.get_position()
            cb_h = pos.height
            ax_cb = fig.add_axes([pos.x0+0.01, pos.y0-0.06, pos.width-0.01, 0.02])
          
            CB1 = plt.colorbar(SHADE1, cax=ax_cb, orientation = 'horizontal')
            CB1.ax.tick_params(labelsize=8)    
      
            # Draw panel title
            ax.set_title(tit_l[idx], size=14, loc = 'center')
      
            # Draw panel caption (e.g., (a), (b)...)
            ax.text(0.01, 1.03, num_l[idx],
                    verticalalignment='bottom', horizontalalignment='left',
                    transform=ax.transAxes, color='k', fontsize=12)
        
            # Draw nvar unit
            ax.text(1.00, 1.01, unit,
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=ax.transAxes, color='k', fontsize=12)
      
          # Draw figure title
          fig.suptitle(tit, fontsize=18)
        
        
          # Set time information
          time_info = "Init: " + itime.strftime('%m/%d/%Y %H:%M:%S UTC') + \
                      "\nFT=" + "{0:.1f}".format(GDIMS["time"][tlev]/3600) + "h"
           
          # Set a footer of figure name
          foot = "_FT" + str(int(GDIMS["time"][tlev])).zfill(6) + "s_i" + \
                 itime.strftime('%Y%m%d%H%M%S') 
      
          if nvar != "CAPE" and nvar != "RAIN":
             time_info += "\nZ=" + "{0:.1f}".format(GDIMS["z"][zlev]/1000) + "km" 
             foot += "z" + "{0:.1f}".format(GDIMS["z"][zlev]/1000) + "km"
      
          # Draw time information
          fig.text(0.99, 0.97, time_info,
                    verticalalignment='top', horizontalalignment='right',
                    color='k', fontsize=10)
        
          ofig = tvar + foot
      
          if not quick: # Store figure
            png_dir = os.path.join(INFO["run_dir"] + "/png",itime.strftime('%Y%m%d'))
            os.makedirs(png_dir, exist_ok=True)
            plt.savefig(os.path.join(png_dir, ofig + ".png"), 
                        bbox_inches="tight", pad_inches = 0.1)
            plt.clf()
            plt.close('all')
          else: # quick view
            print(ofig)
            plt.show()
      
     

###################

# OUTDIR
top = "/work/hp150019/f22013/SCALE-LETKF/scale-5.3.2/OUTPUT"

# Run directory
run_dir = "/work/hp150019/f22013/SCALE-LETKF/scale-5.3.2/scale-5.3.x_LETKF/scale-letkf-SNO/scale/run"

# exp name
exp = "TEST_exp_d2"

# Original SCALE process number
scale_np_org = 256

# SCALE process number after SNO
scale_np = 4

# Forecast initial time (yyyy, mm, dd, hh, nn, ss)
itime = datetime(2019, 1, 30, 0, 0)

# Variable list
nvar_l = ["CAPE", "RAIN"]
#nvar_l = ["RAIN"]

# t level (tlev=0 means FT=0)
tlev = 1
tlev_l = np.arange(1,18+1)
#tlev_l = [1,2]

# z level 
zlev = 30

# Conduct a quick view
quick = True
quick = False

# Member list
mem_l = []

# ensemble size
mems = 8 

# Prepare ensemble list
for m in range(mems):
    mem = str(m+1).zfill(4)
    mem_l.append(mem)
#mem_l.append("mean")
#mem_l.append("mdet")

INFO = {"top":top, "exp":exp, "ctime":itime.strftime('%Y%m%d%H%M%S'),
        "run_dir":run_dir, 
        "SCALE_NP":scale_np, "SCALE_NP_ORG":scale_np_org, "mem_l":mem_l}

# loop for nvar
main(nvar_l, INFO, itime, tlev_l, zlev, quick)
