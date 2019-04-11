import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from datetime import timedelta

import matplotlib.colors as mpc

import os
import sys

def def_cmap(nvar, vmax):

    cmap_s = mpc.ListedColormap(['cyan','dodgerblue',
                                 'lime', 'limegreen','yellow',
                                 'orange', 'red', 'firebrick', 'magenta',
                                 'purple'])
    cmap_s.set_under('w', alpha=0.1)

    if nvar == "QHYD" or nvar == "VORT":
      levs = np.array([0.1, 0.2, 0.4, 0.8, 1.6, 2.4, 3.2, 4.0])
      levs_s = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
      cmap = mpc.ListedColormap(['cyan','dodgerblue',
                                 'lime', 'yellow',
                                 'orange', 'red', 'magenta',
                                 'purple'])

    elif nvar == "DBZ":
      levs = np.array([15, 20, 25, 30, 35, 40, 45, 50, 55, 60])
      levs_s = np.array([5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
      cmap = mpc.ListedColormap(['cyan','dodgerblue',
                                 'lime', 'limegreen','yellow',
                                 'orange', 'red', 'firebrick', 'magenta',
                                 'purple'])

    elif nvar == "CAPE":
      levs = np.array([50, 100, 150, 200, 400, 800, 1200, 1600, 2000, 2400])
      levs_s = np.array([10, 20, 40, 80, 120, 160, 200, 300, 400, 500])
      cmap = mpc.ListedColormap(['cyan','dodgerblue',
                                 'lime', 'limegreen','yellow',
                                 'orange', 'red', 'firebrick', 'magenta',
                                 'purple'])

    elif nvar == "RAIN":
      levs = np.array([0.5, 1, 2, 4, 8, 12, 16, 20, 30, 40])
      #levs = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
      levs_s = np.array([0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9])

      #rat = int(vmax / 10)
      #
      #levs = levs * rat

      cmap = mpc.ListedColormap(['cyan','dodgerblue',
                                 'lime', 'limegreen','yellow',
                                 'orange', 'red', 'firebrick', 'magenta',
                                 'purple'])

    elif nvar == "U" or nvar == "V":
      levs = np.array([-20, -16, -12, -8, -4, 4, 8, 12, 16, 20])
      levs_s = np.array([1, 2, 3, 4, 6, 8, 10, 12, 16, 20])

      rat = int(vmax / 20)

      levs = levs * rat

      cmap = mpc.ListedColormap(['darkblue', 'cornflowerblue','lightskyblue','lightcyan',
                                 'w',
                                 'mistyrose','salmon','firebrick','maroon',
                                 ])
      cmap.set_under('gray', alpha=1.0)
      cmap.set_over('gray', alpha=1.0)

    if nvar != "U" and nvar != "V":
      cmap.set_under('w', alpha=0.1)

    cnorm = mpc.BoundaryNorm(levs, ncolors=len(levs)-1, clip=False)
    cnorm_s = mpc.BoundaryNorm(levs_s, ncolors=len(levs_s)-1, clip=False)

    return(cmap, cmap_s, cnorm, cnorm_s, levs, levs_s)

def readlatlon_split(PLOT_DIMS, GDIMS):
 
    # t,z,y,x
    shape2d = (GDIMS["NYG"], GDIMS["NXG"])

    lon2d = np.zeros(shape2d)
    lat2d = np.zeros(shape2d)

    lon2d = read2d_split("lon", os.path.join(PLOT_DIMS["dir"], PLOT_DIMS["mem_l"][0]),
                         GDIMS, PLOT_DIMS)
    lat2d = read2d_split("lat", os.path.join(PLOT_DIMS["dir"], PLOT_DIMS["mem_l"][0]),
                         GDIMS, PLOT_DIMS)

    return(lon2d, lat2d)
    
def read2d_split(nvar, dir, GDIMS, PLOT_DIMS):
 
    #shape2d = (PLOT_DIMS["jmax_g"] - PLOT_DIMS["jmin_g"], 
    #          PLOT_DIMS["imax_g"] - PLOT_DIMS["imin_g"])
    shape2d = (GDIMS["NYG"], GDIMS["NXG"])
    
    var = np.zeros(shape2d)

    for p in range(GDIMS["SCALE_NP"]):
      fn = dir + "/history.pe" + str(p).zfill(6) + ".nc"

      nc = Dataset(fn, 'r', format='NETCDF4')

      node_x = nc.scale_cartesC_prc_rank_x
      node_y = nc.scale_cartesC_prc_rank_y

      imin = node_x * len(GDIMS["x"])
      imax = imin + len(GDIMS["x"])

      jmin = node_y * len(GDIMS["y"])
      jmax = jmin + len(GDIMS["y"])

      if nvar == "lon" or nvar == "lat":
        var[jmin:jmax,imin:imax] = nc.variables[nvar][:,:]
      elif nvar == "RAIN" or nvar == "CAPE":
        var[jmin:jmax,imin:imax] = nc.variables[nvar][\
                                                PLOT_DIMS["tlev"],:,:]
      else:
        var[jmin:jmax,imin:imax] = nc.variables[nvar][\
                                                PLOT_DIMS["tlev"],
                                                PLOT_DIMS["zlev"],:,:]

    return(var) 

def read_ens2d_split(nvar, GDIMS, PLOT_DIMS):
 
    shape = (len(PLOT_DIMS["mem_l"]),GDIMS["NXG"], GDIMS["NYG"])

    evar = np.zeros(shape)

    for midx, mem in enumerate(PLOT_DIMS["mem_l"]):
       mem_dir = os.path.join(PLOT_DIMS["dir"], mem)

       evar[midx,:,:] = read2d_split(nvar, mem_dir, GDIMS, PLOT_DIMS)

    return(evar)



def get_gdims(scale_np, PLOT_DIMS):

    fn = os.path.join(PLOT_DIMS["dir"], PLOT_DIMS["mem_l"][0], "history.pe000000.nc")

    nc = Dataset(fn, 'r', format='NETCDF4')
    cxg = nc.variables["CXG"][:]
    cyg = nc.variables["CYG"][:]
    z = nc.variables["z"][:]
    x = nc.variables["x"][:]
    y = nc.variables["y"][:]
    time = nc.variables["time"][:]
 
    khalo = nc.scale_atmos_grid_cartesC_index_khalo
    jhalo = nc.scale_atmos_grid_cartesC_index_jhalo
    ihalo = nc.scale_atmos_grid_cartesC_index_ihalo

    unit = nc.variables[PLOT_DIMS["nvar"]].units
  
    nc.close()

    GDIMS = {"CXG":cxg, "CYG":cyg, "z":z, "x":x, "y":y, "time":time,
             "NXG":len(cxg) - 2 * ihalo, "NYG":len(cyg) - 2 * jhalo,
             "KHALO":khalo, "JHALO":jhalo, "IHALO":ihalo,
             "SCALE_NP":scale_np}

    if scale_np > 0:
       lon2d, lat2d = readlatlon_split(PLOT_DIMS, GDIMS)

    GDIMS["lon2d"] = lon2d
    GDIMS["lat2d"] = lat2d
    PLOT_DIMS["unit"] = unit

    return(GDIMS, PLOT_DIMS)

def get_dims(fn, dir):

    nc = Dataset(fn, 'r', format='NETCDF4')

    lon2d = nc.variables["lon"][:]
    lat2d = nc.variables["lat"][:]
    CZ = nc.variables["CZ"][:]
    time = nc.variables["time"][:]
    
    nc.close()

    DIMS = {"lon2d":lon2d, "lat2d":lat2d, "CZ":CZ, "times":time, "dir":dir}

    return(DIMS)

def prep_proj_multi(axs, PLOT_DIMS, quick):

    print("proj")
    from mpl_toolkits.basemap import Basemap
    ddeg = 0.2 # deg

    m_l = []
    cc = "k"
    contc = 'lightgrey'

    pdlat = 5.0
    pdlon = 5.0

    lc = 'k'
    fs = 8
    lw = 0.2

    if quick:
       res = 'c'
    else:
       res = 'i'

    for idx, ax in enumerate(axs):
       m = Basemap(projection="mill",resolution = res,
               llcrnrlon=PLOT_DIMS["lon_ll"], llcrnrlat=PLOT_DIMS["lat_ll"],
               urcrnrlon=PLOT_DIMS["lon_ur"], urcrnrlat=PLOT_DIMS["lat_ur"],
               lat_0=PLOT_DIMS["blat"], lat_1=PLOT_DIMS["lat2"],
               lat_2=PLOT_DIMS["lat2"], lon_0=PLOT_DIMS["blon"],
               ax = ax)
       m_l.append(m)
 
       m.drawcoastlines(linewidth = 0.5, color = cc, zorder=2)
       m.fillcontinents(color=contc,lake_color='w', zorder=0)
       m.drawparallels(np.arange(0,70,pdlat),labels=[1,0,0,0],fontsize=fs,color=lc,linewidth=lw)
       m.drawmeridians(np.arange(0,180,pdlon),labels=[0,0,0,1],fontsize=fs,color=lc,linewidth=lw)

    return(m_l)


def read_ens2d(DIMS, dir, PLOT_DIMS):
 
    odim = (len(PLOT_DIMS["mem_l"]), len(DIMS["lat2d"]), len(DIMS["lon2d"]))

    ovar = np.zeros(odim)

    for midx, mem in enumerate(PLOT_DIMS["mem_l"]):

       fn = os.path.join(DIMS["dir"], "history_" + mem + ".pe000000.nc")
       
       nc = Dataset(fn, 'r', format='NETCDF4')
       if nvar == "RAIN":
          ovar[midx,:,:] = nc.variables[nvar][PLOT_DIMS["tlev"]-1,:,:]
       else:
          ovar[midx,:,:] = nc.variables[nvar][PLOT_DIMS["tlev"]-1,PLOT_DIMS["zlev"]-1,:,:]
     
       nc.close()

    return(ovar)

