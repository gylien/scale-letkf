import numpy as np
from netCDF4 import Dataset

fn = "../tmp_obssim/out/test03km/fcst/Him8_20160601120000_0001_FT01.nc"
fn = "/work/hp150019/share/honda/LMU/experiments/germany_3km_15h_en_fcst_1000m/20160601120000/obssim/Him8_20160601120000_0001_FT01.nc"

nc = Dataset( fn, "r", format="NETCDF4" ) 
lon2d = nc.variables['longitude'][:]
lat2d = nc.variables['latitude'][:]
tbb3d = nc.variables['tbb'][:]
nc.close()

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

print(tbb3d.shape)
print(lon2d.shape, lat2d.shape)
print(np.max(lon2d), np.min(lon2d))
print(np.max(lat2d), np.min(lat2d))

levs_tbb = np.arange(200, 280, 5)
levs_refl = np.arange(0, 1.05, 0.05)

cmap_refl = plt.cm.get_cmap("gist_gray")
cmap_tbb = plt.cm.get_cmap("gist_gray_r")

b = 0
for b in range(16):

    if b <= 5:
       levs = levs_refl
       cmap = cmap_refl
    else:
       levs = levs_tbb
       cmap = cmap_tbb

    norm = BoundaryNorm(levs, ncolors=cmap.N, clip=True)

    fig, (ax1) = plt.subplots(1, 1, figsize=(9,9))
    SHADE = ax1.pcolormesh(lon2d, lat2d,
                           #lon2d[:,:],
                           tbb3d[b,:,:],
                           vmin=np.min(levs),
                           vmax=np.max(levs),
                           cmap=cmap, norm=norm,
                           )
 
    pos = ax1.get_position()
    cb_h = pos.height
    ax_cb = fig.add_axes( [pos.x1+0.005, pos.y0, 0.005, 0.7] )
    cb = plt.colorbar( SHADE, cax=ax_cb, extend='both' )

    tit = "Band: " + str(b+1).zfill(2)
    ax1.set_title( tit, size=12, loc='center')

    plt.show()
