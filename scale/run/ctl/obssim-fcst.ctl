*dset /work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_NATURE/20010101010000/fcst/mean/radar_20010101010000_mean.dat
dset /work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_NATURE/20010101010000/fcst/mean/radar3d_i20010101010000_mean_180_180_t003600.dat

options big_endian
undef -9.99e33
xdef 192 linear 0 2
ydef 192 linear 0 2
zdef 40 linear 0.25 0.5
tdef 1000 linear 00:00z01Jan2001 5mn
vars 4
z 40 99 test
randz 40 99 test
vr 40 99 test
randvr 40 99 test
endvars
