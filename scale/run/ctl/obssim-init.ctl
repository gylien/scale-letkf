dset /work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_NODA_AKSOY/20010101010500/anal/mean/radar_20010101010500_mean.dat
*dset /work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_DA_PAWR_AKSOY/20010101010500/%e/mean/radar_20010101010500_mean.dat
*dset /work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_InSnd_LT_SN14_Mac_0605_PAWR/20000101004030/%e/mean/radar_20000101004030_mean.dat
options template
edef 2 names anal gues
undef -9.99e33
xdef 192 linear 0 2
ydef 192 linear 0 2
zdef 40 linear 0.25 0.5 
tdef 1000 linear 01:05z01Jan2001 1mn
vars 2
z 40 99 test
vr 40 99 test
endvars