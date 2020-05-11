dset /work/hp150019/share/honda/SCALE-LETKF/AIP_D4_VERIFY/DEBUG/D4_500m_TEST_DEFAULT/20190824150000/dafcst/fcst_all_20190824-150030.grd
options big_endian
undef -9.990000e+33
xdef    256 linear 139.02670288085938 0.0045
ydef    256 linear 35.3876838684082 0.00366
zdef    30 levels    110.0000    330.0000    550.0000    788.8249   1055.1472 
           1352.1328   1683.3123   2052.6226   2464.4541   2923.7026  
           3435.8274   4006.9160   4643.7568   5353.9209   6145.8501  
           7028.9585   8013.7441   9111.9131  10336.5205  11702.1240 
          13242.8018  14842.8018  16442.8008  18042.8008  19642.8008 
          21242.8008  22842.8008  24442.8008  26042.8008  27642.8008 
tdef    61 linear 12:00Z30Jan2018 1mn
* Acutual DT is 30 sec but GrADS cannot deal with DT < 1 min
vars 12
u 30 99 U wind
v 30 99 V wind
w 30 99 W wind
tk 30 99 Temperature
p 30 99 Pressure
qv 30 99 Qv
qc 30 99 Qc
qr 30 99 Qr
qi 30 99 Qi
qs 30 99 Qs
qg 30 99 Qg
pw 0 99 Precipitable water
endvars
