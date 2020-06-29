dset /work/hp150019/share/honda/SCALE-LETKF/AIP_D4_VERIFY/DEBUG/D4_500m_TEST_DEFAULT/20190824150000/dafcst/fcst_all_20190824-150030.grd
options big_endian
undef -9.990000e+33
xdef    256 linear 139.02670288085938 0.0045
ydef    256 linear 35.3876838684082 0.00366
zdef    30 levels   55 275 495 727.4924 986.7532  
                    1275.864 1598.262 1957.78 2358.692 2805.763 
                    3304.309 3860.255 4480.21 5171.544  5942.476  
                    6802.168 7760.842 8829.893 10022.03 11351.42  
                   12842.8 14442.8 16042.8 17642.8 19242.8  
                   20842.8 22442.8 24042.8 25642.8 27242.8 
tdef    61 linear 12:00Z30Jan2018 1mn
* Acutual DT is 30 sec but GrADS cannot deal with DT < 1 min
vars 13
u 30 99 U wind (m/s)
v 30 99 V wind (m/s)
w 30 99 W wind (m/s)
tk 30 99 Temperature (K)
p 30 99 Pressure (Pa)
qv 30 99 Qv (kg/kg)
qc 30 99 Qc (kg/kg)
qr 30 99 Qr (kg/kg)
qi 30 99 Qi (kg/kg)
qs 30 99 Qs (kg/kg)
qg 30 99 Qg (kg/kg)
pw 0 99 Precipitable water (g/m2)
prec 0 99 Precipitation amount (mm/h) snapshot
endvars
