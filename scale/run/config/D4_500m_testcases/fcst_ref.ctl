dset /work/hp150019/share/honda/SCALE-LETKF/AIP_D4_VERIFY/D4_500m_TEST_DEFAULT/20190910130000/dafcst/fcst_ref3d_20190910-130500.grd
options big_endian
undef -9.990000e+33
xdef    256 linear 139.02670288085938 0.0045
ydef    256 linear 35.3876838684082 0.00366
zdef    60 levels   55 165 275 385 495 
                    608.08 727.4924 853.592 986.7532 1127.371 
                    1275.864 1432.673 1598.262 1773.125 1957.78 
                    2152.776 2358.692 2576.139 2805.763 3048.247 
                    3304.309 3574.711 3860.255 4161.79 4480.21 
                    4816.462 5171.544 5546.511 5942.476 6360.614 
                    6802.168 7268.449 7760.842 8280.809 8829.893 
                    9409.727 10022.03 10668.62 11351.42 12072.46 
                   12842.8 13642.8 14442.8 15242.8 16042.8 
                   16842.8 17642.8 18442.8 19242.8 20042.8 
                   20842.8 21642.8 22442.8 23242.8 24042.8 
                   24842.8 25642.8 26442.8 27242.8 28042.8
tdef    60 linear 12:00Z30Jan2018 1mn
* Acutual DT is 30 sec but GrADS cannot deal with DT < 1 min
* FT=0 data, which should be the same with the analysis data, is not included
vars 1
dbz     60 99 X-band Radar Reflectivity (dBZ)
endvars
