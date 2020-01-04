dset ^sfc_%y4%m2%d2%h2%n200.grd
undef 9.999E+20
options template big_endian
xdef 720 linear   0.000000 0.500000
ydef 361 linear -90.000000 0.500000
zdef 4 levels 10 40 100 200
tdef 21 linear <STIME> 6hr
VARS 7
MSLETmsl  0 99 mean sea level Pressure Reduced to MSL [Pa]
PRESsfc   0 99 surface Pressure [Pa]
UGRD10m   0 99 10 m above ground U-Component of Wind [m/s]
VGRD10m   0 99 10 m above ground V-Component of Wind [m/s]
TMP2m     0 99 2 m above ground Temperature [K]
RH2m      0 99 2 m above ground Relative Humidity [%]
HGTsfc    0 99 surface Geopotential Height [gpm]
ENDVARS
