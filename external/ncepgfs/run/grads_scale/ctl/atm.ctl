dset ^atm_%y4%m2%d2%h2%n200.grd
undef 9.999E+20
options template big_endian
xdef 720 linear    0.000000 0.500000
ydef 361 linear  -90.000000 0.500000
zdef 31  levels
 1000 975 950 925 900 850 800 750 700 650
  600 550 500 450 400 350 300 250 200 150
  100  70  50  30  20  10   7   5   3   2
    1
tdef 21 linear <STIME> 6hr
VARS 5
HGTprs    31 99  Geopotential Height [gpm]
UGRDprs   31 99  U-Component of Wind [m/s]
VGRDprs   31 99  V-Component of Wind [m/s]
TMPprs    31 99  Temperature [K]
RHprs     31 99  Relative Humidity [%]
ENDVARS
