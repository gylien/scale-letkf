DSET ^history.grd
TITLE SCALE-RM data output
OPTIONS BIG_ENDIAN
UNDEF -0.99999E+31
XDEF   192 LINEAR    138.32     0.0135
YDEF   192 LINEAR     34.81     0.0108
ZDEF     3 LEVELS  1000.0  3000.0  5000.0
TDEF <--NT-->  LINEAR  <--DATE-->   10MN

VARS    18
U10      0   99 velocity u 10m
V10      0   99 velocity v 10m
T2       0   99 Temperature at 2m 
Q2       0   99 Specific humidity at 2m 
MSLP     0   99 Mean Surface Air Pressure
PREC     0   99 Precipitation
SFC_TEMP 0   99 Surface skin temperature
U       3   99 velocity u 10m
V       3   99 velocity v 10m
W       3   99 velocity w 10m
T       3   99 Temperature
RH      3   99 Relative humidity
QV      3   99 mixing ratio of water vapor
QC      3   99 mixing ratio of ice cloud 
QI      3   99 mixing ratio of liquid cloud
QR      3   99 mixing ratio of rain
QS      3   99 mixing ratio of snow
QG      3   99 mixing ratio of graupel 
 
ENDVARS 
