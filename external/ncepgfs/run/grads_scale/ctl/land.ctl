dset ^land_%y4%m2%d2%h2%n200.grd
undef 9.999E+20
options template big_endian
xdef 720 linear   0.000000 0.500000
ydef 361 linear -90.000000 0.500000
zdef 4 levels 0.05 0.25 0.70 1.50
*Levels of soil layer 0-0.1m, 0.1-0.4m, 0.4-1m, 1-2m
tdef 21 linear <STIME> 6hr
VARS 4
LANDsfc  0 99 surface Land Cover (1=land, 0=sea) [Proportion]
TMPsfc   0 99 surface Temperature [K]
TSOIL    4 99 ground Temperature [K]
SOILW    4 99 ground Volumetric Soil Moisture Content [Fraction]
ENDVARS
