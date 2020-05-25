dset ^<--GRADSFILE_P-->
options big_endian
title MSMplev
undef -9.99e+33
xdef 241 linear 120 0.125
ydef 253 linear 22.4 0.1
zdef 16 levels 1000 975 950 925 900 850 800 700
 600 500 400 300 250 200 150 100
tdef 12 linear <--TIMEF--> 180mn
vars 6
z=>z  16  t,z,y,x  geopotential height
w=>w  16  t,z,y,x  vertical velocity in p
u=>u  16  t,z,y,x  eastward component of wind
v=>v  16  t,z,y,x  northward component of wind
temp=>temp  16  t,z,y,x  temperature
rh=>rh  16  t,z,y,x  relative humidity
endvars
