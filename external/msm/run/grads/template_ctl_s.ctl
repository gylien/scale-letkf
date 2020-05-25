dset ^<--GRADSFILE_S-->
options big_endian
title MSMsurf
undef -9.99e+33
xdef 481 linear 120 0.0625
ydef 505 linear 22.4 0.05
zdef 1 linear 0 1
tdef 34 linear <--TIMEF--> 60mn
vars 7
psea=>psea  0  t,y,x  sea level pressure
sp=>sp  0  t,y,x  surface air pressure
u=>u  0  t,y,x  eastward component of wind
v=>v  0  t,y,x  northward component of wind
temp=>temp  0  t,y,x  temperature
rh=>rh  0  t,y,x  relative humidity
r1h=>r1h  0  t,y,x  rainfall in 1 hour
endvars
