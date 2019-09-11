*dset /work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_InSnd_LT_SN14_Mac_0605_PAWR/20000101%h2%n200/%e/sprd/init_grads.dat
*dset /work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_InSnd_LT_SN14_Mac_0605_NODA/20000101%h2%n230/%e/sprd/init_grads.dat

dset /work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_NODA_AKSOY/2001010101%n200/%e/sprd/init_grads.dat

options template
edef 2 names anal gues
undef -9.99e33
xdef 192 linear 0 2
ydef 192 linear 0 2
zdef 40 linear 0.25 0.5 
tdef 1000 linear 01:01z01Jan2001 1mn
vars 17
u 40 99 U wind
v 40 99 V wind
W 40 99 W wind
tk 40 99 Temperature
p 40 99 Pressure
qv 40 99 Water vapor mixing ratio
qc 40 99 Cloud water mixing ratio
qr 40 99 Rain water mixing ratio
qi 40 99 Cloud ice mixing ratio
qs 40 99 Snow mixing ratio
qg 40 99 Graupel mixing ratio
cc 40 99 Charge density (cloud)
cr 40 99 Charge density (rain)
ci 40 99 Charge density (ice)
cs 40 99 Charge density (snow)
cg 40 99 Charge density (graupel)
oep 40 99 Old electric potential 
endvars
