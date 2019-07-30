*dset ^out/test03km/anal/Him8_20160603120000_mean.dat 
dset /volume64/data/ra001011/honda/SCALE-LETKF/VIS_TEST/OUTPUT/test03km/20160603120000/obssim/Him8_20160603120000_%e_FT01.dat 
options big_endian template
*undef -9.990000e+33
undef -9.990000e+33
xdef    499 linear     4.111291     0.023649
ydef    703 linear    45.525261     0.013730
zdef      1 linear 1 1 * B01, B02, B03
edef 5 names mean 0001 0002 0003 0004
tdef     1 linear 13:00Z03Jun2016 60mn
pdef    250    352 lcc    50.400000    10.000000 125.5 176.5    43.000000    57.000000    10.000000  3000.000000  3000.000000
vars 16
refl01     0 99 Reflectance (0.47, 0.51, 0.64 micro meters)
refl02     0 99 Reflectance (0.47, 0.51, 0.64 micro meters)
refl03     0 99 Reflectance (0.47, 0.51, 0.64 micro meters)
refl04     0 99 Reflectance (0.47, 0.51, 0.64 micro meters)
refl05     0 99 Reflectance (0.47, 0.51, 0.64 micro meters)
refl06     0 99 Reflectance (0.47, 0.51, 0.64 micro meters)
tbb07      0 99 brightness temperature (Band 07 of Advanced Himawari Imager)
tbb08      0 99 brightness temperature (Band 08 of Advanced Himawari Imager)
tbb09      0 99 brightness temperature (Band 09 of Advanced Himawari Imager)
tbb10      0 99 brightness temperature (Band 10 of Advanced Himawari Imager)
tbb11      0 99 brightness temperature (Band 11 of Advanced Himawari Imager)
tbb12      0 99 brightness temperature (Band 12 of Advanced Himawari Imager)
tbb13      0 99 brightness temperature (Band 13 of Advanced Himawari Imager)
tbb14      0 99 brightness temperature (Band 14 of Advanced Himawari Imager)
tbb15      0 99 brightness temperature (Band 15 of Advanced Himawari Imager)
tbb16      0 99 brightness temperature (Band 16 of Advanced Himawari Imager)
endvars
