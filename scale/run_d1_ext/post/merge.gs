'open Uprs.ctl'
'open Vprs.ctl'
'open Wprs.ctl'
'open Tprs.ctl'
'open Gprs.ctl'
'open QVprs.ctl'
'open RHprs.ctl'
'open U10.ctl'
'open V10.ctl'
'open T2.ctl'
'open MSLP.ctl'
'open PREC.ctl'
'open SFC_TEMP.ctl'



'set gxout fwrite'
'set undef dfile'
'set fwrite -ap -be history.grd'


tmax=21
zmax=9

t=1
while(t<=tmax)
'set t 't

z=1
while(z<=zmax)
'set z 'z
'd Uprs.1'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd Vprs.2'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd Wprs.3'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd Tprs.4'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd Gprs.5'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd QVprs.6'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd RHprs.7'
z=z+1
endwhile

'set z 1'
'd U10.8'
'd V10.9'
'd T2.10'
'd MSLP.11'
'd PREC.12'
'd SFC_TEMP.13'

t=t+1
endwhile

'disable fwrite'
'quit'
