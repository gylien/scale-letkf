function merge (args)

*'open Uprs.ctl'
*'open Vprs.ctl'
*'open Wprs.ctl'
*'open Tprs.ctl'
*'open Gprs.ctl'
*'open QVprs.ctl'
*'open RHprs.ctl'
'open U10.ctl'
'open V10.ctl'
'open T2.ctl'
'open MSLP.ctl'
'open PREC.ctl'
'open SFC_TEMP.ctl'

'set gxout fwrite'
'set undef dfile'
'set fwrite -ap -be history.grd'


tmax=subwrd(args,1)
zmax=9

t=1
while(t<=tmax)
'set t 't

*z=1
*while(z<=zmax)
*'set z 'z
*'d Uprs.1'
*z=z+1
*endwhile

*z=1
*while(z<=zmax)
*'set z 'z
*'d Vprs.2'
*z=z+1
*endwhile

*z=1
*while(z<=zmax)
*'set z 'z
*'d Wprs.3'
*z=z+1
*endwhile
*
*z=1
*while(z<=zmax)
*'set z 'z
*'d Tprs.4'
*z=z+1
*endwhile
*
*z=1
*while(z<=zmax)
*'set z 'z
*'d Gprs.5'
*z=z+1
*endwhile
*
*z=1
*while(z<=zmax)
*'set z 'z
*'d QVprs.6'
*z=z+1
*endwhile
*
*z=1
*while(z<=zmax)
*'set z 'z
*'d RHprs.7'
*z=z+1
*endwhile         

'set z 1'
'd U10.1'
'd V10.2'
'd T2.3'
'd MSLP.4'
'd PREC.5'
'd SFC_TEMP.6'

t=t+1
endwhile

'disable fwrite'
'quit'
