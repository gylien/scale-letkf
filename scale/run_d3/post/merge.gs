function merge (args)
'set gxout fwrite'
'set undef dfile'
'set fwrite -ap -be history.grd'


tmax=subwrd(args,1)
**zmax=60
zmax=3

t=1
while(t<=tmax)

'open U10.ctl'
'open V10.ctl'
'open T2.ctl'
'open Q2.ctl'
'open MSLP.ctl'
'open PREC.ctl'
'open SFC_TEMP.ctl'

'set t 't

'set z 1'
'd U10.1'
'd V10.2'
'd T2.3'
'd Q2.4'
'd MSLP.5'
'd PREC.6'
'd SFC_TEMP.7'

'close 7'
'close 6'
'close 5'
'close 4'
'close 3'
'close 2'
'close 1'

'open U.ctl'
'open V.ctl'
'open W.ctl'
'open T.ctl'
'open RH.ctl'
'open QV.ctl'
'open QC.ctl'
'open QI.ctl'
'open QR.ctl'
'open QS.ctl'
'open QG.ctl'

'set t 't
z=1
while(z<=zmax)
'set z 'z
'd U.1'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd V.2'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd W.3'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd T.4'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd RH.5'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd QV.6'
z=z+1
endwhile

z=1
while(z<=zmax)
'set z 'z
'd QC.7'
z=z+1
endwhile         

z=1
while(z<=zmax)
'set z 'z
'd QI.8'
z=z+1
endwhile 

z=1
while(z<=zmax)
'set z 'z
'd QR.9'
z=z+1
endwhile 

z=1
while(z<=zmax)
'set z 'z
'd QS.10'
z=z+1
endwhile 

z=1
while(z<=zmax)
'set z 'z
'd QG.11'
z=z+1
endwhile 

'close 11'
'close 10'
'close 9'
'close 8'
'close 7'
'close 6'
'close 5'
'close 4'
'close 3'
'close 2'
'close 1'


t=t+1
endwhile

'disable fwrite'
'quit'
