function convert(args)

ctlfile = subwrd(args,1)

say ctlfile

levs = '1000 925 850 700 500 300 200 100 50 20 10'
nz = 11

var3d = 'ugrdprs vgrdprs TMPprs HGTprs RHprs'
nv3d = 5
*var2d = 'UGRD10m VGRD10m TMP2m PRMSLmsl'
var2d = 'UGRD10m VGRD10m TMP2m MSLETmsl'
nv2d = 4

'reinit'
'open 'ctlfile
'set gxout fwrite'
'set lon 95.5 174.5'
'set lat 12.5 54.5'

v = 1
while (v <= nv3d)
  var = subwrd(var3d, v)

  z = 1
  while (z <= nz)
    llev = subwrd(levs, z)
    'set lev 'llev

    say 'd 'var', z = 'z
    'd 'var

    z = z + 1
  endwhile

  v = v + 1
endwhile

'set lev 1'
v = 1
while (v <= nv2d)
  var = subwrd(var2d, v)
  
  say 'd 'var', z = 1'
  'd 'var
  
  v = v + 1
endwhile

'disable fwrite'
'quit'

return
