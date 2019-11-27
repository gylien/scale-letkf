name = '700_vvel'

convert = 'convert'
pngquant = '/apps/SLES12/opt/pngquant/2.10.1/bin/pngquant'

'settime.gs'
rcl = sublin(result, 1)
tt = subwrd(rcl, 1)
tfilename = subwrd(rcl, 2)
tstring = subwrd(rcl, 3)
trun = subwrd(rcl, 4)
'setrgb.gs'

outf = 'out/'name'_'tfilename


'set mpdset hires'
'set display color white'
'set map 22 1 3'
'set grid on 3 1 1'
'set xlopts 1 4 0.12'
'set ylopts 1 4 0.12'
*'set parea 0 11 0.8 7.8'
'c'
'set grads off'

*'set lon 85 155'
*'set lat 0 45'
'set lev 700'
'set xlint 10'
'set ylint 10'

'set gxout shade2'
'set clevs -1e9 -30 -20 -15 -10 -6 -3 -1.5 1.5 3 6 10 15 20 30 50 70'
'set ccols 17 107 106 105 104 103 102 101 0 121 122 123 124 125 126 127 128 129'
'd const(w*100, -1e10, -u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout contour'
'set clevs -6 6'
'set ccolor 16'
'set cthick 2'
'set clab off'
'd w*100'

'set gxout contour'
'set cint 30'
'set ccolor 1'
'set cthick 4'
'set clab off'
'd z'
'set cint 60'
'set ccolor 1'
'set cthick 4'
'set clab on'
'd z'


'set string 1 c 5'
'set strsiz 0.11 0.15'
'draw string 5.5 7.58 `1700 hPa Vert vel (cm s  ) / Gpt hgt (m)  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

'set string 1 c 3'
'set strsiz 0.07 0.09'
'draw string 3.7 7.64 `1-1'

'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
'!rm -f 'outf'_tmp.png'

