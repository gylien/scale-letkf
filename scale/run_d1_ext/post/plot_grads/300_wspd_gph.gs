name = '300_wspd'

convert = 'convert'
*pngquant = '/apps/SLES12/opt/pngquant/2.10.1/bin/pngquant'

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
'set lev 300'
'set xlint 10'
'set ylint 10'

'set gxout shade2'
'set clevs -1e9 50 60 70 80 90 100 110 120 130 150 160 170 180 190 200 210 220 230'
'set ccols 17 0 41 42 43 44 45 46 47 48 61 62 63 64 65 66 67 68 69 70'
'd const(mag(Uprs,Vprs)*1.94384, -1e10, -u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout contour'
'set clevs 50 100 150 200'
'set ccolor 16'
'set cthick 2'
'set clab off'
'd mag(Uprs,Vprs)*1.94384'

'set gxout contour'
'set cint 120'
'set ccolor 1'
'set cthick 4'
'set clab off'
'd Gprs'
'set cint 240'
'set ccolor 1'
'set cthick 4'
'set clab on'
'd Gprs'


'set string 1 c 5'
'set strsiz 0.11 0.15'
'draw string 5.5 7.58 `1300 hPa Wspd (kts) / Gpt hgt (m)  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
'!rm -f 'outf'_tmp.png'
