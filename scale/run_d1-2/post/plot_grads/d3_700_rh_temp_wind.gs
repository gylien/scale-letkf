name = 'd3_700_temp'

'settime.gs'
rcl = sublin(result, 1)
tt = subwrd(rcl, 1)
tfilename = subwrd(rcl, 2)
tstring = subwrd(rcl, 3)
trun = subwrd(rcl, 4)
'setrgb.gs'

outf = 'out_d3/'name'_'tfilename


'set mpdset hires'
'set display color white'
'set map 1 1 3'
'set grid on 3 1 1'
'set xlopts 1 4 0.12'
'set ylopts 1 4 0.12'
*'set parea 0 11 0.8 7.8'
'c'
'set grads off'

'domain_d3.gs'

*'set lon 85 155'
*'set lat 0 45'
'set lev 700'
*'set xlint 10'
*'set ylint 10'

'set gxout shade2'
'set clevs -1e9 70 80 90'
'set ccols 17 0 211 212 213'
'd const(rhprs, -1e10, -u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout barb'
'set ccolor 4'
'set cthick 3'
'set digsiz 0.042'
'd skip(uprs*1.94384,3,3);vprs*1.94384'

'set gxout contour'
'set clevs 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30'
'set ccolor 201'
'set cthick 4'
'set cstyle 3'
'set clab on'
'd tprs-273.15'
'set clevs 0'
'set ccolor 202'
'set cthick 6'
'set cstyle 3'
'set clab on'
'd tprs-273.15'
'set clevs -18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1'
'set ccolor 203'
'set cthick 4'
'set cstyle 3'
'set clab on'
'd tprs-273.15'

'set string 1 c 5'
'set strsiz 0.09 0.12'
'draw string 5.5 7.58 `1700 hPa RH (%) / Temp (C) / Wind (kts)  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
'!rm -f 'outf'_tmp.png'
