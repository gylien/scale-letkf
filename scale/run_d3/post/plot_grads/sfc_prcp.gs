name = 'sfc_prcp'

convert = 'convert'

'settime.gs'
rcl = sublin(result, 1)
tt = subwrd(rcl, 1)
tfilename = subwrd(rcl, 2)
tstring = subwrd(rcl, 3)
trun = subwrd(rcl, 4)
'setrgb.gs'

outf = 'out/'name'_'tfilename
*'enable print 'outf'.gmf'


'set mpdset hires'
'set display color white'
'set map 1 1 3'
'set grid on 3 1 1'
'set xlopts 1 4 0.12'
'set ylopts 1 4 0.12'
'set parea 0 11 1.15 7.35'
'c'
'set grads off'


'domain.gs'
'set z 1'

'set gxout shaded'
'set clevs -1e9 2 4 8 12 18 25 35 50 70'
'set ccols 17 0 181 182 183 184 185 186 187 188 189'
if (tt = 1)
  'd const(PREC-1e8, -1e10, -u)'
else
  'd const(PREC*3600, -1e10, -u)'
endif
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout contour'
'set cint 2'
'set cstyle 3'
'set ccolor 21'
'set cthick 4'
'set clab off'
'd MSLP*0.01'
'set gxout contour'
'set cint 4'
'set cstyle 1'
'set ccolor 21'
'set cthick 4'
'set clab off'
'd MSLP*0.01'
'set cint 8'
'set ccolor 21'
'set cstyle 1'
'set cthick 4'
'set clab on'
'd MSLP*0.01'


'set string 1 c 5'
'set strsiz 0.09 0.12'
'draw string 5.5 7.5 `1Prev 1h precip (mm) / SLP (hPa)  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
'!rm -f 'outf'_tmp.png'


*'print'
*'disable print'

*'!gxeps -cR -i 'outf'.gmf'
*'!'convert' -density 400x382 -trim +antialias +matte 'outf'.eps 'outf'_tmp.png'
*'!'convert' -resize 25% 'outf'_tmp.png 'outf'_tmp2.png'
*'!'pngquant' --force --quality 60-75 'outf'_tmp2.png -o 'outf'.png'
*'!rm -f 'outf'.gmf'
*'!rm -f 'outf'_tmp.png'
*'!rm -f 'outf'_tmp2.png'

