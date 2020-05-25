function plot (args)

name = '850_theq'

tstring = subwrd(args, 1)
trun = subwrd(args, 2)
tfilename = subwrd(args, 3)

'setrgb.gs'

outf = 'out/'name'_'tfilename

'set mpdset hires'
'set display color white'
'set map 1 1 3'
'set grid on 3 1'
'set xlopts 1 4 0.12'
'set ylopts 1 4 0.12'
'set parea 0 11 1.15 7.35'
'c'
'set grads off'

'domain.gs'

'set lev 850'

'set gxout shaded'
'set clevs -1e9 70 80 90'
'set ccols 17 0 211 212 213'
'd const(rh, -1e10, -u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout barb'
'set ccolor 4'
'set cthick 3'
'set digsiz 0.042'
'd skip(u*1.94384,8,8);skip(v*1.94384,4,4)'

'tlk = 55 + 1 / (1/(temp-55) - (log(const(rh,1e-8,-u)/100)) / 2840 ) '
*** Teten formula
'es = 6.11 * pow(10,(7.5 * (temp-273.15)/(temp-38.2))) '
'qv = 0.01 * rh * (0.622 * es / (lev-es) )'

*** (1000/850)^(rd/cp) := 1.0476

'theq = smth9(temp * 1.0476 * exp(2.675 * 1000 * qv / tlk ))'

'set gxout contour'
'set clevs 336 342 348'
'set ccolor 201'
'set cthick 3'
'set cstyle 1'
'set clab on'
'd theq'

'set clevs 330'
'set ccolor 201'
'set cthick 6'
'set cstyle 1'
'set clab on'
'd theq'

'set clevs 321 327'
'set ccolor 203'
'set cthick 1'
'set cstyle 1'
'set clab off'
*'d theq'

'set clevs 318 324'
'set ccolor 203'
'set cthick 3'
'set cstyle 1'
'set clab on'
'd theq'

'set clevs 315'
'set ccolor 203'
'set cthick 6'
'set cstyle 1'
'set clab off'
'd theq'

'set clevs 291 297 303 309'
'set ccolor 202'
'set cthick 1'
'set cstyle 1'
'set clab off'
*'d theq'

'set clevs 294 306 312'
'set ccolor 202'
'set cthick 3'
'set cstyle 1'
'set clab on'
'd theq'

'set clevs 300'
'set ccolor 202'
'set cthick 6'
'set cstyle 1'
'set clab on'
'd theq'

'set string 1 c 5'
'set strsiz 0.09 0.12'
'draw string 5.5 7.5 `1850 hPa RH (%) / Temp (C) / Wind (kts)  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'


'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
'!rm -f 'outf'_tmp.png'

return
