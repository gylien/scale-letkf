function plot (args)

name = '850_theq'

convert = 'convert'
pngquant = '/work/hp150019/c24140/scale-letkf-rt/external/lib/pngquant/pngquant'

tstring = subwrd(args, 1)
trun = subwrd(args, 2)
tfilename = subwrd(args, 3)

'setrgb.gs'

outf = 'out/'name'_'tfilename
*'enable print 'outf'.gmf'


'set mpdset hires'
'set display color white'
'set map 1 1 3'
'set grid on 3 1'
'set xlopts 1 4 0.12'
'set ylopts 1 4 0.12'
*'set parea 0 11 0.8 7.8'
'c'
'set grads off'

'set lon 95.6726 174.328'
'set lat 12.8555 54.0347'
'set lev 850'
'set xlint 10'
'set ylint 10'

'set lev 700'

'set gxout shade2'
'set clevs -1e9 70 80 90'
'set ccols 17 0 211 212 213'
'd const(RHprs, -1e10, -u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set lev 850'

'set gxout barb'
'set ccolor 1'
'set cthick 3'
'set digsiz 0.042'
'd skip(UGRDprs*1.94384,5,5);VGRDprs*1.94384'

'tlk = 55 + 1 / (1/(TMPprs-55) - (log(const(RHprs,1e-8,-u)/100)) / 2840 ) '

*** (1000/850)^(rd/cp) := 1.0476

'qv = 0.622 * RHprs * 0.01 * (6.11 * pow(10, 7.5*(TMPprs-273.15)/(TMPprs-35.8) ) )  / 850 '

'theq = smth9(TMPprs * 1.0476 * exp(2.675 * 1000 * qv / tlk ))'


'set gxout contour'


'set clevs 342 348'
'set ccolor 201'
'set cthick 3'
'set cstyle 1'
'set clab off'
'd theq'

'set clevs 336'
'set ccolor 201'
'set cthick 6'
'set cstyle 1'
'set clab on'
'd theq'

'set clevs 324 330'
'set ccolor 203'
'set cthick 3'
'set cstyle 1'
'set clab off'
'd theq'

'set clevs 318'
'set ccolor 203'
'set cthick 6'
'set cstyle 1'
'set clab on'
'd theq'

'set clevs 306 312'
'set ccolor 202'
'set cthick 3'
'set cstyle 1'
'set clab off'
'd theq'

'set clevs 300'
'set ccolor 202'
'set cthick 6'
'set cstyle 1'
'set clab on'
'd theq'


'set string 1 c 5'
'set strsiz 0.11 0.15'
'draw string 5.5 7.58 `1700 hPa RH (%) / 850hPa eq.pot.temp (K) / Wind (kts)  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

*'print'
*'disable print'

'gxprint 'outf'.eps'



*'!gxeps -cR -i 'outf'.gmf'
'!'convert' -density 400x382 -trim +antialias +matte 'outf'.eps 'outf'_tmp.png'
'!'convert' -resize 25% 'outf'_tmp.png 'outf'_tmp2.png'
*'!'convert' -resize 25% 'outf'_tmp.png 'outf'.gif'
'!'pngquant' --force --quality 60-75 'outf'_tmp2.png -o 'outf'.png'
'!rm -f 'outf'.gmf'
'!rm -f 'outf'_tmp.png'
'!rm -f 'outf'_tmp2.png'
*'!rm -f 'outf'_tmp.gif'

return
