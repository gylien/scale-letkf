name = 'max_ref'

convert = 'convert'
pngquant = '/work/hp150019/c24140/scale-letkf-rt/external/lib/pngquant/pngquant'

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
'set grid on 3 1'
'set xlopts 1 4 0.12'
'set ylopts 1 4 0.12'
*'set parea 0 11 0.8 7.8'
'c'
'set grads off'

*'set lon 85 155'
*'set lat 0 45'
'set z 1'
'set xlint 10'
'set ylint 10'

'set gxout shaded'
'set clevs -1e9 0 5 10 15 20 25 30 35 40 45 50 55 60'
'set ccols 17 0 240 241 242 243 244 245 246 247 248 249 250 251 252'
'd const(max_dbz, -1e10, -u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'


'set string 1 c 5'
'set strsiz 0.11 0.15'
'draw string 5.5 7.58 `1Max reflec (dBZ) [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

'gxprint 'outf'.eps'


*'print'
*'disable print'

*'!gxeps -cR -i 'outf'.gmf'

'!'convert' -density 400x382 -trim +antialias +matte 'outf'.eps 'outf'_tmp.png'
'!'convert' -resize 25% 'outf'_tmp.png 'outf'_tmp2.png'
*'!'convert' -resize 25% 'outf'_tmp.png 'outf'.gif'
'!'pngquant' --force --quality 60-75 'outf'_tmp2.png -o 'outf'.png'
'!rm -f 'outf'.gmf'
'!rm -f 'outf'_tmp.png'
'!rm -f 'outf'_tmp2.png'
*'!rm -f 'outf'_tmp.gif'

