name = '3000m_uvw'

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
'set map 1 1 3'
'set grid on 3 1 1'
'set xlopts 1 4 0.12'
'set ylopts 1 4 0.12'
'set parea 0 11 1.15 7.35'
'c'
'set grads off'

'domain.gs'
'set lev 3000'

'color_w'
'd w'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'


*'set gxout shade2'
*'set clevs -1e9 2 4 6 8 10 12 14 16 18 20 22 24 27 30'
*'set ccols 17 0 220 221 222 223 224 225 226 227 228 229 230 231 232 233'
*'d const(mag(U10,V10), -1e10, -u)'
*'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout stream'
'set cthick 1'
'set ccolor 18'
'set strmden -2'
'd u;v'

'set string 1 c 5'
'set strsiz 0.09 0.12'
'draw string 5.5 7.5 `13000m wind (m/s)  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'


'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
'!rm -f 'outf'_tmp.png'
