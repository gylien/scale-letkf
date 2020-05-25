function plot (args)

name = 'sfc_wind'

convert = 'convert'
pngquant = '/apps/SLES12/opt/pngquant/2.10.1/bin/pngquant'

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
'set parea 0 11 1.15 7.35'
'c'
'set grads off'

'domain.gs'
'set z 1'


'uabs=mag(u,v)'

'set gxout shaded'
'set clevs -1e9 2 4 6 8 10 12 14 16 18 20 22 24 27 30'
'set ccols 17 0 220 221 222 223 224 225 226 227 228 229 230 231 232 233'
'd const(uabs, -1e10, -u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout stream'
'set cthick 1'
'set ccolor 18'
'set strmden 6'
'd u;v'

'set gxout contour'
'set clab on'
'set ccolor 25'
'set cthick 3'
'set clevs 4'
'd uabs*0.727'
'set clab on'
'set ccolor 25'
'set cthick 3'
'set clevs 7'
'd uabs*0.504'
'set clab on'
'set ccolor 25'
'set cthick 3'
'set clevs 9'
'd uabs*0.433'
'set clab on'
'set ccolor 25'
'set cthick 3'
'set clevs 12'
'd uabs*0.367'


'set gxout shaded'
'set clevs 0'
'set ccols 17 17'
*'d maskout(mask.2(t=1), 0.5-mask.2(t=1)) + lon.2'


'set string 1 c 5'
'set strsiz 0.09 0.12'
'draw string 5.5 7.5 `1Sfc (10m) wind (m/s;                        )  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

'set string 25 c 5'
'set strsiz 0.09 0.12'
'draw string 4.6 7.5 `1contour: Beaufort scale'


'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
'!rm -f 'outf'_tmp.png'

*'print'
*'disable print'

*'!gxeps -cR -i 'outf'.gmf'
*'!'convert' -density 400x382 -trim +antialias +matte 'outf'.eps 'outf'_tmp.png'
*'!'convert' -resize 25% 'outf'_tmp.png 'outf'_tmp2.png'
**'!'convert' -resize 25% 'outf'_tmp.png 'outf'.gif'
*'!'pngquant' --force --quality 60-75 'outf'_tmp2.png -o 'outf'.png'
*'!rm -f 'outf'.gmf'
*'!rm -f 'outf'_tmp.png'
*'!rm -f 'outf'_tmp2.png'
**'!rm -f 'outf'_tmp.gif'


return
