function plot (args)

name = 'sfc_prcp'

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
'set z 1'
'set xlint 10'
'set ylint 10'

'set gxout shaded'
'set clevs -1e9 2 4 8 12 18 25 35 50 70'
'set ccols 17 0 181 182 183 184 185 186 187 188 189'
if (tstring = trun)
  'd const(TMPsfc-1e8, -1e10, -u)'
else
*  'd const(PRATEavesfc*21600, -1e10, -u)'
  'd const(APCPsfc*21600, -1e10, -u)'
endif
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'


'set gxout contour'
'set cint 4'
'set ccolor 21'
'set cthick 4'
'set clab off'
'd PRMSLmsl*0.01'
*'d MSLETmsl*0.01'
'set cint 8'
'set ccolor 21'
'set cthick 4'
'set clab on'
'd PRMSLmsl*0.01'
*'d MSLETmsl*0.01'


'set gxout shaded'
'set clevs 0'
'set ccols 17 17'
'd maskout(mask.2(t=1), 0.5-mask.2(t=1)) + lon.2'


'set string 1 c 5'
'set strsiz 0.11 0.15'
'draw string 5.5 7.58 `1Prev 6h precip (mm) / SLP (hPa)  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

*'!rm -f 'outf'_tmp.png'
*'gxprint 'outf'_tmp.png'
*'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
*'!rm -f 'outf'_tmp.png'

'gxprint 'outf'.eps'

*'print'
*'disable print'

*'!gxeps -cR -i 'outf'.gmf'
'!'convert' -density 400x382 -trim +antialias +matte 'outf'.eps 'outf'_tmp.png'
'!'convert' -resize 25% 'outf'_tmp.png 'outf'_tmp2.png'
**'!'convert' -resize 25% 'outf'_tmp.png 'outf'.gif'
'!'pngquant' --force --quality 60-75 'outf'_tmp2.png -o 'outf'.png'
'!rm -f 'outf'.gmf'
'!rm -f 'outf'_tmp.png'
'!rm -f 'outf'_tmp2.png'
**'!rm -f 'outf'_tmp.gif'

return
