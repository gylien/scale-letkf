function plot (args)

name = '300_wspd'

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
'set map 22 1 3'
'set grid on 3 1'
'set xlopts 1 4 0.12'
'set ylopts 1 4 0.12'
*'set parea 0 11 0.8 7.8'
'c'
'set grads off'

'set lon 95.6726 174.328'
'set lat 12.8555 54.0347'
'set lev 300'
'set xlint 10'
'set ylint 10'

'set gxout shaded'
'set clevs -1e9 50 60 70 80 90 100 110 120 130 150 160 170 180 190 200 210 220 230'
'set ccols 17 0 41 42 43 44 45 46 47 48 61 62 63 64 65 66 67 68 69 70'
'd const(mag(UGRDprs,VGRDprs)*1.94384, -1e10, -u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout contour'
'set clevs 50 100 150 200'
'set ccolor 16'
'set cthick 2'
'set clab off'
'd mag(UGRDprs,VGRDprs)*1.94384'

'set gxout contour'
'set cint 120'
'set ccolor 1'
'set cthick 4'
'set clab off'
'd HGTprs'
'set cint 240'
'set ccolor 1'
'set cthick 4'
'set clab on'
'd HGTprs'


'set gxout shaded'
'set clevs 0'
'set ccols 17 17'
'd maskout(mask.2(t=1), 0.5-mask.2(t=1)) + lon.2'


'set string 1 c 5'
'set strsiz 0.11 0.15'
'draw string 5.5 7.58 `1300 hPa Wspd (kts) / Gpt hgt (m)  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

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

return
