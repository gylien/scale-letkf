function plot (args)

name = '500_vort'

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
'set lev 500'
'set xlint 10'
'set ylint 10'

'set gxout shaded'
'set clevs -1e9 -14 -12 -10 -8 -6 -4 -2 2 4 6 8 10 12 14 16 18 20'
'set ccols 17 107 106 105 104 103 102 101 0 121 122 123 124 125 126 127 128 129 130'
'd const(hcurl(UGRDprs,VGRDprs)*1.e5, -1e10, -u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout contour'
'set clevs -4 4'
'set ccolor 16'
'set cthick 2'
'set clab off'
'd hcurl(UGRDprs,VGRDprs)*1.e5'

'set gxout contour'
'set cint 60'
'set ccolor 1'
'set cthick 4'
'set clab off'
'd HGTprs'
'set cint 120'
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
'draw string 5.5 7.58 `1500 hPa Vort (10  s  ) / Gpt hgt (m)  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'
'set string 1 c 3'
'set strsiz 0.07 0.09'
'draw string 3.38 7.64 `1-5  -1'

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
