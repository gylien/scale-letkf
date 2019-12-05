name = 'olr_ir'

convert = 'convert'
pngquant = '/apps/SLES12/opt/pngquant/2.10.1/bin/pngquant'

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
'set map 23 1 2'
'set grid off'
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

'set gxout shade2'
'set clevs -1e9 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280'
'set ccols 17 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161'
'd const(olr, -1e10, -u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout contour'
'set cint 10'
'set ccolor 1'
'set cstyle 3'
'set cthick 2'
'd lon'
'set cint 10'
'set ccolor 1'
'set cstyle 3'
'set cthick 2'
'd lat'

*'set gxout contour'
*'set clevs 900 908 916 924 932 940 948 956 964 972 980 988 996 1004 1012 1020 1028 1036 1044 1052 1060 1068 1076 1084 1092 1100'
*'set ccolor 24'
*'set cthick 3'
*'set clab off'
*'d slp'
*'set cint 8'
*'set ccolor 24'
*'set cthick 3'
*'set clab masked'
*'d slp'


'set string 1 c 5'
'set strsiz 0.11 0.15'
'draw string 5.5 7.58 `1Out lw rad (W/m  )  [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'
'set string 1 c 3'
'set strsiz 0.07 0.09'
'draw string 4.09 7.64 `1-2'

'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
'!rm -f 'outf'_tmp.png'

