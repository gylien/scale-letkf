function plot (args)

name = 'd3_sfc_prcp'

'settime.gs'
rcl = sublin(result, 1)
tt = subwrd(rcl, 1)
tfilename = subwrd(rcl, 2)
tstring = subwrd(rcl, 3)
trun = subwrd(rcl, 4)
'setrgb.gs'

outf = 'out_d3/'name'_'tfilename

'set mpdset hires'
'set display color white'
'set map 1 1 3'
'set grid on 3 1 1'
'set xlopts 1 4 0.12'
'set ylopts 1 4 0.12'
'set parea 0 11 1.15 7.35'
'c'
'set grads off'

'q dims'
'domain_d3.gs'
'q dims'
'set z 1'
'set gxout shaded'
'set clevs -1e9 2 4 8 12 18 25 35 50 70'
'set ccols 17 0 181 182 183 184 185 186 187 188 189'
if (tstring = trun)
'd const(PREC-1e8, -1e10, -u)'
else
'd const(PREC*3600, -1e10, -u)'
endif
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set gxout shaded'
'set clevs 0'
'set ccols 17 17'

*'d maskout(mask.2(t=1), 0.5-mask.2(t=1)) + lon.2'

'set string 1 c 5'
'set strsiz 0.09 0.12'
'draw string 5.5 7.5 `1Prev 1h precip (mm) [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
'!rm -f 'outf'_tmp.png'

return
