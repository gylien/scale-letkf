function ref5000m (args)

ctlfile = subwrd(args, 1)
tstart = subwrd(args, 2)
tend = subwrd(args, 3)
tskip = subwrd(args, 4)

name = 'dbz5000m_mdet'
msg = 'Usage: plot_driver ctlfile tstart tend tskip'
convert = 'convert'
***pngquant = '/apps/SLES12/opt/pngquant/2.10.1/bin/pngquant'

if (ctlfile = '') ; say msg ; return ; endif
if (tstart = '') ; tstart = 1 ; endif
if (tend = '') ; tend = tstart ; endif
if (tskip = '') ; tskip = 1 ; endif

'open 'ctlfile

it = tstart
while (it <= tend)
'set t 'it 
'settime.gs'
rcl = sublin(result, 1)
tt = subwrd(rcl, 1)
tfilename = subwrd(rcl, 2)
tstring = subwrd(rcl, 3)
trun = subwrd(rcl, 4)

'setrgb.gs'

outf = 'out/'name'_'tfilename

*'set mpdset worldmap'
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

'set lev 27' 

'set gxout shaded'
'set clevs -1e9 0 5 10 15 20 25 30 35 40 45 50 55 60'
'set ccols 17 0 240 241 242 243 244 245 246 247 248 249 250 251 252'

'd const(dbz, -1e10, -u)'

'cbarn.gs 0.8 0 5.5 0.7 0.8 1'

'set string 1 c 5'
'set strsiz 0.09 0.12'
'draw string 5.5 7.5 `1Max reflec (dBZ) [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'_0001.png'
'!rm -f 'outf'_tmp.png'


***'print'
***'disable print'

*'!gxeps -cR -i 'outf'.gmf'
*'!'convert' -density 400x382 -trim +antialias +matte 'outf'.eps 'outf'_tmp.png'
*'!'convert' -resize 25% 'outf'_tmp.png 'outf'_tmp2.png'
*'!'pngquant' --force --quality 60-75 'outf'_tmp2.png -o 'outf'.png'
*'!rm -f 'outf'.gmf'
*'!rm -f 'outf'_tmp.png'
*'!rm -f 'outf'_tmp2.png'

  it = it + tskip
endwhile

'quit'
