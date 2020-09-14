name = '3000m_dbz'

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

*** tori aezu
plev=700

*'d qc'
*'d qi'

*** reflectivity ***

far=2.53e4
fbr=1.84
fas=3.48e4
fbs=1.66
fag=5.54e3
fbg=1.70

'rho= 'plev' * 100 / 287.1 / t '

'zr = 1.0e3 * rho * qr'
'zr = 'far'* pow(zr,'fbr')'
'zs = 1.0e3 * rho * qs'
'zs = 'fas'* pow(zs,'fbs')'
'zg = 1.0e3 * rho * qg'
'zg = 'fag'* pow(zg,'fbg')'
'ref=10*log10(zr+zs+zg+1.01)'

'color_radar_share.gs'
'd const(ref,0,-u)'
'cbarn.gs 0.8 0 5.5 0.7 0.8 1'


******
'srh=smth9(rh)'
'srh=smth9(srh)'
'srh=smth9(srh)'

'set gxout contour'
'set clevs 70 90'
'set ccolor 4'
'set cthick 2'
'set cstyle 1'
'set clab on'
'd srh'


'set string 1 c 5'
'set strsiz 0.09 0.12'
'draw string 5.5 7.5 `13000m radar ref (dbz) [ Run: `0'trun'`1 | VT: `0'tstring'`1 ]'

'!rm -f 'outf'_tmp.png'
'gxprint 'outf'_tmp.png'
'!convert -trim +repage 'outf'_tmp.png 'outf'.png'
'!rm -f 'outf'_tmp.png'
