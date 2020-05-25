function plotDriver (args)

msg = 'Usage: plot_driver ctlfile tfcstbase tfcstplot tfcst'

ctlfile = subwrd(args, 1)
tfcstbase = subwrd(args, 2)
tfcstplot = subwrd(args, 3)
tfcst = subwrd(args, 4)

if (ctlfile = '') ; say msg ; return ; endif
if (tfcstbase = '') ; say msg ; return ; endif
if (tfcstplot = '') ; say msg ; return ; endif
if (tfcst = '') ; say msg ; return ; endif

'open 'ctlfile
'open mask.ctl'

'q dims'
*tline = sublin(result, 5)
*tstring = subwrd(tline, 6)

tstring = tfcstplot

tfilename = 'f'
tth = (tfcst-1)*3600
rc = math_format('%06.0f', tth)
ttf = rc

*tnow=(tfcst-1) * 6 +1

tnow=tfcst

'set t 'tnow


*tsta=tnow-5
*if (tfcst=1)
*'r6h = 0.0'
*else
*'define r6h = sum(r1h,t='tsta',t='tnow')'
*endif



'sfc_prcp_slp.gs 'tstring' 'tfcstbase' f'ttf
'sfc_wind.gs 'tstring' 'tfcstbase' f'ttf
'sfc_2mtemp.gs 'tstring' 'tfcstbase' f'ttf
*'sfc_temp.gs 'tstring' 'tfcstbase' f'ttf


*'olr.gs 'tstring' 'tfcstbase' f'ttf
*'max_ref.gs 'tstring' 'tfcstbase' f'ttf

'quit'

return
