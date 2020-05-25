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
*tth = (tfcst-1)*21600
tth = (tfcst-1)*10800
rc = math_format('%06.0f', tth)
ttf = rc

*tnow=(tfcst-1) * 2 +1
tnow=tfcst
'set t 'tnow

'925_theq_wind.gs 'tstring' 'tfcstbase' f'ttf
'925_rh_temp_wind.gs 'tstring' 'tfcstbase' f'ttf
'850_theq_wind.gs 'tstring' 'tfcstbase' f'ttf
'850_rh_temp_wind.gs 'tstring' 'tfcstbase' f'ttf
'700_rh_temp_wind.gs 'tstring' 'tfcstbase' f'ttf
'500_rh_temp_wind.gs 'tstring' 'tfcstbase' f'ttf

'quit'

return
