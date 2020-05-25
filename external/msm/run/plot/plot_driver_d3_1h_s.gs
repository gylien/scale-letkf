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

tstring = tfcstplot

tfilename = 'f'
tth = (tfcst-1)*3600
rc = math_format('%06.0f', tth)
ttf = rc

tnow=tfcst

'set t 'tnow

'd3_sfc_prcp.gs 'tstring' 'tfcstbase' f'ttf
'd3_sfc_uvt.gs 'tstring' 'tfcstbase' f'ttf

'quit'

return
