function plotDriver (args)

msg = 'Usage: plot_driver ctlfile tfcstbase tfcst'

ctlfile = subwrd(args, 1)
tfcstbase = subwrd(args, 2)
tfcst = subwrd(args, 3)

if (ctlfile = '') ; say msg ; return ; endif
if (tfcstbase = '') ; say msg ; return ; endif
if (tfcst = '') ; say msg ; return ; endif

'open 'ctlfile
'open mask.ctl'

'q dims'
tline = sublin(result, 5)
tstring = subwrd(tline, 6)
tfilename = 'f'
tth = (tfcst-1)*21600
rc = math_format('%06.0f', tth)
ttf = rc

'sfc_prcp_slp.gs 'tstring' 'tfcstbase' f'ttf
'sfc_wind.gs 'tstring' 'tfcstbase' f'ttf
'sfc_2mtemp.gs 'tstring' 'tfcstbase' f'ttf
'sfc_temp.gs 'tstring' 'tfcstbase' f'ttf
'850_theq_wind.gs 'tstring' 'tfcstbase' f'ttf
'850_rh_temp_wind.gs 'tstring' 'tfcstbase' f'ttf
'700_vvel_hgt.gs 'tstring' 'tfcstbase' f'ttf
'500_vort_gph.gs 'tstring' 'tfcstbase' f'ttf
'300_wspd_gph.gs 'tstring' 'tfcstbase' f'ttf

*'olr.gs 'tstring' 'tfcstbase' f'ttf
*'max_ref.gs 'tstring' 'tfcstbase' f'ttf

'quit'

return
