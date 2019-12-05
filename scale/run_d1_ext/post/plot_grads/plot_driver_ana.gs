function plotDriver (args)

msg = 'Usage: plot_driver ctlfile tstart tend tskip'

ctlfile = subwrd(args, 1)
timeana = subwrd(args, 2)

if (ctlfile = '') ; say msg ; return ; endif
if (timeana = '') ; tstart = 1 ; return ; endif

'open 'ctlfile

  'set time 'timeana

  'sfc_prcp_slp.gs'
*  'sfc_wind.gs'
*  'sfc_2mtemp.gs'
*  'sfc_temp.gs'
*  'sfc_accuprcp.gs'
*  'sfc_accusnow.gs'
*  '850_rh_temp_wind.gs'
*  '850_theq_wind.gs'
*  '700_vvel_hgt.gs'
*  '500_vort_gph.gs'
*  '300_wspd_gph.gs'
*  'olr.gs'
*  'max_ref.gs'

'quit'

return
