function plotDriver (args)

msg = 'Usage: plot_driver ctlfile tstart tend tskip'

ctlfile = subwrd(args, 1)
tstart = subwrd(args, 2)
tend = subwrd(args, 3)
tskip = subwrd(args, 4)

if (ctlfile = '') ; say msg ; return ; endif
if (tstart = '') ; tstart = 1 ; endif
if (tend = '') ; tend = tstart ; endif
if (tskip = '') ; tskip = 1 ; endif

'open 'ctlfile


it = tstart
while (it <= tend)
  'set t 'it

  'sfc_prcp_slp.gs'
  'sfc_wind.gs'
  'sfc_2mtemp.gs'
  '925_theq_wind_JMA.gs'
  '925_rh_temp_wind.gs'
  '850_theq_wind_JMA.gs'
  '850_rh_temp_wind.gs'
  '700_rh_temp_wind.gs'
  '500_rh_temp_wind.gs'

  'd3_sfc_uvt.gs'
  'd3_sfc_prcp.gs'

*  'sfc_temp.gs'
*  'sfc_accuprcp.gs'
*  'sfc_accusnow.gs'
*  '700_vvel_hgt.gs'
*  '500_vort_gph.gs'
*  '300_wspd_gph.gs'
*  'olr.gs'
*  'max_ref.gs'

  it = it + tskip
endwhile

'quit'

return
