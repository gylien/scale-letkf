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


  'sfc_prcp.gs'
  'sfc_uvt.gs'
  '1000m_dbz.gs'
  '1000m_uvw.gs'
  '3000m_dbz.gs'
  '3000m_uvw.gs'
  '5000m_dbz.gs'
  '5000m_uvw.gs'


  it = it + tskip
endwhile

'quit'

return
