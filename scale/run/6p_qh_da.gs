
'reinit'

'set display color white'
'c'

'set grads off'
'set tlsupp off'

ctl_nat = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_NATURE/20010101000000/fcst/mean/QHYD_d01z-3d.ctl"
ctl_noda = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_NODA_AKSOY/ctl/mean.ctl"
ctl_da = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_DA_PAWR_AKSOY/ctl/mean.ctl"

radar_lon = 180
radar_lat = 180


var.3 = "qhyd*1.e3"
var.2 = "(qc+qr+qi+qs+qg)*1.e3"
var.1 = "(qc+qr+qi+qs+qg)*1.e3"
tit = "Mixing ratio of hydrometeors (g/kg) & vertical vorticity (`3,`00.5 10`a-2`ns`a-1`n)"

ctit.1 = "NODA"
ctit.2 = "PAWR DA"
ctit.3 = "Nature run"


'defcolor'


level = 4
*level = 2
level = 6
level = 10

gt1.1 = 14
gt1.2 = 15

gt2.1 = 5
gt2.2 = 10

clevs = "0.1 0.5 1 2 3 4 5 6 7"
ccols = "0  242 241 240 230 231 232 233 234 235"

slon = 115
elon = 245
slat = 115
elat = 245

xmin = 0.8
dx = 3
dy = 3
ymax = 7.9
pdx = 0.5
pdy = 0.5

xlint = 20
ylint = 20


imax = 3
i = 1
while(i <= imax)
 xx = (dx + pdx) * (i - 1)

 
 if(i = 3); 'open 'ctl_nat;  endif
 if(i = 2); 'open 'ctl_da; endif
 if(i = 1); 'open 'ctl_noda; endif

 'set mproj off'


 jmax = 1
 j = 1
 while(j <= jmax)
  yy = (dy + pdy) * (j - 1)

  'set parea 'xmin+xx' 'xmin+xx+dx' 'ymax-yy-dy' 'ymax-yy

  'set lon 'slon' 'elon
  'set lat 'slat' 'elat

  'set lev 'level
  clevel = subwrd(result, 4)
  if(i=3); 'set t 'gt1.j; endif
  if(i<=2); 'set t 'gt2.j; endif
  'q time'
  ctime = subwrd(result, 3)
  ft = (gt1.j - 1 ) * 5

  'set xlint 'xlint
  'set ylint 'ylint

  'set clevs 'clevs
  'set ccols 'ccols

  'set xlab on'
  'set ylab on'
  'set grid on'
  'set frame on'

  'set gxout shade2'

  'd 'var.i

  'set xlab off'
  'set ylab off'
  'set grid off'
  'set frame off'


  ptit = "FT="ft"min, z="clevel"km"

  'set string 1 c 3'
  'set strsiz 0.1 0.12'
  'draw string 'xmin+xx+dx/2' 'ymax-yy+0.1' 'ptit

  'q ll2xy 'radar_lon' 'radar_lat
  radar_x = subwrd(result,1)
  radar_y = subwrd(result,2)

  'draw mark 3 'radar_x' 'radar_y' 0.1'

  'q ll2xy 'radar_lon' 'radar_lat+60
  radar_y60 = subwrd(result,2)
  rad = radar_y60 - radar_y
  'draw mark 2 'radar_x' 'radar_y' 'rad*2


  j = j + 1
 endwhile

 'close 1'

  'set string 1 c 4'
  'set strsiz 0.14 0.15'
  'draw string 'xmin+xx+dx/2' 'ymax+0.4' 'ctit.i


 i = i + 1
endwhile





xc1 = xmin
xc2 = xc1 + dx*1.5
yc2 = ymax - yy - dy - 0.3
yc1 = yc2 - 0.15

'xcbar 'xc1' 'xc2' 'yc1' 'yc2' -line -edge triangle -fs 1 -fw 0.1 -fh 0.1'

