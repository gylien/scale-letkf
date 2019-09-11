'reinit'

'set display color white'
'c'

'set grads off'
'set tlsupp off'

ctl = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_NATURE/20010101000000/fcst/mean/QCRG_TOT_d01z-3d.ctl"
uctl = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_NATURE/20010101000000/fcst/mean/U_d01z-3d.ctl"
vctl = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_WK1982_LT_SN14_NATURE/20010101000000/fcst/mean/V_d01z-3d.ctl"



var1 = "qcrg_tot*1.e2"
tit = "Charge density of QHYD (10`a2`n nC m`a-3`n) & vertical vorticity (`3,`00.5 10`a-2`ns`a-1`n)"

'open 'ctl
'open 'uctl
'open 'vctl

'defcolor.gs'


clevs = "-8 -4 -2 -1 -0.5 0.5 1 2 4 8"
ccols = "245 244 242 241 240 0 230 231 233 234 235"


slon = 130 
elon = 250
slat = 130 
elat = 250

* km
level.1 = 6
level.2 = 2

level.1 = 14
level.2 = 10

gt.1 = 7
gt.2 = 13
gt.3 = 19

'set mproj off'

xmin = 0.8
dx = 3
dy = 3
ymax = 7.9
pdx = 0.5
pdy = 0.5

xlint = 20
ylint = 20

jmax = 2
j = 1
while(j <= jmax)
 yy = (dy + pdy) * (j - 1)

 imax = 3
 i = 1
 while(i <= imax)
  xx = (dx + pdx) * (i - 1)

  'set parea 'xmin+xx' 'xmin+xx+dx' 'ymax-yy-dy' 'ymax-yy

  'set lon 'slon-5' 'elon+5
  'set lat 'slat-5' 'elat+5

  'dudy = cdiff(u.2,y)/4000'
  'dvdx = cdiff(v.3,x)/4000'
  'vort = dvdx - dudy'


  'set lon 'slon' 'elon
  'set lat 'slat' 'elat

  'set lev 'level.j
  clevel = subwrd(result, 4)
  'set t 'gt.i
  ft = (gt.i - 1 ) * 5
*  'q time'
*  ctime = subwrd(result, 3)

  'set xlint 'xlint
  'set ylint 'ylint

  'set gxout shade2'

  'set xlab on'
  'set ylab on'
  'set grid on'
  'set frame on'
 
  'set clevs 'clevs 
  'set ccols 'ccols 

  'd 'var1

  'set xlab off'
  'set ylab off'
  'set grid off'
  'set frame off'

  'set gxout contour'
  'set clevs -0.5 0.5'
  'set ccolor 1'
  'set clab off'
  'd vort*1.e2'


  ptit = "FT="ft"min, z="clevel"km"

  'set string 1 c 3'
  'set strsiz 0.1 0.12'
  'draw string 'xmin+xx+dx/2' 'ymax-yy+0.1' 'ptit
  

  i = i + 1
 endwhile

 j = j + 1
endwhile

xc1 = xmin
xc2 = xc1 + dx*1.5
yc2 = ymax - yy - dy - 0.3
yc1 = yc2 - 0.15

'xcbar 'xc1' 'xc2' 'yc1' 'yc2' -line -edge triangle -fs 1 -fw 0.1 -fh 0.1'

'set string 1 l 4'
'set strsiz 0.12 0.14'
'draw string 'xmin' 'ymax+0.4' 'tit


