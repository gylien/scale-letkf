
'reinit'

'set display color white'
'c'

'set grads off'
'set tlsupp off'

top = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT"

exp1 = "2000m_WK1982_NOLT_TOMITA_NORHO"
exp2 = "2000m_WK1982_NOLT_TOMITA_RHO"
exp3 = "2000m_WK1982_NOLT_SN14"

tit.1 = "Tomita"
tit.2 = "Tomita (Roh)"
tit.3 = "SN14"

ctl1 = ""top"/"exp1"/20010101000000/fcst/mean/QHYD_d01z-3d.ctl "
ctl2 = ""top"/"exp2"/20010101000000/fcst/mean/QHYD_d01z-3d.ctl "
ctl3 = ""top"/"exp3"/20010101000000/fcst/mean/QHYD_d01z-3d.ctl "

'open 'ctl1
'open 'ctl2
'open 'ctl3


gt = 25

gts = 1
gte = 25
dgt = 1

gt = gts
while(gt <= gte)

*gz.1 = 6
*gz.2 = 2
*
*gz.1 = 2
*gz.2 = 1
*
*gz.1 = 4
*gz.2 = 3
*

gz.1 = 6
gz.2 = 2
gz.1 = 8
gz.2 = 7


*glev.1 = 17.25
*glev.2 = 4.25

'set mproj off'

'set t 'gt

lons = 105
lone = 295

lats = 105
late = 295

cmin = 0.1
clevs = ""cmin" 0.5 1 2 3 4 5 6 7 8 9 10"

xmin = 0.7
dx = 3
dy = 3
ymax = 7.5
pdx = 0.5
pdy = 0.5

imax = 3
jmax = 2

j = 1
while(j <= jmax)
 yy = (dy + pdy) * (j - 1)

 i = 1
 while(i <= imax)
  xx = (dx + pdx) * (i - 1)
 
  'set parea 'xmin+xx' 'xmin+xx+dx' 'ymax-dy-yy' 'ymax-yy
 
  'set z 'gz.j
*  'set lev 'glev.j
  clev = subwrd(result,4)
 
  'set gxout shade2'
  
  'set lon 'lons' 'lone
  'set lat 'lats' 'late

  'set xlint 20' 
  'set ylint 20' 

  'set grid on'
  'set frame on'
  'set xlab on'
  'set ylab on'
  
  'set clevs 'clevs
  
*  'var = qhyd.'i'*1.e3'
*  'd maskout(var,var-'cmin')'
  'd qhyd.'i'*1.e3'
  
  'set grid off'
  'set frame off'
  'set xlab off'
  'set ylab off'
  
  'set string 1 c 3'
  'set strsiz 0.12 0.14'
 
  ptit = "z="clev"km"
  'draw string 'xmin+xx+dx/2' 'ymax+0.1-yy' 'ptit
 
  if (j=1)
   'draw string 'xmin+xx+dx/2' 'ymax+0.3' 'tit.i
  endif

 i = i + 1
 endwhile

j = j + 1
endwhile

xc1 = xmin
xc2 = xc1 + dx*1.5
yc2 = ymax - yy - dy - 0.3
yc1 = yc2 - 0.15

'xcbar -line 'xc1' 'xc2' 'yc1' 'yc2' -edge triangle -fw 0.1 -fh 0.1'


  
'set string 1 c 3'
'set strsiz 0.16 0.18'
gtt = gt - 1
'draw string 'xmin+dx*1.5 + pdx' 'ymax + 0.7' Hydrotemeteors (g/kg) t='gtt*300's'

gtt3 = gtt
if (gtt < 10)
 gtt3 = "00"gtt""
else
 if (gtt < 100)
  gtt3 = "0"gtt""
 endif
endif
ofig = "png/6p_qh_t"gtt3"s_"gz.1".png"

'gxprint tmp.png'
'!convert -trim tmp.png 'ofig


'reset'
'set display color white'
'c'
'set grads off'
'set tlsupp off'


gt = gt + 1
endwhile

