!==================================================!
module setup
implicit real(a-h,o-z)

integer,parameter::nplot=360 !!! every 5sec for 20min = 240
integer::irecs(nplot)
character*14::ctimes(nplot)
real(4)::dt_obs(nplot)
real(4)::dt_d4_ana(nplot)
real(4)::dt_d4_fcst(nplot)

integer,parameter::nintvsec=5

real(4),parameter::xmin_window=0.0 !!! min
real(4),parameter::xmax_window=30.0 !!! min
real(4),parameter::xxs(nplot)=(/( xmin_window + (xmax_window-xmin_window)*real(i-1)/real(nplot)   ,i=1,nplot)/)


integer,parameter::nmin=30 !!! min
real(4),parameter::xxsh(nmin)=(/(real(i-1), i=1,nmin)/)
character*5::clabels(7)
real(4)::xlabels(5)
integer::nlabel

integer::itime_now
real(4)::xxnow

character*30::cfile='./data_inv/d4.txt'

character*30::ctitle1='Domain 4'
character*30::ctitle2=''

character*30::psfile='./figs/d4'

real(4),parameter::vmin=-20.0
real(4),parameter::vmax= 10.0
real(4),parameter::btic=5.0


end module setup
!==================================================!
program main
use setup

call set_time

call load
call draw

stop
end program main
!==================================================!
subroutine set_time
use setup

character*14::ctime_now

call dateq1(id)
call timeq1(it_now)

call date13(id,iyyyy,imm,idd)


write(ctitle2,'(I4.4,A,I2.2,A,I2.2)') iyyyy,'/',imm,'/',idd



it = (it_now /100) * 100 
itres=mod(it_now,100)

it = it - 90000  - 1200 !!! UTC 
if (it.lt.0.0)then
 it=it+240000
 call datef1(-1,id,id)
end if

call datef1(-1,id,id)
 it=it-nintvsec  

nlabel=0
do iplot=1,nplot
 it=it+nintvsec 
 if (mod(it,100).eq.60) it=it+40
 if (mod(it,10000).eq.6000) it=it+4000
 if (it/10000.eq.24) then
  it=it-240000
  call datef1(1,id,id)
 end if

write(ctimes(iplot),'(I8.8,I6.6)') id,it

if (mod(it,500).eq.0)then !!! label 5min
 nlabel=nlabel+1
 write(clabels(nlabel),'(I2.2)')mod(int(it/100),100)
 if (trim(clabels(nlabel)) .eq."00") write(clabels(nlabel),'(I2.2,A1,I2.2)')it/10000,':',mod(int(it/100),100)
 xlabels(nlabel)=xxs(iplot)
end if
end do


ires_sec = (itres/100) * 60 + mod(itres,100)

xxnow=xxsh(21) + (xxsh(21)-xxsh(20)) * real(ires_sec)/real(3600)


return
end subroutine set_time
!==================================================!
subroutine load
use setup

integer,parameter::nrecmax=300

character*14:: ctrecread, ctobsread, ctd4aread, ctd4fread

irecs=0

open(11,file=trim(cfile),form='formatted')
 do i=1,nrecmax
  read(11,*,iostat=ios) ctrecread, ctobsread, ctd4aread, ctd4fread
 if (ios.ne.0)then
  exit
 else
  iloc=0
 do iplot=1,nplot
  if (ctrecread.eq.ctimes(iplot))then
   iloc=iplot
   irecs(iplot)=1
   call dtime_min (ctobsread, ctrecread, dt_obs(iplot)) 
   call dtime_min (ctd4aread, ctrecread, dt_d4_ana(iplot)) 
   call dtime_min (ctd4fread, ctrecread, dt_d4_fcst(iplot) ) 
  end if
 end do
 if (iloc.eq.0)then
!  write(*,*) 'check :: itrecread',itrecread
!  write(*,*)ctrecread
!  do iplot=1,nplot
!   write(*,*) ctimes(iplot)
!  end do
 exit
!  stop
 end if

 end if

 end do 
close(11)


return
end subroutine load
!==================================================!
subroutine dtime_min(ctime1,ctime2,dtmin)
use setup

character*14::ctime1,ctime2

read(ctime1(1:8),*)id1
read(ctime2(1:8),*)id2

read(ctime1(9:14),*)it1
read(ctime2(9:14),*)it2

!write(*,*) id1,id2
!write(*,*) it1,it2

if (id1.ne.id2) then
 call dateg1(nd,id2,id1)
 it1=it1+240000*nd
end if

idh=it1/10000-it2/10000
idm=mod(it1/100,100)-mod(it2/100,100)
ids=mod(it1,100)-mod(it2,100)

dtmin=real(idh)*60.0 + real(idm) + real(ids)/60.0

return
end subroutine dtime_min
!==================================================!
subroutine draw
use setup

character*1::csgi

! *** quicklook ***

iout=2

! *** general settings ***
      call sgiset ('IFONT',1)
      call swcmll
      call swlset ('LSEP',.FALSE.) ! psfilename numbering
      call swcset ('FNAME',trim(psfile))
      call swiset ('IWIDTH',1000) ! output window size
      call swiset ('IHEIGHT',800)     

      call swiset('IFL',1) !!! PNG
!    call swiset('IFL',2) !!! EPS
!    call swiset('IFL',4) !!! PDF

      call gropn(iout) 
      call sglset ('LFULL',.TRUE.) ! using fullsize
      call slmgn (0.0,0.0,0.0,0.0) ! margin 

      call grfrm 
      call grswnd (xmin_window,xmax_window,vmin,vmax) ! set window

      vpl=0.17; vpr=0.92; vpb=0.23 ; vpt=0.73
      call grsvpt (vpl,vpr,vpb,vpt) ! set viewport

      call grstrn (1) ! linear or log
      call grstrf

! **** Lines & markers ****
      call sglset ('LCLIP',.TRUE.) ! Cliping
     

      do iplot=1,nplot-1
       if (irecs(iplot).eq.1.and.irecs(iplot+1).eq.1)then
          call uulinz (2,xxs(iplot:iplot+1),dt_obs(iplot:iplot+1)+0.05,1,24)
          call uulinz (2,xxs(iplot:iplot+1),dt_d4_ana(iplot:iplot+1)-0.05,1,34)
          call uulinz (2,xxs(iplot:iplot+1),dt_d4_fcst(iplot:iplot+1),1,44)
       end if
      end do

      do il=1,nmin
       xloc=real(il)
       ithck=1
!       if (il.le.nrec.and.cdate(il).ne.'') ithck=3
       call uulinz (2,(/xloc,xloc/),(/vmin,vmax/),3,40+ithck)
      end do

       call uulinz (2,(/xxnow,xxnow/),(/vmin,vmax/),3,20+ithck)

       call uulinz (2,(/xmin_window,xmax_window/),(/0.0,0.0/),3,1)

      
      ! **** x ,y axis ****
      call sglset ('LCLIP',.FALSE.) ! Cliping
     
      call UYSFMT('(I3)')

      call uziset ('INDEXT2',5)
      call uziset ('INDEXL1',5)
      call uziset ('INNER',-1)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)

      call uzlset ('LABELXB',.TRUE.)
      call uzlset ('LABELXT',.FALSE.)

      call uxaxlb ('B',xlabels(1:nlabel),nlabel,xlabels(1:nlabel),clabels(1:nlabel),5,nlabel)
      call uxaxlb ('T',xlabels(1:nlabel),nlabel,xlabels(1:nlabel),clabels(1:nlabel),5,nlabel)

      
!      call uxsttl ('B','JST',0.0)

      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)
      call uyaxdv ('L',btic,btic)
      call uyaxdv ('R',btic,btic)  
      call uziset ('IROTCYL',1)

      call uysttl ('L',csgi(131)//'t (min)',0.0)
      call sgtxzv (0.5*(vpl+vpr),0.035+vpt,trim(ctitle1),0.025,0,0,5) !
      call sgtxzv (vpr-0.02,0.025+vpt,trim(ctitle2),0.016,0,1,3) ! 

      call sgtxzv (vpr-0.01,vpt-0.03,'D4 fcst',0.014,0,1,44) !       
      call sgtxzv (vpr-0.01,vpt-0.28,'D4 analysis',0.014,0,1,34) !       
      call sgtxzv (vpr-0.01,vpt-0.34,'PAWR obs',0.014,0,1,24) !       

      call grcls



return
end subroutine draw
!==================================================!
