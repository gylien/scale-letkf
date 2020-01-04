!==================================================!
module setup
implicit real(a-h,o-z)

integer,parameter::nrec=20
real(4)::bias_a(nrec)
real(4)::rmse_a(nrec)
character*5::cdate(nrec)
integer::ndata(nrec)

real(4)::bias_g(nrec)
real(4)::rmse_g(nrec)


character*30::cfile_a='../data_inv/q_anal.txt'
character*30::cfile_g='../data_inv/q_gues.txt'

character*30::ctitle1='Q'
character*30::ctitle2='(g/kg)'

character*30::psfile='../figs/q_letkf'

real(4),parameter::xmax_window=6.0 !!! day
real(4),parameter::vmin_rmse=0.0
real(4),parameter::vmax_rmse=4.0
real(4),parameter::btic_rmse=1.0
real(4),parameter::vmin_bias=-2.0
real(4),parameter::vmax_bias= 2.0
real(4),parameter::btic_bias=1.0
real(4),parameter::vmin_num= 0.0
real(4),parameter::vmax_num=1500.0
real(4),parameter::btic_num=500.0

real(4),parameter::factor=1000.0


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



return
end subroutine set_time
!==================================================!
subroutine load
use setup

integer,parameter::nrec_max=120

real(4)::rread1(nrec_max)
real(4)::rread2(nrec_max)
character*14::cread(nrec_max)
integer::iread(nrec_max)

character*50::cdumy
character*14::ctime

open(11,file=trim(cfile_a),form='formatted')
 do i=1,nrec
  read(11,*) ctime,bias_a(i),rmse_a(i),idumy
 end do 
close(11)

open(12,file=trim(cfile_g),form='formatted')
 ios=0
 do i=1,nrec
  read(12,*) cread(i),bias_g(i),rmse_g(i),ndata(i)
   if (cread(i)(9:10).eq.'00')then
    write(cdate(i),'(A2,A1,A2)')cread(i)(5:6),'/',cread(i)(7:8)  
   else
    cdate(i)=''
   end if
 end do
close(12)

bias_a(1:nrec)=bias_a(nrec:1:-1) * factor
rmse_a(1:nrec)=rmse_a(nrec:1:-1) * factor
bias_g(1:nrec)=bias_g(nrec:1:-1) * factor
rmse_g(1:nrec)=rmse_g(nrec:1:-1) * factor
ndata(1:nrec) =ndata(nrec:1:-1)
cdate(1:nrec) =cdate(nrec:1:-1)


return
end subroutine load
!==================================================!
subroutine draw
use setup

real(4)::xlocs(nrec)


! *** quicklook ***

xmin=0.0

xlocs=(/( xmin + 0.25*real(i-1) ,i=1,nrec )/)


iout=2
! *** general settings ***
      call sgiset ('IFONT',1)
      call swcmll
      call swlset ('LSEP',.FALSE.) ! psfilename numbering
      call swcset ('FNAME',trim(psfile))
      call swiset ('IWIDTH',1000) ! output window size
      call swiset ('IHEIGHT',800)     

      call swiset('IFL',1) !!! PNG
!   call swiset('IFL',2) !!! EPS
!    call swiset('IFL',4) !!! PDF

      call gropn(iout) 
      call sglset ('LFULL',.TRUE.) ! using fullsize
      call slmgn (0.0,0.0,0.0,0.0) ! margin 

      call grfrm 
      call grswnd (xmin,xmax_window,vmin_rmse,vmax_rmse) ! set window

      vpl=0.17; vpr=0.92; vpb=0.45 ; vpt=0.73
      call grsvpt (vpl,vpr,vpb,vpt) ! set viewport

      call grstrn (1) ! linear or log
      call grstrf

! **** Lines & markers ****
      call sglset ('LCLIP',.TRUE.) ! Cliping
     
      call uulinz (nrec,xlocs,rmse_a,1,5)
      call uulinz (nrec,xlocs,rmse_g,3,5)
      nls=int((xmax_window-xmin) / 0.25) +1

      do il=1,nls
       xloc=real(int(xmin))+0.25* real(il-1)
!       write(*,*) il,xloc
       ithck=1
       if (il.le.nrec.and.cdate(il).ne.'') ithck=3
       call uulinz (2,(/xloc,xloc/),(/vmin_rmse,vmax_rmse/),3,40+ithck)
      end do

!      write(*,*) maxval(bias_a),minval(bias_a)
!      write(*,*) maxval(rmse_a),minval(rmse_a)
      
      ! **** x ,y axis ****
      call UYSFMT('(F5.1)')

      call uziset ('INDEXT2',5)
      call uziset ('INDEXL1',5)
      call uziset ('INNER',-1)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)


      call uzlset ('LABELXB',.FALSE.)
      call uzlset ('LABELXT',.FALSE.)

      call uxaxlb ('B',xlocs,nrec,xlocs,cdate,5,nrec)
      call uxaxlb ('T',xlocs,nrec,xlocs,cdate,5,nrec)
      
!      call uxsttl ('B','JST',0.0)

      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)
      call uyaxdv ('L',btic_rmse,btic_rmse)
      call uyaxdv ('R',btic_rmse,btic_rmse)  
      call uziset ('IROTCYL',1)

      call uysttl ('L','RMSE',0.0)
      call sglset ('LCLIP',.FALSE.) ! Cliping
      call sgtxzv (0.5*(vpl+vpr),0.035+vpt,trim(ctitle1),0.025,0,0,5) !
      call sgtxzv (vpr-0.02,0.025+vpt,trim(ctitle2),0.018,0,1,5) !


vmin=vmin_bias
vmax=vmax_bias
bstics=btic_bias
bmtics=btic_bias
!      call grfrm 
      call grswnd (xmin,xmax_window,vmin,vmax) ! set window

      vpl=0.17; vpr=0.92; vpb=0.13 ; vpt=0.41
      call grsvpt (vpl,vpr,vpb,vpt) ! set viewport

      call grstrn (1) ! linear or log
      call grstrf

! **** Lines & markers ****
      call sglset ('LCLIP',.TRUE.) ! Cliping
     
      call uulinz (nrec,xlocs,bias_a,1,5)
      call uulinz (nrec,xlocs,bias_g,3,5)

      nls=int((xmax_window-xmin) / 0.25) +1

      do il=1,nls
       xloc=real(int(xmin))+0.25* real(il-1)
!       write(*,*) il,xloc
       ithck=1
       if (il.le.nrec.and.cdate(il).ne.'') ithck=3
       call uulinz (2,(/xloc,xloc/),(/vmin,vmax/),3,40+ithck)
      end do

      ! **** x ,y axis ****

      call uzinit
      call UYSFMT('(F5.1)')
      call uziset ('INDEXT2',5)
      call uziset ('INDEXL1',5)
      call uziset ('INNER',-1)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)


      call uzlset ('LABELXB',.TRUE.)
!      call uzlset ('LABELXB',.FALSE.)
      call uzlset ('LABELXT',.FALSE.)

      call uziset ('ICENTXB',0) !!! centering
!      call uxaxdv ('B',astics,amtics)
!      call uxaxdv ('T',astics,amtics)

!      call ucxacl ('B', 20190401, 6)
!      call ucxacl ('T', 20190401, 6)

      call uxaxlb ('B',xlocs,nrec,xlocs,cdate,5,nrec)
      call uxaxlb ('T',xlocs,nrec,xlocs,cdate,5,nrec)
      
!      call uxsttl ('B','JST',0.0)

      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)

      call uyaxdv ('L',bstics,bmtics)
      call uyaxdv ('R',bstics,bmtics)  
      call uziset ('IROTCYL',1)

      call uysttl ('L','BIAS',0.0)
      call sglset ('LCLIP',.FALSE.) ! Cliping

      call grcls




return
end subroutine draw
!==================================================!
