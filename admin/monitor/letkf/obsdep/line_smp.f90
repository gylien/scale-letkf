program line

implicit real(a-h,o-z)

integer,parameter::n=1000
integer,parameter::k=6
real(4),parameter::twpi=6.283
real(4)::y(n+1)
real(4)::axx(n+1)

integer,parameter::nmode=16

character*15::cfile="data.txt"

nw = n/nmode/2 


nl=max(1+(k-1)*nw*2,1)
nr=min(2*nw+(k-1)*nw*2,n)

!nl = max( k*nw*2 -nw+1 -nw,1)
!nr = min( k*nw*2 +nw+1,n)

do i=1,n+1
  axx(i) = real(i-1)/real(n)
end do
do i=1,n+1
!    y(i) = sin(twpi * axx(i)*real(k))* 0.2
!    y(i) = cos(twpi * axx(i)*real(k))* 0.2
 if (i.ge.nl.and.i.le.nr )then
write(*,*)(axx(i)-axx(nl+nw))
    y(i) = exp ( -((axx(i)-axx(nl+nw))/ (real(nw)/real(2*n)) )**2 ) * 0.2
 else
    y(i)=0.0
 end if
end do

! *** quicklook ***

xmin=0.0
xmax=1.0

vmin=-1.0
vmax=1.0
bstics=1.0
bmtics=1.0
astics=1.0
amtics=1.0

base=0.0

iout=2
! *** general settings ***
      call sgiset ('IFONT',1)
      call swcmll
      call swlset ('LSEP',.TRUE.) ! psfilename numbering
      call swcset ('FNAME','figure')
      call swiset ('IWIDTH',1000) ! output window size
      call swiset ('IHEIGHT',800)     

    call swiset('IFL',1) !!! PNG
!   call swiset('IFL',2) !!! EPS
!    call swiset('IFL',4) !!! PDF

      call gropn(iout) 
      call sglset ('LFULL',.TRUE.) ! using fullsize
      call slmgn (0.0,0.0,0.0,0.0) ! margin 

      call grfrm 
      call grswnd (xmin,xmax,vmin,vmax) ! set window

      vpl=0.17; vpr=0.77; vpb=0.23 ; vpt=0.73
      call grsvpt (vpl,vpr,vpb,vpt) ! set viewport

      call grstrn (1) ! linear or log
      call grstrf

! **** Lines & markers ****
      call sglset ('LCLIP',.TRUE.) ! Cliping

      call uuslnt(1)
      call uuslni(3)
      call uulin (n+1,axx(1:n+1),y(1:n+1))

      
      ! **** Box ****
      
!      rsize=0.06*1.0/6.0
      
!      call uvbxaz (n,axtimeh,0.0,prec,999,0)
!       call uvbxaz (n-1,axtime(1:n-1),base(2:n),prec(2:n),4999,999,0.0)!rsize)
      
      ! **** x ,y axis ****
      call UYSFMT('(F5.1)')

      call uziset ('INDEXT2',5)
      call uziset ('INDEXL1',5)
      call uziset ('INNER',-1)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)


!      call uzlset ('LABELXB',.TRUE.)
!      call uzlset ('LABELXT',.FALSE.)

!      call uxaxdv ('B',astics,amtics)
!      call uxaxdv ('T',astics,amtics)

 !     call uxsttl ('B','JST',0.0)

!      call uzlset ('LABELYR',.FALSE.)
!      call uzlset ('LABELYL',.TRUE.)
!      call uyaxdv ('L',bstics,bmtics)
!      call uyaxdv ('R',bstics,bmtics)  
!      call uziset ('IROTCYL',1)

!      call uysttl ('L','Precip (mm/h)',0.0)
      call sglset ('LCLIP',.FALSE.) ! Cliping
!      call sgtxzv (0.47,0.76,'Kyotanabe 2013/07/13',0.025,0,0,5) !


      call grcls


      stop
      
    end program line
