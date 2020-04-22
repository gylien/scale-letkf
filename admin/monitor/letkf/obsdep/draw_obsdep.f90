program draw_obs_dep
  implicit real(a-h,o-z)
!  character(100) :: cfile='./test_data/obsdep.dat'
  character(100) :: cfile='./test_data_ens/obsdep.dat'
  integer :: nobs
  real(4) :: wk(11)
  integer :: n, iunit

  real(4),parameter :: vplabel(4)=(/1000.0,500.0,200.0,100.0/)
  character*4,parameter :: cplabel(4)=(/'1000','500 ','200 ','100 '/)

  integer,parameter::nobsmax=200000

  integer,parameter::nlevmax=200
  real(4)::vprof_lev(nlevmax)
  real(4)::vprof_obs(nlevmax)
  real(4)::vprof_oer(nlevmax)
  real(4)::vprof_fm(nlevmax)
  real(4)::vprof_fer(nlevmax)
  real(4)::vprof_am(nlevmax)

  integer,parameter::nmem=50
  real(4) :: wk_ext(nmem)

  character*2::cvar 
  character*40::psfile 
  character*20::xtitle 
  character*40::title1 
  character*40::title2(2) 

  integer,parameter :: nsta_max=1000
  real(4) :: vlon_stat(nsta_max),vlat_stat(nsta_max)

  cvar='q'

open(21,file='stations.dat',form='formatted')
do ista=1,nsta_max
  read(21,*,iostat=ios) inum,vlon,vlat
  if (ios.ne.0) exit
  vlon_stat(inum)=vlon
  vlat_stat(inum)=vlat
end do
close(21)

nsta=inum

select case (cvar)
  case('u','U')
  ielm_ref=2819
  vmin=-40.0
  vmax=80.0
  amtics=20.0
  vfact=1.0
  xtitle='U (m/s)'
  case('v','V')
  ielm_ref=2820
  vmin=-40.0
  vmax=40.0
  amtics=20.0
  vfact=1.0
  xtitle='V (m/s)'
  case('t','T')
  ielm_ref=3073
  vmin=200.0
  vmax=300.0
  amtics=20.0
  vfact=1.0
  xtitle='Temperature (K)'
  case('tv','Tv','TV')
  ielm_ref=3074
  vmin=200.0
  vmax=300.0
  amtics=20.0
  vfact=1.0
  xtitle='Virtual temp. (K)'
  case('rh','RH')
  ielm_ref=3071
  vmin=200.0
  vmax=300.0
  amtics=20.0
  vfact=1.0
  xtitle='Relative humidity (%)'
  case('q','Q')
  ielm_ref=3330
  vmin=0.0
  vmax=20.0
  amtics=4.0
  vfact=1000.0
  xtitle='Water vapor (g/kg)'
end select

 do ista=1,nsta
  write(*,*) ista
  write(psfile,'(A,A,A,I3.3)') 'prof_',trim(cvar),'_',ista
  write(title2(1),'(A,F6.2)') 'lon: ',vlon_stat(ista)
  write(title2(2),'(A,F6.2)') 'lat: ',vlat_stat(ista)

  vlon_ref=vlon_stat(ista)
  vlat_ref=vlat_stat(ista)

  iunit=92
  open (iunit, file=cfile, form='unformatted', access='sequential', convert='big_endian')


  ilev=0

  do n = 1, nobsmax
    read(iunit,iostat=ios) wk
    read(iunit,iostat=ios) wk_ext
     if (ios.ne.0) exit

    ielm=wk(1)
    vlon=wk(2)
    vlat=wk(3)
    vlev=wk(4)
    vdat=wk(5)
    verr=wk(6)
    ityp=wk(7)
    vdif=wk(8)
    iqc=wk(9)
    vomb=wk(10)
    voma=wk(11)

 
  if (iqc.eq.0.and.ityp.eq.1 .and. &
!  if (ityp.eq.1 .and. &
     abs(vlon-vlon_ref).le.0.01.and.abs(vlat-vlat_ref).le.0.01.and. &
      ielm.eq.ielm_ref ) then
    ilev=ilev+1
    if (ilev.ge.2.and.vlev.ge.vprof_lev(ilev-1)) exit 
!   write(*,'(I5,I5,3F9.2,2E10.2,I4,F9.2,I4,2E10.2)') &
!         n,ielm,vlon,vlat,vlev,vdat,verr,ityp,vdif,iqc,vomb,voma
    vprof_lev(ilev)=vlev
    vprof_obs(ilev)=vdat
    vprof_oer(ilev)=verr
    vprof_fm(ilev)=vdat-vomb
    vprof_fer(ilev)=sqrt(sum(wk_ext**2)/real(nmem-1))
    vprof_am(ilev)=vdat-voma

  end if

  end do
  close (iunit)

  nlev=ilev

   if (nlev.eq.0) then
    write(*,*) 'no data'
    cycle
    stop
   end if

  vprof_obs=vprof_obs*vfact
  vprof_fm=vprof_fm*vfact
  vprof_am=vprof_am*vfact
  vprof_oer=vprof_oer*vfact
  vprof_fer=vprof_fer*vfact
 
! *** quicklook ***

pbot=1050.0
ptop=100.0

astics=amtics


iout=2
! *** general settings ***
      call gliset ('MSGLEV',1)
      call sgiset ('IFONT',1)
      call swcmll
      call swlset ('LSEP',.TRUE.) ! psfilename numbering
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
      call grswnd (vmin,vmax,pbot,ptop) ! set window

      vpl=0.27; vpr=0.77; vpb=0.23 ; vpt=0.73
      call grsvpt (vpl,vpr,vpb,vpt) ! set viewport

      call grstrn (2) ! linear or log
      call grstrf

! **** Lines & markers ****
      call sglset ('LCLIP',.TRUE.) ! Cliping

      call uuslnt(1)
      call uuslni(3)
      call uulin (nlev,vprof_obs(1:nlev),vprof_lev(1:nlev))
      call uuslnt(1)
      call uuslni(1)
!      call uulin (nlev,vprof_obs(1:nlev)+vprof_oer(1:nlev),vprof_lev(1:nlev))
!      call uulin (nlev,vprof_obs(1:nlev)-vprof_oer(1:nlev),vprof_lev(1:nlev))
do ilev=1,nlev
      vleft=vprof_obs(ilev)-vprof_oer(ilev)
      vright=vprof_obs(ilev)+vprof_oer(ilev)
      call uulin (2,(/vleft,vright/),(/vprof_lev(ilev),vprof_lev(ilev)/))
      call uulin (2,(/vleft,vleft/),(/vprof_lev(ilev)*1.01,vprof_lev(ilev)/1.01/))
      call uulin (2,(/vright,vright/),(/vprof_lev(ilev)*1.01,vprof_lev(ilev)/1.01/))
end do


      call uuslnt(1)
      call uuslni(23)
      call uulin (nlev,vprof_fm(1:nlev),vprof_lev(1:nlev))
      call uuslni(21)
do ilev=1,nlev
      vleft=vprof_fm(ilev)-vprof_fer(ilev)
      vright=vprof_fm(ilev)+vprof_fer(ilev)
      call uulin (2,(/vleft,vright/),(/vprof_lev(ilev),vprof_lev(ilev)/))
      call uulin (2,(/vleft,vleft/),(/vprof_lev(ilev)*1.01,vprof_lev(ilev)/1.01/))
      call uulin (2,(/vright,vright/),(/vprof_lev(ilev)*1.01,vprof_lev(ilev)/1.01/))
end do

      call uuslni(43)
      call uulin (nlev,vprof_am(1:nlev),vprof_lev(1:nlev))

      ! **** x ,y axis ****
      call UYSFMT('(F5.1)')

      call uziset ('INDEXT2',5)
      call uziset ('INDEXL1',5)
      call uziset ('INNER',-1)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)


      call uzlset ('LABELXB',.TRUE.)
      call uzlset ('LABELXT',.FALSE.)

      call uxaxdv ('B',astics,amtics)
      call uxaxdv ('T',astics,amtics)

      call uxsttl ('B',trim(xtitle),0.0)
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.TRUE.)
      call uyaxlb ('L',vplabel,3,vplabel,cplabel,4,3)
      call uyaxlb ('R',vplabel,3,vplabel,cplabel,4,3)
      call uziset ('IROTCYL',1)

      call uysttl ('L','Pressure (hPa)',0.0)
      call sglset ('LCLIP',.FALSE.) ! Cliping
      call sgtxzv (0.5*(vpr+vpl),0.020+vpt,trim(title1),0.018,0,0,5) !
      call sgtxzv (vpr-0.02,0.035+vpt,trim(title2(1)),0.014,0,1,3) !
      call sgtxzv (vpr-0.02,0.015+vpt,trim(title2(2)),0.014,0,1,3) !



      call grcls

end do !!! ista

      stop
end program draw_obs_dep


