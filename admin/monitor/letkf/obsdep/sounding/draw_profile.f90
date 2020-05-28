program draw_obs_dep
  implicit real(a-h,o-z)
  character(100) :: cdirbase='/work/hp150019/share/SCALE-LETKF-rt/result/ope/d1/'
  character(100) :: cfile
  integer :: nobs
  real(4) :: wk(11)
  integer :: n, iunit

  real(4),parameter :: vplabel(5)=(/1000.0,500.0,300.0,200.0,100.0/)
  character*4,parameter :: cplabel(5)=(/'1000','500 ','300 ','200 ','100 '/)

  integer,parameter::nobsmax=200000

  integer,parameter::nprofmax=500
  integer,parameter::nlevmax=200
  real(4)::vprof_lev(nlevmax,nprofmax)
  real(4)::vprof_obs(nlevmax,nprofmax)
  real(4)::vprof_oer(nlevmax,nprofmax)
  real(4)::vprof_fm(nlevmax,nprofmax)
  real(4)::vprof_fer(nlevmax,nprofmax)
  real(4)::vprof_am(nlevmax,nprofmax)
  integer::iprof_nlev(nprofmax)
  real(4)::vlons(nprofmax)
  real(4)::vlats(nprofmax)

  integer,parameter::nmem=50
  real(4) :: wk_ext(nmem)

  character*2::cvar 
  character*40::psbase='profile' 
  character*40::psfile 
  character*20::xtitle 
  character*40::title1 
  character*40::title2(2) 
  character*40::title3 

  integer,parameter :: nsta_max=1000
  real(4) :: vlon_stat(nsta_max),vlat_stat(nsta_max)

  character*14::cread
  character*10::ctime
  real(4),parameter::range_lonl=95.0
  real(4),parameter::range_lonr=175.0
  real(4),parameter::range_latl=10.0
  real(4),parameter::range_latr=55.0
  integer,parameter::iclrmap=4

  integer,parameter::nptsh_d1=560
  integer,parameter::npts_d1=nptsh_d1*2+1
  real(4)::vlons_area_d1(npts_d1) 
  real(4)::vlats_area_d1(npts_d1) 

  real(4),parameter::vmin_dtemp=-4.0 
  real(4),parameter::vmax_dtemp=4.0 
  real(4),parameter::tics_dtemp=2.0 



  character*3::cnum3


  nargs=iargc()

  if(nargs.lt.2)then
    write(*,*) 'usage:: ./draw YYYYMMDDHH [u/v/t/tv/q]'
    stop
  end if

  call getarg(1,cread)
  ctime=cread(1:10)
  call getarg(2,cread)
  cvar=trim(cread)

  cfile=trim(cdirbase)//ctime//'0000/obs/obsdep.dat'
!  cfile='../test_data_ens/obsdep.dat'

title1=ctime(5:6)//'/'//ctime(7:8)//' '//ctime(9:10)//'Z'
if (ctime(5:5).eq.'0') title1=ctime(6:6)//'/'//ctime(7:8)//' '//ctime(9:10)//'Z'

!open(21,file='stations.dat',form='formatted')
!do ista=1,nsta_max
!  read(21,*,iostat=ios) inum,vlon,vlat
!  if (ios.ne.0) exit
!  vlon_stat(inum)=vlon
!  vlat_stat(inum)=vlat
!end do
!close(21)

!nsta=inum

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

! do ista=1,nsta
!  write(*,*) ista
!  write(psfile,'(A,A,A,I3.3)') 'prof_',trim(cvar),'_',ista
!  write(title2(1),'(A,F6.2)') 'lon: ',vlon_stat(ista)
!  write(title2(2),'(A,F6.2)') 'lat: ',vlat_stat(ista)

!  vlon_ref=vlon_stat(ista)
!  vlat_ref=vlat_stat(ista)

  iunit=92
  open (iunit, file=trim(cfile), form='unformatted', access='sequential', convert='big_endian')
  ilev=0
  iprof=0

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
!     abs(vlon-vlon_ref).le.0.01.and.abs(vlat-vlat_ref).le.0.01.and. &
      ielm.eq.ielm_ref ) then
    if (iprof.eq.0) then
      iprof=iprof+1
      ilev=ilev+1
    else
      if (vlon.eq.vlons(iprof).and.vlat.eq.vlats(iprof).and.vlev.lt.vprof_lev(ilev,iprof))then
        ilev=ilev+1
      else
        if (ilev.eq.1)then
          continue !!! overwrite (1,iprof)   
        else
          iprof=iprof+1
          ilev=1
        end if
      end if 
    end if
!   write(*,'(I5,I5,3F9.4,2E10.2,I4,F9.2,I4,2E10.2)') &
!         n,ielm,vlon,vlat,vlev,vdat,verr,ityp,vdif,iqc,vomb,voma
    vprof_lev(ilev,iprof)=vlev
    vprof_obs(ilev,iprof)=vdat
    vprof_oer(ilev,iprof)=verr
    vprof_fm(ilev,iprof)=vdat-vomb
    vprof_fer(ilev,iprof)=sqrt(sum(wk_ext**2)/real(nmem-1))
    vprof_am(ilev,iprof)=vdat-voma
    iprof_nlev(iprof)=ilev
    if (ilev.eq.1)then
      vlons(iprof)=vlon
      vlats(iprof)=vlat
    end if
  end if

  end do
  close (iunit)

  nprof=iprof

  write(*,*) '# of profs : ',nprof

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


do iprof=1,nprof
     write(psfile,'(4A,I3.3)')trim(psbase),'_',trim(cvar),'_',iprof
     write(title2(1),'(A,F6.2)') 'lon: ',vlons(iprof)
     write(title2(2),'(A,F6.2)') 'lat: ',vlats(iprof)
nlev=iprof_nlev(iprof)
!
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

if ( ielm_ref.eq.3073.or.ielm_ref.eq.3074) then
      vpl=0.27; vpr=0.62; vpb=0.23 ; vpt=0.73
else
      vpl=0.27; vpr=0.77; vpb=0.23 ; vpt=0.73
end if
      call grsvpt (vpl,vpr,vpb,vpt) ! set viewport

      call grstrn (2) ! linear or log
      call grstrf

! **** Lines & markers ****
      call sglset ('LCLIP',.TRUE.) ! Cliping

      call uuslnt(1)
      call uuslni(3)
      call uulin (nlev,vprof_obs(1:nlev,iprof),vprof_lev(1:nlev,iprof))
      call uuslnt(1)
      call uuslni(1)
!      call uulin (nlev,vprof_obs(1:nlev)+vprof_oer(1:nlev),vprof_lev(1:nlev))
!      call uulin (nlev,vprof_obs(1:nlev)-vprof_oer(1:nlev),vprof_lev(1:nlev))
do ilev=1,nlev
      vleft=vprof_obs(ilev,iprof)-vprof_oer(ilev,iprof)
      vright=vprof_obs(ilev,iprof)+vprof_oer(ilev,iprof)
      call uulin (2,(/vleft,vright/),(/vprof_lev(ilev,iprof),vprof_lev(ilev,iprof)/))
      call uulin (2,(/vleft,vleft/),(/vprof_lev(ilev,iprof)*1.01,vprof_lev(ilev,iprof)/1.01/))
      call uulin (2,(/vright,vright/),(/vprof_lev(ilev,iprof)*1.01,vprof_lev(ilev,iprof)/1.01/))
end do

      call uuslnt(1)
      call uuslni(43)
      call uulin (nlev,vprof_fm(1:nlev,iprof),vprof_lev(1:nlev,iprof))
      call uuslni(41)
do ilev=1,nlev
      vleft=vprof_fm(ilev,iprof)-vprof_fer(ilev,iprof)
      vright=vprof_fm(ilev,iprof)+vprof_fer(ilev,iprof)
      call uulin (2,(/vleft,vright/),(/vprof_lev(ilev,iprof),vprof_lev(ilev,iprof)/))
      call uulin (2,(/vleft,vleft/),(/vprof_lev(ilev,iprof)*1.01,vprof_lev(ilev,iprof)/1.01/))
      call uulin (2,(/vright,vright/),(/vprof_lev(ilev,iprof)*1.01,vprof_lev(ilev,iprof)/1.01/))
end do

      call uuslni(23)
      call uulin (nlev,vprof_am(1:nlev,iprof),vprof_lev(1:nlev,iprof))

      if (vmin*vmax .gt.0.0) call uulinz (2,(/0.0,0.0/),(/pbot,ptop/),3,1)

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
      call uyaxlb ('L',vplabel,5,vplabel,cplabel,4,5)
      call uyaxlb ('R',vplabel,5,vplabel,cplabel,4,5)
      call uziset ('IROTCYL',1)

      call uysttl ('L','Pressure (hPa)',0.0)
      call sglset ('LCLIP',.FALSE.) ! Cliping
      call sgtxzv (0.5*(vpr+vpl),0.020+vpt,trim(title1),0.018,0,0,5) !
      call sgtxzv (vpr,0.035+vpt,trim(title2(1)),0.014,0,1,3) !
      call sgtxzv (vpr,0.015+vpt,trim(title2(2)),0.014,0,1,3) !

 if (ielm_ref.eq.3073.or.ielm_ref.eq.3074)then

      !call grfrm 
      call grswnd (vmin_dtemp,vmax_dtemp,pbot,ptop) ! set window

      vpl=0.65; vpr=0.77; vpb=0.23; vpt=0.73

      call grsvpt (vpl,vpr,vpb,vpt) ! set viewport

      call grstrn (2) ! linear or log
      call grstrf

! **** Lines & markers ****
      call sglset ('LCLIP',.TRUE.) ! Cliping

      call uuslnt(1)
      call uuslni(3)
      call uulin (nlev,vprof_obs(1:nlev,iprof)-vprof_obs(1:nlev,iprof),vprof_lev(1:nlev,iprof))
      call uuslnt(1)
      call uuslni(1)
!      call uulin (nlev,vprof_obs(1:nlev)+vprof_oer(1:nlev),vprof_lev(1:nlev))
!      call uulin (nlev,vprof_obs(1:nlev)-vprof_oer(1:nlev),vprof_lev(1:nlev))
do ilev=1,nlev
      vleft=-vprof_oer(ilev,iprof)
      vright=vprof_oer(ilev,iprof)
      call uulin (2,(/vleft,vright/),(/vprof_lev(ilev,iprof),vprof_lev(ilev,iprof)/))
      call uulin (2,(/vleft,vleft/),(/vprof_lev(ilev,iprof)*1.01,vprof_lev(ilev,iprof)/1.01/))
      call uulin (2,(/vright,vright/),(/vprof_lev(ilev,iprof)*1.01,vprof_lev(ilev,iprof)/1.01/))
end do

      call uuslnt(1)
      call uuslni(43)
      call uulin (nlev,vprof_fm(1:nlev,iprof)-vprof_obs(1:nlev,iprof),vprof_lev(1:nlev,iprof))
      call uuslni(41)
do ilev=1,nlev
      vleft=vprof_fm(ilev,iprof)-vprof_obs(ilev,iprof)-vprof_fer(ilev,iprof)
      vright=vprof_fm(ilev,iprof)-vprof_obs(ilev,iprof)+vprof_fer(ilev,iprof)
      call uulin (2,(/vleft,vright/),(/vprof_lev(ilev,iprof),vprof_lev(ilev,iprof)/))
      call uulin (2,(/vleft,vleft/),(/vprof_lev(ilev,iprof)*1.01,vprof_lev(ilev,iprof)/1.01/))
      call uulin (2,(/vright,vright/),(/vprof_lev(ilev,iprof)*1.01,vprof_lev(ilev,iprof)/1.01/))
end do

      call uuslni(23)
      call uulin (nlev,vprof_am(1:nlev,iprof)-vprof_obs(1:nlev,iprof),vprof_lev(1:nlev,iprof))

      call uulinz (2,(/0.0,0.0/),(/pbot,ptop/),3,1)

      call sglset ('LCLIP',.FALSE.) ! Cliping

      ! **** x ,y axis ****
      call UYSFMT('(F5.1)')

      call uzinit
      call uziset ('INDEXT2',5)
      call uziset ('INDEXL1',5)
      call uziset ('INNER',-1)
      call uzrset ('RSIZEL1',0.018)
      call uzrset ('RSIZEC1',0.018)
      call uzrset ('RSIZET2',0.010)
      call uzrset ('RSIZET1',0.004)


      call uzlset ('LABELXB',.TRUE.)
      call uzlset ('LABELXT',.FALSE.)

      call uxaxdv ('B',tics_dtemp,tics_dtemp)
      call uxaxdv ('T',tics_dtemp,tics_dtemp)

      call uxsttl ('B','diff.',0.0)
      call uzlset ('LABELYR',.FALSE.)
      call uzlset ('LABELYL',.FALSE.)
      call uyaxlb ('L',vplabel,5,vplabel,cplabel,4,5)
      call uyaxlb ('R',vplabel,5,vplabel,cplabel,4,5)
      call uziset ('IROTCYL',1)

!      call uysttl ('L','Pressure (hPa)',0.0)
      call sglset ('LCLIP',.FALSE.) ! Cliping
!      call sgtxzv (0.5*(vpr+vpl),0.020+vpt,trim(title1),0.018,0,0,5) !
!      call sgtxzv (vpr-0.02,0.035+vpt,trim(title2(1)),0.014,0,1,3) !
!      call sgtxzv (vpr-0.02,0.015+vpt,trim(title2(2)),0.014,0,1,3) !
end if

      call grcls

end do !!! ista

  call read_d1(vlons_area_d1,vlats_area_d1)

  psfile='map_'//trim(cvar)

  vpl=0.15
  vpr=0.85
  vpb=0.20
  vpt=0.70

  aratio=cos(3.1415926/180.0 *0.5*(range_latl+range_latr))*(range_lonr-range_lonl) / (range_latr-range_latl)

  vpb=max(vpt-(vpr-vpl)/aratio,0.10)
  if (vpb.eq.0.10)then
   vpl=0.55-0.5*(vpt-vpb)*aratio
   vpr=0.55+0.5*(vpt-vpb)*aratio
  end if

  call sgiset('IFONT',1)
  call swiset('ICLRMAP',iclrmap)
  call swcmll
  call swcset('FNAME',trim(psfile))
  call swlset('LSEP',.false.)
  call swiset('IFL',1) !!! PNG
  call swiset('IWIDTH',1000)
  call swiset('IHEIGHT',800)
  call gropn(iout)
  call sglset('LFULL',.true.)
  call sglset('LCLIP',.true.)
  call slmgn(0.0,0.0,0.0,0.0)
  call grfrm

  call grswnd(range_lonl,range_lonr,range_latl,range_latr)

  call grsvpt(vpl,vpr,vpb,vpt)
  call grstrn(11)

  call umlset('LGLOBE',.false.)
  call umiset('INDEXOUT',31)
!  call umscnt (0.5*(vlonl+vlonr),0.5*(vlatl+vlatr),0.0)

  call umpfit

  call grstrf


  call sglset('LCLIP',.true.)
!  call umpglb
  call umplim
  call umpmap('coast_world')

  call uulinz(npts_d1,vlons_area_d1,vlats_area_d1,1,13)

do iprof=1,nprof
   write(cnum3,'(I3)') iprof
  call sgtxzu(vlons(iprof),vlats(iprof),trim(cnum3),0.008,0,0,21)
end do

open (31,file='location_prof_'//trim(cvar)//'.txt',form='formatted')
do iprof=1,nprof
  call stftrf(vlons(iprof),vlats(iprof),vx,vy)
!  call stfpr2(vx,vy,rx,ry)
!  call stftrf(rx,ry,wx,wy)
!  call swfint(wx,wy,iwx,iwy)
   write(cnum3,'(I3.3)') iprof
   write(31,'(A3,2F10.6)') trim(cnum3),vx,vy
end do
close(31)


  amtics=20.0
  astics=20.0
 

  bmtics=10.0
  bstics=10.0

  call sglset('LCLIP',.false.)
  call uzinit
  call uzlset('LOFFSET',.false.)
  call uziset('INDEXT2',3)
  call uziset('INDEXT1',3)
  call uziset('INNER',-1)
  call uzrset('RSIZEL1',0.016)
  call uzrset('RSIZEC1',0.016)
  call uzrset('RSIZET1',0.006)
  call uzrset('RSIZET2',0.003)
    
  call uxaxdv('B',astics,amtics)
  call uxaxdv('T',astics,amtics)
  call uxsttl('B','Lon',0.0)
  call uyaxdv('L',bstics,bmtics)
  call uyaxdv('R',bstics,bmtics)
  call uysttl('L','Lat',0.0)

  call uzlset('LABELYR',.false.)

  call sgtxzv (0.5*(vpr+vpl),vpt+0.03,trim(title1),0.025,0,0,5) !
  call sgtxzv (vpr-0.01,vpt+0.045,trim(title2(1)),0.016,0,1,3) !
  call sgtxzv (vpr-0.01,vpt+0.020,trim(title2(2)),0.016,0,1,3) !
  call sgtxzv (vpl+0.01,vpt+0.025,trim(title3),0.018,0,-1,3) !
  call grcls 


      stop
end program draw_obs_dep

!==================================================!
subroutine read_d1(varray_lon,varray_lat)

include 'area_d1.h' 

real(4)::varray_lon(npts),varray_lat(npts)

varray_lon = vlons_area_d1
varray_lat = vlats_area_d1

return
end subroutine read_d1
!==================================================!
