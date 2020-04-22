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

  integer,parameter::nprofmax=5000
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

  integer,parameter::nsmpmax=50000
  real(4)::vsmp_obs(nsmpmax)
  real(4)::vsmp_fm(nsmpmax)
  real(4)::vsmp_am(nsmpmax)
  real(4)::vsmp_fer(nsmpmax)
 
  integer,parameter::nmem=50
  real(4) :: wk_ext(nmem)

  character*2::cvar 
  character*40::psbase='scat' 
  character*40::psfile 
  character*20::xtitle 
  character*40::title1 
  character*40::title2(2) 
  character*40::title3 

  integer,parameter :: nsta_max=1000
  real(4) :: vlon_stat(nsta_max),vlat_stat(nsta_max)

  character*14::cread
  character*10::ctimestart
  character*10::ctimeend
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

  character*3::cnum3

  integer,parameter::nbin=20
  real(4)::vref_bin(nbin)
  real(4)::vmean_gues_bin(nbin)
  real(4)::vmean_anal_bin(nbin)
  logical::lmask_bin(nsmpmax)
 
  nargs=iargc()

  if(nargs.lt.3)then
    write(*,*) 'usage:: ./draw YYYYMMDDHHstart  YYYYMMDDHHend  [u/v/t/tv/q]'
    stop
  end if

  call getarg(1,cread)
  ctimestart=cread(1:10)
  call getarg(2,cread)
  ctimeend=cread(1:10)
  call getarg(3,cread)
  cvar=trim(cread)

title1=ctimestart(5:6)//'/'//ctimestart(7:8)//' '//ctimestart(9:10)//'Z'
if (ctime(5:5).eq.'0') title1=ctimestart(6:6)//'/'//ctimestart(7:8)//' '//ctimestart(9:10)//'Z'

if (ctimeend(5:5).eq.'0') then
 title1=trim(title1)//'-'//ctimeend(6:6)//'/'//ctimeend(7:8)//' '//ctimeend(9:10)//'Z'
else
 title1=trim(title1)//'-'//ctimeend(5:6)//'/'//ctimeend(7:8)//' '//ctimeend(9:10)//'Z'
end if

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
  vmin=260.0
  vmax=290.0
  amtics=5.0
  vfact=1.0
  xtitle='Temperature (K)'
  case('tv','Tv','TV')
  ielm_ref=3074
  vmin=275.0
  vmax=305.0
  amtics=5.0
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
  amtics=5.0
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

 ctime=ctimestart
 iprof=0

 do while ( trim(ctime) .ne. trim(ctimeend))
 
 ilev=0
 
 cfile=trim(cdirbase)//ctime//'0000/obs/obsdep.dat'

Write(*,*) cfile
  iunit=92
  open (iunit, file=trim(cfile), form='unformatted', access='sequential', convert='big_endian')
  do n = 1, nobsmax
     read(iunit,iostat=ios) wk
!     read(iunit,iostat=ios) wk_ext
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
      if (vlon.eq.vlons(iprof).and.vlat.eq.vlats(iprof))then
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
!   write(*,'(I5,I5,3F9.2,2E10.2,I4,F9.2,I4,2E10.2)') &
!         n,ielm,vlon,vlat,vlev,vdat,verr,ityp,vdif,iqc,vomb,voma
    vprof_lev(ilev,iprof)=vlev
    vprof_obs(ilev,iprof)=vdat
    vprof_oer(ilev,iprof)=verr
    vprof_fm(ilev,iprof)=vdat-vomb
!    vprof_fer(ilev,iprof)=sqrt(sum(wk_ext**2)/real(nmem-1))
    vprof_fer(ilev,iprof)=0.0
    vprof_am(ilev,iprof)=vdat-voma
    iprof_nlev(iprof)=ilev
    if (ilev.eq.1)then
      vlons(iprof)=vlon
      vlats(iprof)=vlat
    end if
  end if

  end do
  close (iunit)

    select case (ctime(9:10))
      case('00') ; ctime(9:10)='06'
      case('06') ; ctime(9:10)='12'
      case('12') ; ctime(9:10)='18'
      case('18') 
        read(ctime(1:4),*) iyear
        read(ctime(5:6),*) imon
        read(ctime(7:8),*) iday
        call datef3(1,iyear,imon,iday,iyear,imon,iday)
        write(ctime(1:8),'(I4.4,I2.2,I2.2)') iyear,imon,iday
        ctime(9:10)='00'
    end select

  end do !!! ctime

  nprof=iprof

  write(*,*) '# of profs : ',nprof

  vprof_obs=vprof_obs*vfact
  vprof_fm=vprof_fm*vfact
  vprof_am=vprof_am*vfact
  vprof_oer=vprof_oer*vfact
  vprof_fer=vprof_fer*vfact


  vlim_lonl=120.0
  vlim_lonr=150.0
  vlim_latl=20.0
  vlim_latr=60.0
  vlim_levb=1000.0
  vlim_levt=850.0


  ismp=0
  do iprof=1,nprof
   if(vlons(iprof).ge.vlim_lonl .and. vlons(iprof).le.vlim_lonr .and. &
      vlats(iprof).ge.vlim_latl .and. vlats(iprof).le.vlim_latr ) then
   do ilev=1,iprof_nlev(iprof)
     if (vprof_lev(ilev,iprof).ge. vlim_levt .and.vprof_lev(ilev,iprof).le. vlim_levb ) then
       ismp=ismp+1
       vsmp_obs(ismp)=vprof_obs(ilev,iprof)
       vsmp_fm(ismp)=vprof_fm(ilev,iprof)
       vsmp_fer(ismp)=vprof_fer(ilev,iprof)
       vsmp_am(ismp)=vprof_am(ilev,iprof)
     end if
   end do
   end if
  end do
nsmp=ismp

! *** quicklook ***

pbot=1050.0
ptop=100.0

astics=amtics

dvmin=-10.0
dvmax=10.0
bmtics=5.0
bstics=5.0

iout=2


do ibin=1,nbin
 vmin_bin = vmin+real(ibin-1)*(vmax-vmin)/real(nbin)
 vmax_bin = vmin+real(ibin)  *(vmax-vmin)/real(nbin)
 vref_bin(ibin) = 0.5*(vmin_bin+vmax_bin)

! lmask_bin(1:nsmp)= (vsmp_obs(1:nsmp)).gt.vmin_bin .and.&
!                    (vsmp_obs(1:nsmp)).le.vmax_bin  
 lmask_bin(1:nsmp)= (vsmp_fm(1:nsmp)).gt.vmin_bin .and.&
                    (vsmp_fm(1:nsmp)).le.vmax_bin  


 if (count(lmask_bin(1:nsmp)).ne.0)then
 vmean_gues_bin(ibin)=sum( vsmp_obs(1:nsmp)-vsmp_fm(1:nsmp), lmask_bin(1:nsmp)) /  &
                      real(count(lmask_bin(1:nsmp)))
 else
 vmean_gues_bin(ibin)=0.0
 end if 

! lmask_bin(1:nsmp)= (vsmp_obs(1:nsmp)).gt.vmin_bin .and.&
!                    (vsmp_obs(1:nsmp)).le.vmax_bin  
 lmask_bin(1:nsmp)= (vsmp_fm(1:nsmp)).gt.vmin_bin .and.&
                    (vsmp_fm(1:nsmp)).le.vmax_bin  


 if (count(lmask_bin(1:nsmp)).ne.0)then
 vmean_anal_bin(ibin)=sum( vsmp_obs(1:nsmp)-vsmp_am(1:nsmp), lmask_bin(1:nsmp)) /  &
                      real(count(lmask_bin(1:nsmp)))
 else
 vmean_anal_bin(ibin)=0.0
 end if 
 write(*,*) ibin,count(lmask_bin(1:nsmp)),vmean_anal_bin(ibin)
end do 

     write(psfile,'(4A,I3.3)')trim(psbase),'_',trim(cvar)
     write(title2(1),'(A,I6)') 'ndata: ',nsmp
!
! *** general settings ***
      call gliset ('MSGLEV',1)
      call sgiset ('IFONT',1)
      call swiset('ICLRMAP',4)
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
      call grswnd (vmin,vmax,dvmin,dvmax) ! set window

      vpl=0.27; vpr=0.77; vpb=0.23 ; vpt=0.73
      call grsvpt (vpl,vpr,vpb,vpt) ! set viewport

      call grstrn (1) ! linear or log
      call grstrf

! **** Lines & markers ****
      call sglset ('LCLIP',.TRUE.) ! Cliping

!     do ismp=1,nsmp
!        call uulinz(2,(/vsmp_obs(ismp),vsmp_obs(ismp)/),vsmp_obs(ismp)-(/vsmp_fm(ismp),vsmp_am(ismp)/),1,21 )
!        call uulinz(2,(/vsmp_obs(ismp),vsmp_obs(ismp)/),vsmp_obs(ismp)-vsmp_fm(ismp)+(/vsmp_fer(ismp),-vsmp_fer(ismp)/),1,41 )
!     end do

      nskip=(nsmp/8000)+1

!      call uumrkz(nsmp/nskip,vsmp_obs(1:nsmp:nskip),vsmp_obs(1:nsmp:nskip)-vsmp_fm(1:nsmp:nskip),10,41,0.004)
!      call uumrkz(nsmp/nskip,vsmp_obs(1:nsmp:nskip),vsmp_obs(1:nsmp:nskip)-vsmp_am(1:nsmp:nskip),10,21,0.004)
      call uumrkz(nsmp/nskip,vsmp_fm(1:nsmp:nskip),vsmp_obs(1:nsmp:nskip)-vsmp_fm(1:nsmp:nskip),10,41,0.004)
      call uumrkz(nsmp/nskip,vsmp_fm(1:nsmp:nskip),vsmp_obs(1:nsmp:nskip)-vsmp_am(1:nsmp:nskip),10,21,0.004)


      call uulinz (nbin,vref_bin,vmean_gues_bin,1,557)
      call uulinz (nbin,vref_bin,vmean_gues_bin,1,45)
      call uulinz (nbin,vref_bin,vmean_anal_bin,1,557)
      call uulinz (nbin,vref_bin,vmean_anal_bin,1,25)


!      call uulinz (2,(/vmin,vmax/),(/vmin,vmax/),3,1)
      call uulinz (2,(/vmin,vmax/),(/0.0,0.0/),3,1)


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
!      call uyaxdv ('R',astics,amtics)
!      call uyaxdv ('L',astics,amtics)
      call uyaxdv ('R',bstics,bmtics)
      call uyaxdv ('L',bstics,bmtics)
       call uziset ('IROTCYL',1)

      call uysttl ('L','o-b/o-a '//trim(xtitle),0.0)
      call sglset ('LCLIP',.FALSE.) ! Cliping
      call sgtxzv (0.5*(vpr+vpl),0.020+vpt,trim(title1),0.018,0,0,5) !
      call sgtxzv (vpr-0.02,0.04+vpt,trim(title2(1)),0.014,0,1,3) !
      call sgtxzv (vpr-0.02,0.02+vpt,trim(title2(2)),0.014,0,1,3) !



      call grcls


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

