!==================================================!
module setup
  implicit real(a-h,o-z)
  character(100) :: cfile='./test_data/obsdep.dat'
  integer :: nobs
  real(4) :: wk(11)
  integer :: n, iunit


  integer,parameter :: nsta_max=1000
  real(4) :: vlon_sta(nsta_max),vlat_sta(nsta_max)

character*40 ::psfile='map_d1-4'

character*20 ::title1,title2(2),title3

integer,parameter::nptsh_d1=560
integer,parameter::npts_d1=nptsh_d1*2+1
real(4)::vlons_area_d1(npts_d1) 
real(4)::vlats_area_d1(npts_d1) 

integer,parameter::iclrmap=4

real(4),parameter::range_lonl=95.0
real(4),parameter::range_lonr=175.0
real(4),parameter::range_latl=10.0
real(4),parameter::range_latr=55.0

character*3::cnum3

integer::nsta

end module setup
!==================================================!
program main
use setup

  ista=0   

  vlon_sta=0.0
  vlat_sta=0.0

  ielm_ref=3074
  open(21,file='stations.dat',form='formatted')

  open (iunit, file=cfile, form='unformatted', access='sequential', convert='big_endian')
  do n = 1, 200000
    read(iunit,iostat=ios) wk
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
      ielm.eq.ielm_ref ) then
    do itest=1,ista
      if (vlon_sta(itest).eq.vlon.and.vlat_sta(itest).eq.vlat) goto 999
    end do 
    ista=ista+1
    vlon_sta(ista)=vlon
    vlat_sta(ista)=vlat
    write(21,*) ista,vlon,vlat
  end if

  999 continue
  end do
  close (iunit)

  close(21)
  
  nsta=ista


call read_d1(vlons_area_d1,vlats_area_d1)

psfile='map_station'
title1='ADPUPA'
call draw

stop
end program main

!==================================================!
subroutine read_d1(varray_lon,varray_lat)

include 'area_d1.h' 

real(4)::varray_lon(npts),varray_lat(npts)

varray_lon = vlons_area_d1
varray_lat = vlats_area_d1

return
end subroutine read_d1
!==================================================!
subroutine draw 
use setup
  
  iout=2

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

do ista=1,nsta
   write(cnum3,'(I3)') ista
  call sgtxzu(vlon_sta(ista),vlat_sta(ista),trim(cnum3),0.010,0,0,21)
end do

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


return
end subroutine draw

!==================================================!
subroutine dcbar(vpxr,vpyl,dylen)
use setup

integer,parameter::ntpmax=40
integer::itpats(ntpmax)
real(4)::vtpats(ntpmax+1)

real(4)::xm4(4),ym4(4)
real(4)::xm5(5),ym5(5)
real(4)::xm7(7),ym7(7)
character*4::cvar
character*10::cfact


call ueqntl(ntpat)

do itpat=1,ntpat
 call ueqtlv(val1,val2,ipat,itpat)
 itpats(itpat)=ipat 
 vtpats(itpat)=val1
 vtpats(itpat+1)=val2
end do

dyp=dylen/real(ntpat)

call sglset ('LCLIP',.FALSE.) ! Cliping


!!write(*,*) ntpat

!!stop
!!! Color bar

!      do ic=1,ntpat
      do ic=2,ntpat-1
         xm4(1)=vpxr+0.01
         xm4(2)=xm4(1)+0.02
         xm4(3)=xm4(2)
         xm4(4)=xm4(1)
         ym4(1)=vpyl+real(ic-1)*dyp
         ym4(2)=ym4(1)
         ym4(3)=vpyl+real(ic)*dyp
         ym4(4)=ym4(3)
         itpat=itpats(ic)
         call sgtnzv(4,xm4,ym4,itpat)
      end do
      ic=1
         xm4(3)=vpxr+0.01
         xm4(2)=xm4(3)+0.02
         xm4(1)=xm4(3)+0.01
         ym4(3)=vpyl+real(ic)*dyp
         ym4(2)=ym4(3)
         ym4(1)=vpyl+real(ic-1)*dyp
         call sgtnzv(3,xm4(1:3),ym4(1:3),itpats(ic))
      ic=ntpat
         xm4(1)=vpxr+0.01
         xm4(2)=xm4(1)+0.02
         xm4(3)=xm4(1)+0.01
         ym4(1)=vpyl+real(ic-1)*dyp
         ym4(2)=ym4(1)
         ym4(3)=vpyl+real(ic)*dyp
         call sgtnzv(3,xm4(1:3),ym4(1:3),itpats(ic))


!!! Waku
!      xm5=(/vpxr+0.01,vpxr+0.03,vpxr+0.03,vpxr+0.01,vpxr+0.01/)
!      ym5=(/vpyl,vpyl,vpyl+dylen,vpyl+dylen,vpyl/)
!      call sgplzv(5,xm5,ym5,1,1) 
      xm7=(/vpxr+0.01,vpxr+0.02,vpxr+0.03,vpxr+0.03,vpxr+0.02,vpxr+0.01,vpxr+0.01/)
      ym7=(/vpyl+dyp,vpyl,vpyl+dyp,vpyl+dylen-dyp,vpyl+dylen,vpyl+dylen-dyp,vpyl+dyp/)
      call sgplzv(7,xm7,ym7,1,1) 


       do ic=2,ntpat,2
!         if (vtpats(ic).eq.real(int(vtpats(ic)))) then
           write(cvar,'(F4.1)') vtpats(ic) 
           call sgtxzv(xm4(2),vpyl+real(ic-1)*dyp,cvar,0.012,0,-1,3) ! 
!         end if
       end do


!         ixfac=-nint(log(factor)/log(10.0))+6
!         write(cfact,'(A,I2,A)') '*10|',ixfac,'"'
         write(cfact,'(A)') '(hPa)'
         call sgtxzv(xm4(1),vpyl+real(ntpat)*dyp+0.02,cfact,0.015,0,-1,3) ! 


return
end subroutine dcbar

!==============================================================!

