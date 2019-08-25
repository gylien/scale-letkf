!==================================================!
module setup
implicit real(a-h,o-z)

include 'common_d4.h'

! input vars
integer::itimestart
integer::itimeend
integer::itimeskip
character*14::ctimeparent
character*14::ctimebase

!
integer::itimeref

integer,parameter::ioutsec=30
!
integer::izref

!integer,parameter::izref=13 !!! 1598m
!integer,parameter::izref=27 !!! 5171m

real(4)::val_mdet(nlon,nlat)

real(4)::val_plot_c(nlon,nlat)
real(4)::val_plot_s(nlon,nlat)


character*160 ::psfile='figure'
character*160 ::psfile_temp='figure'

character*40::cdummy
character*20 ::title1,title2(2),title3

real,parameter::vcmax=0.1
real,parameter::vcmin=-0.1
real,parameter::vcintv=0.02
real,parameter::vsmax_dif=0.1
real,parameter::vsmin_dif=0.0
real,parameter::vsintv_dif=0.01
real,parameter::vsmax_abs=55.0
real,parameter::vsmin_abs=5.0
real,parameter::vsintv_abs=5.0

real(4)::vsmax
real(4)::vsmin
real(4)::vsintv
integer::iclrmap

character*8::cftime='_fXXXXXX'

character*160::cfile
character*160::cplot
character*100::cfilegrd='dafcst/fcst_ref3d_YYYYMMDD-HHMMSS.grd'

character*8::cYMD
character*6::cHMS

character*160::command

end module setup
!==================================================!
program main
use setup

character*5::cnum5

narg = iargc()
if (narg.lt.6) then
 write(*,*) 'input $ctimebase $icycle $itimestart $itimeend $itimeskip $irefheight'
 stop
else
 call getarg (1,ctimebase)
 call getarg (2,cnum5)
 read(cnum5,*) icycle
 call getarg (3,cnum5)
 read(cnum5,*) itimestart
 call getarg (4,cnum5)
 read(cnum5,*) itimeend
 call getarg (5,cnum5)
 read(cnum5,*) itimeskip
 call getarg (6,cnum5)
 read(cnum5,*) irefheight
end if

izref = iblkge(axz,nz,real(irefheight))

 read(ctimebase(1:8),*) iYMD
 read(ctimebase(9:14),*) iHMS

 call time12(iHMS,iTT)
 iTT=iTT+icycle*isec_cycle 
 if (iTT.eq.86400) then
  iTT=0
  call datef1(1,iYMD,iYMD)
 end if
 call time21(iHMS,iTT)
 write(cYMD,'(I8.8)') iYMD
 write(cHMS,'(I6.6)') iHMS

write(cfilegrd(19:26),'(A8)') cYMD
write(cfilegrd(28:33),'(A6)') cHMS
cfile=trim(cdir_base_fcst)//'/'//ctimebase//'/'//trim(cfilegrd)

write(command,'(A)')'mkdir -p '//trim(cdir_base_plot)//'/'//cYMD//cHMS

call system(trim(command))

open(11,file=trim(cfile),form='unformatted',access='direct',recl=nlon*nlat)

do itimeref=itimestart,itimeend,itimeskip
write(cftime(2:8),'(I7)') 1000000 + ioutsec*(itimeref-1)
write(cftime(1:2),'(A2)')'_f'

 irec=(itimeref-1)*2*nz+izref
 read(11,rec=irec) val_mdet(:,:)

write(title3,'(I4,A)')irefheight,'m radar ref(dbz)'

val_plot_s=val_mdet 
val_plot_c=rmiss
where (abs(val_mdet).gt.1e10) val_plot_s=rmiss

 write(*,*) maxval(val_plot_s), minval(val_plot_s)

write(psfile,'(A,I4,A,A)') trim(cdir_base_plot)//'/'//cYMD//cHMS//'/dbz',irefheight,'m_mdet',cftime

iclrmap=12
title1='mdet'
vsmax=vsmax_abs
vsmin=vsmin_abs
vsintv=vsintv_abs
call draw
write(command,'(A)')'mv '//trim(psfile_temp)//'_0001.png '//trim(psfile)//'_0001.png '

call system(trim(command))


end do

close(11)

stop
end program main

!==================================================!
subroutine draw 
use setup

integer,parameter::nwork=3*(nlon+2)*(nlat+2)/2+1
integer::iwork(nwork)

integer,parameter::npatmax=40
real(4)::vtlevs(npatmax)
integer::itpats(npatmax)

real(4)::vmask(nlon_ext,nlat_ext)

  iout=2


do ilon=1,nlon_ext
do ilat=1,nlat_ext
 if (ilon-nlonadd/2.ge.1.and.ilon-nlonadd/2.lt.nlon .and. &
     ilat-nlatadd/2.ge.1.and.ilat-nlatadd/2.lt.nlat  )then
  if (val_plot_s(ilon-nlonadd/2,ilat-nlatadd/2).ne.rmiss)then
   vmask(ilon,ilat)=0.0
  else
   vmask(ilon,ilat)=rmiss
  end if
 else
  vmask(ilon,ilat)=rmiss
 end if
end do
end do

  vpl=0.15
  vpr=0.85
  vpb=0.20
  vpt=0.70

  aratio=(range_lonr-range_lonl) / (range_latr-range_latl)

  vpb=max(vpt-(vpr-vpl)/aratio,0.10)
  if (vpb.eq.0.10)then
   vpl=0.55-0.5*(vpt-vpb)*aratio
   vpr=0.55+0.5*(vpt-vpb)*aratio
  end if

  call sgiset('IFONT',1)
  call swiset('ICLRMAP',iclrmap)
  call swcmll
  call swcset('FNAME',trim(psfile_temp))
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
  call grstrn(10)

  call umlset('LGLOBE',.false.)
  call umiset('INDEXOUT',31)
!  call umscnt (0.5*(vlonl+vlonr),0.5*(vlatl+vlatr),0.0)

  call umpfit

  call grstrf


  call glrset ('RMISS', rmiss)
  call gllset ('LMISS', .true.)

  call uwsgxa (axlon,nlon)
  call uwsgya (axlatSN,nlat)

!  call uegtla (vsmin,vsmax,vsintv)
!  if (abs(vsmin).eq.abs(vsmax)) call uestlv(-vsintv,vsintv,55999)

  call ueitlv

!  ntpat = int((vsmax-vsmin) / vsintv)
!  vtlevs (1) = -1.0e6
!   do ipat=1,ntpat
!    vtlevs(ipat+1) = vsmin + real(ipat-1) * vsintv
!    vtlevs(ipat+2) = vsmin + real(ipat) * vsintv
!   end do
!  vtlevs (ntpat+3) = 1.0e6

!  if (abs(vsmin).eq.abs(vsmax)) then
!   itpats(1:ntpat+2) = (/ 16, 20,30,40,50, 55,55, 60,70,80,90, 94/) * 1000+999
!  else
!   itpats(1:ntpat+2) = (/ 10, 14,20,30,40, 50,60, 70,78,84,92, 96/) * 1000+999
!  end if


   ntpat = 9

   itpats(1:ntpat+2) = (/ 0, 40,34,50,62,68,74,80,84,92,98/) * 1000+999
   itpats(1) = 0
!   vtlevs(1:ntpat+3) = (/-1.0e6,0.5,1.0,2.0,3.0,5.0,7.0,10.0,15.0,20.0,30.0,40.0,60.0,80.0,1.0e6/)
   vtlevs(2:ntpat+2) = 5.0 + (/( 5.0*real(i), i=1,10 )/)
   vtlevs(1) = -1.0e10
   vtlevs(ntpat+3) = 1.0e10
  call uestln(vtlevs(1:ntpat+3),itpats(1:ntpat+2),ntpat+2)

  call uetone (val_plot_s,nlon,nlon,nlat)

  call dcbar(vpr+0.02,vpb,(vpt-vpb)*0.8)

  call udsfmt ('B')
  call udrset ('RSIZEL',0.012)
  call udlset ('LMSG',.false.)

!  call udgcla (vcmin,vcmax,vcintv)
!  call udcntr (val_plot_c,nlon,nlon,nlat)
!  call udcntz (val_plot_c,nlon,nlon,nlat,iwork,nwork)


!!! masking
  call ueitlv

  call glrset ('RMISS', 0.0)
  call gllset ('LMISS', .false.)
  call sglset('LCLIP',.true.)
  call uwsgxa (axlon_ext,nlon_ext)
  call uwsgya (axlatSN_ext,nlat_ext)
  call uestlv(rmiss-abs(rmiss)*0.01,rmiss+abs(rmiss)*0.01,1602)
  call uetone (vmask,nlon_ext,nlon_ext,nlat_ext)


  call umlset ('LGRIDMJ',.false.)
  call umrset ('DGRIDMN',0.5)
  call umiset ('ITYPEMN',3)
  call umiset ('INDEXMN',1)

  call umpglb
!  call umplim
  call umpmap('coast_japan')

!  call uulinz(npts,vlons_area_d2,vlats_area_d2,3,91)

  amtics=0.5
  astics=0.5

  bmtics=0.5
  bstics=0.5

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
  call uxsfmt ('(F5.1)')
  call uysfmt ('B')

    
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


!       do ic=2,ntpat,2
       do ic=2,ntpat,1
!         if (vtpats(ic).eq.real(int(vtpats(ic)))) then
          if (vtpats(ic).ge.1.0)then
           write(cvar,'(I3)') int(vtpats(ic)) 
          else
           write(cvar,'(F4.1)') vtpats(ic) 
          end if
           call sgtxzv(xm4(2),vpyl+real(ic-1)*dyp,cvar,0.012,0,-1,3) ! 
!         end if
       end do


!         ixfac=-nint(log(factor)/log(10.0))+6
!         write(cfact,'(A,I2,A)') '*10|',ixfac,'"'
         write(cfact,'(A)') '(dbz)'
         call sgtxzv(xm4(1),vpyl+real(ntpat)*dyp+0.02,cfact,0.015,0,-1,3) ! 


return
end subroutine dcbar

!==============================================================!
