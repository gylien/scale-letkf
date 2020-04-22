!==================================================!
module setup
implicit real(a-h,o-z)

! input vars
character*40 ::psfile='map_subdom'
character*40 ::figdir='../figs/map/'

character*40::cdir_data
character*40::cdummy
character*20 ::title1,title2(2),title3

integer,parameter::npts_d1_sub=81
integer,parameter::ndomx=16
integer,parameter::ndomy=12
integer,parameter::ndom=192

real(4)::vlons_area_d1_sub(npts_d1_sub,ndom) 
real(4)::vlats_area_d1_sub(npts_d1_sub,ndom) 
real(4)::val_sub(ndom) 
real(4)::rmse_sub(ndom) 
real(4)::bias_sub(ndom) 
integer::nobs_sub(ndom) 

real(4)::vlons_cen_sub(ndom) 
real(4)::vlats_cen_sub(ndom) 
character*3::cdom_sub(ndom) 


character*14::ctime

real(4),parameter::range_lonl=95.0
real(4),parameter::range_lonr=175.0
real(4),parameter::range_latl=10.0
real(4),parameter::range_latr=55.0

integer,parameter::nargmax=3

character(len=10)::cargs(nargmax)

character(len=2)::citem
character(len=4)::ctype
character(len=4)::cstep

end module setup
!==================================================!
program main
use setup


narg=iargc()
if(narg.lt.2)then
write(*,*) 'usage:: ./draw_map_subdom [u/v/t/q/ps] [rmse/bias/nobs] [anal/gues](optional)'
stop
end if

 call getarg(1,citem)
 call getarg(2,ctype)
cstep='anal'
if (narg .ge.3) call getarg(3,cstep)


select case(trim(cstep))
case('anal');title2(1)='Analysis'
case('gues');title2(1)='First guess'
end select

call read_d1_sub(vlons_area_d1_sub,vlats_area_d1_sub)

cdir_data='../data/'

call readdata

psfile=trim(figdir)//trim(citem)//'_'//trim(cstep)//'_'//trim(ctype)//'_map_d1_sub'

title2(2)=ctime(1:4)//'/'//ctime(5:6)//'/'//ctime(7:8)//' '//ctime(9:10)//'Z'

select case(trim(citem))
case('u'); title3='U (m/s)'
case('v'); title3='V (m/s)'
case('t'); title3='Temperature (K)'
case('q'); title3='Water vapor (g/kg)'
case('ps'); title3='Surf. Pres. (hPa)'
end select

select case(trim(ctype))
case('rmse')
    title1='RMSE'
    val_sub=rmse_sub
case('bias')
    title1='BIAS'
    val_sub=bias_sub
case('nobs')
    title1='Number of obs'
    val_sub=real(nobs_sub)
end select

call draw

stop
end program main

!==================================================!
subroutine readdata
use setup
character*50::cfile

do idom=1,ndom
 write(cfile,'(A,I3.3,5A)')trim(cdir_data),idom-1,'/',trim(citem),'_',trim(cstep),'.txt'
open(21,file=trim(cfile),form='formatted')
  ios=0
  do while (ios.eq.0) 
!  ctime=''
!  do while (trim(ctime).ne.'20200414120000') 
   read(21,*,iostat=ios) cdummy,bias,rmse,nobs  
!   if (ios.ne.0)then 
     rmse_sub(idom) = rmse
     bias_sub(idom) = bias
     nobs_sub(idom) = nobs
     ctime=cdummy
!     write(*,*) idom,ctime, nobs
!     if (ctime.eq.'20200414180000') stop
!   end if
  end do
close(21)
end do

!write(*,*) maxval(nobs_sub)
!write(*,*) maxval(rmse_sub)
!stop


return
end subroutine readdata
!==================================================!
subroutine read_d1_sub(varray_lon,varray_lat)
use setup 
character*40::cdatfile='area_d1_subdom_192.dat'

open(31,file=trim(cdatfile),form='unformatted')
  read(31) vlons_area_d1_sub
  read(31) vlats_area_d1_sub
close(31)


!write(*,*) maxval(vlons_area_d1_sub(:,1))
!write(*,*) maxval(vlats_area_d1_sub(:,1))


do idom = 1,ndom
  vlons_cen_sub(idom)=0.25 *( vlons_area_d1_sub(1, idom) &
                           +  vlons_area_d1_sub(21,idom) &
                           +  vlons_area_d1_sub(41,idom) &
                           +  vlons_area_d1_sub(61,idom) )
  vlats_cen_sub(idom)=0.25 *( vlats_area_d1_sub(1, idom) &
                           +  vlats_area_d1_sub(21,idom) &
                           +  vlats_area_d1_sub(41,idom) &
                           +  vlats_area_d1_sub(61,idom) )
  write(cdom_sub(idom),'(I3.3)') idom-1
end do

return
end subroutine read_d1_sub
!==================================================!

subroutine draw 
use setup

integer,parameter::npatmax=40
real(4)::vtlevs(npatmax)
integer::itpats(npatmax)


ntpat=10

select case(trim(ctype))
case('rmse')
    iclrmap=6
    itpats(1:ntpat+2) = (/ 0, 40,34,30,50,62,68,74,80,84,92,98/) * 1000 + 999
    itpats(1) = 0
    vtlevs(1) = -1.0e10
    vtlevs(ntpat+3) = 1.0e10
    select case(trim(citem))
    case('ps')
      vtlevs(2:ntpat+2) = 0.2 + (/( 0.2*real(i), i=1,11 )/)
      val_sub = val_sub * 0.01
    case('q')
      vtlevs(2:ntpat+2) = 0.15 + (/( 0.15*real(i), i=1,11 )/)
      val_sub = val_sub * 1000.0
    case('u')
      vtlevs(2:ntpat+2) = 0.5 + (/( 0.5*real(i), i=1,11 )/)
    case('v')
      vtlevs(2:ntpat+2) = 0.5 + (/( 0.5*real(i), i=1,11 )/)
    case('t')
      vtlevs(2:ntpat+2) = 0.2 + (/( 0.2*real(i), i=1,11 )/)
     end select 
case('bias')
    iclrmap=4
    itpats(1:ntpat+2) = (/ 0, 40,34,30,50,62,68,74,80,84,92,98/) * 1000 + 999
    itpats(1) = 0
    vtlevs(1) = -1.0e10
    vtlevs(ntpat+3) = 1.0e10
    select case(trim(citem))
    case('ps')
    vtlevs(2:ntpat+2) = 0.15 + (/( 0.15*real(i), i=-6,6 )/)
    val_sub = val_sub * 0.01
    case('q')
    vtlevs(2:ntpat+2) = 0.05 + (/( 0.05*real(i), i=-6,6 )/)
    val_sub = val_sub * 1000.0
    case('u')
    vtlevs(2:ntpat+2) = 0.5 + (/( 0.5*real(i), i=-6,6 )/)
    case('v')
    vtlevs(2:ntpat+2) = 0.5 + (/( 0.5*real(i), i=-6,6 )/)
    case('t')
    vtlevs(2:ntpat+2) = 0.2 + (/( 0.2*real(i), i=-6,6 )/)
    end select 
case('nobs')
    iclrmap=6
    itpats(1:ntpat+2) = (/ 0, 40,34,30,50,62,68,74,80,84,92,98/) * 1000 + 999
    vtlevs(2:ntpat+2) = 1.0e-4+real( (/ 0,1,2,5,10,20,50,100,200,500,1000 /) )
    itpats(1) = 1615
    vtlevs(1) = -1.0e10
    vtlevs(ntpat+3) = 1.0e10
end select


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
  call umiset('INDEXOUT',11)
!  call umscnt (0.5*(vlonl+vlonr),0.5*(vlatl+vlatr),0.0)

  call umpfit

  call grstrf


!  call uulinz(npts_d1,vlons_area_d1,vlats_area_d1,1,13)

do idom=1,ndom
  itpat=iblkge( vtlevs(1:ntpat+3), ntpat+3, val_sub(idom) )
  itone=itpats(itpat)
  if (nobs_sub(idom).le.1) itone=1615
  call sgtnzu(npts_d1_sub,vlons_area_d1_sub(:,idom),vlats_area_d1_sub(:,idom),itone)
  call uulinz(npts_d1_sub,vlons_area_d1_sub(:,idom),vlats_area_d1_sub(:,idom),1,41)
end do

open (31,file='location.txt',form='formatted')
do idom=1,ndom
  call stftrf(vlons_cen_sub(idom),vlats_cen_sub(idom),vx,vy)
  call stfpr2(vx,vy,rx,ry)
  call stftrf(rx,ry,wx,wy)
  call swfint(wx,wy,iwx,iwy)
  write(31,'(A3,2F10.6)') cdom_sub(idom),iwx,iwy
end do
close(31)

  call uestln(vtlevs(1:ntpat+3),itpats(1:ntpat+2),ntpat+2)
  call dcbar(vpr+0.02,vpb,(vpt-vpb)*0.8)
  call sglset('LCLIP',.true.)
!  call umpglb
  call umplim
  call umpmap('coast_world')


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

intv=2
if (ctype.eq.'nobs')intv=1
       do ic=2,ntpat,2
!         if (vtpats(ic).eq.real(int(vtpats(ic)))) then
           write(cvar,'(F4.1)') vtpats(ic) 
           if (ctype.eq.'nobs') write(cvar,'(I4)')int(vtpats(ic))
           call sgtxzv(xm4(2),vpyl+real(ic-1)*dyp,cvar,0.012,0,-1,3) ! 
!         end if
       end do


!         ixfac=-nint(log(factor)/log(10.0))+6
!         write(cfact,'(A,I2,A)') '*10|',ixfac,'"'
!         write(cfact,'(A)') '(hPa)'
!         call sgtxzv(xm4(1),vpyl+real(ntpat)*dyp+0.02,cfact,0.015,0,-1,3) ! 


return
end subroutine dcbar

!==============================================================!
