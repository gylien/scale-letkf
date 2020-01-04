program main

implicit real(a-h,o-z)

!
! conventional observations
!
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_tv_obs=3074
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331!

  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_rain_obs=19999
  INTEGER,PARAMETER :: id_tclon_obs=99991  ! TC vital
  INTEGER,PARAMETER :: id_tclat_obs=99992  ! TC vital
  INTEGER,PARAMETER :: id_tcmip_obs=99993  ! TC vital
!
! radar observations
!
  INTEGER,PARAMETER :: id_radar_ref_obs=4001
  INTEGER,PARAMETER :: id_radar_ref_zero_obs=4004
  INTEGER,PARAMETER :: id_radar_vr_obs=4002
  INTEGER,PARAMETER :: id_radar_prh_obs=4003
!
! Himawari-8 (H08) observations
!
  INTEGER,PARAMETER :: id_H08IR_obs=8800

  INTEGER,PARAMETER :: nobtype=24
!
  CHARACTER(6),PARAMETER :: obtypelist(nobtype)= &
     (/'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR', &
       'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG', &
       'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND', &
       'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW', &
       'TMPAPR', 'PHARAD', 'H08IRB', 'TCVITL'/) ! H08


integer,parameter::nmax=100000
real(4)::vlons(nmax),vlats(nmax)
integer::ivars(nmax)
real(4)::wk(8)
character*120::cdir_in='/work/hp150019/share/SCALE-LETKF-rt/external/ncepobs_gdas_letkf/'
character*120::cfile_in

character*20 ::psfile='figure'

character*40::cdummy
character*20 ::title1,title2(2),title3
character*20 :: carg

call getarg(1,carg)
!
if (len(trim(carg)).ne.10)then
 write(*,*)'input YYYYMMDDHH'
 stop
end if

cfile_in = trim(carg)//'/'//'obs_'//trim(carg)//'0000.dat'
!cfile_in = '/data9/amemiya/realtime_OFP/data/ncepobs_gdas/run/dec_prepbufr/outdir_3/'//trim(carg)//'00.dat'
!cfile_in = '/data9/amemiya/realtime_OFP/data/ncepobs_gdas/run/dec_prepbufr/outdir_2/obs_20190319000000.dat'
!cfile_in = '/data9/amemiya/realtime_OFP/data/ncepobs_gdas_letkf/2019031900/obs_201903190000.dat'

!write(*,*) trim(cdir_in)//trim(cfile_in)
!stop

do iobtype=1,nobtype
 psfile=obtypelist(iobtype)
 title1=obtypelist(iobtype)

 irecs=0
 nsta=0
 vlons=-999.9
 vlats=-999.9
 ivars=0

 ios=0
 open(90,file=trim(cdir_in)//trim(cfile_in),form='unformatted',iostat=ios)
! open(90,file=trim(cfile_in),form='unformatted',iostat=ios)

  do while (ios.eq.0) 
  irecs=irecs+1
   read(90,iostat=ios) wk

!   write(*,*) ios, wk
!   stop

 if (int(wk(7)).eq.iobtype) then
   j=0
   do i=1,nsta
    if (int(vlons(i)*100.0).eq.int(wk(2)*100.0).and.int(vlats(i)*100.0).eq.int(wk(3)*100.0)) j=i
   end do

   if (j.ne.0)then
    ivars(j)=ivars(j)+1
   else
    nsta=nsta+1
    vlons(nsta)=wk(2)
    vlats(nsta)=wk(3)
    ivars(nsta)=ivars(nsta)+1
   end if   

 end if

  end do

close(90)

write(*,*) obtypelist(iobtype), nsta

if (nsta.ne.0)then
write(title2(1),'(I6,A)') nsta,' points'
write(title2(2),'(I6,A)') sum(ivars),'   data'
write(title3,'(A)') trim(carg)

  iout=2

  vlonl=90.0
  vlonr=180.0
  vlatl=0.0
  vlatr=60.0

  vpl=0.25
  vpr=0.85
  vpb=0.30
  vpt=0.70

  aratio=(vlonr-vlonl) / (vlatr-vlatl)

  vpb=max(vpt-(vpr-vpl)/aratio,0.10)
  if (vpb.eq.0.10)then
   vpl=0.55-0.5*(vpt-vpb)*aratio
   vpr=0.55+0.5*(vpt-vpb)*aratio
  end if

  call sgiset('IFONT',1)
  call sgiset('ICLRMAP',8)
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

  call grswnd(vlonl,vlonr,vlatl,vlatr)

  call grsvpt(vpl,vpr,vpb,vpt)
  call grstrn(10)

  call umlset('LGLOBE',.false.)
  call umiset('INDEXOUT',3)
!  call umscnt (0.5*(vlonl+vlonr),0.5*(vlatl+vlatr),0.0)

  call umpfit

  call grstrf

  call uumrkz(nsta, vlons(1:nsta) ,vlats(1:nsta),4,21,0.002)

  call sglset('LCLIP',.true.)
!  call umpglb
  call umplim
  call umpmap('coast_world')

  call uwsgxa((/vlonl,vlonr/),2)
  call uwsgya((/vlatl,vlatr/),2)

  amtics=10.0
  astics=10.0
 

  bmtics=10.0
  bstics=10.0

  call sglset('LCLIP',.false.)
  call uzinit
  call uzlset('LOFFSET',.false.)
  call uziset('INDEXT2',3)
  call uziset('INDEXT1',3)
  call uziset('INNER',-1)
  call uzrset('RSIZEC1',0.025)
  call uzrset('RSIZET1',0.010)
  call uzrset('RSIZET2',0.008)
    
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

end if
  
end do


 
stop
end program main
