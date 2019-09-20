!==================================================!
subroutine plot_dbz_DCL(val_plot_s,land2dgs,psfile,cheight,csec) 
  use common
  use common_scale
  use scale_io, only: &
      H_LONG
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
  use scale_atmos_grid_cartesC, only: &
      DX, DY, &
      GRID_CXG => ATMOS_GRID_CARTESC_CXG, &
      GRID_CYG => ATMOS_GRID_CARTESC_CYG
  use scale_time, only: &
      TIME_NOWDATE

  implicit none

  real(r_sngl),intent(in) :: val_plot_s(nlong,nlatg)
  real(r_sngl),intent(in) :: land2dgs(nlong,nlatg)
  character(*),intent(in) :: psfile
  character(len=5),intent(in) :: cheight
  character(len=4),intent(in) :: csec

  real(r_sngl),allocatable :: val_plot(:,:)

  integer,allocatable :: iwork(:)
  
  integer :: nlonadd 
  integer :: nlatadd 
  real(r_sngl),parameter :: OFFXY = 0.0e3 ! offset (m)
  real(r_sngl),allocatable :: grid_cxg_ext(:)
  real(r_sngl),allocatable :: grid_cyg_ext(:)
  real(r_sngl),allocatable :: vmask(:,:)
  
  real(r_sngl),parameter :: rmiss = -9.99e20
  
  integer,parameter :: npatmax = 40
  real(r_sngl) :: vtlevs(npatmax)
  integer :: itpats(npatmax)  
  character(len=40) :: title1, title2(2), title3

  integer :: iclrmap
  integer :: ilon, ilat, nwork, nlong_ext, nlatg_ext
  real(r_sngl) :: vpr, vpl, vpt, vpb

  real(r_sngl) :: range_lonl, range_lonr
  real(r_sngl) :: range_latl, range_latr
  real(r_sngl) :: aratio
  integer :: i, ntpat
  real(r_sngl) :: amtics, astics
  real(r_sngl) :: bmtics, bstics

  character(len=19) :: ftimelabel



  write(*,*) 'max min val_plot', maxval(val_plot_s),minval(val_plot_s)



  nlonadd = int(OFFXY / DX)
  nlatadd = int(OFFXY / DY)

  write(ftimelabel,'(I4.4,A1,I2.2,A1,I2.2,A1, I2.2,A1,I2.2,A1,I2.2)')&
       TIME_NOWDATE(1), '/', TIME_NOWDATE(2), '/' ,TIME_NOWDATE(3), ' ', &
       TIME_NOWDATE(4), ':', TIME_NOWDATE(5), ':', TIME_NOWDATE(6)

  nlong_ext = nlong + nlonadd * 2
  nlatg_ext = nlatg + nlatadd * 2 

  allocate(val_plot(nlong,nlatg))
  allocate(grid_cxg_ext(nlong_ext))
  allocate(grid_cyg_ext(nlatg_ext))
  allocate(vmask(nlong_ext,nlatg_ext))

  nwork = 3 * (nlong + 2) * (nlatg + 2) / 2 + 1
  allocate(iwork(nwork))

  val_plot = val_plot_s
  where(.not.val_plot > rmiss) val_plot = rmiss

  title1 = trim(ftimelabel) // " UTC (FT=" // csec //"s)"
  if (trim(csec).eq.'anal')  title1 = trim(ftimelabel) // " UTC (Analysis)"
  title2 = (/'',''/)
  title3 = 'radar ref ' // cheight // ' m'
  
  if (cheight(1:1) == '0') title3='Radar reflectivity z=' // cheight(2:5) // ' m'
  if (cheight(1:2) == '00') title3='Radar reflectivity z=' // cheight(3:5) // ' m'


  iclrmap = 12
  
  do ilon = 1, nlong_ext
  do ilat = 1, nlatg_ext
    if (ilon - nlonadd >= 1 .and. ilon - nlonadd < nlong .and. &
       ilat - nlatadd >= 1 .and. ilat - nlatadd < nlatg  )then
      if (val_plot_s(ilon-nlonadd,ilat-nlatadd) /= rmiss)then
        vmask(ilon,ilat) = 0.0
      else
        vmask(ilon,ilat) = rmiss
      end if
    else
      vmask(ilon,ilat) = rmiss
    end if
  end do
  end do

  grid_cxg_ext(nlonadd+1:nlong_ext-nlonadd) = real(grid_cxg(IHALO+1:IHALO+nlong))
  grid_cyg_ext(nlatadd+1:nlatg_ext-nlatadd) = real(grid_cyg(JHALO+1:JHALO+nlatg))

  do ilon = 1, nlonadd
    grid_cxg_ext(ilon) = grid_cxg_ext(nlonadd+1) &
                       - real(nlonadd+1-ilon) * (grid_cxg_ext(nlonadd+2) &
                       - grid_cxg_ext(nlonadd+1)) 
    grid_cxg_ext(nlong_ext-ilon+1) = grid_cxg_ext(nlong_ext-nlonadd) &
                                   + real(nlonadd+1-ilon) * (grid_cxg_ext(nlonadd+2) &
                                   - grid_cxg_ext(nlonadd+1)) 
  end do
  do ilat = 1, nlatadd
    grid_cyg_ext(ilat) = grid_cyg_ext(nlatadd+1) &
                       - real(nlatadd+1-ilat) * (grid_cyg_ext(nlatadd+2) &
                       - grid_cyg_ext(nlatadd+1)) 
    grid_cyg_ext(nlatg_ext-ilat+1) = grid_cyg_ext(nlatg_ext-nlatadd) &
                                   + real(nlatadd+1-ilat) * (grid_cyg_ext(nlatadd+2) &
                                   - grid_cyg_ext(nlatadd+1)) 
  end do

  range_lonl = GRID_CXG(IHALO+1) - OFFXY ! DX(250m) x 20 grids
  range_lonr = GRID_CXG(IHALO+nlong) + OFFXY
  range_latl = GRID_CYG(JHALO+1) - OFFXY
  range_latr = GRID_CYG(JHALO+nlatg) + OFFXY

  vpl = 0.15
  vpr = 0.85
  vpb = 0.20
  vpt = 0.70

  aratio = (range_lonr - range_lonl) / (range_latr - range_latl)

  vpb = max(vpt - (vpr - vpl) / aratio, 0.10)
  if (vpb == 0.10)then
    vpl = 0.55 - 0.5 * (vpt - vpb) * aratio
    vpr = 0.55 + 0.5 * (vpt - vpb) * aratio
  end if

  call gliset('MSGLEV',1)
  call sgiset('IFONT',1)
  call swiset('ICLRMAP',iclrmap)
  call swcmll
  call swcset('FNAME',trim(psfile))
  call swlset('LSEP',.false.)
  call swiset('IFL',1) !!! PNG
  call swiset('IWIDTH',1000)
  call swiset('IHEIGHT',800)
  call gropn(2)
  call sglset('LFULL',.true.)
  call sglset('LCLIP',.true.)
  call slmgn(0.0,0.0,0.0,0.0)
  call grfrm


  call grswnd(range_lonl,range_lonr,range_latl,range_latr)

  call grsvpt(vpl,vpr,vpb,vpt)
  call grstrn(1) !!! Give up map proj

!!!!!  call grstrn(10) !!! Mercator
!!!!  call umlset('LGLOBE',.false.)
!!!!  call umiset('INDEXOUT',31)
!  call umscnt (0.5*(vlonl+vlonr),0.5*(vlatl+vlatr),0.0)
!!!!  call umpfit

  call grstrf

  call glrset ('RMISS', rmiss)
  call gllset ('LMISS', .true.)

  call uwsgxa (grid_cxg_ext(nlonadd+1:nlong_ext-nlonadd),nlong)
  call uwsgya (grid_cyg_ext(nlatadd+1:nlatg_ext-nlatadd),nlatg)

  call ueitlv

!!   vtlevs(1:ntpat+3) = (/-1.0e6,0.5,1.0,2.0,3.0,5.0,7.0,10.0,15.0,20.0,30.0,40.0,60.0,80.0,1.0e6/) !!! rain 
!  ntpat = 9
!  itpats(1:ntpat+2) = (/ 0, 40,34,50,62,68,74,80,84,92,98/) * 1000 + 999
!  vtlevs(2:ntpat+2) = 5.0 + (/( 5.0*real(i), i=1,10 )/)
  ntpat = 10
  itpats(1:ntpat+2) = (/ 0, 40,34,30,50,62,68,74,80,84,92,98/) * 1000 + 999
  vtlevs(2:ntpat+2) = 5.0 + (/( 5.0*real(i), i=1,11 )/)
  itpats(1) = 0
  vtlevs(1) = -1.0e10
  vtlevs(ntpat+3) = 1.0e10
  call uestln(vtlevs(1:ntpat+3),itpats(1:ntpat+2),ntpat+2)

  call uetone (val_plot,nlong,nlong,nlatg)

  call dcbar(vpr+0.02,vpb,(vpt-vpb)*0.8)

  call udlset ('LMSG',.false.)
  call udlset ('LABEL',.false.)

  call udiclv 
  !call udsclv(1.0,31,1,'',-1.0) !!! z > 1.0m  -- approximate coastline 
  call udsclv(0.5,11,1,'',-1.0) !!! z > 1.0m  -- approximate coastline 


  call udcntz (land2dgs,nlong,nlong,nlatg,iwork,nwork)

!!! masking
  call ueitlv

  call glrset ('RMISS', 0.0)
  call gllset ('LMISS', .false.)
  call sglset('LCLIP',.true.)
  call uwsgxa (grid_cxg_ext,nlong_ext)
  call uwsgya (grid_cyg_ext,nlatg_ext)
  call uestlv(rmiss-abs(rmiss)*0.01,rmiss+abs(rmiss)*0.01,1602)
  call uetone (vmask,nlong_ext,nlong_ext,nlatg_ext)

!!! map
 
 ! call umlset ('LGRIDMJ',.false.)
!  call umrset ('DGRIDMN',0.5)
!  call umiset ('ITYPEMN',3)
!  call umiset ('INDEXMN',1)
!
!  call umpglb
!  call umplim
!  call umpmap('coast_japan')

!  call uulinz(npts,vlons_area_d2,vlats_area_d2,3,91)


  
  call uumrkz(1,0.5*(range_lonr+range_lonl),0.5*(range_latr+range_latl),9,21,0.010)


  amtics = 50.0 !! km
  astics = 50.0
  bmtics = 50.0
  bstics = 50.0

  call sglset('LCLIP',.false.)
  call uzinit
  call uzlset('LOFFSET',.true.)
  call uzrset('XOFFSET',-0.001*0.5*(range_lonr+range_lonl))
  call uzrset('YOFFSET',-0.001*0.5*(range_lonr+range_latl))
  call uzrset('XFACT',0.001)
  call uzrset('YFACT',0.001)
  call uziset('INDEXT2',3)
  call uziset('INDEXT1',3)
  call uziset('INNER',-1)
  call uzrset('RSIZEL1',0.016)
  call uzrset('RSIZEC1',0.016)
  call uzrset('RSIZET1',0.006)
  call uzrset('RSIZET2',0.003)
  call uxsfmt ('(I3)')
  call uysfmt ('(I3)')

    
  call uxaxdv('B',astics,amtics)
  call uxaxdv('T',astics,amtics)
  call uxsttl('B','X (km)',0.0)
  call uyaxdv('L',bstics,bmtics)
  call uyaxdv('R',bstics,bmtics)
  call uysttl('L','Y (km)',0.0)

  call uzlset('LABELYR',.false.)

  call sglset('LCLIP',.false.)
 


  call sgtxzv (0.5*(vpr+vpl),vpt+0.02,trim(title1),0.016,0,0,3) !
  call sgtxzv (vpr-0.01,vpt+0.045,trim(title2(1)),0.016,0,1,3) !
  call sgtxzv (vpr-0.01,vpt+0.020,trim(title2(2)),0.016,0,1,3) !
  call sgtxzv (0.5*(vpr+vpl),vpt+0.05,trim(title3),0.017,0,0,4) !

  call grcls 

  deallocate(val_plot,grid_cxg_ext,grid_cyg_ext,vmask,iwork)

return
end subroutine plot_dbz_DCL

!==================================================!
subroutine dcbar(vpxr,vpyl,dylen)
  use common
  implicit none

  integer,parameter :: ntpmax = 40
  integer :: itpats(ntpmax)
  real(r_sngl) :: vtpats(ntpmax+1)
  
  real(r_sngl) :: xm4(4), ym4(4)
  real(r_sngl) :: xm5(5), ym5(5)
  real(r_sngl) :: xm7(7), ym7(7)
  character(len=4) :: cvar
  character(len=10) :: cfact
  
  integer :: ntpat, ipat
  integer :: itpat, ic
  real(r_sngl) :: val1, val2
  real(r_sngl) :: dyp, dylen, vpxr, vpyl
  
  call ueqntl(ntpat)
  
  do itpat = 1, ntpat
    call ueqtlv(val1,val2,ipat,itpat)
    itpats(itpat) = ipat 
    vtpats(itpat) = val1
    vtpats(itpat+1) = val2
  end do
  
  dyp = dylen / real(ntpat)
  
  call sglset ('LCLIP',.FALSE.) ! Cliping
  
  
  !!write(*,*) ntpat
  
  !!stop
  !!! Color bar
  
  !      do ic=1,ntpat
  do ic = 2, ntpat - 1
    xm4(1) = vpxr + 0.01
    xm4(2) = xm4(1) + 0.02
    xm4(3) = xm4(2)
    xm4(4) = xm4(1)
    ym4(1) = vpyl + real(ic-1) * dyp
    ym4(2) = ym4(1)
    ym4(3) = vpyl + real(ic)*dyp
    ym4(4) = ym4(3)
    itpat = itpats(ic)
    call sgtnzv(4,xm4,ym4,itpat)
  end do

  ic = 1
  xm4(3) = vpxr+0.01
  xm4(2) = xm4(3)+0.02
  xm4(1) = xm4(3)+0.01
  ym4(3) = vpyl+real(ic)*dyp
  ym4(2) = ym4(3)
  ym4(1) = vpyl+real(ic-1)*dyp
  call sgtnzv(3,xm4(1:3),ym4(1:3),itpats(ic))

  ic = ntpat
  xm4(1) = vpxr + 0.01
  xm4(2) = xm4(1) + 0.02
  xm4(3) = xm4(1) + 0.01
  ym4(1) = vpyl + real(ic-1) * dyp
  ym4(2) = ym4(1)
  ym4(3) = vpyl + real(ic) * dyp
  call sgtnzv(3,xm4(1:3),ym4(1:3),itpats(ic))
  
  !!! Waku
  !      xm5=(/vpxr+0.01,vpxr+0.03,vpxr+0.03,vpxr+0.01,vpxr+0.01/)
  !      ym5=(/vpyl,vpyl,vpyl+dylen,vpyl+dylen,vpyl/)
  !      call sgplzv(5,xm5,ym5,1,1) 
  xm7 = (/vpxr+0.01,vpxr+0.02,vpxr+0.03,vpxr+0.03,vpxr+0.02,vpxr+0.01,vpxr+0.01/)
  ym7 = (/vpyl+dyp,vpyl,vpyl+dyp,vpyl+dylen-dyp,vpyl+dylen,vpyl+dylen-dyp,vpyl+dyp/)
  call sgplzv(7,xm7,ym7,1,1) 
  
  
  !       do ic=2,ntpat,2
  do ic = 2, ntpat, 1
  !         if (vtpats(ic).eq.real(int(vtpats(ic)))) then
    if (vtpats(ic) >= 1.0)then
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

!==============================================================!==================================================!
