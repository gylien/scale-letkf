!==================================================!
subroutine plot_dbz_DCL_obs(nobs,ze_radar,lon_radar,lat_radar,z_radar,psfile) 
  use common
  use common_scale
  use scale_io, only: &
      H_LONG
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, KHALO
  use scale_atmos_grid_cartesC, only: &
      DX, DY, &
      GRID_CXG => ATMOS_GRID_CARTESC_CXG, &
      GRID_CYG => ATMOS_GRID_CARTESC_CYG, &
      GRID_CZ  => ATMOS_GRID_CARTESC_CZ,  &
      GRID_FZ  => ATMOS_GRID_CARTESC_FZ
  use scale_time, only: &
      TIME_NOWDATE
  use common_nml, only: &
      plot_zlev_min, plot_zlev_max, plot_zlev_intv

  implicit none

  integer, intent(in) :: nobs
  real(r_sngl),intent(in) :: ze_radar(nobs)
  real(r_sngl),intent(in) :: lon_radar(nobs), lat_radar(nobs), z_radar(nobs)
  character(*),intent(in) :: psfile


  real(r_sngl),allocatable :: val_plot(:,:)
  integer,allocatable :: icount_data(:,:)

  integer,allocatable :: iwork(:)
  
  integer :: nlonadd 
  integer :: nlatadd 
  real(r_sngl),parameter :: OFFXY = 0.0e3 ! offset (m)
  real(r_sngl),allocatable :: grid_cxg_ext(:)
  real(r_sngl),allocatable :: grid_cyg_ext(:)
!  real(r_sngl),allocatable :: vmask(:,:)
  
  real(r_sngl),parameter :: rmiss = -9.99e20
  
  integer,parameter :: npatmax = 40
  real(r_sngl) :: vtlevs(npatmax)
  integer :: itpats(npatmax)  
  character(len=40) :: title1, title2(2), title3

  integer :: iclrmap,iobs,iplot_lev
  integer :: ilon, ilat, nwork, nlong_ext, nlatg_ext
  real(r_sngl) :: vpr, vpl, vpt, vpb

  real(r_sngl) :: range_lonl, range_lonr
  real(r_sngl) :: range_latl, range_latr
  real(r_sngl) :: aratio
  integer :: i, ntpat
  real(r_sngl) :: amtics, astics
  real(r_sngl) :: bmtics, bstics

  integer:: iblkle

  character(len=19) :: ftimelabel

  character(len=5)::cheight

  real(r_sngl),parameter::rmiss_radar=-327.68

  integer,parameter::ismth=1

  include 'latlon_d4_grid.h'

  if (nlong_fix.ne.nlong.or.nlatg_fix.ne.nlatg) then
   write(*,*) 'nlong or nlatg does not match'
   stop
  end if

  nlonadd = int(OFFXY / DX)
  nlatadd = int(OFFXY / DY)

  write(ftimelabel,'(I4.4,A1,I2.2,A1,I2.2,A1, I2.2,A1,I2.2,A1,I2.2)')&
       TIME_NOWDATE(1), '/', TIME_NOWDATE(2), '/' ,TIME_NOWDATE(3), ' ', &
       TIME_NOWDATE(4), ':', TIME_NOWDATE(5), ':', TIME_NOWDATE(6)

  nlong_ext = nlong + nlonadd * 2
  nlatg_ext = nlatg + nlatadd * 2 

  allocate(icount_data(nlong,nlatg))
  allocate(val_plot(nlong,nlatg))
  allocate(grid_cxg_ext(nlong_ext))
  allocate(grid_cyg_ext(nlatg_ext))
!  allocate(vmask(nlong_ext,nlatg_ext))

  nwork = 3 * (nlong + 2) * (nlatg + 2) / 2 + 1
  allocate(iwork(nwork))


do iplot_lev=plot_zlev_min,plot_zlev_max,plot_zlev_intv
if ( grid_cz(iplot_lev+KHALO) .le. RADAR_ZMAX) then !!! exclude stratosphere

 !!! def val_plot, cheight
 icount_data=0
 val_plot=0.0
 do iobs=1,nobs
 if (z_radar(iobs).ne.rmiss_radar) then
  if ( z_radar(iobs) .ge. grid_fz(iplot_lev+KHALO-1) .and. & 
       z_radar(iobs) .lt. grid_fz(iplot_lev+KHALO) ) then
  ilon=iblkle(axlon,nlong,lon_radar(iobs))
  ilat=iblkle(axlatSN,nlatg,lat_radar(iobs))
!   if (ilon.lt.nlong.and.abs(lon_radar(iobs)-axlon(ilon)).gt.abs(lon_radar(iobs)-axlon(ilon+1))) ilon=ilon+1
!   if (ilat.lt.nlatg.and.abs(lat_radar(iobs)-axlatSN(ilat)).gt.abs(lat_radar(iobs)-axlatSN(ilat+1))) ilat=ilat+1
   if (ilon.le.nlong-ismth.and.ilon.ge.1+ismth.and.ilat.le.nlatg-ismth.and.ilat.ge.1+ismth)then
   if (ze_radar(iobs).ne.rmiss_radar.and.ze_radar(iobs).gt.9.9999998e-3)then
   icount_data(ilon-ismth:ilon+ismth,ilat-ismth:ilat+ismth) = icount_data(ilon-ismth:ilon+ismth,ilat-ismth:ilat+ismth) + 1
    val_plot(ilon-ismth:ilon+ismth,ilat-ismth:ilat+ismth)  = val_plot(ilon-ismth:ilon+ismth,ilat-ismth:ilat+ismth)  + 10.0 * log10(ze_radar(iobs))
   end if
   end if
  end if
 end if
end do
write(cheight,'(I5.5)')int(grid_cz(iplot_lev+KHALO))

where(icount_data.ne.0) val_plot=val_plot/real(icount_data)

  title1 = trim(ftimelabel) // " UTC (PAWR obs)"
  title2 = (/'',''/)
 write(title3,'(A,I5,A)') 'Radar reflectivity ',int(grid_cz(iplot_lev+KHALO)),' m' 
!  if (cheight(1:1) == '0') title3='Radar reflectivity z=' // cheight(2:5) // ' m'
!  if (cheight(1:2) == '00') title3='Radar reflectivity z=' // cheight(3:5) // ' m'

  iclrmap = 12
  
!  do ilon = 1, nlong_ext
!  do ilat = 1, nlatg_ext
!    if (ilon - nlonadd >= 1 .and. ilon - nlonadd < nlong .and. &
!       ilat - nlatadd >= 1 .and. ilat - nlatadd < nlatg  )then
!      if (val_plot(ilon,ilat) /= rmiss)then
!        vmask(ilon,ilat) = 0.0
!      else
!        vmask(ilon,ilat) = rmiss
!      end if
!    else
!      vmask(ilon,ilat) = rmiss
!    end if
!  end do
!  end do
 
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
  call swcset('FNAME',trim(psfile)//'_z'//trim(cheight)//'m')
  call swlset('LSEP',.false.)
  call swiset('IFL',1) !!! PNG
  call swiset('IWIDTH',1000)
  call swiset('IHEIGHT',800)
  call gropn(2)
  call sglset('LFULL',.true.)
  call sglset('LCLIP',.true.)
  call slmgn(0.0,0.0,0.0,0.0)
  call grfrm


!!! MAP
  call grswnd(axlon(1), axlon(nlong), axlatSN(1),axlatSN(nlatg))

  call grsvpt(vpl,vpr,vpb,vpt)
  call grstrn(10) !!! Mercator
  call umlset('LGLOBE',.false.)
  call umiset('INDEXOUT',31)
!  call umscnt (0.5*(axlon(1)+axlon(nlong)),0.5*(axlatSN(1)+axlatSN(nlatg)),0.0)
  call umpfit

  call grstrf

  call uwsgxa (axlon,nlong)
  call uwsgya (axlatSN,nlatg)

  call sglset('LCLIP',.true.) 
  call umlset ('LGRIDMJ',.false.)
  call umrset ('DGRIDMN',0.5)
!  call umiset ('ITYPEMN',3)
  call umiset ('INDEXMN',1)
  call umpglb
  call umplim
  call umpmap('coast_japan')

!!!!

  call grswnd(range_lonl,range_lonr,range_latl,range_latr)

  call grsvpt(vpl,vpr,vpb,vpt)
  call grstrn(1) !!! Give up map proj

  call grstrf

  call glrset ('RMISS', rmiss)
  call gllset ('LMISS', .true.)

  call uwsgxa (grid_cxg_ext(nlonadd+1:nlong_ext-nlonadd),nlong)
  call uwsgya (grid_cyg_ext(nlatadd+1:nlatg_ext-nlatadd),nlatg)

  call ueitlv

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

!!! Radar location  
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

 end if
 end do !!! iplot_lev

  deallocate(icount_data,val_plot,grid_cxg_ext,grid_cyg_ext,iwork)


return

contains
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

end subroutine plot_dbz_DCL_obs
!==============================================================!==================================================!