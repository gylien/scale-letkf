!==================================================!
subroutine plot_dbz_DCL_obs(nobs, ze_radar, lon_radar, lat_radar, z_radar, nlons, nlats, lons, lats, dlons, dlats, psfile) 
  use common
  use common_scale
  use scale_const, only: &
      D2R => CONST_D2R
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, KHALO
  use scale_atmos_grid_cartesC, only: &
      GRID_CXG => ATMOS_GRID_CARTESC_CXG, &
      GRID_CYG => ATMOS_GRID_CARTESC_CYG, &
      GRID_CZ  => ATMOS_GRID_CARTESC_CZ,  &
      GRID_FZ  => ATMOS_GRID_CARTESC_FZ
  use scale_atmos_grid_cartesC_real, only: &
      BASE_LON => ATMOS_GRID_CARTESC_REAL_BASEPOINT_LON, &
      BASE_LAT => ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  use scale_time, only: &
      TIME_NOWDATE
  use common_nml, only: &
      PLOT_ZLEV_MIN, PLOT_ZLEV_MAX, PLOT_ZLEV_INTV, HORI_LOCAL
  use common_mpi_scale, only: &
      mpi_timer, &
      myrank_o, &
      nprocs_o

  implicit none

  integer, intent(in) :: nobs
  real(r_sngl),intent(in) :: ze_radar(nobs)
  real(r_sngl),intent(in) :: lon_radar(nobs), lat_radar(nobs), z_radar(nobs)
  character(*),intent(in) :: psfile

  integer,intent(in) :: nlons,nlats !!! grid superob
  real(r_sngl),intent(in) :: lons(nlons), lats(nlats)
  real(r_sngl),intent(in) :: dlons, dlats

  real(r_sngl):: val_plot(nlons,nlats)

  integer :: iwork(3 * (nlons + 2) * (nlats + 2) / 2 + 1)
 
  real(r_sngl),parameter :: rmiss = -9.99e20

  integer,parameter :: npatmax = 40
  real(r_sngl) :: vtlevs(npatmax)
  integer :: itpats(npatmax)  
  character(len=40) :: title1, title2(2), title3

  real(RP) :: lon_RP, lat_RP

  integer :: iclrmap,iobs,iplot_lev
  integer :: ilon, ilat
  real(r_sngl) :: vpr, vpl, vpt, vpb

  real(r_sngl) :: range_lonl, range_lonr
  real(r_sngl) :: range_latl, range_latr
  real(r_sngl) :: aratio
  integer :: i, ntpat
  real(r_sngl) :: amtics, astics
  real(r_sngl) :: bmtics, bstics

  integer:: iblkge, itpat

  character(len=19) :: ftimelabel

  character(len=5)::cheight

  real(r_sngl),parameter::rmiss_radar=-327.68

  integer :: pcnt !!! process count

  real(r_sngl)::lon4(4), lat4(4)
  real(r_sngl) :: zmin, zmax

  integer time1, time2, timerate, timemax

  call system_clock(time1, timerate, timemax)

  write(ftimelabel,'(I4.4,A1,I2.2,A1,I2.2,A1, I2.2,A1,I2.2,A1,I2.2)')&
       TIME_NOWDATE(1), '/', TIME_NOWDATE(2), '/' ,TIME_NOWDATE(3), ' ', &
       TIME_NOWDATE(4), ':', TIME_NOWDATE(5), ':', TIME_NOWDATE(6)



  pcnt = 0
  do iplot_lev = PLOT_ZLEV_MIN, PLOT_ZLEV_MAX, PLOT_ZLEV_INTV
    if ( GRID_CZ(iplot_lev+KHALO) > RADAR_ZMAX) cycle !!! exclude stratosphere

    pcnt = pcnt + 1
    if ( mod(pcnt, nprocs_o) /= myrank_o ) cycle
write(*,*) "DEBUG",myrank_o, pcnt, iplot_lev, nprocs_o
    write(cheight,'(I5.5)')int(GRID_CZ(iplot_lev+KHALO))
   
    title1 = trim(ftimelabel) // " UTC (PAWR obs)"
    title2 = (/'',''/)
    write(title3,'(A,I5,A)') 'Radar reflectivity ',int(GRID_CZ(iplot_lev+KHALO)),' m' 
   
    iclrmap = 12
     

    call MAPPROJECTION_xy2lonlat( GRID_CXG(IHALO+1), GRID_CYG(1), lon_RP, lat_RP )
    range_lonl=real(lon_RP/D2R)
    call MAPPROJECTION_xy2lonlat( GRID_CXG(IHALO+nlong), GRID_CYG(1), lon_RP, lat_RP )
    range_lonr=real(lon_RP/D2R)
  
  
    call MAPPROJECTION_xy2lonlat( GRID_CXG(1),GRID_CYG(JHALO+1), lon_RP, lat_RP )
    range_latl=real(lat_RP/D2R)
    call MAPPROJECTION_xy2lonlat( GRID_CXG(1),GRID_CYG(JHALO+nlatg), lon_RP, lat_RP )
    range_latr=real(lat_RP/D2R)
  
    iclrmap = 12
  
  
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
    call grswnd(range_lonl,range_lonr,range_latl,range_latr)
   
    call grsvpt(vpl,vpr,vpb,vpt)
    call grstrn(10) !!! Mercator
    call umlset('LGLOBE',.false.)
    call umiset('INDEXOUT',31)
    call umpfit
   
    call grstrf
   
 
    call glrset ('RMISS', rmiss)
    call gllset ('LMISS', .true.)
     
    call ueitlv
   
    ntpat = 10
    itpats(1:ntpat+2) = (/ 0, 40,34,30,50,62,68,74,80,84,92,98/) * 1000 + 999
    vtlevs(2:ntpat+2) = 5.0 + (/( 5.0*real(i), i=1,11 )/)
    itpats(1) = 0
    vtlevs(1) = -1.0e10
    vtlevs(ntpat+3) = 1.0e10
    call uestln(vtlevs(1:ntpat+3),itpats(1:ntpat+2),ntpat+2)
 

    val_plot=0.0

    zmin = real( GRID_FZ(iplot_lev+KHALO-1) - RADAR_SO_SIZE_VERT )
    zmax = real( GRID_FZ(iplot_lev+KHALO-1) + RADAR_SO_SIZE_VERT )

    call system_clock(time2, timerate, timemax)
    if (myrank_o == 1 ) write(*, *) "plot setting", (time2 - time1) / dble(timerate), myrank_o
    time1 = time2

    do iobs = 1, nobs
      if ( z_radar(iobs) >= zmin .and. & 
          z_radar(iobs) < zmax  ) then
          ilon = int((lon_radar(iobs)-lons(1)-0.5*dlons)/dlons)+1
          ilat = int((lat_radar(iobs)-lats(1)-0.5*dlats)/dlats)+1
          if (ilon.ge.1.and.ilon.le.nlons.and.ilat.ge.1.and.ilat.le.nlats) then
           itpat=iblkge( vtlevs(1:ntpat+3),ntpat+3,10.0*log10(max(ze_radar(iobs),1.0e-10)))
            if (itpat >= 2 .and. itpat <= ntpat+1 )then
             call sgtnzu(4, &
                         (/ lons(ilon)-0.5*dlons,lons(ilon)+0.5*dlons,lons(ilon)+0.5*dlons,lons(ilon)-0.5*dlons/), &
                         (/ lats(ilat)-0.5*dlats,lats(ilat)-0.5*dlats,lats(ilat)+0.5*dlats,lats(ilat)+0.5*dlats/), &
                         itpats(itpat))
            end if
          end if
     end if
    end do
    call system_clock(time2, timerate, timemax)
    if (myrank_o == 1 ) write(*, *) "plot obs loop", (time2 - time1) / dble(timerate), myrank_o
    time1 = time2

    call uwsgxa (lons,nlons)
    call uwsgya (lats,nlats)
!    call uetone (val_plot,nlons,nlons,nlats)
   
    call dcbar(vpr+0.02,vpb,(vpt-vpb)*0.8)
   
    call udlset ('LMSG',.false.)
    call udlset ('LABEL',.false.)
   
    call system_clock(time2, timerate, timemax)
    if (myrank_o == 1 ) write(*, *) "plot chk1", (time2 - time1) / dble(timerate), myrank_o
    time1 = time2
!!! map  
 
    call sglset('LCLIP',.true.)
    call umlset ('LGRIDMJ',.false.)
    call umrset ('DGRIDMN',0.5)
    call umiset ('ITYPEMN',3)
    call umiset ('INDEXMN',1)
  
    call umpglb
    call umplim
    call umpmap('coast_japan')
  
    call uumrkz(1,real(BASE_LON/D2R),real(BASE_LAT/D2R),9,21,0.010) !!! RADAR location
  
    call system_clock(time2, timerate, timemax)
    if (myrank_o == 1 ) write(*, *) "plot chk2", (time2 - time1) / dble(timerate), myrank_o
    time1 = time2

    amtics = 0.5 !! deg
    astics = 0.5
    bmtics = 0.5
    bstics = 0.5
  
     
    call sglset('LCLIP',.false.)
    call uzinit
    call uziset('INDEXT2',3)
    call uziset('INDEXT1',3)
    call uziset('INNER',-1)
    call uzrset('RSIZEL1',0.016)
    call uzrset('RSIZEC1',0.016)
    call uzrset('RSIZET1',0.006)
    call uzrset('RSIZET2',0.003)
    call uxsfmt ('(F5.1)')
    call uysfmt ('(F5.1)')
   
    call system_clock(time2, timerate, timemax)
    if (myrank_o == 1 ) write(*, *) "plot chk3", (time2 - time1) / dble(timerate), myrank_o
    time1 = time2
       
    call uxaxdv('B',astics,amtics)
    call uxaxdv('T',astics,amtics)
    call uxsttl('B','Longitude',0.0)
    call uyaxdv('L',bstics,bmtics)
    call uyaxdv('R',bstics,bmtics)
    call uysttl('L','Latitude',0.0)
  
    call system_clock(time2, timerate, timemax)
    if (myrank_o == 1 ) write(*, *) "plot chk4", (time2 - time1) / dble(timerate), myrank_o
    time1 = time2

    call uzlset('LABELYR',.false.)
   
    call sglset('LCLIP',.false.)
    
    call sgtxzv (0.5*(vpr+vpl),vpt+0.02,trim(title1),0.016,0,0,3) !
    call sgtxzv (vpr-0.01,vpt+0.045,trim(title2(1)),0.016,0,1,3) !
    call sgtxzv (vpr-0.01,vpt+0.020,trim(title2(2)),0.016,0,1,3) !
    call sgtxzv (0.5*(vpr+vpl),vpt+0.05,trim(title3),0.017,0,0,4) !
   
    call system_clock(time2, timerate, timemax)
    if (myrank_o == 1 ) write(*, *) "plot chk5", (time2 - time1) / dble(timerate), myrank_o
    time1 = time2

    call grcls 
   
    call system_clock(time2, timerate, timemax)
    if (myrank_o == 1 ) write(*, *) "plot chk6", (time2 - time1) / dble(timerate), myrank_o
    time1 = time2

  end do !!! iplot_lev


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
