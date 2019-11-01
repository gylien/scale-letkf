module radar_obs
!=======================================================================
!
! [PURPOSE:] Radar observations
!
!=======================================================================
!$USE OMP_LIB
  use common
  use common_nml
  use common_scale
  use common_obs_scale
  use common_mpi_scale, only: &
    MPI_COMM_o, &
    myrank_a, myrank_o, &
    mpi_timer!, &
    !pawr_toshiba_scattv_mpi, &
    !pawr_3dvar_allreduce

  implicit none
  public
  real(r_size), allocatable, save :: radlon(:, :, :), radlat(:, :, :), radz(:, :, :)
  integer, save :: utime_obs(6) = (/-1,-1,-1,-1,-1,-1/)

contains

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY_radar(elm,radar_lon,radar_lat,radar_z,ri,rj,rk,lon,lat,lev,v3d,v2d,yobs,qc,stggrd)
  use scale_mapprojection, only: &
      MAPPROJECTION_rotcoef
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk,radar_lon,radar_lat,radar_z !!!!! Use only, ri, rj, rk eventually... (radar_lon,lat,z in ri,rj,rk)
  REAL(r_size),INTENT(IN) :: lon,lat,lev
  REAL(r_size),INTENT(IN) :: v3d(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  REAL(r_size),INTENT(OUT) :: yobs
  INTEGER,INTENT(OUT) :: qc
  INTEGER,INTENT(IN),OPTIONAL :: stggrd
  INTEGER :: stggrd_ = 0

  REAL(r_size) :: qvr,qcr,qrr,qir,qsr,qgr,ur,vr,wr,tr,pr !,rhr
  REAL(r_size) :: dist , dlon , dlat , az , elev , radar_ref,radar_rv

  REAL(r_size) :: utmp, vtmp
  REAL(RP) :: rotc(1,1,2)
  real(RP) :: lon_tmp(1,1),lat_tmp(1,1)

!  integer :: ierr
!  REAL(r_dble) :: rrtimer00,rrtimer
!  rrtimer00 = MPI_WTIME()


  if (present(stggrd)) stggrd_ = stggrd


  yobs = undef
  qc = iqc_good

  if (stggrd_ == 1) then
    CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri-0.5_r_size,rj,ur)  !###### should modity itpl_3d to prevent '1.0' problem....??
    CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj-0.5_r_size,vr)  !######
    CALL itpl_3d(v3d(:,:,:,iv3dd_w),rk-0.5_r_size,ri,rj,wr)  !######
  else
    CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri,rj,ur)
    CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj,vr)
    CALL itpl_3d(v3d(:,:,:,iv3dd_w),rk,ri,rj,wr)
  end if
  CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,tr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_p),rk,ri,rj,pr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,qvr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qc),rk,ri,rj,qcr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qr),rk,ri,rj,qrr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qi),rk,ri,rj,qir)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qs),rk,ri,rj,qsr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qg),rk,ri,rj,qgr)


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:itpl_3d:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !Compute az and elevation for the current observation.
  !Simple approach (TODO: implement a more robust computation)

  !Azimuth
  dlon=lon-radar_lon
  dlat=lat-radar_lat
  IF ( dlon == 0.0d0 .and. dlat == 0.0d0  )THEN
!      WRITE(6,*)'OOPS',dlon,dlat,lon,lat,radar_lon,radar_lat
    qc = iqc_out_h
    RETURN
  ELSE
    az = rad2deg*atan2(dlon*cos(radar_lat*deg2rad),dlat)
  ENDIF
  !WRITE(6,*)dlon,dlat,lon,lat,radar_lon(ityp),radar_lat(ityp),dlon*cos(radar_lat(ityp))
  !IF( abs(dlon) > maxdlon )maxdlon=abs(dlon)
  !IF( abs(dlat) > maxdlat )maxdlat=abs(dlat)
  !WRITE(6,*)maxdlon,maxdlat
  IF( az < 0) az = 360.0d0 + az
  !elevation
  CALL com_distll_1(lon,lat,radar_lon,radar_lat,dist)
  elev=rad2deg*atan2(lev-radar_z,dist)

  lon_tmp(1,1) = lon*deg2rad
  lat_tmp(1,1) = lat*deg2rad
  call MAPPROJECTION_rotcoef(1, 1, 1, 1, 1, 1, &
                             lon_tmp(1,1),lat_tmp(1,1),rotc)

  utmp = ur
  vtmp = vr
  ur = utmp * rotc(1,1,1) - vtmp * rotc(1,1,2)
  vr = utmp * rotc(1,1,2) + vtmp * rotc(1,1,1)

!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:radar_coord:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !Check that the azimuth and elevation angles are within the expected range.
  !Some grid points may be at the radar location.
  !IF( .NOT. ( az .GT. 0.0d0 .AND. az .LT. 360.0d0 ) )RETURN
  !IF( .NOT. ( elev .GT. 0.0d0 .AND. elev .LT. 90.0d0 ))RETURN



  !DEBUG---------------------------------------------------------------
  !WRITE(6,*)'RADAR lat ',radar_lat(ityp),' RADAR lon ',radar_lon(ityp)
  !WRITE(6,*)'OBS   lat ',dlat      ,' OBS   lon ',dlon
  !WRITE(6,*)'AZIMUTH   ',az
  !WRITE(6,*)'RADAR z   ',radar_z ,'  OBS   z   ',lev
  !WRITE(6,*)'ELEVATION ',elev
  !DEGUB---------------------------------------------------------------

  !WRITE(6,*)'BCRV',dlon,dlat,az,elev

  CALL calc_ref_vr(qvr,qcr,qrr,qir,qsr,qgr,ur,vr,wr,tr,pr,az,elev,radar_ref,radar_rv)

!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:calc_ref_vr:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  SELECT CASE (elm)
  CASE(id_radar_ref_obs,id_radar_ref_zero_obs)
!!!!    if (radar_ref < MIN_RADAR_REF) then
!!!!      !In this case we will replace the observation by -RH
!!!!      !This allows us to use pseudo rh observations in some cases.
!!!!      !Later we will take the decision on what to do with these cases...

!!!!      CALL itpl_3d(v3d(:,:,:,iv3dd_rh),rk,ri,rj,rhr)
!!!!!      qc = 
!!!!      yobs = -rhr

!!!!    else                      !!!!!! --------- Pesudo RH: TO BE DONE...
    if (radar_ref < MIN_RADAR_REF) then
      qc = iqc_ref_low
      yobs = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT  !!! even if the above qc is bad, still return the value
    else
      yobs = 10.0d0 * log10(radar_ref)
    end if
!!!!    end if
  CASE(id_radar_vr_obs)
    if (radar_ref < MIN_RADAR_REF) then
      qc = iqc_ref_low
    end if
    yobs = radar_rv  !!! even if the above qc is bad, still return the value
  CASE DEFAULT
    qc = iqc_otype
  END SELECT


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:conversion:',rrtimer-rrtimer00
!  rrtimer00=rrtimer

  RETURN
END SUBROUTINE Trans_XtoY_radar

!-----------------------------------------------------------------------
! Compute radar reflectivity and radial wind.
! Radial wind computations for certain methods depend on model reflectivity
! so both functions has been merged into a single one.
! First reflectivity is computed, and the the radial velocity is computed.
!-----------------------------------------------------------------------
SUBROUTINE calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,az,elev,ref,vr)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: qv        !Water vapor
  REAL(r_size), INTENT(IN) :: qc,qr     !Cloud and rain water
  REAL(r_size), INTENT(IN) :: qci,qs,qg !Cloud ice, snow and graupel
  REAL(r_size), INTENT(IN) :: t,p       !Temperature and pressure.
  REAL(r_size), INTENT(IN) :: u,v,w     !velocities with respecto to earth.
  REAL(r_size), INTENT(INOUT) :: ref      !Reflectivity.
  REAL(r_size)              :: ro
  REAL(r_size), INTENT(IN) :: az     !Azimuth respect to the radar.
  REAL(r_size), INTENT(IN) :: elev   !Elevation angle respect to the surface.
  REAL(r_size), INTENT(INOUT) :: vr    !Radial velocity.
  REAL(r_size)  :: qms , qmg !Melting species concentration (METHOD_REF_CALC 3)
  REAL(r_size)  :: qt        !Total condensate mixing ratio (METHOD_REF_CALC 1)
  REAL(r_size)  :: zr , zs , zg !Rain, snow and graupel's reflectivities.
  REAL(r_size)  :: wr , ws , wg !Rain, snow and graupel's mean terminal velocities.
  REAL(r_size)  :: wt           !Total mean terminal velocity.
  REAL(r_size)  :: nor, nos, nog !Rain, snow and graupel's intercepting parameters.
  REAL(r_size)  :: ror, ros, rog , roi !Rain, snow and graupel, ice densities.
  REAL(r_size)  :: a,b,c,d,Cd    !Constant for fall speed computations.
  REAL(r_size)  :: cf, pip , roo
  REAL(r_size)  :: ki2 , kr2
  REAL(r_size)  :: lr , ls , lg
  REAL(r_size)  :: tmp_factor , rofactor
  REAL(r_size)  :: p0
  REAL(r_size)  :: Fs, Fg , zms , zmg , fws , fwg !METHOD_REF_CALC 3
  REAL(r_size)  :: qrp , qsp , qgp
  REAL(r_size)  :: maxf                     !Maximum mixture relative concentration. (METHOD_REF_CALC 3)

  !Note: While equivalent reflectivity is assumed to be independent of the radar, in 
  !practice short wavelengths as those associated with K band radars migh frequently
  !experience Mie scattering. In that case, the equivalent reflectivity is not longer
  !radar independent and an appropiate relationship between the forecasted concentrations
  !and the reflectivity should be used.


  !REAL(r_size)  :: trqr,trqs,trqg
  ![P] Pa
  ![T] K
  ![q...] dimensionless
  ![ref] mm^6/m^3
  ![refdb] refdb
  ![U,V,W] m/s
  ![az, elev] degree
  ![vr] m/s

  !This model reflectivity won't be lower than this value.

  !Initialize reflectivities
  zr=0.0d0
  zs=0.0d0
  zg=0.0d0
  zms=0.0d0
  zmg=0.0d0
  ref=0.0d0

  !Compute air density (all methods use this)

  ro =  p / (rd * t)

  !Begin computation of reflectivity and vr

  if (METHOD_REF_CALC == 1) then

    !WRF method: See for example Sugimoto et al. Evaluation of the Performance of Ra
    !dial Velocity Assimilation with WRF-3DVAR System and Simulated Multiple-Doppler
    !Radar Data
    !Only rain is used to estimate the terminal velocity of hidrometeors.
    !Only rain is used to compute equivalent reflectivity.
    !Marshall-Palmer distributions are assumed to find the relationship bestween
    !concentration and equivalent reflectivity.
    ! Sun and Crook 1997 , 1998.
    !Derived for C-band radars.
    !Attenuation is not computed.

    !Reflectivity
    nor=8.0d6      ![m^-4]
    ror=1000.0d0   ![Kg/m3]
    pip=pi ** 1.75 !factor
    cf =10.0d18 * 72 !factor
    p0=1.0d5            !Reference pressure.

    qt=qr + qs + qg  !Assume that the total condensate is rain water
                     !But ignore cloud ice and cloud water

    IF( qt .GT. 0.0d0 )THEN
    ref = cf * ( ( ro * qt )**1.75 )
    ref = ref / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
   !ref= 2.04d4 *( ( ro * qt * 1.0d3 ) ** 1.75 ) !Original Sun and Crook expresion.
    ELSE
    ref=0.0d0
    ENDIF

    !Radial wind

    IF ( qt .GT. 0.0d0 )THEN
    a=(p0/p)**0.4
    wt = 5.40d0 * a * ( qt ** 0.125 )
    ELSE
    wt=0d0
    ENDIF
   !WRITE(6,*)qr,qs,qg


  else if (METHOD_REF_CALC == 2) then
    !Observation operator from Tong and Xue 2006, 2008 a and b.
    !Based on Smith et al 1975.
    !It includes reflectivity contribution by all the microphisical species.
    !is assumes Marshall and Palmer distributions.
    !Based on C band radars.
    !Attenuation is not computed.
    nor=8.0d6      ![m^-4]
    nos=3.0d6      ![m^-4]
    nog=4.0d4      ![m^-4]
    ror=1000.0d0   ![Kg/m3]
    ros=100.0d0    ![Kg/m3]
    rog=913.0d0    ![Kg/m3] 
    roi=917.0d0    ![Kg/m3]
    roo=1.0d0      ![Kg/m3] Surface air density.
    ki2=0.176d0    !Dielectric factor for ice.
    kr2=0.930d0    !Dielectric factor for water.
    pip=pi ** 1.75 !factor
    cf =1.0d18 * 720 !factor 

    IF( qr .GT. 0.0d0 )THEN
    zr= cf * ( ( ro * qr )**1.75 )
    zr= zr / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
    ENDIF
    !The contribution of snow depends on temperature (bright band effect)
    IF( qs .GT. 0.0d0 )THEN
    IF ( t <= 273.16 )THEN
     zs = cf * ki2 * ( ros ** 0.25 ) * ( ( ro * qs ) ** 1.75 )
     zs = zs / ( pip * kr2 * ( nos ** 0.75  ) * ( roi ** 2 ) )
    ELSE
    !WARNING: This formulation has to be checked the paper says that 
     !ros instead of roi should be used in thes expresion, but that 
     !leads to unrealistic jumps in the model derived reflectivity.
     zs = cf * ( ( ro * qs ) ** 1.75 )
     zs = zs / ( pip * ( nos ** 0.75 ) * ( roi ** 1.75 ) )
    ENDIF
    ENDIF

    !Only dry graupel contribution is ussed.
    IF( qg .GT. 0.0d0 )THEN
    zg= ( cf / ( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) ) ** 0.95
    zg= zg * ( ( ro * qg ) ** 1.6625 )

    !zg=(ki2/kr2)*( cf / ( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) ) 
    !zg= zg * ( ( ro * qg ) ** 1.75 )

    !if( qg * ro * 1000 .GT. 1.0d0 )THEN
    !  WRITE(6,*)qg, zg, ro, rd , t , p, 10*log10(zg)
    !endif
    ENDIF

    ref = zr + zs + zg

    !Compute reflectivity weigthed terminal velocity.
    !Lin et al 1983.
    IF( ref > 0.0d0 )THEN
    !There are hidrometeors, compute their terminal velocity.
    !Change units to be consistent with Lin et al 1983 and
    !to obtain wt in m/s
    nor=nor*1e-3      ![cm^-4]
    nos=nos*1e-3      ![cm^-4]
    nog=nog*1e-3      ![cm^-4]
    ror=ror*1e-3        ![g/cm3]
    ros=ros*1e-3        ![g/cm3]
    rog=rog*1e-3      ![g/cm3] 
    roo=roo*1e-3      ![g/cm3] Surface air density.
    ro= ro*1e-3

    a=2115d0   ![cm**1-b / s]
    b=0.8d0
    c=152.93d0 ![cm**1-b / s]
    d=0.25d0
    Cd=0.6d0

    rofactor= ( roo / ro  ) ** 0.25
    if(qr > 0.0d0)then
      CALL com_gamma( 4.0_r_size + b , tmp_factor )
      lr= ( pi * ror * nor / ( ro * qr ) ) ** 0.25
      wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
      wr= 1.0d-2*wr * rofactor
    else
      wr = 0.0d0
    endif

    if(qs > 0.0d0)then
      CALL com_gamma( 4.0_r_size + d , tmp_factor )
      ls= ( pi * ros * nos / ( ro * qs ) ) ** 0.25
      ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
      ws= 1.0d-2*ws * rofactor
    else
      ws = 0.0d0
    endif
 
    if(qg > 0.0d0)then
      CALL com_gamma( 4.5_r_size , tmp_factor )
      lg= ( pi * rog * nog / ( ro * qg ) ) ** 0.25
      wg= tmp_factor * ( ( ( 4.0d0 * gg * 100.0d0 * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
      wg= 1.0d-2*wg / ( 6.0d0 * ( lg ** 0.5 ) )
    else
      wg = 0.0d0
    endif

    !Reflectivity weighted terminal velocity. 
    wt = ( wr * zr + ws * zs + wg * zg )/ ( zr + zs + zg )

    ELSE

    wt=0.0d0

    ENDIF

  else if (METHOD_REF_CALC == 3) then
    !Observation operator from Xue et al 2007
    !Asumes power law between Z and mass concentration of different 
    !hydrometeor categories. 
    !Derived for X-Band radars tacking into account Mie scattering.
    !Includes a computation of the attenuation. (k)

    MAXF=0.5d0

    !First we need to compute the mixtures between rain, snow and hail
    !Following Jung et al 2007 eq (2) and (3)
    Fg=0.0d0
    Fs=0.0d0
    fwg=0.0d0
    fws=0.0d0
    IF( qr .GT. 0.0d0 .AND. qg .GT. 0.0d0)THEN
      Fg=MAXF * ( min( qr/qg , qg/qr ) )**(1.0d0/3.0d0)
      fwg= qr / ( qr + qg )
    ENDIF
    IF( qr .GT. 0.0d0 .AND. qs .GT. 0.0d0)THEN
      Fs=MAXF * ( min( qr/qs , qs/qr ) )**(1.0d0/3.0d0)
      fws= qr / ( qr + qs )
    ENDIF


    !Correct the rain, snow and hail mixing ratios assuming
    !that we have a mixture due to melting.

    qrp=(1.0d0-Fs-Fg)*qr

    qsp=(1.0d0-Fs)*qs

    qgp=(1.0d0-Fg)*qg

    !Compute the new species resulting from melting.

    qms=Fs * (qr + qs) !Melting snow concentration.

    qmg=Fg * (qr + qg) !Melting hail concentration.

    !Compute reflectivities for each species including the melting species.

    IF( qrp .GT. 0.0d0)THEN
    zr= 2.53d4 * ( ro * qrp * 1.0d3 )**1.84
    ENDIF
    IF( qsp .GT. 0.0d0)THEN
    zs= 3.48d3 * ( ro * qsp * 1.0d3 )**1.66
    ENDIF
    IF( qgp .GT. 0.0d0)THEN
!!!    zg= 8.18d4 * ( ro * qgp * 1.0d3 )**1.50 !!! hail (Xue et al. 2009)
    zg= 5.54d3 * ( ro * qgp * 1.0d3 )**1.70  !!! graupel (A.Amemiya 2019)
    ENDIF
    IF( qms .GT. 0.0d0 )THEN
    zms=( 0.00491 + 5.75*fws - 5.588*(fws**2) )*1.0d5
    zms= zms * ( ro * qms * 1.0d3 )**( 1.67 - 0.202*fws + 0.398*(fws**2) )

    ENDIF
    IF( qmg .GT. 0.0d0 )THEN
    zmg=( 0.809 + 10.13*fwg -5.98*(fwg**2) )*1.0d5
    zmg= zmg * ( ro * qmg * 1.0d3 )**( 1.48 + 0.0448*fwg - 0.0313*(fwg**2) )
    ENDIF

    ref = zr +  zg  + zs + zms + zmg

    !Compute reflectivity weigthed terminal velocity.
    !Lin et al 1983. (The distribution parameters are 
    !consistent with the work of Jung et al 2007)

    IF( ref > 0.0d0 )THEN
      !There are hidrometeors, compute their terminal velocity.
      !Units according to Lin et al 1983.
      nor=8.0d-2      ![cm^-4]
      nos=3.0d-2      ![cm^-4]
      nog=4.0d-4      ![cm^-4]
      ror=1.0d0        ![g/cm3]
      ros=0.1d0        ![g/cm3]
      rog=0.917d0      ![g/cm3] 
      roo=0.001d0      ![g/cm3] Surface air density.
      ro=1.0d-3 * ro
      a=2115d0   ![cm**1-b / s]
      b=0.8d0
      c=152.93d0 ![cm**1-b / s]
      d=0.25d0
      Cd=0.6d0

      rofactor= ( roo / ro  ) ** 0.5

      IF ( qr .GT. 0.0d0 )THEN
      CALL com_gamma( 4.0_r_size + b , tmp_factor )
      lr= ( pi * ror * nor / ( ro * qr ) ) ** 0.25
      wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
      wr= 1.0d-2 * wr * rofactor
      ELSE
      wr=0.0d0
      ENDIF

      IF( qs .GT. 0.0d0 )THEN
      ls= ( pi * ros * nos / ( ro * qs ) ) ** 0.25
      CALL com_gamma( 4.0_r_size + d , tmp_factor )
      ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
      ws= 1.0d-2 * ws * rofactor
      ELSE
      ws=0.0d0
      ENDIF

      IF ( qg .GT. 0.0d0 )THEN
      lg= ( pi * rog * nog / ( ro * qg ) ) ** 0.25
      CALL com_gamma( 4.5_r_size , tmp_factor )
      wg= tmp_factor * ( ( ( 4.0d0 * gg * 100.0d0 * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
      wg= 1.0d-2 * wg / ( 6.0d0 * ( lg ** 0.5 ) )
      ELSE
      wg=0.0d0
      ENDIF

      !Reflectivity weighted terminal velocity. 
      !The melting species are assumed to fail as their non-melting counterparts.
      !however this might not be the case for melting snow.
      wt = ( wr * zr + ws *  zs +  ws * zms + wg *  zg + wg * zmg ) / ( zr + zs + zg + zms + zmg )

    ELSE

      wt=0.0d0

    ENDIF

  else  !IF OVER DIFFERENT OPTIONS

    WRITE(6,*)'[Error] Not recognized method for radar reflectivity and wind computation'
    STOP

  end if ! [METHOD_REF_CALC == ?]


  !Compute radial velocity
  !WRITE(6,*)'ICRV',u,v,w,wt,az,elev
  vr = u * cos(elev*deg2rad) * sin(az*deg2rad)
  vr = vr + v * cos(elev*deg2rad) * cos(az*deg2rad)
  IF( USE_TERMINAL_VELOCITY )THEN
    vr = vr + (w - wt)*sin(elev*deg2rad)
  ELSE
    vr = vr + (w)*sin(elev*deg2rad)
  ENDIF

  !WRITE(6,*) u , v , w
  !WRITE(6,*) wt , vr
  !WRITE(6,*) elev , az , deg2rad


  RETURN
END SUBROUTINE calc_ref_vr

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
SUBROUTINE get_nobs_radar(cfile,nn,radarlon,radarlat,radarz)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn
!  REAL(r_sngl) :: wk(8)
  INTEGER :: nrec
  REAL(r_sngl) :: tmp
  INTEGER :: ios
!  INTEGER :: ir,iv
  INTEGER :: iunit
  LOGICAL :: ex
  INTEGER :: sz
  REAL(r_size),INTENT(OUT) :: radarlon,radarlat,radarz

  IF(RADAR_OBS_4D) THEN
    nrec = 8
  ELSE
    nrec = 7
  END IF
  nn = 0
!  iv = 0
!  ir = 0
  iunit=91

  radarlon=0.0d0
  radarlat=0.0d0
  radarz  =0.0d0

  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarlon=REAL(tmp,r_size)
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarlat=REAL(tmp,r_size)
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarz=REAL(tmp,r_size)

! get file size by reading through the entire file... too slow for big files
!-----------------------------
!    DO
!      READ(iunit,IOSTAT=ios) wk(1:nrec)
!      IF(ios /= 0) EXIT
!!      SELECT CASE(NINT(wk(1)))
!!      CASE(id_radar_ref_obs,id_radar_ref_zero_obs)
!!        ir = ir + 1
!!      CASE(id_radar_vr_obs)
!!        iv = iv + 1
!!      END SELECT
!      nn = nn + 1
!    END DO
!-----------------------------

! get file size by INQUIRE statement... may not work for some older fortran compilers
!-----------------------------
    INQUIRE(UNIT=iunit, SIZE=sz)
    sz = sz - r_sngl * (1+2) * 3 ! substract the radar data header
    IF (MOD(sz, r_sngl * (nrec+2)) /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    nn = sz / (r_sngl * (nrec+2))
!-----------------------------

    WRITE(6,*)' RADAR FILE ', cfile
    WRITE(6,*)' RADAR LON = ',radarlon
    WRITE(6,*)' RADAR LAT = ',radarlat
    WRITE(6,*)' RADAR Z   = ',radarz
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
!    WRITE(6,'(A12,I10)') '   REFLECTIVITY:',ir
!    WRITE(6,'(A12,I10)') 'RADIAL VELOCITY:',iv
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs_radar

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
SUBROUTINE read_obs_radar(cfile,obs)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(INOUT) :: obs
  REAL(r_sngl) :: wk(8)
  INTEGER :: nrec
  REAL(r_sngl) :: tmp
  INTEGER :: n,iunit,ios

  IF(RADAR_OBS_4D) THEN
    nrec = 8
  ELSE
    nrec = 7
  END IF
  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  DO n=1,obs%nobs
    READ(iunit) wk(1:nrec)
    obs%elm(n) = NINT(wk(1))
    obs%lon(n) = REAL(wk(2),r_size)
    obs%lat(n) = REAL(wk(3),r_size)
    obs%lev(n) = REAL(wk(4),r_size)
    obs%dat(n) = REAL(wk(5),r_size)
    obs%err(n) = REAL(wk(6),r_size)
!    obs%typ(n) = NINT(wk(7))
    obs%typ(n) = 22
    IF(RADAR_OBS_4D) THEN
      obs%dif(n) = REAL(wk(8),r_size)
    ELSE
      obs%dif(n) = 0.0d0
    END IF
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_radar

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
SUBROUTINE write_obs_radar(cfile,obs,append,missing)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(IN) :: obs
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL,INTENT(IN),OPTIONAL :: missing
  LOGICAL :: append_
  LOGICAL :: missing_
  REAL(r_sngl) :: wk(8)
  INTEGER :: nrec
  INTEGER :: n,iunit

  IF(RADAR_OBS_4D) THEN
    nrec = 8
  ELSE
    nrec = 7
  END IF
  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  missing_ = .true.
  IF(present(missing)) missing_ = missing

  IF(append_) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  END IF
  WRITE(iunit) REAL(obs%meta(1),r_sngl)
  WRITE(iunit) REAL(obs%meta(2),r_sngl)
  WRITE(iunit) REAL(obs%meta(3),r_sngl)
  DO n=1,obs%nobs
    if (missing_ .or. abs(obs%dat(n) - undef) > tiny(wk(5))) then
      wk(1) = REAL(obs%elm(n),r_sngl)
      wk(2) = REAL(obs%lon(n),r_sngl)
      wk(3) = REAL(obs%lat(n),r_sngl)
      wk(4) = REAL(obs%lev(n),r_sngl)
      wk(5) = REAL(obs%dat(n),r_sngl)
      wk(6) = REAL(obs%err(n),r_sngl)
      wk(7) = REAL(obs%typ(n),r_sngl)
      IF(RADAR_OBS_4D) THEN
        wk(8) = REAL(obs%dif(n),r_sngl)
      END IF
      WRITE(iunit) wk(1:nrec)
    end if
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_radar

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
subroutine read_obs_radar_toshiba(cfile, obs)
  use iso_c_binding
  use read_toshiba_f
#ifdef JITDT
  use jitdt_read_toshiba_f
#endif
  use radar_tools
  use scale_atmos_grid_cartesC, only: &
      DX, DY
  use scale_time, only: &
      TIME_gettimelabel, &
      TIME_NOWDATE
#ifdef PLOT_DCL
  use common_mpi_scale, only: &
      myrank_o, &
      myrank_da, &
      pawr_toshiba_hd_mpi
#endif
  implicit none

  character(len=*), intent(in) :: cfile
  type(obs_info), intent(out) :: obs
!  REAL(r_sngl) :: wk(8)
!  INTEGER :: nrec
!  REAL(r_sngl) :: tmp
!  INTEGER :: n,iunit,ios

  integer, parameter :: n_type = 3
  character(len=1024) :: jitdt_place
  character(len=4), parameter :: file_type_sfx(n_type) = &
    (/'.ze', '.vr', '.qcf'/)
  logical, parameter :: input_is_dbz = .true.

  type(c_pawr_header) :: hd(n_type)
!  real(kind=c_float) :: az(AZDIM, ELDIM, n_type)
!  real(kind=c_float) :: el(AZDIM, ELDIM, n_type)
!  real(kind=c_float) :: rtdat(RDIM, AZDIM, ELDIM, n_type)
  real(kind=c_float), allocatable :: rtdat(:, :, :, :)
  real(kind=c_float), allocatable :: az(:, :, :)
  real(kind=c_float), allocatable :: el(:, :, :)
  integer j, ierr
  character(len=3) :: fname
  integer, save::i=0

  real(r_size), allocatable :: ze(:, :, :), vr(:, :, :), qcflag(:, :, :), attenuation(:, :, :), rrange(:)
  real(r_size), allocatable :: lon(:), lat(:), z(:)
  integer(8), allocatable :: grid_index(:), grid_count_ze(:), grid_count_vr(:)
  real(r_size), allocatable :: grid_ze(:), grid_vr(:)
  real(r_size), allocatable :: grid_lon_ze(:), grid_lat_ze(:), grid_z_ze(:)
  real(r_size), allocatable :: grid_lon_vr(:),  grid_lat_vr(:),  grid_z_vr(:)

  character(len=1024) :: input_fname(n_type)
  integer na, nr, ne, ia, ir, ie
  real(r_size) :: lon0, lat0, z0
  real(r_size) :: dlon, dlat
  real(r_size) :: missing
  integer :: nlon , nlat , nlev
  integer(8) nobs_sp


  real(r_size) :: max_obs_ze , min_obs_ze , max_obs_vr , min_obs_vr 
  integer :: nobs_ze, nobs_vr
  integer(8) :: idx, n
  integer :: pos
  integer, parameter :: int1 = selected_int_kind(1) !1-BYTE INT
  integer(kind = int1) :: tmp_qcf, valid_qcf

!  integer,parameter :: qcf_mask(8)=(/ 1, 0, 0, 0, 0, 0, 0, 0 /) !! valid, shadow, clutter possible, clutter certain, interference, range sidelobe /
  integer,parameter :: qcf_mask(8)=(/ 0, 1, 1, 1, 1, 0, 0, 0 /) !! valid, shadow, clutter possible, clutter certain, interference, range sidelobe /
!!!  integer,parameter :: qcf_mask(8)=(/ 0, 0, 0, 0, 0, 0, 0, 0 /) !! valid, shadow, clutter possible, clutter certain, interference, range sidelobe /

  integer::qcf_count(0:255)

#ifdef PLOT_DCL
  character(len=8)  :: date
  character(len=10) :: time
  character(len=90) :: plotname
  character(len=19) :: timelabel
#endif

  integer :: range_res

!  real(kind=r_sngl), allocatable :: rtdat_l(:,:,:,:)
!  real(kind=r_sngl), allocatable :: az_l(:,:,:), el_l(:,:,:)
!  real(r_size), allocatable :: ze_l(:, :, :), vr_l(:, :, :), qcflag_l(:, :, :), attenuation_l(:, :, :)
!  real(r_size), allocatable :: radlon_l(:, :, :), radlat_l(:, :, :), radz_l(:, :, :)
!  integer :: ne_lmax

  call mpi_timer('', 3)

  if ( OBS_JITDT_CHECK_RADAR_TIME .and. minval(utime_obs) >= 0 ) then
    if ( .not. obs_da_same_time(utime_obs) ) then
      if (myrank_o == 0 ) then
        write(6,'(a)') "Obs & analysis times are different!"
        write(6,'(a,i4.4,i2.2,i2.2,1x,i2.2,1a,i2.2,1a,i2.2)') "SCALE-LETKF:",TIME_NOWDATE(1),TIME_NOWDATE(2),TIME_NOWDATE(3),&
                                                                     TIME_NOWDATE(4),":",TIME_NOWDATE(5),":",TIME_NOWDATE(6)
        write(6,'(a,i4.4,i2.2,i2.2,1x,i2.2,1a,i2.2,1a,i2.2)') "PAWR OBS (previous):",utime_obs(1),utime_obs(2),utime_obs(3),&
                                                                             utime_obs(4),":",utime_obs(5),":",utime_obs(6)
        obs%nobs = 0
      endif

      return
    endif

  endif

  allocate(rtdat(RDIM, AZDIM, ELDIM, n_type))
  allocate(az(AZDIM, ELDIM, n_type))
  allocate(el(AZDIM, ELDIM, n_type))

  RADAR_SO_SIZE_HORI = max(real(DX,kind=r_size),RADAR_SO_SIZE_HORI)
  RADAR_SO_SIZE_HORI = max(real(DY,kind=r_size),RADAR_SO_SIZE_HORI)

  if (LOG_LEVEL >= 3 .and. myrank_o == 0) then
    write(*, *) RDIM, AZDIM, ELDIM
    write(*, *) "dx = ", RADAR_SO_SIZE_HORI
    write(*, *) "dy = ", RADAR_SO_SIZE_HORI
    write(*, *) "dz = ", RADAR_SO_SIZE_VERT
  endif

  if (myrank_o == 0) then
#ifdef JITDT
    if (OBS_USE_JITDT) then
  !    jitdt_place = trim(OBS_JITDT_DATADIR) !// '/'
      jitdt_place = trim(OBS_JITDT_IP)
      write(*, *) "jitdt_place = ", trim(jitdt_place)
  
      ierr = jitdt_read_toshiba(n_type, jitdt_place, hd, az, el, rtdat)
  
      call mpi_timer('read_obs_radar_toshiba:jitdt_read_toshiba:', 2)
      if (ierr /= 0) then
        obs%nobs = 0
        return
      endif
    else
#endif
      do j = 1, n_type
        input_fname(j) = trim(cfile)
        call str_replace(input_fname(j), '<type>', trim(file_type_sfx(j)), pos)
        if (pos == 0) then
          write (6, '(5A)') "[Error] Keyword '<type>' is not found in '", trim(cfile), "'."
          stop 1
        end if
      end do
  
      if (LOG_LEVEL >= 3) then
        write(*, *) "file1 = ", trim(input_fname(1))
        write(*, *) "file2 = ", trim(input_fname(2))
        write(*, *) "file3 = ", trim(input_fname(3))
      endif
  
      do j = 1, n_type
        ierr = read_toshiba(input_fname(j), hd(j), az(:, :, j), el(:, :, j), rtdat(:, :, :, j))
        if (LOG_LEVEL >= 3) then
          write(*, *) "return code = ", ierr
        endif
        if (ierr /= 0) then
          obs%nobs = 0
          return
        endif
      end do
  
#ifdef JITDT
    end if
#endif
  end if  ! myrank_o == 0
  call mpi_timer('read_obs_radar_toshiba:read_toshiba:', 2, barrier=MPI_COMM_o)

  ! Set obs information
  if (myrank_o == 0) then
    lon0 = hd(1)%longitude
    lat0 = hd(1)%latitude
    z0 = hd(1)%altitude
    missing = real(hd(1)%mesh_offset, r_size)
    range_res = hd(1)%range_res

    na = hd(1)%sector_num
    nr = hd(1)%range_num
    ne = hd(1)%el_num

    call jst2utc(hd(1)%s_yr, hd(1)%s_mn, hd(1)%s_dy, hd(1)%s_hr, hd(1)%s_mi, hd(1)%s_sc, 0.0_DP, utime_obs)

  endif

  call  pawr_toshiba_hd_mpi(lon0, lat0, z0, missing,&
                           range_res, na, nr, ne, &
                           AZDIM, ELDIM, n_type, RDIM, &
                           az, el, rtdat, utime_obs)

!  ! Set local index for elavation (ne_lmax)
!  ie = mod(ELDIM, nprocs_o) 
!  ne_lmax = ((ELDIM - ie ) / nprocs_o) + 1

!  allocate(rtdat_l(1:RDIM, 1:AZDIM, ne_lmax, n_type))

!  call pawr_toshiba_scattv_mpi(RDIM, AZDIM, ELDIM, n_type, ne_lmax, &
!                              rtdat, rtdat_l)

  call mpi_timer('read_obs_radar_toshiba:comm:', 2, barrier=MPI_COMM_o)

  if (LOG_LEVEL >= 2 .and. myrank_o == 0) then
    write(*, '(I4.4, "-", I2.2, "-", I2.2, "T", I2.2, ":", I2.2, ":", I2.2, &
         &     " -> ", I4.4, "-", I2.2, "-", I2.2, "T", I2.2, ":", I2.2, ":", I2.2)') &
         & hd(1)%s_yr, hd(1)%s_mn, hd(1)%s_dy, hd(1)%s_hr, hd(1)%s_mi, hd(1)%s_sc, &
         & hd(1)%e_yr, hd(1)%e_mn, hd(1)%e_dy, hd(1)%e_hr, hd(1)%e_mi, hd(1)%e_sc
    write(*, *) lon0, lat0, z0
    write(*, *) hd(1)%range_num, hd(1)%sector_num, hd(1)%el_num
    write(*, *) "missing = ", missing
  endif

  if ( OBS_JITDT_CHECK_RADAR_TIME ) then
    if ( .not. obs_da_same_time(utime_obs) ) then
      if (myrank_o == 0 ) then
        write(*,*) "Obs & analysis times are different!"
        write(6,'(a,i4.4,i2.2,i2.2,1x,i2.2,1a,i2.2,1a,i2.2)') "SCALE-LETKF:",TIME_NOWDATE(1),TIME_NOWDATE(2),TIME_NOWDATE(3),&
                                                                     TIME_NOWDATE(4),":",TIME_NOWDATE(5),":",TIME_NOWDATE(6)
        write(6,'(a,i4.4,i2.2,i2.2,1x,i2.2,1a,i2.2,1a,i2.2)') "PAWR OBS:",utime_obs(1),utime_obs(2),utime_obs(3),&
                                                                          utime_obs(4),":",utime_obs(5),":",utime_obs(6)
        obs%nobs = 0
      endif

      return
    endif
  endif

!  i = 1
!! OUTPUT SPHERICAL COORDINATE DATA FOR DEBUG !!!
!  write(fname, '(I03.3)') i
!  open(1, file = trim(fname) // ".bin", access = "stream", form = "unformatted")
!  write(1) rtdat(1:hd(1)%range_num, 1:hd(1)%sector_num, 1:hd(1)%el_num, :)
!  close(1)

!  i=0

  allocate(ze(na, nr, ne), vr(na, nr, ne), qcflag(na, nr, ne), attenuation(na, nr, ne))

!  allocate(ze_l(na, nr, ne_lmax), vr_l(na, nr, ne_lmax))
!  allocate(qcflag_l(na, nr, ne_lmax), attenuation_l(na, nr, ne_lmax))

  valid_qcf = 0
  do j = 1, 8  
    if(qcf_mask(j) > 0) valid_qcf = ibset(valid_qcf, j - 1) 
  end do

!
!qcf_count=0

!$omp parallel do private(ia, ir, ie)
  do ie = 1, ne
     do ir = 1, nr
        do ia = 1, na
           ze(ia, ir, ie) = rtdat(ir, ia, ie, 1)
           vr(ia, ir, ie) = rtdat(ir, ia, ie, 2)                                                                      
           tmp_qcf = int(rtdat(ir, ia, ie, 3), int1)

 !          qcf_count(tmp_qcf)=qcf_count(tmp_qcf)+1
           
           if(iand(valid_qcf, tmp_qcf) == 0) then      
              qcflag(ia, ir, ie) = 0.0d0 !valid
           else
              qcflag(ia, ir, ie) = 1000.0d0 !invalid
           end if
           if(vr(ia, ir, ie) > RADAR_MAX_ABS_VR .or. vr(ia, ir, ie) < -RADAR_MAX_ABS_VR) vr(ia, ir, ie) = missing
           attenuation(ia, ir, ie) = 1.0d0 !not implemented yet
        end do
     end do
  end do
!$omp end parallel do

  deallocate(rtdat)

!do j=0,255
! write(*,*) j,qcf_count(j)
!end do
!stop

  call mpi_timer('read_obs_radar_toshiba:preliminary_qc:', 2, barrier=MPI_COMM_o)

  allocate(rrange(nr))
!$omp parallel do private(ir)
  do ir = 1, nr
     rrange(ir) = (dble(ir) - 0.5d0) * range_res
  end do
!$omp end parallel do

  if ((.not. allocated(radlon)) .and. (.not. allocated(radlat)) .and. (.not. allocated(radz)) ) then 
    allocate(radlon(na, nr, ne), radlat(na, nr, ne), radz(na, nr, ne))
    call radar_georeference(lon0, lat0, z0, na, nr, ne, &                                   ! input
         &                  real(az(:, 1, 1), r_size), rrange, real(el(1, :, 1), r_size), & ! input (assume ordinary scan strategy)
         &                  radlon, radlat, radz)                                     ! output

    call mpi_timer('read_obs_radar_toshiba:radar_georeference:', 2, barrier=MPI_COMM_o)
  endif

  deallocate(az, el)

!  write(*, *) "call define_grid"
  call define_grid(lon0, lat0, nr, rrange, rrange(nr), RADAR_ZMAX, RADAR_SO_SIZE_HORI, RADAR_SO_SIZE_HORI, RADAR_SO_SIZE_VERT, & ! input
       &           dlon, dlat, nlon, nlat, nlev, lon, lat, z)              ! output

  call mpi_timer('read_obs_radar_toshiba:define_grid:', 2, barrier=MPI_COMM_o)

!  write(*, *) "call radar_superobing"
  call radar_superobing(na, nr, ne, radlon, radlat, radz, ze, vr, &    ! input spherical
       &                qcflag, attenuation, &                                         ! input spherical
       &                nlon, nlat, nlev, lon, lat, z, dlon, dlat, RADAR_SO_SIZE_VERT, & ! input cartesian
       &                missing, input_is_dbz, &                                       ! input param
       &                lon0, lat0, &
       &                nobs_sp, grid_index, &                                         ! output array info
       &                grid_ze, grid_lon_ze, grid_lat_ze, grid_z_ze, grid_count_ze, & ! output ze
       &                grid_vr, grid_lon_vr, grid_lat_vr, grid_z_vr, grid_count_vr)   ! output vr
!  write(*, *) "done"

  if(allocated(ze)) deallocate(ze)
  if(allocated(vr)) deallocate(vr)
  if(allocated(qcflag)) deallocate(qcflag)
  if(allocated(attenuation)) deallocate(attenuation)
  if(allocated(rrange)) deallocate(rrange)


  call mpi_timer('read_obs_radar_toshiba:radar_superobing:', 2, barrier=MPI_COMM_o)

#ifdef PLOT_DCL
  if (PLOT_OBS)then
!    call date_and_time(date=date, time=time)
!    write (6, '(2a,1x,a,1x,a)') '[Info:plot] obs start plotting: ', date, time, trim(timelabel(1:15))

    call TIME_gettimelabel(timelabel)
    plotname = "obs_dbz_"//trim(timelabel(1:15))
    call plot_dbz_DCL_obs (nobs_sp,real(grid_ze),real(grid_lon_ze),real(grid_lat_ze),real(grid_z_ze),trim(plotname))

!    call date_and_time(date=date, time=time)
!    write (6, '(2a,1x,a,1x,a)') '[Info:plot] obs finish plotting: ', date, time, trim(timelabel(1:15))
  endif
  call mpi_timer('read_obs_radar_toshiba:plot_obs:', 2, barrier=MPI_COMM_o)
#endif

  if (myrank_o /= 0) return 


!!!!! check
!write(*,*) nlon,nlat,nlev,nobs_sp




  obs%meta(1) = lon0
  obs%meta(2) = lat0
  obs%meta(3) = z0

  obs%nobs = 0
  do idx = 1, nobs_sp
    if (grid_count_ze(idx) > 0) then
      obs%nobs = obs%nobs + 1
    end if
    if (grid_count_vr(idx) > 0) then
      obs%nobs = obs%nobs + 1
    end if
  end do
  call obs_info_allocate(obs, extended=.true.)

  n = 0
  nobs_ze = 0
  nobs_vr = 0
  min_obs_ze = huge(1.0d0)
  max_obs_ze = -huge(1.0d0)
  min_obs_vr = huge(1.0d0)
  max_obs_vr = -huge(1.0d0)
  do idx = 1, nobs_sp
    if (grid_count_ze(idx) > 0) then
      n = n + 1
      obs%elm(n) = id_radar_ref_obs
      obs%lon(n) = grid_lon_ze(idx)
      obs%lat(n) = grid_lat_ze(idx)
      obs%lev(n) = grid_z_ze(idx)
      obs%dat(n) = grid_ze(idx)
      obs%err(n) = OBSERR_RADAR_REF
      obs%typ(n) = 22
      obs%dif(n) = 0.0d0
      nobs_ze = nobs_ze + 1
      if (grid_ze(idx) > max_obs_ze) max_obs_ze = grid_ze(idx)
      if (grid_ze(idx) < min_obs_ze) min_obs_ze = grid_ze(idx)
    end if

    if (grid_count_vr(idx) > 0) then
      n = n + 1
      obs%elm(n) = id_radar_vr_obs
      obs%lon(n) = grid_lon_ze(idx)
      obs%lat(n) = grid_lat_ze(idx)
      obs%lev(n) = grid_z_ze(idx)
      obs%dat(n) = grid_vr(idx)
      obs%err(n) = OBSERR_RADAR_VR
      obs%typ(n) = 22
      obs%dif(n) = 0.0d0
      nobs_vr = nobs_vr + 1
      if (grid_vr(idx) > max_obs_vr) max_obs_vr = grid_vr(idx)
      if (grid_vr(idx) < min_obs_vr) min_obs_vr = grid_vr(idx)
    end if
  end do

  if (LOG_LEVEL >= 2) then
    write (6, *) "Reflectivity obs. range = ", min_obs_ze, " to ", max_obs_ze
    write (6, *) "Radial vel. obs. range  = ", min_obs_vr, " to ", max_obs_vr
    write (6, *) "ze: ", nobs_ze, ", vr: ", nobs_vr
  endif

  call mpi_timer('read_obs_radar_toshiba:save_obs_info:', 2)




!  write(*, *) "writing LETKF data ..."
!  call output_letkf_obs(lon0, lat0, z0, 1, nobs_sp, &
!       &                grid_ze, grid_lon_ze, grid_lat_ze, grid_z_ze, grid_count_ze, &
!       &                grid_vr, grid_count_vr, &
!       &                OBSERR_RADAR_REF, OBSERR_RADAR_VR)
!  write(*, *) "done"

!  !!! OUTPUT CARTESIAN COORDINATE DATE FOR DEBUG !!!
!if (nobs_ze.gt.0) then 
!  write(*, *) "writing cartesian data ..."
!  i=i+1
!  write(fname, '(I03.3)') i
!  call output_grads_obs("super_" // fname // ".grd", nlon, nlat, nlev, nobs_sp, grid_index, &
!       &                grid_ze, grid_lon_ze, grid_lat_ze, grid_z_ze, grid_count_ze, &
!       &                grid_vr, grid_lon_vr, grid_lat_vr, grid_z_vr, grid_count_vr, missing)
!  write(*, *) "done"
! endif



! Assume radar geographical information does not change during DA cyclying
!  if(allocated(radlon)) deallocate(radlon)
!  if(allocated(radlat)) deallocate(radlat)
!  if(allocated(radz)) deallocate(radz)

  if(allocated(grid_index)) deallocate(grid_index)
  if(allocated(grid_ze)) deallocate(grid_ze)
  if(allocated(grid_lon_ze)) deallocate(grid_lon_ze)
  if(allocated(grid_lat_ze)) deallocate(grid_lat_ze)
  if(allocated(grid_z_ze)) deallocate(grid_z_ze)
  if(allocated(grid_count_ze)) deallocate(grid_count_ze)
  if(allocated(grid_vr)) deallocate(grid_vr)
  if(allocated(grid_lon_vr)) deallocate(grid_lon_vr)
  if(allocated(grid_lat_vr)) deallocate(grid_lat_vr)
  if(allocated(grid_z_vr)) deallocate(grid_z_vr)
  if(allocated(grid_count_vr)) deallocate(grid_count_vr)

  call mpi_timer('read_obs_radar_toshiba:deallocate_vars:', 2)

  return
end subroutine read_obs_radar_toshiba

subroutine read_obs_radar_jrc(cfile, obs)
  use netcdf
  use common_ncio
  implicit none

  character(*), intent(in) :: cfile
  type(obs_info), intent(inout) :: obs
  integer :: ncid, varid
  integer :: status
  integer :: xdim, ydim, zdim
  integer :: i, j, k
  integer :: n, nobs_vr, nobs_zh

  real(r_sngl), allocatable :: vr3d(:,:,:), zh3d(:,:,:)
  real(r_sngl), allocatable :: lon1d(:), lat1d(:), z1d(:)
  real(r_sngl) :: sf_vr, sf_zh ! scale factor for vr & Zh
  real(r_sngl) :: fill_vr, fill_zh ! fill values for vr & Zh
  real(r_size) :: radar_lon, radar_lat
#ifdef SINGLELETKF
  real(8) :: radar_lon_r8, radar_lat_r8
#endif

  real(r_sngl) :: max_obs_zh, min_obs_zh
  real(r_sngl) :: max_obs_vr, min_obs_vr

  call mpi_timer('', 3)

  status = nf90_open(trim(cfile), NF90_NOWRITE, ncid)
  if (status /= nf90_noerr) then
    obs%nobs = 0
    return
  endif

  call ncio_read_dim(ncid, "lon", xdim)
  call ncio_read_dim(ncid, "lat", ydim)
  call ncio_read_dim(ncid, "alt", zdim)

  allocate(vr3d(xdim,ydim,zdim), zh3d(xdim,ydim,zdim))
  allocate(lon1d(xdim), lat1d(ydim), z1d(zdim)) 

  ! get lon/lat/height information
  call ncio_check(nf90_inq_varid(ncid, "lon",  varid))
  call ncio_check(nf90_get_var(ncid, varid, lon1d(:), &
                               start = (/ 1/),    &
                               count = (/ xdim/)))

  call ncio_check(nf90_inq_varid(ncid, "lat",  varid))
  call ncio_check(nf90_get_var(ncid, varid, lat1d(:), &
                               start = (/ 1/),    &
                               count = (/ ydim/)))

  call ncio_check(nf90_inq_varid(ncid, "alt",  varid))
  call ncio_check(nf90_get_var(ncid, varid, z1d(:), &
                               start = (/ 1/),    &
                               count = (/ zdim/)))



  call ncio_check(nf90_inq_varid(ncid, "Vel",  varid))
  call ncio_check(nf90_get_var(ncid, varid, vr3d(:,:,:), &
                               start = (/ 1, 1, 1/),    &
                               count = (/ xdim, ydim, zdim /)))

  status = nf90_get_att(ncid,varid,"scale_factor",sf_vr)
  if (status /= nf90_noerr) then
    sf_vr = 1.0
  endif
  status = nf90_get_att(ncid,varid,"_FillValue",fill_vr)
  if (status /= nf90_noerr) then
    fill_vr = -32768.0
  endif

  call ncio_check(nf90_inq_varid(ncid, "Zh MTI",  varid))
  call ncio_check(nf90_get_var(ncid, varid, zh3d(:,:,:), &
                               start = (/ 1, 1, 1/),    &
                               count = (/ xdim, ydim, zdim /)))

  status = nf90_get_att(ncid,varid,"scale_factor",sf_zh)
  if (status /= nf90_noerr) then
    sf_zh = 1.0
  endif
  status = nf90_get_att(ncid,varid,"_FillValue",fill_zh)
  if (status /= nf90_noerr) then
    fill_zh = -32768.0
  endif

#ifdef SINGLELETKF
  call ncio_read_gattr_r8(ncid, "site_positions_center_latitude",  radar_lon_r8)
  call ncio_read_gattr_r8(ncid, "site_positions_center_longitude", radar_lat_r8)
 radar_lon=real(radar_lon_r8)
 radar_lat=real(radar_lat_r8)
#else
  call ncio_read_gattr_r8(ncid, "site_positions_center_latitude",  radar_lon)
  call ncio_read_gattr_r8(ncid, "site_positions_center_longitude", radar_lat)
#endif

  call ncio_close(ncid)


  obs%nobs = 0

  obs%meta(1) = real(radar_lon, kind=r_size)
  obs%meta(2) = real(radar_lat, kind=r_size)
  obs%meta(3) = 82.0_r_size

  ! count number of obs 
  do k = 1, zdim
    do j = 1, ydim
      do i = 1, xdim
        if (zh3d(i,j,k) /= fill_zh .and. vr3d(i,j,k) /= fill_vr)then
          obs%nobs = obs%nobs + 2
        endif
      enddo ! i
    enddo ! j
  enddo ! k

  call obs_info_allocate(obs, extended=.true.)

  ! substitute obs info
  n = 0
  nobs_vr = 0
  nobs_zh = 0
  min_obs_zh = huge(1.0)
  max_obs_zh = -huge(1.0)
  min_obs_vr = huge(1.0)
  max_obs_vr = -huge(1.0)

  do k = 1, zdim
    do j = 1, ydim
      do i = 1, xdim
        if (zh3d(i,j,k) /= fill_zh .and. vr3d(i,j,k) /= fill_vr)then
          ! zh
          n = n + 1
          obs%elm(n) = id_radar_ref_obs
          obs%lon(n) = real(lon1d(i), kind=r_size)
          obs%lat(n) = real(lat1d(j), kind=r_size)
          obs%lev(n) = real(z1d(k), kind=r_size)
          obs%dat(n) = real(zh3d(i,j,k) * sf_zh, kind=r_size)
          obs%err(n) = OBSERR_RADAR_REF
          obs%typ(n) = 22
          obs%dif(n) = 0.0_r_size
          nobs_zh = nobs_zh + 1

          if (zh3d(i,j,k) > max_obs_zh) max_obs_zh = zh3d(i,j,k)
          if (zh3d(i,j,k) < min_obs_zh) min_obs_zh = zh3d(i,j,k)

          ! vr
          n = n + 1
          obs%elm(n) = id_radar_vr_obs
          obs%lon(n) = real(lon1d(i), kind=r_size)
          obs%lat(n) = real(lat1d(j), kind=r_size)
          obs%lev(n) = real(z1d(k), kind=r_size)
          obs%dat(n) = real(vr3d(i,j,k) * sf_vr, kind=r_size)
          obs%err(n) = OBSERR_RADAR_VR
          obs%typ(n) = 22
          obs%dif(n) = 0.0_r_size
          nobs_vr = nobs_vr + 1

          if (vr3d(i,j,k) > max_obs_vr) max_obs_vr = vr3d(i,j,k)
          if (vr3d(i,j,k) < min_obs_vr) min_obs_vr = vr3d(i,j,k)
        endif
      enddo ! i
    enddo ! j
  enddo ! k

  write (6, *) "Reflectivity obs. range = ", min_obs_zh * sf_zh, &
               " to ", max_obs_zh * sf_zh
  write (6, *) "Radial vel. obs. range  = ", min_obs_vr * sf_vr, &
               " to ", max_obs_vr * sf_vr
  write (6, *) "zh: ", nobs_zh, ", vr: ", nobs_vr

  deallocate(vr3d, zh3d)
  deallocate(lon1d, lat1d, z1d)

  call mpi_timer('read_obs_radar_jrc:read_obs_radar_jrc:', 2)

  return
end subroutine read_obs_radar_jrc

function obs_da_same_time(utime_obs)
  use scale_time, only: &
      TIME_NOWDATE
  implicit none

  integer :: utime_obs(6)
  integer :: i
  logical :: obs_da_same_time

  obs_da_same_time = .true.
  do i = 1, 6
    if (utime_obs(i) /= TIME_NOWDATE(i)) then
      obs_da_same_time = .false.
      exit
    endif
  enddo

end function obs_da_same_time

!=======================================================================
end module radar_obs
