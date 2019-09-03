MODULE common_obs_scale
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   10/04/2012 Guo-Yuan Lien     modified for GFS model
!   07/25/2014 Guo-Yuan Lien     modified for SCALE model
!   .......... See git history for the following revisions
!
!=======================================================================
!
! [LETKF observation format]
!   (In files, all variables are stored in single-precision float)
!
!  column  description
!     (1)  variable type (1..nid_obs; see 'id_*_obs' parameters)
!     (2)  longitude (degree)
!     (3)  latitude (degree)
!     (4)  level/height
!            u,v,t,tv,q,rh: level (hPa)
!            ps: station elevation (m)
!     (5)  observation value
!            wind (m/s)
!            temperature (K)
!            specific humidity (kg/kg)
!            relative humidity (%)
!            surface pressure (hPa)
!     (6)  observation error
!            unit same as observation value
!     (7)  observation platform type (1..nobtype+1; see 'obtypelist' array)
!     (8)  observation time relative to analysis time (sec)
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_nml
  USE common_scale
  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: nid_obs_varlocal=9 !Him8
!
! conventional observations
!
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_tv_obs=3074
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
!
! surface observations codes > 9999
!
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

  INTEGER,PARAMETER :: elem_uid(nid_obs)= &
     (/id_u_obs, id_v_obs, id_t_obs, id_tv_obs, id_q_obs, id_rh_obs, &
       id_ps_obs, id_rain_obs, id_radar_ref_obs, id_radar_ref_zero_obs, id_radar_vr_obs, id_radar_prh_obs, &
       id_H08IR_obs, id_tclon_obs, id_tclat_obs, id_tcmip_obs/)

  CHARACTER(3),PARAMETER :: obelmlist(nid_obs)= &
     (/'  U', '  V', '  T', ' Tv', '  Q', ' RH', ' PS', 'PRC', 'REF', 'RE0', ' Vr', 'PRH',&
       'H08', 'TCX', 'TCY', 'TCP'/)

  CHARACTER(3),PARAMETER :: obelmlist_varlocal(nid_obs_varlocal)= &
     (/'WND', '  T', 'MOI', ' PS', 'PRC', 'TCV', 'REF', ' Vr', 'H08'/)

  ! Parameter 'nobtype' is set in common_nml.f90
  CHARACTER(6),PARAMETER :: obtypelist(nobtype)= &
     (/'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR', &
       'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG', &
       'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND', &
       'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW', &
       'TMPAPR', 'PHARAD', 'H08IRB', 'TCVITL'/) ! Him8

  INTEGER,PARAMETER :: max_obs_info_meta = 3 ! maximum array size for type(obs_info)%meta

  TYPE obs_info
    INTEGER :: nobs = 0
    INTEGER,ALLOCATABLE :: elm(:)
    REAL(r_size),ALLOCATABLE :: lon(:)
    REAL(r_size),ALLOCATABLE :: lat(:)
    REAL(r_size),ALLOCATABLE :: lev(:)
    REAL(r_size),ALLOCATABLE :: dat(:)
    REAL(r_size),ALLOCATABLE :: err(:)
    INTEGER,ALLOCATABLE :: typ(:)
    REAL(r_size),ALLOCATABLE :: dif(:)
    REAL(r_size) :: meta(max_obs_info_meta) = undef
    REAL(r_size),ALLOCATABLE :: ri(:)
    REAL(r_size),ALLOCATABLE :: rj(:)
    INTEGER,ALLOCATABLE :: rank(:)
  END TYPE obs_info

  TYPE obs_da_value
    INTEGER :: nobs = 0
    INTEGER :: nobs_in_key = 0
    INTEGER,ALLOCATABLE :: set(:)
    INTEGER,ALLOCATABLE :: idx(:)
    INTEGER,ALLOCATABLE :: key(:)
    REAL(r_size),ALLOCATABLE :: val(:)
    !
    ! obsda%lev array is used only for Himawari-8 assimilation.
    ! This array preserves the most sensitive height derived from transmittance outputs from RTTOV.
    ! For Himawari-8 assimilation, LETKF uses obsda%lev instead of obs%lev.
    ! 
#ifdef H08
    REAL(r_size),ALLOCATABLE :: lev(:) ! Him8
    REAL(r_size),ALLOCATABLE :: val2(:) ! Him8 sigma_o for AOEI (not CA)
    REAL(r_size),ALLOCATABLE :: sprd(:) ! background spread
!    REAL(r_size),ALLOCATABLE :: pred1(:) ! Him8 bias correction predictor 1 (nobs)
!    REAL(r_size),ALLOCATABLE :: pred2(:) ! Him8 bias correction predictor 1 (nobs)
#endif
    REAL(r_size),ALLOCATABLE :: ensval(:,:)
    INTEGER,ALLOCATABLE :: qc(:)
  END TYPE obs_da_value

  character(obsformatlenmax), parameter :: obsfmt_prepbufr = 'PREPBUFR'
  character(obsformatlenmax), parameter :: obsfmt_radar    = 'RADAR'
  character(obsformatlenmax), parameter :: obsfmt_h08      = 'HIMAWARI8'
!  integer, parameter :: nobsformats = 3
!  character(obsformatlenmax), parameter :: obsformat(nobsformats) = &
!    (/obsfmt_prepbufr, obsfmt_radar, obsfmt_h08/)

  INTEGER,PARAMETER :: iqc_good=0
  INTEGER,PARAMETER :: iqc_gross_err=5
  INTEGER,PARAMETER :: iqc_ps_ter=10
  INTEGER,PARAMETER :: iqc_ref_low=11
  INTEGER,PARAMETER :: iqc_ref_mem=12
  INTEGER,PARAMETER :: iqc_radar_vhi=19
  INTEGER,PARAMETER :: iqc_out_vhi=20
  INTEGER,PARAMETER :: iqc_out_vlo=21
  INTEGER,PARAMETER :: iqc_obs_bad=50
  INTEGER,PARAMETER :: iqc_otype=90
  INTEGER,PARAMETER :: iqc_time=97
  INTEGER,PARAMETER :: iqc_out_h=98
  INTEGER,PARAMETER :: iqc_undef=99

  type(obs_info),allocatable,save :: obs(:) ! observation information
  type(obs_da_value),save :: obsda_sort     ! sorted obsda

  integer, save :: obsdep_nobs                     ! obsdep information
  integer, allocatable, save :: obsdep_set(:)      ! 
  integer, allocatable, save :: obsdep_idx(:)      ! 
  integer, allocatable, save :: obsdep_qc(:)       ! 
  real(r_size), allocatable, save :: obsdep_omb(:) ! 
  real(r_size), allocatable, save :: obsdep_oma(:) ! 

  REAL(r_size),SAVE :: MIN_RADAR_REF
  REAL(r_size),SAVE :: RADAR_REF_THRES

CONTAINS

!-----------------------------------------------------------------------
! Convert a raw obsID to a sequential obsID (1 - nid_obs)
!-----------------------------------------------------------------------
function uid_obs(id_obs)
  implicit none
  integer :: id_obs
  integer :: uid_obs

  select case(id_obs)
  case(id_u_obs)
    uid_obs = 1
  case(id_v_obs)
    uid_obs = 2
  case(id_t_obs)
    uid_obs = 3
  case(id_tv_obs)
    uid_obs = 4
  case(id_q_obs)
    uid_obs = 5
  case(id_rh_obs)
    uid_obs = 6
  case(id_ps_obs)
    uid_obs = 7
  case(id_rain_obs)
    uid_obs = 8
  case(id_radar_ref_obs)
    uid_obs = 9
  case(id_radar_ref_zero_obs)
    uid_obs = 10
  case(id_radar_vr_obs)
    uid_obs = 11
  case(id_radar_prh_obs)
    uid_obs = 12
  case(id_h08ir_obs) ! H08
    uid_obs = 13     ! H08
  case(id_tclon_obs)
    uid_obs = 14
  case(id_tclat_obs)
    uid_obs = 15
  case(id_tcmip_obs)
    uid_obs = 16
  case default
    uid_obs = -1     ! error
  end select
end function uid_obs
!-----------------------------------------------------------------------
! Convert a raw obsID to a sequential obsID for variable localization (1 - nid_obs_verlocal)
!-----------------------------------------------------------------------
function uid_obs_varlocal(id_obs)
  implicit none
  integer :: id_obs
  integer :: uid_obs_varlocal

  select case(id_obs)
  case(id_u_obs, id_v_obs)
    uid_obs_varlocal = 1
  case(id_t_obs, id_tv_obs)
    uid_obs_varlocal = 2
  case(id_q_obs, id_rh_obs)
    uid_obs_varlocal = 3
  case(id_ps_obs)
    uid_obs_varlocal = 4
  case(id_rain_obs)
    uid_obs_varlocal = 5
  case(id_tclon_obs, id_tclat_obs, id_tcmip_obs)
    uid_obs_varlocal = 6
  case(id_radar_ref_obs, id_radar_ref_zero_obs, id_radar_prh_obs)
    uid_obs_varlocal = 7
  case(id_radar_vr_obs)
    uid_obs_varlocal = 8
  case(id_h08ir_obs)      ! H08
    uid_obs_varlocal = 9  ! H08
  case default
    uid_obs_varlocal = -1 ! error
  end select
end function uid_obs_varlocal

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine set_common_obs_scale
  implicit none

  MIN_RADAR_REF = 10.0d0 ** (MIN_RADAR_REF_DBZ/10.0d0)
  RADAR_REF_THRES = 10.0d0 ** (RADAR_REF_THRES_DBZ/10.0d0)

  return
end subroutine set_common_obs_scale

!-----------------------------------------------------------------------
! Transformation from model variables to an observation
!
! stggrd: grid type of u and v
!  0: non-staggered grid
!  1: staggered grid
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY(elm,ri,rj,rk,lon,lat,v3d,v2d,yobs,qc,stggrd)
  use scale_mapprojection, only: &
      MAPPROJECTION_rotcoef
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_size),INTENT(IN) :: lon,lat
  REAL(r_size),INTENT(IN) :: v3d(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  REAL(r_size),INTENT(OUT) :: yobs
  INTEGER,INTENT(OUT) :: qc
  INTEGER,INTENT(IN),OPTIONAL :: stggrd
  REAL(r_size) :: u,v,t,q,topo
  REAL(RP) :: rotc(1,1,2)
  real(r_size) :: lon_tmp(1,1),lat_tmp(1,1)

  INTEGER :: stggrd_ = 0
  if (present(stggrd)) stggrd_ = stggrd

  yobs = undef
  qc = iqc_good

  SELECT CASE (elm)
  CASE(id_u_obs,id_v_obs)  ! U,V
    if (stggrd_ == 1) then
      CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri-0.5_r_size,rj,u)  !###### should modity itpl_3d to prevent '1.0' problem....??
      CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj-0.5_r_size,v)  !######
    else
      CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri,rj,u)
      CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj,v)
    end if
    lon_tmp(1,1) = lon*deg2rad
    lat_tmp(1,1) = lat*deg2rad
    call MAPPROJECTION_rotcoef(1, 1, 1, 1, 1, 1, &
                               lon_tmp(1,1),lat_tmp(1,1),rotc)
    if (elm == id_u_obs) then
      yobs = u * rotc(1,1,1) - v * rotc(1,1,2)
    else
      yobs = u * rotc(1,1,2) + v * rotc(1,1,1)
    end if
  CASE(id_t_obs)  ! T
    CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,yobs)
  CASE(id_tv_obs)  ! Tv
    CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,yobs)
    CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,q)
    yobs = yobs * (1.0d0 + fvirt * q)
  CASE(id_q_obs)  ! Q
    CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,yobs)
  CASE(id_ps_obs) ! PS
    CALL itpl_2d(v2d(:,:,iv2dd_t2m),ri,rj,t)
    CALL itpl_2d(v2d(:,:,iv2dd_q2m),ri,rj,q)
    CALL itpl_2d(v2d(:,:,iv2dd_topo),ri,rj,topo)
    CALL itpl_2d(v2d(:,:,iv2dd_ps),ri,rj,yobs)
    call prsadj(yobs,rk-topo,t,q)
    if (abs(rk-topo) > PS_ADJUST_THRES) then
      if (LOG_LEVEL >= 2) then
        write (6,'(A,F6.1)') '[Warning] PS observation height adjustment exceeds the threshold. dz=', abs(rk-topo)
      end if
      qc = iqc_ps_ter
    end if
!  CASE(id_rain_obs) ! RAIN                        ############# (not finished)
!    CALL itpl_2d(v2d(:,:,iv2dd_rain),ri,rj,yobs) !#############
  CASE(id_rh_obs) ! RH
    CALL itpl_3d(v3d(:,:,:,iv3dd_rh),rk,ri,rj,yobs)
!  CASE(id_tclon_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(1)
!  CASE(id_tclat_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(2)
!  CASE(id_tcmip_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(3)
  CASE DEFAULT
    qc = iqc_otype
  END SELECT

  RETURN
END SUBROUTINE Trans_XtoY
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY_radar(elm,radar_lon,radar_lat,radar_z,ri,rj,rk,lon,lat,lev,v3d,v2d,yobs,qc,stggrd)
  use scale_mapprojection, only: &
      MAPPROJECTION_rotcoef
!  USE common_mpi
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

  real(r_size) :: utmp, vtmp
  REAL(RP) :: rotc(1,1,2)
  real(RP) :: lon_tmp(1,1),lat_tmp(1,1)

!  integer :: ierr
!  REAL(r_dble) :: rrtimer00,rrtimer
!  rrtimer00 = MPI_WTIME()


  if (present(stggrd)) stggrd_ = stggrd


  yobs = undef
  qc = iqc_good

  if (stggrd_ == 1) then
    CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri-0.5_r_size,rj,utmp)  !###### should modity itpl_3d to prevent '1.0' problem....??
    CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj-0.5_r_size,vtmp)  !######
    CALL itpl_3d(v3d(:,:,:,iv3dd_w),rk-0.5_r_size,ri,rj,wr)  !######
  else
    CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri,rj,utmp)
    CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj,vtmp)
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
!

  lon_tmp(1,1) = lon*deg2rad
  lat_tmp(1,1) = lat*deg2rad
  call MAPPROJECTION_rotcoef(1, 1, 1, 1, 1, 1, &
                             lon_tmp(1,1),lat_tmp(1,1),rotc)

  ur = utmp * rotc(1,1,1) - vtmp * rotc(1,1,2)
  vr = utmp * rotc(1,1,2) + vtmp * rotc(1,1,1)


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

  !WRITE(6,*)'BCRV',dlon,dlat,az,eltmpv

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
! TC center search
!  [AUTHORS:] T. Miyoshi and M. Kunii
!-----------------------------------------------------------------------
!SUBROUTINE tctrk(ps,t2,ri,rj,trk)
!  IMPLICIT NONE
!  INTEGER,PARAMETER :: isearch = 10 !search radius [grid points]
!  REAL(r_size),INTENT(IN) :: ps(nlon,nlat)
!  REAL(r_size),INTENT(IN) :: t2(nlon,nlat)
!  REAL(r_size),INTENT(IN) :: ri,rj
!  REAL(r_size),INTENT(OUT) :: trk(3) !1:lon, 2:lat, 3:minp
!  REAL(r_size) :: wk(100,3)
!  REAL(r_size) :: slp(nlon,nlat)
!  REAL(r_size) :: p1,p2,p3,p4,p5,a,c,d,xx,yy
!  INTEGER :: i,j,i0,i1,j0,j1,n

!  i0 = MAX(1,FLOOR(ri)-isearch)
!  i1 = MIN(nlon,CEILING(ri)+isearch)
!  j0 = MAX(1,FLOOR(rj)-isearch)
!  j1 = MIN(nlat,CEILING(rj)+isearch)
!  trk = undef

!  DO j=j0,j1
!    DO i=i0,i1
!      slp(i,j) = ps(i,j) * (1.0d0 - 0.0065d0 * phi0(i,j) / &
!        & (t2(i,j) + 0.0065d0 * phi0(i,j))) ** -5.257d0
!    END DO
!  END DO

!  n=0
!  DO j=j0+1,j1-1
!    DO i=i0+1,i1-1
!      IF(slp(i,j) > slp(i  ,j-1)) CYCLE
!      IF(slp(i,j) > slp(i  ,j+1)) CYCLE
!      IF(slp(i,j) > slp(i-1,j  )) CYCLE
!      IF(slp(i,j) > slp(i+1,j  )) CYCLE
!      IF(slp(i,j) > slp(i-1,j-1)) CYCLE
!      IF(slp(i,j) > slp(i+1,j-1)) CYCLE
!      IF(slp(i,j) > slp(i-1,j+1)) CYCLE
!      IF(slp(i,j) > slp(i+1,j+1)) CYCLE
!      p1 = slp(i,j)
!      p2 = slp(i-1,j)
!      p3 = slp(i+1,j)
!      p4 = slp(i,j-1)
!      p5 = slp(i,j+1)
!      c = (p3-p2)*0.5d0
!      d = (p5-p4)*0.5d0
!      a = (p2+p3+p4+p5)*0.25d0 - p1
!      IF(a == 0.0d0) CYCLE
!      xx = -0.5d0 * c / a
!      yy = -0.5d0 * d / a
!      n = n+1
!      wk(n,3) = p1 - a*(xx*xx + yy*yy)
!      wk(n,2) = lat(i,j) * (1.0d0 - yy) + lat(i,j+1) * yy
!      wk(n,1) = lon(i,j) * (1.0d0 - xx) + lon(i+1,j) * xx
!    END DO
!  END DO

!  j=1
!  IF(n > 1) THEN
!    a = wk(1,3)
!    DO i=2,n
!      IF(wk(i,3) < a) THEN
!        a = wk(i,3)
!        j = i
!      END IF
!    END DO
!  END IF
!  trk = wk(j,:)

!END SUBROUTINE tctrk
!-----------------------------------------------------------------------
! Compute relative humidity (RH)
!-----------------------------------------------------------------------
SUBROUTINE calc_rh(t,q,p,rh)
  IMPLICIT NONE
  REAL(r_size),PARAMETER :: t0=273.15d0
  REAL(r_size),PARAMETER :: e0c=6.11d0
  REAL(r_size),PARAMETER :: al=17.3d0
  REAL(r_size),PARAMETER :: bl=237.3d0
  REAL(r_size),PARAMETER :: e0i=6.1121d0
  REAL(r_size),PARAMETER :: ai=22.587d0
  REAL(r_size),PARAMETER :: bi=273.86d0
  REAL(r_size),INTENT(IN) :: t,q,p
  REAL(r_size),INTENT(OUT) :: rh
  REAL(r_size) :: e,es,tc

  e = q * p * 0.01d0 / (0.378d0 * q + 0.622d0)

  tc = t-t0
  IF(tc >= 0.0d0) THEN
    es = e0c * exp(al*tc/(bl+tc))
  ELSE IF(tc <= -15.d0) THEN
    es = e0i * exp(ai*tc/(bi+tc))
  ELSE
    es = e0c * exp(al*tc/(bl+tc)) * (15.0d0+tc)/15.0d0 &
       + e0i * exp(ai*tc/(bi+tc)) * (-tc) / 15.0d0
  END IF

  rh = e/es

  RETURN
END SUBROUTINE calc_rh
!-----------------------------------------------------------------------
! Pressure adjustment for a different height level
!-----------------------------------------------------------------------
SUBROUTINE prsadj(p,dz,t,q)
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: p
  REAL(r_size),INTENT(IN) :: dz ! height difference (target - original) [m]
  REAL(r_size),INTENT(IN) :: t  ! temperature [K] at original level
  REAL(r_size),INTENT(IN) :: q  ! humidity [kg/kg] at original level
  REAL(r_size),PARAMETER :: gamma=5.0d-3 ! lapse rate [K/m]
  REAL(r_size) :: tv

  IF(dz /= 0) THEN
    tv = t * (1.0d0 + 0.608d0 * q)
    p = p * ((-gamma*dz+tv)/tv)**(gg/(gamma*rd)) !tv is at original level
!    p = p * (tv/(tv+gamma*dz))**(gg/(gamma*rd)) !tv is at target level
  END IF

  RETURN
END SUBROUTINE prsadj



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
!!!    zg= 8.18d4 * ( ro * qgp * 1.0d3 )**1.50  !!! hail
    zg= 5.54d3 * ( ro * qgp * 1.0d3 )**1.70   !!! graupel (A. Amemiya 2019.5)
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
! Coordinate conversion (find rk in pressure levels)
!
! rk = 0.0d0  : surface observation
!-----------------------------------------------------------------------
SUBROUTINE phys2ijk(p_full,elem,ri,rj,rlev,rk,qc)
  use scale_atmos_grid_cartesC_index, only: &
      KHALO
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: p_full(nlevh,nlonh,nlath)
  INTEGER,INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels (for 3D variable only)
  REAL(r_size),INTENT(OUT) :: rk
  INTEGER,INTENT(OUT) :: qc
  REAL(r_size) :: ak
!  REAL(r_size) :: lnps(nlonh,nlath)
  REAL(r_size) :: lnps(nlevh,nlonh,nlath)
  REAL(r_size) :: plev(nlevh)
  REAL(r_size) :: ptmp
  INTEGER :: i,j,k, ii, jj, ks

  qc = iqc_good
!
! rlev -> rk
!
  if (ri < 1.0d0 .or. ri > nlonh .or. rj < 1.0d0 .or. rj > nlath) then
    if (LOG_LEVEL >= 1) then
      write (6,'(A)') '[Warning] observation is outside of the horizontal domain'
    end if
    rk = undef
    qc = iqc_out_h
    return
  end if
  !
  IF(elem > 9999) THEN ! surface observation
    rk = rlev
  ELSE
    !
    ! horizontal interpolation
    !
    i = CEILING(ri)
    j = CEILING(rj)
    !
    ! Find the lowest valid level
    !
    ks = 1+KHALO
    do jj = j-1, j
      do ii = i-1, i


!print *, p_full(:,ii,jj)

        DO k=1+KHALO,nlev+KHALO
          if (p_full(k,ii,jj) >= 0.0d0) exit
        END DO
        if (k > ks) ks = k
      end do
    end do
!    DO k=1+KHALO,nlev+KHALO
!!      IF(i <= nlon+IHALO) THEN
!        lnps(i-1:i,j-1:j) = LOG(p_full(k,i-1:i,j-1:j))
!!      ELSE
!!        lnps(i-1,j-1:j) = LOG(p_full(k,i-1,j-1:j))
!!        lnps(1,j-1:j) = LOG(p_full(k,1,j-1:j))
!!      END IF
!      CALL itpl_2d(lnps,ri,rj,plev(k))
!    END DO

    lnps(:,i-1:i,j-1:j) = LOG(p_full(:,i-1:i,j-1:j))
    call itpl_2d_column(lnps,ri,rj,plev)

    !
    ! Log pressure
    !
    rk = LOG(rlev)
    !
    ! determine if rk is within bound.
    !
    IF(rk < plev(nlev+KHALO)) THEN
      call itpl_2d(p_full(nlev+KHALO,:,:),ri,rj,ptmp)
      if (LOG_LEVEL >= 2) then
        write(6,'(A,F8.1,A,F8.1,A,I5)') '[Warning] observation is too high: ptop=', ptmp, ', lev=', rlev, ', elem=', elem
      end if
      rk = undef
      qc = iqc_out_vhi
      RETURN
    END IF
    IF(rk > plev(ks)) THEN
      call itpl_2d(p_full(ks,:,:),ri,rj,ptmp)
!print *, ks, rk, plev(ks)
      if (LOG_LEVEL >= 2) then
        write(6,'(A,F8.1,A,F8.1,A,I5)') '[Warning] observation is too low: pbottom=', ptmp, ', lev=', rlev, ', elem=', elem
      end if
      rk = undef
      qc = iqc_out_vlo

!print *, plev
!print *, elem, ri, rj, rlev, rk, qc


      RETURN
    END IF
    !
    ! find rk
    !
    DO k=ks+1,nlev+KHALO
      IF(plev(k) < rk) EXIT ! assuming descending order of plev
    END DO
    ak = (rk - plev(k-1)) / (plev(k) - plev(k-1))
    rk = REAL(k-1,r_size) + ak
  END IF

  RETURN
END SUBROUTINE phys2ijk
!-----------------------------------------------------------------------
! Coordinate conversion (find rk in height levels)
!
! rk = 0.0d0  : surface observation
!-----------------------------------------------------------------------
SUBROUTINE phys2ijkz(z_full,ri,rj,rlev,rk,qc)
  use scale_atmos_grid_cartesC_index, only: &
      KHALO
!  use common_mpi
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: z_full(nlevh,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev ! height levels
  REAL(r_size),INTENT(OUT) :: rk
  INTEGER,INTENT(OUT) :: qc
  REAL(r_size) :: ak
  REAL(r_size) :: zlev(nlevh)
  REAL(r_size) :: ztmp
  INTEGER :: i,j,k, ii, jj, ks

!  integer :: ierr
!  REAL(r_dble) :: rrtimer00,rrtimer
!  rrtimer00 = MPI_WTIME()


  qc = iqc_good
!
! rlev -> rk
!
  if (ri < 1.0d0 .or. ri > nlonh .or. rj < 1.0d0 .or. rj > nlath) then
    if (LOG_LEVEL >= 1) then
      write (6,'(A)') '[Warning] observation is outside of the horizontal domain'
    end if
    rk = undef
    qc = iqc_out_h
    return
  end if


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:check_domain:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !
  ! horizontal interpolation
  !
  i = CEILING(ri)
  j = CEILING(rj)
  !
  ! Find the lowest valid level
  !
  ks = 1+KHALO
  do jj = j-1, j
    do ii = i-1, i
      DO k=1+KHALO,nlev+KHALO
        if (z_full(k,ii,jj) > -300.0d0 .and. z_full(k,ii,jj) < 10000.0d0) exit
      END DO
      if (k > ks) ks = k
    end do
  end do


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:find_lowesr_level:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


!  DO k=1+KHALO,nlev+KHALO
!    CALL itpl_2d(z_full(k,:,:),ri,rj,zlev(k))
!  END DO

  call itpl_2d_column(z_full,ri,rj,zlev)


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:itpl_2d:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !
  ! determine if rlev is within bound.
  !
  IF(rlev > zlev(nlev+KHALO)) THEN
    call itpl_2d(z_full(nlev+KHALO,:,:),ri,rj,ztmp)
    if (LOG_LEVEL >= 2) then
      write(6,'(A,F8.1,A,F8.1)') '[Warning] observation is too high: ztop=', ztmp, ', lev=', rlev
    end if
    rk = undef
    qc = iqc_out_vhi
    RETURN
  END IF
  IF(rlev < zlev(ks)) THEN
    call itpl_2d(z_full(ks,:,:),ri,rj,ztmp)
    if (LOG_LEVEL >= 2) then
      write(6,'(A,F8.1,A,F8.1)') '[Warning] observation is too low: zbottom=', ztmp, ', lev=', rlev
    end if
    rk = undef
    qc = iqc_out_vlo
    RETURN
  END IF


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:check_height:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !
  ! find rk
  !
  DO k=ks+1,nlev+KHALO
    IF(zlev(k) > rlev) EXIT ! assuming ascending order of zlev
  END DO
  ak = (rlev - zlev(k-1)) / (zlev(k) - zlev(k-1))
  rk = REAL(k-1,r_size) + ak


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:calc_rk:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  RETURN
END SUBROUTINE phys2ijkz
!-----------------------------------------------------------------------
! Coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE phys2ij(rlon,rlat,rig,rjg)
  use scale_atmos_grid_cartesC, only: &
      CXG => ATMOS_GRID_CARTESC_CXG, &
      CYG => ATMOS_GRID_CARTESC_CYG, &
      DX, &
      DY
  use scale_mapprojection, only: &
      MAPPROJECTION_lonlat2xy
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(OUT) :: rig
  REAL(r_size),INTENT(OUT) :: rjg
!
! rlon,rlat -> ri,rj
!
  call MAPPROJECTION_lonlat2xy(rlon*pi/180.0_r_size,rlat*pi/180.0_r_size,rig,rjg)
  rig = (rig - CXG(1)) / DX + 1.0d0
  rjg = (rjg - CYG(1)) / DY + 1.0d0

  RETURN
END SUBROUTINE phys2ij

SUBROUTINE ij2phys(rig,rjg,rlon,rlat)
  use scale_atmos_grid_cartesC, only: &
      CXG => ATMOS_GRID_CARTESC_CXG, &
      CYG => ATMOS_GRID_CARTESC_CYG, &
      DX, &
      DY
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rig
  REAL(r_size),INTENT(IN) :: rjg
  REAL(r_size),INTENT(OUT) :: rlon ! (deg)
  REAL(r_size),INTENT(OUT) :: rlat ! (deg)
  REAL(r_size) :: x, y ! (m)
!
! ri,rj -> rlon,rlat
!
  x = (rig - 1.0d0) * DX + CXG(1) 
  y = (rjg - 1.0d0) * DY + CYG(1) 

  call MAPPROJECTION_xy2lonlat(x,y,rlon,rlat)

  rlon = rlon * rad2deg
  rlat = rlat * rad2deg

  RETURN
END SUBROUTINE ij2phys
!
!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
     & + var(i  ,j-1) *    ai  * (1-aj) &
     & + var(i-1,j  ) * (1-ai) *    aj  &
     & + var(i  ,j  ) *    ai  *    aj

  RETURN
END SUBROUTINE itpl_2d

SUBROUTINE itpl_2d_column(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlevh,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5(nlevh)
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  var5(:) = var(:,i-1,j-1) * (1-ai) * (1-aj) &
        & + var(:,i  ,j-1) *    ai  * (1-aj) &
        & + var(:,i-1,j  ) * (1-ai) *    aj  &
        & + var(:,i  ,j  ) *    ai  *    aj

  RETURN
END SUBROUTINE itpl_2d_column

SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlevh,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rk
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj,ak
  INTEGER :: i,j,k

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)
  k = CEILING(rk)
  ak = rk - REAL(k-1,r_size)

  var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
     & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
     & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
     & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
     & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
     & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
     & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
     & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak

  RETURN
END SUBROUTINE itpl_3d
!-----------------------------------------------------------------------
! Monitor observation departure by giving the v3dg,v2dg data
!-----------------------------------------------------------------------
subroutine monit_obs(v3dg,v2dg,topo,nobs,bias,rmse,monit_type,use_key,&
                     nobs_H08,bias_H08,rmse_H08,yobs_H08,step)!,bias_H08_bc,rmse_H08_bc,&
                     !aH08,bH08,vbcf,step)
  use scale_prc, only: &
      PRC_myrank
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO

  implicit none

  REAL(RP),intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),intent(in) :: v2dg(nlon,nlat,nv2d)
  real(r_size),intent(in) :: topo(nlon,nlat)
  INTEGER,INTENT(OUT) :: nobs(nid_obs)
  REAL(r_size),INTENT(OUT) :: bias(nid_obs)
  REAL(r_size),INTENT(OUT) :: rmse(nid_obs)
  LOGICAL,INTENT(OUT) :: monit_type(nid_obs)
  logical,intent(in) :: use_key
  integer,intent(in) :: step

  REAL(r_size) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size) :: v2dgh(nlonh,nlath,nv2dd)
  integer :: nnobs
  integer :: n,nn
  integer :: iset,iidx
  real(r_size) :: ril,rjl,rk,rkz

  real(r_size),allocatable :: oelm(:)
  real(r_size),allocatable :: ohx(:)
  integer,allocatable :: oqc(:)

  integer :: OMP_GET_NUM_THREADS, omp_chunk

!  REAL(r_size) :: timer
!  INTEGER :: ierr

  INTEGER,INTENT(OUT) :: nobs_H08(NIRB_HIM8)
  REAL(r_size),INTENT(OUT) :: bias_H08(NIRB_HIM8)
  REAL(r_size),INTENT(OUT) :: rmse_H08(NIRB_HIM8)
!  REAL(r_size),INTENT(OUT) :: bias_H08_bc(NIRB_HIM8)
!  REAL(r_size),INTENT(OUT) :: rmse_H08_bc(NIRB_HIM8)

  real(r_size),intent(out) :: yobs_H08(nlon,nlat,NIRB_HIM8)
  real(r_size) :: yobs_H08_monit(nlon,nlat,NIRB_HIM8)
  real(r_size) :: plev_obs_H08(nlon,nlat,NIRB_HIM8)
!  real(r_size) :: yobs_H08_bc(nlon,nlat,NIRB_HIM8)
  real(r_size) :: yobs_H08_clr(nlon,nlat,NIRB_HIM8)
  integer :: qc_H08(nlon,nlat,NIRB_HIM8)
  real(r_size) :: zangle_H08(nlon,nlat)

  integer :: i8, j8, b8
  integer :: ch

!  ! bias correction
!  integer,parameter :: nmin = 400 ! parameter from Miyoshi et al. (2010) & Sato (2007)
!  integer :: nobs_b
!  integer :: i, j, didx, npr 
!  real(r_size), intent(in) :: vbcf(H08_NPRED,NIRB_HIM8)
!  real(r_size), intent(out) :: aH08(H08_NPRED,H08_NPRED,NIRB_HIM8)
!  real(r_size), intent(out) :: bH08(H08_NPRED,NIRB_HIM8) ! bias correction
!  real(r_size), allocatable :: pred(:,:)
!  real(r_size) :: tmp
!  real(r_size) :: pbeta, predt

#ifdef TCV
! Multiple TCs are not considered (04/14/2017)
!  real(r_size) :: TC_rij(2) = -1.0d0
!  integer :: obs_n_TCX, obs_n_TCY, obs_n_TCP ! TCX, TCY, TCP
!  real(r_size) :: bTC(3) = 9.99d33
!  integer :: k
#endif

  call state_to_history(v3dg, v2dg, topo, v3dgh, v2dgh)

  if (use_key) then
    nnobs = obsda_sort%nobs_in_key
  else
    nnobs = obsda_sort%nobs
  end if

  allocate (oelm(nnobs))
  allocate (ohx(nnobs))
  allocate (oqc(nnobs))

#ifdef DEBUG
  if (step < 0 .or. step > 2) then
    write (6, *) '[Error] monit_obs: step should be 0, 1, or 2.'
    stop
  end if
#endif
  if (step == 1) then
    obsdep_nobs = nnobs
    allocate (obsdep_set(obsdep_nobs))
    allocate (obsdep_idx(obsdep_nobs))
    allocate (obsdep_qc (obsdep_nobs))
    allocate (obsdep_omb(obsdep_nobs))
    allocate (obsdep_oma(obsdep_nobs))
  end if

  oqc = -1

!  obs_idx_TCX = -1
!  obs_idx_TCY = -1
!  obs_idx_TCP = -1

!  yobs_H08 = -1.0d0
!  yobs_H08_monit = -1.0d0
!!  if (USE_HIM8) then 
! Always calculate Him8 radiances

write(6,'(a,i9)'),"Hello from monit_obs before Trans",obsdep_nobs
  call Trans_XtoY_H08_allg(v3dgh,v2dgh,yobs_H08,yobs_H08_clr,&
                           plev_obs_H08,qc_H08,zangle_H08)
write(6,'(a)'),"Hello from monit_obs afte Trans"
  !
  ! Initialize qc flag (set iqc is "bad")
  ! This will be overwritten by obsda_sort%qc
  qc_H08 = iqc_obs_bad
  yobs_H08_monit = yobs_H08

!  endif

!##!$OMP PARALLEL PRIVATE(n,nn,iset,iidx,ril,rjl,rk,rkz,i8,j8,b8,ch)
!##  omp_chunk = min(4, max(1, (nnobs-1) / OMP_GET_NUM_THREADS() + 1))
!##!$OMP DO SCHEDULE(DYNAMIC,omp_chunk)
  do n = 1, nnobs

    if (use_key) then
      nn = obsda_sort%key(n)
    else
      nn = n
    end if

    iset = obsda_sort%set(nn)
    iidx = obsda_sort%idx(nn)

    if (step == 1) then
      obsdep_set(n) = iset
      obsdep_idx(n) = iidx
#ifdef DEBUG
    else if (step == 2) then
      if (obsdep_set(n) /= iset) then
        write (6, *) "[Error] 'set' for y_b and y_a are inconsistent!"
        stop
      end if
      if (obsdep_idx(n) /= iidx) then
        write (6, *) "[Error] 'idx' for y_b and y_a are inconsistent!"
        stop
      end if
#endif
    end if

#ifdef DEBUG
    if (obsda_sort%qc(nn) /= iqc_good) then
      write (6, *) "[Error] The QC value of this observation provided for monitoring is not good: ", &
                   obsda_sort%qc(nn)
      stop
    end if
#endif

    oelm(n) = obs(iset)%elm(iidx)


    call rij_g2l(PRC_myrank, obs(iset)%ri(iidx), obs(iset)%rj(iidx), ril, rjl)
#ifdef DEBUG
    if (PRC_myrank /= obs(iset)%rank(iidx) .or. obs(iset)%rank(iidx) == -1) then
      write (6, *) "[Error] This observation provided for monitoring does not reside in my rank: ", &
                   PRC_myrank, obs(iset)%rank(iidx), obs(iset)%ri(iidx), obs(iset)%rj(iidx), ril, rjl
      stop
    end if
#endif

    if (DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
        abs(obs(iset)%dif(iidx)) <= DEPARTURE_STAT_T_RANGE) then

      oqc(n) = iqc_otype

      select case (OBS_IN_FORMAT(iset))
      !=========================================================================
      case (obsfmt_prepbufr)
      !-------------------------------------------------------------------------
        call phys2ijk(v3dgh(:,:,:,iv3dd_p),obs(iset)%elm(iidx), &
                      ril,rjl,obs(iset)%lev(iidx),rk,oqc(n))
        if (oqc(n) == iqc_good) then
          call Trans_XtoY(obs(iset)%elm(iidx),ril,rjl,rk, &
                          obs(iset)%lon(iidx),obs(iset)%lat(iidx), &
                          v3dgh,v2dgh,ohx(n),oqc(n),stggrd=1)
        end if
      !=========================================================================
      case (obsfmt_radar)
      !-------------------------------------------------------------------------
        if (DEPARTURE_STAT_RADAR) then
          call phys2ijkz(v3dgh(:,:,:,iv3dd_hgt),ril,rjl,obs(iset)%lev(iidx),rkz,oqc(n))
          if (oqc(n) == iqc_good) then
            call Trans_XtoY_radar(obs(iset)%elm(iidx),obs(iset)%meta(1), &
                                  obs(iset)%meta(2),obs(iset)%meta(3),ril,rjl,rkz, &
                                  obs(iset)%lon(iidx),obs(iset)%lat(iidx), &
                                  obs(iset)%lev(iidx),v3dgh,v2dgh,ohx(n),oqc(n),stggrd=1)
            if (oqc(n) == iqc_ref_low) oqc(n) = iqc_good ! when process the observation operator, we don't care if reflectivity is too small
          end if
        end if
      !=========================================================================
      case (obsfmt_h08)
      !-------------------------------------------------------------------------
        i8 = nint(ril-IHALO)
        j8 = nint(rjl-JHALO)
        b8 = nint(obs(iset)%lev(iidx))

        ohx(n) = yobs_H08(i8,j8,b8-6) ! Obs - B/A
        oqc(n) = obsda_sort%qc(nn) !qc_H08(i8,j8,b8-6) ! QC

        yobs_H08_monit(i8,j8,b8-6) = obs(iset)%dat(iidx) - yobs_H08(i8,j8,b8-6) 
        qc_H08(i8,j8,b8-6) = obsda_sort%qc(nn)

        if(sum(H08_BAND_USE) /= NIRB_HIM8)then
          do ch = 1, NIRB_HIM8            
            if (H08_BAND_USE(ch) /= 1) then
              yobs_H08_monit(i8,j8,ch) = obs(iset)%dat(iidx-(b8-6)+ch) - yobs_H08(i8,j8,ch) 
              qc_H08(i8,j8,ch) = oqc(n)
            endif
          enddo
        endif

!      if(H08_VBC_USE)then
!        pbeta = 0.0d0
!        do npr = 1, H08_NPRED
!          if(npr == 1) predt = zangle_H08(np)
!          if(npr == 2) predt = yobs_H08_bc(ch,np)
!          pbeta = pbeta + predt * vbcf(npr,ch)
!        enddo
!        yobs_H08_bc(ch,np) = yobs_H08_bc(ch,np) - pbeta
!       
!        if (ch == nint(obs(iset)%dat(iidx))) then
!          ohx(n) = ohx(n) - pbeta
!        endif
!      endif
!
!    do np = 1, nprof
!      do ch = 1, NIRB_HIM8
!        yobs_H08(ch,np) = obs(iset)%dat(prof2B07(np)+ch-1) - yobs_H08(ch,np)
!        yobs_H08_bc(ch,np) = obs(iset)%dat(prof2B07(np)+ch-1) - yobs_H08_bc(ch,np)
!      enddo ! [ch = 1, NIRB_HIM8]
!    enddo ! [np = 1, nprof]


      !=========================================================================
      end select

      if (oqc(n) == iqc_good) then
        ohx(n) = obs(iset)%dat(iidx) - ohx(n)
      else
        ohx(n) = undef
      end if

      if (step == 1) then
        obsdep_qc(n) = oqc(n)
        obsdep_omb(n) = ohx(n)
      else if (step == 2) then
        if (obsdep_qc(n) == iqc_good) then ! Use the QC value of y_a only if the QC of y_b is good
          obsdep_qc(n) = oqc(n)            !
        end if                             !
        obsdep_oma(n) = ohx(n)
      end if

      if (LOG_LEVEL >= 3) then
        write (6, '(2I6,2F8.2,4F12.4,I3)') &
              obs(iset)%elm(iidx), &
              obs(iset)%typ(iidx), &
              obs(iset)%lon(iidx), &
              obs(iset)%lat(iidx), &
              obs(iset)%lev(iidx), &
              obs(iset)%dat(iidx), &
              obs(iset)%err(iidx), &
              ohx(n), &
              oqc(n)
      end if

    end if ! [ DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
           !   abs(obs(iset)%dif(iidx)) <= DEPARTURE_STAT_T_RANGE ]

  end do ! [ n = 1, nnobs ]
!##!$OMP END DO
!##!$OMP END PARALLEL

!   -- TC vital DA -- 
!   -- End of TC vital DA -- 

write(6,'(a)')"DEBUG before monit dep"
  call monit_dep(nnobs,oelm,ohx,oqc,nobs,bias,rmse)
write(6,'(a)')"DEBUG after monit dep"

!  aH08 = 0.0d0
!  bH08 = 0.0d0
  if (USE_HIM8) then
!    call monit_dep_H08(yobs_H08,yobs_H08_bc,qc_H08,nobs_H08,bias_H08,rmse_H08)!,bias_H08_bc,rmse_H08_bc)
    call monit_dep_H08(yobs_H08_monit,qc_H08,nobs_H08,bias_H08,rmse_H08)
!
!    ! bias correction
!    aH08 = 0.0d0
!    bH08 = 0.0d0
!
!    if(step == 2 .and. H08_VBC_USE) then
!
!      do ch = 1, NIRB_HIM8
!        if(H08_BAND_USE(ch) /= 1) cycle
!        nobs_b = nobs_H08(ch)
!    allocate(bTC(3,0:nprocs_d-1))
!    allocate(bufr(3,0:nproces_d-1))
!
!        if(nobs_b > 0)then
!
!          if(allocated(pred)) deallocate(pred)
!          allocate(pred(H08_NPRED,nobs_b)) ! pr 
!
!          do np = 1, nprof ! profile/obs loop
!            if(qc_H08(ch,np) /= iqc_good) cycle
!            n = prof2nda(np)
!            if (use_key) then
!              nn = obsda_sort%key(n)
!            else
!              nn = n
!            end if
!            iset = obsda_sort%set(nn)
!            iidx = obsda_sort%idx(nn)
!
!            didx = (prof2B07(np) + ch - 1) - iidx
!
!            if(didx /= 0)then
!              nn = nn + didx
!              iset = obsda_sort%set(nn)
!              iidx = obsda_sort%idx(nn)
!            endif
!     
!            !
!            ! p R^-1 p^T
!            !
!            do j = 1, H08_NPRED
!              if(j == 1) pred(j,np) = zangle_H08(np)
!              if(j == 2) pred(j,np) = obs(iset)%dat(iidx)
!
!              do i = 1, H08_NPRED
!                if(i == 1) pred(i,np) = zangle_H08(np)
!                if(i == 2) pred(i,np) = obs(iset)%dat(iidx)
!
!                aH08(i,j,ch) = aH08(i,j,ch) + pred(i,np) * pred(j,np) / (obsda_sort%val2(nn)**2)
!              enddo ! i
!            enddo ! j
!            !
!            ! Add B_beta^-1
!            !
!            ! Sato (2007)
!            if(nobs_b < nmin)then
!              tmp = real(nmin, kind=r_size) / (obsda_sort%val2(nn)**2)
!            else
!              tmp = real(nobs_b, kind=r_size) / dlog10(real(nobs_b,kind=r_size)/real(nmin,kind=r_size)) &
!                      / (obsda_sort%val2(nn)**2)
!            endif
!            do i = 1, H08_NPRED
!              aH08(i,i,ch) = aH08(i,i,ch) + tmp
!            enddo
!            !
!            ! p R^-1 d
!            !
!            do i = 1, H08_NPRED
!              bH08(i,ch) = bH08(i,ch) + pred(i,np)  / (obsda_sort%val2(nn)**2) &
!                             & * yobs_H08(ch,np)
!            enddo 
!          enddo ! np
!        endif ! [nobs_b > 0]
!
!        deallocate(pred)
!
!      enddo ! ch
!
!    endif ! [step == 2]
!
!    !deallocate(yobs_H08, qc_H08)
  else 
    nobs_H08 = 0
    bias_H08 = 0.0d0
    rmse_H08 = 0.0d0
!    bias_H08_bc = 0.0d0
!    rmse_H08_bc = 0.0d0
  endif ! [USE_HIM8]
  
write(6,'(a)')"DEBUG after USE_HIM8"

  monit_type = .false.
  monit_type(uid_obs(id_u_obs)) = .true.
  monit_type(uid_obs(id_v_obs)) = .true.
  monit_type(uid_obs(id_t_obs)) = .true.
  monit_type(uid_obs(id_tv_obs)) = .true.
  monit_type(uid_obs(id_q_obs)) = .true.
!  monit_type(uid_obs(id_rh_obs)) = .true.
  monit_type(uid_obs(id_ps_obs)) = .true.
  if (DEPARTURE_STAT_RADAR) then
    monit_type(uid_obs(id_radar_ref_obs)) = .true.
    monit_type(uid_obs(id_radar_ref_zero_obs)) = .true.
    monit_type(uid_obs(id_radar_vr_obs)) = .true.
!    monit_type(uid_obs(id_radar_prh_obs)) = .true.
  end if
  if (DEPARTURE_STAT_H08) then
    monit_type(uid_obs(id_H08IR_obs)) = .true.
  end if

#ifdef TCV
!  monit_type(uid_obs(id_tclon_obs)) = .true.
!  monit_type(uid_obs(id_tclat_obs)) = .true.
!  monit_type(uid_obs(id_tcmip_obs)) = .true.
#endif

  deallocate (oelm)
  deallocate (ohx)
  deallocate (oqc)

write(6,'(a)')"DEBUG after USE_HIM82"

  return
end subroutine monit_obs
!-----------------------------------------------------------------------
! Monitor departure
!  ofmt: output format
!    0: U,V,T(Tv),Q,RH,PS (default)
!    1: U,V,T(Tv),Q,RH,PS,RAIN
!-----------------------------------------------------------------------
SUBROUTINE monit_dep(nn,elm,dep,qc,nobs,bias,rmse)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elm(nn)
  REAL(r_size),INTENT(IN) :: dep(nn)
  INTEGER,INTENT(IN) :: qc(nn)
  INTEGER,INTENT(OUT) :: nobs(nid_obs)
  REAL(r_size),INTENT(OUT) :: bias(nid_obs)
  REAL(r_size),INTENT(OUT) :: rmse(nid_obs)
  INTEGER :: n,i,ielm

  nobs = 0
  bias = 0.0d0
  rmse = 0.0d0

!!!!!!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i,ielm) REDUCTION(+:nobs,bias,rmse)
  DO n=1,nn
    IF(qc(n) /= iqc_good) CYCLE

    ielm = NINT(elm(n))
    if (ielm == id_tv_obs) then ! compute Tv as T
      ielm = id_t_obs
    end if
    if (ielm == id_radar_ref_zero_obs) then ! compute RE0 as REF
      ielm = id_radar_ref_obs
    end if

    i = uid_obs(ielm)
    nobs(i) = nobs(i) + 1
    bias(i) = bias(i) + dep(n)
    rmse(i) = rmse(i) + dep(n)**2
  END DO

  DO i = 1, nid_obs
    IF(nobs(i) == 0) THEN
      bias(i) = undef
      rmse(i) = undef
    ELSE
      bias(i) = bias(i) / REAL(nobs(i),r_size)
      rmse(i) = SQRT(rmse(i) / REAL(nobs(i),r_size))
    END IF
  END DO

  RETURN
END SUBROUTINE monit_dep
!-----------------------------------------------------------------------
! Monitor departure
!-----------------------------------------------------------------------
SUBROUTINE monit_print(nobs,bias,rmse,monit_type)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nobs(nid_obs)
  REAL(r_size),INTENT(IN) :: bias(nid_obs)
  REAL(r_size),INTENT(IN) :: rmse(nid_obs)
  LOGICAL,INTENT(IN),OPTIONAL :: monit_type(nid_obs)

  character(12) :: var_show(nid_obs)
  character(12) :: nobs_show(nid_obs)
  character(12) :: bias_show(nid_obs)
  character(12) :: rmse_show(nid_obs)

  integer :: i, n
  character(4) :: nstr

  logical :: monit_type_(nid_obs)

  monit_type_ = .true.
  if (present(monit_type)) monit_type_ = monit_type

  n = 0
  do i = 1, nid_obs
!!!!!!
!    if (monit_type_(i) .and. i /= uid_obs(id_tv_obs)) then
    if (monit_type_(i) .and. i /= uid_obs(id_tv_obs) .and. i /= uid_obs(id_radar_ref_zero_obs)) then
!!!!!!
      n = n + 1
      write(var_show(n),'(A12)') obelmlist(i)
      write(nobs_show(n),'(I12)') nobs(i)
      if (nobs(i) > 0) then
        write(bias_show(n),'(ES12.3)') bias(i)
        write(rmse_show(n),'(ES12.3)') rmse(i)
      else
        write(bias_show(n),'(A12)') 'N/A'
        write(rmse_show(n),'(A12)') 'N/A'
      end if
    end if
  end do
  write(nstr, '(I4)') n

  WRITE(6,'(A,' // trim(nstr) // "('============'))") '======'
  WRITE(6,'(6x,' // trim(nstr) // 'A)')          var_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // "('------------'))") '------'
  WRITE(6,'(A,' // trim(nstr) // 'A)') 'BIAS  ', bias_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') 'RMSE  ', rmse_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') 'NUMBER', nobs_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // "('============'))") '======'

  RETURN
END SUBROUTINE monit_print
!
! monitor for Himawari-8 IR observations --
!SUBROUTINE monit_dep_H08(dep,dep_bc,qc,nobs,bias,rmse)!,bias_bc,rmse_bc)
SUBROUTINE monit_dep_H08(dep,qc,nobs,bias,rmse)!,bias_bc,rmse_bc)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: dep(nlon,nlat,NIRB_HIM8)
!  REAL(r_size),INTENT(IN) :: dep_bc(nlon,nlat,NIRB_HIM8)
  INTEGER,INTENT(IN) :: qc(nlon,nlat,NIRB_HIM8)
  INTEGER,INTENT(OUT) :: nobs(NIRB_HIM8)
  REAL(r_size),INTENT(OUT) :: bias(NIRB_HIM8)
  REAL(r_size),INTENT(OUT) :: rmse(NIRB_HIM8)
!  REAL(r_size),INTENT(OUT) :: bias_bc(NIRB_HIM8)
!  REAL(r_size),INTENT(OUT) :: rmse_bc(NIRB_HIM8)
!  INTEGER :: n
  integer i, j, ch  

  nobs = 0
  rmse = 0.0d0
  bias = 0.0d0
!  rmse_bc = 0.0d0
!  bias_bc = 0.0d0
    
  do ch = 1, NIRB_HIM8
    do j = 1, nlat
    do i = 1, nlon

      if(qc(i,j,ch) /= iqc_good) cycle
      if(dep(i,j,ch) /= dep(i,j,ch)) cycle ! NaN QC

      nobs(ch) = nobs(ch) + 1
      bias(ch) = bias(ch) + dep(i,j,ch)
      rmse(ch) = rmse(ch) + dep(i,j,ch)**2
!      bias_bc(ch) = bias_bc(ch) + dep_bc(i,j,ch)
!      rmse_bc(ch) = rmse_bc(ch) + dep_bc(i,j,ch)**2

    enddo
    enddo
  enddo

  DO ch = 1, NIRB_HIM8
    IF(nobs(ch) == 0) THEN
      bias(ch) = undef
      rmse(ch) = undef
!      bias_bc(ch) = undef
!      rmse_bc(ch) = undef
    ELSE
      bias(ch) = bias(ch) / REAL(nobs(ch),r_size)
      rmse(ch) = SQRT(rmse(ch) / REAL(nobs(ch),r_size))
!      bias_bc(ch) = bias_bc(ch) / REAL(nobs(ch),r_size)
!      rmse_bc(ch) = SQRT(rmse_bc(ch) / REAL(nobs(ch),r_size))
    END IF
  END DO

  RETURN
END SUBROUTINE monit_dep_H08
!
SUBROUTINE monit_print_H08(nobs,bias,rmse,monit_type)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nobs(NIRB_HIM8)
  REAL(r_size),INTENT(IN) :: bias(NIRB_HIM8)
  REAL(r_size),INTENT(IN) :: rmse(NIRB_HIM8)
  LOGICAL,INTENT(IN),OPTIONAL :: monit_type(NIRB_HIM8)
  
  character(12) :: var_show(NIRB_HIM8)
  character(12) :: nobs_show(NIRB_HIM8)
  character(12) :: bias_show(NIRB_HIM8)
  character(12) :: rmse_show(NIRB_HIM8)
  character(12) :: flag_show(NIRB_HIM8)

  integer :: i, n
  character(4) :: nstr
  character(12) :: tmpstr(NIRB_HIM8)
  character(12) :: tmpstr2(NIRB_HIM8)

  character(3) :: B3(NIRB_HIM8)

  logical :: monit_type_(NIRB_HIM8)

  monit_type_ = .true.
  if (present(monit_type)) monit_type_ = monit_type

  n = 0
  do i = 1, NIRB_HIM8
    n = n + 1

    B3(i) = ch2BB_Him8(i)

    if(H08_BAND_USE(i) == 1)then
      write(flag_show(n),'(A12)') "YES"
    else
      write(flag_show(n),'(A12)') " NO"
    endif
    write(var_show(n),'(A12)') B3(i)
    write(nobs_show(n),'(I12)') nobs(i)
    if (nobs(i) > 0) then
      write(bias_show(n),'(ES12.3)') bias(i)
      write(rmse_show(n),'(ES12.3)') rmse(i)
    else
      write(bias_show(n),'(A12)') 'N/A'
      write(rmse_show(n),'(A12)') 'N/A'
    end if
  end do
  write(nstr, '(I4)') n
  tmpstr(1:n) = '============'
  tmpstr2(1:n) = '------------'

  WRITE(6,'(A,' // trim(nstr) // 'A)') '======', tmpstr(1:n)
  WRITE(6,'(6x,' // trim(nstr) // 'A)')          var_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') '------', tmpstr2(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') 'BIAS  ', bias_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') 'RMSE  ', rmse_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') 'NUMBER', nobs_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') 'USED? ', flag_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') '======', tmpstr(1:n)

  RETURN
END SUBROUTINE monit_print_H08
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_info_allocate(obs, extended)
  IMPLICIT NONE
  TYPE(obs_info),INTENT(INOUT) :: obs
  logical, optional, intent(in) :: extended

  call obs_info_deallocate(obs)

  ALLOCATE( obs%elm (obs%nobs) )
  ALLOCATE( obs%lon (obs%nobs) )
  ALLOCATE( obs%lat (obs%nobs) )
  ALLOCATE( obs%lev (obs%nobs) )
  ALLOCATE( obs%dat (obs%nobs) )
  ALLOCATE( obs%err (obs%nobs) )
  ALLOCATE( obs%typ (obs%nobs) )
  ALLOCATE( obs%dif (obs%nobs) )

  obs%elm = 0
  obs%lon = 0.0d0
  obs%lat = 0.0d0
  obs%lev = 0.0d0
  obs%dat = 0.0d0
  obs%err = 0.0d0
  obs%typ = 0
  obs%dif = 0.0d0

  if (present(extended)) then
    if (extended) then
      allocate( obs%ri (obs%nobs) )
      allocate( obs%rj (obs%nobs) )
      allocate( obs%rank (obs%nobs) )

      obs%ri = 0.0d0
      obs%rj = 0.0d0
      obs%rank = -1
    end if
  end if

  RETURN
END SUBROUTINE obs_info_allocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_info_deallocate(obs)
  IMPLICIT NONE
  TYPE(obs_info),INTENT(INOUT) :: obs

  IF(ALLOCATED(obs%elm)) DEALLOCATE(obs%elm)
  IF(ALLOCATED(obs%lon)) DEALLOCATE(obs%lon)
  IF(ALLOCATED(obs%lat)) DEALLOCATE(obs%lat)
  IF(ALLOCATED(obs%lev)) DEALLOCATE(obs%lev)
  IF(ALLOCATED(obs%dat)) DEALLOCATE(obs%dat)
  IF(ALLOCATED(obs%err)) DEALLOCATE(obs%err)
  IF(ALLOCATED(obs%typ)) DEALLOCATE(obs%typ)
  IF(ALLOCATED(obs%dif)) DEALLOCATE(obs%dif)

  RETURN
END SUBROUTINE obs_info_deallocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_da_value_allocate(obsda,member)
  IMPLICIT NONE
  TYPE(obs_da_value),INTENT(INOUT) :: obsda
  INTEGER,INTENT(IN) :: member

  call obs_da_value_deallocate(obsda)

  ALLOCATE( obsda%set (obsda%nobs) )
  ALLOCATE( obsda%idx (obsda%nobs) )
  ALLOCATE( obsda%key (obsda%nobs) )
  ALLOCATE( obsda%val (obsda%nobs) )
#ifdef H08
  ALLOCATE( obsda%lev (obsda%nobs) ) ! Him8
  ALLOCATE( obsda%val2 (obsda%nobs) ) ! Him8
  ALLOCATE( obsda%sprd (obsda%nobs) ) ! Him8
!  ALLOCATE( obsda%pred1 (obsda%nobs) ) ! Him8
!  ALLOCATE( obsda%pred2 (obsda%nobs) ) ! Him8
#endif
  ALLOCATE( obsda%qc  (obsda%nobs) )

  obsda%nobs_in_key = 0
  obsda%idx = 0
  obsda%key = 0
  obsda%val = 0.0d0
#ifdef H08
  obsda%lev = 0.0d0 ! Him8
  obsda%val2 = 0.0d0 ! Him8
  obsda%sprd = 0.0d0 ! Him8
!  obsda%pred1 = 0.0d0 ! Him8
!  obsda%pred2 = 0.0d0 ! Him8
#endif
  obsda%qc = 0

  if (member >= 0) then
    ALLOCATE( obsda%ensval (member,obsda%nobs) )
    obsda%ensval = 0.0d0
  end if

  RETURN
END SUBROUTINE obs_da_value_allocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_da_value_deallocate(obsda)
  IMPLICIT NONE
  TYPE(obs_da_value),INTENT(INOUT) :: obsda

  obsda%nobs_in_key = 0

  IF(ALLOCATED(obsda%set   )) DEALLOCATE(obsda%set   )
  IF(ALLOCATED(obsda%idx   )) DEALLOCATE(obsda%idx   )
  IF(ALLOCATED(obsda%key   )) DEALLOCATE(obsda%key   )
  IF(ALLOCATED(obsda%val   )) DEALLOCATE(obsda%val   )
#ifdef H08
  IF(ALLOCATED(obsda%lev   )) DEALLOCATE(obsda%lev   ) ! Him8
  IF(ALLOCATED(obsda%val2   )) DEALLOCATE(obsda%val2   ) ! Him8
  IF(ALLOCATED(obsda%sprd   )) DEALLOCATE(obsda%sprd   ) ! Him8
!  IF(ALLOCATED(obsda%pred1   )) DEALLOCATE(obsda%pred1   ) ! Him8
!  IF(ALLOCATED(obsda%pred2   )) DEALLOCATE(obsda%pred2   ) ! Him8
#endif
  IF(ALLOCATED(obsda%ensval)) DEALLOCATE(obsda%ensval)
  IF(ALLOCATED(obsda%qc    )) DEALLOCATE(obsda%qc    )

  RETURN
END SUBROUTINE obs_da_value_deallocate
!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nrec,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nrec
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl),ALLOCATABLE :: wk(:)
!  INTEGER :: ios
!  INTEGER :: iu,iv,it,iq,irh,ips,itc
  INTEGER :: iunit
  LOGICAL :: ex
  INTEGER :: sz

  ALLOCATE(wk(nrec))
  nn = 0
!  iu = 0
!  iv = 0
!  it = 0
!  iq = 0
!  irh = 0
!  ips = 0
!  itc = 0
  iunit=91
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

! get file size by reading through the entire file... too slow for big files
!-----------------------------
!    DO
!      READ(iunit,IOSTAT=ios) wk
!      IF(ios /= 0) EXIT
!!      SELECT CASE(NINT(wk(1)))
!!      CASE(id_u_obs)
!!        iu = iu + 1
!!      CASE(id_v_obs)
!!        iv = iv + 1
!!      CASE(id_t_obs,id_tv_obs)
!!        it = it + 1
!!      CASE(id_q_obs)
!!        iq = iq + 1
!!      CASE(id_rh_obs)
!!        irh = irh + 1
!!      CASE(id_ps_obs)
!!        ips = ips + 1
!!      CASE(id_tclon_obs)
!!        itc = itc + 1
!!      END SELECT
!      nn = nn + 1
!    END DO

! get file size by INQUIRE statement... may not work for some older fortran compilers
!-----------------------------
    INQUIRE(UNIT=iunit, SIZE=sz)
    IF (MOD(sz, r_sngl * (nrec+2)) /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    nn = sz / (r_sngl * (nrec+2))
!-----------------------------

    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
!    WRITE(6,'(A12,I10)') '          U:',iu
!    WRITE(6,'(A12,I10)') '          V:',iv
!    WRITE(6,'(A12,I10)') '      T(Tv):',it
!    WRITE(6,'(A12,I10)') '          Q:',iq
!    WRITE(6,'(A12,I10)') '         RH:',irh
!    WRITE(6,'(A12,I10)') '         Ps:',ips
!    WRITE(6,'(A12,I10)') '   TC TRACK:',itc
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF
  DEALLOCATE(wk)

  RETURN
END SUBROUTINE get_nobs

SUBROUTINE read_obs(cfile,obs)
  use scale_mapprojection, only: &
      MAPPROJECTION_lonlat2xy
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(INOUT) :: obs
  REAL(r_sngl) :: wk(8)
  REAL(r_size) :: x, y
  INTEGER :: n,iunit

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,obs%nobs
    READ(iunit) wk
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_v_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_t_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_tv_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_q_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_ps_obs)
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = wk(6) * 100.0 ! hPa -> Pa
    CASE(id_rh_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = wk(5) * 0.01 ! percent input
      wk(6) = wk(6) * 0.01 ! percent input
    CASE(id_tcmip_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = real(OBSERR_TCP,kind=r_sngl)
    CASE(id_tclon_obs)
      call MAPPROJECTION_lonlat2xy(REAL(wk(2),kind=r_size)*pi/180.0_r_size,&
                                   REAL(wk(3),kind=r_size)*pi/180.0_r_size,&
                                   x,y)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = real(x,kind=r_sngl)
      wk(6) = real(OBSERR_TCXY,kind=r_sngl)
    CASE(id_tclat_obs)
      call MAPPROJECTION_lonlat2xy(REAL(wk(2),kind=r_size)*pi/180.0_r_size,&
                                   REAL(wk(3),kind=r_size)*pi/180.0_r_size,&
                                   x,y)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = real(y,kind=r_sngl)
      wk(6) = real(OBSERR_TCXY,kind=r_sngl)
    END SELECT
    obs%elm(n) = NINT(wk(1))
    obs%lon(n) = REAL(wk(2),r_size)
    obs%lat(n) = REAL(wk(3),r_size)
    obs%lev(n) = REAL(wk(4),r_size)
    obs%dat(n) = REAL(wk(5),r_size)
    obs%err(n) = REAL(wk(6),r_size)
    obs%typ(n) = NINT(wk(7))
    obs%dif(n) = REAL(wk(8),r_size)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs

SUBROUTINE write_obs(cfile,obs,append,missing)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(IN) :: obs
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL,INTENT(IN),OPTIONAL :: missing
  LOGICAL :: append_
  LOGICAL :: missing_
  REAL(r_sngl) :: wk(8)
  INTEGER :: n,iunit

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
  DO n=1,obs%nobs
    if (missing_ .or. abs(obs%dat(n) - undef) > tiny(wk(5))) then
      wk(1) = REAL(obs%elm(n),r_sngl)
      wk(2) = REAL(obs%lon(n),r_sngl)
      wk(3) = REAL(obs%lat(n),r_sngl)
      wk(4) = REAL(obs%lev(n),r_sngl)
      wk(5) = REAL(obs%dat(n),r_sngl)
      wk(6) = REAL(obs%err(n),r_sngl)
      wk(7) = REAL(obs%typ(n),r_sngl)
      wk(8) = REAL(obs%dif(n),r_sngl)
      SELECT CASE(NINT(wk(1)))
      CASE(id_u_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_v_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_t_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_tv_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_q_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_ps_obs)
        wk(5) = wk(5) * 0.01 ! Pa -> hPa
        wk(6) = wk(6) * 0.01 ! Pa -> hPa
      CASE(id_rh_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
        wk(5) = wk(5) * 100.0 ! percent output
        wk(6) = wk(6) * 100.0 ! percent output
      CASE(id_tcmip_obs)
        wk(5) = wk(5) * 0.01 ! Pa -> hPa
        wk(6) = wk(6) * 0.01 ! Pa -> hPa
      END SELECT
      WRITE(iunit) wk
    end if
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs

SUBROUTINE read_obs_da(cfile,obsda,im)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_da_value),INTENT(INOUT) :: obsda
  INTEGER,INTENT(IN) :: im
#ifdef H08
  REAL(r_sngl) :: wk(8) ! H08
#else
  REAL(r_sngl) :: wk(4) ! H08
#endif
  INTEGER :: n,iunit

!  call obs_da_value_allocate(obsda)

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,obsda%nobs
    READ(iunit) wk
    obsda%set(n) = NINT(wk(1))
    obsda%idx(n) = NINT(wk(2))  !!!!!! will overflow......
    if (im == 0) then
      obsda%val(n) = REAL(wk(3),r_size)
    else
      obsda%ensval(im,n) = REAL(wk(3),r_size)
    end if
    obsda%qc(n) = NINT(wk(4))
#ifdef H08
    obsda%lev(n) = REAL(wk(5),r_size) ! Him8
    obsda%val2(n) = REAL(wk(6),r_size) ! Him8
!    obsda%pred1(n) = REAL(wk(7),r_size) ! Him8
!    obsda%pred2(n) = REAL(wk(8),r_size) ! Him8
#endif
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_da

SUBROUTINE write_obs_da(cfile,obsda,im,append)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_da_value),INTENT(IN) :: obsda
  INTEGER,INTENT(IN) :: im
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL :: append_
#ifdef H08
!  REAL(r_sngl) :: wk(8) ! H08
  REAL(r_sngl) :: wk(6) ! H08
#else
  REAL(r_sngl) :: wk(4) 
#endif
  INTEGER :: n,iunit

  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  IF(append_) THEN
    IF(obsda%nobs <= 0) RETURN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append',STATUS='replace')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',STATUS='replace')
  END IF
  DO n=1,obsda%nobs
    wk(1) = REAL(obsda%set(n),r_sngl)
    wk(2) = REAL(obsda%idx(n),r_sngl)  !!!!!! will overflow......
    if (im == 0) then
      wk(3) = REAL(obsda%val(n),r_sngl)
    else
      wk(3) = REAL(obsda%ensval(im,n),r_sngl)
    end if
    wk(4) = REAL(obsda%qc(n),r_sngl)
#ifdef H08
    wk(5) = REAL(obsda%lev(n),r_sngl) ! Him8
    wk(6) = REAL(obsda%val2(n),r_sngl) ! Him8
!    wk(7) = REAL(obsda%pred1(n),r_sngl) ! Him8
!    wk(8) = REAL(obsda%pred2(n),r_sngl) ! Him8
#endif
    WRITE(iunit) wk
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_da

subroutine write_obs_dep(cfile, nobs, set, idx, qc, omb, oma)
  implicit none
  character(*), intent(in) :: cfile
  integer, intent(in) :: nobs
  integer, intent(in) :: set(nobs)
  integer, intent(in) :: idx(nobs)
  integer, intent(in) :: qc(nobs)
  real(r_size), intent(in) :: omb(nobs)
  real(r_size), intent(in) :: oma(nobs)

  real(r_sngl) :: wk(11)
  integer :: n, iunit

  iunit=92

  open (iunit, file=cfile, form='unformatted', access='sequential')
  do n = 1, nobs
    wk(1) = real(obs(set(n))%elm(idx(n)), r_sngl)
    wk(2) = real(obs(set(n))%lon(idx(n)), r_sngl)
    wk(3) = real(obs(set(n))%lat(idx(n)), r_sngl)
    wk(4) = real(obs(set(n))%lev(idx(n)), r_sngl)
    wk(5) = real(obs(set(n))%dat(idx(n)), r_sngl)
    wk(6) = real(obs(set(n))%err(idx(n)), r_sngl)
    wk(7) = real(obs(set(n))%typ(idx(n)), r_sngl)
    wk(8) = real(obs(set(n))%dif(idx(n)), r_sngl)
    wk(9) = real(qc(n), r_sngl)
    wk(10) = real(omb(n), r_sngl)
    wk(11) = real(oma(n), r_sngl)
    select case (nint(wk(1)))
    case (id_u_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    case (id_v_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    case (id_t_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    case (id_tv_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    case (id_q_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    case (id_ps_obs)
      wk(5) = wk(5) * 0.01 ! Pa -> hPa
      wk(6) = wk(6) * 0.01 ! Pa -> hPa
    case (id_rh_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
      wk(5) = wk(5) * 100.0 ! percent output
      wk(6) = wk(6) * 100.0 ! percent output
    case (id_tcmip_obs)
      wk(5) = wk(5) * 0.01 ! Pa -> hPa
      wk(6) = wk(6) * 0.01 ! Pa -> hPa
    end select
    write (iunit) wk
  end do
  close (iunit)

  return
end subroutine write_obs_dep

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

subroutine read_obs_all(obs)
  implicit none

  type(obs_info), intent(out) :: obs(OBS_IN_NUM)
  integer :: iof
  logical :: ex

  do iof = 1, OBS_IN_NUM
    inquire (file=trim(OBS_IN_NAME(iof)), exist=ex)
    if (.not. ex) then
      write(6,*) '[Warning] FILE ',trim(OBS_IN_NAME(iof)),' NOT FOUND'


      obs(iof)%nobs = 0
      call obs_info_allocate(obs(iof), extended=.true.) !!! check why this is necessary !!!


      cycle
    end if

    select case (OBS_IN_FORMAT(iof))
    case (obsfmt_prepbufr)
      call get_nobs(trim(OBS_IN_NAME(iof)),8,obs(iof)%nobs)
    case (obsfmt_radar)
      call get_nobs_radar(trim(OBS_IN_NAME(iof)), obs(iof)%nobs, obs(iof)%meta(1), obs(iof)%meta(2), obs(iof)%meta(3))
    case (obsfmt_h08)
      if(H08_FORMAT_NC)then
        call get_nobs_allgHim8(obs(iof)%nobs) 
      else
        call get_nobs_H08(trim(OBS_IN_NAME(iof)),obs(iof)%nobs) 
      endif
    case default
      write(6,*) '[Error] Unsupported observation file format!'
      stop
    end select

    write(6,'(5A,I9,A)') 'OBS FILE [', trim(OBS_IN_NAME(iof)), '] (FORMAT ', &
                         trim(OBS_IN_FORMAT(iof)), '): TOTAL ', &
                         obs(iof)%nobs, ' OBSERVATIONS'

    call obs_info_allocate(obs(iof), extended=.true.)

    select case (OBS_IN_FORMAT(iof))
    case (obsfmt_prepbufr)
      call read_obs(trim(OBS_IN_NAME(iof)),obs(iof))
    case (obsfmt_radar)
      call read_obs_radar(trim(OBS_IN_NAME(iof)),obs(iof))
    case (obsfmt_h08)
      if(.not. H08_FORMAT_NC)then
        call read_obs_H08(trim(OBS_IN_NAME(iof)),obs(iof)) 
      endif
    end select
  end do ! [ iof = 1, OBS_IN_NUM ]

  return
end subroutine read_obs_all

subroutine write_obs_all(obs, missing, file_suffix)
  implicit none

  type(obs_info), intent(in) :: obs(OBS_IN_NUM)
  logical, intent(in), optional :: missing
  character(len=*), intent(in), optional :: file_suffix
  logical :: missing_
  integer :: iof, strlen1, strlen2
  character(200) :: filestr

  missing_ = .true.
  IF(present(missing)) missing_ = missing

  do iof = 1, OBS_IN_NUM
    if (present(file_suffix)) then
      strlen1 = len(trim(OBS_IN_NAME(iof)))
      strlen2 = len(trim(file_suffix))
      write (filestr(1:strlen1),'(A)') trim(OBS_IN_NAME(iof))
      write (filestr(strlen1+1:strlen1+strlen2),'(A)') trim(file_suffix)
    else
      filestr = OBS_IN_NAME(iof)
    end if
    select case (OBS_IN_FORMAT(iof))
    case (obsfmt_prepbufr)
      call write_obs(trim(filestr),obs(iof),missing=missing_)
    case (obsfmt_radar)
      call write_obs_radar(trim(filestr),obs(iof),missing=missing_)
    case (obsfmt_h08)
      call write_obs_H08(trim(filestr),obs(iof),missing=missing_) ! H08
    end select
  end do ! [ iof = 1, OBS_IN_NUM ]

  return
end subroutine write_obs_all
!
!-----------------------------------------------------------------------
!   TC vital obs subroutines by T. Honda (03/28/2016)
!-----------------------------------------------------------------------
!
SUBROUTINE search_tc_subdom(ritc,rjtc,v2d,yobs_tcx,yobs_tcy,yobs_mslp)
  use scale_atmos_grid_cartesC, only: &
      CXG => ATMOS_GRID_CARTESC_CXG, &
      CYG => ATMOS_GRID_CARTESC_CYG, &
      DX, &
      DY
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
  use scale_prc, only: &
      PRC_myrank

  IMPLICIT NONE
  INTEGER :: il, jl, ig, jg
  REAL(r_size) :: xdis, ydis, rdis
  REAL(r_size),INTENT(IN) :: ritc, rjtc
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  REAL(r_size),INTENT(OUT) :: yobs_mslp !(Pa)
!  REAL(r_size),INTENT(OUT) :: yobs_lon, yobs_lat !(deg)
  REAL(r_size),INTENT(OUT) :: yobs_tcx, yobs_tcy !(m)

  REAL(r_size) :: slp2d(nlonh,nlath)
  REAL(r_size) :: dz, t, q, var5

  yobs_mslp = 9.99d33
  yobs_tcx = 9.99d33
  yobs_tcy = 9.99d33

  DO jl = 1, nlat 
  DO il = 1, nlon 
    t = v2d(il,jl,iv2dd_t2m)
    q = v2d(il,jl,iv2dd_q2m)
    dz = -1.0d0 * v2d(il,jl,iv2dd_topo)
    slp2d(il,jl) = v2d(il,jl,iv2dd_ps)
    call prsadj(slp2d(il,jl),dz,t,q)
  ENDDO
  ENDDO

  DO jl = JHALO + 1, nlat - JHALO
  DO il = IHALO + 1, nlon - IHALO
    call ij_l2g(PRC_myrank, il, jl, ig, jg)
    xdis = abs(real(ig,kind=r_size) - ritc) * DX
    ydis = abs(real(jg,kind=r_size) - rjtc) * DY
    rdis = sqrt(xdis*xdis + ydis*ydis)

    IF(rdis > TC_SEARCH_DIS)CYCLE

    IF(IHALO >= 2 .and. JHALO >= 2)THEN
      call wgt_ave2d(slp2d(:,:),il,jl,var5)
    ELSE
      var5 = slp2d(il,jl)
    ENDIF

    if(var5 < yobs_mslp)then
      yobs_mslp = var5
      yobs_tcx = (real(ig,kind=r_size) - 1.0d0) * DX + CXG(1)
      yobs_tcy = (real(jg,kind=r_size) - 1.0d0) * DY + CYG(1)
    endif
  ENDDO
  ENDDO

  RETURN
END SUBROUTINE search_tc_subdom

!-- 25 points weighted average (tentative)--
! 2D weight is...
!     1 1 1 1 1
!     1 3 3 3 1
!     1 3 5 3 1
!     1 3 3 3 1
!     1 1 1 1 1
!
SUBROUTINE wgt_ave2d(var,i,j,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlonh,nlath)
  INTEGER,INTENT(IN) :: i,j
  REAL(r_size),INTENT(OUT) :: var5

  var5 = ((var(i,j) * 5.0d0 + &
         (sum(var(i-1:i+1,j-1:j+1)) - var(i,j)) * 3.0d0 + &
         (sum(var(i-2:i+2,j-2:j+2)) - sum(var(i-1:i+1,j-1:j+1))) * 1.0d0)) / 45.0d0

  RETURN
END SUBROUTINE wgt_ave2d


!
!-----------------------------------------------------------------------
!   Himawari-8 obs subroutines by T. Honda (10/29/2015)
!-----------------------------------------------------------------------
! --
!
SUBROUTINE Trans_XtoY_H08(nprof,ri,rj,lon,lat,v3d,v2d,yobs,yobs_clr,mwgt_plev,qc,zenith1d,stggrd)
  use scale_mapprojection, only: &
      MAPPROJECTION_rotcoef
  use scale_H08_fwd12
  use scale_atmos_grid_cartesC_index, only: &
      KHALO, IHALO, JHALO, &
      KS, KE, KA, KMAX
  use scale_atmos_grid_cartesC, only: &
      CZ => ATMOS_GRID_CARTESC_CZ, &
      FZ => ATMOS_GRID_CARTESC_FZ
  use scale_const, only: &
      CONST_D2R
  use scale_atmos_phy_rd_profile, only: &
      ATMOS_PHY_RD_PROFILE_read, &
      ATMOS_PHY_RD_PROFILE_setup_zgrid
  use scale_atmos_hydrometeor, only: &
      N_HYD
  use scale_atmos_aerosol, only: &
      N_AE

  IMPLICIT NONE
  INTEGER :: np, ch
  REAL(r_size),PARAMETER :: HIM8_LON = 140.7d0

  INTEGER,INTENT(IN) :: nprof ! Num of Brightness Temp "Loc" observed by Himawari-8
                              ! NOTE: multiple channels (obs) on each grid point !!
  REAL(r_size),INTENT(IN) :: ri(nprof),rj(nprof)
  REAL(r_size),INTENT(IN) :: lon(nprof),lat(nprof)
  REAL(r_size),INTENT(IN) :: v3d(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  INTEGER,INTENT(IN),OPTIONAL :: stggrd

  INTEGER :: stggrd_ = 0

! -- 2D (nlevh,nbtobs) or 1D (nbtobs) profiles for RTTOV --  
  REAL(r_size) :: prs2d(nlevh,nprof)
  REAL(r_size) :: tk2d(nlevh,nprof)
  REAL(r_size) :: qv2d(nlevh,nprof)
  REAL(r_size) :: qliq2d(nlevh,nprof)
  REAL(r_size) :: qice2d(nlevh,nprof)

  REAL(r_size) :: tsfc1d(nprof)
  REAL(r_size) :: qsfc1d(nprof)
  REAL(r_size) :: psfc1d(nprof)
  REAL(r_size) :: usfc1d(nprof)
  REAL(r_size) :: vsfc1d(nprof)
  REAL(r_size) :: lon1d(nprof)
  REAL(r_size) :: lat1d(nprof)
  REAL(r_size) :: topo1d(nprof)
  REAL(r_size) :: lsmask1d(nprof)
  REAL(r_size),INTENT(OUT) :: zenith1d(nprof) ! predictor for bias correction

! -- brightness temp from RTTOV
  REAL(r_size) :: btall_out(NIRB_HIM8,nprof) ! NOTE: RTTOV always calculates all (10) channels!!
  REAL(r_size) :: btclr_out(NIRB_HIM8,nprof) ! NOTE: RTTOV always calculates all (10) channels!!
! -- cloud top height
  REAL(r_size) :: ctop_out(nprof) 

  REAL(r_size),INTENT(OUT) :: yobs(NIRB_HIM8,nprof)
  REAL(r_size),INTENT(OUT) :: yobs_clr(NIRB_HIM8,nprof)
  INTEGER,INTENT(OUT) :: qc(NIRB_HIM8,nprof)
  REAL(r_size),INTENT(OUT) :: mwgt_plev(NIRB_HIM8,nprof)

  INTEGER :: slev, elev

  REAL(r_size) :: utmp, vtmp ! U10m & V10m tmp for rotation
  REAL(r_size),PARAMETER :: btmax = 400.0d0
  REAL(r_size),PARAMETER :: btmin = 100.0d0
  REAL(RP) :: rotc(1,1,2)
  real(r_size) :: lon_tmp(1,1),lat_tmp(1,1)

  real(r_size) :: blon, blat ! lat/lon at the domain center

  real(RP), parameter:: RD_TOA  = 100.0_RP !< top of atmosphere [km]
  integer :: RD_KMAX ! # of computational cells: z for radiation scheme

  integer, parameter :: MSTRN_ngas     =  7 !< # of gas species ! MSTRNX
  integer, parameter :: MSTRN_ncfc     = 28 !< # of CFC species ! MSTRNX

  integer, parameter :: ngas = MSTRN_ngas
  integer, parameter :: ncfc = MSTRN_ncfc
  integer, parameter :: RD_naero      = N_HYD + N_AE ! # of cloud/aerosol species

  real(RP), allocatable :: RD_zh          (:)   ! altitude    at the interface [km]
  real(RP), allocatable :: RD_z           (:)   ! altitude    at the center [km]
  real(RP), allocatable :: RD_rhodz       (:)   ! density * delta z [kg/m2]
  real(RP), allocatable :: RD_pres        (:)   ! pressure    at the center [hPa]
  real(RP), allocatable :: RD_presh       (:)   ! pressure    at the interface [hPa]
  real(RP), allocatable :: RD_temp        (:)   ! temperature at the center [K]
  real(RP), allocatable :: RD_temph       (:)   ! temperature at the interface [K]
  real(RP), allocatable :: RD_gas         (:,:) ! gas species   volume mixing ratio [ppmv]
  real(RP), allocatable :: RD_cfc         (:,:) ! CFCs          volume mixing ratio [ppmv]
  real(RP), allocatable :: RD_aerosol_conc(:,:) ! cloud/aerosol volume mixing ratio [ppmv]
  real(RP), allocatable :: RD_aerosol_radi(:,:) ! cloud/aerosol effective radius [cm]
  real(RP), allocatable :: RD_cldfrac     (:)   ! cloud fraction (0-1)


  !
  ! Extrapolate input profiles by using climatology (MIPAS)
  ! Based on "scalelib/src/atmos-physics/scale_atmos_phy_rd_mstrnx.F90"
  !

  ! Get basepoint lat/lon
  call ij2phys(real(nlong/2+IHALO, kind=r_size),&
               real(nlatg/2+JHALO, kind=r_size),&
               blon, blat)

  ! --- setup MSTRN parameter
  !call RD_MSTRN_setup( ngas, & ! [OUT]
  !                     ncfc  ) ! [OUT]

  !--- setup climatological profile
  !    Done from common_mpi_scale

  !--- setup climatological profile
  !    Done from common_mpi_scale

  RD_KMAX      = KMAX + H08_RTTOV_KADD

  !--- allocate arrays
  ! input
  allocate( RD_zh   (RD_KMAX+1) )
  allocate( RD_z    (RD_KMAX  ) )

  allocate( RD_rhodz(RD_KMAX  ) )
  allocate( RD_pres (RD_KMAX  ) )
  allocate( RD_presh(RD_KMAX+1) )
  allocate( RD_temp (RD_KMAX  ) )
  allocate( RD_temph(RD_KMAX+1) )

  allocate( RD_gas         (RD_KMAX,ngas    ) )
  allocate( RD_cfc         (RD_KMAX,ncfc    ) )
  allocate( RD_aerosol_conc(RD_KMAX,RD_naero) )
  allocate( RD_aerosol_radi(RD_KMAX,RD_naero) )
  allocate( RD_cldfrac     (RD_KMAX         ) )

  !--- setup vartical grid for radiation (larger TOA than Model domain)
  call ATMOS_PHY_RD_PROFILE_setup_zgrid( KA, KS, KE, RD_KMAX, H08_RTTOV_KADD, & ! [IN]
                                         RD_TOA, CZ, FZ, & ! [IN]
                                         RD_zh(:), RD_z(:)         ) ! [INOUT]

  !--- read climatological profile
  call ATMOS_PHY_RD_PROFILE_read( RD_KMAX,                & ! [IN]
                                  ngas,                   & ! [IN]
                                  ncfc,                   & ! [IN]
                                  RD_naero,               & ! [IN]
                                  blat*CONST_D2R,         & ! [IN]
                                  H08_NOWDATE    (:),     & ! [IN]
                                  RD_zh          (:),     & ! [IN]
                                  RD_z           (:),     & ! [IN]
                                  RD_rhodz       (:),     & ! [OUT]
                                  RD_pres        (:),     & ! [OUT]
                                  RD_presh       (:),     & ! [OUT]
                                  RD_temp        (:),     & ! [OUT]
                                  RD_temph       (:),     & ! [OUT]
                                  RD_gas         (:,:),   & ! [OUT]
                                  RD_cfc         (:,:),   & ! [OUT]
                                  RD_aerosol_conc(:,:),   & ! [OUT]
                                  RD_aerosol_radi(:,:),   & ! [OUT]
                                  RD_cldfrac     (:)      ) ! [OUT]

!  kidx_rlx = 1
!  do k = 1, RD_KMAX + 1
!    if(RD_zh(k)*1.0d3 < H08_RTTOV_RLX_HGT)then
!      kidx_rlx = k - RD_KADD
!      exit
!    endif
!  enddo


  if (present(stggrd)) stggrd_ = stggrd

  slev = 1 + KHALO
  elev = nlevh - KHALO

! -- make profile arrays for RTTOV --
  DO np = 1, nprof ! -- make profiles

    lon1d(np) = lon(np)
    lat1d(np) = lat(np)

    CALL zenith_geosat(HIM8_LON,lon(np),lat(np),zenith1d(np))
    CALL itpl_2d(v2d(:,:,iv2dd_skint),ri(np),rj(np),tsfc1d(np)) ! T2 is better??
!    CALL itpl_2d(v2d(:,:,iv2dd_t2m),ri(np),rj(np),tsfc1d(np))
    CALL itpl_2d(v2d(:,:,iv2dd_q2m),ri(np),rj(np),qsfc1d(np))
!    CALL itpl_2d(v3d(KHALO+1,:,:,iv3dd_hgt),ri(np),rj(np),topo1d(np))
    CALL itpl_2d(v2d(:,:,iv2dd_topo),ri(np),rj(np),topo1d(np))
    CALL itpl_2d(v2d(:,:,iv2dd_lsmask),ri(np),rj(np),lsmask1d(np))
    CALL itpl_2d(v2d(:,:,iv2dd_ps),ri(np),rj(np),psfc1d(np))

!    call prsadj(yobs,rk-topo,t,q)
!    if (abs(rk-topo) > PS_ADJUST_THRES) then
!      write (6,'(A,F6.1)') '[Warning] PS observation height adjustment exceeds the threshold. dz=', abs(rk-topo)
!      qc = iqc_ps_ter
!    end if

    if (stggrd_ == 1) then
      CALL itpl_2d(v2d(:,:,iv2dd_u10m),ri(np)-0.5,rj(np),utmp)  !###### should modity itpl_3d to prevent '1.0' problem....??
      CALL itpl_2d(v2d(:,:,iv2dd_v10m),ri(np),rj(np)-0.5,vtmp)  !######
    else
      CALL itpl_2d(v2d(:,:,iv2dd_u10m),ri(np),rj(np),utmp)
      CALL itpl_2d(v2d(:,:,iv2dd_v10m),ri(np),rj(np),vtmp)
    end if

    lon_tmp(1,1) = lon(np)*deg2rad
    lat_tmp(1,1) = lat(np)*deg2rad
    call MAPPROJECTION_rotcoef(1, 1, 1, 1, 1, 1, &
                               lon_tmp(1,1),lat_tmp(1,1),rotc)
    usfc1d(np) = utmp * rotc(1,1,1) - vtmp * rotc(1,1,2)
    vsfc1d(np) = utmp * rotc(1,1,2) + vtmp * rotc(1,1,1)

    CALL itpl_2d_column(v3d(:,:,:,iv3dd_p),ri(np),rj(np),prs2d(:,np))
    CALL itpl_2d_column(v3d(:,:,:,iv3dd_t),ri(np),rj(np),tk2d(:,np))
    CALL itpl_2d_column(v3d(:,:,:,iv3dd_q),ri(np),rj(np),qv2d(:,np))
    CALL itpl_2d_column(v3d(:,:,:,iv3dd_qc),ri(np),rj(np),qliq2d(:,np))
    CALL itpl_2d_column((v3d(:,:,:,iv3dd_qi) &
                       + v3d(:,:,:,iv3dd_qs) &
                       + v3d(:,:,:,iv3dd_qg)),ri(np),rj(np),qice2d(:,np))

  ENDDO ! -- make profiles


!
! -- NOTE: The channel number for RTTOV is always 10, because it should be the same
!          with that in Himawari-8 RTTOV coef files.
!
!        : Satellite zenith angles are computed within SCALE_RTTOV_fwd using (lon,lat).
!

  CALL SCALE_RTTOV_fwd12(NIRB_HIM8, & ! num of channels
                       nlev,& ! num of levels
                       nprof,& ! num of profs
                       prs2d(elev:slev:-1,1:nprof),& ! (Pa)
                       tk2d(elev:slev:-1,1:nprof),& ! (K)
                       qv2d(elev:slev:-1,1:nprof),& ! (kg/kg)
                       qliq2d(elev:slev:-1,1:nprof),& ! (kg/kg)
                       qice2d(elev:slev:-1,1:nprof),& ! (kg/kg)
                       tsfc1d(1:nprof),& ! (K)
                       qsfc1d(1:nprof),& ! (kg/kg)
                       psfc1d(1:nprof),& ! (Pa)
                       usfc1d(1:nprof),& ! (m/s)
                       vsfc1d(1:nprof),& ! (m/s)
                       topo1d(1:nprof),& ! (m)
                       lon1d(1:nprof),& ! (deg)
                       lat1d(1:nprof),& ! (deg)
                       lsmask1d(1:nprof),& ! (0-1)
                       zenith1d(1:nprof), & ! (deg) 
                       RD_presh(1:RD_KMAX+1), & ! (hPa) 
                       RD_temph(1:RD_KMAX+1), & ! (K) 
                       btall_out(1:NIRB_HIM8,1:nprof),& ! (K)
                       btclr_out(1:NIRB_HIM8,1:nprof),& ! (K)
                       mwgt_plev(1:NIRB_HIM8,1:nprof),& ! (Pa)
                       ctop_out(1:nprof))


  deallocate(RD_zh,RD_z,RD_rhodz,RD_pres,RD_presh)
  deallocate(RD_temp,RD_temph,RD_gas,RD_cfc,RD_cldfrac)
  deallocate(RD_aerosol_conc,RD_aerosol_radi)

!
! -- btall_out is substituted into yobs
!
  do np = 1, nprof
  do ch = 1, NIRB_HIM8

    yobs(ch,np) = btall_out(ch,np)
    qc(ch,np) = iqc_good

    ! Band QC will be done in lekt_obs
    !if (H08_BAND_USE(ch) /= 1) then
    !  qc(ch,np) = iqc_obs_bad
    !end if

    if(H08_VLOCAL_CTOP)then
      if((ctop_out(np) > 0.0d0) .and. (ctop_out(np) < mwgt_plev(ch,np)) .and. &
         (mwgt_plev(ch,np)>H08_LIMIT_LEV)) then
        mwgt_plev(ch,np) = (ctop_out(np) + mwgt_plev(ch,np))*0.5d0
      endif
    endif

    if(H08_REJECT_LAND .and. (lsmask1d(np) > 0.5d0))then
      qc(ch,np) = iqc_obs_bad
    endif

    ! QC
    if(yobs(ch,np) > btmax .or. yobs(ch,np) < btmin .or. yobs(ch,np) /= yobs(ch,np))then
      qc(ch,np) = iqc_obs_bad
    endif

    yobs_clr(ch,np) = btclr_out(ch,np)

  enddo ! ch
  enddo ! np

  return
END SUBROUTINE Trans_XtoY_H08

!
SUBROUTINE Trans_XtoY_H08_allg(v3d,v2d,yobs,yobs_clr,mwgt_plev2d,qc,zenith1d,stggrd)
  use scale_mapprojection, only: &
      MAPPROJECTION_rotcoef, &
      MAPPROJECTION_xy2lonlat
  use scale_H08_fwd12
  use scale_atmos_grid_cartesC_index, only: &
      KHALO, IHALO, JHALO, &
      KS, KE, KA, KMAX
  use scale_atmos_grid_cartesC, only: &
      CZ  => ATMOS_GRID_CARTESC_CZ, &
      FZ  => ATMOS_GRID_CARTESC_FZ, &
      CX => ATMOS_GRID_CARTESC_CX, &
      CY => ATMOS_GRID_CARTESC_CY, &
      DX, DY
  use scale_const, only: &
      CONST_D2R
  use scale_atmos_phy_rd_profile, only: &
      ATMOS_PHY_RD_PROFILE_read, &
      ATMOS_PHY_RD_PROFILE_setup_zgrid
  use scale_atmos_hydrometeor, only: &
      N_HYD
  use scale_atmos_aerosol, only: &
      N_AE

  IMPLICIT NONE
  INTEGER :: np, ch
  REAL(r_size),PARAMETER :: HIM8_LON = 140.7d0

  REAL(r_size),INTENT(IN) :: v3d(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  INTEGER,INTENT(IN),OPTIONAL :: stggrd
  REAL(RP) :: rotc(1,1,2)

  INTEGER :: stggrd_ = 0

! -- 2D (nlevh,nbtobs) or 1D (nbtobs) profiles for RTTOV --  
  REAL(r_size) :: prs2d(nlev,nlon*nlat)
  REAL(r_size) :: tk2d(nlev,nlon*nlat)
  REAL(r_size) :: qv2d(nlev,nlon*nlat)
  REAL(r_size) :: qliq2d(nlev,nlon*nlat)
  REAL(r_size) :: qice2d(nlev,nlon*nlat)

  REAL(r_size) :: tsfc1d(nlon*nlat)
  REAL(r_size) :: qsfc1d(nlon*nlat)
  REAL(r_size) :: psfc1d(nlon*nlat)
  REAL(r_size) :: usfc1d(nlon*nlat)
  REAL(r_size) :: vsfc1d(nlon*nlat)
  REAL(r_size) :: lon1d(nlon*nlat)
  REAL(r_size) :: lat1d(nlon*nlat)
  REAL(r_size) :: topo1d(nlon*nlat)
  REAL(r_size) :: lsmask1d(nlon*nlat)
  REAL(r_size),INTENT(OUT) :: zenith1d(nlon*nlat) ! predictor for bias correction

! -- brightness temp from RTTOV
  REAL(r_size) :: btall_out(NIRB_HIM8,nlon*nlat) ! NOTE: RTTOV always calculates all (10) channels!!
  REAL(r_size) :: btclr_out(NIRB_HIM8,nlon*nlat) ! NOTE: RTTOV always calculates all (10) channels!!
! -- cloud top height
  REAL(r_size) :: ctop_out1d(nlon*nlat) 

  REAL(r_size),INTENT(OUT) :: yobs(nlon,nlat,NIRB_HIM8)
  REAL(r_size),INTENT(OUT) :: yobs_clr(nlon,nlat,NIRB_HIM8)
  REAL(r_size),INTENT(OUT) :: mwgt_plev2d(nlon,nlat,NIRB_HIM8)
  INTEGER,INTENT(OUT) :: qc(nlon,nlat,NIRB_HIM8)
  REAL(r_size) :: mwgt_plev1d(NIRB_HIM8,nlon*nlat)

  REAL(r_size) :: utmp, vtmp ! U10m & V10m tmp for rotation
  real(r_size) :: lon_tmp(1,1),lat_tmp(1,1)
  REAL(r_size),PARAMETER :: btmax = 400.0d0
  REAL(r_size),PARAMETER :: btmin = 100.0d0

  real(r_size) :: blon, blat ! lat/lon at the domain center
  integer :: k

  real(RP), parameter:: RD_TOA  = 50.0_RP !< top of atmosphere [km]

  integer, parameter :: MSTRN_ngas     =  7 !< # of gas species ! MSTRNX
  integer, parameter :: MSTRN_ncfc     = 28 !< # of CFC species ! MSTRNX

  integer, parameter :: ngas = MSTRN_ngas
  integer, parameter :: ncfc = MSTRN_ncfc
  integer, parameter :: RD_naero      = N_HYD + N_AE ! # of cloud/aerosol species

  integer :: RD_KMAX ! # of computational cells: z for radiation scheme
  real(RP) :: RD_zh          (KMAX + H08_RTTOV_KADD + 1)   ! altitude    at the interface [km]
  real(RP) :: RD_z           (KMAX + H08_RTTOV_KADD)   ! altitude    at the center [km]
  real(RP) :: RD_rhodz       (KMAX + H08_RTTOV_KADD)   ! density * delta z [kg/m2]
  real(RP) :: RD_pres        (KMAX + H08_RTTOV_KADD)   ! pressure    at the center [hPa]
  real(RP) :: RD_presh       (KMAX + H08_RTTOV_KADD + 1)   ! pressure    at the interface [hPa]
  real(RP) :: RD_temp        (KMAX + H08_RTTOV_KADD)   ! temperature at the center [K]
  real(RP) :: RD_temph       (KMAX + H08_RTTOV_KADD + 1)   ! temperature at the interface [K]
  real(RP) :: RD_gas         (KMAX + H08_RTTOV_KADD,ngas) ! gas species   volume mixing ratio [ppmv]
  real(RP) :: RD_cfc         (KMAX + H08_RTTOV_KADD,ncfc) ! CFCs          volume mixing ratio [ppmv]
  real(RP) :: RD_aerosol_conc(KMAX + H08_RTTOV_KADD,RD_naero) ! cloud/aerosol volume mixing ratio [ppmv]
  real(RP) :: RD_aerosol_radi(KMAX + H08_RTTOV_KADD,RD_naero) ! cloud/aerosol effective radius [cm]
  real(RP) :: RD_cldfrac     (KMAX + H08_RTTOV_KADD)   ! cloud fraction (0-1)

  integer :: i, j
  real(r_size) :: ri, rj

  !
  ! Extrapolate input profiles by using climatology (MIPAS)
  ! Based on "scalelib/src/atmos-physics/scale_atmos_phy_rd_mstrnx.F90"
  !

  ! Get basepoint lat/lon
  call ij2phys(real(nlong/2+IHALO, kind=r_size),&
               real(nlatg/2+JHALO, kind=r_size),&
               blon, blat)

  RD_KMAX = KMAX + H08_RTTOV_KADD

  !--- setup vartical grid for radiation (larger TOA than Model domain)
  call ATMOS_PHY_RD_PROFILE_setup_zgrid( KA, KS, KE, RD_KMAX, H08_RTTOV_KADD, & ! [IN]
                                         RD_TOA, CZ, FZ, & ! [IN]
                                         RD_zh(:), RD_z(:)         ) ! [INOUT]

  !--- read climatological profile
  call ATMOS_PHY_RD_PROFILE_read( KMAX + H08_RTTOV_KADD,  & ! [IN]
                                  ngas,                   & ! [IN]
                                  ncfc,                   & ! [IN]
                                  RD_naero,               & ! [IN]
                                  blat*CONST_D2R,         & ! [IN]
                                  H08_NOWDATE    (:),     & ! [IN]
                                  RD_zh          (:),     & ! [IN]
                                  RD_z           (:),     & ! [IN]
                                  RD_rhodz       (:),     & ! [OUT]
                                  RD_pres        (:),     & ! [OUT]
                                  RD_presh       (:),     & ! [OUT]
                                  RD_temp        (:),     & ! [OUT]
                                  RD_temph       (:),     & ! [OUT]
                                  RD_gas         (:,:),   & ! [OUT]
                                  RD_cfc         (:,:),   & ! [OUT]
                                  RD_aerosol_conc(:,:),   & ! [OUT]
                                  RD_aerosol_radi(:,:),   & ! [OUT]
                                  RD_cldfrac     (:)      ) ! [OUT]

  if (present(stggrd)) stggrd_ = stggrd

! -- make profile arrays for RTTOV --
  do j = 1, nlat
  do i = 1, nlon
    np = (j - 1) * nlon + i

    ri = real(i + IHALO, r_size)
    rj = real(j + JHALO, r_size)
    call MAPPROJECTION_xy2lonlat((ri-1.0_r_size) * DX + CX(1), &
                                 (rj-1.0_r_size) * DY + CY(1),&
                                 lon1d(np), lat1d(np))

    lon1d(np) = lon1d(np) * rad2deg
    lat1d(np) = lat1d(np) * rad2deg

    CALL zenith_geosat(HIM8_LON,lon1d(np),lat1d(np),zenith1d(np))
    tsfc1d(np) = v2d(i+IHALO,j+JHALO,iv2dd_skint)
    qsfc1d(np) = v2d(i+IHALO,j+JHALO,iv2dd_q2m)
    topo1d(np) = v2d(i+IHALO,j+JHALO,iv2dd_topo)
    lsmask1d(np) = v2d(i+IHALO,j+JHALO,iv2dd_lsmask)
    psfc1d(np) = v2d(i+IHALO,j+JHALO,iv2dd_ps)

    ! assume not staggerd grid
    utmp = v2d(i+IHALO,j+JHALO,iv2dd_u10m)
    vtmp = v2d(i+IHALO,j+JHALO,iv2dd_v10m)
    call MAPPROJECTION_rotcoef(1, 1, 1, 1, 1, 1, &
                               lon1d(np),lat1d(np),rotc)
    usfc1d(np) = utmp * rotc(1,1,1) - vtmp * rotc(1,1,2)
    vsfc1d(np) = utmp * rotc(1,1,2) + vtmp * rotc(1,1,1)

    do k = 1, KMAX
      prs2d(k,np) = v3d(KHALO+KMAX-k+1,i+IHALO,j+JHALO,iv3dd_p)
      tk2d(k,np) = v3d(KHALO+KMAX-k+1,i+IHALO,j+JHALO,iv3dd_t)
      qv2d(k,np) = v3d(KHALO+KMAX-k+1,i+IHALO,j+JHALO,iv3dd_q)
      qliq2d(k,np) = v3d(KHALO+KMAX-k+1,i+IHALO,j+JHALO,iv3dd_qc)
      qice2d(k,np) = v3d(KHALO+KMAX-k+1,i+IHALO,j+JHALO,iv3dd_qi) & 
                   + v3d(KHALO+KMAX-k+1,i+IHALO,j+JHALO,iv3dd_qs) &
                   + v3d(KHALO+KMAX-k+1,i+IHALO,j+JHALO,iv3dd_qg)  

    enddo


  enddo ! i
  enddo ! j

!
! -- NOTE: The channel number for RTTOV is always 10, because it should be the same
!          with that in Himawari-8 RTTOV coef files.
!
!        : Satellite zenith angles are computed within SCALE_RTTOV_fwd using (lon,lat).
!

write(6,'(a)') "DEBUG before RTTOV"

  CALL SCALE_RTTOV_fwd12(NIRB_HIM8, & ! num of channels
                       KMAX,& ! num of levels
                       nlon*nlat,& ! num of profs
                       prs2d(:,:),& ! (Pa)
                       tk2d(:,:),& ! (K)
                       qv2d(:,:),& ! (kg/kg)
                       qliq2d(:,:),& ! (kg/kg)
                       qice2d(:,:),& ! (kg/kg)
                       tsfc1d(:),& ! (K)
                       qsfc1d(:),& ! (kg/kg)
                       psfc1d(:),& ! (Pa)
                       usfc1d(:),& ! (m/s)
                       vsfc1d(:),& ! (m/s)
                       topo1d(:),& ! (m)
                       lon1d(:),& ! (deg)
                       lat1d(:),& ! (deg)
                       lsmask1d(:),& ! (0-1)
                       zenith1d(:), & ! (deg) 
                       RD_presh(:), & ! (hPa) 
                       RD_temph(:), & ! (K) 
                       btall_out(:,:),& ! (K)
                       btclr_out(:,:),& ! (K)
                       mwgt_plev1d(:,:),& ! (Pa)
                       ctop_out1d(:))

write(6,'(a)') "DEBUG after RTTOV"

!
! -- btall_out is substituted into yobs
!

  do j = 1, nlat
  do i = 1, nlon
    np = (j - 1) * nlon + i

    do ch = 1, NIRB_HIM8
      qc(i,j,ch) = iqc_good
      yobs(i,j,ch) = btall_out(ch,np)
      yobs_clr(i,j,ch) = btclr_out(ch,np)
      mwgt_plev2d(i,j,ch) = mwgt_plev1d(ch,np)

      if(H08_VLOCAL_CTOP)then
        if((ctop_out1d(np) > 0.0d0) .and. (ctop_out1d(np) < mwgt_plev1d(ch,np)) .and. &
           (mwgt_plev1d(ch,np)>H08_LIMIT_LEV)) then
          mwgt_plev1d(ch,np) = (ctop_out1d(np) + mwgt_plev1d(ch,np))*0.5d0
        endif
      endif

      ! QC
      if(H08_REJECT_LAND .and. (lsmask1d(np) > 0.5d0))then
        qc(i,j,ch) = iqc_obs_bad
      endif

      if(yobs(i,j,ch) > btmax .or. yobs(i,j,ch) < btmin .or. yobs(i,j,ch) /= yobs(i,j,ch))then
        qc(i,j,np) = iqc_obs_bad
      endif

    enddo ! ch

  enddo ! i
  enddo ! j

  return
END SUBROUTINE Trans_XtoY_H08_allg

SUBROUTINE write_vbc_Him8(vbca,ANAL)
  implicit none

  real(r_size), intent(in) :: vbca(H08_NPRED,NIRB_HIM8)
  logical :: ANAL
  integer :: npr, ich
  character(255) :: OUTFILE
  integer,parameter :: iunit = 99

  if(ANAL) then
    OUTFILE = trim(H08_VBC_PATH)//'/Him8_vbca.dat'
  else
    OUTFILE = trim(H08_VBC_PATH)//'/Him8_vbcf.dat'
  endif


  open(iunit,file=trim(OUTFILE),status='unknown',access='direct',&
       form='unformatted',recl=NIRB_HIM8*8)

  do npr = 1, H08_NPRED    
    write(iunit,rec=npr)(vbca(npr,ich),ich=1,NIRB_HIM8)
  enddo

  do npr = 1, H08_NPRED    
    write(6,'(a,f15.10)')"DEBUG VBCA",vbca(npr,3)
  enddo

  close(iunit)

  return
END SUBROUTINE write_vbc_Him8


SUBROUTINE read_vbc_Him8(vbc)
  implicit none

  real(r_size), intent(out) :: vbc(H08_NPRED,NIRB_HIM8)
  integer :: npr, ich
  integer,parameter :: iunit = 99
  logical :: ex

  inquire(file=trim(H08_VBC_PATH)//'/Him8_vbcf.dat', exist=ex)

  if(ex)then
    open(iunit,file=trim(H08_VBC_PATH)//'/Him8_vbcf.dat',status='unknown',access='direct',&
         form='unformatted',recl=NIRB_HIM8*8)

    do npr = 1, H08_NPRED
      read(iunit,rec=npr)(vbc(npr,ich),ich=1,NIRB_HIM8)
    enddo
    close(iunit)
  else
    write(6,'(a)')' xxx Failed to open vbc_Him8'
    vbc = 0.0d0
  endif

  return
END SUBROUTINE read_vbc_Him8


SUBROUTINE zenith_geosat(sat_lon,lon,lat,z_angle)
! 
! Compute geostatinoary-satelitte zenith angle from lat/lon information
!
! -- Note: Computation of the zenith angle in each obs point (P) is based on the
! formula in
!          LRIT/HRIT Global Specification.
!          http://www.cgms-info.org/index_.php/cgms/page?cat=publications&page=technical+publications
! 
  USE scale_const, ONLY: &
      Deg2Rad => CONST_D2R

  IMPLICIT NONE 

  REAL(r_size),INTENT(IN) :: sat_lon ! longitude of Himawari-8 satellite
  REAL(r_size),INTENT(IN) :: lon, lat ! (degree)
  REAL(r_size),INTENT(OUT) :: z_angle ! zenith angle

  REAL(r_size),PARAMETER :: Rpol = 6356.7523d3 ! a polar radius of Earth (m) 
  REAL(r_size) :: Rl ! a local radius of Earth
  REAL(r_size) :: rlon, rlat ! (Radian)
!
!
! Vector components for a satellite coordinate frame
!
!
  REAL(r_size) :: rnps, rnep, c_lat ! auxiliary variables
  REAL(r_size) :: r1, r2, r3       ! components of location vector for point P 
  REAL(r_size) :: r1ps, r2ps, r3ps ! components of the vector from P to the satellite 
  REAL(r_size) :: r1ep, r2ep, r3ep  ! components of the vector from the center of Earth to P

! sattelite zenith angle 

  rlat = lat * Deg2Rad
  rlon = lon * Deg2Rad

  c_lat = datan(0.993305616d0 * dtan(rlat))
  Rl = Rpol / dsqrt(1.0d0 - 0.00669438444d0 * dcos(c_lat)*dcos(c_lat))
  r1 = 42164.0d3 - Rl * dcos(c_lat) * dcos(rlon - sat_lon*Deg2Rad)
  r2 = -Rl * dcos(c_lat) * dsin(rlon - sat_lon*Deg2Rad)
  r3 = Rl * dsin(c_lat)
  rnps = dsqrt(r1*r1+r2*r2+r3*r3)

  r1ps = r1 * (-1.0d0)
  r2ps = r2 * (-1.0d0)
  r3ps = r3 * (-1.0d0)

  r1ep = r1 - 42164.0d3
  r2ep = r2
  r3ep = r3

  rnep = dsqrt(r1ep*r1ep+r2ep*r2ep+r3ep*r3ep)

  z_angle = r1ps * r1ep + r2ps * r2ep + r3ps * r3ep ! internal product 
  z_angle = dacos(z_angle/(rnps*rnep))/Deg2Rad


  RETURN
END SUBROUTINE zenith_geosat

SUBROUTINE get_nobs_H08(cfile,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn ! num of all H08 obs
!  REAL(r_sngl) :: wk(4+NIRB_HIM8)
!  INTEGER :: ios 
  INTEGER :: iprof
  INTEGER :: iunit
  LOGICAL :: ex
  INTEGER :: sz

  nn = 0 
  iprof = 0
  iunit=91
      
      
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    
!    READ(iunit,IOSTAT=ios)wk
!    IF(ios /= 0) THEN 
!      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
!      RETURN
!    END IF
    
! get file size by reading through the entire file... too slow for big files
!-----------------------------
!    DO
!      READ(iunit,IOSTAT=ios) wk
!      IF(ios /= 0) EXIT
!      iprof = iprof + 1
!      nn = nn + NIRB_HIM8
!    END DO
!-----------------------------

! get file size by INQUIRE statement... may not work for some older fortran compilers
!-----------------------------
    INQUIRE(UNIT=iunit, SIZE=sz)
    IF (MOD(sz, r_sngl * (H08_OBS_RECL+2)) /= 0) THEN
      WRITE(6,'(2A)') cfile,': Reading error -- skipped'
      RETURN
    END IF
    iprof = sz / (r_sngl * (H08_OBS_RECL+2))
    nn = iprof * NIRB_HIM8
!-----------------------------

    WRITE(6,*)' H08 FILE ', cfile
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    WRITE(6,'(A12,I10)') '   num of prof:',iprof
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs_H08

SUBROUTINE read_obs_H08(cfile,obs)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(INOUT) :: obs
  REAL(r_sngl) :: wk(H08_OBS_RECL)
  INTEGER :: n,iunit

  INTEGER :: nprof, np, ch
  INTEGER :: istd, idif

  if(H08_OBS_STD)then
    istd = 4 + NIRB_HIM8 + 1
    idif = 4 + NIRB_HIM8 + 2
  else
    idif = 4 + NIRB_HIM8 + 1
  endif


  nprof = obs%nobs / NIRB_HIM8

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

  n = 0
  DO np=1,nprof
    READ(iunit) wk

    DO ch = 1, NIRB_HIM8
      n = n + 1

      obs%elm(n) = NINT(wk(1))
      obs%typ(n) = NINT(wk(2))
      obs%lon(n) = REAL(wk(3),r_size)
      obs%lat(n) = REAL(wk(4),r_size)
      obs%dat(n) = REAL(wk(4+ch),r_size)
      obs%lev(n) = ch + 6.0 ! substitute channnel number instead of the obs level
      !obs%err(n) = REAL(OBSERR_H08(ch),r_size)

      if(H08_OBS_STD)then
        obs%err(n) = -REAL(wk(istd),r_size)
      else
        obs%err(n) = REAL(OBSERR_H08(ch),r_size)
      endif

      if(H08_OBS_4D)then
        obs%dif(n) = -REAL(wk(idif),r_size)
      else
        obs%dif(n) = 0.0d0
      endif
    END DO
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_H08

SUBROUTINE write_obs_H08(cfile,obs,append,missing)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(IN) :: obs
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL,INTENT(IN),OPTIONAL :: missing
  LOGICAL :: append_
  LOGICAL :: missing_
  REAL(r_sngl) :: wk(4+NIRB_HIM8)
  INTEGER :: n,iunit
  INTEGER :: iprof, ns, ne

  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  missing_ = .true.
  IF(present(missing)) missing_ = missing

  iprof = obs%nobs / NIRB_HIM8


  IF(append_) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  END IF

  DO n=1,iprof
    ns = (n-1)*NIRB_HIM8 + 1
    ne = n*NIRB_HIM8

    wk(1) = REAL(obs%elm(ns),r_sngl)
    wk(2) = REAL(obs%typ(ns),r_sngl)
    wk(3) = REAL(obs%lon(ns),r_sngl)
    wk(4) = REAL(obs%lat(ns),r_sngl)
    wk(5:5+NIRB_HIM8-1) = REAL(obs%dat(ns:ne),r_size)
    WRITE(iunit) wk

  ENDDO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_H08

subroutine get_dim_Him8_nc(filename,imax_him8,jmax_him8)
  use netcdf
  use common_ncio
  implicit none

  character(*),intent(in) :: filename
!  character(len=12) :: filesuffix = '.nc'

  integer :: ncid, varid
!  integer, intent(out) :: nobs
  integer,intent(out) :: imax_him8, jmax_him8

  call ncio_open(trim(filename), NF90_NOWRITE, ncid)

  call ncio_read_dim(ncid,'longitude',imax_him8)
  call ncio_read_dim(ncid,'latitude',jmax_him8)

  call ncio_close(ncid)

  return
end subroutine get_dim_Him8_nc

subroutine read_Him8_nc(filename,imax_him8,jmax_him8,lon_him8,lat_him8,tbb)
  use netcdf
  use common_ncio
  implicit none

  character(*),intent(in) :: filename

  integer :: ncid, varid
  integer,intent(in) :: imax_him8, jmax_him8

  real(r_sngl),intent(out) :: lon_him8(imax_him8)
  real(r_sngl),intent(out) :: lat_him8(jmax_him8)
  real(r_sngl),intent(out) :: tbb(imax_him8,jmax_him8,NIRB_HIM8)

  call ncio_open(trim(filename), NF90_NOWRITE, ncid)

  call ncio_check(nf90_inq_varid(ncid, 'longitude', varid))
  call ncio_check(nf90_get_var(ncid, varid, lon_him8, &
                               start = (/ 1 /), count = (/ imax_him8 /)))

  call ncio_check(nf90_inq_varid(ncid, 'latitude', varid))
  call ncio_check(nf90_get_var(ncid, varid, lat_him8, &
                               start = (/ 1 /), count = (/ jmax_him8 /)))

  call ncio_check(nf90_inq_varid(ncid, 'tbb', varid))
  call ncio_check(nf90_get_var(ncid, varid, tbb, &
                               start = (/ 1, 1, 1 /), count = (/ imax_him8, jmax_him8, NIRB_HIM8 /)))

  call ncio_close(ncid)

  return
end subroutine read_Him8_nc

subroutine sobs_Him8(imax_him8,jmax_him8,lon_him8,lat_him8,tbb_org,tbb_sobs)
  use scale_atmos_grid_cartesC, only: &
      CX => ATMOS_GRID_CARTESC_CX, &
      CY => ATMOS_GRID_CARTESC_CY, &
      DX, DY
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  implicit none

  integer,intent(in) :: imax_him8, jmax_him8

  real(r_sngl),intent(in) :: lon_him8(imax_him8)
  real(r_sngl),intent(in) :: lat_him8(jmax_him8)
  real(r_sngl),intent(in) :: tbb_org(imax_him8,jmax_him8,NIRB_HIM8)

  real(r_size),intent(out) :: tbb_sobs(nlon,nlat,NIRB_HIM8)

  real(r_size) :: ri, rj  
  real(r_size) :: ri_tmp(2), rj_tmp(2)  
  real(r_size) :: rlon_tmp(2), rlat_tmp(2)
  real(r_size) :: lon2d(nlon,nlat), lat2d(nlon,nlat)

  integer :: i, j, k, ii, jj
  integer :: is, ie, js, je
  integer :: dix, diy, cnt
  integer :: i_him8, j_him8

  ! Assumte that Himawari-8 obs is based on a uniform lat-lon coordinate
  !
  ! Use a distance btw. the center of a subdomain & an adjacent grid point
  do i = 1, 2
    ri_tmp(i) = real(int(nlon/2) + i - 1 + IHALO, r_size) 
    rj_tmp(i) = real(int(nlat/2) + i - 1 + JHALO, r_size)
    call MAPPROJECTION_xy2lonlat((ri_tmp(i)-1.0_r_size) * DX + CX(1), &
                                 (rj_tmp(i)-1.0_r_size) * DY + CY(1),&
                                  rlon_tmp(i), rlat_tmp(i))
    rlon_tmp(i) = rlon_tmp(i) * rad2deg
    rlat_tmp(i) = rlat_tmp(i) * rad2deg
  enddo

  dix = max(nint(abs(rlon_tmp(2) - rlon_tmp(1)) * 0.5d0 / abs(lon_him8(2)-lon_him8(1))),1)
  diy = max(nint(abs(rlat_tmp(2) - rlat_tmp(1)) * 0.5d0 / abs(lat_him8(2)-lat_him8(1))),1)

  tbb_sobs = 0.0d0

  do j = 1, nlat
  do i = 1, nlon

    ri = real(i + IHALO, r_size)
    rj = real(j + JHALO, r_size)
    call MAPPROJECTION_xy2lonlat((ri-1.0_r_size) * DX + CX(1), &
                                 (rj-1.0_r_size) * DY + CY(1),&
                                  lon2d(i,j), lat2d(i,j))
    lon2d(i,j) = lon2d(i,j) * rad2deg
    lat2d(i,j) = lat2d(i,j) * rad2deg

    call phys2ij_Him8(imax_him8,jmax_him8,lon_him8,lat_him8,lon2d(i,j),lat2d(i,j),i_him8,j_him8)


    is = max(i_him8 - dix,1)
    ie = min(i_him8 + dix, imax_him8)

    js = max(j_him8 - diy,1)
    je = min(j_him8 + diy, jmax_him8)

    cnt = 0
    do jj = js, je
    do ii = is, ie
      if (minval(tbb_org(ii,jj,:)) < 0.0) cycle ! undef
      cnt = cnt + 1

      do k = 1, NIRB_HIM8
        tbb_sobs(i,j,k) = tbb_sobs(i,j,k) + real(tbb_org(ii,jj,k),kind=r_size)
      enddo
    enddo
    enddo
   
    if (cnt > 0) then
      do k = 1, NIRB_HIM8
        tbb_sobs(i,j,k) = tbb_sobs(i,j,k) / real(cnt,kind=r_size)
      enddo
    else
      tbb_sobs(i,j,:) = -1.0d0
    endif
 
  enddo ! i
  enddo ! j

  return
end subroutine sobs_Him8

subroutine phys2ij_Him8(imax_him8,jmax_him8,lon_him8,lat_him8,rlon,rlat,ig,jg)
  implicit none

  integer, intent(in) :: imax_him8, jmax_him8
  real(r_sngl),intent(in) :: lon_him8(imax_him8)
  real(r_sngl),intent(in) :: lat_him8(jmax_him8)
  real(r_size),intent(in) :: rlon
  real(r_size),intent(in) :: rlat
  integer,intent(out) :: ig
  integer,intent(out) :: jg

  real(r_sngl) :: dlon_him8, dlat_him8

!
! rlon,rlat -> ri,rj in Himawari 8 
!

  dlon_him8 = (maxval(lon_him8) - minval(lon_him8)) / real(imax_him8-1,kind=r_sngl)
  dlat_him8 = (maxval(lat_him8) - minval(lat_him8)) / real(jmax_him8-1,kind=r_sngl)

  ig = nint((rlon - minval(lon_him8)) / dlon_him8) + 1
  jg = nint((rlat - minval(lat_him8)) / dlat_him8) + 1

  if(ig > imax_him8 .or. ig < 1)then
    ig = -1
  endif
  if(jg > jmax_him8 .or. jg < 1)then
    jg = -1
  endif

  return
end subroutine phys2ij_Him8

subroutine get_nobs_allgHim8(nobs)
  implicit none

  integer, intent(out) :: nobs

  nobs = nlong * nlatg * NIRB_HIM8

  return
end subroutine get_nobs_allgHim8

subroutine allgHim82obs(tbb_allg,tbb_allg_prep,qc_allg_prep,obsdat,obslon,obslat,obslev,obserr)
  use scale_atmos_grid_cartesC, only: &
      CXG => ATMOS_GRID_CARTESC_CXG, &
      CYG => ATMOS_GRID_CARTESC_CYG, &
      DX, DY
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  implicit none

  real(r_size),intent(in) :: tbb_allg(nlong,nlatg,NIRB_HIM8)
  real(r_size),intent(out) :: tbb_allg_prep(nlong,nlatg,NIRB_HIM8)

  integer,intent(out),optional :: qc_allg_prep(nlong,nlatg,NIRB_HIM8)

  real(r_size),intent(out),optional :: obsdat(nlong*nlatg*NIRB_HIM8)
  real(r_size),intent(out),optional :: obslon(nlong*nlatg*NIRB_HIM8)
  real(r_size),intent(out),optional :: obslat(nlong*nlatg*NIRB_HIM8)
  real(r_size),intent(out),optional :: obslev(nlong*nlatg*NIRB_HIM8)
  real(r_size),intent(out),optional :: obserr(nlong*nlatg*NIRB_HIM8)

  real(r_size) :: ril, rjl
  real(r_size) :: lon, lat

  integer :: i, j
  integer :: ch
  integer :: n
  integer :: is, ie, js, je
  integer :: ii, jj
  integer :: ave_ng

  tbb_allg_prep(:,:,:) = undef

  if (present(obsdat) .and. present(obslon) .and. present(obslat) .and. present(obslev) .and. present(obserr)) then
    obsdat(:) = undef
    obslon(:) = 0.0
    obslat(:) = 0.0
    obslev(:) = undef
    obserr(:) = undef
  endif


  if (present(qc_allg_prep)) then
    qc_allg_prep = iqc_obs_bad
  endif

  ave_ng = 2 * H08_OBS_AVE_NG + 1

  do j = 1, nlatg
  do i = 1, nlong
    if (present(obslon) .and. present(obslat) .and. present(obslev) .and. present(obserr)) then
      ril = real(i+IHALO,kind=r_size)
      rjl = real(j+JHALO,kind=r_size)

      call MAPPROJECTION_xy2lonlat((ril - 1.0_r_size) * DX + CXG(1), &
                                   (rjl - 1.0_r_size) * DY + CYG(1), lon, lat)
    endif

    do ch = 1, NIRB_HIM8
      n = ((j - 1) * nlong + i - 1) * NIRB_HIM8 + ch

      if (present(obslon) .and. present(obslat) .and. present(obslev) .and. present(obserr)) then
        obslon(n) = lon * rad2deg
        obslat(n) = lat * rad2deg
        obslev(n) = ch + 6.0
        obserr(n) = REAL(OBSERR_H08(ch),r_size)
      endif

      select case(H08_OBS_METHOD)
      case(1) ! simple thinning
        tbb_allg_prep(i,j,ch) = tbb_allg(i,j,ch)

      case(2) ! averaging adjacent grids
        if (i <= H08_OBS_AVE_NG .or. (nlong - i) <= H08_OBS_AVE_NG .or.&
            j <= H08_OBS_AVE_NG .or. (nlatg - j) <= H08_OBS_AVE_NG) cycle

        is = i - H08_OBS_AVE_NG       
        ie = i + H08_OBS_AVE_NG       
        js = j - H08_OBS_AVE_NG       
        je = j + H08_OBS_AVE_NG       
   
        tbb_allg_prep(i,j,ch) = 0.0d0
        do jj = js, je
        do ii = is, ie
          tbb_allg_prep(i,j,ch) = tbb_allg_prep(i,j,ch) + tbb_allg(ii,jj,ch)
        enddo ! ii
        enddo ! jj
        tbb_allg_prep(i,j,ch) = tbb_allg_prep(i,j,ch) / (ave_ng**2)

      case(3) ! take a difference btw two bands
        tbb_allg_prep(i,j,ch) = tbb_allg(i,j,ch) -  tbb_allg(i,j,H08_OBS_SWD_B-6) 
      end select

      if (H08_OBS_THIN_LEV > 1) then
        if ((mod(i, H08_OBS_THIN_LEV) /= 0) .or. (mod(j, H08_OBS_THIN_LEV) /= 0)) then
          tbb_allg_prep(i,j,ch) = abs(tbb_allg_prep(i,j,ch)) * (-1.0d10) 
        endif
      endif


      if (present(obsdat)) then
        obsdat(n) = tbb_allg_prep(i,j,ch)
      endif

      if (present(qc_allg_prep)) then
        qc_allg_prep(i,j,ch) = iqc_good

        ! tbb_allg_prep can be negative when [H08_OBS_METHOD == 3]:
        ! take a difference btw two bands
        if (tbb_allg_prep(i,j,ch) < -200.0d0) then
          qc_allg_prep(i,j,ch) = iqc_obs_bad
        endif
      endif


    enddo ! ch
  enddo ! i
  enddo ! j

  return
end subroutine allgHim82obs

END MODULE common_obs_scale
