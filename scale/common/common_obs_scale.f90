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

  INTEGER,PARAMETER :: nid_obs_varlocal=9 !H08
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
       'TMPAPR', 'PHARAD', 'H08IRB', 'TCVITL'/) ! H08

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
    REAL(r_size),ALLOCATABLE :: lev(:) ! H08
    REAL(r_size),ALLOCATABLE :: val2(:) ! H08 ! clear sky BT (gues)
#endif
    REAL(r_size),ALLOCATABLE :: ensval(:,:)
    INTEGER,ALLOCATABLE :: qc(:)
  END TYPE obs_da_value

  character(obsformatlenmax), parameter :: obsfmt_prepbufr     = 'PREPBUFR'
  character(obsformatlenmax), parameter :: obsfmt_radar        = 'RADAR'
  character(obsformatlenmax), parameter :: obsfmt_pawr_toshiba = 'PAWR_TOSHIBA'
  character(obsformatlenmax), parameter :: obsfmt_pawr_jrc     = 'PAWR_JRC'
  character(obsformatlenmax), parameter :: obsfmt_h08          = 'HIMAWARI8'
!  integer, parameter :: nobsformats = 4
!  character(obsformatlenmax), parameter :: obsformat(nobsformats) = &
!    (/obsfmt_prepbufr, obsfmt_radar, obsfmt_pawr_toshiba, obsfmt_h08/)

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
  real(RP) :: lon_tmp(1,1),lat_tmp(1,1)

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
      if (LOG_LEVEL >= 4) then
        write(6,'(A,F8.1,A,F8.1,A,I5)') '[Warning] observation is too high: ptop=', ptmp, ', lev=', rlev, ', elem=', elem
      end if
      rk = undef
      qc = iqc_out_vhi
      RETURN
    END IF
    IF(rk > plev(ks)) THEN
      call itpl_2d(p_full(ks,:,:),ri,rj,ptmp)
!print *, ks, rk, plev(ks)
      if (LOG_LEVEL >= 4) then
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
    if (LOG_LEVEL >= 4) then
      write(6,'(A,F8.1,A,F8.1)') '[Warning] observation is too high: ztop=', ztmp, ', lev=', rlev
    end if
    rk = undef
    qc = iqc_out_vhi
    RETURN
  END IF
  IF(rlev < zlev(ks)) THEN
    call itpl_2d(z_full(ks,:,:),ri,rj,ztmp)
    if (LOG_LEVEL >= 4) then
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
      ATMOS_GRID_CARTESC_CXG, &
      ATMOS_GRID_CARTESC_CYG, &
      DX, &
      DY
  use scale_mapprojection, only: &
      MAPPROJECTION_lonlat2xy
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(OUT) :: rig
  REAL(r_size),INTENT(OUT) :: rjg
  real(RP) :: rig_RP, rjg_RP
!
! rlon,rlat -> ri,rj
!
  call MAPPROJECTION_lonlat2xy(real(rlon*pi/180.0_r_size,RP),real(rlat*pi/180.0_r_size,RP),rig_RP,rjg_RP)
  rig = (real(rig_RP,r_size) - ATMOS_GRID_CARTESC_CXG(1)) / DX + 1.0d0
  rjg = (real(rjg_RP,r_size) - ATMOS_GRID_CARTESC_CYG(1)) / DY + 1.0d0

  RETURN
END SUBROUTINE phys2ij

SUBROUTINE ij2phys(rig,rjg,rlon,rlat)
  use scale_atmos_grid_cartesC, only: &
      ATMOS_GRID_CARTESC_CXG, &
      ATMOS_GRID_CARTESC_CYG, &
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
  real(RP) :: rlon_RP, rlat_RP
!
! ri,rj -> rlon,rlat
!
  x = (rig - 1.0d0) * DX + ATMOS_GRID_CARTESC_CXG(1) 
  y = (rjg - 1.0d0) * DY + ATMOS_GRID_CARTESC_CYG(1) 

  call MAPPROJECTION_xy2lonlat(real(x,RP),real(y,RP),rlon_RP,rlat_RP)

  rlon = real(rlon_RP,r_size) * rad2deg
  rlat = real(rlat_RP,r_size) * rad2deg

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
!
! monit_obs is ported into obs/obs_tools.f90
!
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
  IF(ALLOCATED(obs%ri)) DEALLOCATE(obs%ri)
  IF(ALLOCATED(obs%rj)) DEALLOCATE(obs%rj)
  IF(ALLOCATED(obs%rank)) DEALLOCATE(obs%rank)

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
  ALLOCATE( obsda%lev (obsda%nobs) ) ! H08
  ALLOCATE( obsda%val2(obsda%nobs) ) ! H08
#endif
  ALLOCATE( obsda%qc  (obsda%nobs) )

  obsda%nobs_in_key = 0
  obsda%idx = 0
  obsda%key = 0
  obsda%val = 0.0d0
#ifdef H08
  obsda%lev = 0.0d0 ! H08
  obsda%val2 = 0.0d0 ! H08
#endif
  obsda%qc = 0

  if (member > 0) then
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
  IF(ALLOCATED(obsda%lev   )) DEALLOCATE(obsda%lev   ) ! H08
  IF(ALLOCATED(obsda%val2  )) DEALLOCATE(obsda%val2  ) ! H08
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
  REAL(RP) :: x, y
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
      call MAPPROJECTION_lonlat2xy(real(wk(2)*pi/180.0_r_size,RP),&
                                   real(wk(3)*pi/180.0_r_size,RP),&
                                   x,y)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = real(x,kind=r_sngl)
      wk(6) = real(OBSERR_TCX,kind=r_sngl)
    CASE(id_tclat_obs)
      call MAPPROJECTION_lonlat2xy(real(wk(2)*pi/180.0_r_size,RP),&
                                   real(wk(3)*pi/180.0_r_size,RP),&
                                   x,y)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = real(y,kind=r_sngl)
      wk(6) = real(OBSERR_TCY,kind=r_sngl)
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
  REAL(r_sngl) :: wk(6) ! H08
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
    obsda%lev(n) = REAL(wk(5),r_size) ! H08
    obsda%val2(n) = REAL(wk(6),r_size) ! H08
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
    wk(5) = REAL(obsda%lev(n),r_sngl) ! H08
    wk(6) = REAL(obsda%val2(n),r_sngl) ! H08
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
!
!-----------------------------------------------------------------------
!   TC vital obs subroutines by T. Honda (03/28/2016)
!-----------------------------------------------------------------------
!
SUBROUTINE search_tc_subdom(ritc,rjtc,v2d,yobs_tcx,yobs_tcy,yobs_mslp)
  use scale_atmos_grid_cartesC, only: &
      ATMOS_GRID_CARTESC_CXG, &
      ATMOS_GRID_CARTESC_CYG, &
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
      yobs_tcx = (real(ig,kind=r_size) - 1.0d0) * DX + ATMOS_GRID_CARTESC_CXG(1)
      yobs_tcy = (real(jg,kind=r_size) - 1.0d0) * DY + ATMOS_GRID_CARTESC_CYG(1)
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
#ifdef H08
!
SUBROUTINE Trans_XtoY_H08(nprof,ri,rj,lon,lat,v3d,v2d,yobs,plev_obs,qc,stggrd,yobs_H08_clr)
  use scale_mapprojection, only: &
      MAPPROJECTION_rotcoef
  use scale_H08_fwd
  use scale_atmos_grid_cartesC_index, only: &
    KHALO

  IMPLICIT NONE
  INTEGER :: n, np, k, ch
  INTEGER,INTENT(IN) :: nprof ! Num of Brightness Temp "Loc" observed by Himawari-8
                              ! NOTE: multiple channels (obs) on each grid point !!
  REAL(r_size),INTENT(IN) :: ri(nprof),rj(nprof)
  REAL(r_size),INTENT(IN) :: lon(nprof),lat(nprof)
  REAL(r_size),INTENT(IN) :: v3d(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  INTEGER,INTENT(IN),OPTIONAL :: stggrd
  REAL(RP) :: rotc(2)

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

! -- brightness temp from RTTOV
  REAL(r_size) :: btall_out(nch,nprof) ! NOTE: RTTOV always calculates all (10) channels!!
  REAL(r_size) :: btclr_out(nch,nprof) ! NOTE: RTTOV always calculates all (10) channels!!
! -- transmittance from RTTOV
  REAL(r_size) :: trans_out(nlev,nch,nprof)
 
  REAL(r_size) :: max_weight(nch,nprof)
  REAL(r_size) :: tmp_weight

  REAL(r_size),INTENT(OUT) :: yobs(nprof*nch)
  REAL(r_size),OPTIONAL,INTENT(OUT) :: yobs_H08_clr(nprof*nch)
  INTEGER,INTENT(OUT) :: qc(nprof*nch)
  REAL(r_size),INTENT(OUT) :: plev_obs(nch*nprof)

  REAL(r_size) :: rdp ! delta p
  INTEGER :: slev, elev

  REAL(r_size) :: utmp, vtmp ! U10m & V10m tmp for rotation
  REAL(RP) :: rotc(1,1,2)
  real(r_size) :: lon_tmp(1,1),lat_tmp(1,1)

  if (present(stggrd)) stggrd_ = stggrd

  yobs = undef
  qc = iqc_good

  lon1d(:) = lon(:)
  lat1d(:) = lat(:)

! -- make profile arrays for RTTOV --
  DO np = 1, nprof ! -- make profiles

    CALL itpl_2d(v2d(:,:,iv2dd_skint),ri(np),rj(np),tsfc1d(np)) ! T2 is better??
!    CALL itpl_2d(v2d(:,:,iv2dd_t2m),ri(np),rj(np),tsfc1d(np))
    CALL itpl_2d(v2d(:,:,iv2dd_q2m),ri(np),rj(np),qsfc1d(np))
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

  slev = 1 + KHALO
  elev = nlevh - KHALO

  CALL SCALE_RTTOV_fwd(nch, & ! num of channels
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
                       btall_out(1:nch,1:nprof),& ! (K)
                       btclr_out(1:nch,1:nprof),& ! (K)
                       trans_out(nlev:1:-1,1:nch,1:nprof))
!
! -- Compute max weight level using trans_out 
! -- (Transmittance from each user pressure level to Top Of the Atmosphere)
! -- btall_out is substituted into yobs

  n = 0
  DO np = 1, nprof
  DO ch = 1, nch
    n = n + 1

    rdp = 1.0d0 / (prs2d(slev+1,np) - prs2d(slev,np))
    max_weight(ch,np) = abs((trans_out(2,ch,np) - trans_out(1,ch,np)) * rdp )


    plev_obs(n) = (prs2d(slev+1,np) + prs2d(slev,np)) * 0.5d0 ! (Pa)

    DO k = 2, nlev-1

      rdp = 1.0d0 / abs(prs2d(slev+k,np) - prs2d(slev+k-1,np))
      tmp_weight = (trans_out(k+1,ch,np) - trans_out(k,ch,np)) * rdp 
      if(tmp_weight > max_weight(ch,np))then
        max_weight(ch,np) = tmp_weight
        plev_obs(n) = (prs2d(slev+k,np) + prs2d(slev+k-1,np)) * 0.5d0 ! (Pa)
      endif
    ENDDO

    yobs(n) = btall_out(ch,np)
!
! ## comment out by T.Honda (02/09/2016)
! -- tentative QC here --
!    IF(plev_obs(n) >= H08_LIMIT_LEV)THEN
!      qc(n) = iqc_good
!    ELSE
!      qc(n) = iqc_obs_bad
!    ENDIF

    SELECT CASE(H08_CH_USE(ch))
    CASE(1)
      qc(n) = iqc_good
    CASE DEFAULT
      qc(n) = iqc_obs_bad
    END SELECT

    IF(H08_REJECT_LAND .and. (lsmask1d(np) > 0.5d0))THEN
      qc(n) = iqc_obs_bad
    ENDIF

    IF(abs(btall_out(ch,np) - btclr_out(ch,np)) > H08_CLDSKY_THRS)THEN
! Cloudy sky
      yobs(n) = yobs(n) * (-1.0d0)
    ELSE
! Clear sky
      yobs(n) = yobs(n) * 1.0d0
    ENDIF
   
    yobs_H08_clr(n) = btclr_out(ch,np)

  ENDDO ! ch
  ENDDO ! np

  RETURN
END SUBROUTINE Trans_XtoY_H08
#endif

SUBROUTINE get_nobs_H08(cfile,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn ! num of all H08 obs
!  REAL(r_sngl) :: wk(4+nch)
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
!      nn = nn + nch
!    END DO
!-----------------------------

! get file size by INQUIRE statement... may not work for some older fortran compilers
!-----------------------------
    INQUIRE(UNIT=iunit, SIZE=sz)
    IF (MOD(sz, r_sngl * (4+nch+2)) /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    iprof = sz / (r_sngl * (4+nch+2))
    nn = iprof * nch
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
  REAL(r_sngl) :: wk(4+nch)
!  REAL(r_sngl) :: tmp
  INTEGER :: n,iunit

  INTEGER :: nprof, np, ch

  nprof = obs%nobs / nch

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

  n = 0
  DO np=1,nprof
    READ(iunit) wk

    DO ch = 1, nch
      n = n + 1

      obs%elm(n) = NINT(wk(1))
      obs%typ(n) = NINT(wk(2))
      obs%lon(n) = REAL(wk(3),r_size)
      obs%lat(n) = REAL(wk(4),r_size)
      obs%dat(n) = REAL(wk(4+ch),r_size)
      obs%dif(n) = 0.0d0
      obs%lev(n) = ch + 6.0 ! substitute channnel number instead of the obs level
      obs%err(n) = REAL(OBSERR_H08(ch),r_size)
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
  REAL(r_sngl) :: wk(4+nch)
  INTEGER :: n,iunit
  INTEGER :: iprof, ns, ne

  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  missing_ = .true.
  IF(present(missing)) missing_ = missing

  iprof = obs%nobs / nch


  IF(append_) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  END IF

  DO n=1,iprof
    ns = (n-1)*nch + 1
    ne = n*nch

    wk(1) = REAL(obs%elm(ns),r_sngl)
    wk(2) = REAL(obs%typ(ns),r_sngl)
    wk(3) = REAL(obs%lon(ns),r_sngl)
    wk(4) = REAL(obs%lat(ns),r_sngl)
    wk(5:5+nch-1) = REAL(obs%dat(ns:ne),r_size)
    WRITE(iunit) wk

  ENDDO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_H08

END MODULE common_obs_scale
