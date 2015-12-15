MODULE common_nml
!=======================================================================
!
! [PURPOSE:] Read namelist
!
! [HISTORY:]
!   November 2014   Guo-Yuan Lien     created
!
!=======================================================================
  use common, only: r_size
  use scale_stdio, only: IO_FID_CONF

  implicit none
  public

  !----
  integer, parameter :: nvarmax = 100
  integer, parameter :: nch = 10 ! Num of Himawari-8 (IR) channels 

  !--- PARAM_ENSEMBLE
  integer :: MEMBER = 3      ! ensemble size
  integer :: MEMBER_RUN = 1  !
  integer :: MEMBER_ITER = 0 !

  !--- PARAM_LETKF
  integer :: SLOT_START = 1
  integer :: SLOT_END = 1
  integer :: SLOT_BASE = 1
  real(r_size) :: SLOT_TINTERVAL = 3600.0d0

  real(r_size) :: SIGMA_OBS = 500.0d3
  real(r_size) :: SIGMA_OBS_RAIN = -1.0d0  ! < 0: same as SIGMA_OBS
  real(r_size) :: SIGMA_OBS_RADAR = -1.0d0 ! < 0: same as SIGMA_OBS
  real(r_size) :: SIGMA_OBS_H08 = -1.0d0 ! < 0: same as SIGMA_OBS ! H08
  real(r_size) :: SIGMA_OBSV = 0.4d0
  real(r_size) :: SIGMA_OBSV_RAIN = -1.0d0 ! < 0: same as SIGMA_OBSV
  real(r_size) :: SIGMA_OBSZ_RADAR = 1000.0d0
  real(r_size) :: SIGMA_OBSV_H08 = -1.0d0 ! < 0: same as SIGMA_OBSV ! H08
  real(r_size) :: SIGMA_OBST = 3.0d0
  real(r_size) :: BASE_OBSV_RAIN = 85000.0d0

  real(r_size) :: COV_INFL_MUL = 1.0d0    ! > 0: globally constant covariance inflation
                                          ! < 0: 3D inflation values input from a GPV file "infl_mul.grd"
  real(r_size) :: MIN_INFL_MUL = 0.0d0    ! minimum inlfation factor
  logical :: ADAPTIVE_INFL_INIT = .false.
  real(r_size) :: RELAX_ALPHA = 0.0d0        ! RTPP relaxation parameter
  real(r_size) :: RELAX_ALPHA_SPREAD = 0.0d0 ! RTPS relaxation parameter
  real(r_size) :: SP_INFL_ADD = 0.0d0        ! additive inflation

  real(r_size) :: GROSS_ERROR = 5.0d0
  real(r_size) :: GROSS_ERROR_RAIN = -1.0d0      ! < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_RADAR_REF = -1.0d0 ! < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_RADAR_VR = -1.0d0  ! < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_RADAR_PRH = -1.0d0 ! < 0: same as GROSS_ERROR

  integer :: LEV_UPDATE_Q = 100000        ! q and qc are only updated below and equal to this model level
  real(r_size) :: Q_SPRD_MAX = 0.5        ! maximum q (ensemble spread)/(ensemble mean)

  real(r_size) :: BOUNDARY_TAPER_WIDTH = 0.0d0

  !--- PARAM_LETKF_PRC
  integer :: NNODES = 1
  integer :: PPN = 1
  integer :: MEM_NODES = 1
  integer :: MEM_NP = 1
!  integer :: PRC_NUM_X_LETKF = 1
!  integer :: PRC_NUM_Y_LETKF = 1

  !--- PARAM_LETKF_VAR_LOCAL
  real(r_size) :: VAR_LOCAL_UV(nvarmax)        = 1.0d0
  real(r_size) :: VAR_LOCAL_T(nvarmax)         = 1.0d0
  real(r_size) :: VAR_LOCAL_Q(nvarmax)         = 1.0d0
  real(r_size) :: VAR_LOCAL_PS(nvarmax)        = 1.0d0
  real(r_size) :: VAR_LOCAL_RAIN(nvarmax)      = 1.0d0
  real(r_size) :: VAR_LOCAL_TC(nvarmax)        = 1.0d0
  real(r_size) :: VAR_LOCAL_RADAR_REF(nvarmax) = 1.0d0
  real(r_size) :: VAR_LOCAL_RADAR_VR(nvarmax)  = 1.0d0
  real(r_size) :: VAR_LOCAL_H08(nvarmax)       = 1.0d0 ! H08

  !--- PARAM_LETKF_OBSERR
  real(r_size) :: OBSERR_U = 1.0d0
  real(r_size) :: OBSERR_V = 1.0d0
  real(r_size) :: OBSERR_T = 1.0d0
  real(r_size) :: OBSERR_Q = 0.001d0
  real(r_size) :: OBSERR_RH = 10.0d0
  real(r_size) :: OBSERR_PS = 100.0d0
  real(r_size) :: OBSERR_RADAR_REF = 5.0d0
  real(r_size) :: OBSERR_RADAR_VR = 3.0d0
  real(r_size) :: OBSERR_H08 = 0.3d0 ! H08  
  logical :: USE_OBSERR_RADAR_REF = .false.
  logical :: USE_OBSERR_RADAR_VR = .false.

  !--- PARAM_LETKF_MONITOR
  logical :: DEPARTURE_STAT = .true.
  logical :: DEPARTURE_STAT_RADAR = .false.
  logical :: DEPARTURE_STAT_H08 = .false.
  real(r_size) :: DEPARTURE_STAT_T_RANGE = 0.0d0 ! time range within which observations are considered in the departure statistics.
                                                 ! 0.0d0: no limit

  LOGICAL :: OMB_OUTPUT = .true.
  LOGICAL :: OMA_OUTPUT = .true.
  LOGICAL :: OBSGUES_OUTPUT = .false.
  LOGICAL :: OBSANAL_OUTPUT = .false.

  !--- PARAM_LETKF_RADAR
  INTEGER :: MIN_RADAR_REF_MEMBER = 1          !Ensemble members with reflectivity greather than 0.
  REAL(r_size) :: MIN_RADAR_REF_DBZ = 0.0d0    !Minimum reflectivity
  REAL(r_size) :: RADAR_REF_THRES_DBZ = 15.0d0 !Threshold of rain/no rain

  REAL(r_size) :: RADAR_PRH_ERROR = 0.1d0      !Obserational error for pseudo RH observations.

  logical :: USE_RADAR_PSEUDO_RH = .false.

  real(r_size) :: RADAR_EDGE_TAPER_WIDTH = 0.0d0
  real(r_size) :: RADAR_RANGE = 0.0d0              !!!!!! should not use this and should save this information in obs files !!!!!!

  !These 2 flags affects the computation of model reflectivity and radial velocity. 
  INTEGER :: INTERPOLATION_TECHNIQUE = 1
  INTEGER :: METHOD_REF_CALC = 3

  LOGICAL :: USE_TERMINAL_VELOCITY = .false.

  ! PARAMETERS FOR RADAR DATA ASSIMILATION
  INTEGER :: NRADARTYPE = 1  !Currently PAWR (1) and LIDAR (2) ... not used?

  !---PARAM_LETKF_H08
  real(r_size) :: H08_LIMIT_LEV = 20000.0d0 ! (Pa) Upper limit level of the sensitive height for Himawari-8 IR
  integer :: H08_CH_USE(nch) = (/0,1,1,1,0,0,0,0,0,0/)
                        !! ch = (1,2,3,4,5,6,7,8,9,10)
                        !! (B07,B08,B09,B10,B11,B12,B13,B14,B15,B16)
                        !! ==1: Assimilate
                        !! ==0: NOT assimilate (rejected by QC in trans_XtoY_H08)
                        !! It is better to reject B11(ch=5) & B12(ch=6) obs because these bands are 
                        !! sensitive to chemicals.

contains
!-----------------------------------------------------------------------
! PARAM_ENSEMBLE
!-----------------------------------------------------------------------
subroutine read_nml_ensemble
  implicit none
  integer :: ierr
  
  namelist /PARAM_ENSEMBLE/ &
    MEMBER, &
    MEMBER_RUN, &
    MEMBER_ITER

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_ENSEMBLE,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Error: /PARAM_ENSEMBLE/ is not found in namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_ENSEMBLE. Check!'
    stop
  endif

  write(6, nml=PARAM_ENSEMBLE)

  return
end subroutine read_nml_ensemble

!-----------------------------------------------------------------------
! PARAM_LETKF
!-----------------------------------------------------------------------
subroutine read_nml_letkf
  implicit none
  integer :: ierr
  
  namelist /PARAM_LETKF/ &
    SLOT_START, &
    SLOT_END, &
    SLOT_BASE, &
    SLOT_TINTERVAL, &
    SIGMA_OBS, &
    SIGMA_OBS_RAIN, &
    SIGMA_OBS_RADAR, &
    SIGMA_OBS_H08, & ! H08
    SIGMA_OBSV, &
    SIGMA_OBSV_RAIN, &
    SIGMA_OBSV_H08, & ! H08
    SIGMA_OBSZ_RADAR, &
    SIGMA_OBST, &
    BASE_OBSV_RAIN, &
    COV_INFL_MUL, &
    MIN_INFL_MUL, &
    ADAPTIVE_INFL_INIT, &
    RELAX_ALPHA, &
    RELAX_ALPHA_SPREAD, &
    SP_INFL_ADD, &
    GROSS_ERROR, &
    GROSS_ERROR_RAIN, &
    GROSS_ERROR_RADAR_REF, &
    GROSS_ERROR_RADAR_VR, &
    GROSS_ERROR_RADAR_PRH, &
    LEV_UPDATE_Q, &
    Q_SPRD_MAX, &
    BOUNDARY_TAPER_WIDTH

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Error: /PARAM_LETKF/ is not found in namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF. Check!'
    stop
  endif

  if (GROSS_ERROR_RAIN < 0.0d0) then
    GROSS_ERROR_RAIN = GROSS_ERROR
  end if
  if (GROSS_ERROR_RADAR_REF < 0.0d0) then
    GROSS_ERROR_RADAR_REF = GROSS_ERROR
  end if
  if (GROSS_ERROR_RADAR_VR < 0.0d0) then
    GROSS_ERROR_RADAR_VR = GROSS_ERROR
  end if
  if (GROSS_ERROR_RADAR_PRH < 0.0d0) then
    GROSS_ERROR_RADAR_PRH = GROSS_ERROR
  end if
  if (SIGMA_OBS_RAIN < 0.0d0) then
    SIGMA_OBS_RAIN = SIGMA_OBS
  end if
  if (SIGMA_OBS_RADAR < 0.0d0) then
    SIGMA_OBS_RADAR = SIGMA_OBS
  end if
  if (SIGMA_OBSV_RAIN < 0.0d0) then
    SIGMA_OBSV_RAIN = SIGMA_OBSV
  end if

  if (SIGMA_OBS_H08 < 0.0d0) then ! H08
    SIGMA_OBS_H08 = SIGMA_OBS
  end if
  if (SIGMA_OBSV_H08 < 0.0d0) then ! H08
    SIGMA_OBSV_H08 = SIGMA_OBSV
  end if

  write(6, nml=PARAM_LETKF)

  return
end subroutine read_nml_letkf

!-----------------------------------------------------------------------
! PARAM_LETKF_PRC
!-----------------------------------------------------------------------
subroutine read_nml_letkf_prc
  implicit none
  integer :: ierr
  
  namelist /PARAM_LETKF_PRC/ &
    NNODES, &
    PPN, &
    MEM_NODES, &
    MEM_NP
!    PRC_NUM_X_LETKF, &
!    PRC_NUM_Y_LETKF

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_PRC,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_LETKF_PRC/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_PRC. Check!'
    stop
  endif

  write(6, nml=PARAM_LETKF_PRC)

  return
end subroutine read_nml_letkf_prc

!-----------------------------------------------------------------------
! PARAM_LETKF_VAR_LOCAL
!-----------------------------------------------------------------------
subroutine read_nml_letkf_var_local
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_VAR_LOCAL/ &
    VAR_LOCAL_UV, &
    VAR_LOCAL_T, &
    VAR_LOCAL_Q, &
    VAR_LOCAL_PS, &
    VAR_LOCAL_RAIN, &
    VAR_LOCAL_TC, &
    VAR_LOCAL_RADAR_REF, &
    VAR_LOCAL_RADAR_VR, &
    VAR_LOCAL_H08 ! H08

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_VAR_LOCAL,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_LETKF_VAR_LOCAL/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_VAR_LOCAL. Check!'
    stop
  endif

  write(6, nml=PARAM_LETKF_VAR_LOCAL)

  return
end subroutine read_nml_letkf_var_local

!-----------------------------------------------------------------------
! PARAM_LETKF_OBSERR
!-----------------------------------------------------------------------
subroutine read_nml_letkf_obserr
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_OBSERR/ &
    OBSERR_U, &
    OBSERR_V, &
    OBSERR_T, &
    OBSERR_Q, &
    OBSERR_RH, &
    OBSERR_PS, &
    OBSERR_RADAR_REF, &
    OBSERR_RADAR_VR, &
    OBSERR_H08, & ! H08
    USE_OBSERR_RADAR_REF, &
    USE_OBSERR_RADAR_VR

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_OBSERR,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_LETKF_OBSERR/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_OBSERR. Check!'
    stop
  endif

  write(6, nml=PARAM_LETKF_OBSERR)

  return
end subroutine read_nml_letkf_obserr

!-----------------------------------------------------------------------
! PARAM_LETKF_MONITOR
!-----------------------------------------------------------------------
subroutine read_nml_letkf_monitor
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_MONITOR/ &
    DEPARTURE_STAT, &
    DEPARTURE_STAT_RADAR, &
    DEPARTURE_STAT_H08, &
    DEPARTURE_STAT_T_RANGE, &
    OMB_OUTPUT, &
    OMA_OUTPUT, &
    OBSGUES_OUTPUT, &
    OBSANAL_OUTPUT

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_MONITOR,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_LETKF_MONITOR/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_MONITOR. Check!'
    stop
  endif

  write(6, nml=PARAM_LETKF_MONITOR)

  return
end subroutine read_nml_letkf_monitor

!-----------------------------------------------------------------------
! PARAM_LETKF_RADAR
!-----------------------------------------------------------------------
subroutine read_nml_letkf_radar
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_RADAR/ &
    MIN_RADAR_REF_MEMBER, &
    MIN_RADAR_REF_DBZ, &
    RADAR_REF_THRES_DBZ, &
    RADAR_PRH_ERROR, &
    USE_RADAR_PSEUDO_RH, &
    RADAR_EDGE_TAPER_WIDTH, &
    RADAR_RANGE, &               !!!!!! should not use this and should save this information in obs files !!!!!!
    INTERPOLATION_TECHNIQUE, &
    METHOD_REF_CALC, &
    USE_TERMINAL_VELOCITY, &
    NRADARTYPE

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_RADAR,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_LETKF_RADAR/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_RADAR. Check!'
    stop
  endif

  if (RADAR_REF_THRES_DBZ < MIN_RADAR_REF_DBZ) then
    RADAR_REF_THRES_DBZ = MIN_RADAR_REF_DBZ
  end if

  write(6, nml=PARAM_LETKF_RADAR)

  return
end subroutine read_nml_letkf_radar

!-----------------------------------------------------------------------
! PARAM_LETKF_RADAR
!-----------------------------------------------------------------------
subroutine read_nml_letkf_h08
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_H08/ &
    H08_LIMIT_LEV, &
    H08_CH_USE

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_H08,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_LETKF_H08/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_H08. Check!'
    stop
  endif

  write(6, nml=PARAM_LETKF_H08)

  return
end subroutine read_nml_letkf_h08

end module common_nml
