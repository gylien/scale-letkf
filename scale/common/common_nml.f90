MODULE common_nml
!===============================================================================
!
! [PURPOSE:] Read namelist
!
! [HISTORY:]
!   November 2014   Guo-Yuan Lien     created
!   .............   See git history for the following revisions
!
!===============================================================================
  use common, only: r_size
  use scale_stdio, only: IO_FID_CONF

  implicit none
  public

  !----
  integer, parameter :: nv3d = 11    ! number of 3D state variables (in SCALE restart files)
  integer, parameter :: nv2d = 0     ! number of 2D state variables (in SCALE restart files)
  integer, parameter :: nid_obs = 16 ! number of variable types
  integer, parameter :: nobtype = 24 ! number of observation report types
  integer, parameter :: nch = 10     ! H08 Num of Himawari-8 (IR) channels

  integer, parameter :: nobsfilemax = 10
  integer, parameter :: obsformatlenmax = 10
  integer, parameter :: filelenmax = 256

  integer, parameter :: memflen = 4                           ! Length of formatted member strings
  character(len=memflen), parameter :: memf_notation = '@@@@' ! Notation of the member string
  character(len=memflen), parameter :: memf_mean = 'mean'
  character(len=memflen), parameter :: memf_mdet = 'mdet'
  character(len=memflen), parameter :: memf_sprd = 'sprd'

  !--- PARAM_ENSEMBLE
  integer :: MEMBER = 3      ! ensemble size
  integer :: MEMBER_RUN = 1  !
  integer :: MEMBER_ITER = 0 !
  character(filelenmax) :: CONF_FILES = 'run.@@@@.conf'
  logical :: CONF_FILES_SEQNUM = .false.

  logical :: DET_RUN = .false.
  logical :: DET_RUN_CYCLED = .true.

  !--- PARAM_MODEL
  character(len=10) :: MODEL = 'scale-rm'
  logical :: VERIFY_COORD = .false.

!  !--- PARAM_IO
!  integer :: IO_AGGREGATE = .false.

  !--- PARAM_LOG
  integer :: LOG_LEVEL = 2                        ! Log message output level:
                                                  !  0: Minimum log output
                                                  !  1: Reduced log output
                                                  !  2: Normal  log output
                                                  !  3: Verbose log output
  logical :: USE_MPI_BARRIER = .true.             ! Whether enabling some MPI_Barrier for better timing measurement?

  !--- PARAM_OBSOPE
  integer               :: OBS_IN_NUM = 1
  character(filelenmax) :: OBS_IN_NAME(nobsfilemax) = 'obs.dat'
  character(obsformatlenmax) :: OBS_IN_FORMAT(nobsfilemax) = 'PREPBUFR'
  logical               :: OBSDA_RUN(nobsfilemax) = .true.
  logical               :: OBSDA_OUT = .false.
  character(filelenmax) :: OBSDA_OUT_BASENAME = 'obsda.@@@@'
  character(filelenmax) :: OBSDA_MEAN_OUT_BASENAME = ''
  character(filelenmax) :: OBSDA_MDET_OUT_BASENAME = ''

  character(filelenmax) :: HISTORY_IN_BASENAME = 'hist.@@@@'
  character(filelenmax) :: HISTORY_MEAN_IN_BASENAME = ''
  character(filelenmax) :: HISTORY_MDET_IN_BASENAME = ''

  integer               :: SLOT_START = 1
  integer               :: SLOT_END = 1
  integer               :: SLOT_BASE = 1
  real(r_size)          :: SLOT_TINTERVAL = 3600.0d0

  !--- PARAM_LETKF
  logical               :: OBSDA_IN = .false.
  character(filelenmax) :: OBSDA_IN_BASENAME = 'obsda.@@@@'
  character(filelenmax) :: OBSDA_MEAN_IN_BASENAME = ''
  character(filelenmax) :: OBSDA_MDET_IN_BASENAME = ''
  character(filelenmax) :: GUES_IN_BASENAME = 'gues.@@@@'
  character(filelenmax) :: GUES_MEAN_INOUT_BASENAME = ''
  character(filelenmax) :: GUES_MDET_IN_BASENAME = ''
  logical               :: GUES_SPRD_OUT = .true.
  character(filelenmax) :: GUES_SPRD_OUT_BASENAME = ''
  character(filelenmax) :: ANAL_OUT_BASENAME = 'anal.@@@@'
  character(filelenmax) :: ANAL_MEAN_OUT_BASENAME = ''
  character(filelenmax) :: ANAL_MDET_OUT_BASENAME = ''
  logical               :: ANAL_SPRD_OUT = .true.
  character(filelenmax) :: ANAL_SPRD_OUT_BASENAME = ''
  character(filelenmax) :: LETKF_TOPO_IN_BASENAME = 'topo'  !!!!!! -- directly use the SCALE namelist --???? !!!!!!

  real(r_size) :: INFL_MUL = 1.0d0           ! >  0: globally constant covariance inflation
                                             ! <= 0: use 3D inflation field from 'INFL_MUL_IN_BASENAME' file
  real(r_size) :: INFL_MUL_MIN = -1.0d0      ! minimum inlfation factor (<= 0: not used)
  logical :: INFL_MUL_ADAPTIVE = .false.     ! if true, outout adaptively estimated 3D inlfation field to 'INFL_MUL_OUT_BASENAME' file
  character(filelenmax) :: INFL_MUL_IN_BASENAME = 'infl'
  character(filelenmax) :: INFL_MUL_OUT_BASENAME = 'infl'

  real(r_size) :: INFL_ADD = 0.0d0           ! additive inflation
  character(filelenmax) :: INFL_ADD_IN_BASENAME = 'addi.@@@@'
  logical :: INFL_ADD_SHUFFLE = .false.      ! shuffle the additive inflation members?
  logical :: INFL_ADD_Q_RATIO = .false.
  logical :: INFL_ADD_REF_ONLY = .false.

  real(r_size) :: RELAX_ALPHA = 0.0d0        ! RTPP relaxation parameter
  real(r_size) :: RELAX_ALPHA_SPREAD = 0.0d0 ! RTPS relaxation parameter
  logical :: RELAX_TO_INFLATED_PRIOR = .false. ! .true. : relaxation to multiplicatively inflated prior
                                               ! .false.: relaxation to original prior
  logical :: RELAX_SPREAD_OUT = .false.
  character(filelenmax) :: RELAX_SPREAD_OUT_BASENAME = 'rtps'

  real(r_size) :: GROSS_ERROR = 5.0d0
  real(r_size) :: GROSS_ERROR_RAIN = -1.0d0      ! < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_RADAR_REF = -1.0d0 ! < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_RADAR_VR = -1.0d0  ! < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_RADAR_PRH = -1.0d0 ! < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_H08 = -1.0d0      ! < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_TCX = -1.0d0 ! debug ! < 0: same as GROSS_ERROR 
  real(r_size) :: GROSS_ERROR_TCY = -1.0d0 ! debug ! < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_TCP = -1.0d0 ! debug ! < 0: same as GROSS_ERROR

  real(r_size) :: Q_UPDATE_TOP = 0.0d0     ! water vapor and hydrometeors are updated only below this pressure level (Pa)
  real(r_size) :: Q_SPRD_MAX = -1.0D0      ! maximum q (ensemble spread)/(ensemble mean) (only effective when > 0)

  real(r_size) :: BOUNDARY_BUFFER_WIDTH = 0.0d0

  logical :: POSITIVE_DEFINITE_Q = .false.
  logical :: POSITIVE_DEFINITE_QHYD = .false.
  real(r_size) :: TC_SEARCH_DIS = 200.0d3 ! (m) ! tentative! Should be modify !!

  real(r_size) :: PS_ADJUST_THRES = 100.d0

  logical :: NOBS_OUT = .false.
  character(filelenmax) :: NOBS_OUT_BASENAME = 'nobs'

  !*** for backward compatibility ***
  real(r_size) :: COV_INFL_MUL = 1.0d0
  real(r_size) :: MIN_INFL_MUL = 0.0d0
  logical :: ADAPTIVE_INFL_INIT = .false.
  real(r_size) :: BOUNDARY_TAPER_WIDTH = 0.0d0

  !--- PARAM_LETKF_PRC
  integer :: NNODES = 1
  integer :: PPN = 1
  integer :: MEM_NODES = 1
  integer :: MEM_NP = 1
!  integer :: PRC_NUM_X_LETKF = 1
!  integer :: PRC_NUM_Y_LETKF = 1

  !--- PARAM_LETKF_OBS
  logical :: USE_OBS(nobtype) = .true.

  ! >0: localization length scale (m)
  !  0: no localization XXX not implemented yet XXX
  ! <0: same as HORI_LOCAL(1)
  real(r_size) :: HORI_LOCAL(nobtype) = &
    (/500.0d3, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
       -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
       -1.0d0, -1.0d0, -1.0d0, -1.0d0/)

  ! >0: localization length scale [ln(p) or m depends on obstype]
  !  0: no localization
  ! <0: same as VERT_LOCAL(1)
  real(r_size) :: VERT_LOCAL(nobtype) = &
    (/ 0.4d0,   -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0,   -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, 1000.0d0, -1.0d0, -1.0d0/)
!      -1.0d0, 1000.0d0, -1.0d0,  0.0d0/)

  ! >0: localization length scale (sec) XXX not implemented yet XXX
  !  0: no localization
  ! <0: same as TIME_LOCAL(1)
  real(r_size) :: TIME_LOCAL(nobtype) = &
    (/ 0.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0/)

  real(r_size) :: HORI_LOCAL_RADAR_OBSNOREF = -1.0d0 ! <0: same as HORI_LOCAL(22=PHARAD)
  real(r_size) :: HORI_LOCAL_RADAR_VR = -1.0d0       ! <0: same as HORI_LOCAL(22=PHARAD)
  real(r_size) :: VERT_LOCAL_RADAR_VR = -1.0d0       ! <0: same as VERT_LOCAL(22=PHARAD)
  real(r_size) :: VERT_LOCAL_RAIN_BASE = 85000.0d0

  ! >0: observation number limit
  !  0: do not limit observation numbers
  ! <0: same as MAX_NOBS_PER_GRID(1)
  integer :: MAX_NOBS_PER_GRID(nobtype) = &
    (/ 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
      -1, -1, -1, -1/)

  integer :: MAX_NOBS_PER_GRID_CRITERION = 1 ! 1: normalized 3D distance (from closest)
                                             ! 2: localization weight (from largest)
                                             ! 3: weighted observation error variance (from smallest)

  ! >0: typical minimum spacing of the obsetvation types in the densest observed area (not tuned carefully yet)
  !     *this is only used for automatically determine OBS_SORT_GRID_SPACING. if using pre-set OBS_SORT_GRID_SPACING, this has no effect.
  ! <=0: same as OBS_MIN_SPACING(1)
  real(r_size) :: OBS_MIN_SPACING(nobtype) = &
    (/300.0d3, 100.0d3, 100.0d3, 150.0d3, 300.0d3, 150.0d3, 150.0d3, 100.0d3, 150.0d3, 150.0d3, &
      150.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 300.0d3, 150.0d3, 150.0d3, 150.0d3, &
      150.0d3,   1.0d3,  15.0d3,1000.0d3/)

  ! >0: optimal grid spacing for bucket sorting of observations
  !  0: automatically determined based on HORI_LOCAL, MAX_NOBS_PER_GRID, and OBS_MIN_SPACING
  ! <0: same as OBS_SORT_GRID_SPACING(1)
  real(r_size) :: OBS_SORT_GRID_SPACING(nobtype) = &
    (/ 0.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0/)

  !--- PARAM_LETKF_VAR_LOCAL
  real(r_size) :: VAR_LOCAL_UV(nv3d+nv2d)        = 1.0d0
  real(r_size) :: VAR_LOCAL_T(nv3d+nv2d)         = 1.0d0
  real(r_size) :: VAR_LOCAL_Q(nv3d+nv2d)         = 1.0d0
  real(r_size) :: VAR_LOCAL_PS(nv3d+nv2d)        = 1.0d0
  real(r_size) :: VAR_LOCAL_RAIN(nv3d+nv2d)      = 1.0d0
  real(r_size) :: VAR_LOCAL_TC(nv3d+nv2d)        = 1.0d0
  real(r_size) :: VAR_LOCAL_RADAR_REF(nv3d+nv2d) = 1.0d0
  real(r_size) :: VAR_LOCAL_RADAR_VR(nv3d+nv2d)  = 1.0d0
  real(r_size) :: VAR_LOCAL_H08(nv3d+nv2d)       = 1.0d0 ! H08

  !--- PARAM_LETKF_MONITOR
  logical :: DEPARTURE_STAT = .true.
  logical :: DEPARTURE_STAT_RADAR = .false.
  logical :: DEPARTURE_STAT_H08 = .false.
  real(r_size) :: DEPARTURE_STAT_T_RANGE = 0.0d0   ! time range within which observations are considered in the departure statistics.
                                                   ! 0: no limit
  logical :: DEPARTURE_STAT_ALL_PROCESSES = .true. ! print the departure statistics by all processes?
                                                   ! if set to .false., the statistics are only printed by the ensemble mean group, which may save time

  LOGICAL               :: OBSDEP_OUT = .true.
  character(filelenmax) :: OBSDEP_OUT_BASENAME = 'obsdep'
  LOGICAL               :: OBSGUES_OUT = .false.                  !XXX not implemented yet...
  character(filelenmax) :: OBSGUES_OUT_BASENAME = 'obsgues.@@@@'  !XXX not implemented yet...
  LOGICAL               :: OBSANAL_OUT = .false.                  !XXX not implemented yet...
  character(filelenmax) :: OBSANAL_OUT_BASENAME = 'obsanal.@@@@'  !XXX not implemented yet...

  !--- PARAM_LETKF_RADAR
  logical :: USE_RADAR_REF       = .true.
  logical :: USE_RADAR_VR        = .true.
  logical :: USE_RADAR_PSEUDO_RH = .false.

  logical :: USE_OBSERR_RADAR_REF = .false.
  logical :: USE_OBSERR_RADAR_VR = .false.

  logical :: RADAR_OBS_4D = .false.

  REAL(r_size) :: RADAR_REF_THRES_DBZ = 15.0d0 !Threshold of rain/no rain
  INTEGER :: MIN_RADAR_REF_MEMBER = 1          !Ensemble members with reflectivity greather than RADAR_REF_THRES_DBZ
  INTEGER :: MIN_RADAR_REF_MEMBER_OBSREF = 1   !Ensemble members with

  REAL(r_size) :: MIN_RADAR_REF_DBZ = 0.0d0    !Minimum reflectivity
  REAL(r_size) :: LOW_REF_SHIFT = 0.0d0

  real(r_size) :: RADAR_ZMAX = 99.0d3          !Height limit of radar data to be used

  REAL(r_size) :: RADAR_PRH_ERROR = 0.1d0      !Obserational error for pseudo RH observations.

  !These 2 flags affects the computation of model reflectivity and radial velocity. 
  INTEGER :: INTERPOLATION_TECHNIQUE = 1
  INTEGER :: METHOD_REF_CALC = 3

  LOGICAL :: USE_TERMINAL_VELOCITY = .false.

  ! PARAMETERS FOR RADAR DATA ASSIMILATION
  INTEGER :: NRADARTYPE = 1  !Currently PAWR (1) and LIDAR (2) ... not used?

  !---PARAM_LETKF_H08
  logical :: H08_REJECT_LAND = .false. ! true: reject Himawari-8 radiance over the land
  logical :: H08_RTTOV_CLD = .true. ! true: all-sky, false: CSR in RTTOV fwd model
  real(r_size) :: H08_RTTOV_MINQ = 0.10d0 ! Threshold of water/ice contents for diagnosing cloud fraction (g m-3)
  real(r_size) :: H08_LIMIT_LEV = 20000.0d0 ! (Pa) Upper limit level of the sensitive height for Himawari-8 IR
  real(r_size) :: H08_RTTOV_CFRAC_CNST = 0.10d0 ! Denominator constant for diagnosing SEQUENTIAL(0-1) cloud fraction (g m-3)
                                                ! Negative values indicate DISCRETE (0/1) cloud fraction 
  real(r_size) :: H08_BT_MIN = 0.0d0 ! Lower limit of the BT for Himawari-8 IR
  real(r_size) :: H08_CLDSKY_THRS = -5.0d0 ! Threshold for diagnosing the sky condition using [BT(all-sky) - BT(clr)].
                                           ! Negative values: turn off
  integer :: H08_MIN_CLD_MEMBER = 1       ! If the number of the cloudy members is larger than H08_MIN_CLD_MEMBER,
                                           ! the first guess is diagnosed as cloudy. ! Not finished yet!
  integer :: H08_CH_USE(nch) = (/0,0,1,0,0,0,0,0,0,0/)
                        !! ch = (1,2,3,4,5,6,7,8,9,10)
                        !! (B07,B08,B09,B10,B11,B12,B13,B14,B15,B16)
                        !! ==1: Assimilate
                        !! ==0: NOT assimilate (rejected by QC in trans_XtoY_H08)
                        !! It is better to reject B11(ch=5) & B12(ch=6) obs because these bands are 
                        !! sensitive to chemicals.

  !--- PARAM_OBS_ERROR
  real(r_size) :: OBSERR_U = 1.0d0
  real(r_size) :: OBSERR_V = 1.0d0
  real(r_size) :: OBSERR_T = 1.0d0
  real(r_size) :: OBSERR_Q = 0.001d0
  real(r_size) :: OBSERR_RH = 10.0d0
  real(r_size) :: OBSERR_PS = 100.0d0
  real(r_size) :: OBSERR_RADAR_REF = 5.0d0
  real(r_size) :: OBSERR_RADAR_VR = 3.0d0
  real(r_size) :: OBSERR_TCX = 50.0d3 ! (m)
  real(r_size) :: OBSERR_TCY = 50.0d3 ! (m)
  real(r_size) :: OBSERR_TCP = 5.0d2 ! (Pa)
  real(r_size) :: OBSERR_H08(nch) = (/5.0d0,5.0d0,5.0d0,5.0d0,5.0d0,&
                                      5.0d0,5.0d0,5.0d0,5.0d0,5.0d0/) ! H08

  !--- PARAM_OBSSIM
  character(filelenmax) :: OBSSIM_IN_TYPE = 'history'
  character(filelenmax) :: OBSSIM_RESTART_IN_BASENAME = 'restart'
  character(filelenmax) :: OBSSIM_HISTORY_IN_BASENAME = 'history'
  character(filelenmax) :: OBSSIM_TOPO_IN_BASENAME = 'topo'
  integer               :: OBSSIM_TIME_START = 1
  integer               :: OBSSIM_TIME_END = 1
  character(filelenmax) :: OBSSIM_GRADS_OUT_NAME = ''
  integer               :: OBSSIM_NUM_3D_VARS = 0
  integer               :: OBSSIM_3D_VARS_LIST(nid_obs) = 0
  integer               :: OBSSIM_NUM_2D_VARS = 0
  integer               :: OBSSIM_2D_VARS_LIST(nid_obs) = 0
  real(r_size)          :: OBSSIM_RADAR_LON = 0.0d0
  real(r_size)          :: OBSSIM_RADAR_LAT = 0.0d0
  real(r_size)          :: OBSSIM_RADAR_Z = 0.0d0

contains
!-------------------------------------------------------------------------------
! PARAM_ENSEMBLE
!-------------------------------------------------------------------------------
subroutine read_nml_ensemble
  implicit none
  integer :: ierr
  
  namelist /PARAM_ENSEMBLE/ &
    MEMBER, &
    MEMBER_RUN, &
    MEMBER_ITER, &
    CONF_FILES, &
    CONF_FILES_SEQNUM, &
    DET_RUN, &
    DET_RUN_CYCLED

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_ENSEMBLE,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Error] /PARAM_ENSEMBLE/ is not found in namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_ENSEMBLE. Check!'
    stop
  endif

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_ENSEMBLE)
  end if

  return
end subroutine read_nml_ensemble

!-------------------------------------------------------------------------------
! PARAM_MODEL
!-------------------------------------------------------------------------------
subroutine read_nml_model
  implicit none
  integer :: ierr

  namelist /PARAM_MODEL/ &
    MODEL, &
    VERIFY_COORD

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_MODEL,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Warning] /PARAM_MODEL/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_MODEL. Check!'
    stop
  endif

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_MODEL)
  end if

  return
end subroutine read_nml_model

!-------------------------------------------------------------------------------
! PARAM_IO
!-------------------------------------------------------------------------------
!subroutine read_nml_io
!  implicit none
!  integer :: ierr

!  namelist /PARAM_IO/ &
!    IO_AGGREGATE

!  rewind(IO_FID_CONF)
!  read(IO_FID_CONF,nml=PARAM_IO,iostat=ierr)
!  if (ierr < 0) then !--- missing
!    write(6,*) '[Warning] /PARAM_IO/ is not found in namelist.'
!!    stop
!  elseif (ierr > 0) then !--- fatal error
!    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_IO. Check!'
!    stop
!  endif

!  if (LOG_LEVEL >= 2) then
!    write(6, nml=PARAM_IO)
!  end if

!  return
!end subroutine read_nml_io

!-------------------------------------------------------------------------------
! PARAM_LOG
!-------------------------------------------------------------------------------
subroutine read_nml_log
  implicit none
  integer :: ierr

  namelist /PARAM_LOG/ &
    LOG_LEVEL, &
    USE_MPI_BARRIER

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LOG,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Warning] /PARAM_LOG/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_LOG. Check!'
    stop
  endif

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_LOG)
  end if

  return
end subroutine read_nml_log

!-------------------------------------------------------------------------------
! PARAM_OBSOPE
!-------------------------------------------------------------------------------
subroutine read_nml_obsope
  implicit none
  integer :: ierr
  
  namelist /PARAM_OBSOPE/ &
    OBS_IN_NUM, &
    OBS_IN_NAME, &
    OBS_IN_FORMAT, &
    OBSDA_RUN, &
    OBSDA_OUT, &
    OBSDA_OUT_BASENAME, &
    OBSDA_MEAN_OUT_BASENAME, &
    OBSDA_MDET_OUT_BASENAME, &
    HISTORY_IN_BASENAME, &
    HISTORY_MEAN_IN_BASENAME, &
    HISTORY_MDET_IN_BASENAME, &
    SLOT_START, &
    SLOT_END, &
    SLOT_BASE, &
    SLOT_TINTERVAL

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_OBSOPE,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Error] /PARAM_OBSOPE/ is not found in namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_OBSOPE. Check!'
    stop
  endif

  if (trim(OBSDA_MEAN_OUT_BASENAME) == '') then
    call file_member_replace(0, OBSDA_OUT_BASENAME, OBSDA_MEAN_OUT_BASENAME, memf_mean)
  end if
  if (trim(OBSDA_MDET_OUT_BASENAME) == '') then
    call file_member_replace(0, OBSDA_OUT_BASENAME, OBSDA_MDET_OUT_BASENAME, memf_mdet)
  end if

  if (trim(HISTORY_MEAN_IN_BASENAME) == '') then
    call file_member_replace(0, HISTORY_IN_BASENAME, HISTORY_MEAN_IN_BASENAME, memf_mean)
  end if
  if (trim(HISTORY_MDET_IN_BASENAME) == '') then
    call file_member_replace(0, HISTORY_IN_BASENAME, HISTORY_MDET_IN_BASENAME, memf_mdet)
  end if

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_OBSOPE)
  end if

  return
end subroutine read_nml_obsope

!-------------------------------------------------------------------------------
! PARAM_LETKF
!-------------------------------------------------------------------------------
subroutine read_nml_letkf
  implicit none
  integer :: ierr
  
  namelist /PARAM_LETKF/ &
    OBSDA_IN, &
    OBSDA_IN_BASENAME, &
    OBSDA_MEAN_IN_BASENAME, &
    OBSDA_MDET_IN_BASENAME, &
    GUES_IN_BASENAME, &
    GUES_MEAN_INOUT_BASENAME, &
    GUES_MDET_IN_BASENAME, &
    GUES_SPRD_OUT, &
    GUES_SPRD_OUT_BASENAME, &
    ANAL_OUT_BASENAME, &
    ANAL_MEAN_OUT_BASENAME, &
    ANAL_MDET_OUT_BASENAME, &
    ANAL_SPRD_OUT, &
    ANAL_SPRD_OUT_BASENAME, &
    LETKF_TOPO_IN_BASENAME, &
    INFL_MUL, &
    INFL_MUL_MIN, &
    INFL_MUL_ADAPTIVE, &
    INFL_MUL_IN_BASENAME, &
    INFL_MUL_OUT_BASENAME, &
    INFL_ADD, &
    INFL_ADD_IN_BASENAME, &
    INFL_ADD_SHUFFLE, &
    INFL_ADD_Q_RATIO, &
    INFL_ADD_REF_ONLY, &
    RELAX_ALPHA, &
    RELAX_ALPHA_SPREAD, &
    RELAX_TO_INFLATED_PRIOR, &
    RELAX_SPREAD_OUT, &
    RELAX_SPREAD_OUT_BASENAME, &
    GROSS_ERROR, &
    GROSS_ERROR_RAIN, &
    GROSS_ERROR_RADAR_REF, &
    GROSS_ERROR_RADAR_VR, &
    GROSS_ERROR_RADAR_PRH, &
    GROSS_ERROR_H08, &
    GROSS_ERROR_TCX, &
    GROSS_ERROR_TCY, &
    GROSS_ERROR_TCP, &
    Q_UPDATE_TOP, &
    Q_SPRD_MAX, &
    BOUNDARY_BUFFER_WIDTH, &
    POSITIVE_DEFINITE_Q, &
    POSITIVE_DEFINITE_QHYD, &
    TC_SEARCH_DIS, &
    PS_ADJUST_THRES, &
    NOBS_OUT, &
    NOBS_OUT_BASENAME, &
    !*** for backward compatibility ***
    COV_INFL_MUL, &
    MIN_INFL_MUL, &
    ADAPTIVE_INFL_INIT, &
    BOUNDARY_TAPER_WIDTH

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Error] /PARAM_LETKF/ is not found in namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_LETKF. Check!'
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
  if (GROSS_ERROR_H08 < 0.0d0) then ! H08
    GROSS_ERROR_H08 = GROSS_ERROR
  end if
  if (GROSS_ERROR_TCX < 0.0d0) then
    GROSS_ERROR_TCX = GROSS_ERROR
  end if
  if (GROSS_ERROR_TCY < 0.0d0) then
    GROSS_ERROR_TCY = GROSS_ERROR
  end if
  if (GROSS_ERROR_TCP < 0.0d0) then
    GROSS_ERROR_TCP = GROSS_ERROR
  end if

  if (trim(OBSDA_MEAN_IN_BASENAME) == '') then
    call file_member_replace(0, OBSDA_IN_BASENAME, OBSDA_MEAN_IN_BASENAME, memf_mean)
  end if
  if (trim(OBSDA_MDET_IN_BASENAME) == '') then
    call file_member_replace(0, OBSDA_IN_BASENAME, OBSDA_MDET_IN_BASENAME, memf_mdet)
  end if

  if (trim(GUES_MEAN_INOUT_BASENAME) == '') then
    call file_member_replace(0, GUES_IN_BASENAME, GUES_MEAN_INOUT_BASENAME, memf_mean)
  end if
  if (trim(GUES_MDET_IN_BASENAME) == '') then
    call file_member_replace(0, GUES_IN_BASENAME, GUES_MDET_IN_BASENAME, memf_mdet)
  end if
  if (trim(GUES_SPRD_OUT_BASENAME) == '') then
    call file_member_replace(0, GUES_IN_BASENAME, GUES_SPRD_OUT_BASENAME, memf_sprd)
  end if
  if (trim(ANAL_MEAN_OUT_BASENAME) == '') then
    call file_member_replace(0, ANAL_OUT_BASENAME, ANAL_MEAN_OUT_BASENAME, memf_mean)
  end if
  if (trim(ANAL_MDET_OUT_BASENAME) == '') then
    call file_member_replace(0, ANAL_OUT_BASENAME, ANAL_MDET_OUT_BASENAME, memf_mdet)
  end if
  if (trim(ANAL_SPRD_OUT_BASENAME) == '') then
    call file_member_replace(0, ANAL_OUT_BASENAME, ANAL_SPRD_OUT_BASENAME, memf_sprd)
  end if

  if (trim(INFL_MUL_OUT_BASENAME) == '') then
    INFL_MUL_ADAPTIVE = .false.
  end if
  if (trim(INFL_ADD_IN_BASENAME) == '') then
    INFL_ADD = 0.0d0
  end if
  if (trim(RELAX_SPREAD_OUT_BASENAME) == '') then
    RELAX_SPREAD_OUT = .false.
  end if
  if (trim(NOBS_OUT_BASENAME) == '') then
    NOBS_OUT = .false.
  end if

  !*** for backward compatibility ***
  if (COV_INFL_MUL /= 1.0d0 .and. INFL_MUL == 1.0d0) then
    INFL_MUL = COV_INFL_MUL
  end if
  if (MIN_INFL_MUL /= 0.0d0 .and. INFL_MUL_MIN == 0.0d0) then
    INFL_MUL_MIN = MIN_INFL_MUL
  end if
  if (ADAPTIVE_INFL_INIT .and. (.not. INFL_MUL_ADAPTIVE)) then
    INFL_MUL_ADAPTIVE = ADAPTIVE_INFL_INIT
  end if
  if (BOUNDARY_TAPER_WIDTH /= 0.0d0 .and. BOUNDARY_BUFFER_WIDTH == 0.0d0) then
    BOUNDARY_BUFFER_WIDTH = BOUNDARY_TAPER_WIDTH
  end if

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_LETKF)
  end if

  return
end subroutine read_nml_letkf

!-------------------------------------------------------------------------------
! PARAM_LETKF_PRC
!-------------------------------------------------------------------------------
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
    write(6,*) '[Warning] /PARAM_LETKF_PRC/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_LETKF_PRC. Check!'
    stop
  endif

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_LETKF_PRC)
  end if

  return
end subroutine read_nml_letkf_prc

!-------------------------------------------------------------------------------
! PARAM_LETKF_OBS
!-------------------------------------------------------------------------------
subroutine read_nml_letkf_obs
  implicit none
  integer :: itype
  integer :: ierr

  namelist /PARAM_LETKF_OBS/ &
    USE_OBS, &
    HORI_LOCAL, &
    VERT_LOCAL, &
    TIME_LOCAL, &
    HORI_LOCAL_RADAR_OBSNOREF, &
    HORI_LOCAL_RADAR_VR, &
    VERT_LOCAL_RADAR_VR, &
    VERT_LOCAL_RAIN_BASE, &
    MAX_NOBS_PER_GRID, &
    MAX_NOBS_PER_GRID_CRITERION, &
    OBS_MIN_SPACING, &
    OBS_SORT_GRID_SPACING

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_OBS,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Warning] /PARAM_LETKF_OBS/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_LETKF_OBS. Check!'
    stop
  endif

  do itype = 2, nobtype
    if (HORI_LOCAL(itype) < 0.0d0) then
      HORI_LOCAL(itype) = HORI_LOCAL(1)
    end if
    if (VERT_LOCAL(itype) < 0.0d0) then
      VERT_LOCAL(itype) = VERT_LOCAL(1)
    end if
    if (TIME_LOCAL(itype) < 0.0d0) then
      TIME_LOCAL(itype) = TIME_LOCAL(1)
    end if

    if (MAX_NOBS_PER_GRID(itype) < 0) then
      MAX_NOBS_PER_GRID(itype) = MAX_NOBS_PER_GRID(1)
    end if

    if (MAX_NOBS_PER_GRID_CRITERION < 1 .or. MAX_NOBS_PER_GRID_CRITERION > 3) then
      write (6, '(A,I4)') "[Error] Unsupported 'MAX_NOBS_PER_GRID_CRITERION':", MAX_NOBS_PER_GRID_CRITERION
      stop 99
    end if

    if (OBS_MIN_SPACING(itype) <= 0.0d0) then
      OBS_MIN_SPACING(itype) = OBS_MIN_SPACING(1)
    end if
    if (OBS_SORT_GRID_SPACING(itype) < 0.0d0) then
      OBS_SORT_GRID_SPACING(itype) = OBS_SORT_GRID_SPACING(1)
    end if
  end do

  if (HORI_LOCAL_RADAR_OBSNOREF < 0.0d0) then
    HORI_LOCAL_RADAR_OBSNOREF = HORI_LOCAL(22) !PHARAD
  end if
  if (HORI_LOCAL_RADAR_VR < 0.0d0) then
    HORI_LOCAL_RADAR_VR = HORI_LOCAL(22) !PHARAD
  end if
  if (VERT_LOCAL_RADAR_VR < 0.0d0) then
    VERT_LOCAL_RADAR_VR = VERT_LOCAL(22) !PHARAD
  end if

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_LETKF_OBS)
  end if

  return
end subroutine read_nml_letkf_obs

!-------------------------------------------------------------------------------
! PARAM_LETKF_VAR_LOCAL
!-------------------------------------------------------------------------------
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
    write(6,*) '[Warning] /PARAM_LETKF_VAR_LOCAL/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_LETKF_VAR_LOCAL. Check!'
    stop
  endif

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_LETKF_VAR_LOCAL)
  end if

  return
end subroutine read_nml_letkf_var_local

!-------------------------------------------------------------------------------
! PARAM_LETKF_MONITOR
!-------------------------------------------------------------------------------
subroutine read_nml_letkf_monitor
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_MONITOR/ &
    DEPARTURE_STAT, &
    DEPARTURE_STAT_RADAR, &
    DEPARTURE_STAT_H08, &
    DEPARTURE_STAT_T_RANGE, &
    DEPARTURE_STAT_ALL_PROCESSES, &
    OBSDEP_OUT, &
    OBSDEP_OUT_BASENAME, &
    OBSGUES_OUT, &
    OBSGUES_OUT_BASENAME, &
    OBSANAL_OUT, &
    OBSANAL_OUT_BASENAME

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_MONITOR,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Warning] /PARAM_LETKF_MONITOR/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_LETKF_MONITOR. Check!'
    stop
  endif

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_LETKF_MONITOR)
  end if

  return
end subroutine read_nml_letkf_monitor

!-------------------------------------------------------------------------------
! PARAM_LETKF_RADAR
!-------------------------------------------------------------------------------
subroutine read_nml_letkf_radar
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_RADAR/ &
    USE_RADAR_REF, &
    USE_RADAR_VR, &
    USE_RADAR_PSEUDO_RH, &
    USE_OBSERR_RADAR_REF, &
    USE_OBSERR_RADAR_VR, &
    RADAR_OBS_4D, &
    RADAR_REF_THRES_DBZ, &
    MIN_RADAR_REF_MEMBER, &
    MIN_RADAR_REF_MEMBER_OBSREF, &
    MIN_RADAR_REF_DBZ, &
    LOW_REF_SHIFT, &
    RADAR_ZMAX, &
    RADAR_PRH_ERROR, &
    INTERPOLATION_TECHNIQUE, &
    METHOD_REF_CALC, &
    USE_TERMINAL_VELOCITY, &
    NRADARTYPE

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_RADAR,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Warning] /PARAM_LETKF_RADAR/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_LETKF_RADAR. Check!'
    stop
  endif

  if (RADAR_REF_THRES_DBZ < MIN_RADAR_REF_DBZ) then
    RADAR_REF_THRES_DBZ = MIN_RADAR_REF_DBZ
  end if

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_LETKF_RADAR)
  end if

  return
end subroutine read_nml_letkf_radar

!-------------------------------------------------------------------------------
! PARAM_LETKF_H08
!-------------------------------------------------------------------------------
subroutine read_nml_letkf_h08
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_H08/ &
    H08_REJECT_LAND, &
    H08_RTTOV_CLD, &
    H08_MIN_CLD_MEMBER, &
    H08_CLDSKY_THRS, &
    H08_RTTOV_MINQ, &
    H08_RTTOV_CFRAC_CNST, &
    H08_LIMIT_LEV, &
    H08_BT_MIN, &
    H08_CH_USE

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_H08,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Warning] /PARAM_LETKF_H08/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_LETKF_H08. Check!'
    stop
  endif

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_LETKF_H08)
  end if

  return
end subroutine read_nml_letkf_h08

!-------------------------------------------------------------------------------
! PARAM_OBS_ERROR
!-------------------------------------------------------------------------------
subroutine read_nml_obs_error
  implicit none
  integer :: ierr

  namelist /PARAM_OBS_ERROR/ &
    OBSERR_U, &
    OBSERR_V, &
    OBSERR_T, &
    OBSERR_Q, &
    OBSERR_RH, &
    OBSERR_PS, &
    OBSERR_RADAR_REF, &
    OBSERR_RADAR_VR, &
    OBSERR_TCX, &
    OBSERR_TCY, &
    OBSERR_TCP, &
    OBSERR_H08    ! H08

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_OBS_ERROR,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Warning] /PARAM_OBS_ERROR/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_OBS_ERROR. Check!'
    stop
  endif

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_OBS_ERROR)
  end if

  return
end subroutine read_nml_obs_error

!-------------------------------------------------------------------------------
! PARAM_OBSSIM
!-------------------------------------------------------------------------------
subroutine read_nml_obssim
  implicit none
  integer :: ierr

  namelist /PARAM_OBSSIM/ &
    OBSSIM_IN_TYPE, &
    OBSSIM_RESTART_IN_BASENAME, &
    OBSSIM_HISTORY_IN_BASENAME, &
    OBSSIM_TOPO_IN_BASENAME, &
    OBSSIM_TIME_START, &
    OBSSIM_TIME_END, &
    OBSSIM_GRADS_OUT_NAME, &
    OBSSIM_NUM_3D_VARS, &
    OBSSIM_3D_VARS_LIST, &
    OBSSIM_NUM_2D_VARS, &
    OBSSIM_2D_VARS_LIST, &
    OBSSIM_RADAR_LON, &
    OBSSIM_RADAR_LAT, &
    OBSSIM_RADAR_Z

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_OBSSIM,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) '[Warning] /PARAM_OBSSIM/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) '[Error] xxx Not appropriate names in namelist PARAM_OBSSIM. Check!'
    stop
  endif

  if (trim(OBSSIM_GRADS_OUT_NAME) == '') then
    if (trim(OBSSIM_IN_TYPE) == 'restart') then
      OBSSIM_GRADS_OUT_NAME = trim(OBSSIM_RESTART_IN_BASENAME) // '.grd'
    else if (trim(OBSSIM_IN_TYPE) == 'history') then
      OBSSIM_GRADS_OUT_NAME = trim(OBSSIM_HISTORY_IN_BASENAME) // '.grd'
    end if
  end if

  if (LOG_LEVEL >= 2) then
    write(6, nml=PARAM_OBSSIM)
  end if

  return
end subroutine read_nml_obssim

!-------------------------------------------------------------------------------
! Replace the member notation by the formatted member string
! * will be wrong if memflen /= 4
!-------------------------------------------------------------------------------
! [INPUT]
!   mem          : member number
!   filename     : input filename string
!   str          : (optional) use this formatted member string if mem <= 0
! [OUTPUT]
!   filename_out : output filename string with the member notation replaced
!-------------------------------------------------------------------------------
subroutine file_member_replace(mem, filename, filename_out, memfstr)
  implicit none
  integer, intent(in) :: mem
  character(len=*), intent(in) :: filename
  character(len=*), intent(out) :: filename_out
  character(len=memflen), intent(in), optional :: memfstr
  integer :: s, is

  s = 0
  filename_out = filename
  do is = 1, len(filename)-memflen+1
    if (filename(is:is+memflen-1) == memf_notation) then
      if (mem >= 1) then
        write (filename_out(is:is+memflen-1), '(I4.4)') mem
      else if (present(memfstr)) then
        write (filename_out(is:is+memflen-1), '(A4)') memfstr
      end if
      s = is
      exit
    end if
  end do

  if (s == 0) then
    write (6, '(3A)') "[Warning] Keyword '@@@@' not found in '", filename, "'"
    stop 99
  end if

  return
end subroutine file_member_replace

!===============================================================================
end module common_nml
