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
  integer, parameter :: nv3d = 22    ! number of 3D state variables (in SCALE restart files)
  integer, parameter :: nv2d = 0     ! number of 2D state variables (in SCALE restart files)
  integer, parameter :: nid_obs = 18 ! number of variable types
  integer, parameter :: nobtype = 25 ! number of observation report types
  integer, parameter :: nch = 10     ! H08 Num of Himawari-8 (IR) channels
  integer, parameter :: NIRB_HIM8 = 10     ! H08 Num of Himawari-8 (IR) bands

  integer, parameter :: nobsfilemax = 10
  integer, parameter :: filelenmax = 256
  integer, parameter :: memberflen = 4 ! Length of member # in filename

  !--- PARAM_ENSEMBLE
  integer :: MEMBER = 3      ! ensemble size
  integer :: MEMBER_RUN = 1  !
  integer :: MEMBER_ITER = 0 !

  !--- PARAM_OBSOPE
  integer               :: OBS_IN_NUM = 1
  character(filelenmax) :: OBS_IN_NAME(nobsfilemax) = 'obs.dat'
  integer               :: OBS_IN_FORMAT(nobsfilemax) = 1
  logical               :: OBSDA_RUN(nobsfilemax) = .true.
  logical               :: OBSDA_OUT = .true.
  character(filelenmax) :: OBSDA_OUT_BASENAME = 'obsda.@@@@'

  character(filelenmax) :: HISTORY_IN_BASENAME = 'hist.@@@@'

  integer               :: SLOT_START = 1
  integer               :: SLOT_END = 1
  integer               :: SLOT_BASE = 1
  real(r_size)          :: SLOT_TINTERVAL = 3600.0d0

  !--- PARAM_LETKF
  integer               :: MAX_UPDATE_ZNUM = 100
  logical               :: QC_SIGB = .false.
  logical               :: WRITE_GRADS_SPRD = .true.
  logical               :: WRITE_GRADS_MEAN = .true.
  logical               :: OBSDA_IN = .false.
  character(filelenmax) :: OBSDA_IN_BASENAME = 'obsda.@@@@'
  character(filelenmax) :: GUES_IN_BASENAME = 'gues.@@@@'
  character(filelenmax) :: GUES_OUT_MEAN_BASENAME = 'gues.mean'
  character(filelenmax) :: GUES_OUT_SPRD_BASENAME = 'gues.sprd'
  character(filelenmax) :: ANAL_OUT_BASENAME = 'anal.@@@@'
  character(filelenmax) :: ANAL_OUT_MEAN_BASENAME = 'anal.mean'
  character(filelenmax) :: ANAL_OUT_SPRD_BASENAME = 'anal.sprd'
  character(filelenmax) :: LETKF_TOPO_IN_BASENAME = 'topo'  !!!!!! -- directly use the SCALE namelist --???? !!!!!!

  real(r_size) :: INFL_MUL = 1.0d0           ! >  0: globally constant covariance inflation
                                             ! <= 0: use 3D inflation field from 'INFL_MUL_IN_BASENAME' file
  real(r_size) :: INFL_MUL_MIN = -1.0d0      ! minimum inlfation factor (<= 0: not used)
  logical :: INFL_MUL_ADAPTIVE = .false.     ! if true, outout adaptively estimated 3D inlfation field to 'INFL_MUL_OUT_BASENAME' file
  character(filelenmax) :: INFL_MUL_IN_BASENAME = 'infl'
  character(filelenmax) :: INFL_MUL_OUT_BASENAME = 'infl'

  real(r_size) :: INFL_ADD = 0.0d0           ! additive inflation
  character(filelenmax) :: INFL_ADD_IN_BASENAME = 'addi.@@@@'

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
  real(r_size) :: GROSS_ERROR_TCX = -1.0d0 !  < 0: same as GROSS_ERROR 
  real(r_size) :: GROSS_ERROR_TCY = -1.0d0 !  < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_TCP = -1.0d0 !  < 0: same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_LT = -1.0d0 ! < 0: same as GROSS_ERROR

  real(r_size) :: Q_UPDATE_TOP = 0.0d0     ! water vapor and hydrometeors are updated only below this pressure level (Pa)
  real(r_size) :: Q_UPDATE_LOW_HEIGHT = -99.0d0     ! water vapor is updated only above this pressure level (m)
  real(r_size) :: Q_SPRD_MAX = -1.0D0      ! maximum q (ensemble spread)/(ensemble mean) (only effective when > 0)

  real(r_size) :: BOUNDARY_BUFFER_WIDTH = 0.0d0

  logical :: POSITIVE_DEFINITE_Q = .false.
  logical :: POSITIVE_DEFINITE_QHYD = .false.
  logical :: POSITIVE_DEFINITE_QHYD_QCRG = .true. ! Modify qcharge if qhyd is zero
  logical :: N_LOG_TRANS = .false. ! Log transformation for Nx
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
  ! <0: same as HORI_LOCAL_SIGMA(1)
  real(r_size) :: HORI_LOCAL(nobtype) = &
    (/500.0d3, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
       -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
       -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0/)

  ! >0: localization length scale [ln(p) or m depends on obstype]
  !  0: no localization
  ! <0: same as VERTIME_LOCAL(1)
  real(r_size) :: VERT_LOCAL(nobtype) = &
    (/ 0.4d0,   -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0,   -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, 1000.0d0, -1.0d0, -1.0d0, -1.0d0/)
!      -1.0d0, 1000.0d0, -1.0d0,  0.0d0/)

  ! >0: localization length scale (sec) XXX not implemented yet XXX
  !  0: no localization
  ! <0: same as TIME_LOCAL_SIGMA(1)
  real(r_size) :: TIME_LOCAL(nobtype) = &
    (/ 0.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0/)

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
      -1, -1, -1, -1, -1/)

  integer :: MAX_NOBS_PER_GRID_CRITERION = 1 ! 1: normalized 3D distance (from closest)
                                             ! 2: localization weight (from largest)
                                             ! 3: weighted observation error variance (from smallest)

  ! >0: typical minimum spacing of the obsetvation types in the densest observed area (not tuned carefully yet)
  !     *this is only used for automatically determine OBS_SORT_GRID_SPACING. if using pre-set OBS_SORT_GRID_SPACING, this has no effect.
  ! <=0: same as OBS_MIN_SPACING(1)
  real(r_size) :: OBS_MIN_SPACING(nobtype) = &
    (/300.0d3, 100.0d3, 100.0d3, 150.0d3, 300.0d3, 150.0d3, 150.0d3, 100.0d3, 150.0d3, 150.0d3, &
      150.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 300.0d3, 150.0d3, 150.0d3, 150.0d3, &
      150.0d3,   1.0d3,  15.0d3,1000.0d3, 5.0d2/)

  ! >0: optimal grid spacing for bucket sorting of observations
  !  0: automatically determined based on HORI_LOCAL, MAX_NOBS_PER_GRID, and OBS_MIN_SPACING
  ! <0: same as OBS_SORT_GRID_SPACING(1)
  real(r_size) :: OBS_SORT_GRID_SPACING(nobtype) = &
    (/ 0.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0/)

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
  real(r_size) :: VAR_LOCAL_LT(nv3d+nv2d)       = 1.0d0 ! LT

  !--- PARAM_LETKF_MONITOR
  logical :: DEPARTURE_STAT = .true.
  logical :: DEPARTURE_STAT_RADAR = .false.
  logical :: DEPARTURE_STAT_H08 = .false.
  real(r_size) :: DEPARTURE_STAT_T_RANGE = 0.0d0 ! time range within which observations are considered in the departure statistics.
                                                 ! 0: no limit

  LOGICAL :: OMB_OUTPUT = .true.
  LOGICAL :: OMA_OUTPUT = .true.
  LOGICAL :: OBSGUES_OUTPUT = .false.
  LOGICAL :: OBSANAL_OUTPUT = .false.

  !--- PARAM_LETKF_LT
  character(4) :: LT_OBS_NANE = 'FP3D'
  logical :: LT_NOB_ONLY = .false. ! Use only if O-B < 0
  logical :: LT_POB_ONLY = .false. ! Use only if O-B < 0
  logical :: LT_TEST_SINGLE = .false. ! Sigle-observation experiment? 
  integer :: LT_TEST_SINGLE_I = 1 ! grid index for single obs test
  integer :: LT_TEST_SINGLE_J = 1 ! gird index for single obs test
  real :: LT_TEST_SINGLE_LON = 175.0 
  real :: LT_TEST_SINGLE_LAT = 183.0
  real(r_size) :: LT_OBSERR_GROSS = 1.0d0 ! obs error for gross-error check
  integer :: MIN_LT_MEMBER_OBSON = 1
  integer :: MIN_LT_MEMBER_OBSOFF = 1
  real(r_size) :: LT_ZMAX = 12.0d3 ! !Height limit of lightning data to be used
  integer :: XY_THINNING_LT = 1 ! Horizontal thinning interval (1: no thinning)
  integer :: Z_THINNING_LT = 1 ! Vertical thinning interval (1: no thinning)
  logical :: LT_LOG = .false. ! Log transformation
  real(r_size) :: LT_LOG_CONST = 1.0d0 ! constant for log transformation
  real(r_size) :: LT_LOG_OERR = 1.0d0 ! obs error for log transformation
  logical :: LT_2DLOC = .false. ! Assimilate maximum loc (2D)
  real(r_size) :: LT_2DLOC_OERR = 40.0d3 ! obs error for 2d location (m)
  real(r_size) :: LT_ON_THRS = 0.0d0 ! threashold for flash on/off
  logical :: USE_GT = .false. ! Gaussian transformation
  character(filelenmax) :: CDF_FP_FILENAME = '' ! CDF for flash point
  real(r_size) :: LT_GT_OERR = 0.5d0 ! obs err for lt with gaussian transform

  !--- PARAM_LETKF_RADAR
  logical :: USE_RADAR_REF       = .true.
  logical :: USE_RADAR_VR        = .true.
  logical :: USE_RADAR_PSEUDO_RH = .false.

  logical :: USE_OBSERR_RADAR_REF = .false.
  logical :: USE_OBSERR_RADAR_VR = .false.

  REAL(r_size) :: RADAR_REF_THRES_DBZ = 15.0d0 !Threshold of rain/no rain
  INTEGER :: MIN_RADAR_REF_MEMBER = 1          !Ensemble members with reflectivity greather than RADAR_REF_THRES_DBZ
  INTEGER :: MIN_RADAR_REF_MEMBER_OBSREF = 1   !Ensemble members with

  REAL(r_size) :: MIN_RADAR_REF_DBZ = 0.0d0    !Minimum reflectivity
  REAL(r_size) :: LOW_REF_SHIFT = 0.0d0

  real(r_size) :: RADAR_ZMAX = 99.0d3          !Height limit of radar data to be used

  REAL(r_size) :: RADAR_PRH_ERROR = 0.1d0      !Obserational error for pseudo RH observations.
  logical :: USE_RADAR_METHOD3_MELT = .true.

  !These 2 flags affects the computation of model reflectivity and radial velocity. 
  INTEGER :: INTERPOLATION_TECHNIQUE = 1
  INTEGER :: METHOD_REF_CALC = 3

  LOGICAL :: USE_TERMINAL_VELOCITY = .false.
  logical :: RADAR_IDEAL = .true.
  logical :: ADD_OBSERR_TRUTH = .true. ! Add noise to the truth?

  ! PARAMETERS FOR RADAR DATA ASSIMILATION
  INTEGER :: NRADARTYPE = 1  !Currently PAWR (1) and LIDAR (2) ... not used?

  !---PARAM_LETKF_H08
  character(filelenmax) :: H08_RTTOV_COEF_PATH = '.'
  real(r_size) :: H08_RTTOV_MINQ_CTOP = 0.10d0 ! Threshold of water/ice contents for diagnosing the cloud top (g m-3)

  logical :: H08_RTTOV_PROF_SHIFT = .false. ! true: shift the climatological profile above the model top 
                                            !       (equivalent to extrapolate
                                            !       by using the climatological
                                            !       lapse rate)
                                            ! false: relax the original (model)
                                            ! profiles above [H08_RTTOV_RLX_HGT]
                                            ! m back to the climatological
                                            ! profile 
  integer :: H08_RTTOV_KADD = 0
  integer :: H08_RTTOV_CFRAC =  1 ! cloud fraction diagnosis 
                                  ! 0: using H08_RTTOV_CFRAC_CNST following
                                  ! Honda et al. (2017a,b)
                                  ! 1: SCALE method as of 11/15/2017 with a
                                  ! minor modification (excluding qr)
                                  ! 2: Tompkins and Janiskova (2004QJRMS) method
                                  ! (as in Okamoto 2017QJRMS)

  integer :: H08_THIN_NG = 10 ! 20 km
  real(r_size) :: H08_OBSERR_TRUE = 3.0 ! true obs error for OSSE
  real(r_size) :: H08_OBSERR_RUN_CLR = 3.0 ! obs error for OSSE
  real(r_size) :: H08_OBSERR_RUN_CLD = 6.0 ! obs error for OSSE
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
  real(r_size) :: OBSERR_LT3D = 1.0d0 ! tentative!!
  real(r_size) :: OBSERR_LT2D = 1.0d0 ! tentative!!
  real(r_size) :: OBSERR_FP_RAT = 0.5d0 ! 50% of observed value: Lien et al. (2013)
  real(r_size) :: OBSERR_FP = 1.0d0 
  real(r_size) :: OBSERR_FP_TRUE = 1.0d0
  real(r_size) :: OBSERR_FP_OBSON = 1.0d0  ! obs err for LETKF when obs has lightning
  real(r_size) :: OBSERR_FP_OBSOFF = 0.1d0  ! obs err for LETKF when obs has no lightning
  real(r_size) :: OBSERR_FP_LOC2D_MAX = 1.0d0  ! intensity obs
  real(r_size) :: OBSERR_FP_LOC2D_LL = 40.0d3  ! 2d location obs

  !--- PARAM_OBSSIM
  logical               :: OBSSIM_OBSOUT = .false.
  character(filelenmax) :: OBSSIM_IN_TYPE = 'history'
  character(filelenmax) :: OBSSIM_RESTART_IN_BASENAME = 'restart'
  character(filelenmax) :: OBSSIM_HISTORY_IN_BASENAME = 'history'
  character(filelenmax) :: OBSSIM_TOPO_IN_BASENAME = 'topo'
  integer               :: OBSSIM_TIME_START = 1
  integer               :: OBSSIM_TIME_END = 1
  integer               :: OBSSIM_TIME_INT = 30
  character(filelenmax) :: OBSSIM_GRADS_OUT_NAME = ''
  integer               :: OBSSIM_NUM_3D_VARS = 0
  integer               :: OBSSIM_3D_VARS_LIST(100) = 0
  integer               :: OBSSIM_NUM_2D_VARS = 0
  integer               :: OBSSIM_2D_VARS_LIST(100) = 0
  real(r_size)          :: OBSSIM_RADAR_LON = 0.0d0
  real(r_size)          :: OBSSIM_RADAR_LAT = 0.0d0
  real(r_size)          :: OBSSIM_RADAR_Z = 0.0d0
  character(3) :: OBSSIM_RADAR_LONc = ''
  character(3) :: OBSSIM_RADAR_LATc = ''
  character(3) :: OBSSIM_RADAR_Zc = ''
  logical               :: OBSSIM_RADAR_ERR_10 = .true.
  integer               :: OBSSIM_RADAR_CLR_THIN = 1
  integer               :: OBSSIM_RADAR_CLR_ZTHIN = 1
  integer               :: OBSSIM_RADAR_RAIN_THIN = 1
  integer               :: OBSSIM_RADAR_RAIN_ZTHIN = 1
  integer               :: OBSSIM_RADAR_VR_THIN = -1
  integer               :: OBSSIM_RADAR_VR_ZTHIN = -1
  real(r_size)          :: OBSSIM_RADAR_RANGE = 60.0d3
  character(filelenmax) :: OBSSIM_OBSOUT_FNAME = ''

  !--- PARAM_SC
  ! inital perturbation will be added within a SC_DIST_ECHO [m] circle centered at the first echo in nature run
  real(r_size) :: SC_DIST_ECHO = 20.0d3 ! [m] 
  real(r_size) :: SC_PERT_COEF = 1.0d0 ! coefficient for perturbation
  real(r_size) :: SC_FECHO_X = 40.0d3 ! [m]  ! first echo location
  real(r_size) :: SC_FECHO_Y = 40.0d3 ! [m]  ! first echo location
  character(filelenmax) :: SC_NATURE_IN_BASENAME = ''
  character(filelenmax) :: SC_OUT_BASENAME = ''

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
! PARAM_ENSEMBLE
!-----------------------------------------------------------------------
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
    HISTORY_IN_BASENAME, &
    SLOT_START, &
    SLOT_END, &
    SLOT_BASE, &
    SLOT_TINTERVAL

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_OBSOPE,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Error: /PARAM_OBSOPE/ is not found in namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_OBSOPE. Check!'
    stop
  endif

  write(6, nml=PARAM_OBSOPE)

  return
end subroutine read_nml_obsope

!-----------------------------------------------------------------------
! PARAM_LETKF
!-----------------------------------------------------------------------
subroutine read_nml_letkf
  implicit none
  integer :: ierr
  
  namelist /PARAM_LETKF/ &
    MAX_UPDATE_ZNUM, &
    QC_SIGB, &
    WRITE_GRADS_MEAN, &
    WRITE_GRADS_SPRD, &
    OBSDA_IN, &
    OBSDA_IN_BASENAME, &
    GUES_IN_BASENAME, &
    GUES_OUT_MEAN_BASENAME, &
    GUES_OUT_SPRD_BASENAME, &
    ANAL_OUT_BASENAME, &
    ANAL_OUT_MEAN_BASENAME, &
    ANAL_OUT_SPRD_BASENAME, &
    LETKF_TOPO_IN_BASENAME, &
    INFL_MUL, &
    INFL_MUL_MIN, &
    INFL_MUL_ADAPTIVE, &
    INFL_MUL_IN_BASENAME, &
    INFL_MUL_OUT_BASENAME, &
    INFL_ADD, &
    INFL_ADD_IN_BASENAME, &
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
    GROSS_ERROR_LT, &
    Q_UPDATE_TOP, &
    Q_UPDATE_LOW_HEIGHT, &
    Q_SPRD_MAX, &
    BOUNDARY_BUFFER_WIDTH, &
    POSITIVE_DEFINITE_Q, &
    POSITIVE_DEFINITE_QHYD, &
    POSITIVE_DEFINITE_QHYD_QCRG, &
    N_LOG_TRANS, &
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
  if (GROSS_ERROR_LT < 0.0d0) then
    GROSS_ERROR_LT = GROSS_ERROR
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
! PARAM_LETKF_OBS
!-----------------------------------------------------------------------
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
    write(6,*) 'Warning: /PARAM_LETKF_OBS/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_OBS. Check!'
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

  write(6, nml=PARAM_LETKF_OBS)

  return
end subroutine read_nml_letkf_obs

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
    VAR_LOCAL_H08, & ! H08
    VAR_LOCAL_LT ! LT

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
    USE_RADAR_REF, &
    USE_RADAR_VR, &
    USE_RADAR_PSEUDO_RH, &
    USE_OBSERR_RADAR_REF, &
    USE_OBSERR_RADAR_VR, &
    RADAR_REF_THRES_DBZ, &
    MIN_RADAR_REF_MEMBER, &
    MIN_RADAR_REF_MEMBER_OBSREF, &
    MIN_RADAR_REF_DBZ, &
    LOW_REF_SHIFT, &
    RADAR_ZMAX, &
    USE_RADAR_METHOD3_MELT, &
    RADAR_PRH_ERROR, &
    INTERPOLATION_TECHNIQUE, &
    METHOD_REF_CALC, &
    USE_TERMINAL_VELOCITY, &
    RADAR_IDEAL, &
    ADD_OBSERR_TRUTH, &
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
! PARAM_LETKF_H08
!-----------------------------------------------------------------------
subroutine read_nml_letkf_h08
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_H08/ &
    H08_REJECT_LAND, &
    H08_OBSERR_TRUE, &
    H08_OBSERR_RUN_CLR, &
    H08_OBSERR_RUN_CLD, &
    H08_THIN_NG, &
    H08_RTTOV_CLD, &
    H08_MIN_CLD_MEMBER, &
    H08_CLDSKY_THRS, &
    H08_RTTOV_MINQ, &
    H08_RTTOV_CFRAC, &
    H08_RTTOV_CFRAC_CNST, &
    H08_LIMIT_LEV, &
    H08_BT_MIN, &
    H08_RTTOV_COEF_PATH, &
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

!-----------------------------------------------------------------------
! PARAM_OBS_ERROR
!-----------------------------------------------------------------------
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
    OBSERR_H08, &    ! H08
    OBSERR_LT3D, &
    OBSERR_LT2D, &
    OBSERR_FP_RAT, &
    OBSERR_FP_OBSOFF, &
    OBSERR_FP_OBSON, &
    OBSERR_FP_TRUE, &
    OBSERR_FP_LOC2D_MAX, &
    OBSERR_FP_LOC2D_LL,  &
    OBSERR_FP

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_OBS_ERROR,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_OBS_ERROR/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_OBS_ERROR. Check!'
    stop
  endif

  write(6, nml=PARAM_OBS_ERROR)

  return
end subroutine read_nml_obs_error

!-----------------------------------------------------------------------
! PARAM_OBSSIM
!-----------------------------------------------------------------------
subroutine read_nml_obssim
  implicit none
  integer :: ierr

  namelist /PARAM_OBSSIM/ &
    OBSSIM_OBSOUT, &
    OBSSIM_IN_TYPE, &
    OBSSIM_RESTART_IN_BASENAME, &
    OBSSIM_HISTORY_IN_BASENAME, &
    OBSSIM_TOPO_IN_BASENAME, &
    OBSSIM_TIME_START, &
    OBSSIM_TIME_END, &
    OBSSIM_TIME_INT, &
    OBSSIM_GRADS_OUT_NAME, &
    OBSSIM_NUM_3D_VARS, &
    OBSSIM_3D_VARS_LIST, &
    OBSSIM_NUM_2D_VARS, &
    OBSSIM_2D_VARS_LIST, &
    OBSSIM_RADAR_LON, &
    OBSSIM_RADAR_LAT, &
    OBSSIM_RADAR_Z,   &
    OBSSIM_RADAR_ERR_10, &
    OBSSIM_RADAR_CLR_THIN, &
    OBSSIM_RADAR_RAIN_THIN, &
    OBSSIM_RADAR_CLR_ZTHIN, &
    OBSSIM_RADAR_VR_THIN, &
    OBSSIM_RADAR_VR_ZTHIN, &
    OBSSIM_RADAR_RAIN_ZTHIN, &
    OBSSIM_RADAR_RANGE, &
    OBSSIM_OBSOUT_FNAME

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_OBSSIM,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_OBSSIM/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_OBSSIM. Check!'
    stop
  endif

  if (trim(OBSSIM_GRADS_OUT_NAME) == '') then
    if (trim(OBSSIM_IN_TYPE) == 'restart') then
      OBSSIM_GRADS_OUT_NAME = trim(OBSSIM_RESTART_IN_BASENAME) // '.grd'
    else if (trim(OBSSIM_IN_TYPE) == 'history') then
      OBSSIM_GRADS_OUT_NAME = trim(OBSSIM_HISTORY_IN_BASENAME) // '.grd'
    end if
  end if

  if ( OBSSIM_RADAR_VR_THIN < 0 ) then
    OBSSIM_RADAR_VR_THIN = OBSSIM_RADAR_RAIN_THIN
  endif

  if ( OBSSIM_RADAR_VR_ZTHIN < 0 ) then
    OBSSIM_RADAR_VR_ZTHIN = OBSSIM_RADAR_RAIN_ZTHIN
  endif

  write(OBSSIM_RADAR_LONc,'(i3.3)')int(OBSSIM_RADAR_LON / 1000.0d0)
  write(OBSSIM_RADAR_LATc,'(i3.3)')int(OBSSIM_RADAR_LAT / 1000.0d0)
  write(OBSSIM_RADAR_Zc,'(i3.3)')int(OBSSIM_RADAR_Z / 1000.0d0)

  write(6, nml=PARAM_OBSSIM)

  return
end subroutine read_nml_obssim

!-----------------------------------------------------------------------
! PARAM_SC
!-----------------------------------------------------------------------
subroutine read_nml_sc
  implicit none
  integer :: ierr

  real(r_size) :: SC_DIST_ECHO = 20.0d3 ! [m] 
  real(r_size) :: SC_PERT_COEF = 1.0d0 ! coefficient for perturbation
  real(r_size) :: SC_FECHO_X = 40.0d3 ! [m]  ! first echo location
  real(r_size) :: SC_FECHO_Y = 40.0d3 ! [m]  ! first echo location
  character(filelenmax) :: SC_NATURE_IN_BASENAME = ''
  character(filelenmax) :: SC_OUT_BASENAME = ''
  namelist /PARAM_SC/ &
    SC_DIST_ECHO, &
    SC_PERT_COEF, &
    SC_FECHO_X, &
    SC_FECHO_Y, &
    SC_NATURE_IN_BASENAME, &
    SC_OUT_BASENAME

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_SC,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_SC/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_SC. Check!'
    stop
  endif

  write(6, nml=PARAM_SC)

  return
end subroutine read_nml_sc

!-----------------------------------------------------------------------
! PARAM_LETKF_LT
!-----------------------------------------------------------------------
subroutine read_nml_letkf_lt
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_LT/ &
    LT_NOB_ONLY, &
    LT_POB_ONLY, &
    LT_TEST_SINGLE, &
    LT_TEST_SINGLE_I, &
    LT_TEST_SINGLE_J, &
    LT_TEST_SINGLE_LON, &
    LT_TEST_SINGLE_LAT, &
    LT_OBS_NANE, &
    LT_OBSERR_GROSS, &
    MIN_LT_MEMBER_OBSON, &
    MIN_LT_MEMBER_OBSOFF, &
    XY_THINNING_LT, &
    Z_THINNING_LT, &
    LT_ZMAX, &
    LT_LOG, &
    LT_LOG_CONST, &
    LT_LOG_OERR, &
    LT_2DLOC, &
    LT_2DLOC_OERR, &
    USE_GT, &
    LT_GT_OERR, &
    CDF_FP_FILENAME

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_LT,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_LETKF_LT/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_LT. Check!'
    stop
  endif

  if ( LT_LOG ) then
    LT_ON_THRS = dlog( real(LT_LOG_CONST) + 0.0d0 )
  endif

  write(6, nml=PARAM_LETKF_LT)

  return
end subroutine read_nml_letkf_lt

!-----------------------------------------------------------------------
! file_member_replace
!-----------------------------------------------------------------------
subroutine file_member_replace(mem, filename, filename_out)
  implicit none
  integer, intent(in) :: mem
  character(*), intent(in) :: filename
  character(filelenmax), intent(out) :: filename_out

  character(memberflen) :: memberfstr = '@@@@'
  integer :: s, is

  s = 0
  filename_out = filename
  do is = 1, len(filename)-memberflen+1
    if (filename(is:is+memberflen-1) == memberfstr) then
      if (mem <= MEMBER) then
        write (filename_out(is:is+memberflen-1), '(I4.4)') mem
      else if (mem == MEMBER+1) then
        write (filename_out(is:is+memberflen-1), '(A4)') 'mean'  !!!!!! will be wrong if memberflen != 4 !!!!!!
      else if (mem == MEMBER+2) then
        write (filename_out(is:is+memberflen-1), '(A4)') 'sprd'  !!!!!! will be wrong if memberflen != 4 !!!!!!
      end if
      s = is
      exit
    end if
  end do

  if (s == 0) then
    write (6, '(3A)') "[Warning] Keyword '@@@@' not found in '", filename, "'"
    stop 1
  end if

  return
end subroutine file_member_replace

end module common_nml
