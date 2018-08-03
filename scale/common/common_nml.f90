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
  integer, parameter :: nv3d = 11                 ! Number of 3D state variables (in SCALE restart files)
  integer, parameter :: nv2d = 0                  ! Number of 2D state variables (in SCALE restart files)
  integer, parameter :: nid_obs = 16              ! Number of variable types
  integer, parameter :: nobtype = 24              ! Number of observation report types
  integer, parameter :: nch = 10                  ! H08 Num of Himawari-8 (IR) channels

  integer, parameter :: nobsfilemax = 10
  integer, parameter :: filelenmax = 256

  integer, parameter :: memflen = 4               ! Length of a member string
  character(len=memflen), parameter :: memf_notation = '@@@@' ! Special notation to be replaced with a member string
  character(len=memflen), parameter :: memf_mean = 'mean' ! Member string for the ensemble mean
  character(len=memflen), parameter :: memf_mdet = 'mdet' ! Member string for the deterministic run
  character(len=memflen), parameter :: memf_sprd = 'sprd' ! Member string for the ensemble spread

  !--- &PARAM_ENSEMBLE
  ! Ensemble settings
  integer :: MEMBER = 3                           ! Ensemble size
  integer :: MEMBER_RUN = 1                       ! Actual number of members to be run
  integer :: MEMBER_ITER = 0                      ! Iteration number determining the range of ensemble members to be run
                                                  ! - 0:  Run all iteration for all members
  character(filelenmax) :: CONF_FILES = 'run.@@@@.conf' ! Filename pattern of the configuration files ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used)
  logical :: CONF_FILES_SEQNUM = .false.          ! [DEPRECATED]

  logical :: DET_RUN = .false.                    ! Enable the deterministic run (Schraff et al. 2016 QJRMS)?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  logical :: DET_RUN_CYCLED = .true.              ! Cycle the deterministic run?
                                                  ! - .false.:  No - Use the forecast from the analysis ensemble mean (in the previous cycle) as the first guess for the deterministic analysis
                                                  ! - .true.:  Yes - Use the forecast from the deterministic analysis (in the previous cycle) as the first guess for the deterministic analysis

  !--- &PARAM_MODEL
  ! Model-related settings
  character(len=10) :: MODEL = 'scale-rm'         ! Model name (always set to 'scale-rm')
  logical :: VERIFY_COORD = .false.               ! Verify the vertical coordinate settings with the vertical coordinate values in the input file?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes

!  !--- &PARAM_IO
!  ! IO settings
!  integer :: IO_AGGREGATE = .false.

  !--- &PARAM_OBSOPE
  ! Observation operator settings
  integer               :: OBS_IN_NUM = 1         ! Number of input observation data files
  character(filelenmax) :: OBS_IN_NAME(nobsfilemax) = 'obs.dat' ! Array of filenames of each observation file
  integer               :: OBS_IN_FORMAT(nobsfilemax) = 1 ! Array of data format of each observation file
                                                          ! - 1:  Conventional data in LETKF format
                                                          ! - 2:  Radar data in LETKF format
                                                          ! - 3:  Himawari-8 data in LETKF format
  logical               :: OBSDA_RUN(nobsfilemax) = .true. ! Array setting whether to run observation operator for each observation file
                                                           ! - .false.:  No
                                                           ! - .true.:  Yes
  logical               :: OBSDA_OUT = .false.    ! Output observation operator results [i.e., H(x^b)] to files?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes - Can be used as a separate observation operator program
  character(filelenmax) :: OBSDA_OUT_BASENAME = 'obsda.@@@@' ! Base filename pattern of observation operator result outputs ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used)
  character(filelenmax) :: OBSDA_MEAN_OUT_BASENAME = '' ! Base filename of observation operator result outputs for the ensemble mean
  character(filelenmax) :: OBSDA_MDET_OUT_BASENAME = '' ! Base filename of observation operator result outputs for the deterministic run

  character(filelenmax) :: HISTORY_IN_BASENAME = 'hist.@@@@' ! Base filename pattern of input history files for observation operator calculation ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used)
  character(filelenmax) :: HISTORY_MEAN_IN_BASENAME = '' ! Base filename of input history files for the ensemble mean for observation operator calculation
  character(filelenmax) :: HISTORY_MDET_IN_BASENAME = '' ! Base filename of input history files for the deterministic run for observation operator calculation

  integer               :: SLOT_START = 1         ! Start time-slot number for 4D-LETKF
  integer               :: SLOT_END = 1           ! End time-slot number for 4D-LETKF
  integer               :: SLOT_BASE = 1          ! Base time-slot number for 4D-LETKF
  real(r_size)          :: SLOT_TINTERVAL = 3600.0d0 ! Time-slot interval for 4D-LETKF

  !--- &PARAM_LETKF
  ! General LETKF settings
  logical               :: OBSDA_IN = .false.     ! Input observation operator results [i.e., H(x^b)] from a separate observation operator program
  character(filelenmax) :: OBSDA_IN_BASENAME = 'obsda.@@@@' ! Base filename pattern of input observation operator results from a separate observation operator program ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used)
  character(filelenmax) :: OBSDA_MEAN_IN_BASENAME = '' ! Base filename of input observation operator results for the ensemble mean
  character(filelenmax) :: OBSDA_MDET_IN_BASENAME = '' ! Base filename of input observation operator results for the deterministic run
  character(filelenmax) :: GUES_IN_BASENAME = 'gues.@@@@' ! Base filename pattern of input first-guess files ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used)
  character(filelenmax) :: GUES_MEAN_INOUT_BASENAME = '' ! Base filename of first-guess files for the ensemble mean (may be used as both input and output)
  character(filelenmax) :: GUES_MDET_IN_BASENAME = '' ! Base filename of first-guess files for the deterministic run
  logical               :: GUES_SPRD_OUT = .true. ! Output first-guess ensemble spread?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  character(filelenmax) :: GUES_SPRD_OUT_BASENAME = '' ! Base filename of output first-guess ensemble spread
  character(filelenmax) :: ANAL_OUT_BASENAME = 'anal.@@@@' ! Base filename pattern of output analysis files ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used)
  character(filelenmax) :: ANAL_MEAN_OUT_BASENAME = '' ! Base filename of output analysis files for the ensemble mean
  character(filelenmax) :: ANAL_MDET_OUT_BASENAME = '' ! Base filename of output analysis files for the deterministic run
  logical               :: ANAL_SPRD_OUT = .true. ! Output analysis ensemble spread?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  character(filelenmax) :: ANAL_SPRD_OUT_BASENAME = '' ! Base filename of output analysis spread
  character(filelenmax) :: LETKF_TOPO_IN_BASENAME = 'topo' ! Base filename of the input topographic file
                                                           !!!!!!! directly use the SCALE namelist ?? !!!!!!

  real(r_size) :: INFL_MUL = 1.0d0                ! Multiplicative covariance inflation parameter
                                                  ! - 1:  Disable the multiplicative inflation
                                                  ! - > 0:  Use a global constant inflation parameter
                                                  ! - <= 0:  Use a 3D inflation field from INFL_MUL_IN_BASENAME file
  real(r_size) :: INFL_MUL_MIN = -1.0d0           ! Minimum multiplicative inlfation parameter
                                                  ! - <= 0:  No minimum setting
  logical :: INFL_MUL_ADAPTIVE = .false.          ! Output adaptively estimated 3D inflation field to INFL_MUL_OUT_BASENAME file?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  character(filelenmax) :: INFL_MUL_IN_BASENAME = 'infl' ! Base filename of input 3D inflation field
  character(filelenmax) :: INFL_MUL_OUT_BASENAME = 'infl' ! Base filename of output (adaptively estimated) 3D inflation field

  real(r_size) :: INFL_ADD = 0.0d0                ! Additive covariance inflation parameter; this value will be multiplied to the input additive inflation field when using additive inflation
                                                  ! - < 0:  Disable the additive inflation
  character(filelenmax) :: INFL_ADD_IN_BASENAME = 'addi.@@@@' ! Base filename pattern of the input additive inflation field ([special notation](#special-notations-that-can-be-used-in-some-namelist-variables-for-filename-patterns) can be used)
  logical :: INFL_ADD_SHUFFLE = .false.           ! Shuffle the ensemble members for additive inflation field?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  logical :: INFL_ADD_Q_RATIO = .false.           ! For moisture field, further multiply the additive inflation field by the first-guess ensemble mean values (i.e., use the additive inflation field as ratio)?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  logical :: INFL_ADD_REF_ONLY = .false.          ! Apply the additive inflation only around where raining reflectivity (> RADAR_REF_THRES_DBZ) observations exist?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes

  real(r_size) :: RELAX_ALPHA = 0.0d0             ! Relaxation-to-prior-perturbation (RTPP) parameter (Zhang et al. 2004 MWR)
  real(r_size) :: RELAX_ALPHA_SPREAD = 0.0d0      ! Relaxation-to-prior-spread (RTPS) parameter (Whitaker and Hamill 2012 MWR)
  logical :: RELAX_TO_INFLATED_PRIOR = .false.    ! Choice of using covariance relaxation and multiplicative inflation together
                                                  ! - .true.:  Relaxation to the prior after the multiplicative inflation
                                                  ! - .false.:  Relaxation to the original prior before the multiplicative inflation
  logical :: RELAX_SPREAD_OUT = .false.           ! Output the equivalent multiplicative inflation field when using RTPS?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  character(filelenmax) :: RELAX_SPREAD_OUT_BASENAME = 'rtps' ! Base filename of equivalent multiplicative inflation field output (when using RTPS)

  real(r_size) :: GROSS_ERROR = 5.0d0             ! Threshold of gross error check (times of observation errors)
  real(r_size) :: GROSS_ERROR_RAIN = -1.0d0       ! Threshold of gross error check for precipitation data
                                                  ! - 0:  Same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_RADAR_REF = -1.0d0  ! Threshold of gross error check for radar reflectivity data
                                                  ! - 0:  Same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_RADAR_VR = -1.0d0   ! Threshold of gross error check for radar radial velocity data
                                                  ! - 0:  Same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_RADAR_PRH = -1.0d0  ! [NOT IMPLEMENTED] Threshold of gross error check for radar pseudo-RH data
                                                  ! - 0:  Same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_H08 = -1.0d0        ! Threshold of gross error check for Himawari-8 radiance data
                                                  ! - 0:  Same as GROSS_ERROR
  real(r_size) :: GROSS_ERROR_TCX = -1.0d0        ! [TENTATIVE]
  real(r_size) :: GROSS_ERROR_TCY = -1.0d0        ! [TENTATIVE]
  real(r_size) :: GROSS_ERROR_TCP = -1.0d0        ! [TENTATIVE]
  real(r_size) :: Q_UPDATE_TOP = 0.0d0            ! Pressure level (Pa) only below which water vapor and hydrometeors are updated
  real(r_size) :: Q_SPRD_MAX = -1.0D0             ! Maximum ratio of ensemble spread to ensemble mean for mositure in the analysis; if the analysis ensemble spread is greater than this ratio, scale the ensemble perturbations to reduce the spread to this ratio 
                                                  ! - <= 0:  Disabled

  real(r_size) :: BOUNDARY_BUFFER_WIDTH = 0.0d0   ! Width (m) of the buffer area along the lateral boundary where the analysis increment is gradually reduced to zero

  logical :: POSITIVE_DEFINITE_Q = .false.        ! Force setting the negative values in the analysis water vapor field to zero?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  logical :: POSITIVE_DEFINITE_QHYD = .false.     ! Force setting the negative values in the analysis hydrometeor fields to zero?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  real(r_size) :: TC_SEARCH_DIS = 200.0d3         ! [TENTATIVE] (m)

  real(r_size) :: PS_ADJUST_THRES = 100.d0        ! Threshold of elevation difference (m) between the station report and the model topography
                                                  ! Within the threshold surface pressure observations are assimilated (height adjustment will be performed to compensate this difference); beyond this threshold the surface pressure observations are discarded

  logical :: NOBS_OUT = .false.                   ! Output the field of actual observation numbers assimilated in each grid when the observation number limit (Hamrud et al. 2015 MWR) is enabled?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  character(filelenmax) :: NOBS_OUT_BASENAME = 'nobs' ! Base filename of the field of actual observation numbers assimilated in each grid when the observation number limit is enabled

  real(r_size) :: COV_INFL_MUL = 1.0d0            ! [FOR BACKWARD COMPATIBILITY] = INFL_MUL
  real(r_size) :: MIN_INFL_MUL = 0.0d0            ! [FOR BACKWARD COMPATIBILITY] = INFL_MUL_MIN
  logical :: ADAPTIVE_INFL_INIT = .false.         ! [FOR BACKWARD COMPATIBILITY] = INFL_MUL_ADAPTIVE
  real(r_size) :: BOUNDARY_TAPER_WIDTH = 0.0d0    ! [FOR BACKWARD COMPATIBILITY] = BOUNDARY_BUFFER_WIDTH

  !--- &PARAM_LETKF_PRC
  ! Parallelization settings for LETKF
  integer :: NNODES = 1                           ! Total number of nodes used for the LETKF 
  integer :: PPN = 1                              ! Number of MPI processes used per nodes for the LETKF
  integer :: MEM_NODES = 1                        ! Number of nodes used for one ensemble member in the LETKF
  integer :: MEM_NP = 1                           ! Number of MPI processes used for one ensemble member in the LETKF

!  integer :: PRC_NUM_X_LETKF = 1
!  integer :: PRC_NUM_Y_LETKF = 1

  !--- &PARAM_LETKF_OBS
  ! Observation-specific settings
  logical :: USE_OBS(nobtype) = .true.            ! Array setting whether each observation report type is used?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes

  real(r_size) :: HORI_LOCAL(nobtype) = &
    (/500.0d3, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
       -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
       -1.0d0, -1.0d0, -1.0d0, -1.0d0/)           ! Array of horizontal covariance localization length scale for each observation report type
                                                  ! - > 0:  Horizontal localization length scale (m)
                                                  ! - 0:  [NOT IMPLEMENTED] No horizontal localization
                                                  ! - < 0:  Same setting as HORI_LOCAL(1)

  real(r_size) :: VERT_LOCAL(nobtype) = &
    (/ 0.4d0,   -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0,   -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, 1000.0d0, -1.0d0, -1.0d0/)          ! Array of vertical covariance localization length scale for each observation report type
                                                  ! - > 0:  Vertical localization length scale [ln(p) or m depending on the report type]
                                                  ! - 0:  No vertical localization
                                                  ! - < 0:  Same setting as VERT_LOCAL(1)

  real(r_size) :: TIME_LOCAL(nobtype) = &
    (/ 0.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0/)            ! [NOT IMPLEMENTED] Array of temporal covariabce localization interval for each observation report type
                                                  ! - > 0:  [NOT IMPLEMENTED] Temporal localization interval (sec)
                                                  ! - 0:  No temporal localization
                                                  ! - < 0:  Same setting as TIME_LOCAL(1)

  real(r_size) :: HORI_LOCAL_RADAR_OBSNOREF = -1.0d0 ! Horizontal covariance localization length scale (m) for clear-sky radar reflectivity data (<= RADAR_REF_THRES_DBZ)
                                                     ! - < 0:  Same setting as HORI_LOCAL(22) for all radar data
  real(r_size) :: HORI_LOCAL_RADAR_VR = -1.0d0    ! Horizontal covariance localization length scale (m) for radar radial velocity data
                                                  ! - < 0:  Same setting as HORI_LOCAL(22) for all radar data
  real(r_size) :: VERT_LOCAL_RADAR_VR = -1.0d0    ! Vertical covariance localization length scale (m) for radar radial velocity data
                                                  ! - < 0:  Same setting as HORI_LOCAL(22) for all radar data
  real(r_size) :: VERT_LOCAL_RAIN_BASE = 85000.0d0 ! [NOT IMPLEMENTED] Base level (Pa) for the vertical covariance localization for rain data

  integer :: MAX_NOBS_PER_GRID(nobtype) = &
    (/ 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
      -1, -1, -1, -1/)                            ! Observation number limit: maximum number of observations for each observation report type and each variable assimilated in a grid (Hamrud et al. 2015 MWR)
                                                  ! - > 0:  Enable the observation number limit
                                                  ! - 0:  No observation number limit
                                                  ! - < 0:  Same setting as MAX_NOBS_PER_GRID(1)

  integer :: MAX_NOBS_PER_GRID_CRITERION = 1      ! Criterion to select observations for the observation number limit
                                                  ! - 1:  Normalized 3D distance (from closest)
                                                  ! - 2:  Covariance localization weight (from largest)
                                                  ! - 3:  Weighted observation error variance (from smallest)

  real(r_size) :: OBS_MIN_SPACING(nobtype) = &
    (/300.0d3, 100.0d3, 100.0d3, 150.0d3, 300.0d3, 150.0d3, 150.0d3, 100.0d3, 150.0d3, 150.0d3, &
      150.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 150.0d3, 300.0d3, 150.0d3, 150.0d3, 150.0d3, &
      150.0d3,   1.0d3,  15.0d3,1000.0d3/)        ! Array of estimates of a typical minimum horizontal observation spacing (m) (in the densest observed area) for each obsetvation report type.
                                                  ! * This setting only affects the computational speed but not the analysis results
                                                  ! * This setting is used to automatically determine OBS_SORT_GRID_SPACING, effective only when OBS_SORT_GRID_SPACING = 0
                                                  ! - <= 0:  Same setting as OBS_MIN_SPACING(1)

  real(r_size) :: OBS_SORT_GRID_SPACING(nobtype) = &
    (/ 0.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
      -1.0d0, -1.0d0, -1.0d0, -1.0d0/)            ! Array of optimal grid spacing (m) for observation bucket sorting for each observation report type
                                                  ! * This setting only affects the computational speed but not the analysis results
                                                  ! - 0:  Automatically determined based on HORI_LOCAL, MAX_NOBS_PER_GRID, and OBS_MIN_SPACING
                                                  ! - < 0:  Same setting as OBS_SORT_GRID_SPACING(1)

  !--- &PARAM_LETKF_VAR_LOCAL
  ! Variable localization settings
  real(r_size) :: VAR_LOCAL_UV(nv3d+nv2d)        = 1.0d0 ! Array of variable covariance localization factors between u- and v-wind observations and each state variable
                                                         ! - 1:  No variable localization
  real(r_size) :: VAR_LOCAL_T(nv3d+nv2d)         = 1.0d0 ! Same as above, but for temperature observations
  real(r_size) :: VAR_LOCAL_Q(nv3d+nv2d)         = 1.0d0 ! Same as above, but for water vapor observations
  real(r_size) :: VAR_LOCAL_PS(nv3d+nv2d)        = 1.0d0 ! Same as above, but for surface pressure observations
  real(r_size) :: VAR_LOCAL_RAIN(nv3d+nv2d)      = 1.0d0 ! [NOT IMPLEMENTED] Same as above, but for rain observations
  real(r_size) :: VAR_LOCAL_TC(nv3d+nv2d)        = 1.0d0 ! [NOT IMPLEMENTED] Same as above, but for TC vital observations
  real(r_size) :: VAR_LOCAL_RADAR_REF(nv3d+nv2d) = 1.0d0 ! Same as above, but for radar reflectivity observations
  real(r_size) :: VAR_LOCAL_RADAR_VR(nv3d+nv2d)  = 1.0d0 ! Same as above, but for radar radial velocity observations
  real(r_size) :: VAR_LOCAL_H08(nv3d+nv2d)       = 1.0d0 ! [TENTATIVE] Same as above, but for Himawari-8 observations

  !--- &PARAM_LETKF_MONITOR
  ! Observation diagnostic settings
  logical :: DEPARTURE_STAT = .true.              ! Output observation departure statistics (O-B and O-A) for conventional observations?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  logical :: DEPARTURE_STAT_RADAR = .false.       ! Output observation departure statistics (O-B and O-A) for radar observations?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  logical :: DEPARTURE_STAT_H08 = .false.         ! Output observation departure statistics (O-B and O-A) for Himawari-8 observations?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  real(r_size) :: DEPARTURE_STAT_T_RANGE = 0.0d0  ! Range of time difference to the analysis time, within which observations are considered in the departure statistics
                                                  ! - 0:  No time range restriction
  logical :: DEPARTURE_STAT_ALL_PROCESSES = .true. ! Print the departure statistics by all processes?
                                                   ! - .false.:  No - The statistics are only printed by the ensemble mean group, which may save computational time
                                                   ! - .true.:  Yes - The same statistics are printed by all processes

  LOGICAL               :: OBSDEP_OUT = .true.    ! Output observation departure (innovation) data for all observations into a binary file?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  character(filelenmax) :: OBSDEP_OUT_BASENAME = 'obsdep' ! Filename of observation departure (innovation) data output
  LOGICAL               :: OBSGUES_OUT = .false.                 ! [NOT IMPLEMENTED]
  character(filelenmax) :: OBSGUES_OUT_BASENAME = 'obsgues.@@@@' ! [NOT IMPLEMENTED]
  LOGICAL               :: OBSANAL_OUT = .false.                 ! [NOT IMPLEMENTED]
  character(filelenmax) :: OBSANAL_OUT_BASENAME = 'obsanal.@@@@' ! [NOT IMPLEMENTED]

  !--- &PARAM_LETKF_RADAR
  ! Settings for radar data
  logical :: USE_RADAR_REF       = .true.         ! Assimilate radar reflectivity observations?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  logical :: USE_RADAR_VR        = .true.         ! Assimilate radar radial velocity observations?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes
  logical :: USE_RADAR_PSEUDO_RH = .false.        ! [NOT IMPLEMENTED] Assimilate pseudo-relative himidity observations?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes

  logical :: USE_OBSERR_RADAR_REF = .false.       ! Use OBSERR_RADAR_REF for the observation error of radar reflectivity observations instead of that provided in the observation files
  logical :: USE_OBSERR_RADAR_VR = .false.        ! Use OBSERR_RADAR_VR for the observation error of radar radial velocity observations instead of that provided in the observation files

  logical :: RADAR_OBS_4D = .false.               ! Radar observation data file is in a new format with the "time difference to the analysis time" column, allowing for 4D LETKF?
                                                  ! - .false.:  No - Old-format file (without the "time difference to the analysis time" column)
                                                  ! - .true.:  Yes - New-format file (with the "time difference to the analysis time" column)

  REAL(r_size) :: RADAR_REF_THRES_DBZ = 15.0d0    ! Threshold of raining and clear-sky radar reflectivity observations (dBZ)
  INTEGER :: MIN_RADAR_REF_MEMBER = 1             ! Threshold of number of first-guess ensemble members with raining reflectivity to assimilate the *__clear-sky__* radar reflectivity data
                                                  ! * The observation data are assimilated only when the number of raining (> RADAR_REF_THRES_DBZ) members is above this threshold. 
  INTEGER :: MIN_RADAR_REF_MEMBER_OBSREF = 1      ! Same as above, but the threshold to assimilate the *__raining__* reflectivity observations

  REAL(r_size) :: MIN_RADAR_REF_DBZ = 0.0d0       ! Minimum useful radar reflectivity value (dBZ); all reflectivity data below this value are re-assigned to a constant depending on LOW_REF_SHIFT
  REAL(r_size) :: LOW_REF_SHIFT = 0.0d0           ! Shift of the constant refelectivity value for those data smaller than MIN_RADAR_REF_DBZ (dBZ); all reflectivity data below MIN_RADAR_REF_DBZ are set to (MIN_RADAR_REF_DBZ + LOW_REF_SHIFT)
                                                  ! * This setting should be zero or negative

  real(r_size) :: RADAR_ZMAX = 99.0d3             ! Maximum height level (m) of radar data to be assimilated

  REAL(r_size) :: RADAR_PRH_ERROR = 0.1d0         ! [NOT IMPLEMENTED]

  !These 2 flags affects the computation of model reflectivity and radial velocity. 
  INTEGER :: INTERPOLATION_TECHNIQUE = 1          ! [TENTATIVE]
  INTEGER :: METHOD_REF_CALC = 3                  ! Method to compute the radar reflectivity in the radar observation operator

  LOGICAL :: USE_TERMINAL_VELOCITY = .false.      ! Consider the terminal velocity of the hydrometeors in the radar observation operator?
                                                  ! - .false.:  No
                                                  ! - .true.:  Yes

  INTEGER :: NRADARTYPE = 1                       ! [NOT IMPLEMENTED]

  !--- &PARAM_LETKF_H08
  ! Settings for Himawari-8 data
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

  !--- &PARAM_OBS_ERROR
  ! Observation error settings
  ! __* Ususally the observation errors provided in the observation data files, instead of the values here, are used in the assimilation, unless some special options are enabled__
  real(r_size) :: OBSERR_U = 1.0d0                ! Observaiton error of u-wind observations (m/s)
  real(r_size) :: OBSERR_V = 1.0d0                ! Observation error of v-wind observations (m/s)
  real(r_size) :: OBSERR_T = 1.0d0                ! Observation error of temperature observations (K)
  real(r_size) :: OBSERR_Q = 0.001d0              ! Observation error of water vapor observations (kg/kg)
  real(r_size) :: OBSERR_RH = 10.0d0              ! Observation error of relative humidity observations (%)
  real(r_size) :: OBSERR_PS = 100.0d0             ! Observation error of surface pressure observations (Pa)
  real(r_size) :: OBSERR_RADAR_REF = 5.0d0        ! Observation error of radar reflectivity observations (dBZ)
  real(r_size) :: OBSERR_RADAR_VR = 3.0d0         ! Observation error of radar radial velocity observations (m/s)
  real(r_size) :: OBSERR_TCX = 50.0d3             ! [TENTATIVE] Observation error of TC x-center-position observations (m)
  real(r_size) :: OBSERR_TCY = 50.0d3             ! [TENTATIVE] Observation error of TC y-center-position observations (m)
  real(r_size) :: OBSERR_TCP = 5.0d2              ! [TENTATIVE] Observation error of TC minimum sea level pressure observations (Pa)
  real(r_size) :: OBSERR_H08(nch) = (/5.0d0,5.0d0,5.0d0,5.0d0,5.0d0,&
                                      5.0d0,5.0d0,5.0d0,5.0d0,5.0d0/) ! [TENTATIVE] Observation error of Himawari-8 observations

  !--- &PARAM_OBSSIM
  ! Settings for the "observation simulator" (obssim) program
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

  !--- &

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
    write(6,*) 'Error: /PARAM_ENSEMBLE/ is not found in namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_ENSEMBLE. Check!'
    stop
  endif

  write(6, nml=PARAM_ENSEMBLE)

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
    write(6,*) 'Warning: /PARAM_MODEL/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_MODEL. Check!'
    stop
  endif

  write(6, nml=PARAM_MODEL)

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
!    write(6,*) 'Warning: /PARAM_IO/ is not found in namelist.'
!!    stop
!  elseif (ierr > 0) then !--- fatal error
!    write(6,*) 'xxx Not appropriate names in namelist PARAM_IO. Check!'
!    stop
!  endif

!  write(6, nml=PARAM_IO)

!  return
!end subroutine read_nml_io

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
    write(6,*) 'Error: /PARAM_OBSOPE/ is not found in namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_OBSOPE. Check!'
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

  write(6, nml=PARAM_OBSOPE)

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

  write(6, nml=PARAM_LETKF)

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
    write(6,*) 'Warning: /PARAM_LETKF_PRC/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_PRC. Check!'
    stop
  endif

  write(6, nml=PARAM_LETKF_PRC)

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
    write(6,*) 'Warning: /PARAM_LETKF_VAR_LOCAL/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_VAR_LOCAL. Check!'
    stop
  endif

  write(6, nml=PARAM_LETKF_VAR_LOCAL)

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
    write(6,*) 'Warning: /PARAM_LETKF_MONITOR/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_MONITOR. Check!'
    stop
  endif

  write(6, nml=PARAM_LETKF_MONITOR)

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
    write(6,*) 'Warning: /PARAM_LETKF_H08/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_H08. Check!'
    stop
  endif

  write(6, nml=PARAM_LETKF_H08)

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
    write(6,*) 'Warning: /PARAM_OBS_ERROR/ is not found in namelist.'
!    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_OBS_ERROR. Check!'
    stop
  endif

  write(6, nml=PARAM_OBS_ERROR)

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

  write(6, nml=PARAM_OBSSIM)

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
