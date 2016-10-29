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
  integer, parameter :: nobsfilemax = 10
  integer, parameter :: filelenmax = 256
  integer, parameter :: memberflen = 4 ! Length of member # in filename
  integer, parameter :: nch = 10 ! Num of Himawari-8 (IR) channels 
  integer,parameter  :: PNUM_TOMITA = 6 ! Number of the Tomita (2008)'s parameters.
                                        ! [Cr, Cs, drag_g, beta_saut, gamma_saut,
                                        ! gamma_sacr]

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
  logical               :: OBSDA_IN = .false.
  character(filelenmax) :: OBSDA_IN_BASENAME = 'obsda.@@@@'
  character(filelenmax) :: GUES_IN_BASENAME = 'gues.@@@@'
  character(filelenmax) :: GUES_OUT_MEAN_BASENAME = 'gues.mean'
  character(filelenmax) :: GUES_OUT_SPRD_BASENAME = 'gues.sprd'
  character(filelenmax) :: ANAL_OUT_BASENAME = 'anal.@@@@'
  character(filelenmax) :: ANAL_OUT_MEAN_BASENAME = 'anal.mean'
  character(filelenmax) :: ANAL_OUT_SPRD_BASENAME = 'anal.sprd'
  character(filelenmax) :: LETKF_TOPO_IN_BASENAME = 'topo'  !!!!!! -- directly use the SCALE namelist --???? !!!!!!

  real(r_size) :: SIGMA_OBS = 500.0d3
  real(r_size) :: SIGMA_OBS_RAIN = -1.0d0  ! < 0: same as SIGMA_OBS
  real(r_size) :: SIGMA_OBS_RADAR = -1.0d0 ! < 0: same as SIGMA_OBS
  real(r_size) :: SIGMA_OBS_RADAR_OBSNOREF = -1.0d0 ! < 0: same as SIGMA_OBS_RADAR
  real(r_size) :: SIGMA_OBS_H08 = -1.0d0 ! < 0: same as SIGMA_OBS ! H08
  real(r_size) :: SIGMA_OBS_PEST = 1000.0d33 ! < 0: same as SIGMA_OBS 
  real(r_size) :: SIGMA_OBS_TC = -1.0d0 ! < 0: same as SIGMA_OBS 
  real(r_size) :: SIGMA_OBSV = 0.4d0
  real(r_size) :: SIGMA_OBSV_RAIN = -1.0d0 ! < 0: same as SIGMA_OBSV
  real(r_size) :: SIGMA_OBSZ_RADAR = 1000.0d0
  real(r_size) :: SIGMA_OBSV_H08 = -1.0d0 ! < 0: same as SIGMA_OBSV ! H08
  real(r_size) :: SIGMA_OBSV_TC = -1.0d0 ! < 0: same as SIGMA_OBSV 
  real(r_size) :: SIGMA_OBST = 3.0d0
  real(r_size) :: BASE_OBSV_RAIN = 85000.0d0

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

  integer :: MAX_NOBS_PER_GRID = 0   ! observation number limit; <= 0: Do not use
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
  real(r_size) :: OBSERR_TCX = 50.0d3 ! (m)
  real(r_size) :: OBSERR_TCY = 50.0d3 ! (m)
  real(r_size) :: OBSERR_TCP = 5.0d2 ! (Pa)
  logical :: USE_OBSERR_RADAR_REF = .false.
  logical :: USE_OBSERR_RADAR_VR = .false.
!
! 
  real(r_size) :: OBSERR_H08(nch) = (/5.0d0,5.0d0,5.0d0,5.0d0,5.0d0,&
                                      5.0d0,5.0d0,5.0d0,5.0d0,5.0d0/) ! H08

  real(r_size) :: OBSERR_H08_MAX = 15.0d0
  real(r_size) :: OBSERR_H08_MIN = 0.5d0

  !--- PARAM_LETKF_MONITOR
  logical :: DEPARTURE_STAT = .true.
  logical :: DEPARTURE_STAT_RADAR = .false.
  logical :: DEPARTURE_STAT_H08 = .false.
  logical :: DEPARTURE_STAT_H08_ALL = .false.
  real(r_size) :: DEPARTURE_STAT_T_RANGE = 0.0d0 ! time range within which observations are considered in the departure statistics.
                                                 ! 0.0d0: no limit

  LOGICAL :: OMB_OUTPUT = .true.
  LOGICAL :: OMA_OUTPUT = .true.
  LOGICAL :: OBSGUES_OUTPUT = .false.
  LOGICAL :: OBSANAL_OUTPUT = .false.

  !--- PARAM_LETKF_RADAR
  logical :: USE_RADAR_REF       = .true.
  logical :: USE_RADAR_VR        = .true.
  logical :: USE_RADAR_PSEUDO_RH = .false.

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
  logical :: H08_OB_OBSERR = .false. ! Obs err for Him8 is assigned by using abs(O-B).
  logical :: H08_CLD_OBSERR = .false. ! Cloud dependent obs error for Him8. If this is true, obs error depending on CA is assigned in letkf
  real(r_size) :: H08_CLD_OBSERR_WTH = 1.0d0 ! Bin width of CA for cloud dependent obs error.
  integer :: H08_CLD_OBSERR_NBIN = 51 ! Number of bins for CA.
  integer :: H08_CLD_OBSERR_MTIME = 18 ! Max number of analysis time that is used to diagnose cloud dependent obserr  

  real(r_size) :: H08_CLD_OBSERR_GROSS_ERR = 20.0d0
  integer :: H08_CLD_OBSERR_MIN_SAMPLE = 500
  logical :: H08_DEBIAS_AMEAN = .false.
  logical :: H08_DEBIAS_CA = .false.
  integer :: H08_DEBIAS_CA_MIN_SAMPLE = 500
  logical :: H08_DEBIAS_CA_CLR = .false.
  logical :: H08_CLD_OBSERR_OB2 = .false. ! Maximum threshold for Him8 obs err using sqrt([O-B]**2)
  logical :: H08_CLD_OBSERR_BSPRD2 = .false. ! Output background spread**2 in obs space as a function of CA ! ?? not used
  logical :: H08_RTTOV_EXTRA_US76 = .false.

  integer :: H08_CH_USE(nch) = (/0,0,1,0,0,0,0,0,0,0/)
                        !! ch = (1,2,3,4,5,6,7,8,9,10)
                        !! (B07,B08,B09,B10,B11,B12,B13,B14,B15,B16)
                        !! ==1: Assimilate
                        !! ==0: NOT assimilate (rejected by QC in trans_XtoY_H08)
                        !! It is better to reject B11(ch=5) & B12(ch=6) obs because these bands are 
                        !! sensitive to chemicals.

#ifdef PEST_TOMITA
  !---PARAM_LETKF_PEST_TOMITA 
  !-- Parameter estimation for Tomita (2008)'s namelist parameters. 
  !-- !!
  !-- !! Parameters are treated as global constants.
  !-- !! If you wish to estimate as a local parameter, you should set
  !-- !! PEST_TOMITA_LOCAL2D = T
  !-- !!
  !-- !! Following Kotsuki et al. (20xx), parameters are converted into tanh space.
  !-- !! If you wish to treat in physical space, you should set
  !-- !! PEST_TOMITA_ARCTAN = F
  !
  integer  :: EPNUM_TOMITA = 0 ! Number of  Tomita (2008) parameters to be estimated.
                               ! Depending on  PEST_TOMITA_FLAG.
  CHARACTER(10) :: PNAME_TOMITA(PNUM_TOMITA) = (/'Cr        ', 'Cs        ', &
                                                 'drag_g    ', 'beta_saut ', &
                                                 'gamma_saut', 'gamma_sacr'/)
  LOGICAL :: PEST_TOMITA_FLAG(PNUM_TOMITA) = (/.false.,.false.,.false.,.false.,.false.,.false./)
  REAL(r_size) :: DEF_PARAM_TOMITA(PNUM_TOMITA) = (/130.D0, 4.84D0, &
                                                    0.6D0,  6.0D-3, & 
                                                   60.0D-3, 25.0D-3/) ! Default parameters in Tomita (2008).
  REAL(r_size) :: EPARAM_TOMITA_CR_LIMIT(2) = (/10.0D0,180.0D0/) ! Max/min threshold for Tomita (2008) parameter (Cr)
  REAL(r_size) :: EPARAM_TOMITA_CS_LIMIT(2) = (/0.1D0,8.0D0/) ! Max/min threshold for Tomita (2008) parameter (Cs)
  REAL(r_size) :: EPARAM_TOMITA_DRAGG_LIMIT(2) = (/0.1D0,5.0D0/) ! Max/min threshold for Tomita (2008) parameter (drag_g)
  REAL(r_size) :: EPARAM_TOMITA_BETAS_LIMIT(2) = (/0.1D-3,10.0D-3/) ! Max/min threshold for Tomita (2008) parameter (beta_saut)
  REAL(r_size) :: EPARAM_TOMITA_GAMMAT_LIMIT(2) = (/1.0D-3,7.0D-2/) ! Max/min threshold for Tomita (2008) parameter (gamma_saut)
  REAL(r_size) :: EPARAM_TOMITA_GAMMAR_LIMIT(2) = (/0.1D-3,7.0D-2/) ! Max/min threshold for Tomita (2008) parameter (gamma_sacr)
  REAL(r_size),ALLOCATABLE :: EPARAM_TOMITA_LIMIT(:,:) ! Max/min threshold for Tomita (2008) parameters.
  REAL(r_size) :: PEST_RELAX_ALPHA_SPREAD = 0.0d0
  LOGICAL :: PEST_TOMITA_ARCTAN = .true.
  LOGICAL :: PEST_TOMITA_LOCAL2D = .false.
#endif


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
    OBSDA_IN, &
    OBSDA_IN_BASENAME, &
    GUES_IN_BASENAME, &
    GUES_OUT_MEAN_BASENAME, &
    GUES_OUT_SPRD_BASENAME, &
    ANAL_OUT_BASENAME, &
    ANAL_OUT_MEAN_BASENAME, &
    ANAL_OUT_SPRD_BASENAME, &
    LETKF_TOPO_IN_BASENAME, &
    SIGMA_OBS, &
    SIGMA_OBS_RAIN, &
    SIGMA_OBS_RADAR, &
    SIGMA_OBS_RADAR_OBSNOREF, &
    SIGMA_OBS_H08, & ! H08
    SIGMA_OBS_TC, & 
    SIGMA_OBSV, &
    SIGMA_OBSV_RAIN, &
    SIGMA_OBSV_H08, & ! H08
    SIGMA_OBSV_TC, & 
    SIGMA_OBSZ_RADAR, &
    SIGMA_OBST, &
    BASE_OBSV_RAIN, &
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
    Q_UPDATE_TOP, &
    Q_SPRD_MAX, &
    BOUNDARY_BUFFER_WIDTH, &
    POSITIVE_DEFINITE_Q, &
    POSITIVE_DEFINITE_QHYD, &
    TC_SEARCH_DIS, &
    PS_ADJUST_THRES, &
    MAX_NOBS_PER_GRID, &
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
  if (SIGMA_OBS_RAIN < 0.0d0) then
    SIGMA_OBS_RAIN = SIGMA_OBS
  end if
  if (SIGMA_OBS_RADAR < 0.0d0) then
    SIGMA_OBS_RADAR = SIGMA_OBS
  end if
  if (SIGMA_OBS_RADAR_OBSNOREF < 0.0d0) then
    SIGMA_OBS_RADAR_OBSNOREF = SIGMA_OBS_RADAR
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
  if (SIGMA_OBS_TC < 0.0d0) then 
    SIGMA_OBS_TC = SIGMA_OBS
  end if
  if (SIGMA_OBSV_TC < 0.0d0) then 
    SIGMA_OBSV_TC = SIGMA_OBSV
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
    OBSERR_TCX, &
    OBSERR_TCY, &
    OBSERR_TCP, &
    OBSERR_H08, & ! H08
    OBSERR_H08_MAX, & ! H08
    OBSERR_H08_MIN, & ! H08
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
    DEPARTURE_STAT_H08_ALL, &
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

!-----------------------------------------------------------------------
! PARAM_LETKF_RADAR
!-----------------------------------------------------------------------
subroutine read_nml_letkf_h08
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_H08/   &
    H08_REJECT_LAND,           &
    H08_RTTOV_CLD,             &
    H08_RTTOV_MINQ,            &
    H08_RTTOV_CFRAC_CNST,      &
    H08_LIMIT_LEV,             &
    H08_BT_MIN,                &
    H08_RTTOV_EXTRA_US76,      &
    H08_OB_OBSERR,             &
    H08_CLD_OBSERR,            &
    H08_CLD_OBSERR_WTH,        &
    H08_CLD_OBSERR_NBIN,       &
    H08_CLD_OBSERR_GROSS_ERR,  &
    H08_CLD_OBSERR_MIN_SAMPLE, &
    H08_CLD_OBSERR_MTIME,      &
    H08_DEBIAS_AMEAN,          &
    H08_DEBIAS_CA,             &
    H08_DEBIAS_CA_CLR,         &
    H08_DEBIAS_CA_MIN_SAMPLE,  &
    H08_CLD_OBSERR_OB2,        &
    H08_CLD_OBSERR_BSPRD2,     &
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

#ifdef PEST_TOMITA
!-----------------------------------------------------------------------
! PARAM_LETKF_PEST_TOMITA
!-----------------------------------------------------------------------
subroutine read_nml_letkf_pest_tomita
  implicit none
  integer :: ierr
!  integer :: pr1

  namelist /PARAM_LETKF_PEST_TOMITA/ &
    PEST_TOMITA_FLAG,  &
    EPNUM_TOMITA,      &
    DEF_PARAM_TOMITA,  &
    EPARAM_TOMITA_CR_LIMIT, &
    EPARAM_TOMITA_CS_LIMIT, &
    EPARAM_TOMITA_DRAGG_LIMIT, &
    EPARAM_TOMITA_BETAS_LIMIT, &
    EPARAM_TOMITA_GAMMAT_LIMIT, &
    EPARAM_TOMITA_GAMMAR_LIMIT, &
    PEST_RELAX_ALPHA_SPREAD,&
    PEST_TOMITA_ARCTAN, &
    PEST_TOMITA_LOCAL2D

  if(.not.allocated(EPARAM_TOMITA_LIMIT)) allocate(EPARAM_TOMITA_LIMIT(2,1:PNUM_TOMITA))

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_PEST_TOMITA,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'Warning: /PARAM_LETKF_PEST_TOMITA/ is not found in namelist.' 
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist PARAM_LETKF_PEST_TOMITA. Check!'
    stop
  endif

  EPARAM_TOMITA_LIMIT(:,1) = EPARAM_TOMITA_CR_LIMIT(1:2)
  EPARAM_TOMITA_LIMIT(:,2) = EPARAM_TOMITA_CS_LIMIT(1:2)
  EPARAM_TOMITA_LIMIT(:,3) = EPARAM_TOMITA_DRAGG_LIMIT(1:2)
  EPARAM_TOMITA_LIMIT(:,4) = EPARAM_TOMITA_BETAS_LIMIT(1:2)
  EPARAM_TOMITA_LIMIT(:,5) = EPARAM_TOMITA_GAMMAT_LIMIT(1:2)
  EPARAM_TOMITA_LIMIT(:,6) = EPARAM_TOMITA_GAMMAR_LIMIT(1:2)

  EPNUM_TOMITA = COUNT(PEST_TOMITA_FLAG)

  write(6, nml=PARAM_LETKF_PEST_TOMITA)

  return
end subroutine read_nml_letkf_pest_tomita
#endif
!-----------------------------------------------------------------------

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
