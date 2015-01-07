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

  !--- PARAM_LETKF
  integer :: MEMBER = 3      ! ensemble size

  integer :: SLOT_START = 1
  integer :: SLOT_END = 1
  integer :: SLOT_BASE = 1
  real(r_size) :: SLOT_TINTERVAL = 3600.0d0

  real(r_size) :: SIGMA_OBS = 500.0d3
!  real(r_size) :: SIGMA_OBS_RAIN = 350.0d3
  real(r_size) :: SIGMA_OBSV = 0.4d0
!  real(r_size) :: SIGMA_OBSV_RAIN = 0.4d0
!  real(r_size) :: BASE_OBSV_RAIN = 85000.0d0
  real(r_size) :: SIGMA_OBST = 3.0d0
  real(r_size) :: VAR_LOCAL_UV(nvarmax)   = 1.0d0
  real(r_size) :: VAR_LOCAL_T(nvarmax)    = 1.0d0
  real(r_size) :: VAR_LOCAL_Q(nvarmax)    = 1.0d0
  real(r_size) :: VAR_LOCAL_PS(nvarmax)   = 1.0d0
  real(r_size) :: VAR_LOCAL_RAIN(nvarmax) = 1.0d0
  real(r_size) :: VAR_LOCAL_TC(nvarmax)   = 1.0d0

!RESHAPE( (/ &
!!       U    V    W    T    P    Q   QC   QR   QI   QS   QG
!   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! U,V
!   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! T,Tv
!   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! Q,RH
!   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! PS
!   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! RAIN
!   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0 & ! TC
!   & /),(/nv3d+nv2d,nvarlocal/))

  real(r_size) :: COV_INFL_MUL = 1.0d0   ! > 0: globally constant covariance inflation
                                          ! < 0: 3D inflation values input from a GPV file "infl_mul.grd"
  real(r_size) :: MIN_INFL_MUL = 0.0d0    ! minimum inlfation factor
  logical :: ADAPTIVE_INFL_INIT = .false.
  real(r_size) :: RELAX_ALPHA = 0.0d0     ! relaxation parameter
  real(r_size) :: SP_INFL_ADD = 0.0d0     ! additive inflation

  integer :: LEV_UPDATE_Q = 100000        ! q and qc are only updated below and equal to this model level
  real(r_size) :: Q_SPRD_MAX = 0.5        ! maximum q (ensemble spread)/(ensemble mean)

  !--- PARAM_LETKF_PRC
  integer :: NNODES = 1
  integer :: PPN = 1
  integer :: MEM_NODES = 1
  integer :: MEM_NP = 1
!  integer :: PRC_NUM_X_LETKF = 1
!  integer :: PRC_NUM_Y_LETKF = 1

  !--- PARAM_LETKF_OBS
  LOGICAL :: OMB_OUTPUT = .true.
  LOGICAL :: OMA_OUTPUT = .true.
  LOGICAL :: OBSGUES_OUTPUT = .false.
  LOGICAL :: OBSANAL_OUTPUT = .false.

contains
!-----------------------------------------------------------------------
! PARAM_LETKF
!-----------------------------------------------------------------------
subroutine read_nml_letkf
  use common_mpi, only: nprocs
  implicit none
  integer :: ierr
  
  namelist /PARAM_LETKF/ &
    MEMBER, &
    SLOT_START, &
    SLOT_END, &
    SLOT_BASE, &
    SLOT_TINTERVAL, &
    SIGMA_OBS, &
!    SIGMA_OBS_RAIN, &
    SIGMA_OBSV, &
!    SIGMA_OBSV_RAIN, &
!    BASE_OBSV_RAIN, &
    SIGMA_OBST, &
    VAR_LOCAL_UV, &
    VAR_LOCAL_T, &
    VAR_LOCAL_Q, &
    VAR_LOCAL_PS, &
    VAR_LOCAL_RAIN, &
    VAR_LOCAL_TC, &
    COV_INFL_MUL, &
    MIN_INFL_MUL, &
    ADAPTIVE_INFL_INIT, &
    RELAX_ALPHA, &
    SP_INFL_ADD, &
    LEV_UPDATE_Q, &
    Q_SPRD_MAX

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'xxx Not found namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist LETKF_PARAM_PRC. Check!'
    stop
  endif

  return
end subroutine read_nml_letkf

!-----------------------------------------------------------------------
! PARAM_LETKF_PRC
!-----------------------------------------------------------------------
subroutine read_nml_letkf_prc
  use common_mpi, only: nprocs
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
    write(6,*) 'xxx Not found namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist LETKF_PARAM_PRC. Check!'
    stop
  endif

  return
end subroutine read_nml_letkf_prc

!-----------------------------------------------------------------------
! PARAM_LETKF_OBS
!-----------------------------------------------------------------------
subroutine read_nml_letkf_obs
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_OBS/ &
    OMB_OUTPUT, &
    OMA_OUTPUT, &
    OBSGUES_OUTPUT, &
    OBSANAL_OUTPUT

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_OBS,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'xxx Not found namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist LETKF_PARAM_OBS. Check!'
    stop
  endif

  return
end subroutine read_nml_letkf_obs


end module common_nml
