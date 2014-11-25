module common_scalelib
!=======================================================================
!
! [PURPOSE:] Use the SCALE library
!
! [HISTORY:]
!   Novermber 2014  Guo-Yuan Lien  created
!
!=======================================================================
!$USE OMP_LIB
  use common
  use common_mpi
  use common_scale
  use common_mpi_scale, only: &
    nitmax, &
    proc2mem

!  use common_letkf, only: nbv

  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use dc_log, only: &
    loginit
!  use gtool_file, only: &
!!     fileread, &
!     filecloseall
!  use scale_grid_index, only: &
!    KHALO, IHALO, JHALO

  implicit none
  public

!  integer,parameter :: mpibufsize=1000000
!  integer,save :: nij1
!  integer,save :: nij1max
!  integer,allocatable,save :: nij1node(:)
!!  real(r_size),allocatable,save :: phi1(:)
!  real(r_size),allocatable,save :: lon1(:),lat1(:)
!  real(r_size),allocatable,save :: lonu1(:),latu1(:)
!  real(r_size),allocatable,save :: lonv1(:),latv1(:)
!  real(r_size),allocatable,save :: ri1(:),rj1(:)
!!  real(r_size),allocatable,save :: wg1(:)

!  integer,save :: nitmax ! maximum number of model files processed by a process
!  integer,allocatable,save :: procs(:)
!  integer,allocatable,save :: mem2node(:,:)
!  integer,allocatable,save :: mem2proc(:,:)
!  integer,allocatable,save :: proc2mem(:,:,:)
!  integer,save :: n_mem
!  integer,save :: n_mempn

contains
!-----------------------------------------------------------------------
! Start using SCALE library
!-----------------------------------------------------------------------
!subroutine set_scalelib(mem_np)
!subroutine set_scalelib(nitmax, proc2mem)
subroutine set_scalelib
  use common_nml, only: &
    MEM_NP

  use gtool_history, only: &
    historyinit
  use scale_process, only: &
    PRC_setup,    &
    PRC_MPIstart, &
!      PRC_mpifinish, &
    PRC_master, &
    PRC_myrank, &
    PRC_myrank_world, &
    PRC_2Drank, &
    PRC_NUM_X, &
    PRC_NUM_Y
!    prc_nu, &
  use scale_const, only: &
    CONST_setup
  use scale_calendar, only: &
    CALENDAR_setup
  use scale_random, only: &
    RANDOM_setup
  use scale_time, only: &
    TIME_setup
  use scale_grid, only: &
    GRID_setup, &
    GRID_DOMAIN_CENTER_X, &
    GRID_DOMAIN_CENTER_Y
!  use scale_grid_nest, only: &
!    NEST_setup
!  use scale_land_grid_index, only: &
!    LAND_GRID_INDEX_setup
!  use scale_land_grid, only: &
!    LAND_GRID_setup
!  use scale_urban_grid_index, only: &
!    URBAN_GRID_INDEX_setup
!  use scale_urban_grid, only: &
!    URBAN_GRID_setup
!  use scale_tracer, only: &
!    TRACER_setup
  use scale_fileio, only: &
     FILEIO_setup
  use scale_comm, only: &
    COMM_setup
!  use scale_topography, only: &
!    TOPO_setup
!  use scale_landuse, only: &
!    LANDUSE_setup
!  use scale_grid_real, only: &
!    REAL_setup
  use scale_mapproj, only: &
    MPRJ_setup
  implicit none

!  integer,intent(in) :: nitmax ! maximum number of model files processed by a process
!  integer,intent(in) :: proc2mem(2,nitmax,nprocs)
!  integer,intent(in) :: proc2mem(:,:,:)

!    integer,intent(in) :: mem_np

!    character(len=H_MID) :: DATATYPE = 'DEFAULT' !< REAL4 or REAL8
  integer :: rankidx(2)

  !-----------------------------------------------------------------------------

  ! start SCALE MPI
  call PRC_MPIstart(MEM_NP, nitmax, nprocs, proc2mem)
!  call PRC_MPIstart(nbv, MEM_NP, nitmax, nprocs, proc2mem)

  ! setup process
  call PRC_setup

  ! setup Log
  call LogInit(IO_FID_CONF, IO_FID_LOG, IO_L)

  ! setup constants
  call CONST_setup

  ! setup time
  call TIME_setup( setup_TimeIntegration = .false. )

  call PROF_rapstart('Initialize')

  ! setup horizontal/vertical grid coordinates
  call GRID_INDEX_setup
  call GRID_setup

  ! check if the namelist seetings are consistent
  if (MEM_NP /= PRC_NUM_X * PRC_NUM_Y) then
    write(6,*) 'MEM_NP should be equal to PRC_NUM_X * PRC_NUM_Y.'
    stop
  else if (IMAX /= nlonsub) then
    write(6,*) 'IMAX should be equal to nlonsub.'
    stop
  else if (JMAX /= nlatsub) then
    write(6,*) 'JMAX should be equal to nlatsub.'
    stop
  else if (KMAX /= nlev) then
    write(6,*) 'KMAX should be equal to nlev.'
    stop
  end if

!  call LAND_GRID_INDEX_setup
!  call LAND_GRID_setup

!  call URBAN_GRID_INDEX_setup
!  call URBAN_GRID_setup

  ! setup file I/O
  call FILEIO_setup

  ! setup mpi communication
  call COMM_setup

  ! setup topography
!  call TOPO_setup
  ! setup land use category index/fraction
!  call LANDUSE_setup

  ! setup grid coordinates (real world)
!  call REAL_setup

  ! setup map projection
  call MPRJ_setup( GRID_DOMAIN_CENTER_X, GRID_DOMAIN_CENTER_Y )

  ! setup history file I/O
  rankidx(1) = PRC_2Drank(PRC_myrank, 1)
  rankidx(2) = PRC_2Drank(PRC_myrank, 2)
  call HistoryInit('','','',IMAX*JMAX*KMAX,PRC_master,PRC_myrank,rankidx)

  call PROF_rapend('Initialize')

  call PROF_rapstart('Main')

  return
end subroutine set_scalelib

!-----------------------------------------------------------------------
! Finish using SCALE library
!-----------------------------------------------------------------------
subroutine unset_scalelib
  use scale_process, only: &
    PRC_MPIfinish
  use gtool_file, only: &
    FileCloseAll
  implicit none

  call PROF_rapend('Main')

  call PROF_rapreport

  call FileCloseAll

  ! stop SCALE MPI
  call PRC_MPIfinish

  return
end subroutine unset_scalelib
!=======================================================================

END MODULE common_scalelib
