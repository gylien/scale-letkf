module common_scalerm
!===============================================================================
!
! [PURPOSE:] Common tools for SCALE-RM
!
!===============================================================================
!!$USE OMP_LIB
!  use mpi
!  use common
  use common_mpi
  use common_nml
  use common_scale, only: &
    modelname
  use common_mpi_scale, only: &
    mpi_timer, &
    nprocs_m, &
    nitmax, &
    myrank_to_mem, &
    myrank_to_pe, &
    myrank_use, &
    mydom, &
    MPI_COMM_u, nprocs_u, myrank_u, &
    MPI_COMM_a, nprocs_a, myrank_a, &
    MPI_COMM_d, nprocs_d, myrank_d
  implicit none
  public

contains

!-------------------------------------------------------------------------------
! Setup SCALE-RM
!-------------------------------------------------------------------------------
subroutine scalerm_setup(execname)
  use dc_log, only: &
    LogInit
  use scale_stdio, only: &
    IO_setup, &
    IO_LOG_setup, &
    IO_FID_CONF, &
    IO_FID_LOG, &
    IO_L, &
    H_LONG
  use scale_prof, only: &
    PROF_setup, &
    PROF_setprefx, &
    PROF_rapstart, &
    PROF_rapend

  use scale_process, only: &
    PRC_LOCAL_setup, &
!    PRC_mpi_alive, &
!    PRC_MPIstart, &
!    PRC_UNIVERSAL_setup, &
    PRC_MPIsplit, &
    PRC_GLOBAL_setup, &
    PRC_UNIVERSAL_IsMaster, &
    PRC_nprocs, &
    PRC_myrank, &
    PRC_masterrank, &
    PRC_DOMAIN_nlim
  use scale_rm_process, only: &
    PRC_setup, &
    PRC_2Drank
  use scale_const, only: &
    CONST_setup
  use scale_calendar, only: &
    CALENDAR_setup
  use scale_random, only: &
    RANDOM_setup
  use scale_grid_index, only: &
    GRID_INDEX_setup, &
    IMAX, &
    JMAX, &
    KMAX
  use scale_grid, only: &
    GRID_setup, &
    DX, &
    DY, &
    GRID_DOMAIN_CENTER_X, &
    GRID_DOMAIN_CENTER_Y
  use scale_grid_nest, only: &
    NEST_setup
  use scale_land_grid_index, only: &
    LAND_GRID_INDEX_setup
  use scale_land_grid, only: &
    LAND_GRID_setup
  use scale_urban_grid_index, only: &
    URBAN_GRID_INDEX_setup
  use scale_urban_grid, only: &
    URBAN_GRID_setup
  use scale_fileio, only: &
    FILEIO_setup
!    FILEIO_cleanup
  use scale_comm, only: &
    COMM_setup
!    COMM_cleanup
  use scale_topography, only: &
    TOPO_setup
  use scale_landuse, only: &
    LANDUSE_setup
  use scale_grid_real, only: &
    REAL_setup
  use scale_mapproj, only: &
    MPRJ_setup
  use scale_gridtrans, only: &
    GTRANS_setup
  use scale_interpolation, only: &
    INTERP_setup
  use scale_rm_statistics, only: &
    STAT_setup
  use scale_history, only: &
    HIST_setup
!    HIST_write
  use gtool_history, only: &
    HistoryInit
  use scale_monitor, only: &
    MONIT_setup
!    MONIT_write, &
!    MONIT_finalize
  use scale_external_input, only: &
    EXTIN_setup
  use scale_atmos_hydrostatic, only: &
    ATMOS_HYDROSTATIC_setup
  use scale_atmos_thermodyn, only: &
    ATMOS_THERMODYN_setup
  use scale_atmos_hydrometeor, only: &
    ATMOS_HYDROMETEOR_setup
  use scale_atmos_saturation, only: &
    ATMOS_SATURATION_setup
  use scale_bulkflux, only: &
    BULKFLUX_setup
  use scale_roughness, only: &
    ROUGHNESS_setup

  use mod_atmos_driver, only: &
    ATMOS_driver_config
  use scale_atmos_phy_mp, only: &
    ATMOS_PHY_MP_config
  use mod_admin_restart, only: &
    ADMIN_restart_setup
!    ADMIN_restart_write
  use mod_admin_time, only: &
    ADMIN_TIME_setup
!    ADMIN_TIME_checkstate, &
!    ADMIN_TIME_advance,    &
!    TIME_DOATMOS_step,     &
!    TIME_DOLAND_step,      &
!    TIME_DOURBAN_step,     &
!    TIME_DOOCEAN_step,     &
!    TIME_DOresume,         &
!    TIME_DOend
  use mod_atmos_admin, only: &
    ATMOS_admin_setup
!    ATMOS_do
  use mod_atmos_vars, only: &
    ATMOS_vars_setup
!    ATMOS_sw_check => ATMOS_RESTART_CHECK,    &
!    ATMOS_vars_restart_check
  use mod_atmos_driver, only: &
    ATMOS_driver_setup
!    ATMOS_driver,           &
!    ATMOS_driver_finalize
  use mod_ocean_admin, only: &
    OCEAN_admin_setup
!    OCEAN_do
  use mod_ocean_vars, only: &
    OCEAN_vars_setup
  use mod_ocean_driver, only: &
    OCEAN_driver_setup
!    OCEAN_driver
  use mod_land_admin, only: &
    LAND_admin_setup
!    LAND_do
  use mod_land_vars, only: &
    LAND_vars_setup
  use mod_land_driver, only: &
    LAND_driver_setup
!    LAND_driver
  use mod_urban_admin, only: &
    URBAN_admin_setup
!    URBAN_do
  use mod_urban_vars, only: &
    URBAN_vars_setup
  use mod_urban_driver, only: &
    URBAN_driver_setup
!    URBAN_driver
  use mod_cpl_admin, only: &
    CPL_admin_setup
  use mod_cpl_vars, only: &
    CPL_vars_setup
  use mod_user, only: &
    USER_config, &
    USER_setup
!    USER_step
  implicit none

  character(len=*), intent(in), optional :: execname

!  integer :: universal_comm
!  integer :: universal_nprocs
!  logical :: universal_master
  integer :: global_comm
  integer :: local_comm
  integer :: local_myrank
  logical :: local_ismaster
  integer :: intercomm_parent
  integer :: intercomm_child
  character(len=H_LONG) :: confname_domains(PRC_DOMAIN_nlim)
  character(len=H_LONG) :: confname_mydom
  character(len=H_LONG) :: confname

  integer :: color, key, it, its, ite, im, idom, ierr

  integer :: rankidx(2)
  integer :: HIST_item_limit    ! dummy
  integer :: HIST_variant_limit ! dummy

  character(len=7) :: execname_ = 'LETKF  '

  if (present(execname)) execname_ = execname

  call mpi_timer('', 2, barrier=MPI_COMM_WORLD)

  ! Communicator for all processes used
  !-----------------------------------------------------------------------------

  if (myrank_use) then
    color = 0
    key   = (myrank_to_mem(1) - 1) * nprocs_m + myrank_to_pe
!    key   = myrank
  else
    color = MPI_UNDEFINED
    key   = MPI_UNDEFINED
  end if

  call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, MPI_COMM_u, ierr)

  if (.not. myrank_use) then
    write (6, '(A,I6.6,A)') 'MYRANK=', myrank, ': This process is not used!'
    return
  end if

  call MPI_COMM_SIZE(MPI_COMM_u, nprocs_u, ierr)
  call MPI_COMM_RANK(MPI_COMM_u, myrank_u, ierr)

  call mpi_timer('scalerm_setup:mpi_comm_split_u:', 2)

  ! Communicator for all domains of single members
  !-----------------------------------------------------------------------------

  ! start SCALE MPI
!  call PRC_MPIstart( universal_comm ) ! [OUT]

!  PRC_mpi_alive = .true.
!  universal_comm = MPI_COMM_u

!  call PRC_UNIVERSAL_setup( universal_comm,   & ! [IN]
!                            universal_nprocs, & ! [OUT]
!                            universal_master  ) ! [OUT]

!  if (myrank_to_mem(1) >= 1) then
    color = myrank_to_mem(1) - 1
    key   = myrank_to_pe
!  else
!    color = MPI_UNDEFINED
!    key   = MPI_UNDEFINED
!  endif

  call MPI_COMM_SPLIT(MPI_COMM_u, color, key, global_comm, ierr)

  call PRC_GLOBAL_setup( .false.,    & ! [IN]
                         global_comm ) ! [IN]

  call mpi_timer('scalerm_setup:mpi_comm_split_d_global:', 2)

  ! Communicator for one domain
  !-----------------------------------------------------------------------------

  do idom = 1, NUM_DOMAIN
    confname_domains(idom) = trim(CONF_FILES)
    call filename_replace_dom(confname_domains(idom), idom)
  end do

  !--- split for nesting
  ! communicator split for nesting domains
  call PRC_MPIsplit( global_comm,      & ! [IN]
                     NUM_DOMAIN,       & ! [IN]
                     PRC_DOMAINS(:),   & ! [IN]
                     confname_domains(:), & ! [IN]
                     .false.,          & ! [IN]
                     .false.,          & ! [IN] flag bulk_split
                     COLOR_REORDER,    & ! [IN]
                     local_comm,       & ! [OUT]
                     intercomm_parent, & ! [OUT]
                     intercomm_child,  & ! [OUT]
                     confname_mydom    ) ! [OUT]

  MPI_COMM_d = local_comm

  do idom = 1, NUM_DOMAIN
    if (trim(confname_mydom) == trim(confname_domains(idom))) then
      mydom = idom
      exit
    end if
  end do

#ifdef DEBUG
  if (mydom <= 0) then
    write (6, '(A)'), '[Error] Cannot determine my domain ID.'
    stop
  end if
#endif

  ! Setup standard I/O and read LETKF namelists
  !-----------------------------------------------------------------------------

  if (execname_ == 'SCALERM') then

    if (MEMBER_ITER == 0) then
      its = 1
      ite = nitmax
    else
      its = MEMBER_ITER
      ite = MEMBER_ITER
    end if

    do it = its, ite
      im = myrank_to_mem(it)
      if (im >= 1 .and. im <= MEMBER_RUN) then
        confname = confname_mydom
        if (CONF_FILES_SEQNUM) then
          call filename_replace_mem(confname, im)
        else
          if (im <= MEMBER) then
            call filename_replace_mem(confname, im)
          else if (im == MEMBER+1) then
            call filename_replace_mem(confname, memf_mean)
          else if (im == MEMBER+2) then
            call filename_replace_mem(confname, memf_mdet)
          end if
        end if
        write (6,'(A,I6.6,2A)') 'MYRANK ', myrank, ' is running a model with configuration file: ', trim(confname)
      end if
    end do ! [ it = its, ite ]

    ! setup standard I/O
    call IO_setup( modelname, .true., trim(confname) )
  else
    if (mydom >= 2) then ! In d01, keep using the original launcher config file; skip re-opening config files here
      ! setup standard I/O (for inner domains)
      call IO_setup( modelname, .true., confname_mydom )

!      call read_nml_log
!      call read_nml_model
!      call read_nml_ensemble
!      call read_nml_process
    end if

    select case (execname_)
    case ('LETKF  ')
      call read_nml_obs_error
      call read_nml_obsope
      call read_nml_letkf
      call read_nml_letkf_obs
      call read_nml_letkf_var_local
      call read_nml_letkf_monitor
      call read_nml_letkf_radar
      call read_nml_letkf_h08
    case ('OBSOPE ', 'OBSMAKE')
      call read_nml_obs_error
      call read_nml_obsope
      call read_nml_letkf_radar
      call read_nml_letkf_h08
    case ('OBSSIM ')
      call read_nml_obssim
      call read_nml_letkf_radar
      call read_nml_letkf_h08
    end select
  end if

  !-----------------------------------------------------------------------------

  ! setup MPI
  call PRC_LOCAL_setup( local_comm,    & ! [IN]
                        local_myrank,  & ! [OUT]
                        local_ismaster ) ! [OUT]

!  call MPI_COMM_SIZE(MPI_COMM_d, nprocs_d, ierr)
  nprocs_d = PRC_nprocs
!  call MPI_COMM_RANK(MPI_COMM_d, myrank_d, ierr)
!  myrank_d = PRC_myrank
  myrank_d = local_myrank

  call mpi_timer('scalerm_setup:mpi_comm_split_d_local:', 2)

  ! Communicator for all processes for single domains
  !-----------------------------------------------------------------------------

!  if (mydom > 0) then
    color = mydom - 1
    key   = (myrank_to_mem(1) - 1) * nprocs_m + myrank_to_pe
!    key   = myrank
!  else
!    color = MPI_UNDEFINED
!    key   = MPI_UNDEFINED
!  end if

  call MPI_COMM_SPLIT(MPI_COMM_u, color, key, MPI_COMM_a, ierr)

  call MPI_COMM_SIZE(MPI_COMM_a, nprocs_a, ierr)
  call MPI_COMM_RANK(MPI_COMM_a, myrank_a, ierr)

  call mpi_timer('scalerm_setup:mpi_comm_split_a:', 2)

  ! Setup scalelib LOG output (only for the universal master rank)
  !-----------------------------------------------------------------------------

  ! setup Log
  if (execname_ == 'SCALERM') then
    call IO_LOG_setup( local_myrank, local_ismaster )
  else
    call IO_LOG_setup( local_myrank, PRC_UNIVERSAL_IsMaster )
  end if
  call LogInit( IO_FID_CONF, IO_FID_LOG, IO_L )

  call mpi_timer('scalerm_setup:log_setup_init:', 2)

  ! Other scalelib setups
  !-----------------------------------------------------------------------------

  ! setup process
  call PRC_setup

  if (execname_ == 'SCALERM') then
    ! setup PROF
    call PROF_setup

    ! profiler start
    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize', 0)
  end if

  ! setup constants
  call CONST_setup

  if (execname_ == 'SCALERM') then
    ! setup calendar
    call CALENDAR_setup

    ! setup random number
    call RANDOM_setup
  end if

  ! setup horizontal/vertical grid coordinates (cartesian,idealized)
  call GRID_INDEX_setup
  call GRID_setup

#ifdef PNETCDF
  call LAND_GRID_INDEX_setup
  if (execname_ == 'SCALERM') then
    call LAND_GRID_setup
  end if

  call URBAN_GRID_INDEX_setup
  if (execname_ == 'SCALERM') then
    call URBAN_GRID_setup
  end if
#else
  if (execname_ == 'SCALERM') then
    call LAND_GRID_INDEX_setup
    call LAND_GRID_setup

    call URBAN_GRID_INDEX_setup
    call URBAN_GRID_setup
  end if
#endif

  if (execname_ == 'SCALERM') then
    ! setup submodel administrator
    call ATMOS_admin_setup
    call OCEAN_admin_setup
    call LAND_admin_setup
    call URBAN_admin_setup
    call CPL_admin_setup
  end if

  ! setup tracer index
  call ATMOS_HYDROMETEOR_setup

  if (execname_ == 'SCALERM') then
    call ATMOS_driver_config
    call USER_config
  else
!   call ATMOS_driver_config -->
!     call ATMOS_PHY_MP_driver_config -->
!       if ( ATMOS_sw_phy_mp ) then
          call ATMOS_PHY_MP_config('TOMITA08') !!!!!!!!!!!!!!! tentative
!       end if
!     <-- ATMOS_PHY_MP_driver_config
!   <-- ATMOS_driver_config
  end if

  ! setup file I/O
  call FILEIO_setup

  ! setup mpi communication
  call COMM_setup

  if (execname_ == 'SCALERM') then
    ! setup topography
    call TOPO_setup

    ! setup land use category index/fraction
    call LANDUSE_setup
  end if

  ! setup grid coordinates (real world)
  if (execname_ == 'SCALERM') then
    call REAL_setup
  else
!   call REAL_setup -->
      ! setup map projection
      call MPRJ_setup( GRID_DOMAIN_CENTER_X, GRID_DOMAIN_CENTER_Y )
!   <-- REAL_setup
  end if

  if (execname_ == 'SCALERM') then
    ! setup grid transfer metrics (uses in ATMOS_dynamics)
    call GTRANS_setup

    ! setup Z-ZS interpolation factor (uses in History)
    call INTERP_setup

    ! setup restart
    call ADMIN_restart_setup

    ! setup time
    call ADMIN_TIME_setup( setup_TimeIntegration = .true. )

    ! setup statistics
    call STAT_setup
  end if

  ! setup history I/O
  if (execname_ == 'SCALERM') then
    call HIST_setup
  else
!   call HIST_setup -->
      ! setup history file I/O
      rankidx(1) = PRC_2Drank(PRC_myrank, 1)
      rankidx(2) = PRC_2Drank(PRC_myrank, 2)

      call HistoryInit( HIST_item_limit,                  & ! [OUT]
                        HIST_variant_limit,               & ! [OUT]
                        IMAX, JMAX, KMAX,                 & ! [IN]
                        PRC_masterrank,                   & ! [IN]
                        PRC_myrank,                       & ! [IN]
                        rankidx,                          & ! [IN]
                        '',                               & ! [IN]
                        '',                               & ! [IN]
                        '',                               & ! [IN]
                        0.0d0,                            & ! [IN]
                        1.0d0,                            & ! [IN]
                        default_basename='history',       & ! [IN]
                        default_zcoord = 'model',         & ! [IN]
                        default_tinterval = 1.0d0,        & ! [IN]
                        namelist_fid=IO_FID_CONF          ) ! [IN]
!   <-- HIST_setup
  end if

  if (execname_ == 'SCALERM') then
    ! setup monitor I/O
    call MONIT_setup

    ! setup external in
    call EXTIN_setup( 'RM' )

    ! setup nesting grid
    call NEST_setup ( intercomm_parent, intercomm_child )

    ! setup common tools
    call ATMOS_HYDROSTATIC_setup
    call ATMOS_THERMODYN_setup
    call ATMOS_SATURATION_setup

    call BULKFLUX_setup( sqrt(DX**2+DY**2) )
    call ROUGHNESS_setup

    ! setup variable container
    call ATMOS_vars_setup
    call OCEAN_vars_setup
    call LAND_vars_setup
    call URBAN_vars_setup
    call CPL_vars_setup

    ! setup submodel driver
    call ATMOS_driver_setup
    call OCEAN_driver_setup
    call LAND_driver_setup
    call URBAN_driver_setup

    call USER_setup

    call PROF_rapend('Initialize', 0)
  end if

  call mpi_timer('scalerm_setup:other_setup:', 2)

  return
end subroutine scalerm_setup

!-------------------------------------------------------------------------------
! Finalize SCALE-RM
!-------------------------------------------------------------------------------
subroutine scalerm_finalize(execname)
!  use scale_stdio, only: &
!    IO_FID_CONF, &
!    IO_FID_LOG, &
!    IO_L, &
!    IO_FID_STDOUT
  use scale_prof, only: &
!    PROF_setprefx, &
    PROF_rapstart, &
    PROF_rapend, &
    PROF_rapreport
  use scale_process, only: &
    PRC_MPIfinish
  use scale_fileio, only: &
    FILEIO_cleanup
  use gtool_file, only: &
    FileCloseAll
  use scale_comm, only: &
    COMM_cleanup
  use scale_monitor, only: &
    MONIT_finalize
  use mod_atmos_vars, only: &
    ATMOS_sw_check => ATMOS_RESTART_CHECK, &
    ATMOS_vars_restart_check
  implicit none

  character(len=*), intent(in), optional :: execname
  character(len=7) :: execname_ = 'LETKF  '
  integer :: ierr

  if (present(execname)) execname_ = execname

  if (myrank_use) then
    if (execname_ == 'SCALERM') then
      ! check data
      if( ATMOS_sw_check ) call ATMOS_vars_restart_check

      call PROF_rapstart('Monit', 2)
      call MONIT_finalize
      call PROF_rapend  ('Monit', 2)

      call PROF_rapstart('File', 2)
      ! clean up resource allocated for I/O
      call FILEIO_cleanup

      call COMM_cleanup
    end if
  end if

  call FileCloseAll

  if (execname_ == 'SCALERM') then
    call PROF_rapend  ('File', 2)

    call PROF_rapend  ('All', 1)

    call PROF_rapreport
  end if

  if (myrank_use) then
    call MPI_COMM_FREE(MPI_COMM_d, ierr)
    call MPI_COMM_FREE(MPI_COMM_a, ierr)
    call MPI_COMM_FREE(MPI_COMM_u, ierr)
  end if

  ! stop MPI
!  call PRC_MPIfinish

!  ! Close logfile, configfile
!  if ( IO_L ) then
!    if( IO_FID_LOG /= IO_FID_STDOUT ) close(IO_FID_LOG)
!  endif
!  close(IO_FID_CONF)

  return
end subroutine scalerm_finalize

!-------------------------------------------------------------------------------
! resume_state
!-------------------------------------------------------------------------------
subroutine resume_state
  use mod_atmos_driver, only: &
    ATMOS_driver_resume1, &
    ATMOS_driver_resume2, &
    ATMOS_SURFACE_SET
  use mod_ocean_driver, only: &
    OCEAN_driver_resume, &
    OCEAN_SURFACE_SET
  use mod_land_driver, only: &
    LAND_driver_resume, &
    LAND_SURFACE_SET
  use mod_urban_driver, only: &
    URBAN_driver_resume, &
    URBAN_SURFACE_SET
  use mod_atmos_vars, only: &
    ATMOS_vars_diagnostics,     &
    ATMOS_vars_history_setpres, &
    ATMOS_vars_restart_read
  use mod_ocean_vars, only: &
    OCEAN_vars_restart_read
  use mod_land_vars, only: &
    LAND_vars_restart_read
  use mod_urban_vars, only: &
    URBAN_vars_restart_read
  use mod_user, only: &
    USER_resume0, &
    USER_resume
  use mod_atmos_admin, only: &
    ATMOS_do
  use mod_ocean_admin, only: &
    OCEAN_do
  use mod_land_admin, only: &
    LAND_do
  use mod_urban_admin, only: &
    URBAN_do
  use mod_admin_restart, only: &
    ADMIN_restart_read
  implicit none
  !---------------------------------------------------------------------------

  ! read restart data
  call ADMIN_restart_read

  ! setup user-defined procedure before setup of other components
  call USER_resume0

  if ( ATMOS_do ) then
    ! calc diagnostics
    call ATMOS_vars_diagnostics
    call ATMOS_vars_history_setpres
  endif

  ! setup surface condition
  if( ATMOS_do ) call ATMOS_SURFACE_SET( countup=.false. )
  if( OCEAN_do ) call OCEAN_SURFACE_SET( countup=.false. )
  if( LAND_do  ) call LAND_SURFACE_SET ( countup=.false. )
  if( URBAN_do ) call URBAN_SURFACE_SET( countup=.false. )

  ! setup submodel driver
  if( ATMOS_do ) call ATMOS_driver_resume1
  if( OCEAN_do ) call OCEAN_driver_resume
  if( LAND_do  ) call LAND_driver_resume
  if( URBAN_do ) call URBAN_driver_resume
  if( ATMOS_do ) call ATMOS_driver_resume2

  ! setup user-defined procedure
  call USER_resume

  return
end subroutine resume_state

!===============================================================================
end module common_scalerm
