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
    modelname, &
    confname
  use common_mpi_scale, only: &
    mpi_timer, &
    nprocs_m, &
    nitmax, &
    myrank_to_mem, &
    myrank_to_pe, &
    myrank_use, &
    myrank_use_da, &
    mydom, &
    MPI_COMM_u, nprocs_u, myrank_u, &
    MPI_COMM_a, nprocs_a, myrank_a, &
    MPI_COMM_da, nprocs_da, myrank_da, &
    MPI_COMM_d, nprocs_d, myrank_d
  implicit none
  public

  integer, save                :: scalerm_mem = -1
  character(len=memflen), save :: scalerm_memf = '????'
  logical, save                :: scalerm_run = .false.

contains

!-------------------------------------------------------------------------------
! Setup SCALE-RM
!-------------------------------------------------------------------------------
subroutine scalerm_setup(execname)
  use scale_io, only: &
    IO_setup, &
    IO_LOG_setup, &
#ifdef SCALEUV
    H_LONG, &
    IO_filename_replace_setup
#else
    H_LONG
#endif
  use scale_prof, only: &
    PROF_setup, &
    PROF_setprefx, &
    PROF_rapstart, &
    PROF_rapend

!  use scale_file, only: &
!    FILE_Close_All
  use scale_prc, only: &
!    PRC_mpi_alive, &
!    PRC_MPIstart, &
!    PRC_UNIVERSAL_setup, &
    PRC_MPIsplit, &
    PRC_GLOBAL_setup, &
    PRC_LOCAL_setup, &
    PRC_UNIVERSAL_IsMaster, &
    PRC_nprocs, &
    PRC_myrank, &
    PRC_DOMAIN_nlim
  use scale_prc_cartesC, only: &
    PRC_CARTESC_setup
  use scale_const, only: &
    CONST_setup
  use scale_calendar, only: &
    CALENDAR_setup
  use scale_random, only: &
    RANDOM_setup
  use scale_atmos_grid_cartesC_index, only: &
    ATMOS_GRID_CARTESC_INDEX_setup
  use scale_atmos_grid_cartesC, only: &
    ATMOS_GRID_CARTESC_setup, &
    ATMOS_GRID_CARTESC_DOMAIN_CENTER_X, &
    ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
    DX, &
    DY
  use scale_comm_cartesC_nest, only: &
    COMM_CARTESC_NEST_setup
  use scale_ocean_grid_cartesC_index, only: &
    OCEAN_GRID_CARTESC_INDEX_setup
  use scale_ocean_grid_cartesC, only: &
    OCEAN_GRID_CARTESC_setup
  use scale_ocean_grid_cartesC_real, only: &
    OCEAN_GRID_CARTESC_REAL_setup
  use scale_land_grid_cartesC_index, only: &
    LAND_GRID_CARTESC_INDEX_setup
  use scale_land_grid_cartesC, only: &
    LAND_GRID_CARTESC_setup
  use scale_land_grid_cartesC_real, only: &
    LAND_GRID_CARTESC_REAL_setup
  use scale_urban_grid_cartesC_index, only: &
    URBAN_GRID_CARTESC_INDEX_setup
  use scale_urban_grid_cartesC, only: &
    URBAN_GRID_CARTESC_setup
  use scale_urban_grid_cartesC_real, only: &
     URBAN_GRID_CARTESC_REAL_setup
  use scale_file_cartesC, only: &
    FILE_CARTESC_setup
!    FILE_CARTESC_cleanup
  use scale_file, only: &
    FILE_setup
  use scale_comm_cartesC, only: &
    COMM_setup
!    COMM_cleanup
  use scale_topography, only: &
    TOPO_setup
  use scale_landuse, only: &
    LANDUSE_setup
  use scale_atmos_grid_cartesC_real, only: &
    ATMOS_GRID_CARTESC_REAL_setup
  use scale_mapprojection, only: &
    MAPPROJECTION_setup
  use scale_atmos_grid_cartesC_metric, only: &
    ATMOS_GRID_CARTESC_METRIC_setup
  use scale_statistics, only: &
    STATISTICS_setup
  use scale_time, only: &
    TIME_NOWDATE, &
    TIME_NOWMS,   &
    TIME_NOWSTEP, &
    TIME_DTSEC
!  use scale_file_history, only: &
!    FILE_HISTORY_write, &
!    FILE_HISTORY_set_nowdate, &
!    FILE_HISTORY_finalize
  use scale_file_history_cartesC, only: &
    FILE_HISTORY_CARTESC_setup
  use scale_monitor_cartesC, only: &
    MONITOR_CARTESC_setup
!    MONIT_write, &
!    MONIT_finalize
  use scale_file_external_input_cartesC, only: &
    FILE_EXTERNAL_INPUT_CARTESC_setup
  use scale_atmos_hydrostatic, only: &
    ATMOS_HYDROSTATIC_setup
  use scale_atmos_thermodyn, only: &
    ATMOS_THERMODYN_setup
  use scale_atmos_hydrometeor, only: &
    ATMOS_HYDROMETEOR_setup!, &
!    ATMOS_HYDROMETEOR_regist
!    I_QC, &
!    I_QR, &
!    I_QI, &
!    I_QS, &
!    I_QG
!  use scale_atmos_phy_mp, only: &
!    ATMOS_PHY_MP_config
!  use scale_atmos_phy_mp_tomita08, only: &
!!    ATMOS_PHY_MP_TOMITA08_ntracers, &
!    ATMOS_PHY_MP_TOMITA08_nwaters, &
!    ATMOS_PHY_MP_TOMITA08_nices, &
!    ATMOS_PHY_MP_TOMITA08_tracer_names, &
!    ATMOS_PHY_MP_TOMITA08_tracer_descriptions, &
!    ATMOS_PHY_MP_TOMITA08_tracer_units
!!  use mod_atmos_phy_mp_driver, only: &
!!    ATMOS_PHY_MP_driver_tracer_setup
!!  use mod_atmos_admin, only: &
!!    ATMOS_PHY_MP_TYPE, &
!!    ATMOS_sw_phy_mp
!!  use mod_atmos_phy_mp_vars, only: &
!!!    QA_MP, &
!!    QS_MP, &
!!!    QE_MP
  use scale_atmos_saturation, only: &
    ATMOS_SATURATION_setup
  use scale_bulkflux, only: &
    BULKFLUX_setup
!  use scale_roughness, only: &
!    ROUGHNESS_setup
  use mod_atmos_driver, only: &
    ATMOS_driver_tracer_setup
  use mod_atmos_phy_mp_vars, only: &
    QA_MP
  use mod_admin_restart, only: &
    ADMIN_restart_setup
!    ADMIN_restart_write
  use mod_admin_time, only: &
    ADMIN_TIME_setup
!    ADMIN_TIME_checkstate, &
!    ADMIN_TIME_advance, &
!    TIME_DOATMOS_step, &
!    TIME_DOLAND_step, &
!    TIME_DOURBAN_step, &
!    TIME_DOOCEAN_step, &
!    TIME_DOresume, &
!    TIME_DOend
  use mod_atmos_admin, only: &
    ATMOS_admin_setup, &
    ATMOS_do, &
    ATMOS_PHY_MP_TYPE
  use mod_atmos_vars, only: &
    ATMOS_vars_setup
!    ATMOS_sw_check => ATMOS_RESTART_CHECK, &
!    ATMOS_vars_restart_check
  use mod_atmos_driver, only: &
    ATMOS_driver_setup
!    ATMOS_driver, &
!    ATMOS_driver_finalize
  use mod_ocean_admin, only: &
    OCEAN_admin_setup, &
    OCEAN_do
  use mod_ocean_vars, only: &
    OCEAN_vars_setup
  use mod_ocean_driver, only: &
    OCEAN_driver_setup
!    OCEAN_driver
  use mod_land_admin, only: &
    LAND_admin_setup, &
    LAND_do
  use mod_land_vars, only: &
    LAND_vars_setup
  use mod_land_driver, only: &
    LAND_driver_setup
!    LAND_driver
  use mod_urban_admin, only: &
    URBAN_admin_setup, &
    URBAN_do,          &
    URBAN_land
  use mod_urban_vars, only: &
    URBAN_vars_setup
  use mod_urban_driver, only: &
    URBAN_driver_setup
!    URBAN_driver
  use mod_lake_admin, only: &
    LAKE_admin_setup, &
    LAKE_do
  use mod_cpl_admin, only: &
    CPL_admin_setup, &
    CPL_sw
  use mod_cpl_vars, only: &
    CPL_vars_setup
  use mod_user, only: &
    USER_tracer_setup,  &
!    USER_config, &
    USER_setup
!    USER_step
  use mod_convert, only: &
    CONVERT_setup
!    CONVERT
  use mod_mktopo, only: &
    MKTOPO_setup
!    MKTOPO
  use mod_mkinit, only: &
    MKINIT_setup
!    MKINIT
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
  character(len=H_LONG) :: confname_new
  character(len=H_LONG) :: confname_new2
  integer :: qs_mp_dummy

  integer :: color, key, idom, ierr
  integer :: color_da, key_da, mem_da
  character(len=2) :: fmttmp

  logical :: exec_model
  logical :: exec_modelonly

  character(len=7) :: execname_ = 'LETKF  '

  if (present(execname)) execname_ = execname

  exec_model = .false.
  if (execname_ == 'SCALERM' .or. execname_ == 'RMPREP ' .or. execname_ == 'DACYCLE') then
    exec_model = .true.
  end if
  exec_modelonly = .false.
  if (execname_ == 'SCALERM' .or. execname_ == 'RMPREP ') then
    exec_modelonly = .true.
  end if

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

  ! Input fake confname, unique for each domains, to PRC_MPIsplit,
  ! to easily determine 'my domain' later
  confname_domains(1:NUM_DOMAIN) = domf_notation
  do idom = 1, NUM_DOMAIN
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

  ! Determine my domain
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

  confname_new = confname

  ! Set real confname for my domain
  if (trim(CONF_FILES) /= '') then
    if (exec_model .or. mydom >= 2) then
      confname_new = trim(CONF_FILES)
      call filename_replace_dom(confname_new, mydom)
    end if
  end if

  ! Setup standard I/O
  !-----------------------------------------------------------------------------

  if (exec_model) then
    if (MEMBER_ITER == 0) then
      write (6, '(A)') '[Warning] Currently this code can only run a single member iteration; reset MEMBER_ITER to 1!'
      MEMBER_ITER = 1
    end if

    if (myrank_to_mem(MEMBER_ITER) >= 1 .and. myrank_to_mem(MEMBER_ITER) <= MEMBER_RUN) then
      scalerm_run = .true.
      if (MEMBER_SEQ(1) == -1) then
        scalerm_mem = myrank_to_mem(MEMBER_ITER)
      else
        scalerm_mem = MEMBER_SEQ(myrank_to_mem(MEMBER_ITER))
      end if

      mem_da = MEMBER
      if (ENS_WITH_MEAN .and. ENS_WITH_MDET) then
        mem_da = MEMBER + 2
      elseif (ENS_WITH_MEAN .or. ENS_WITH_MDET) then
        mem_da = MEMBER + 1
      endif

      myrank_use_da = .true.

      if (scalerm_mem >= 1 .and. scalerm_mem <= MEMBER) then
        write (fmttmp, '(I2)') memflen
        write (scalerm_memf, '(I'//trim(fmttmp)//'.'//trim(fmttmp)//')') scalerm_mem
      else if (scalerm_mem == MEMBER+1 .and. ENS_WITH_MEAN) then
        scalerm_memf = memf_mean
      else if (scalerm_mem == MEMBER+2 .and. ENS_WITH_MDET) then
        scalerm_memf = memf_mdet
      else if (scalerm_mem <= MEMBER+2+MAX_DACYCLE_RUN_FCST .and. scalerm_mem > mem_da) then
        scalerm_memf = memf_mean
        myrank_use_da = .false.
      else
        write (6, '(A,I7)') '[Error] Invalid member number for this rank:', scalerm_mem
        write (6, '(A,I7)') '        MEMBER =', MEMBER
        stop 1
      end if

#ifdef SCALEUV
      call IO_filename_replace_setup(memf_notation, scalerm_memf)
#endif

      if (trim(CONF_FILES) /= '') then
        confname_new2 = confname_new
        call filename_replace_mem(confname_new2, scalerm_memf)
        if (mydom == 1) then
          if (trim(confname_new2) /= trim(confname_new)) then ! In domain #1, reset config file only when
            confname_new = confname_new2                      ! <member> keyword exists in CONF_FILES
          else
            confname_new = confname
          end if
        else
          confname_new = confname_new2 ! In other domains, always reset config file as long as CONF_FILES is set
        end if
      end if
    end if
  end if ! [ exec_model ]

  if ((.not. exec_model) .or. scalerm_run) then
    ! setup standard I/O: Re-open the new config file; change of IO_LOG_BASENAME is effective in this step
    confname = confname_new
    call IO_setup( modelname, trim(confname) )

  !  call read_nml_log
  !  call read_nml_model
  !  call read_nml_ensemble
  !  call read_nml_process

    if (myrank == 0) then
      if (scalerm_run) then
        write (6, '(A,I6.6,2A)') '[Info] MYRANK = ', myrank, ' is running SCALE using configuration file: ', trim(confname)
      else
        write (6, '(A,I6.6,2A)') '[Info] MYRANK = ', myrank, ' is using configuration file: ', trim(confname)
      end if
    endif
  else
    write (6, '(A,I6.6,A)') '[Info] MYRANK = ', myrank, ' is not used for SCALE!'
  end if

  call mpi_timer('scalerm_setup:standard I/O:', 2)

  ! Read LETKF namelists
  !-----------------------------------------------------------------------------

  select case (execname_)
  case ('DACYCLE')
    call read_nml_dacycle
    call read_nml_obs_error
    call read_nml_obsope
    call read_nml_letkf
    call read_nml_letkf_obs
    call read_nml_letkf_var_local
    call read_nml_letkf_monitor
    call read_nml_letkf_radar
    call read_nml_letkf_h08
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

  call mpi_timer('scalerm_setup:read_nml:', 2)

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

  if (mydom > 0) then
    color = mydom - 1
    key   = (myrank_to_mem(1) - 1) * nprocs_m + myrank_to_pe
  else
    color = MPI_UNDEFINED
    key   = MPI_UNDEFINED
  end if

  call MPI_COMM_SPLIT(MPI_COMM_u, color, key, MPI_COMM_a, ierr)

  call MPI_COMM_SIZE(MPI_COMM_a, nprocs_a, ierr)
  call MPI_COMM_RANK(MPI_COMM_a, myrank_a, ierr)

  call mpi_timer('scalerm_setup:mpi_comm_split_a:', 2)

! Define another MPI communicator in which processes used for DA
!  color_da = 0:used for DA, 1: not used for DA
  color_da = 0
  key_da = key
  if (.not. myrank_use_da) then
    color_da = 1
    key_da = myrank_a - mem_da * nprocs_m
  endif

  call MPI_COMM_SPLIT(MPI_COMM_a, color_da, key_da, MPI_COMM_da, ierr)

  call MPI_COMM_SIZE(MPI_COMM_da, nprocs_da, ierr)
  call MPI_COMM_RANK(MPI_COMM_da, myrank_da, ierr)

  call mpi_timer('scalerm_setup:mpi_comm_split_da:', 2)

  if (exec_modelonly .and. (.not. scalerm_run)) then
    return
  end if

  ! Setup scalelib LOG output (only for the universal master rank)
  !-----------------------------------------------------------------------------

  ! setup Log
  if (exec_model) then
    call IO_LOG_setup( local_myrank, local_ismaster )
  else
    call IO_LOG_setup( local_myrank, PRC_UNIVERSAL_IsMaster )
  end if

  call mpi_timer('scalerm_setup:log_setup_init:', 2)

  ! Other scalelib setups
  !-----------------------------------------------------------------------------

  ! setup process
  call PRC_CARTESC_setup

  if (exec_model .and. scalerm_run) then
    ! setup PROF
    call PROF_setup

    ! profiler start
    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize', 0)
  end if

  ! setup constants
  call CONST_setup

  if (exec_model .and. scalerm_run) then
    ! setup calendar
    call CALENDAR_setup

    ! setup random number
    call RANDOM_setup

    ! setup submodel administrator
    call ATMOS_admin_setup
    call OCEAN_admin_setup
    call LAND_admin_setup
    call URBAN_admin_setup
    call LAKE_admin_setup
    call CPL_admin_setup
  end if

  ! setup horizontal/vertical grid coordinates (cartesian,idealized)
  call ATMOS_GRID_CARTESC_INDEX_setup
  call ATMOS_GRID_CARTESC_setup

  if (exec_model .and. scalerm_run) then
    if ( OCEAN_do ) then
       call OCEAN_GRID_CARTESC_INDEX_setup
       call OCEAN_GRID_CARTESC_setup
    endif

    if ( LAND_do ) then
       call LAND_GRID_CARTESC_INDEX_setup
       call LAND_GRID_CARTESC_setup
    endif

    if ( URBAN_do ) then
       call URBAN_GRID_CARTESC_INDEX_setup
       call URBAN_GRID_CARTESC_setup
    endif

  end if

!#ifdef PNETCDF
!  call LAND_GRID_CARTESC_INDEX_setup
!  if (exec_model .and. scalerm_run) then
!    call LAND_GRID_CARTESC_setup
!  end if
!
!  call URBAN_GRID_CARTESC_INDEX_setup
!  if (exec_model .and. scalerm_run) then
!    call URBAN_GRID_CARTESC_setup
!  end if
!#else
!  if (exec_model .and. scalerm_run) then
!    call LAND_GRID_CARTESC_INDEX_setup
!    call LAND_GRID_CARTESC_setup
!
!    call URBAN_GRID_CARTESC_INDEX_setup
!    call URBAN_GRID_CARTESC_setup
!  end if
!#endif

  ! setup tracer index
  call ATMOS_HYDROMETEOR_setup

  if (exec_model) then
    if (scalerm_run) then
      call ATMOS_driver_tracer_setup
      call USER_tracer_setup
    end if
  else
!!   call ATMOS_driver_config -->
!!     call ATMOS_PHY_MP_driver_tracer_setup -->
!!        if ( ATMOS_sw_phy_mp ) then
!!          select case ( ATMOS_PHY_MP_TYPE )
!!          case ( 'TOMITA08' )
!            call ATMOS_HYDROMETEOR_regist( &
!!                 QS_MP,                                        & ! [OUT]
!                 qs_mp_dummy,                                  & ! [OUT]
!                 ATMOS_PHY_MP_TOMITA08_nwaters,                & ! [IN]
!                 ATMOS_PHY_MP_TOMITA08_nices,                  & ! [IN]
!                 ATMOS_PHY_MP_TOMITA08_tracer_names(:),        & ! [IN]
!                 ATMOS_PHY_MP_TOMITA08_tracer_descriptions(:), & ! [IN]
!                 ATMOS_PHY_MP_TOMITA08_tracer_units(:)         ) ! [IN]
!!            QA_MP = ATMOS_PHY_MP_TOMITA08_ntracers
!!            I_QC = QS_MP+1
!!            I_QR = QS_MP+2
!!            I_QI = QS_MP+3
!!            I_QS = QS_MP+4
!!            I_QG = QS_MP+5
!!          case default
!!            call ATMOS_PHY_MP_config( ATMOS_PHY_MP_TYPE )
!!          end select
!!        end if
!!        QE_MP = QS_MP + QA_MP - 1
!!     <-- ATMOS_PHY_MP_driver_tracer_setup
!!   <-- ATMOS_driver_config
  end if

  ! setup file I/O
  if (exec_model) then
    if (scalerm_run) then
      call FILE_CARTESC_setup
    end if
  else
!   call FILE_CARTESC_setup -->
      call FILE_setup( PRC_myrank )
!   <-- FILE_CARTESC_setup
  end if

  ! setup mpi communication
  call COMM_setup

  if (exec_model .and. scalerm_run) then
    ! setup topography
    call TOPO_setup

    ! setup land use category index/fraction
    call LANDUSE_setup( OCEAN_do, (.not. URBAN_land), LAKE_do )
  end if

  ! setup grid coordinates (real world)
  if (exec_model) then
    if (scalerm_run) then
#ifdef SCALEUV
      call ATMOS_GRID_CARTESC_REAL_setup( catalogue_output = (myrank_to_mem(1) == 1) ) ! Only output catalogue file in the first member of this execution
#else
      call ATMOS_GRID_CARTESC_REAL_setup
#endif

      ! setup grid transfer metrics (uses in ATMOS_dynamics)
      call ATMOS_GRID_CARTESC_METRIC_setup

      if ( OCEAN_do ) call OCEAN_GRID_CARTESC_REAL_setup
      if ( LAND_do  ) call LAND_GRID_CARTESC_REAL_setup
      if ( URBAN_do ) call URBAN_GRID_CARTESC_REAL_setup

    end if
  else
!   call ATMOS_GRID_CARTESC_REAL_setup -->
      ! setup map projection
      call MAPPROJECTION_setup( ATMOS_GRID_CARTESC_DOMAIN_CENTER_X, ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y )
!   <-- ATMOS_GRID_CARTESC_REAL_setup
  end if

  if (exec_model .and. scalerm_run) then

    ! setup restart
#ifdef SCALEUV
    call ADMIN_restart_setup( member = scalerm_mem )
#else
    call ADMIN_restart_setup
#endif
  end if

  ! setup time
  if (scalerm_run) then
    if (execname_ == 'SCALERM' .or. execname_ == 'DACYCLE') then
      call ADMIN_TIME_setup( setup_TimeIntegration = .true. )
    else if (execname_ == 'RMPREP ') then
      call ADMIN_TIME_setup( setup_TimeIntegration = .false. )
    end if
  end if

  if (exec_model .and. scalerm_run) then
    ! setup statistics
    call STATISTICS_setup
  end if

  if ((execname_ == 'SCALERM' .or. execname_ == 'DACYCLE') .and. scalerm_run) then
    ! setup history I/O
    call FILE_HISTORY_CARTESC_setup

    ! setup monitor I/O
    call MONITOR_CARTESC_setup( TIME_DTSEC, ATMOS_do, OCEAN_do, LAND_do, URBAN_do )

    ! setup external in
    call FILE_EXTERNAL_INPUT_CARTESC_setup
  end if

  if (exec_model .and. scalerm_run) then
    ! setup nesting grid
    call COMM_CARTESC_NEST_setup ( QA_MP, ATMOS_PHY_MP_TYPE, intercomm_parent, intercomm_child )

    ! setup common tools
    call ATMOS_HYDROSTATIC_setup
    call ATMOS_THERMODYN_setup
    call ATMOS_SATURATION_setup
  end if

  if ((execname_ == 'SCALERM' .or. execname_ == 'DACYCLE') .and. scalerm_run) then
    call BULKFLUX_setup( sqrt(DX**2+DY**2) )
  end if

  if (exec_model .and. scalerm_run) then
    ! setup variable container
    if ( ATMOS_do ) call ATMOS_vars_setup
    if ( OCEAN_do ) call OCEAN_vars_setup
    if ( LAND_do  ) call LAND_vars_setup
    if ( URBAN_do ) call URBAN_vars_setup
    if ( CPL_sw   ) call CPL_vars_setup
  end if

  if (scalerm_run) then
    if (execname_ == 'SCALERM' .or. execname_ == 'DACYCLE') then
      ! setup driver
      if ( ATMOS_do ) call ATMOS_driver_setup
      if ( OCEAN_do ) call OCEAN_driver_setup
      if ( LAND_do  ) call LAND_driver_setup
      if ( URBAN_do ) call URBAN_driver_setup

      call USER_setup
    else if (execname_ == 'RMPREP ') then
      ! setup preprocess converter
      call CONVERT_setup

      ! setup mktopo
      call MKTOPO_setup

      ! setup mkinit
      call MKINIT_setup
    end if
  end if

  if (exec_model .and. scalerm_run) then
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
!  use scale_process, only: &
!    PRC_MPIfinish
  use scale_file_cartesC, only: &
    FILE_CARTESC_cleanup
  use scale_file, only: &
    FILE_Close_All
  use scale_comm_cartesC, only: &
    COMM_cleanup
  use scale_file_history, only: &
    FILE_HISTORY_finalize
  use scale_monitor, only: &
    MONITOR_finalize
  use mod_atmos_vars, only: &
    ATMOS_sw_check => ATMOS_RESTART_CHECK, &
    ATMOS_vars_restart_check
  implicit none

  character(len=*), intent(in), optional :: execname
  integer :: ierr
  logical :: exec_model

  character(len=7) :: execname_ = 'LETKF  '

  if (present(execname)) execname_ = execname
  exec_model = .false.
  if (execname_ == 'SCALERM' .or. execname_ == 'RMPREP ' .or. execname_ == 'DACYCLE') then
    exec_model = .true.
  end if

  if (myrank_use .and. scalerm_run) then
    if (execname_ == 'SCALERM' .or. execname_ == 'DACYCLE') then
      ! check data
      if( ATMOS_sw_check ) call ATMOS_vars_restart_check

      call PROF_rapstart('Monit', 2)
      call MONITOR_finalize
      call PROF_rapend  ('Monit', 2)

      call PROF_rapstart('File', 2)

      call FILE_HISTORY_finalize
    end if

    if (exec_model) then
      ! clean up resource allocated for I/O
      call FILE_CARTESC_cleanup
    end if

    if (execname_ == 'SCALERM' .or. execname_ == 'DACYCLE') then
      call COMM_cleanup
    end if
  end if

  call FILE_Close_All

  if (myrank_use .and. scalerm_run) then
    if (execname_ == 'SCALERM' .or. execname_ == 'DACYCLE') then
      call PROF_rapend  ('File', 2)

      call PROF_rapend  ('All', 1)
    end if

    if (exec_model) then
      call PROF_rapreport
    end if
  end if

  if (myrank_use) then
    if (NUM_DOMAIN <= 1) then ! When NUM_DOMAIN >= 2, 'PRC_MPIfinish' in SCALE library will free the communicator
      call MPI_COMM_FREE(MPI_COMM_d, ierr)
    end if
    call MPI_COMM_FREE(MPI_COMM_a, ierr)
    call MPI_COMM_FREE(MPI_COMM_da, ierr)
    call MPI_COMM_FREE(MPI_COMM_u, ierr)
  end if

  ! stop MPI
!  call PRC_MPIfinish

!  ! Close logfile, configfile
!  if ( IO_L ) then
!    if( IO_FID_LOG /= IO_FID_STDOUT ) close(IO_FID_LOG)
!  endif
!  close(IO_FID_CONF)

!  if( universal_master ) write(*,*) '*** End   Launch System for SCALE-RM'

  return
end subroutine scalerm_finalize

!-------------------------------------------------------------------------------
! resume_state
!-------------------------------------------------------------------------------
subroutine resume_state(do_restart_read)
  use scale_atmos_grid_cartesC_index
  use scale_atmos_grid_cartesC, only: &
     CZ   => ATMOS_GRID_CARTESC_CZ,  &
     FZ   => ATMOS_GRID_CARTESC_FZ,  &
     FDZ  => ATMOS_GRID_CARTESC_FDZ, &
     RCDZ => ATMOS_GRID_CARTESC_RCDZ
  use scale_atmos_grid_cartesC_real, only: &
     REAL_CZ  => ATMOS_GRID_CARTESC_REAL_CZ,  &
     REAL_FZ  => ATMOS_GRID_CARTESC_REAL_FZ,  &
     REAL_PHI => ATMOS_GRID_CARTESC_REAL_PHI, &
     AREA     => ATMOS_GRID_CARTESC_REAL_AREA
  use scale_time, only: &
     TIME_NOWDAYSEC
  use mod_admin_restart, only: &
     ADMIN_restart_read
  use mod_atmos_admin, only: &
     ATMOS_do
  use mod_atmos_driver, only: &
     ATMOS_driver_calc_tendency, &
     ATMOS_driver_calc_tendency_from_sflux, &
     ATMOS_SURFACE_SET
  use mod_atmos_vars, only: &
     ATMOS_vars_calc_diagnostics, &
     ATMOS_vars_history_setpres,  &
     ATMOS_vars_history,          &
     ATMOS_vars_monitor,          &
     DENS,                        &
     POTT,                        &
     TEMP,                        &
     PRES,                        &
     QV
  use mod_atmos_bnd_driver, only: &
     ATMOS_BOUNDARY_driver_set
  use scale_atmos_refstate, only: &
     ATMOS_REFSTATE_UPDATE
  use mod_ocean_admin, only: &
     OCEAN_do
  use mod_ocean_driver, only: &
     OCEAN_driver_calc_tendency, &
     OCEAN_SURFACE_SET
  use mod_ocean_vars, only: &
     OCEAN_vars_history
  use mod_land_admin, only: &
     LAND_do
  use mod_land_driver, only: &
     LAND_driver_calc_tendency, &
     LAND_SURFACE_SET
  use mod_land_vars, only: &
     LAND_vars_history
  use mod_urban_admin, only: &
     URBAN_do
  use mod_urban_driver, only: &
     URBAN_driver_calc_tendency, &
     URBAN_SURFACE_SET
  use mod_urban_vars, only: &
     URBAN_vars_history
  use mod_cpl_admin, only: &
     CPL_sw
  use mod_user, only: &
     USER_calc_tendency

  implicit none
  logical, intent(in), optional :: do_restart_read

  logical :: do_restart_read_ = .true.
  !---------------------------------------------------------------------------

  if (present(do_restart_read)) do_restart_read_ = do_restart_read

  ! read restart data
  if (do_restart_read_) then
    call ADMIN_restart_read
  end if

  if ( ATMOS_do ) then
    ! calc diagnostics
    call ATMOS_vars_calc_diagnostics
    call ATMOS_REFSTATE_update( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                DENS(:,:,:), POTT(:,:,:), TEMP(:,:,:),  PRES(:,:,:), QV(:,:,:), & ! [IN]
                                CZ(:), FZ(:), FDZ(:), RCDZ(:), & ! [IN]
                                REAL_CZ(:,:,:), REAL_FZ(:,:,:), REAL_PHI(:,:,:), AREA(:,:),    & ! [IN]
                                TIME_NOWDAYSEC, & ! [IN]
                                force = .true.)
    call ATMOS_BOUNDARY_driver_set
    call ATMOS_vars_history_setpres
  endif

  ! setup surface condition
  if( ATMOS_do ) call ATMOS_SURFACE_SET( countup=.false. )
  if( OCEAN_do ) call OCEAN_SURFACE_SET( countup=.false. )
  if( LAND_do  ) call LAND_SURFACE_SET ( countup=.false. )
  if( URBAN_do ) call URBAN_SURFACE_SET( countup=.false. )

  ! calc tendencies
  if( ATMOS_do ) call ATMOS_driver_calc_tendency           ( force=.true. )
  if( OCEAN_do ) call OCEAN_driver_calc_tendency           ( force=.true. )
  if( LAND_do  ) call LAND_driver_calc_tendency            ( force=.true. )
  if( URBAN_do ) call URBAN_driver_calc_tendency           ( force=.true. )
  if( CPL_sw   ) call ATMOS_driver_calc_tendency_from_sflux( force=.true. )
!                 call USER_calc_tendency

  !########## History & Monitor ##########
  if( ATMOS_do ) call ATMOS_vars_history
  if( OCEAN_do ) call OCEAN_vars_history
  if( LAND_do  ) call LAND_vars_history
  if( URBAN_do ) call URBAN_vars_history

  call ATMOS_vars_monitor

  return
end subroutine resume_state

!===============================================================================
end module common_scalerm
