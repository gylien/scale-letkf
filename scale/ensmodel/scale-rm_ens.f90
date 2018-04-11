program scale_rm_ens
!=======================================================================
!
! [PURPOSE:] Ensemble forecasts with SCALE-RM
!
!=======================================================================
!!$USE OMP_LIB
!  use mpi
!  use common
  use common_mpi, only: &
    nprocs, &
    myrank
  use common_nml
  use common_scale, only: &
    set_common_conf
  use common_scalerm
  use common_mpi_scale, only: &
    mpi_timer, &
    myrank_use, &
    set_mem_node_proc

  use scale_stdio, only: &
    IO_L, &
    IO_FID_LOG
  use scale_prof, only: &
    PROF_setprefx, &
    PROF_rapstart, &
    PROF_rapend
  use scale_process, only: &
    PRC_MPIstart, &
    PRC_mpi_alive, &
    PRC_UNIVERSAL_setup, &
    PRC_MPIfinish, &
    PRC_UNIVERSAL_myrank
  use scale_history, only: &
    HIST_write
  use scale_monitor, only: &
    MONIT_write
  use mod_admin_restart, only: &
#ifdef SCALEUV
    ADMIN_restart_write, &
    ADMIN_restart_write_additional
#else
    ADMIN_restart_write
#endif
  use mod_admin_time, only: &
    ADMIN_TIME_checkstate, &
    ADMIN_TIME_advance, &
    TIME_DOATMOS_step, &
    TIME_DOLAND_step, &
    TIME_DOURBAN_step, &
    TIME_DOOCEAN_step, &
    TIME_DOresume, &
    TIME_DOend
  use mod_atmos_admin, only: &
    ATMOS_do
  use mod_atmos_driver, only: &
    ATMOS_driver, &
    ATMOS_driver_finalize
  use mod_ocean_admin, only: &
    OCEAN_do
  use mod_ocean_driver, only: &
    OCEAN_driver
  use mod_land_admin, only: &
    LAND_do
  use mod_land_driver, only: &
    LAND_driver
  use mod_urban_admin, only: &
    URBAN_do
  use mod_urban_driver, only: &
    URBAN_driver
  use mod_user, only: &
    USER_step
  implicit none

  character(7) :: stdoutf = '-000000'
  integer :: member_iter_cmd = -1

  integer :: universal_comm
  integer :: universal_nprocs
  logical :: universal_master
  integer :: universal_myrank

  character(len=6400) :: cmd1, cmd2, icmd
  character(len=10) :: myranks
  integer :: iarg

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  ! start MPI
  call PRC_MPIstart( universal_comm ) ! [OUT]

  PRC_mpi_alive = .true.

  call mpi_timer('', 1)

  call PRC_UNIVERSAL_setup( universal_comm,   & ! [IN]
                            universal_nprocs, & ! [OUT]
                            universal_master  ) ! [OUT]
  universal_myrank = PRC_UNIVERSAL_myrank
  nprocs = universal_nprocs
  myrank = universal_myrank

  WRITE(6,'(A,I6.6,A,I6.6)') 'Hello from MYRANK ',universal_myrank,'/',universal_nprocs-1

  if (command_argument_count() >= 2) then
    call get_command_argument(2, icmd)
    if (trim(icmd) /= '' .and. trim(icmd) /= '-') then
      read (icmd, '(I)') member_iter_cmd
    end if
  end if

  if (command_argument_count() >= 3) then
    call get_command_argument(3, icmd)
    if (trim(icmd) /= '' .and. trim(icmd) /= '-') then
      WRITE(stdoutf(2:7), '(I6.6)') universal_myrank
!      WRITE(6,'(3A,I6.6)') 'STDOUT goes to ',trim(icmd)//stdoutf,' for MYRANK ', universal_myrank
      OPEN(6,FILE=trim(icmd)//stdoutf)
      WRITE(6,'(A,I6.6,2A)') 'MYRANK=',universal_myrank,', STDOUTF=',trim(icmd)//stdoutf
    end if
  end if

  if (command_argument_count() >= 4) then
    write (myranks, '(I10)') universal_myrank
    call get_command_argument(4, icmd)
    cmd1 = 'bash ' // trim(icmd) // ' ensfcst_1' // ' ' // trim(myranks)
    cmd2 = 'bash ' // trim(icmd) // ' ensfcst_2' // ' ' // trim(myranks)
    do iarg = 5, command_argument_count()
      call get_command_argument(iarg, icmd)
      cmd1 = trim(cmd1) // ' ' // trim(icmd)
      cmd2 = trim(cmd2) // ' ' // trim(icmd)
    end do
  end if

!-----------------------------------------------------------------------
! Pre-processing scripts
!-----------------------------------------------------------------------

  if (command_argument_count() >= 4) then
    write (6,'(A)') 'Run pre-processing scripts'
    write (6,'(A,I6.6,3A)') 'MYRANK ',universal_myrank,' is running a script: [', trim(cmd1), ']'
    call system(trim(cmd1))
  end if

  call mpi_timer('PRE_SCRIPT', 1, barrier=universal_comm)

!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------

  call set_common_conf(universal_nprocs)

  if (member_iter_cmd >= 0) then
    MEMBER_ITER = member_iter_cmd
    write (6, '(A,I4,A)') '[Info] Reset MEMBER_ITER = ', MEMBER_ITER, ' based on the commannd line input'
  end if

  call set_mem_node_proc(MEMBER_RUN)

  call scalerm_setup('SCALERM')

  call mpi_timer('INITIALIZE', 1, barrier=universal_comm)

  if (myrank_use .and. scalerm_run) then

!-----------------------------------------------------------------------
! Main loop
!-----------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START TIMESTEP ++++++'
    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_Loop', 0)

    do

      ! report current time
      call ADMIN_TIME_checkstate

      if ( TIME_DOresume ) then
        ! resume state from restart files
        call resume_state

        ! history&monitor file output
        call MONIT_write('MAIN')
        call HIST_write ! if needed
      end if

      ! time advance
      call ADMIN_TIME_advance

      ! user-defined procedure
      call USER_step

      ! change to next state
      if( OCEAN_do .AND. TIME_DOOCEAN_step ) call OCEAN_driver
      if( LAND_do  .AND. TIME_DOLAND_step  ) call LAND_driver
      if( URBAN_do .AND. TIME_DOURBAN_step ) call URBAN_driver
      if( ATMOS_do .AND. TIME_DOATMOS_step ) call ATMOS_driver

      ! history&monitor file output
      call MONIT_write('MAIN')
      call HIST_write

      ! restart output
      call ADMIN_restart_write
#ifdef SCALEUV
      call ADMIN_restart_write_additional
#endif

      if( TIME_DOend ) exit

      if( IO_L ) call flush(IO_FID_LOG)

    end do

    call PROF_rapend('Main_Loop', 0)

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END TIMESTEP ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    call PROF_setprefx('FIN')

    call PROF_rapstart('All', 1)

    if( ATMOS_do ) call ATMOS_driver_finalize

  end if ! [ myrank_use .and. scalerm_run ]

  call scalerm_finalize('SCALERM')

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_WORLD)

!-----------------------------------------------------------------------
! Post-processing scripts
!-----------------------------------------------------------------------

  if (command_argument_count() >= 4) then
    write (6,'(A)') 'Run post-processing scripts'
    write (6,'(A,I6.6,3A)') 'MYRANK ',universal_myrank,' is running a script: [', trim(cmd2), ']'
    call system(trim(cmd2))
  end if

  call mpi_timer('POST_SCRIPT', 1, barrier=universal_comm)

!-----------------------------------------------------------------------

  call PRC_MPIfinish

  stop
!=======================================================================
end program scale_rm_ens
