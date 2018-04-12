program scale_rm_pp_ens
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
  use scale_atmos_grid_cartesC_real, only: &
    ATMOS_GRID_CARTESC_REAL_update_Z
  use mod_convert, only: &
    CONVERT
  use mod_mktopo, only: &
    MKTOPO
  use mod_mkinit, only: &
    MKINIT
  implicit none

  character(7) :: stdoutf = '-000000'

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
    if (trim(icmd) /= '') then
      WRITE(stdoutf(2:7), '(I6.6)') universal_myrank
!      WRITE(6,'(3A,I6.6)') 'STDOUT goes to ',trim(icmd)//stdoutf,' for MYRANK ', universal_myrank
      OPEN(6,FILE=trim(icmd)//stdoutf)
      WRITE(6,'(A,I6.6,2A)') 'MYRANK=',universal_myrank,', STDOUTF=',trim(icmd)//stdoutf
    end if
  end if

  if (command_argument_count() >= 3) then
    write (myranks, '(I10)') universal_myrank
    call get_command_argument(3, icmd)
    cmd1 = 'bash ' // trim(icmd) // ' enspp_1' // ' ' // trim(myranks)
    cmd2 = 'bash ' // trim(icmd) // ' enspp_2' // ' ' // trim(myranks)
    do iarg = 4, command_argument_count()
      call get_command_argument(iarg, icmd)
      cmd1 = trim(cmd1) // ' ' // trim(icmd)
      cmd2 = trim(cmd2) // ' ' // trim(icmd)
    end do
  end if

!-----------------------------------------------------------------------
! Pre-processing scripts
!-----------------------------------------------------------------------

  if (command_argument_count() >= 3) then
    write (6,'(A)') 'Run pre-processing scripts'
    write (6,'(A,I6.6,3A)') 'MYRANK ',universal_myrank,' is running a script: [', trim(cmd1), ']'
    call system(trim(cmd1))
  end if

  call mpi_timer('PRE_SCRIPT', 1, barrier=universal_comm)

!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------

  call set_common_conf(universal_nprocs)

  call set_mem_node_proc(MEMBER_RUN)

  call scalerm_setup('RMPREP')

  call mpi_timer('INITIALIZE', 1, barrier=universal_comm)

  if (myrank_use .and. scalerm_run) then

!-----------------------------------------------------------------------
! Main
!-----------------------------------------------------------------------

    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_prep', 0)

    ! execute preprocess
    call PROF_rapstart('Convert', 1)
    call CONVERT
    call PROF_rapend  ('Convert', 1)

    ! execute mktopo
    call PROF_rapstart('MkTopo', 1)
    call MKTOPO
    call PROF_rapend  ('MkTopo', 1)

    ! re-setup
    call ATMOS_GRID_CARTESC_REAL_update_Z

    ! execute mkinit
    call PROF_rapstart('MkInit', 1)
    call MKINIT
    call PROF_rapend  ('MkInit', 1)

    call PROF_rapend('Main_prep', 0)

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  end if ! [ myrank_use .and. scalerm_run ]

  call scalerm_finalize('RMPREP')

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_WORLD)

!-----------------------------------------------------------------------
! Post-processing scripts
!-----------------------------------------------------------------------

  if (command_argument_count() >= 3) then
    write (6,'(A)') 'Run post-processing scripts'
    write (6,'(A,I6.6,3A)') 'MYRANK ',universal_myrank,' is running a script: [', trim(cmd2), ']'
    call system(trim(cmd2))
  end if

  call mpi_timer('POST_SCRIPT', 1, barrier=universal_comm)

!-----------------------------------------------------------------------

  call PRC_MPIfinish

  stop
!=======================================================================
end program scale_rm_pp_ens
