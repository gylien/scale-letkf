PROGRAM obsmake
!=======================================================================
!
! [PURPOSE:] Main program of synthetic observation generator
!
! [HISTORY:]
!   November 2014  Guo-Yuan Lien     Created
!   .............  See git history for the following revisions
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
  USE common_obs_scale
  USE common_nml
  USE obsope_tools
  IMPLICIT NONE

  character(len=7) :: stdoutf = '-000000'
  character(len=6400) :: icmd

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  call initialize_mpi_scale
  call mpi_timer('', 1)

  if (command_argument_count() >= 2) then
    call get_command_argument(2, icmd)
    if (trim(icmd) /= '') then
      write (stdoutf(2:7), '(I6.6)') myrank
!      write (6,'(3A,I6.6)') 'STDOUT goes to ', trim(icmd)//stdoutf, ' for MYRANK ', myrank
      open (6, file=trim(icmd)//stdoutf)
      write (6,'(A,I6.6,2A)') 'MYRANK=', myrank, ', STDOUTF=', trim(icmd)//stdoutf
    end if
  end if

!-----------------------------------------------------------------------

  call set_common_conf

  call set_mem_node_proc(1)
  call set_scalelib('OBSMAKE')

  if (myrank_use) then

    call set_common_scale
    call set_common_mpi_scale
    call set_common_obs_scale

    call mpi_timer('INITIALIZE', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Read observations
!-----------------------------------------------------------------------

    allocate(obs(OBS_IN_NUM))
    call read_obs_all(obs)

    call mpi_timer('READ_OBS', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Generate observations
!-----------------------------------------------------------------------

    call obsmake_cal(obs)

    call mpi_timer('OBSMAKE', 1, barrier=MPI_COMM_a)

    deallocate(obs)

    call unset_common_mpi_scale

  end if ! [ myrank_use ]

  call unset_scalelib

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_WORLD)

  call finalize_mpi_scale

  STOP
END PROGRAM obsmake
