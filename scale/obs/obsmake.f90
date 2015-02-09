PROGRAM obsmake
!=======================================================================
!
! [PURPOSE:] Main program of observation operator
!
! [HISTORY:]
!   November 2014  Guo-Yuan Lien     Created
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

  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(11) :: stdoutf='NOUT-000000'

  type(obs_info) :: obs(nobsfiles)
  real(r_size) :: radarlon, radarlat, radarz

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  CALL initialize_mpi
  rtimer00 = MPI_WTIME()

  WRITE(stdoutf(6:11), '(I6.6)') myrank
  WRITE(6,'(3A,I6.6)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I6.6,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf

!-----------------------------------------------------------------------

  CALL set_common_scale(1)
  CALL set_common_obs_scale

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  rtimer = MPI_WTIME()
  WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

!-----------------------------------------------------------------------
! Read observations
!-----------------------------------------------------------------------

  if (scale_IO_mygroup == 1) then ! only run at the first group
    call read_obs_all(obs, radarlon, radarlat, radarz)
  end if

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  rtimer = MPI_WTIME()
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_OBS):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

!-----------------------------------------------------------------------
! Generate observations
!-----------------------------------------------------------------------

  if (scale_IO_mygroup == 1) then ! only run at the first group
    CALL obsmake_cal(obs, radarlon, radarlat, radarz)
  end if

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  rtimer = MPI_WTIME()
  WRITE(6,'(A,2F10.2)') '### TIMER(OBSMAKE):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  CALL unset_common_scale
  CALL finalize_mpi

  STOP
END PROGRAM obsmake
