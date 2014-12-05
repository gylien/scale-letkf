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

  TYPE(obs_info) :: obs
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(8) :: stdoutf='NOUT-000000'

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi
!
  WRITE(stdoutf(6:11), '(I6.6)') myrank
  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!
  CALL set_common_scale

!-----------------------------------------------------------------------

  ! setup standard I/O
  call IO_setup( MODELNAME )

  call read_nml_obsmake

  if (nprocs /= NNODES * PPN) then
    write(6,*) 'Number of MPI processes should be equal to NNODES * PPN.'
    stop
  else if (nprocs /= MEM_NP) then
    write(6,*) 'Number of MPI processes should be equal to MEM_NP.'
    stop
  end if

  CALL set_mem_node_proc(1,NNODES,PPN,MEM_NODES,MEM_NP)

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

!-----------------------------------------------------------------------

  CALL get_nobs(obsfile,8,obs%nobs)
  WRITE(6,'(A,I9,A)') 'TOTAL: ', obs%nobs, ' OBSERVATIONS'

  CALL obs_info_allocate(obs)

  CALL read_obs(obsfile,obs)

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_OBS):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

!-----------------------------------------------------------------------

  if (scale_IO_group_n == 1) then ! only run at the first group
    CALL obsmake_cal(obs)
  end if

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(OBSMAKE):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

  STOP
END PROGRAM obsmake
