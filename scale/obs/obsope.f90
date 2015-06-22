PROGRAM obsope
!=======================================================================
!
! [PURPOSE:] Main program of observation operator
!
! [HISTORY:]
!   11/12/2014 Guo-Yuan Lien     Created
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
  REAL(r_dble) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(11) :: stdoutf='NOUT-000000'
  CHARACTER(11) :: timer_fmt='(A30,F10.2)'

  type(obs_info) :: obs(nobsfiles)
  real(r_size) :: radarlon, radarlat, radarz

  character(len=6400) :: cmd, icmd
  character(len=10) :: myranks

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  if (COMMAND_ARGUMENT_COUNT() >= 5) then
    call chdir('run/obsope')
  end if

  CALL initialize_mpi
  rtimer00 = MPI_WTIME()

  WRITE(stdoutf(6:11), '(I6.6)') myrank
  WRITE(6,'(3A,I6.6)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I6.6,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf

!-----------------------------------------------------------------------
! Pre-processing scripts
!-----------------------------------------------------------------------

  if (COMMAND_ARGUMENT_COUNT() >= 5) then
    write(6,'(A)') 'Run pre-processing scripts'
    cmd = 'bash'
    call get_command_argument(2, icmd)
    cmd = trim(cmd) // ' ' // trim(icmd)
    call get_command_argument(3, icmd)
    cmd = trim(cmd) // ' ' // trim(icmd) // '_pre'
    call get_command_argument(4, icmd)
    cmd = trim(cmd) // ' ' // trim(icmd)
    call get_command_argument(5, icmd)
    cmd = trim(cmd) // ' ' // trim(icmd)
    write (myranks, '(I10)') myrank
    cmd = trim(cmd) // ' ' // trim(myranks)

    WRITE(6,'(A,I6.6,3A)') 'MYRANK ',myrank,' is running a script: [', trim(cmd), ']'
    call system(trim(cmd))
  end if

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  rtimer = MPI_WTIME()
  WRITE(6,timer_fmt) '### TIMER(PRE_SCRIPT):',rtimer-rtimer00
  rtimer00=rtimer

!-----------------------------------------------------------------------

  call set_common_conf

  call read_nml_letkf_obs
  call read_nml_letkf_obserr
  call read_nml_letkf_obs_radar

  call set_mem_node_proc(MEMBER+1,NNODES,PPN,MEM_NODES,MEM_NP)

  if (myrank_use) then

    call set_scalelib

    call set_common_scale
    CALL set_common_mpi_scale
    call set_common_obs_scale

    CALL MPI_BARRIER(MPI_COMM_a,ierr)
    rtimer = MPI_WTIME()
    WRITE(6,timer_fmt) '### TIMER(INITIALIZE):',rtimer-rtimer00
    rtimer00=rtimer

!-----------------------------------------------------------------------
! Read observations
!-----------------------------------------------------------------------

    if (myrank_mem_use) then
      call read_obs_all(obs, radarlon, radarlat, radarz)
    end if

    CALL MPI_BARRIER(MPI_COMM_a,ierr)
    rtimer = MPI_WTIME()
    WRITE(6,timer_fmt) '### TIMER(READ_OBS):',rtimer-rtimer00
    rtimer00=rtimer

!-----------------------------------------------------------------------
! Observation operator
!-----------------------------------------------------------------------

    if (myrank_mem_use) then
      call obsope_cal(obs, radarlon, radarlat, radarz)
    end if

    CALL MPI_BARRIER(MPI_COMM_a,ierr)
    rtimer = MPI_WTIME()
    WRITE(6,timer_fmt) '### TIMER(OBS_OPERATOR):',rtimer-rtimer00
    rtimer00=rtimer

    CALL unset_common_mpi_scale

    call unset_scalelib

  else ! [ myrank_use ]

    write (6, '(A,I6.6,A)') 'MYRANK=',myrank,': This process is not used!'

  end if ! [ myrank_use ]

!-----------------------------------------------------------------------

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  rtimer = MPI_WTIME()
  WRITE(6,timer_fmt) '### TIMER(FINALIZE):',rtimer-rtimer00
  rtimer00=rtimer

!-----------------------------------------------------------------------
! Post-processing scripts
!-----------------------------------------------------------------------

  if (COMMAND_ARGUMENT_COUNT() >= 5) then
    write(6,'(A)') 'Run post-processing scripts'
    cmd = 'bash'
    call get_command_argument(2, icmd)
    cmd = trim(cmd) // ' ' // trim(icmd)
    call get_command_argument(3, icmd)
    cmd = trim(cmd) // ' ' // trim(icmd) // '_post'
    call get_command_argument(4, icmd)
    cmd = trim(cmd) // ' ' // trim(icmd)
    call get_command_argument(5, icmd)
    cmd = trim(cmd) // ' ' // trim(icmd)
!    write (myranks, '(I10)') myrank
    cmd = trim(cmd) // ' ' // trim(myranks)

    WRITE(6,'(A,I6.6,3A)') 'MYRANK ',myrank,' is running a script: [', trim(cmd), ']'
    call system(trim(cmd))
  end if

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  rtimer = MPI_WTIME()
  WRITE(6,timer_fmt) '### TIMER(POST_SCRIPT):',rtimer-rtimer00
  rtimer00=rtimer

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  CALL finalize_mpi

  STOP
END PROGRAM obsope
