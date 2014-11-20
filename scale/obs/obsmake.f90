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

!  use common_letkf, only: nbv

  use common_nml

  use obsope_tools

!  use scale_process, only: &
!       prc_myrank, &
!       prc_myrank_world

!  use scale_grid_index, only: &
!    KHALO, IHALO, JHALO


  IMPLICIT NONE
!  REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:)
!  REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)
!  REAL(r_size),ALLOCATABLE :: anal3d(:,:,:,:)
!  REAL(r_size),ALLOCATABLE :: anal2d(:,:,:)


  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(8) :: stdoutf='NOUT-000000'



  TYPE(obs_info) :: obs
!  TYPE(obs_ensval) :: obsval


!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi
!
  WRITE(stdoutf(6:11), '(I6.6)') myrank
!  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
!  OPEN(6,FILE=stdoutf)
!  WRITE(6,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!
  CALL set_common_scale

!  ALLOCATE(gues3d(nij1,nlev,nbv,nv3d))
!  ALLOCATE(gues2d(nij1,nbv,nv2d))
!  ALLOCATE(anal3d(nij1,nlev,nbv,nv3d))
!  ALLOCATE(anal2d(nij1,nbv,nv2d))
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

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


!-----------------------------------------------------------------------


  CALL get_nobs(obsfile,8,obs%nobs)
  WRITE(6,'(A,I9,A)') 'TOTAL: ', obs%nobs, ' OBSERVATIONS'

  !IF(obs%nobs == 0) ...

  CALL obs_info_allocate(obs)
  CALL read_obs(obsfile,obs)

!  print *, obs%nobs
!  print *, obs%lat, obs%lon, obs%lev, obs%elm



  call set_common_mpi_scale(1,NNODES,PPN,MEM_NODES,MEM_NP)


  call obsmake_cal(obs)




!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

  STOP
END PROGRAM obsmake
