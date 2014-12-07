PROGRAM letkf
!=======================================================================
!
! [PURPOSE:] Main program of LETKF
!
! [HISTORY:]
!   01/16/2009 Takemasa Miyoshi  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
!  USE common_obs_scale

  USE common_letkf, only: nbv

  use common_nml

  USE letkf_obs
  USE letkf_tools

  IMPLICIT NONE
  REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: anal3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: anal2d(:,:,:)
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(11) :: stdoutf='NOUT-000000'


!  TYPE(obs_info) :: obs
!  TYPE(obs_da_value) :: obsval



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
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '  LOCAL ENSEMBLE TRANSFORM KALMAN FILTERING  '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '   LL      EEEEEE  TTTTTT  KK  KK  FFFFFF    '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(6,'(A)') '   LL      EEEEE     TT    KKK     FFFFF     '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(6,'(A)') '   LLLLLL  EEEEEE    TT    KK  KK  FF        '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '             WITHOUT LOCAL PATCH             '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '          Coded by Takemasa Miyoshi          '
  WRITE(6,'(A)') '  Based on Ott et al (2004) and Hunt (2005)  '
  WRITE(6,'(A)') '  Tested by Miyoshi and Yamane (2006)        '
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '              LETKF PARAMETERS               '
  WRITE(6,'(A)') ' ------------------------------------------- '
  WRITE(6,'(A,I15)')   '  nbv          :',nbv
  WRITE(6,'(A,F15.2)') '  sigma_obs    :',sigma_obs
  WRITE(6,'(A,F15.2)') '  sigma_obsv   :',sigma_obsv
  WRITE(6,'(A,F15.2)') '  sigma_obst   :',sigma_obst
  WRITE(6,'(A)') '============================================='
  CALL set_common_scale

!-----------------------------------------------------------------------

  ! setup standard I/O
  call IO_setup( MODELNAME )

  call read_nml_letkf

  if (nprocs /= NNODES * PPN) then
    write(6,*) 'Number of MPI processes should be equal to NNODES * PPN.'
    stop
  end if

!
!  CALL set_mem_node_proc(nbv+1,NNODES,PPN,MEM_NODES,MEM_NP)
!

  CALL set_common_mpi_scale(nbv+1,NNODES,PPN,MEM_NODES,MEM_NP)


!!!
  if (scale_IO_group_n >= 1) then
!!!
  ALLOCATE(gues3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(gues2d(nij1,nbv,nv2d))
  ALLOCATE(anal3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(anal2d(nij1,nbv,nv2d))
!!!
  end if
!!!

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

!-----------------------------------------------------------------------

!!!
  if (scale_IO_group_n >= 1) then
!!!

  CALL get_nobs(obsfile,8,obs%nobs)
  WRITE(6,'(A,I9,A)') 'TOTAL: ', obs%nobs, ' OBSERVATIONS'

  CALL obs_info_allocate(obs)

  CALL read_obs(obsfile,obs)

!!!
  end if
!!!

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_OBS):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

!!-----------------------------------------------------------------------
!! Observations
!!-----------------------------------------------------------------------
!  !
!  ! CONVENTIONAL OBS
!  !
  CALL set_letkf_obs
!

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_OBS):',rtimer,rtimer-rtimer00
  rtimer00=rtimer



!  if (scale_IO_group_n >= 1) then
!!  write (6,*) obsda%idx
!!  write (6,*) obsda%val(3)
!!  write (6,*) obsda%ensval(:,3)
!!  write (6,*) obsda%qc(3)
!!  write (6,*) obsda%ri(3)
!!  write (6,*) obsda%rj(3)
!!  write (6,*) obsda2%idx
!!  write (6,*) obsda2%val(3)
!!  write (6,*) obsda2%ensval(:,3)
!!  write (6,*) obsda2%qc(3)
!!  write (6,*) obsda2%ri(3)
!!  write (6,*) obsda2%rj(3)
!  write (6,*) obsda2%ri
!  write (6,*) obsda2%rj
!  end if


!-----------------------------------------------------------------------
! First guess ensemble
!-----------------------------------------------------------------------
  !
  ! READ GUES
  !
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!!!
  if (scale_IO_group_n >= 1) then
!!!
  call read_ens_mpi('gues',gues3d,gues2d)
!!!
  end if
!!!



!  if (scale_IO_group_n >= 1) then
!  write (6,*) gues3d(20,:,3,iv3d_t)
!!  write (6,*) gues2d
!  end if


  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_GUES):',rtimer,rtimer-rtimer00
  rtimer00=rtimer


  !
  ! WRITE ENS MEAN and SPRD
  !
!!!
  if (scale_IO_group_n >= 1) then
!!!
  CALL write_ensmspr_mpi('gues',gues3d,gues2d)
!!!
  end if
!!!
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(GUES_MEAN):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!!-----------------------------------------------------------------------
!! Data Assimilation
!!-----------------------------------------------------------------------
!  !
!  ! LETKF
!  !

  anal3d = gues3d
  anal2d = gues2d

!  CALL das_letkf(gues3d,gues2d,anal3d,anal2d)
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(DAS_LETKF):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Analysis ensemble
!-----------------------------------------------------------------------
  !
  ! WRITE ANAL
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!!
  if (scale_IO_group_n >= 1) then
!!!
  CALL write_ens_mpi('anal',anal3d,anal2d)
!!!
  end if
!!!

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(WRITE_ANAL):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
  !
  ! WRITE ENS MEAN and SPRD
  !
!!!
  if (scale_IO_group_n >= 1) then
!!!
  CALL write_ensmspr_mpi('anal',anal3d,anal2d)
!!!
  end if
!!!
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(ANAL_MEAN):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!!-----------------------------------------------------------------------
!! Monitor
!!-----------------------------------------------------------------------
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  CALL monit_obs
!!
!  CALL CPU_TIME(rtimer)
!  WRITE(6,'(A,2F10.2)') '### TIMER(MONIT_MEAN):',rtimer,rtimer-rtimer00
!  rtimer00=rtimer

  CALL unset_common_mpi_scale


!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

  STOP
END PROGRAM letkf
