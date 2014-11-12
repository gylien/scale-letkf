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

  IMPLICIT NONE
  INTEGER :: ierr
  CHARACTER(8) :: stdoutf='NOUT-000'

  INTEGER,ALLOCATABLE :: ranks(:)
  INTEGER :: MPI_G, MPI_G_WORLD, MPI_COMM_WORLD_ALL
  INTEGER :: myrank2, myrank3, myrank4
  INTEGER :: nprocs2, nprocs3, nprocs4

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL initialize_mpi
!
  WRITE(stdoutf(6:8), '(I3.3)') myrank
  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
!  WRITE(6,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!

  ALLOCATE(ranks(4))




  call MPI_Comm_group(MPI_COMM_WORLD, MPI_G_WORLD, ierr)

  MPI_COMM_WORLD_ALL = MPI_COMM_WORLD 

  if (myrank < 4) then
    ranks(1) = 0
    ranks(2) = 1
    ranks(3) = 2
    ranks(4) = 3
    call MPI_GROUP_INCL(MPI_G_WORLD,4,ranks,MPI_G,ierr)
  else
    ranks(1) = 4
    ranks(2) = 5
    ranks(3) = 6
    ranks(4) = 7
    call MPI_GROUP_INCL(MPI_G_WORLD,4,ranks,MPI_G,ierr)
  end if

!  call MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_G,MPI_C,ierr)
  call MPI_COMM_CREATE(MPI_COMM_WORLD_ALL,MPI_G,MPI_COMM_WORLD,ierr)

!  CALL MPI_COMM_SIZE(MPI_C,nprocs2,ierr)
!  CALL MPI_COMM_RANK(MPI_C,myrank2,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs2,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank2,ierr)
!  call MPI_Group_rank (MPI_G, myrank3, ierr)


  call MPI_Comm_free(MPI_COMM_WORLD, ierr)
!  MPI_COMM_WORLD = MPI_COMM_WORLD_ALL

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs4,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank4,ierr)

  write(*,*) nprocs, nprocs2, nprocs4, myrank, myrank2, myrank4

!printf("rank= %d newrank= %d recvbuf= %d\n",rank,new_rank,recvbuf); 

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

  STOP
END PROGRAM letkf
