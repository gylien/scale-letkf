program pdbash2
!=======================================================================
!
! [PURPOSE:] Execute parallel distributed bash script
!
!=======================================================================
  USE common
  USE common_mpi
  USE common_scale, only: set_mem_node_proc
!  USE common_mpi_scale

  implicit none

  integer, parameter :: maxlen = 6400
!  character(len=maxlen) :: cmd
!  character(len=10) :: myranks

  integer :: MEMBER = 1
  integer :: NNODES = 1
  integer :: PPN = 1
  integer :: MEM_NODES = 1
  integer :: MEM_NP = 1



  integer :: i, cnt, ilen, ierr
  CHARACTER(len=32) :: arg

  cnt = command_argument_count()
  write (*,*) 'number of command arguments = ', cnt

  do i = 1, cnt
     call get_command_argument(i, arg, ilen, ierr)
     if (ierr /= 0) then
       write (*,*) 'get_command_argument failed: status = ', ierr, ' arg = ', arg
       stop
     end if
     write (*,*) 'command arg ', i, ' = ', arg(1:ilen)
  end do

  write (*,*) 'command line processed'



  CALL initialize_mpi

!  CALL set_common_scale(-1)



  ! setup standard I/O
!  call IO_setup( MODELNAME, .false.)


!  call read_nml_letkf
!  call read_nml_letkf_prc
!  call read_nml_letkf_obs !!!!!!!!!!!!!!!!!!!!!! move outside of subroutine????
!  call read_nml_letkf_obserr !!!!!!!!!!!!!!!!!!! move outside of subroutine????
!  call read_nml_letkf_obs_radar !!!!!!!!!!!!!!!! move outside of subroutine????

!  if (nprocs /= NNODES * PPN) then
!    write(6,'(A,I10)') 'Number of MPI processes = ', nprocs
!    write(6,'(A,I10)') 'NNODES = ', NNODES
!    write(6,'(A,I10)') 'PPN    = ', PPN
!    write(6,'(A)') 'Number of MPI processes should be equal to NNODES * PPN.'
!    stop
!  end if

  !
  ! Set up node and process distribution
  !
  call set_mem_node_proc(MEMBER,NNODES,PPN,MEM_NODES,MEM_NP)

!! print process distribution!!!


!  if (scale_IO_mygroup > 0) then

!    call set_scalelib

!    CALL set_common_mpi_scale

!    CALL unset_common_mpi_scale

!    call unset_scalelib

!  else
!    write (6, '(A,I6.6,A)') 'MYRANK=',myrank,': This process is not used!'
!  end if

  CALL finalize_mpi









!  call MPI_INIT(ierr)
!  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
!  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

!  call get_command(cmd)
!  pos = index(cmd, ' ')
!  cmd = cmd(pos+1:maxlen)
!  pos = index(cmd, ' ')

!  write (myranks, '(I10)') myrank

!!     print *, 'bash ' // trim(cmd(1:pos-1)) // trim(myranks) // trim(cmd(pos:maxlen))
!  call system('bash ' // trim(cmd(1:pos-1)) // trim(myranks) // trim(cmd(pos:maxlen)))

!  call MPI_FINALIZE(ierr)

end program pdbash2
