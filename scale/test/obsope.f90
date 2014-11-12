!-------------------------------------------------------------------------------
!> Program SCALE-LES file I/O sample
!!
!! @par Description
!!          File I/O sample
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-10-08 (H.Yashiro)  [add] New
!!
!<
!-------------------------------------------------------------------------------
program scaleles_fileiosample
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  USE common
  USE common_mpi

  use mod_ioapi, only: &
     IOAPI

  IMPLICIT NONE
  INTEGER :: ierr
  CHARACTER(8) :: stdoutf='NOUT-000'

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL initialize_mpi
!
!  WRITE(stdoutf(6:8), '(I3.3)') myrank
!  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
!  OPEN(6,FILE=stdoutf)

  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

  call IOAPI





  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi


  stop
  !=============================================================================
end program scaleles_fileiosample
