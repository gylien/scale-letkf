MODULE letkf_namelist
!=======================================================================
!
! [PURPOSE:] Read namelist
!
! [HISTORY:]
!   November 2014   Guo-Yuan Lien     created
!
!=======================================================================
!  use common, only :: r_size
  use scale_stdio, only: IO_FID_CONF

  implicit none
  public

  !--- PARAM_LETKF_PRC
  integer :: NNODES = 1
  integer :: PPN = 1
  integer :: MEM_NODES = 1
  integer :: MEM_NP = 1
  integer :: PRC_NUM_X_LETKF = 1
  integer :: PRC_NUM_Y_LETKF = 1

contains
!-----------------------------------------------------------------------
! Read namelist for obsope
!-----------------------------------------------------------------------
subroutine read_nml_obsope
  implicit none

  call read_nml_letkf_prc

  return
end subroutine read_nml_obsope

!-----------------------------------------------------------------------
! PARAM_LETKF_PRC
!-----------------------------------------------------------------------
subroutine read_nml_letkf_prc
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_PRC/ &
    NNODES, &
    PPN, &
    MEM_NODES, &
    MEM_NP, &
    PRC_NUM_X_LETKF, &
    PRC_NUM_Y_LETKF

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=LETKF_PARAM_PRC,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'xxx Not found namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist LETKF_PARAM_PRC. Check!'
    stop
  endif

  if (nprocs /= NNODES * PPN) then
    write(6,*) 'Number of MPI processes should be equal to NNODES * PPN.'
    stop
  else if (MEM_NP /= PRC_NUM_X * PRC_NUM_Y) then
    write(6,*) 'MEM_NP should be equal to PRC_NUM_X * PRC_NUM_Y.'
    stop
  end if

  return
end subroutine read_nml_letkf_prc

end module letkf_namelist
