MODULE common_nml
!=======================================================================
!
! [PURPOSE:] Read namelist
!
! [HISTORY:]
!   November 2014   Guo-Yuan Lien     created
!
!=======================================================================
  use common, only: r_size
  use scale_stdio, only: IO_FID_CONF

  implicit none
  public

  !--- PARAM_LETKF_PRC
  integer :: NNODES = 1
  integer :: PPN = 1
  integer :: MEM_NODES = 1
  integer :: MEM_NP = 1
!  integer :: PRC_NUM_X_LETKF = 1
!  integer :: PRC_NUM_Y_LETKF = 1

  !--- PARAM_LETKF_OBS
  integer :: SLOT_START = 1
  integer :: SLOT_END = 1
  integer :: SLOT_BASE = 1
  real(r_size) :: SLOT_TINTERVAL = 3600.0d0

contains
!-----------------------------------------------------------------------
! PARAM_LETKF_PRC
!-----------------------------------------------------------------------
subroutine read_nml_letkf_prc
  use common_mpi, only: nprocs
  implicit none
  integer :: ierr
  
  namelist /PARAM_LETKF_PRC/ &
    NNODES, &
    PPN, &
    MEM_NODES, &
    MEM_NP
!    PRC_NUM_X_LETKF, &
!    PRC_NUM_Y_LETKF

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_PRC,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'xxx Not found namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist LETKF_PARAM_PRC. Check!'
    stop
  endif

!  if (nprocs /= NNODES * PPN) then
!    write(6,*) 'Number of MPI processes should be equal to NNODES * PPN.'
!    stop
!  else if (MEM_NP /= PRC_NUM_X * PRC_NUM_Y) then
!    write(6,*) 'MEM_NP should be equal to PRC_NUM_X * PRC_NUM_Y.'
!    stop
!  else if (IMAX /= nlon) then
!    write(6,*) 'IMAX should be equal to nlon.'
!    stop
!  else if (JMAX /= nlat) then
!    write(6,*) 'JMAX should be equal to nlat.'
!    stop
!  else if (KMAX /= nlev) then
!    write(6,*) 'KMAX should be equal to nlev.'
!    stop
!  end if

  return
end subroutine read_nml_letkf_prc

!-----------------------------------------------------------------------
! PARAM_LETKF_OBS
!-----------------------------------------------------------------------
subroutine read_nml_letkf_obs
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_OBS/ &
    SLOT_START, &
    SLOT_END, &
    SLOT_BASE, &
    SLOT_TINTERVAL

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_OBS,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'xxx Not found namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist LETKF_PARAM_OBS. Check!'
    stop
  endif

  return
end subroutine read_nml_letkf_obs


end module common_nml
