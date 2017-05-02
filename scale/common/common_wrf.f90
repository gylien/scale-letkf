MODULE common_wrf
!===============================================================================
!
! [PURPOSE:] Common Information for WRF
!
!===============================================================================
!$USE OMP_LIB
  USE common
  use common_nml

  use scale_precision, only: RP, SP
  use scale_prof

  IMPLICIT NONE
  PUBLIC
!-------------------------------------------------------------------------------
! General parameters
!-------------------------------------------------------------------------------

  INTEGER,PARAMETER :: vname_max = 10

  ! 
  !--- 3D, 2D state variables (in SCALE restart files)
  ! 
  INTEGER,PARAMETER :: nv3d = 12
  INTEGER,PARAMETER :: nv2d = 1
  INTEGER,PARAMETER :: iv3d_theta_pert=4 !-- in restart files
  INTEGER,PARAMETER :: iv3d_p_pert=5     !
  INTEGER,PARAMETER :: iv3d_ph_pert=6    !
  INTEGER,PARAMETER :: iv2d_mu_pert=1    !
  INTEGER,PARAMETER :: iv3d_u=1          !-- for LETKF
  INTEGER,PARAMETER :: iv3d_v=2          !
  INTEGER,PARAMETER :: iv3d_w=3          !
  INTEGER,PARAMETER :: iv3d_t=4          !
  INTEGER,PARAMETER :: iv3d_p=5          !
  INTEGER,PARAMETER :: iv3d_ph=6         !
  INTEGER,PARAMETER :: iv3d_q=7          !
  INTEGER,PARAMETER :: iv3d_qc=8         !
  INTEGER,PARAMETER :: iv3d_qr=9         !
  INTEGER,PARAMETER :: iv3d_qi=10        !
  INTEGER,PARAMETER :: iv3d_qs=11        !
  INTEGER,PARAMETER :: iv3d_qg=12        !
  INTEGER,PARAMETER :: iv2d_mu=1         !
  CHARACTER(vname_max),PARAMETER :: v3d_name(nv3d) = &
     (/'U', 'V', 'W', 'T', 'P', 'PH', &
       'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP'/)
  CHARACTER(vname_max),PARAMETER :: v2d_name(nv2d) = &
     (/'MU'/)

  ! 
  !--- 3D, 2D diagnostic variables (in SCALE history files)
  ! 
  INTEGER,PARAMETER :: nv3dd=12
  INTEGER,PARAMETER :: nv2dd=6
  INTEGER,PARAMETER :: iv3dd_theta_pert=4 !-- in history files
  INTEGER,PARAMETER :: iv3dd_p_pert=5     !
  INTEGER,PARAMETER :: iv3dd_ph_pert=6    !
  INTEGER,PARAMETER :: iv3dd_u=1
  INTEGER,PARAMETER :: iv3dd_v=2
  INTEGER,PARAMETER :: iv3dd_w=3
  INTEGER,PARAMETER :: iv3dd_t=4
  INTEGER,PARAMETER :: iv3dd_p=5
  INTEGER,PARAMETER :: iv3dd_hgt=6
  INTEGER,PARAMETER :: iv3dd_q=7
  INTEGER,PARAMETER :: iv3dd_qc=8
  INTEGER,PARAMETER :: iv3dd_qr=9
  INTEGER,PARAMETER :: iv3dd_qi=10
  INTEGER,PARAMETER :: iv3dd_qs=11
  INTEGER,PARAMETER :: iv3dd_qg=12
  INTEGER,PARAMETER :: iv3dd_rh=-999 ! cannot use
  INTEGER,PARAMETER :: iv2dd_topo=1
  INTEGER,PARAMETER :: iv2dd_ps=2
  INTEGER,PARAMETER :: iv2dd_u10m=3
  INTEGER,PARAMETER :: iv2dd_v10m=4
  INTEGER,PARAMETER :: iv2dd_t2m=5
  INTEGER,PARAMETER :: iv2dd_q2m=6
  INTEGER,PARAMETER :: iv2dd_rain=-999 ! connot use
  CHARACTER(vname_max),PARAMETER :: v3dd_name(nv3dd) = &
     (/'U', 'V', 'W', 'T', 'P', 'PH', &
       'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP'/)
  CHARACTER(vname_max),PARAMETER :: v2dd_name(nv2dd) = &
     (/'HGT', 'PSFC', 'U10', 'V10', 'T2', 'Q2'/)

  ! 
  !--- Variables for model coordinates
  ! 
  REAL(r_size),ALLOCATABLE,SAVE :: p_base3d(:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: ph_base3d(:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: ph_pert3d(:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: mu_base2d(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lon2d(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lat2d(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: topo2d(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: eta_stag(:)
  REAL(r_size),ALLOCATABLE,SAVE :: eta(:)
  REAL(r_size),SAVE :: t00 = undef
  REAL(r_size),SAVE :: p00 = undef
  REAL(r_size),SAVE :: ptop = undef
  CHARACTER(vname_max),PARAMETER :: p_base3d_name = 'PB'
  CHARACTER(vname_max),PARAMETER :: ph_base3d_name = 'PHB'
  CHARACTER(vname_max),PARAMETER :: ph_pert3d_name = 'PH'
  CHARACTER(vname_max),PARAMETER :: mu_base2d_name = 'MUB'
  CHARACTER(vname_max),PARAMETER :: lon2d_name = 'XLONG'
  CHARACTER(vname_max),PARAMETER :: lat2d_name = 'XLAT'
  CHARACTER(vname_max),PARAMETER :: topo2d_name = 'HGT'
  CHARACTER(vname_max),PARAMETER :: eta_stag_name = 'ZNW'
  CHARACTER(vname_max),PARAMETER :: eta_name = 'ZNU'
  CHARACTER(vname_max),PARAMETER :: t00_name = 'T00'
  CHARACTER(vname_max),PARAMETER :: p00_name = 'P00'
  CHARACTER(vname_max),PARAMETER :: ptop_name = 'P_TOP'

  ! 
  !--- grid settings
  ! 
  INTEGER,SAVE :: nlon  ! # grids in I-direction [subdomain]
  INTEGER,SAVE :: nlat  ! # grids in J-direction [subdomain]
  INTEGER,SAVE :: nlev  ! # grids in K-direction
  INTEGER,SAVE :: nlong ! # grids in I-direction [global domain]
  INTEGER,SAVE :: nlatg ! # grids in J-direction [global domain]
  INTEGER,SAVE :: nlonh ! # grids in I-direction [subdomain with halo]
  INTEGER,SAVE :: nlath ! # grids in J-direction [subdomain with halo]
  INTEGER,SAVE :: nlevh ! # grids in K-direction [with halo]
  INTEGER,SAVE :: nij0
  INTEGER,SAVE :: nlevall
  INTEGER,SAVE :: nlevalld
  INTEGER,SAVE :: ngpv
  INTEGER,SAVE :: ngpvd

CONTAINS

!-------------------------------------------------------------------------------
! [File I/O] Read WRF restart files <PnetCDF>
!-------------------------------------------------------------------------------
#ifdef PNETCDF
subroutine read_restart_par(filename,v3dg,v2dg,comm,trans,verify_p)
  use scale_process, only: &
      PRC_myrank
  use scale_rm_process, only: &
      PRC_2Drank
  use scale_grid_index, only: &
      IMAX, JMAX, KMAX
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL, MPI_COMM_WORLD
  use common_mpi, only: myrank
  use pnetcdf
  implicit none

  character(*),intent(in) :: filename
  real(r_size),intent(out) :: v3dg(nlev,nlon,nlat,nv3d)
  real(r_size),intent(out) :: v2dg(nlon,nlat,nv2d)
  integer,intent(in) :: comm
  logical,intent(in),optional :: trans
  logical,intent(in),optional :: verify_p
  logical :: trans_
  logical :: verify_p_

  integer :: i,j,k,iv3d,iv2d
  real(SP) :: var3D(nlon,nlat,nlev)
  real(SP) :: var2D(nlon,nlat)
  real(r_size), allocatable :: var_p(:,:,:)
  real(r_size) :: kappa, kapdiv, Rd_wrf, RvRd_wrf, p00_kappa_inv
  real(r_size) :: tmprho, tmptheta_p00_kappa_inv

  integer :: err, ncid, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(4), count(4)

  trans_ = .true.
  if (present(trans)) trans_ = trans
  verify_p_ = .false.
  if (present(verify_p)) verify_p_ = (trans_ .and. verify_p)

  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename)

  err = nfmpi_open(comm, trim(filename), NF_NOWRITE, MPI_INFO_NULL, ncid)

  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//nfmpi_strerror(err)

  start(4) = 1
  count(1) = IMAX
  count(2) = JMAX
  count(3) = KMAX
  count(4) = 1

  ! 3D variables
  !-------------
  do iv3d = 1, nv3d
    select case (iv3d)
    case (iv3d_u)
      start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 2
      start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
      start(3) = 1
    case (iv3d_v)
      start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
      start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 2
      start(3) = 1
    case (iv3d_w, iv3d_ph_pert)
      start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
      start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
      start(3) = 2
    case default
      start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
      start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
      start(3) = 1
    end select

    write(6,'(1x,A,A15,A,8I6)') '*** Read 3D var: ', trim(v3d_name(iv3d)), ' >> PnetCDF start(4), count(4) =', start, count

    err = nfmpi_inq_varid(ncid, trim(v3d_name(iv3d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start, count, var3D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

    forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg(k,i,j,iv3d) = real(var3D(i,j,k), r_size)
  end do

  start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  start(3) = 1
  count(1) = IMAX
  count(2) = JMAX
  count(3) = 1

  ! 2D variables
  !-------------
  do iv2d = 1, nv2d
    write(6,'(1x,A,A15,A,6I6)') '*** Read 2D var: ', trim(v2d_name(iv2d)), ' >> PnetCDF start(3), count(3) =', start(1:3), count(1:3)

    err = nfmpi_inq_varid(ncid, trim(v2d_name(iv2d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:3), count(1:3), var2D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

    v2dg(:,:,iv2d) = real(var2D(:,:), r_size)
  end do

  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  ! transform
  !----------
  if (t00 == undef .or. p00 == undef .or. &
      (.not. allocated(p_base3d)) .or. (.not. allocated(ph_base3d)) .or. (.not. allocated(mu_base2d)) .or. (.not. allocated(eta_stag))) then
    write (6, '(A)') '[Error] p_base3d, ph_base3d, mu_base2d, eta_stag, t00, p00 should be initialized before.'
    stop 1
  end if

  if (trans_) then
    kappa = 2.0_r_size / 7.0_r_size
    kapdiv = 1.0_r_size / (1.0_r_size - kappa)
    Rd_wrf = 287.0_r_size
    RvRd_wrf = 461.6_r_size / Rd_wrf
    p00_kappa_inv = p00 ** (- kappa)
    if (verify_p_) then
      allocate (var_p(nlev,nlon,nlat))
    end if

    do j = 1, nlat
      do i = 1, nlon
        v2dg(i,j,iv2d_mu) = mu_base2d(i,j) + v2dg(i,j,iv2d_mu_pert)

        do k = 1, nlev
          v3dg(k,i,j,iv3d_ph) = ph_base3d(k+1,i,j) + v3dg(k,i,j,iv3d_ph_pert)
          v3dg(k,i,j,iv3d_p) = p_base3d(k,i,j) + v3dg(k,i,j,iv3d_p_pert)
          tmptheta_p00_kappa_inv = (t00 + v3dg(k,i,j,iv3d_theta_pert)) * p00_kappa_inv
          v3dg(k,i,j,iv3d_t) = tmptheta_p00_kappa_inv * v3dg(k,i,j,iv3d_p) ** kappa

          if (verify_p_) then
            if (k == 1) then
              tmprho = v2dg(i,j,iv2d_mu) * (eta_stag(k) - eta_stag(k+1)) / (v3dg(k,i,j,iv3d_ph) - ph_base3d(k,i,j))
            else
              tmprho = v2dg(i,j,iv2d_mu) * (eta_stag(k) - eta_stag(k+1)) / (v3dg(k,i,j,iv3d_ph) - v3dg(k-1,i,j,iv3d_ph))
            end if
            var_p(k,i,j) = (tmprho * Rd_wrf * tmptheta_p00_kappa_inv * (1.0_r_size + v3dg(k,i,j,iv3d_q) * RvRd_wrf)) ** kapdiv
          end if
        end do
      end do
    end do

    if (verify_p_) then
      if (maxval(abs(v3dg(:,:,:,iv3d_p) - var_p)) > 20d0) then
        write (6, '(A,F15.7)') '[Error] Pressure calculation is incorrect! -- maxdiff(p) = ', &
                               maxval(abs(v3dg(:,:,:,iv3d_p) - var_p))
        stop
      else
        write (6, '(A)') 'VERIFY_COORD: Pressure calculation is verified!'
      end if
      deallocate (var_p)
    end if
  end if ! [ trans_ ]

  return
end subroutine read_restart_par
#endif

!-------------------------------------------------------------------------------
! [File I/O] Write WRF restart files <PnetCDF>
!-------------------------------------------------------------------------------
! mode  1: Use PH analysis; output P analysis as well as it is (default)
!       2: Use P analysis; output re-constructed PH field
!-------------------------------------------------------------------------------
#ifdef PNETCDF
subroutine write_restart_par(filename,v3dg,v2dg,comm,trans,mode)
  use scale_process, only: &
      PRC_myrank
  use scale_rm_process, only: &
      PRC_2Drank
  use scale_grid_index, only: &
      IMAX, JMAX, KMAX
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL, MPI_COMM_WORLD
  use common_mpi, only: myrank
  use pnetcdf
  implicit none

  character(*),intent(in) :: filename
  real(r_size),intent(inout) :: v3dg(nlev,nlon,nlat,nv3d)
  real(r_size),intent(inout) :: v2dg(nlon,nlat,nv2d)
  integer,intent(in) :: comm
  logical,intent(in),optional :: trans
  integer,intent(in),optional :: mode
  logical :: trans_
  integer :: mode_

  integer :: i,j,k,iv3d,iv2d
  real(SP) :: var3D(nlon,nlat,nlev)
  real(SP) :: var2D(nlon,nlat)
  real(r_size) :: kappa, Rd_wrf, RvRd_wrf, p00_kappa
  real(r_size) :: tmprho, tmptheta

  integer :: err, ncid, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(4), count(4)

  trans_ = .true.
  if (present(trans)) trans_ = trans
  mode_ = 1
  if (present(mode)) mode_ = mode

  ! transform
  !----------
  if (t00 == undef .or. p00 == undef .or. &
      (.not. allocated(p_base3d)) .or. (.not. allocated(ph_base3d)) .or. (.not. allocated(mu_base2d)) .or. (.not. allocated(eta_stag))) then
    write (6, '(A)') '[Error] p_base3d, ph_base3d, mu_base2d, eta_stag, t00, p00 should be initialized before.'
    stop 1
  end if

  if (trans_) then
    kappa = 2.0_r_size / 7.0_r_size
    Rd_wrf = 287.0_r_size
    RvRd_wrf = 461.6_r_size / Rd_wrf
    p00_kappa = p00 ** kappa

    if (mode_ == 2) then
      write (6,'(1x,A)') '*** Re-construct the PH analysis field based on P, T, Q, MU analyses'
    end if

    do j = 1, nlat
      do i = 1, nlon
        do k = 1, nlev
          tmptheta = v3dg(k,i,j,iv3d_t) * p00_kappa * v3dg(k,i,j,iv3d_p) ** (- kappa)
          tmprho = v3dg(k,i,j,iv3d_p) / (v3dg(k,i,j,iv3d_t) * (1.0_r_size + v3dg(k,i,j,iv3d_q) * RvRd_wrf) * Rd_wrf)
!          tmprho = (p00_kappa * v3dg(k,i,j,iv3d_p) ** (1.0_r_size - kappa)) / (tmptheta * (1.0_r_size + v3dg(k,i,j,iv3d_q) * RvRd_wrf) * Rd_wrf)
          v3dg(k,i,j,iv3d_theta_pert) = tmptheta - t00
          if (mode_ == 2) then
            if (k == 1) then
              v3dg(k,i,j,iv3d_ph) = ph_base3d(k,i,j) + v2dg(i,j,iv2d_mu) * (eta_stag(k) - eta_stag(k+1)) / tmprho
            else
              v3dg(k,i,j,iv3d_ph) = v3dg(k-1,i,j,iv3d_ph) + v2dg(i,j,iv2d_mu) * (eta_stag(k) - eta_stag(k+1)) / tmprho
            end if
          end if
          v3dg(k,i,j,iv3d_p_pert) = v3dg(k,i,j,iv3d_p) - p_base3d(k,i,j)
          v3dg(k,i,j,iv3d_ph_pert) = v3dg(k,i,j,iv3d_ph) - ph_base3d(k+1,i,j)
        end do

        v2dg(i,j,iv2d_mu_pert) = v2dg(i,j,iv2d_mu) - mu_base2d(i,j)
      end do
    end do
  end if ! [ trans_ ]

  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is writing a file ',trim(filename)

  err = nfmpi_open(comm, trim(filename), NF_WRITE, MPI_INFO_NULL, ncid)

  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//nfmpi_strerror(err)

  start(4) = 1
  count(1) = IMAX
  count(2) = JMAX
  count(3) = KMAX
  count(4) = 1

  ! 3D variables
  !-------------
  do iv3d = 1, nv3d
    select case (iv3d)
    case (iv3d_u)
      start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 2
      start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
      start(3) = 1
    case (iv3d_v)
      start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
      start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 2
      start(3) = 1
    case (iv3d_w, iv3d_ph_pert)
      start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
      start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
      start(3) = 2
    case default
      start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
      start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
      start(3) = 1
    end select

    write(6,'(1x,A,A15,A,8I6)') '*** Write 3D var: ', trim(v3d_name(iv3d)), ' >> PnetCDF start(4), count(4) =', start, count

    forall (i=1:nlon, j=1:nlat, k=1:nlev) var3D(i,j,k) = real(v3dg(k,i,j,iv3d), SP)

    err = nfmpi_inq_varid(ncid, trim(v3d_name(iv3d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iput_vara_real(ncid, varid, start, count, var3D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iput_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)
  end do

  start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  start(3) = 1
  count(1) = IMAX
  count(2) = JMAX
  count(3) = 1

  ! 2D variables
  !-------------
  do iv2d = 1, nv2d
    write(6,'(1x,A,A15,A,6I6)') '*** Write 2D var: ', trim(v2d_name(iv2d)), ' >> PnetCDF start(3), count(3) =', start(1:3), count(1:3)

    var2D(:,:) = real(v2dg(:,:,iv2d), SP)

    err = nfmpi_inq_varid(ncid, trim(v2d_name(iv2d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iput_vara_real(ncid, varid, start(1:3), count(1:3), var2D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iput_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)
  end do

  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  return
end subroutine write_restart_par
#endif

!-------------------------------------------------------------------------------
! [File I/O] Read WRF restart files for model coordinates <PnetCDF>
!-------------------------------------------------------------------------------
#ifdef PNETCDF
subroutine read_coor_par(filename,comm,t00,p00,ptop,p_base,ph_base,ph_pert,mu_base,lon,lat,topo,eta_stag,eta)
  use netcdf, only: NF90_NOWRITE
  use scale_process, only: &
    PRC_myrank
  use scale_rm_process, only: &
    PRC_2Drank
  use scale_grid_index, only: &
    IMAX, JMAX, KMAX
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL
  use common_mpi, only: myrank
  use pnetcdf
  implicit none

  character(*), intent(in) :: filename
  integer,      intent(in) :: comm
  real(r_size), intent(out), optional :: t00
  real(r_size), intent(out), optional :: p00
  real(r_size), intent(out), optional :: ptop
  real(r_size), intent(out), optional :: p_base(nlev,nlon,nlat)
  real(r_size), intent(out), optional :: ph_base(nlev+1,nlon,nlat)
  real(r_size), intent(out), optional :: ph_pert(nlev+1,nlon,nlat)
  real(r_size), intent(out), optional :: mu_base(nlon,nlat)
  real(r_size), intent(out), optional :: lon(nlon,nlat)
  real(r_size), intent(out), optional :: lat(nlon,nlat)
  real(r_size), intent(out), optional :: topo(nlon,nlat)
  real(r_size), intent(out), optional :: eta_stag(nlev+1)
  real(r_size), intent(out), optional :: eta(nlev)

  integer :: i,j,k
  real(SP) :: var3D(nlon,nlat,nlev)
  real(SP) :: var3D_stag(nlon,nlat,nlev+1)
  real(SP) :: var2D(nlon,nlat)
  real(SP) :: var1DZ(nlev)
  real(SP) :: var1DZ_stag(nlev+1)
  real(SP) :: var0D(1)

  integer :: err, ncid, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(4), count(4)

  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename)

  err = nfmpi_open(comm, trim(filename), NF_NOWRITE, MPI_INFO_NULL, ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//nfmpi_strerror(err)

  ! 3D variables
  !-------------
  start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  start(3) = 1
  start(4) = 1
  count(1) = IMAX
  count(2) = JMAX
  count(4) = 1

  !--- p_base
  if (present(p_base)) then
    count(3) = KMAX

    write(6,'(1x,A,A15,A,8I6)') '*** Read 3D var: ', trim(p_base3d_name), ' >> PnetCDF start(4), count(4) =', start, count

    err = nfmpi_inq_varid(ncid, trim(p_base3d_name), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start, count, var3D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    if (.not. allocated(p_base)) allocate(p_base(nlev,nlon,nlat))
    forall (i=1:nlon, j=1:nlat, k=1:nlev) p_base(k,i,j) = real(var3D(i,j,k), r_size)
  end if

  !--- ph_base
  if (present(ph_base)) then
    count(3) = KMAX + 1

    write(6,'(1x,A,A15,A,8I6)') '*** Read 3D var: ', trim(ph_base3d_name), ' >> PnetCDF start(4), count(4) =', start, count

    err = nfmpi_inq_varid(ncid, trim(ph_base3d_name), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start, count, var3D_stag, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    if (.not. allocated(ph_base)) allocate(ph_base(nlev+1,nlon,nlat))
    forall (i=1:nlon, j=1:nlat, k=1:nlev+1) ph_base(k,i,j) = real(var3D_stag(i,j,k), r_size)
  end if

  !--- ph_pert
  if (present(ph_pert)) then
    count(3) = KMAX + 1

    write(6,'(1x,A,A15,A,8I6)') '*** Read 3D var: ', trim(ph_pert3d_name), ' >> PnetCDF start(4), count(4) =', start, count

    err = nfmpi_inq_varid(ncid, trim(ph_pert3d_name), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start, count, var3D_stag, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    if (.not. allocated(ph_pert)) allocate(ph_pert(nlev+1,nlon,nlat))
    forall (i=1:nlon, j=1:nlat, k=1:nlev+1) ph_pert(k,i,j) = real(var3D_stag(i,j,k), r_size)
  end if

  ! 2D variables
  !-------------
  start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  start(3) = 1
  count(1) = IMAX
  count(2) = JMAX
  count(3) = 1

  !--- mu_base
  if (present(mu_base)) then
    write(6,'(1x,A,A15,A,6I6)') '*** Read 2D var: ', trim(mu_base2d_name), ' >> PnetCDF start(3), count(3) =', start(1:3), count(1:3)

    err = nfmpi_inq_varid(ncid, trim(mu_base2d_name), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:3), count(1:3), var2D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    if (.not. allocated(mu_base)) allocate(mu_base(nlon,nlat))
    mu_base(:,:) = real(var2D(:,:), r_size)
  end if

  !--- lon
  if (present(lon)) then
    write(6,'(1x,A,A15,A,6I6)') '*** Read 2D var: ', trim(lon2d_name), ' >> PnetCDF start(3), count(3) =', start(1:3), count(1:3)

    err = nfmpi_inq_varid(ncid, trim(lon2d_name), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:3), count(1:3), var2D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    if (.not. allocated(lon)) allocate(lon(nlon,nlat))
    lon(:,:) = real(var2D(:,:), r_size)
  end if

  !--- lat
  if (present(lat)) then
    write(6,'(1x,A,A15,A,6I6)') '*** Read 2D var: ', trim(lat2d_name), ' >> PnetCDF start(3), count(3) =', start(1:3), count(1:3)

    err = nfmpi_inq_varid(ncid, trim(lat2d_name), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:3), count(1:3), var2D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    if (.not. allocated(lat)) allocate(lat(nlon,nlat))
    lat(:,:) = real(var2D(:,:), r_size)
  end if

  !--- topo
  if (present(topo)) then
    write(6,'(1x,A,A15,A,6I6)') '*** Read 2D var: ', trim(topo2d_name), ' >> PnetCDF start(3), count(3) =', start(1:3), count(1:3)

    err = nfmpi_inq_varid(ncid, trim(topo2d_name), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:3), count(1:3), var2D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    if (.not. allocated(topo)) allocate(topo(nlon,nlat))
    topo(:,:) = real(var2D(:,:), r_size)
  end if

  ! 1D variables
  !-------------
  start(1) = 1
  start(2) = 1
  count(2) = 1

  !--- eta
  if (present(eta)) then
    count(1) = KMAX

    write(6,'(1x,A,A15,A,4I6)') '*** Read 2D var: ', trim(eta_name), ' >> PnetCDF start(2), count(2) =', start(1:2), count(1:2)

    err = nfmpi_inq_varid(ncid, trim(eta_name), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:2), count(1:2), var1DZ, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    if (.not. allocated(eta)) allocate(eta(nlev))
    eta(:) = real(var1DZ(:), r_size)
  end if

  !--- eta_stag
  if (present(eta_stag)) then
    count(1) = KMAX + 1

    write(6,'(1x,A,A15,A,4I6)') '*** Read 2D var: ', trim(eta_stag_name), ' >> PnetCDF start(2), count(2) =', start(1:2), count(1:2)

    err = nfmpi_inq_varid(ncid, trim(eta_stag_name), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:2), count(1:2), var1DZ_stag, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    if (.not. allocated(eta_stag)) allocate(eta_stag(nlev+1))
    eta_stag(:) = real(var1DZ_stag(:), r_size)
  end if


  ! Scalar variables
  !-----------------
  start(1) = 1
  count(1) = 1

  !--- t00
  if (present(t00)) then
    t00 = 300.0_r_size
  end if

  !--- p00
  if (present(p00)) then
    p00 = 100000.0_r_size
  end if

  !--- ptop
  if (present(ptop)) then
    write(6,'(1x,A,A15,A,2I6)') '*** Read 2D var: ', trim(ptop_name), ' >> PnetCDF start(1), count(1) =', start(1:1), count(1:1)

    err = nfmpi_inq_varid(ncid, trim(ptop_name), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:1), count(1:1), var0D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

    ptop = real(var0D(1), r_size)
  end if

  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  return
end subroutine read_coor_par
#endif

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE history files <PnetCDF>
!-------------------------------------------------------------------------------
! If v3dg_state and v2dg_state present, return the state variables as well
!-------------------------------------------------------------------------------
#ifdef PNETCDF
subroutine read_history_par(filename,step,v3dg,v2dg,comm,trans,v3dg_state,v2dg_state)
  use scale_process, only: &
      PRC_myrank
  use scale_rm_process, only: &
    PRC_2Drank
  use scale_grid_index, only: &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA, &
      IMAX, JMAX, KMAX
  use scale_comm, only: &
      COMM_vars8, &
      COMM_wait
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL, MPI_COMM_WORLD
  use common_mpi, only: myrank
  use pnetcdf
  implicit none

  character(*),intent(in) :: filename
  integer,intent(in) :: step
  real(r_size),intent(out) :: v3dg(nlevh,nlonh,nlath,nv3dd)
  real(r_size),intent(out) :: v2dg(nlonh,nlath,nv2dd)
  integer,intent(in) :: comm
  logical,intent(in),optional :: trans
  real(r_size),intent(out),optional :: v3dg_state(nlev,nlon,nlat,nv3d)
  real(r_size),intent(out),optional :: v2dg_state(nlon,nlat,nv2d)
  logical :: trans_
  logical :: return_state

  integer :: i,j,k,iv3d,iv2d
  real(SP) :: var3D(nlon,nlat,nlev)
  real(SP) :: var2D(nlon,nlat)
  real(SP) :: var_u(nlon+1,nlat,nlev)
  real(SP) :: var_v(nlon,nlat+1,nlev)
  real(SP) :: var_w(nlon,nlat,nlev+1)
  real(SP) :: var_ph(nlon,nlat,nlev+1)
  real(r_size) :: kappa, p00_kappa_inv

  integer :: err, ncid, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(4), count(4)

  trans_ = .true.
  if (present(trans)) trans_ = trans
  return_state = .false.
  if (present(v3dg_state) .and. present(v2dg_state)) return_state = .true.

  write (6,'(A,I6.6,3A,I5)') 'MYRANK ',myrank,' is reading a file ',trim(filename),', t =',step

  if (return_state) then
    write (6,'(1x,A)') '*** Also return the state variables from this history file'
  end if

  err = nfmpi_open(comm, trim(filename), NF_NOWRITE, MPI_INFO_NULL, ncid)

  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//nfmpi_strerror(err)

  ! 3D variables
  !-------------
  start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  start(3) = 1
  start(4) = step
  count(4) = 1

  do iv3d = 1, nv3dd
    select case (iv3d)
    case (iv3dd_u)
      count(1) = IMAX + 1
      count(2) = JMAX
      count(3) = KMAX

      write(6,'(1x,A,A15,A,8I6)') '*** Read 3D var: ', trim(v3dd_name(iv3d)), ' >> PnetCDF start(4), count(4) =', start, count

      err = nfmpi_inq_varid(ncid, trim(v3dd_name(iv3d)), varid)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
      err = nfmpi_iget_vara_real(ncid, varid, start, count, var_u, req)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
      err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

      forall (i=1:nlon+1, j=1:nlat, k=1:nlev) v3dg(k+KHALO,i+IHALO-1,j+JHALO,iv3d) = real(var_u(i,j,k), r_size)
    case (iv3dd_v)
      count(1) = IMAX
      count(2) = JMAX + 1
      count(3) = KMAX

      write(6,'(1x,A,A15,A,8I6)') '*** Read 3D var: ', trim(v3dd_name(iv3d)), ' >> PnetCDF start(4), count(4) =', start, count

      err = nfmpi_inq_varid(ncid, trim(v3dd_name(iv3d)), varid)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
      err = nfmpi_iget_vara_real(ncid, varid, start, count, var_v, req)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
      err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

      forall (i=1:nlon, j=1:nlat+1, k=1:nlev) v3dg(k+KHALO,i+IHALO,j+JHALO-1,iv3d) = real(var_v(i,j,k), r_size)
    case (iv3dd_w)
      count(1) = IMAX
      count(2) = JMAX
      count(3) = KMAX + 1

      write(6,'(1x,A,A15,A,8I6)') '*** Read 3D var: ', trim(v3dd_name(iv3d)), ' >> PnetCDF start(4), count(4) =', start, count

      err = nfmpi_inq_varid(ncid, trim(v3dd_name(iv3d)), varid)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
      err = nfmpi_iget_vara_real(ncid, varid, start, count, var_w, req)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
      err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

      forall (i=1:nlon, j=1:nlat, k=1:nlev+1) v3dg(k+KHALO-1,i+IHALO,j+JHALO,iv3d) = real(var_w(i,j,k), r_size)
    case (iv3dd_ph_pert)
      count(1) = IMAX
      count(2) = JMAX
      count(3) = KMAX + 1

      write(6,'(1x,A,A15,A,8I6)') '*** Read 3D var: ', trim(v3dd_name(iv3d)), ' >> PnetCDF start(4), count(4) =', start, count

      err = nfmpi_inq_varid(ncid, trim(v3dd_name(iv3d)), varid)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
      err = nfmpi_iget_vara_real(ncid, varid, start, count, var_ph, req)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
      err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

      forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg(k+KHALO,i+IHALO,j+JHALO,iv3d) = 0.5_r_size * real(var_ph(i,j,k) + var_ph(i,j,k+1), r_size)
    case default
      count(1) = IMAX
      count(2) = JMAX
      count(3) = KMAX

      write(6,'(1x,A,A15,A,8I6)') '*** Read 3D var: ', trim(v3dd_name(iv3d)), ' >> PnetCDF start(4), count(4) =', start, count

      err = nfmpi_inq_varid(ncid, trim(v3dd_name(iv3d)), varid)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
      err = nfmpi_iget_vara_real(ncid, varid, start, count, var3D, req)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
      err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
      if ( err .NE. NF_NOERR ) &
         write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

      forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg(k+KHALO,i+IHALO,j+JHALO,iv3d) = real(var3D(i,j,k), r_size)
    end select
  end do

  ! 2D variables
  !-------------
  start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  start(3) = step
  count(1) = IMAX
  count(2) = JMAX
  count(3) = 1

  do iv2d = 1, nv2dd
    write(6,'(1x,A,A15,A,6I6)') '*** Read 2D var: ', trim(v2dd_name(iv2d)), ' >> PnetCDF start(3), count(3) =', start(1:3), count(1:3)

    err = nfmpi_inq_varid(ncid, trim(v2dd_name(iv2d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:3), count(1:3), var2D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

    v2dg(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2d) = real(var2D(:,:), r_size)
  end do

  if (return_state) then
    !--- MU
    write(6,'(1x,A,A15,A,6I6)') '*** Read 2D var: ', trim(v2d_name(iv2d_mu)), ' >> PnetCDF start(3), count(3) =', start(1:3), count(1:3)

    err = nfmpi_inq_varid(ncid, trim(v2d_name(iv2d_mu)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    err = nfmpi_iget_vara_real(ncid, varid, start(1:3), count(1:3), var2D, req)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)
    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

    v2dg_state(:,:,iv2d_mu) = real(var2D(:,:), r_size)
  end if

  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  ! additionally return state variables
  !------------------------------------
  if (return_state) then
    v3dg_state(:,:,:,iv3d_u ) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_u )
    v3dg_state(:,:,:,iv3d_v ) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_v )
    v3dg_state(:,:,:,iv3d_w ) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_w )
    v3dg_state(:,:,:,iv3d_q ) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_q )
    v3dg_state(:,:,:,iv3d_qc) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qc)
    v3dg_state(:,:,:,iv3d_qr) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qr)
    v3dg_state(:,:,:,iv3d_qi) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qi)
    v3dg_state(:,:,:,iv3d_qs) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qs)
    v3dg_state(:,:,:,iv3d_qg) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qg)
  end if

  ! transform
  !----------
  if (t00 == undef .or. p00 == undef .or. &
      (.not. allocated(p_base3d)) .or. (.not. allocated(ph_base3d)) .or. (.not. allocated(mu_base2d))) then
    write (6, '(A)') '[Error] p_base3d, ph_base3d, mu_base2d, t00, p00 should be initialized before.'
    stop 1
  end if

  if (trans_) then
    kappa = 2.0_r_size / 7.0_r_size
    p00_kappa_inv = p00 ** (- kappa)

    do j = 1+JHALO, nlat+JHALO
      do i = 1+IHALO, nlon+IHALO
        do k = 1+KHALO, nlev+KHALO
          v3dg(k,i,j,iv3dd_p) = p_base3d(k-KHALO,i-IHALO,j-JHALO) + v3dg(k,i,j,iv3dd_p_pert)
          v3dg(k,i,j,iv3dd_t) = (t00 + v3dg(k,i,j,iv3dd_theta_pert)) * p00_kappa_inv * v3dg(k,i,j,iv3dd_p) ** kappa
          v3dg(k,i,j,iv3dd_hgt) = (0.5_r_size * (ph_base3d(k-KHALO,i-IHALO,j-JHALO) + ph_base3d(k-KHALO+1,i-IHALO,j-JHALO)) + &
                                   v3dg(k,i,j,iv3dd_ph_pert)) / gg
        end do
      end do
    end do

    if (return_state) then
      v3dg_state(:,:,:,iv3d_t) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_t)
      v3dg_state(:,:,:,iv3d_p) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_p)
      forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg_state(k,i,j,iv3d_ph) = ph_base3d(k+1,i,j) + real(var_ph(i,j,k+1), r_size)

      do j = 1, nlat
        do i = 1, nlon
          v2dg_state(i,j,iv2d_mu) = mu_base2d(i,j) + v2dg_state(i,j,iv2d_mu_pert)
        end do
      end do
    end if
  else ! [ trans_ ]
    if (return_state) then
      v3dg_state(:,:,:,iv3d_theta_pert) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_theta_pert)
      v3dg_state(:,:,:,iv3d_p_pert    ) = v3dg(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_p_pert    )
      forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg_state(k,i,j,iv3d_ph_pert) = real(var_ph(i,j,k+1), r_size)
    end if
  end if ! [ trans_ ]

  ! communicate halo
  !-----------------
  do iv3d = 1, nv3dd
    call COMM_vars8( v3dg(:,:,:,iv3d), iv3d )
  end do
  do iv3d = 1, nv3dd
    call COMM_wait ( v3dg(:,:,:,iv3d), iv3d )
  end do

  do iv3d = 1, nv3dd
    select case (iv3d)
    case (iv3dd_w)
      do j = 1, nlath
        do i = 1, nlonh
          v3dg(   1:KS-2,i,j,iv3d) = v3dg(KS-1,i,j,iv3d)
          v3dg(KE+1:KA,  i,j,iv3d) = v3dg(KE,i,j,iv3d)
        end do
      end do
    case default
      do j = 1, nlath
        do i = 1, nlonh
          v3dg(   1:KS-1,i,j,iv3d) = v3dg(KS,i,j,iv3d)
          v3dg(KE+1:KA,  i,j,iv3d) = v3dg(KE,i,j,iv3d)
        end do
      end do
    end select
  end do

  do iv2d = 1, nv2dd
    call COMM_vars8( v2dg(:,:,iv2d), iv2d )
  end do
  do iv2d = 1, nv2dd
    call COMM_wait ( v2dg(:,:,iv2d), iv2d )
  end do

  return
end subroutine read_history_par
#endif

!-------------------------------------------------------------------------------
! Transform the LETKF state variables to the variables in SCALE history files
! (with HALO), so that they can be used for observation operator calculation
!-------------------------------------------------------------------------------
! [INPUT]
!   v3dg, v2dg   : 3D, 2D state variables
!   topo         : topography
! [OUTPUT]
!   v3dgh, v2dgh : 3D, 2D SCALE history variables
!-------------------------------------------------------------------------------
subroutine state_to_history(v3dg, v2dg, ph_base, eta_stag, ptop, v3dgh, v2dgh, mode)
  use scale_grid_index, only: &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA
  use scale_grid, only: &
      GRID_CZ, &
      GRID_FZ
  use scale_comm, only: &
      COMM_vars8, &
      COMM_wait
  implicit none

  real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP), intent(in) :: v2dg(nlon,nlat,nv2d)
  real(r_size), intent(in) :: ph_base(nlev+1,nlon,nlat)
  real(r_size), intent(in) :: eta_stag(nlev+1)
  real(r_size), intent(in) :: ptop
  real(r_size), intent(out) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  real(r_size), intent(out) :: v2dgh(nlonh,nlath,nv2dd)
  integer, intent(in), optional :: mode
  integer :: mode_

  integer :: i, j, k, iv3d, iv2d
  real(r_size) :: Rd_wrf, RvRd_wrf, tmprho, moist_to_dry
  real(r_size) :: tmpph(nlev)

  mode_ = 1
  if (present(mode)) mode_ = mode

  ! Variables that can be directly copied
  !---------------------------------------------------------

  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_u) = v3dg(:,:,:,iv3d_u)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_v) = v3dg(:,:,:,iv3d_v)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_w) = v3dg(:,:,:,iv3d_w)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_t) = v3dg(:,:,:,iv3d_t)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_p) = v3dg(:,:,:,iv3d_p)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_q) = v3dg(:,:,:,iv3d_q)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qc) = v3dg(:,:,:,iv3d_qc)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qr) = v3dg(:,:,:,iv3d_qr)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qi) = v3dg(:,:,:,iv3d_qi)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qs) = v3dg(:,:,:,iv3d_qs)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qg) = v3dg(:,:,:,iv3d_qg)

  ! 3D height: calculate from ph; re-constructed from p, t, q, mu if mode = 2
  !---------------------------------------------------------

  Rd_wrf = 287.0_r_size
  RvRd_wrf = 461.6_r_size / Rd_wrf

!$OMP PARALLEL DO PRIVATE(i,j,k,tmprho,tmpph) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      if (mode_ == 2) then
        do k = 1, nlev
          tmprho = v3dg(k,i,j,iv3d_p) / (v3dg(k,i,j,iv3d_t) * (1.0_r_size + v3dg(k,i,j,iv3d_q) * RvRd_wrf) * Rd_wrf)
          if (k == 1) then
            tmpph(k) = ph_base(k,i,j) + v2dg(i,j,iv2d_mu) * (eta_stag(k) - eta_stag(k+1)) / tmprho
          else
            tmpph(k) = v3dg(k-1,i,j,iv3d_ph) + v2dg(i,j,iv2d_mu) * (eta_stag(k) - eta_stag(k+1)) / tmprho
          end if
        end do
        v3dgh(1+KHALO,i+IHALO,j+JHALO,iv3dd_hgt) = 0.5_r_size * (ph_base(1,i,j) + tmpph(1)) / gg
        do k = 2, nlev
          v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_hgt) = 0.5_r_size * (tmpph(k-1) + tmpph(k)) / gg
        end do
      else
        v3dgh(1+KHALO,i+IHALO,j+JHALO,iv3dd_hgt) = 0.5_r_size * (ph_base(1,i,j) + v3dg(1,i,j,iv3d_ph)) / gg
        do k = 2, nlev
          v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_hgt) = 0.5_r_size * (v3dg(k-1,i,j,iv3d_ph) + v3dg(k,i,j,iv3d_ph)) / gg
        end do
      end if
    end do
  end do
!$OMP END PARALLEL DO

  ! RH
  !---------------------------------------------------------

!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_rh) = [[RH calculator]]

  ! Topographic height: calculate from ph_base
  !---------------------------------------------------------

  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_topo) = ph_base(1,:,:) / gg

  ! Surface pressure: calculate from mu, q*, eta_stag, ptop
  !---------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(i,j,k,moist_to_dry) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      moist_to_dry = 1.0_r_size
      do k = 1, nlev
        moist_to_dry = moist_to_dry + (v3dg(k,i,j,iv3d_q)  + v3dg(k,i,j,iv3d_qc) + v3dg(k,i,j,iv3d_qr) + &
                                       v3dg(k,i,j,iv3d_qi) + v3dg(k,i,j,iv3d_qs) + v3dg(k,i,j,iv3d_qg)) * (eta_stag(k) - eta_stag(k+1))
      end do
      v2dgh(i+IHALO,j+IHALO,iv2dd_ps) = v2dg(i,j,iv2d_mu) * moist_to_dry + ptop
    end do
  end do
!$OMP END PARALLEL DO

  ! Surface variables: use the 1st level as the surface (although it is not)
  !---------------------------------------------------------

  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_u10m) = v3dg(1,:,:,iv3d_u)
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_v10m) = v3dg(1,:,:,iv3d_v)
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_t2m) = v3dg(1,:,:,iv3d_t)
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_q2m) = v3dg(1,:,:,iv3d_q)

!#ifdef H08
!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_skint) = v3dg(1,:,:,iv3d_t)

!  ! Assume the point where terrain height is less than 10 m is the ocean. T.Honda (02/09/2016)
!  !---------------------------------------------------------

!!$OMP PARALLEL DO PRIVATE(j,i)
!  do j = 1, nlat
!    do i = 1, nlon
!      v2dgh(i+IHALO,j+JHALO,iv2dd_lsmask) = min(max(topo(i,j) - 10.0d0, 0.0d0), 1.0d0)
!    enddo
!  enddo
!!$OMP END PARALLEL DO
!#endif

  ! Pad the upper and lower halo areas
  !---------------------------------------------------------

  do iv3d = 1, nv3dd
    do j  = JS, JE
      do i  = IS, IE
        v3dgh(   1:KS-1,i,j,iv3d) = v3dgh(KS,i,j,iv3d)
        v3dgh(KE+1:KA,  i,j,iv3d) = v3dgh(KE,i,j,iv3d)
      end do
    end do
  end do

  ! Communicate the lateral halo areas
  !---------------------------------------------------------

  do iv3d = 1, nv3dd
    call COMM_vars8( v3dgh(:,:,:,iv3d), iv3d )
  end do
  do iv3d = 1, nv3dd
    call COMM_wait ( v3dgh(:,:,:,iv3d), iv3d )
  end do

  do iv2d = 1, nv2dd
    call COMM_vars8( v2dgh(:,:,iv2d), iv2d )
  end do
  do iv2d = 1, nv2dd
    call COMM_wait ( v2dgh(:,:,iv2d), iv2d )
  end do

!write (6, *) '###### U',    v3dgh(:,10,10,iv3dd_u)
!write (6, *) '###### V',    v3dgh(:,10,10,iv3dd_v)
!write (6, *) '###### W',    v3dgh(:,10,10,iv3dd_w)
!write (6, *) '###### T',    v3dgh(:,10,10,iv3dd_t)
!write (6, *) '###### P',    v3dgh(:,10,10,iv3dd_p)
!write (6, *) '###### HGT',  v3dgh(:,10,10,iv3dd_hgt)
!write (6, *) '###### Q',    v3dgh(:,10,10,iv3dd_q)
!write (6, *) '###### QC',   v3dgh(:,10,10,iv3dd_qc)
!write (6, *) '###### QR',   v3dgh(:,10,10,iv3dd_qr)
!write (6, *) '###### QI',   v3dgh(:,10,10,iv3dd_qi)
!write (6, *) '###### QS',   v3dgh(:,10,10,iv3dd_qs)
!write (6, *) '###### QG',   v3dgh(:,10,10,iv3dd_qg)
!write (6, *) '###### TOPO', v2dgh(10,10,iv2dd_topo)
!write (6, *) '###### PS',   v2dgh(10,10,iv2dd_ps)
!write (6, *) '###### U10M', v2dgh(10,10,iv2dd_u10m)
!write (6, *) '###### V10M', v2dgh(10,10,iv2dd_v10m)
!write (6, *) '###### T2M',  v2dgh(10,10,iv2dd_t2m)
!write (6, *) '###### Q2M',  v2dgh(10,10,iv2dd_q2m)

  return
end subroutine state_to_history

!===============================================================================
END MODULE common_wrf
