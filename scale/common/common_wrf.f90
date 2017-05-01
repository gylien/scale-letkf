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
  REAL(r_size),SAVE :: t00
  REAL(r_size),SAVE :: p00
  REAL(r_size),SAVE :: ptop
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

!!-------------------------------------------------------------------------------
!! Initialize standard I/O and read common namelist of SCALE-LETKF
!!-------------------------------------------------------------------------------
!subroutine set_common_conf(nprocs)
!  use scale_stdio

!  implicit none
!  integer, intent(in) :: nprocs
!  character(len=H_MID), parameter :: MODELNAME = "WRF-LETKF"

!  ! setup standard I/O
!  call IO_setup( MODELNAME, .false.)

!  call read_nml_ensemble
!  call read_nml_letkf_prc

!  if (nprocs /= NNODES * PPN) then
!    write(6,'(A,I10)') 'Number of MPI processes = ', nprocs
!    write(6,'(A,I10)') 'NNODES = ', NNODES
!    write(6,'(A,I10)') 'PPN    = ', PPN
!    write(6,'(A)') 'Number of MPI processes should be equal to NNODES * PPN.'
!    stop
!  end if

!  return
!end subroutine set_common_conf

!!-------------------------------------------------------------------------------
!! Set the parameters related to the SCALE model
!!-------------------------------------------------------------------------------
!SUBROUTINE set_common_scale
!  use scale_rm_process, only: &
!    PRC_NUM_X, &
!    PRC_NUM_Y
!  use scale_grid_index, only: &
!    IMAX, &
!    JMAX, &
!    KMAX, &
!    IHALO, &
!    JHALO, &
!    KHALO

!  IMPLICIT NONE
!!  REAL(r_sngl) :: slat(nlat), wlat(nlat)
!!  REAL(r_size) :: totalwg, wgtmp, latm1, latm2
!!  INTEGER :: i,j

!  WRITE(6,'(A)') 'Hello from set_common_scale'

!  !
!  ! Set up node and process distribution
!  !
!  if (MEM_NP /= PRC_NUM_X * PRC_NUM_Y) then
!    write(6,'(A,I10)') 'MEM_NP    = ', MEM_NP
!    write(6,'(A,I10)') 'PRC_NUM_X = ', PRC_NUM_X
!    write(6,'(A,I10)') 'PRC_NUM_Y = ', PRC_NUM_Y
!    write(6,'(A)') 'MEM_NP should be equal to PRC_NUM_X * PRC_NUM_Y.'
!    stop
!  end if

!  nlon = IMAX
!  nlat = JMAX
!  nlev = KMAX
!  nlong = nlon * PRC_NUM_X
!  nlatg = nlat * PRC_NUM_Y
!  nlonh = nlon + IHALO * 2
!  nlath = nlat + JHALO * 2
!  nlevh = nlev + KHALO * 2

!  nij0 = nlon * nlat
!  nlevall  = nlev * nv3d  + nv2d
!  nlevalld = nlev * nv3dd + nv2dd
!  ngpv  = nij0 * nlevall
!  ngpvd = nij0 * nlevalld

!  t00 = undef
!  p00 = undef
!  ptop = undef

!!  !
!!  ! Lon, Lat
!!  !
!!!$OMP PARALLEL DO PRIVATE(i)
!!  DO i=1,nlon
!!    lon(i) = 360.d0/nlon*(i-1)
!!  END DO
!!!$OMP END PARALLEL DO
!!  CALL SPLAT(idrt,nlat,slat,wlat)
!!  do j=1,nlat
!!    lat(j) = 180.d0/pi*asin(slat(nlat-j+1))
!!  end do
!!  !
!!  ! dx and dy
!!  !
!!!$OMP PARALLEL
!!!$OMP WORKSHARE
!!  dx(:) = 2.0d0 * pi * re * cos(lat(:) * pi / 180.0d0) / REAL(nlon,r_size)
!!!$OMP END WORKSHARE

!!!$OMP DO
!!  DO i=1,nlat-1
!!    dy(i) = 2.0d0 * pi * re * (lat(i+1) - lat(i)) / 360.0d0
!!  END DO
!!!$OMP END DO
!!!$OMP END PARALLEL
!!  dy(nlat) = 2.0d0 * pi * re * (90.0d0 - lat(nlat)) / 180.0d0

!!!$OMP PARALLEL DO
!!  DO i=2,nlat
!!    dy2(i) = (dy(i-1) + dy(i)) * 0.5d0
!!  END DO
!!!$OMP END PARALLEL DO
!!  dy2(1) = (dy(nlat) + dy(1)) * 0.5d0
!!  !
!!  ! Corioris parameter
!!  !
!!!$OMP PARALLEL WORKSHARE
!!  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)
!!!$OMP END PARALLEL WORKSHARE
!!  !
!!  ! Weight for global average
!!  !
!!  totalwg = 0.0_r_size
!!  DO j=1,nlat
!!    if (j == 1) then
!!      latm1 = -0.5d0*pi !-90 degree
!!    else
!!      latm1 = 0.5d0*(lat(j-1) + lat(j))*pi/180.0d0
!!    end if
!!    if (j == nlat) then
!!      latm2 = 0.5d0*pi !90 degree
!!    else
!!      latm2 = 0.5d0*(lat(j) + lat(j+1))*pi/180.0d0
!!    end if
!!    wgtmp = abs(sin(latm2) - sin(latm1))
!!    wg(:,j) = wgtmp
!!    totalwg = totalwg + wgtmp * nlon
!!  END DO
!!  totalwg = 1.0_r_size / totalwg
!!  wg(:,:) = sqrt(wg(:,:) * totalwg)

!  RETURN
!END SUBROUTINE set_common_scale

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
    if (.not. allocated(p_base3d)) allocate(p_base3d(nlev,nlon,nlat))
    if (.not. allocated(ph_base3d)) allocate(ph_base3d(nlev+1,nlon,nlat))
    if (.not. allocated(mu_base2d)) allocate(mu_base2d(nlon,nlat))
    if (.not. allocated(eta_stag)) allocate(eta_stag(nlev+1))
    call read_coor_par(filename, comm, t00=t00, p00=p00, p_base=p_base3d, ph_base=ph_base3d, mu_base=mu_base2d, eta_stag=eta_stag)
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
    if (.not. allocated(p_base3d)) allocate(p_base3d(nlev,nlon,nlat))
    if (.not. allocated(ph_base3d)) allocate(ph_base3d(nlev+1,nlon,nlat))
    if (.not. allocated(mu_base2d)) allocate(mu_base2d(nlon,nlat))
    if (.not. allocated(eta_stag)) allocate(eta_stag(nlev+1))
    call read_coor_par(filename, comm, t00=t00, p00=p00, p_base=p_base3d, ph_base=ph_base3d, mu_base=mu_base2d, eta_stag=eta_stag)
  end if

  if (trans_) then
    kappa = 2.0_r_size / 7.0_r_size
    Rd_wrf = 287.0_r_size
    RvRd_wrf = 461.6_r_size / Rd_wrf
    p00_kappa = p00 ** kappa

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

  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename)

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
    if (.not. allocated(p_base3d)) allocate(p_base3d(nlev,nlon,nlat))
    if (.not. allocated(ph_base3d)) allocate(ph_base3d(nlev+1,nlon,nlat))
    if (.not. allocated(mu_base2d)) allocate(mu_base2d(nlon,nlat))
    if (.not. allocated(eta_stag)) allocate(eta_stag(nlev+1))
!    call read_coor_par(filename, comm, t00=t00, p00=p00, p_base=p_base3d, ph_base=ph_base3d, mu_base=mu_base2d)
    call read_coor_par(filename, comm, t00=t00, p00=p00, p_base=p_base3d, ph_base=ph_base3d, mu_base=mu_base2d, eta_stag=eta_stag)
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
subroutine state_to_history(v3dg, v2dg, topo, v3dgh, v2dgh)
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
  real(RP), intent(in) :: topo(nlon,nlat)
  real(r_size), intent(out) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  real(r_size), intent(out) :: v2dgh(nlonh,nlath,nv2dd)

!  real(r_size) :: height(nlev,nlon,nlat)
!  integer :: i, j, k, iv3d, iv2d

!  ! Variables that can be directly copied
!  !---------------------------------------------------------

!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_u) = v3dg(:,:,:,iv3d_u)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_v) = v3dg(:,:,:,iv3d_v)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_w) = v3dg(:,:,:,iv3d_w)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_t) = v3dg(:,:,:,iv3d_t)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_p) = v3dg(:,:,:,iv3d_p)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_q) = v3dg(:,:,:,iv3d_q)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qc) = v3dg(:,:,:,iv3d_qc)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qr) = v3dg(:,:,:,iv3d_qr)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qi) = v3dg(:,:,:,iv3d_qi)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qs) = v3dg(:,:,:,iv3d_qs)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qg) = v3dg(:,:,:,iv3d_qg)

!  ! RH
!  !---------------------------------------------------------

!!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_rh) = [[RH calculator]]

!  ! Calculate height based the the topography and vertical coordinate
!  !---------------------------------------------------------

!  call scale_calc_z(topo, height)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_hgt) = height

!  ! Surface variables: use the 1st level as the surface (although it is not)
!  !---------------------------------------------------------

!  v2dgh(:,:,iv2dd_topo) = v3dgh(1+KHALO,:,:,iv3dd_hgt)                ! Use the first model level as topography (is this good?)
!!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_topo) = topo(:,:) ! Use the real topography

!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_ps) = v3dg(1,:,:,iv3d_p)
!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_u10m) = v3dg(1,:,:,iv3d_u)
!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_v10m) = v3dg(1,:,:,iv3d_v)
!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_t2m) = v3dg(1,:,:,iv3d_t)
!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_q2m) = v3dg(1,:,:,iv3d_q)

!!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_rain) = [[No way]]

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

!  ! Pad the upper and lower halo areas
!  !---------------------------------------------------------

!  do iv3d = 1, nv3dd
!    do j  = JS, JE
!      do i  = IS, IE
!        v3dgh(   1:KS-1,i,j,iv3d) = v3dgh(KS,i,j,iv3d)
!        v3dgh(KE+1:KA,  i,j,iv3d) = v3dgh(KE,i,j,iv3d)
!      end do
!    end do
!  end do

!  ! Communicate the lateral halo areas
!  !---------------------------------------------------------

!  do iv3d = 1, nv3dd
!    call COMM_vars8( v3dgh(:,:,:,iv3d), iv3d )
!  end do
!  do iv3d = 1, nv3dd
!    call COMM_wait ( v3dgh(:,:,:,iv3d), iv3d )
!  end do

!  do iv2d = 1, nv2dd
!    call COMM_vars8( v2dgh(:,:,iv2d), iv2d )
!  end do
!  do iv2d = 1, nv2dd
!    call COMM_wait ( v2dgh(:,:,iv2d), iv2d )
!  end do

  return
end subroutine state_to_history

!-------------------------------------------------------------------------------
! Monitor state variables
!-------------------------------------------------------------------------------
!SUBROUTINE monit_grd(v3d,v2d)
!  IMPLICIT NONE
!  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
!  INTEGER :: k,n

!  DO k=1,nlev
!    WRITE(6,'(I2,A)') k,'th level'
!    DO n=1,nv3d
!      WRITE(6,'(A,2ES10.2)') v3dname(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
!    END DO
!  END DO

!  DO n=1,nv2d
!    WRITE(6,'(A,2ES10.2)') v2dname(n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
!  END DO

!  RETURN
!END SUBROUTINE monit_grd

subroutine scale_calc_z(phb, ph, z)
  implicit none

  real(r_size), intent(in) :: phb(nlev+1,nlon,nlat)
  real(r_size), intent(in) :: ph(nlev+1,nlon,nlat)
  real(r_size), intent(out) :: z(nlev,nlon,nlat)
  integer :: i, j, k

!$OMP PARALLEL DO PRIVATE(j,i,k) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev
        z(k, i, j) = 0.5_r_size * (phb(k,i,j) + phb(k+1,i,j) + ph(k,i,j) + ph(k+1,i,j)) / gg
      end do
    enddo
  enddo
!$OMP END PARALLEL DO

  return
end subroutine scale_calc_z

!!-------------------------------------------------------------------------------
!! Calculate 3D height coordinate given the topography height (on scattered grids)
!!-------------------------------------------------------------------------------
!! [INPUT]
!!   nij         : scattered grid numbers
!!   topo(nij)   : topography height (on scattered grids)
!! [OUTPUT]
!!   z(nij,nlev) : 3D height coordinate (on scattered grids)
!!-------------------------------------------------------------------------------
!subroutine scale_calc_z_grd(nij, topo, z)
!  use scale_grid, only: &
!     GRID_CZ, &
!     GRID_FZ
!  use scale_grid_index, only: &
!     KHALO, KS, KE
!  implicit none

!  integer, intent(in) :: nij
!  real(r_size), intent(in) :: topo(nij)
!  real(r_size), intent(out) :: z(nij,nlev)
!  real(r_size) :: ztop
!  integer :: k, i

!  ztop = GRID_FZ(KE) - GRID_FZ(KS-1)
!!$OMP PARALLEL DO PRIVATE(i,k)
!  do k = 1, nlev
!    do i = 1, nij
!      z(i,k) = (ztop - topo(i)) / ztop * GRID_CZ(k+KHALO) + topo(i)
!    end do
!  end do
!!$OMP END PARALLEL DO

!  return
!end subroutine scale_calc_z_grd

!!-------------------------------------------------------------------------------
!! Calculate ensemble mean (on scattered grids)
!!-------------------------------------------------------------------------------
!! [INPUT]
!!   mem                     : ensemble size
!!   nens                    : ensemble demension of state variables
!!   nij                     : scattered grid numbers
!!   v3d(nij,nlev,nens,nv3d) : 3D ensemble state variables (on scattered grids)
!!                             inputted by (:,:,1..mem,:)
!!   v2d(nij,     nens,nv3d) : 2D ensemble state variables (on scattered grids)
!!                             inputted by (:,  1..mem,:)
!! [OUTPUT]
!!   v3d(nij,nlev,nens,nv3d) : ensemble mean of 3D state variables (on scattered grids)
!!                             outputted by (:,:,mem+1,:)
!!   v2d(nij,     nens,nv3d) : ensemble mean of 2D state variables (on scattered grids)
!!                             outputted by (:  ,mem+1,:)
!!-------------------------------------------------------------------------------
!subroutine ensmean_grd(mem, nens, nij, v3d, v2d)
!  implicit none
!  integer, intent(in) :: mem
!  integer, intent(in) :: nens
!  integer, intent(in) :: nij
!  real(r_size), intent(inout) :: v3d(nij,nlev,nens,nv3d)
!  real(r_size), intent(inout) :: v2d(nij,nens,nv2d)
!  integer :: i, k, m, n, mmean

!  mmean = mem + 1

!!$OMP PARALLEL DO PRIVATE(i,k,m,n) COLLAPSE(3)
!  do n = 1, nv3d
!    do k = 1, nlev
!      do i = 1, nij
!        v3d(i,k,mmean,n) = v3d(i,k,1,n)
!        do m = 2, mem
!          v3d(i,k,mmean,n) = v3d(i,k,mmean,n) + v3d(i,k,m,n)
!        end do
!        v3d(i,k,mmean,n) = v3d(i,k,mmean,n) / real(mem, r_size)
!      end do
!    end do
!  end do
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE(i,m,n) COLLAPSE(2)
!  do n = 1, nv2d
!    do i = 1, nij
!      v2d(i,mmean,n) = v2d(i,1,n)
!      do m = 2, mem
!        v2d(i,mmean,n) = v2d(i,mmean,n) + v2d(i,m,n)
!      end do
!      v2d(i,mmean,n) = v2d(i,mmean,n) / real(mem, r_size)
!    end do
!  end do
!!$OMP END PARALLEL DO

!  return
!end subroutine ensmean_grd

!!-------------------------------------------------------------------------------
!! Calculate ensemble spread (on scattered grids)
!! * the ensemble mean has to be already calculated
!!-------------------------------------------------------------------------------
!! [INPUT]
!!   mem                     : ensemble size
!!   nens                    : ensemble demension of state variables
!!   nij                     : scattered grid numbers
!!   v3d(nij,nlev,nens,nv3d) : 3D ensemble state variables with mean (on scattered grids)
!!                             inputted by (:,:,1..mem+1,:)
!!   v2d(nij,     nens,nv3d) : 2D ensemble state variables with mean (on scattered grids)
!!                             inputted by (:,  1..mem+1,:)
!! [OUTPUT]
!!   v3ds(nij,nlev,nv3d)     : ensemble spread of 3D state variables (on scattered grids)
!!   v2ds(nij,     nv3d)     : ensemble spread of 2D state variables (on scattered grids)
!!-------------------------------------------------------------------------------
!subroutine enssprd_grd(mem, nens, nij, v3d, v2d, v3ds, v2ds)
!  implicit none
!  integer, intent(in) :: mem
!  integer, intent(in) :: nens
!  integer, intent(in) :: nij
!  real(r_size), intent(in) :: v3d(nij,nlev,nens,nv3d)
!  real(r_size), intent(in) :: v2d(nij,nens,nv2d)
!  real(r_size), intent(out) :: v3ds(nij,nlev,nv3d)
!  real(r_size), intent(out) :: v2ds(nij,nv2d)
!  integer :: i, k, m, n, mmean

!  mmean = mem + 1

!!$OMP PARALLEL DO PRIVATE(i,k,m,n) COLLAPSE(3)
!  do n = 1, nv3d
!    do k = 1, nlev
!      do i = 1, nij
!        v3ds(i,k,n) = (v3d(i,k,1,n) - v3d(i,k,mmean,n)) ** 2
!        do m = 2, mem
!          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n) - v3d(i,k,mmean,n)) ** 2
!        end do
!        v3ds(i,k,n) = sqrt(v3ds(i,k,n) / real(mem-1, r_size))
!      end do
!    end do
!  end do
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE(i,m,n) COLLAPSE(2)
!  do n = 1, nv2d
!    do i = 1, nij
!      v2ds(i,n) = (v2d(i,1,n) - v2d(i,mmean,n)) ** 2
!      do m = 2, mem
!        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n) - v2d(i,mmean,n)) ** 2
!      end do
!      v2ds(i,n) = sqrt(v2ds(i,n) / real(mem-1, r_size))
!    end do
!  end do
!!$OMP END PARALLEL DO

!  return
!end subroutine enssprd_grd

!!-------------------------------------------------------------------------------
!! Convert 1D rank of process to 2D rank
!!-------------------------------------------------------------------------------
!subroutine rank_1d_2d(proc, iproc, jproc)
!  use scale_rm_process, only: PRC_2Drank
!  implicit none
!  integer, intent(in) :: proc
!  integer, intent(out) :: iproc, jproc

!  iproc = PRC_2Drank(proc,1)
!  jproc = PRC_2Drank(proc,2)

!  return  
!end subroutine rank_1d_2d

!!-------------------------------------------------------------------------------
!! Convert 2D rank of process to 1D rank
!!-------------------------------------------------------------------------------
!subroutine rank_2d_1d(iproc, jproc, proc)
!  use scale_rm_process, only: PRC_NUM_X
!  implicit none
!  integer, intent(in) :: iproc, jproc
!  integer, intent(out) :: proc

!  proc = jproc * PRC_NUM_X + iproc

!  return  
!end subroutine rank_2d_1d

!!-------------------------------------------------------------------------------
!! Convert <integer> global grid coordinates (i,j) to local given the 1D rank of process
!!-------------------------------------------------------------------------------
!subroutine ij_g2l(proc, ig, jg, il, jl)
!  implicit none
!  integer, intent(in) :: proc
!  integer, intent(in) :: ig
!  integer, intent(in) :: jg
!  integer, intent(out) :: il
!  integer, intent(out) :: jl
!  integer :: iproc, jproc

!  call rank_1d_2d(proc, iproc, jproc)
!  il = ig - iproc * nlon
!  jl = jg - jproc * nlat

!  return  
!end subroutine ij_g2l

!!-------------------------------------------------------------------------------
!! Convert <integer> local grid coordinates (i,j) to global given the 1D rank of process
!!-------------------------------------------------------------------------------
!subroutine ij_l2g(proc, il, jl, ig, jg)
!  implicit none
!  integer, intent(in) :: proc
!  integer, intent(in) :: il
!  integer, intent(in) :: jl
!  integer, intent(out) :: ig
!  integer, intent(out) :: jg
!  integer :: iproc, jproc

!  call rank_1d_2d(proc, iproc, jproc)
!  ig = il + iproc * nlon
!  jg = jl + jproc * nlat

!  return  
!end subroutine ij_l2g

!!-------------------------------------------------------------------------------
!! Convert <real> global grid coordinates (i,j) to local given the 1D rank of process
!!-------------------------------------------------------------------------------
!subroutine rij_g2l(proc, ig, jg, il, jl)
!  implicit none
!  integer, intent(in) :: proc
!  real(r_size), intent(in) :: ig
!  real(r_size), intent(in) :: jg
!  real(r_size), intent(out) :: il
!  real(r_size), intent(out) :: jl
!  integer :: iproc, jproc

!  call rank_1d_2d(proc, iproc, jproc)
!  il = ig - real(iproc * nlon,r_size)
!  jl = jg - real(jproc * nlat,r_size)

!  return  
!end subroutine rij_g2l

!!-------------------------------------------------------------------------------
!! Convert <real> local grid coordinates (i,j) to global given the 1D rank of process
!!-------------------------------------------------------------------------------
!subroutine rij_l2g(proc, il, jl, ig, jg)
!  implicit none
!  integer, intent(in) :: proc
!  real(r_size), intent(in) :: il
!  real(r_size), intent(in) :: jl
!  real(r_size), intent(out) :: ig
!  real(r_size), intent(out) :: jg
!  integer :: iproc, jproc

!  call rank_1d_2d(proc, iproc, jproc)
!  ig = il + real(iproc * nlon,r_size)
!  jg = jl + real(jproc * nlat,r_size)

!  return  
!end subroutine rij_l2g

!!-------------------------------------------------------------------------------
!! Convert <real> global grid coordinates (i,j) to local where the grid resides
!! * HALO grids are used
!!-------------------------------------------------------------------------------
!! [INPUT]
!!   ig, jg : global grid coordinates
!! [OUTPUT]
!!   proc   : the 1D rank of process where the grid resides;
!!            * return -1 if the grid is outside of the global domain
!!   il, jl : local grid coordinates
!!-------------------------------------------------------------------------------
!subroutine rij_g2l_auto(proc,ig,jg,il,jl)
!  use scale_rm_process, only: &
!      PRC_NUM_X, PRC_NUM_Y
!#ifdef DEBUG
!  use scale_grid_index, only: &
!      IHALO, JHALO, &
!      IA, JA
!  use scale_process, only: &
!      PRC_myrank
!  use scale_grid, only: &
!      GRID_CX, &
!      GRID_CY, &
!      GRID_CXG, &
!      GRID_CYG, &
!      DX, &
!      DY
!#else
!  use scale_grid_index, only: &
!      IHALO, JHALO
!#endif
!  implicit none
!  integer, intent(out) :: proc
!  real(r_size), intent(in) :: ig
!  real(r_size), intent(in) :: jg
!  real(r_size), intent(out) :: il
!  real(r_size), intent(out) :: jl
!  integer :: iproc, jproc

!  if (ig < real(1+IHALO,r_size) .or. ig > real(nlon*PRC_NUM_X+IHALO,r_size) .or. &
!      jg < real(1+JHALO,r_size) .or. jg > real(nlat*PRC_NUM_Y+JHALO,r_size)) then
!    il = -1.0d0
!    jl = -1.0d0
!    proc = -1
!    return
!  end if

!  iproc = ceiling((ig-real(IHALO,r_size)-0.5d0) / real(nlon,r_size)) - 1
!  jproc = ceiling((jg-real(JHALO,r_size)-0.5d0) / real(nlat,r_size)) - 1
!  il = ig - real(iproc * nlon,r_size)
!  jl = jg - real(jproc * nlat,r_size)
!  call rank_2d_1d(iproc,jproc,proc)

!#ifdef DEBUG
!  if (PRC_myrank == proc) then
!    if (ig < (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0 .or. &
!        ig > (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0 .or. &
!        jg < (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0 .or. &
!        jg > (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0) then
!      write (6,'(A)') '[Error] Process assignment fails!'
!      write (6,'(3F10.2)') ig, (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0, (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0
!      write (6,'(3F10.2)') jg, (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0, (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0
!      stop
!    end if
!  end if
!#endif

!  return
!end subroutine rij_g2l_auto

!===============================================================================
END MODULE common_wrf
