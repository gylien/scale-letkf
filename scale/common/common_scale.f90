MODULE common_scale
!===============================================================================
!
! [PURPOSE:] Common Information for SCALE
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   10/03/2012 Guo-Yuan Lien     modified for GFS model
!   07/24/2014 Guo-Yuan Lien     modified for SCALE model
!   .......... See git history for the following revisions
!
!===============================================================================
!$USE OMP_LIB
  USE common
  use common_nml

  use scale_precision, only: RP, SP, DP
  use scale_io, only: H_MID, H_LONG
  use scale_prof

  IMPLICIT NONE
  PUBLIC
!-------------------------------------------------------------------------------
! General parameters
!-------------------------------------------------------------------------------

  character(len=H_MID), parameter :: modelname = "SCALE-LETKF"
  INTEGER,PARAMETER :: vname_max = 10

  character(len=H_LONG), save :: confname

  ! 
  !--- 3D, 2D state variables (in SCALE restart files)
  ! 
  ! Parameter 'nv3d' is set in common_nml.f90
  ! Parameter 'nv2d' is set in common_nml.f90
  INTEGER,PARAMETER :: iv3d_rho=1  !-- State in restart files
  INTEGER,PARAMETER :: iv3d_rhou=2 !
  INTEGER,PARAMETER :: iv3d_rhov=3 !
  INTEGER,PARAMETER :: iv3d_rhow=4 !
  INTEGER,PARAMETER :: iv3d_rhot=5 !
  INTEGER,PARAMETER :: iv3d_u=1    !-- State for LETKF
  INTEGER,PARAMETER :: iv3d_v=2    !
  INTEGER,PARAMETER :: iv3d_w=3    !
  INTEGER,PARAMETER :: iv3d_t=4    !
  INTEGER,PARAMETER :: iv3d_p=5    !
  INTEGER,PARAMETER :: iv3d_q=6    !
  INTEGER,PARAMETER :: iv3d_qc=7   !
  INTEGER,PARAMETER :: iv3d_qr=8   !
  INTEGER,PARAMETER :: iv3d_qi=9   !
  INTEGER,PARAMETER :: iv3d_qs=10  !
  INTEGER,PARAMETER :: iv3d_qg=11  !
  CHARACTER(vname_max),PARAMETER :: v3d_name(nv3d) = &
     (/'DENS', 'MOMX', 'MOMY', 'MOMZ', 'RHOT', &
       'QV', 'QC', 'QR', 'QI', 'QS', 'QG'/)
  CHARACTER(vname_max) :: v2d_name(nv2d)

  ! 
  !--- 3D, 2D diagnostic variables (in SCALE history files)
  ! 
  INTEGER,PARAMETER :: nv3dd=13
  INTEGER,PARAMETER :: nv2dd=7
  INTEGER,PARAMETER :: iv3dd_u=1
  INTEGER,PARAMETER :: iv3dd_v=2
  INTEGER,PARAMETER :: iv3dd_w=3
  INTEGER,PARAMETER :: iv3dd_t=4
  INTEGER,PARAMETER :: iv3dd_p=5
  INTEGER,PARAMETER :: iv3dd_q=6
  INTEGER,PARAMETER :: iv3dd_qc=7
  INTEGER,PARAMETER :: iv3dd_qr=8
  INTEGER,PARAMETER :: iv3dd_qi=9
  INTEGER,PARAMETER :: iv3dd_qs=10
  INTEGER,PARAMETER :: iv3dd_qg=11
  INTEGER,PARAMETER :: iv3dd_rh=12
  INTEGER,PARAMETER :: iv3dd_hgt=13
  INTEGER,PARAMETER :: iv2dd_topo=1
  INTEGER,PARAMETER :: iv2dd_ps=2
  INTEGER,PARAMETER :: iv2dd_rain=3
  INTEGER,PARAMETER :: iv2dd_u10m=4
  INTEGER,PARAMETER :: iv2dd_v10m=5
  INTEGER,PARAMETER :: iv2dd_t2m=6
  INTEGER,PARAMETER :: iv2dd_q2m=7
  CHARACTER(vname_max),PARAMETER :: v3dd_name(nv3dd) = &
     (/'U', 'V', 'W', 'T', 'PRES', &
       'QV', 'QC', 'QR', 'QI', 'QS', 'QG', 'RH', 'height'/)
  LOGICAL,PARAMETER :: v3dd_hastime(nv3dd) = &
     (/.true., .true., .true., .true., .true., &
       .true., .true., .true., .true., .true., .true., .true., .false./)
  CHARACTER(vname_max),PARAMETER :: v2dd_name(nv2dd) = &
     (/'topo', 'SFC_PRES', 'PREC', 'U10', 'V10', 'T2', 'Q2'/)
  LOGICAL,PARAMETER :: v2dd_hastime(nv2dd) = &
     (/.false., .true., .true., .true., .true., .true., .true./)

  ! 
  !--- Variables for model coordinates
  ! 
  REAL(r_size),ALLOCATABLE,SAVE :: height3d(:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lon2d(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lat2d(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: topo2d(:,:)
  reaL(r_sngl), allocatable, save :: lsmask2dgs(:,:)
  CHARACTER(vname_max),PARAMETER :: height3d_name = 'height' ! (in SCALE restart files)
  CHARACTER(vname_max),PARAMETER :: lon2d_name = 'lon'       ! (in SCALE restart files)
  CHARACTER(vname_max),PARAMETER :: lat2d_name = 'lat'       ! (in SCALE restart files)
  CHARACTER(vname_max),PARAMETER :: topo2d_name = 'TOPO'     ! (in SCALE topo files)

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

!  INTEGER,PARAMETER :: nlonsub=200
!  INTEGER,PARAMETER :: nlatsub=200
!  INTEGER,PARAMETER :: nlonns=6
!  INTEGER,PARAMETER :: nlatns=6
!  INTEGER,PARAMETER :: nlon=nlonsub*nlonns
!  INTEGER,PARAMETER :: nlat=nlatsub*nlatns
!  INTEGER,PARAMETER :: nlev=60
!  integer,parameter :: nlonhalo=nlonsub+4
!  integer,parameter :: nlathalo=nlatsub+4
!  integer,parameter :: nlevhalo=nlev+4

!  REAL(r_size),SAVE :: lon(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: lat(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: lonu(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: latu(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: lonv(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: latv(nlonsub,nlatsub)
!!  REAL(r_size),SAVE :: dx(nlatsub)
!!  REAL(r_size),SAVE :: dy(nlatsub)
!!  REAL(r_size),SAVE :: dy2(nlatsub)
!  REAL(r_size),SAVE :: fcori(nlatsub)
!!  REAL(r_size),SAVE :: wg(nlonsub,nlatsub)

  ! 
  !--- time labels for files
  ! 
  character(len=24), save :: timelabel_obs  = ''
  character(len=20), save :: timelabel_hist = ''
  character(len=20), save :: timelabel_anal = ''

CONTAINS

!-------------------------------------------------------------------------------
! Initialize standard I/O and read common namelist of SCALE-LETKF
!-------------------------------------------------------------------------------
subroutine set_common_conf(nprocs)
  use scale_io, only: &
    IO_setup, &
    IO_ARG_getfname

  implicit none
  integer, intent(in) :: nprocs

  ! setup standard I/O
  call IO_setup( modelname )
  confname = IO_ARG_getfname( is_master=.true. )

  call read_nml_log
  call read_nml_model
  call read_nml_ensemble
  call read_nml_process

  return
end subroutine set_common_conf

!-------------------------------------------------------------------------------
! Set the parameters related to the SCALE model
!-------------------------------------------------------------------------------
SUBROUTINE set_common_scale
  use scale_prc_cartesC, only: &
    PRC_NUM_X, &
    PRC_NUM_Y
  use scale_atmos_grid_cartesC_index, only: &
    IMAX, &
    JMAX, &
    KMAX, &
    IHALO, &
    JHALO, &
    KHALO

  IMPLICIT NONE
!  REAL(r_sngl) :: slat(nlat), wlat(nlat)
!  REAL(r_size) :: totalwg, wgtmp, latm1, latm2
!  INTEGER :: i,j

!  WRITE(6,'(A)') 'Hello from set_common_scale'

  !
  ! Set up node and process distribution
  !
  nlon = IMAX
  nlat = JMAX
  nlev = KMAX
  nlong = nlon * PRC_NUM_X
  nlatg = nlat * PRC_NUM_Y
  nlonh = nlon + IHALO * 2
  nlath = nlat + JHALO * 2
  nlevh = nlev + KHALO * 2

  nij0 = nlon * nlat
  nlevall  = nlev * nv3d  + nv2d
  nlevalld = nlev * nv3dd + nv2dd
  ngpv  = nij0 * nlevall
  ngpvd = nij0 * nlevalld

  call set_lonlat2d

!  !
!  ! Lon, Lat
!  !
!!$OMP PARALLEL DO PRIVATE(i)
!  DO i=1,nlon
!    lon(i) = 360.d0/nlon*(i-1)
!  END DO
!!$OMP END PARALLEL DO
!  CALL SPLAT(idrt,nlat,slat,wlat)
!  do j=1,nlat
!    lat(j) = 180.d0/pi*asin(slat(nlat-j+1))
!  end do
!  !
!  ! dx and dy
!  !
!!$OMP PARALLEL
!!$OMP WORKSHARE
!  dx(:) = 2.0d0 * pi * re * cos(lat(:) * pi / 180.0d0) / REAL(nlon,r_size)
!!$OMP END WORKSHARE

!!$OMP DO
!  DO i=1,nlat-1
!    dy(i) = 2.0d0 * pi * re * (lat(i+1) - lat(i)) / 360.0d0
!  END DO
!!$OMP END DO
!!$OMP END PARALLEL
!  dy(nlat) = 2.0d0 * pi * re * (90.0d0 - lat(nlat)) / 180.0d0

!!$OMP PARALLEL DO
!  DO i=2,nlat
!    dy2(i) = (dy(i-1) + dy(i)) * 0.5d0
!  END DO
!!$OMP END PARALLEL DO
!  dy2(1) = (dy(nlat) + dy(1)) * 0.5d0
!  !
!  ! Corioris parameter
!  !
!!$OMP PARALLEL WORKSHARE
!  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)
!!$OMP END PARALLEL WORKSHARE
!  !
!  ! Weight for global average
!  !
!  totalwg = 0.0_r_size
!  DO j=1,nlat
!    if (j == 1) then
!      latm1 = -0.5d0*pi !-90 degree
!    else
!      latm1 = 0.5d0*(lat(j-1) + lat(j))*pi/180.0d0
!    end if
!    if (j == nlat) then
!      latm2 = 0.5d0*pi !90 degree
!    else
!      latm2 = 0.5d0*(lat(j) + lat(j+1))*pi/180.0d0
!    end if
!    wgtmp = abs(sin(latm2) - sin(latm1))
!    wg(:,j) = wgtmp
!    totalwg = totalwg + wgtmp * nlon
!  END DO
!  totalwg = 1.0_r_size / totalwg
!  wg(:,:) = sqrt(wg(:,:) * totalwg)

  RETURN
END SUBROUTINE set_common_scale

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE restart files
!-------------------------------------------------------------------------------
!SUBROUTINE read_restart(filename,v3dg,v2dg)
!  use gtool_file, only: FileRead, FileCloseAll
!  use scale_prc, only: PRC_myrank
!  use common_mpi, only: myrank
!  IMPLICIT NONE

!  CHARACTER(*),INTENT(IN) :: filename
!  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
!  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
!  INTEGER :: iv3d,iv2d

!  WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',PRC_myrank,'.nc'

!  do iv3d = 1, nv3d
!    if (LOG_LEVEL >= 1) then
!      write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
!    end if
!    call FileRead(v3dg(:,:,:,iv3d), filename, trim(v3d_name(iv3d)), 1, PRC_myrank)
!  end do

!  do iv2d = 1, nv2d
!    if (LOG_LEVEL >= 1) then
!      write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
!    end if
!    call FileRead(v2dg(:,:,iv2d), filename, trim(v2d_name(iv2d)), 1, PRC_myrank)
!  end do

!  call FileCloseAll

!  RETURN
!END SUBROUTINE read_restart

SUBROUTINE read_restart(filename,v3dg,v2dg)
  use netcdf
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_HAS_W,  &
    PRC_HAS_S
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use common_mpi, only: myrank
  use common_ncio
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: iv3d, iv2d, ncid, varid
  integer :: is, js, ks

  is = 1
  js = 1
  if (.not. PRC_HAS_W) then
    is = is + IHALO
  end if
  if (.not. PRC_HAS_S) then
    js = js + JHALO
  end if

  if (LOG_LEVEL >= 3) then
    write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',PRC_myrank,'.nc'
  endif

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
!  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)

  do iv3d = 1, nv3d
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
    end if

    if (iv3d == iv3d_rhow) then
      ks = 2 ! ignore zh=1 (the surface where MOMZ = 0.0)
    else
      ks = 1
    endif

    call ncio_check(nf90_inq_varid(ncid, trim(v3d_name(iv3d)), varid))
    call ncio_check(nf90_get_var(ncid, varid, v3dg(:,:,:,iv3d), &
                                 start = (/ ks, is, js, 1 /),    &
                                 count = (/ KMAX, IMAX, JMAX, 1 /)))
  end do

  do iv2d = 1, nv2d
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
    end if
    call ncio_check(nf90_inq_varid(ncid, trim(v2d_name(iv2d)), varid))
    call ncio_check(nf90_get_var(ncid, varid, v2dg(:,:,iv2d), &
                                 start = (/ is, js, 1 /),     &
                                 count = (/ IMAX, JMAX, 1 /)))
  end do

  call ncio_close(ncid)

  RETURN
END SUBROUTINE read_restart

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE restart files <PnetCDF>
!-------------------------------------------------------------------------------
#ifdef PNETCDF
SUBROUTINE read_restart_par(filename,v3dg,v2dg,comm)
!  use common_mpi_scale, only: &
!    MPI_COMM_a
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_2Drank
!    PRC_PERIODIC_X, PRC_PERIODIC_Y, &
!    PRC_HAS_W,  &
!    PRC_HAS_E,  &
!    PRC_HAS_S,  &
!    PRC_HAS_N
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL
  use common_mpi, only: myrank
  use pnetcdf
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  integer,intent(in) :: comm
  integer :: iv3d,iv2d,ncid

  integer :: err, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(3), count(3)

  ! calculate subarray's start() and count() to the global variables
  start(1) = 1
  start(2) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(3) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  count(1) = KMAX
  count(2) = IMAX
  count(3) = JMAX
  start(2) = start(2) + IHALO
  start(3) = start(3) + JHALO
!  if (.NOT. PRC_PERIODIC_X) start(2) = start(2) + IHALO
!  if (.NOT. PRC_PERIODIC_Y) start(3) = start(3) + JHALO

!  write (6,'(A,I6.6,3A,6I6)') 'MYRANK ',myrank,' is reading a file ',trim(filename)//'.nc', ' >> PnetCDF start(3), count(3) =', start, count

  err = nfmpi_open(comm, trim(filename)//".nc", NF_NOWRITE, MPI_INFO_NULL, ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//'.nc '//nfmpi_strerror(err)

  do iv3d = 1, nv3d
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
    end if
    err = nfmpi_inq_varid(ncid, trim(v3d_name(iv3d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
#ifdef SINGLE
    err = nfmpi_iget_vara_real(ncid, varid, start, count, v3dg(:,:,:,iv3d), req)
#else
    err = nfmpi_iget_vara_double(ncid, varid, start, count, v3dg(:,:,:,iv3d), req)
#endif
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_get_vara_double_all '//' '//nfmpi_strerror(err)
  end do

  do iv2d = 1, nv2d
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
    end if
    err = nfmpi_inq_varid(ncid, trim(v2d_name(iv2d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
#ifdef SINGLE
    err = nfmpi_iget_vara_real(ncid, varid, start(2:3), count(2:3), v2dg(:,:,iv2d), req)
#else
    err = nfmpi_iget_vara_double(ncid, varid, start(2:3), count(2:3), v2dg(:,:,iv2d), req)
#endif
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_get_vara_double_all '//' '//nfmpi_strerror(err)
  end do

  err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  RETURN
END SUBROUTINE read_restart_par
#endif

!-------------------------------------------------------------------------------
! [Direct transfer] Read SCALE restart files
!-------------------------------------------------------------------------------
subroutine read_restart_direct(v3dg,v2dg)
  use mod_atmos_vars, only: &
    DENS, &
    MOMX, &
    MOMY, &
    MOMZ, &
    RHOT, &
    QV, QC, QR, &
    QI, QS, QG
  use scale_atmos_grid_cartesC_index, only: &
    IS, IE, JS, JE, KS, KE
  implicit none

  real(RP), intent(out) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP), intent(out) :: v2dg(nlon,nlat,nv2d)
  integer :: iv3d, iv2d

  do iv3d = 1, nv3d
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 3D var [direct transfer]: ', trim(v3d_name(iv3d))
    end if
    select case (iv3d)
    case (iv3d_rho)
      v3dg(:,:,:,iv3d) = DENS(KS:KE,IS:IE,JS:JE)
    case (iv3d_rhou)
      v3dg(:,:,:,iv3d) = MOMX(KS:KE,IS:IE,JS:JE)
    case (iv3d_rhov)
      v3dg(:,:,:,iv3d) = MOMY(KS:KE,IS:IE,JS:JE)
    case (iv3d_rhow)
      v3dg(:,:,:,iv3d) = MOMZ(KS:KE,IS:IE,JS:JE)
    case (iv3d_rhot)
      v3dg(:,:,:,iv3d) = RHOT(KS:KE,IS:IE,JS:JE)
    case (iv3d_q)
      v3dg(:,:,:,iv3d) = QV(KS:KE,IS:IE,JS:JE)
    case (iv3d_qc)
      v3dg(:,:,:,iv3d) = QC(KS:KE,IS:IE,JS:JE)
    case (iv3d_qr)
      v3dg(:,:,:,iv3d) = QR(KS:KE,IS:IE,JS:JE)
    case (iv3d_qi)
      v3dg(:,:,:,iv3d) = QI(KS:KE,IS:IE,JS:JE)
    case (iv3d_qs)
      v3dg(:,:,:,iv3d) = QS(KS:KE,IS:IE,JS:JE)
    case (iv3d_qg)
      v3dg(:,:,:,iv3d) = QG(KS:KE,IS:IE,JS:JE)
    case default
      write (6, '(3A)') "[Error] Variable '", trim(v3d_name(iv3d)), "' is not recognized."
      stop
    end select
  end do

  do iv2d = 1, nv2d
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 2D var [direct transfer]: ', trim(v2d_name(iv2d))
    end if
!    select case (iv2d)
!    case default
!      write (6, '(3A)') "[Error] Variable '", trim(v2d_name(iv2d)), "' is not recognized."
!      stop
!    end select
  end do

  return
end subroutine read_restart_direct

!-------------------------------------------------------------------------------
! [File I/O] Write SCALE restart files
!-------------------------------------------------------------------------------
!SUBROUTINE write_restart(filename,v3dg,v2dg)
!  use netcdf, only: NF90_WRITE

!!  use gtool_file, only: FileOpen, FileClose, FileWrite
!!  use gtool_file_h

!  use scale_prc, only: PRC_myrank
!  use common_mpi, only: myrank
!  use common_ncio
!  implicit none

!  CHARACTER(*),INTENT(IN) :: filename
!  REAL(RP),INTENT(IN) :: v3dg(nlev,nlon,nlat,nv3d)
!  REAL(RP),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
!  REAL(RP) :: rhotmp(nlev,nlon,nlat)
!  character(len=12) :: filesuffix = '.pe000000.nc'
!  integer :: iv3d,iv2d,ncid

!  REAL(RP) :: v3dgtmp(nlev,nlon,nlat)
!  REAL(RP) :: v2dgtmp(nlon,nlat)

!!  integer :: vid,error

!  WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',PRC_myrank,'.nc'

!  write (filesuffix(4:9),'(I6.6)') PRC_myrank
!  call ncio_open(filename // filesuffix, NF90_WRITE, ncid)

!!  call file_open(ncid, filename // filesuffix, File_FWRITE)
!!!  call FileOpen(ncid, filename, 2)

!!print *, '######### ncid:', ncid

!  DO iv3d = 1, nv3d
!    if (LOG_LEVEL >= 1) then
!      write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
!    end if
!    call ncio_write(ncid, trim(v3d_name(iv3d)), nlev, nlon, nlat, 1, v3dg(:,:,:,iv3d))
!    call ncio_read (ncid, trim(v3d_name(iv3d)), nlev, nlon, nlat, 1, v3dgtmp)
!    call ncio_write(ncid, trim(v3d_name(iv3d)), nlev, nlon, nlat, 1, v3dgtmp)
!  END DO

!  DO iv2d = 1, nv2d
!    if (LOG_LEVEL >= 1) then
!      write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
!    end if
!    call ncio_write(ncid, trim(v2d_name(iv2d)), nlon, nlat, 1, v2dg(:,:,iv2d))
!    call ncio_read (ncid, trim(v2d_name(iv2d)), nlon, nlat, 1, v2dgtmp)
!    call ncio_write(ncid, trim(v2d_name(iv2d)), nlon, nlat, 1, v2dgtmp)
!  END DO

!!  DO iv3d = 1, nv3d
!!    if (LOG_LEVEL >= 1) then
!!      write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
!!    end if
!!    call ncio_check(nf90_inq_varid(ncid, trim(v3d_name(iv3d)), vid))

!!!print *, '######### vid:', myrank, vid

!!!    call FileWrite(vid, v3dg(:,:,:,iv3d), -1.0d0, -1.0d0)
!!    call file_write_data( vid, v3dg(:,:,:,iv3d), -1.0d0, -1.0d0, RP, & ! (in)
!!         error                                     ) ! (out)
!!  END DO

!!  DO iv2d = 1, nv2d
!!    if (LOG_LEVEL >= 1) then
!!      write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
!!    end if
!!    call ncio_check(nf90_inq_varid(ncid, trim(v2d_name(iv2d)), vid))
!!!    call FileWrite(vid, v2dg(:,:,iv2d), 1.0d0, 1.0d0)
!!    call file_write_data( vid, v2dg(:,:,iv2d), -1.0d0, -1.0d0, RP, & ! (in)
!!         error                                     ) ! (out)
!!  END DO

!  call ncio_close(ncid)

!!  call file_close(ncid)
!!!  call FileClose(ncid)

!  RETURN
!END SUBROUTINE write_restart

SUBROUTINE write_restart(filename,v3dg,v2dg)
  use netcdf
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_HAS_W,  &
    PRC_HAS_S
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use common_mpi, only: myrank
  use common_ncio
  implicit none

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: iv3d, iv2d, ncid, varid
  integer :: is, js, ks

  is = 1
  js = 1
  if (.not. PRC_HAS_W) then
    is = is + IHALO
  end if
  if (.not. PRC_HAS_S) then
    js = js + JHALO
  end if

!  write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',PRC_myrank,'.nc'

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  if (LOG_LEVEL >= 5) then
    write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is writing a file ',trim(filename) // filesuffix
  endif
  call ncio_open(trim(filename) // filesuffix, NF90_WRITE, ncid)

  do iv3d = 1, nv3d
    if (LOG_LEVEL >= 5) then
      write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
    end if

    if (iv3d == iv3d_rhow) then
      ks = 2 ! ignore zh=1 (the surface where MOMZ = 0.0)
    else
      ks = 1
    endif

    call ncio_check(nf90_inq_varid(ncid, trim(v3d_name(iv3d)), varid))
    call ncio_check(nf90_put_var(ncid, varid, v3dg(:,:,:,iv3d), &
                                 start = (/ ks, is, js, 1 /),    &
                                 count = (/ KMAX, IMAX, JMAX, 1 /)))
  end do

  do iv2d = 1, nv2d
    if (LOG_LEVEL >= 5) then
      write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
    end if
    call ncio_check(nf90_inq_varid(ncid, trim(v2d_name(iv2d)), varid))
    call ncio_check(nf90_put_var(ncid, varid, v2dg(:,:,iv2d), &
                                 start = (/ is, js, 1 /),     &
                                 count = (/ IMAX, JMAX, 1 /)))
  end do

  call ncio_close(ncid)

  RETURN
END SUBROUTINE write_restart

!-------------------------------------------------------------------------------
! [File I/O] Write SCALE restart files <PnetCDF>
!-------------------------------------------------------------------------------
#ifdef PNETCDF
SUBROUTINE write_restart_par(filename,v3dg,v2dg,comm)
!  use common_mpi_scale, only: &
!    MPI_COMM_a
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_2Drank
!    PRC_PERIODIC_X, PRC_PERIODIC_Y, &
!    PRC_HAS_W,  &
!    PRC_HAS_E,  &
!    PRC_HAS_S,  &
!    PRC_HAS_N
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL
  use common_mpi, only: myrank
  use pnetcdf
  implicit none

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(INOUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(INOUT) :: v2dg(nlon,nlat,nv2d)
  integer,intent(in) :: comm
  integer :: iv3d,iv2d,ncid

  integer :: err, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(3), count(3)

  ! calculate subarray's start() and count() to the global variables
  start(1) = 1
  start(2) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(3) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  count(1) = KMAX
  count(2) = IMAX
  count(3) = JMAX
  start(2) = start(2) + IHALO
  start(3) = start(3) + JHALO
!  if (.NOT. PRC_PERIODIC_X) start(2) = start(2) + IHALO
!  if (.NOT. PRC_PERIODIC_Y) start(3) = start(3) + JHALO

  write (6,'(A,I6.6,3A,6I6)') 'MYRANK ',myrank,' is writing a file ',trim(filename)//'.nc', ' >> PnetCDF start(3), count(3) =', start, count

  err = nfmpi_open(comm, trim(filename)//".nc", NF_WRITE, MPI_INFO_NULL, ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//'.nc '//nfmpi_strerror(err)

  do iv3d = 1, nv3d
    if (LOG_LEVEL >= 5) then
      write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
    end if
    err = nfmpi_inq_varid(ncid, trim(v3d_name(iv3d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
#ifdef SINGLE
    err = nfmpi_iput_vara_real(ncid, varid, start, count, v3dg(:,:,:,iv3d), req)
#else
    err = nfmpi_iput_vara_double(ncid, varid, start, count, v3dg(:,:,:,iv3d), req)
#endif
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iput_vara_double '//' '//nfmpi_strerror(err)
  end do

  do iv2d = 1, nv2d
    if (LOG_LEVEL >= 5) then
      write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
    end if
    err = nfmpi_inq_varid(ncid, trim(v2d_name(iv2d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
#ifdef SINGLE
    err = nfmpi_iput_vara_real(ncid, varid, start(2:3), count(2:3), v2dg(:,:,iv2d), req)
#else
    err = nfmpi_iput_vara_double(ncid, varid, start(2:3), count(2:3), v2dg(:,:,iv2d), req)
#endif
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iput_vara_double '//' '//nfmpi_strerror(err)
  end do

  err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  RETURN
END SUBROUTINE write_restart_par
#endif

!-------------------------------------------------------------------------------
! [Direct transfer] Write SCALE restart files
!-------------------------------------------------------------------------------
subroutine write_restart_direct(v3dg,v2dg)
  use mod_atmos_vars, only: &
    DENS, &
    MOMX, &
    MOMY, &
    MOMZ, &
    RHOT, &
    QTRC, &
    QV, Qe, &
    ATMOS_vars_calc_diagnostics, &
    ATMOS_vars_fillhalo
  use mod_atmos_phy_mp_driver, only: &
    ATMOS_PHY_MP_driver_qhyd2qtrc
  use scale_atmos_hydrometeor, only: &
    I_QV, I_HC, I_HR, I_HI, I_HS, I_HG
  use mod_atmos_phy_mp_vars, only: &
    QS_MP, QE_MP
  use scale_atmos_grid_cartesC_index, only: &
    IS, IE, IA, JS, JE, JA, KS, KE, KA
  implicit none

  real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP), intent(in) :: v2dg(nlon,nlat,nv2d)
  integer :: iv3d, iv2d

  do iv3d = 1, nv3d
    if (LOG_LEVEL >= 5) then
      write(6,'(1x,A,A15)') '*** Write 3D var [direct transfer]: ', trim(v3d_name(iv3d))
    end if
    select case (iv3d)
    case (iv3d_rho)
      DENS(KS:KE,IS:IE,JS:JE) = v3dg(:,:,:,iv3d)
    case (iv3d_rhou)
      MOMX(KS:KE,IS:IE,JS:JE) = v3dg(:,:,:,iv3d)
    case (iv3d_rhov)
      MOMY(KS:KE,IS:IE,JS:JE) = v3dg(:,:,:,iv3d)
    case (iv3d_rhow)
      MOMZ(KS:KE,IS:IE,JS:JE) = v3dg(:,:,:,iv3d)
    case (iv3d_rhot)
      RHOT(KS:KE,IS:IE,JS:JE) = v3dg(:,:,:,iv3d)
    case (iv3d_q)
      QV(KS:KE,IS:IE,JS:JE) = v3dg(:,:,:,iv3d)
    case (iv3d_qc)
      Qe(KS:KE,IS:IE,JS:JE,I_HC) = v3dg(:,:,:,iv3d)
    case (iv3d_qr)
      Qe(KS:KE,IS:IE,JS:JE,I_HR) = v3dg(:,:,:,iv3d)
    case (iv3d_qi)
      Qe(KS:KE,IS:IE,JS:JE,I_HI) = v3dg(:,:,:,iv3d)
    case (iv3d_qs)
      Qe(KS:KE,IS:IE,JS:JE,I_HS) = v3dg(:,:,:,iv3d)
    case (iv3d_qg)
      Qe(KS:KE,IS:IE,JS:JE,I_HG) = v3dg(:,:,:,iv3d)
    case default
      write (6, '(3A)') "[Error] Variable '", trim(v3d_name(iv3d)), "' is not recognized."
      stop
    end select
  end do

  ! Assume Tomita08
  call ATMOS_PHY_MP_driver_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, & 
                                      QV, Qe, &                    ! [IN]
                                      QTRC(:,:,:,QS_MP:QE_MP)    ) ! [OUT] 
  call ATMOS_vars_fillhalo 

  call ATMOS_vars_calc_diagnostics 

  do iv2d = 1, nv2d
    if (LOG_LEVEL >= 5) then
      write(6,'(1x,A,A15)') '*** Write 2D var [direct transfer]: ', trim(v2d_name(iv2d))
    end if
!    select case (iv2d)
!    case default
!      write (6, '(3A)') "[Error] Variable '", trim(v2d_name(iv2d)), "' is not recognized."
!      stop
!    end select
  end do

  return
end subroutine write_restart_direct

subroutine write_dafcst_nc( filename, step, nlev_plot, ref3d ) 
  use netcdf
  use common_ncio
  use scale_const, only: &
    D2R => CONST_D2R
  use scale_atmos_grid_cartesC, only: &
    CZ => ATMOS_GRID_CARTESC_CZ, &
    CXG => ATMOS_GRID_CARTESC_CXG, &
    CYG => ATMOS_GRID_CARTESC_CYG
  use scale_atmos_grid_cartesC_index, only: &
    KHALO, IHALO, JHALO, KMAX
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  implicit none

  character(*), intent(in) :: filename
  integer, intent(in) :: step
  integer, intent(in) :: nlev_plot
  real(r_sngl) :: ref3d(nlev_plot, nlong, nlatg)
  integer :: ncid, varid

  integer :: tlev

  character(len=*), parameter :: lat_name = "Latitude"
  character(len=*), parameter :: lon_name = "Longitude"
  character(len=*), parameter :: z_name   = "Height"
  character(len=*), parameter :: t_name = "time"
  character(len=*), parameter :: ref_name = "Reflectivity"
  character(len=*), parameter :: lon_unit = "degree_east"
  character(len=*), parameter :: lat_unit = "degree_north"
  character(len=*), parameter :: z_unit = "m (above ground level)"
  character(len=*), parameter :: t_unit = "seconds"
  character(len=*), parameter :: ref_unit = "dBZ"
  integer :: z_dimid, lon_dimid, lat_dimid, t_dimid
  integer :: lon_varid, lat_varid, z_varid, t_varid
  integer :: ref_varid

  integer :: dimids(4)
  integer :: start(4), count(4)

  real :: lons(nlong), lats(nlatg), zlevs(nlev_plot)
  real(RP) :: lon_RP, lat_RP
  real, allocatable :: tlevs(:)
  integer :: i, j, k

  if ( step == 0 ) then

    tlev = int( DACYCLE_RUN_FCST_TIME / OUT_DAFCST_DSEC ) + 1
    allocate( tlevs(tlev) )

    ! Create the file. 
    call ncio_create( trim(filename), nf90_clobber, ncid )

    ! Define the dimensions. 
    call ncio_check( nf90_def_dim( ncid, trim(lon_name), nlong, lon_dimid) )
    call ncio_check( nf90_def_dim( ncid, trim(lat_name), nlatg, lat_dimid) )
    call ncio_check( nf90_def_dim( ncid, trim(z_name), nlev_plot, z_dimid) )
    call ncio_check( nf90_def_dim( ncid, trim(t_name), tlev, t_dimid) )

    ! Define the coordinate variables. 
    call ncio_check( nf90_def_var(ncid, trim(lon_name), nf90_real, lon_dimid, lon_varid) )
    call ncio_check( nf90_def_var(ncid, trim(lat_name), nf90_real, lat_dimid, lat_varid) )
    call ncio_check( nf90_def_var(ncid, trim(z_name), nf90_real, z_dimid, z_varid) )
    call ncio_check( nf90_def_var(ncid, trim(t_name), nf90_real, t_dimid, t_varid) )

    ! Assume Mercator projection
    do i = 1, nlong
      call MAPPROJECTION_xy2lonlat( CXG(IHALO+i), CYG(JHALO+1), lon_RP, lat_RP )
      lons(i) = real( lon_RP/D2R, kind=r_sngl )
    enddo

    do j = 1, nlatg
      call MAPPROJECTION_xy2lonlat( CXG(IHALO+1), CYG(JHALO+j), lon_RP, lat_RP )
      lats(j) = real( lat_RP/D2R, kind=r_sngl )
    enddo

    do i = 1, tlev
      tlevs(i) = real( i - 1 ) * real( OUT_DAFCST_DSEC, kind=r_sngl )
    enddo


    ! Assign units attributes to coordinate variables.
    call ncio_check( nf90_put_att(ncid, lon_varid, "units", trim(lon_unit) ) )
    call ncio_check( nf90_put_att(ncid, lat_varid, "units", trim(lat_unit) ) )
    call ncio_check( nf90_put_att(ncid, z_varid, "units", trim(z_unit) ) )
    call ncio_check( nf90_put_att(ncid, t_varid, "units", trim(t_unit) ) )
    dimids = (/ z_dimid, lon_dimid, lat_dimid, t_dimid /)

    ! Define the netCDF variables
    call ncio_check( nf90_def_var(ncid, trim(ref_name), nf90_real, dimids, ref_varid) )
    call ncio_check( nf90_put_att(ncid, ref_varid, "units", trim(ref_unit) ) )

    ! End define mode
    call ncio_check( nf90_enddef(ncid) )

    ! Write the coordinate variable data. 
    call ncio_check( nf90_put_var(ncid, lat_varid, lats) )
    call ncio_check( nf90_put_var(ncid, lon_varid, lons) )
    call ncio_check( nf90_put_var(ncid, t_varid, tlevs) )

    j = 0
    do k = max( 1, OUT_NETCDF_ZLEV_MIN), min( nlev, OUT_NETCDF_ZLEV_MAX), OUT_NETCDF_ZLEV_INTV
       j = j + 1
       zlevs(j) = real( CZ(KHALO+k), kind=r_sngl )

    enddo

    call ncio_check( nf90_put_var(ncid, z_varid, zlevs ) )

    deallocate( tlevs )

  else

    call ncio_open( trim(filename), nf90_write, ncid )
    call ncio_check( nf90_inq_varid( ncid, trim(ref_name), ref_varid) )

  endif

  count = (/ nlev_plot, nlong, nlatg, 1 /)
  start = (/ 1, 1, 1, step+1 /)

  ! Write the data.
  call ncio_check( nf90_put_var(ncid, ref_varid, ref3d, start = start, &
                   count = count) )

  call ncio_close( ncid )

  return
end subroutine write_dafcst_nc

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE restart files for model coordinates
!-------------------------------------------------------------------------------
SUBROUTINE read_restart_coor(filename,lon,lat,height)
  use netcdf
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_HAS_W,  &
    PRC_HAS_S
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use common_mpi, only: myrank
  use common_ncio
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: lon(nlon,nlat)
  REAL(RP),INTENT(OUT) :: lat(nlon,nlat)
  REAL(RP),INTENT(OUT) :: height(nlev,nlon,nlat)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: ncid, varid
  integer :: is, js

  is = 1
  js = 1
  if (.not. PRC_HAS_W) then
    is = is + IHALO
  end if
  if (.not. PRC_HAS_S) then
    js = js + JHALO
  end if

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
!  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)

!!! restart files do not contain 3D height variable before SCALE v5.1
!  if (LOG_LEVEL >= 1) then
!    write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(height3d_name)
!  end if
!  call ncio_check(nf90_inq_varid(ncid, trim(height3d_name), varid))
!  call ncio_check(nf90_get_var(ncid, varid, height,        &
!                               start = (/ 1, is, js, 1 /), &
!                               count = (/ KMAX, IMAX, JMAX, 1 /)))

  if (LOG_LEVEL >= 4) then
    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(lon2d_name)
  end if
  call ncio_check(nf90_inq_varid(ncid, trim(lon2d_name), varid))
  call ncio_check(nf90_get_var(ncid, varid, lon,        &
                               start = (/ is, js, 1 /), &
                               count = (/ IMAX, JMAX, 1 /)))

  if (LOG_LEVEL >= 4) then
    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(lat2d_name)
  end if
  call ncio_check(nf90_inq_varid(ncid, trim(lat2d_name), varid))
  call ncio_check(nf90_get_var(ncid, varid, lat,        &
                               start = (/ is, js, 1 /), &
                               count = (/ IMAX, JMAX, 1 /)))

  call ncio_close(ncid)

  RETURN
END SUBROUTINE read_restart_coor

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE topography files
!-------------------------------------------------------------------------------
SUBROUTINE read_topo(filename,topo)
  use netcdf
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_HAS_W,  &
    PRC_HAS_S
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX
  use common_mpi, only: myrank
  use common_ncio
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: topo(nlon,nlat)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: ncid, varid
  integer :: is, js

  is = 1
  js = 1
  if (.not. PRC_HAS_W) then
    is = is + IHALO
  end if
  if (.not. PRC_HAS_S) then
    js = js + JHALO
  end if

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  if (LOG_LEVEL >= 3) then
    write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix
  endif
  call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)

  if (LOG_LEVEL >= 4) then
    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(topo2d_name)
  end if
  call ncio_check(nf90_inq_varid(ncid, trim(topo2d_name), varid))
  call ncio_check(nf90_get_var(ncid, varid, topo,       &
                               start = (/ is, js, 1 /), &
                               count = (/ IMAX, JMAX, 1 /)))

  call ncio_close(ncid)

  RETURN
END SUBROUTINE read_topo

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE topography files <PnetCDF>
!-------------------------------------------------------------------------------
#ifdef PNETCDF
SUBROUTINE read_topo_par(filename,topo,comm)
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_2Drank
!    PRC_PERIODIC_X, PRC_PERIODIC_Y
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL
  use common_mpi, only: myrank
  use pnetcdf
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: topo(nlon,nlat)
  real(RP) :: topo_RP(nlon,nlat)
  integer,intent(in) :: comm

  integer :: ncid

  integer :: err, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(2), count(2)

  ! calculate subarray's start() and count() to the global variables
  start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  count(1) = IMAX
  count(2) = JMAX
  start(1) = start(1) + IHALO
  start(2) = start(2) + JHALO
!  if (.NOT. PRC_PERIODIC_X) start(2) = start(2) + IHALO
!  if (.NOT. PRC_PERIODIC_Y) start(3) = start(3) + JHALO

!  write (6,'(A,I6.6,3A,4I6)') 'MYRANK ',myrank,' is reading a file ',trim(filename)//'.nc', ' >> PnetCDF start(2), count(2) =', start, count

  err = nfmpi_open(comm, trim(filename)//".nc", NF_NOWRITE, MPI_INFO_NULL, ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//'.nc '//nfmpi_strerror(err)

  if (LOG_LEVEL >= 4) then
    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(topo2d_name)
  end if
  err = nfmpi_inq_varid(ncid, trim(topo2d_name), varid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
#ifdef SINGLE
  err = nfmpi_iget_vara_real(ncid, varid, start, count, topo_RP, req)
#else
  err = nfmpi_iget_vara_double(ncid, varid, start, count, topo_RP, req)
#endif
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_get_vara_double_all '//' '//nfmpi_strerror(err)

  err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  topo = real(topo_RP, kind=r_size)

  RETURN
END SUBROUTINE read_topo_par
#endif

!-------------------------------------------------------------------------------
! [Direct transfer] Read SCALE topography files
!-------------------------------------------------------------------------------
subroutine read_topo_direct(topo)
  use scale_topography, only: &
    TOPO_Zsfc
  use scale_atmos_grid_cartesC_index, only: &
    IS, IE, JS, JE
  implicit none

  real(r_size), intent(out) :: topo(nlon,nlat)

  if (LOG_LEVEL >= 4) then
    write(6,'(1x,A,A15)') '*** Read 2D var [direct transfer]: ', trim(topo2d_name)
  end if
  topo(:,:) = TOPO_Zsfc(IS:IE,JS:JE)

  return
end subroutine read_topo_direct

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE history files
!-------------------------------------------------------------------------------
subroutine read_history(filename,step,v3dg,v2dg)
  use scale_prc, only: &
      PRC_myrank
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA
  use scale_file, only: &
      FILE_read
  use scale_comm_cartesC, only: &
      COMM_vars8, &
      COMM_wait
  use common_mpi, only: myrank
  implicit none

  character(*),intent(in) :: filename
  integer,intent(in) :: step
  real(r_size),intent(out) :: v3dg(nlevh,nlonh,nlath,nv3dd)
  real(r_size),intent(out) :: v2dg(nlonh,nlath,nv2dd)
  real(RP) :: v3dg_RP(nlevh,nlonh,nlath,nv3dd)
  real(RP) :: v2dg_RP(nlonh,nlath,nv2dd)
  integer :: i,j,k,iv3d,iv2d
  character(len=12) :: filesuffix = '.pe000000.nc'
  real(RP) :: var3D(nlon,nlat,nlev)
  real(RP) :: var2D(nlon,nlat)

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
!  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix

!  comm = MPI_COMM_NULL
!  if ( FILE_AGGREGATE ) then
!    comm = PRC_LOCAL_COMM_WORLD
!  end if
!  call FILE_open( basename,           & ! [IN]
!                  fid,                ) ! [OUT]
!!                  mpi_comm = comm,    ) ! [IN]

  ! 3D variables
  !-------------
  do iv3d = 1, nv3dd
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 3D hist var: ', trim(v3dd_name(iv3d))
    end if
    if (v3dd_hastime(iv3d)) then
      call FILE_read( filename,              & ! [IN]
                      trim(v3dd_name(iv3d)), & ! [IN]
                      var3D,                 & ! [OUT]
                      step=step              ) ! [IN]
    else
      call FILE_read( filename,              & ! [IN]
                      trim(v3dd_name(iv3d)), & ! [IN]
                      var3D                  ) ! [OUT]
    end if
    forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg(k+KHALO,i+IHALO,j+JHALO,iv3d) = var3D(i,j,k) ! use FORALL to change order of dimensions
  end do

  ! 2D variables
  !-------------
  do iv2d = 1, nv2dd
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 2D hist var: ', trim(v2dd_name(iv2d))
    end if
    if (v2dd_hastime(iv2d)) then
      call FILE_read( filename,              & ! [IN]
                      trim(v2dd_name(iv2d)), & ! [IN]
                      var2D,                 & ! [OUT]
                      step=step              ) ! [IN]
    else
      call FILE_read( filename,              & ! [IN]
                      trim(v2dd_name(iv2d)), & ! [IN]
                      var2D                  ) ! [OUT]
    end if
    v2dg(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2d) = var2D(:,:)
  end do

  ! Communicate halo
  !-------------
!$OMP PARALLEL DO PRIVATE(i,j,iv3d) SCHEDULE(STATIC) COLLAPSE(2)
  do iv3d = 1, nv3dd
    do j = JS, JE
      do i = IS, IE
        v3dg_RP(   1:KS-1,i,j,iv3d) = v3dg_RP(KS,i,j,iv3d)
        v3dg_RP(KE+1:KA,  i,j,iv3d) = v3dg_RP(KE,i,j,iv3d)
      end do
    end do
  end do
!$OMP END PARALLEL DO

  do iv3d = 1, nv3dd
    call COMM_vars8( v3dg_RP(:,:,:,iv3d), iv3d )
  end do
  do iv3d = 1, nv3dd
    call COMM_wait ( v3dg_RP(:,:,:,iv3d), iv3d )
  end do

  do iv2d = 1, nv2dd
    call COMM_vars8( v2dg_RP(:,:,iv2d), iv2d )
  end do
  do iv2d = 1, nv2dd
    call COMM_wait ( v2dg_RP(:,:,iv2d), iv2d )
  end do

  v3dg = real(v3dg_RP, kind=r_size)
  v2dg = real(v2dg_RP, kind=r_size)

  ! Save topo for later use
  !-------------
  if (.not. allocated(topo2d)) then
    allocate (topo2d(nlon,nlat))
    topo2d = v2dg(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_topo)
  end if

  return
end subroutine read_history

!-------------------------------------------------------------------------------
! [File I/O] Read SCALE history files <PnetCDF>
!-------------------------------------------------------------------------------
#ifdef PNETCDF
subroutine read_history_par(filename,step,v3dg,v2dg,comm)
  use scale_prc, only: &
      PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_2Drank
!    PRC_PERIODIC_X, PRC_PERIODIC_Y
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA, &
      IMAX, JMAX, KMAX
!  use gtool_history, only: &
!      HistoryGet
  use scale_comm_cartesC, only: &
      COMM_vars8, &
      COMM_wait
  use mpi, only: MPI_OFFSET_KIND, MPI_INFO_NULL
  use common_mpi, only: myrank
  use pnetcdf
  implicit none

  character(*),intent(in) :: filename
  integer,intent(in) :: step
  real(r_size),intent(out) :: v3dg(nlevh,nlonh,nlath,nv3dd)
  real(r_size),intent(out) :: v2dg(nlonh,nlath,nv2dd)
  real(RP) :: v3dg_RP(nlevh,nlonh,nlath,nv3dd)
  real(RP) :: v2dg_RP(nlonh,nlath,nv2dd)
  integer,intent(in) :: comm
  integer :: i,j,k,iv3d,iv2d
  real(SP) :: var3D(nlon,nlat,nlev,nv3dd)
  real(SP) :: var2D(nlon,nlat,nv2dd)

!  integer :: fid
  integer :: ncid

  integer :: err, varid, req, reqs(1), sts(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(4), count(4)

  ! calculate subarray's start() and count() to the global variables
  start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1
  start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1
  start(3) = 1
  start(4) = step
  count(1) = IMAX
  count(2) = JMAX
  count(3) = KMAX
  count(4) = 1
!  start(1) = start(1) + IHALO   ! History files always have no halo
!  start(2) = start(2) + JHALO   !
!  if (.NOT. PRC_PERIODIC_X) start(1) = start(1) + IHALO
!  if (.NOT. PRC_PERIODIC_Y) start(2) = start(2) + JHALO

!  write (6,'(A,I6.6,3A,8I6)') 'MYRANK ',myrank,' is reading a file ',trim(filename)//'.nc', ' >> PnetCDF start(4), count(4) =', start, count

!  call FILEIO_open( fid, trim(filename) )
  err = nfmpi_open(comm, trim(filename)//".nc", NF_NOWRITE, MPI_INFO_NULL, ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_open '//trim(filename)//'.nc '//nfmpi_strerror(err)

  ! 3D variables
  !-------------
  do iv3d = 1, nv3dd
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 3D hist var: ', trim(v3dd_name(iv3d))
    end if

!--- neither of these work now ---
!    call FILEIO_read( var3D(:,:,:,iv3d),    & ! [OUT]
!                      fid, trim(v3dd_name(iv3d)), 'XYZ', step=step ) ! [IN]   !!! 'XYZ' is not supported.
!    call HistoryGet( var3D(:,:,:,iv3d),     & ! [OUT]                         !!! 'HistoryGet' does not support PNETCDF
!                     filename,              & ! [IN]
!                     trim(v3dd_name(iv3d)), & ! [IN]
!                     step                   ) ! [IN]
    err = nfmpi_inq_varid(ncid, trim(v3dd_name(iv3d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    if (v3dd_hastime(iv3d)) then
      err = nfmpi_iget_vara_real(ncid, varid, start, count, var3D(:,:,:,iv3d), req)
    else
      err = nfmpi_iget_vara_real(ncid, varid, start(1:3), count(1:3), var3D(:,:,:,iv3d), req)
    end if
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)

!!    call FILEIO_flush( fid )
!    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
!    if ( err .NE. NF_NOERR ) &
!       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg_RP(k+KHALO,i+IHALO,j+JHALO,iv3d) = real(var3D(i,j,k,iv3d), r_size) ! use FORALL to change order of dimensions
  end do

  start(3) = step
  count(3) = 1

  ! 2D variables
  !-------------
  do iv2d = 1, nv2dd
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 2D hist var: ', trim(v2dd_name(iv2d))
    end if

!--- neither of these work now ---
!    call FILEIO_read( var2D(:,:,iv2d),      & ! [OUT]
!                      fid, trim(v2dd_name(iv2d)), 'XY', step=step ) ! [IN]
!    call HistoryGet( var2D(:,:,iv2d),       & ! [OUT]                        !!! 'HistoryGet' does not support PNETCDF
!                     filename,              & ! [IN]
!                     trim(v2dd_name(iv2d)), & ! [IN]
!                     step                   ) ! [IN]

    err = nfmpi_inq_varid(ncid, trim(v2dd_name(iv2d)), varid)
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_inq_varid '//' '//nfmpi_strerror(err)
    if (v2dd_hastime(iv2d)) then
      err = nfmpi_iget_vara_real(ncid, varid, start(1:3), count(1:3), var2D(:,:,iv2d), req)
    else
      err = nfmpi_iget_vara_real(ncid, varid, start(1:2), count(1:2), var2D(:,:,iv2d), req)
    end if
    if ( err .NE. NF_NOERR ) &
       write (6,'(A)') 'failed nfmpi_iget_vara_real '//' '//nfmpi_strerror(err)

!!    call FILEIO_flush( fid )
!    err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
!    if ( err .NE. NF_NOERR ) &
!       write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!    v2dg_RP(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2d) = real(var2D(:,:,iv2d), r_size)
  end do

!  call FILEIO_flush( fid )
  err = nfmpi_wait_all(ncid, NF_REQ_ALL, reqs, sts)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_wait_all '//' '//nfmpi_strerror(err)

!  call FILEIO_close( fid )
  err = nfmpi_close(ncid)
  if ( err .NE. NF_NOERR ) &
     write (6,'(A)') 'failed nfmpi_close '//' '//nfmpi_strerror(err)

  ! Copy data buffer
  !-------------
!$OMP PARALLEL PRIVATE(i,j,k,iv3d,iv2d)
!$OMP DO SCHEDULE(STATIC)
  do iv3d = 1, nv3dd
    forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg_RP(k+KHALO,i+IHALO,j+JHALO,iv3d) = real(var3D(i,j,k,iv3d), r_size) ! use FORALL to change order of dimensions
  end do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  do iv2d = 1, nv2dd
    v2dg_RP(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2d) = real(var2D(:,:,iv2d), r_size)
  end do
!$OMP END DO
!$OMP END PARALLEL

  ! Communicate halo
  !-------------
!$OMP PARALLEL DO PRIVATE(i,j,iv3d) SCHEDULE(STATIC) COLLAPSE(2)
  do iv3d = 1, nv3dd
    do j = JS, JE
      do i = IS, IE
        v3dg_RP(   1:KS-1,i,j,iv3d) = v3dg_RP(KS,i,j,iv3d)
        v3dg_RP(KE+1:KA,  i,j,iv3d) = v3dg_RP(KE,i,j,iv3d)
      end do
    end do
  end do
!$OMP END PARALLEL DO

  do iv3d = 1, nv3dd
    call COMM_vars8( v3dg_RP(:,:,:,iv3d), iv3d )
  end do
  do iv3d = 1, nv3dd
    call COMM_wait ( v3dg_RP(:,:,:,iv3d), iv3d )
  end do

  do iv2d = 1, nv2dd
    call COMM_vars8( v2dg_RP(:,:,iv2d), iv2d )
  end do
  do iv2d = 1, nv2dd
    call COMM_wait ( v2dg_RP(:,:,iv2d), iv2d )
  end do

  v3dg = real(v3dg_RP, kind=r_size)
  v2dg = real(v2dg_RP, kind=r_size)

  ! Save topo for later use
  !-------------
  if (.not. allocated(topo2d)) then
    allocate (topo2d(nlon,nlat))
    topo2d = v2dg(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_topo)
  end if

  return
end subroutine read_history_par
#endif

!-------------------------------------------------------------------------------
! [Direct transfer] Read SCALE history files
!-------------------------------------------------------------------------------
subroutine read_history_direct(v3dg, v2dg)
  use mod_atmos_vars, only: &
    ATMOS_vars_calc_diagnostics, &
    ATMOS_vars_get_diagnostic, &
    QV, QC, QR, &
    QI, QS, QG, &
    U, V, W, &
    PRES, TEMP
  use mod_atmos_phy_sf_vars, only: &
    ATMOS_PHY_SF_SFC_PRES, &
    ATMOS_PHY_SF_U10, &
    ATMOS_PHY_SF_V10, &
    ATMOS_PHY_SF_T2, &
    ATMOS_PHY_SF_Q2
  use scale_topography, only: &
    TOPO_Zsfc
  use scale_atmos_grid_cartesC_real, only: &
    ATMOS_GRID_CARTESC_REAL_CZ
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IS, IE, JS, JE, KS, KE, KA
  use scale_comm_cartesC, only: &
    COMM_vars8, &
    COMM_wait
  implicit none

  real(r_size),intent(out) :: v3dg(nlevh,nlonh,nlath,nv3dd)
  real(r_size),intent(out) :: v2dg(nlonh,nlath,nv2dd)
  real(RP) :: v3dg_RP(nlevh,nlonh,nlath,nv3dd)
  real(RP) :: v2dg_RP(nlonh,nlath,nv2dd)
  integer :: i, j, iv3d, iv2d

  ! 3D variables
  !-------------

  call ATMOS_vars_calc_diagnostics 

  do iv3d = 1, nv3dd
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 3D hist var [direct transfer]: ', trim(v3dd_name(iv3d))
    end if
    select case (iv3d)
    case (iv3dd_u)
      v3dg_RP(:,:,:,iv3d) = U(:,:,:)
    case (iv3dd_v)
      v3dg_RP(:,:,:,iv3d) = V(:,:,:)
    case (iv3dd_w)
      v3dg_RP(:,:,:,iv3d) = W(:,:,:)
    case (iv3dd_t)
      v3dg_RP(:,:,:,iv3d) = TEMP(:,:,:)
    case (iv3dd_rh)
      ! RH relative to liquid
      ! Not used as of 12/11/2019
      call ATMOS_vars_get_diagnostic(trim(v3dd_name(iv3d)), v3dg_RP(:,:,:,iv3d))
    case (iv3dd_p)
      v3dg_RP(:,:,:,iv3d) = PRES(:,:,:)
    case (iv3d_q)
      v3dg_RP(:,:,:,iv3d) = QV(:,:,:)
    case (iv3d_qc)
      v3dg_RP(:,:,:,iv3d) = QC(:,:,:)
    case (iv3d_qr)
      v3dg_RP(:,:,:,iv3d) = QR(:,:,:)
    case (iv3d_qi)
      v3dg_RP(:,:,:,iv3d) = QI(:,:,:)
    case (iv3d_qs)
      v3dg_RP(:,:,:,iv3d) = QS(:,:,:)
    case (iv3d_qg)
      v3dg_RP(:,:,:,iv3d) = QG(:,:,:)
    case (iv3dd_hgt)
      v3dg_RP(:,:,:,iv3d) = ATMOS_GRID_CARTESC_REAL_CZ(:,:,:)
    case default
      write (6, '(3A)') "[Error] Variable '", trim(v3dd_name(iv3d)), "' is not recognized."
      stop
    end select
  end do

  ! 2D variables
  !-------------
  do iv2d = 1, nv2dd
    if (LOG_LEVEL >= 4) then
      write(6,'(1x,A,A15)') '*** Read 2D hist var [direct transfer]: ', trim(v2dd_name(iv2d))
    end if
    select case (iv2d)
    case (iv2dd_rain)
      call ATMOS_vars_get_diagnostic(trim(v2dd_name(iv2d)), v2dg_RP(:,:,iv2d))
    case (iv2dd_topo)
      v2dg_RP(:,:,iv2d) = TOPO_Zsfc(:,:)
    case (iv2dd_ps)
      v2dg_RP(:,:,iv2d) = ATMOS_PHY_SF_SFC_PRES(:,:)
    case (iv2dd_u10m)
      v2dg_RP(:,:,iv2d) = ATMOS_PHY_SF_U10(:,:)
    case (iv2dd_v10m)
      v2dg_RP(:,:,iv2d) = ATMOS_PHY_SF_V10(:,:)
    case (iv2dd_t2m)
      v2dg_RP(:,:,iv2d) = ATMOS_PHY_SF_T2(:,:)
    case (iv2dd_q2m)
      v2dg_RP(:,:,iv2d) = ATMOS_PHY_SF_Q2(:,:)
    case default
      write (6, '(3A)') "[Error] Variable '", trim(v2dd_name(iv2d)), "' is not recognized."
      stop
    end select
  end do

  ! Communicate halo
  !-------------
!$OMP PARALLEL DO PRIVATE(i,j,iv3d) SCHEDULE(STATIC) COLLAPSE(2)
  do iv3d = 1, nv3dd
    do j = JS, JE
      do i = IS, IE
        v3dg_RP(   1:KS-1,i,j,iv3d) = v3dg_RP(KS,i,j,iv3d)
        v3dg_RP(KE+1:KA,  i,j,iv3d) = v3dg_RP(KE,i,j,iv3d)
      end do
    end do
  end do
!$OMP END PARALLEL DO

  do iv3d = 1, nv3dd
    call COMM_vars8( v3dg_RP(:,:,:,iv3d), iv3d )
  end do
  do iv3d = 1, nv3dd
    call COMM_wait ( v3dg_RP(:,:,:,iv3d), iv3d )
  end do

  do iv2d = 1, nv2dd
    call COMM_vars8( v2dg_RP(:,:,iv2d), iv2d )
  end do
  do iv2d = 1, nv2dd
    call COMM_wait ( v2dg_RP(:,:,iv2d), iv2d )
  end do

  v3dg = real(v3dg_RP, kind=r_size)
  v2dg = real(v2dg_RP, kind=r_size)

  ! Save topo for later use
  !-------------
  if (.not. allocated(topo2d)) then
    allocate (topo2d(nlon,nlat))
    topo2d = v2dg(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_topo)
  end if

  return
end subroutine read_history_direct

!-------------------------------------------------------------------------------
! Transform the SCALE restart variables to the LETKF state variables
!-------------------------------------------------------------------------------
subroutine state_trans(v3dg)
  use scale_tracer, only: TRACER_CV
    use scale_const, only: &
       Rdry  => CONST_Rdry, &
       Rvap  => CONST_Rvap, &
       CVdry => CONST_CVdry, &
       PRE00 => CONST_PRE00
  implicit none

  real(RP), intent(inout) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: rho,pres,temp
  real(RP) :: qdry,CVtot,Rtot,CPovCV
  integer :: i,j,k,iv3d

!$OMP PARALLEL DO PRIVATE(i,j,k,iv3d,qdry,CVtot,Rtot,CPovCV,rho,pres,temp) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev
       qdry  = 1.0d0
       CVtot = 0.0d0
       do iv3d = iv3d_q, nv3d ! loop over all moisture variables
         qdry  = qdry - v3dg(k,i,j,iv3d)
         CVtot = CVtot + v3dg(k,i,j,iv3d) * TRACER_CV(iv3d-iv3d_q+1)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * v3dg(k,i,j,iv3d_q)
       CPovCV = ( CVtot + Rtot ) / CVtot

       rho = v3dg(k,i,j,iv3d_rho)
       pres = PRE00 * ( v3dg(k,i,j,iv3d_rhot) * Rtot / PRE00 )**CPovCV
       temp = pres / ( rho * Rtot )

       v3dg(k,i,j,iv3d_u) = v3dg(k,i,j,iv3d_rhou) / rho !!!!!! inaccurate! do not consider staggered grid !!!!!!
       v3dg(k,i,j,iv3d_v) = v3dg(k,i,j,iv3d_rhov) / rho !!!!!!
       v3dg(k,i,j,iv3d_w) = v3dg(k,i,j,iv3d_rhow) / rho !!!!!!
       v3dg(k,i,j,iv3d_t) = temp
       v3dg(k,i,j,iv3d_p) = pres
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

  return
end subroutine state_trans

!-------------------------------------------------------------------------------
! Inversely transform the LETKF state variables to the SCALE restart variables
!-------------------------------------------------------------------------------
subroutine state_trans_inv(v3dg)
  use scale_tracer, only: TRACER_CV
    use scale_const, only: &
       Rdry  => CONST_Rdry, &
       Rvap  => CONST_Rvap, &
       CVdry => CONST_CVdry, &
       PRE00 => CONST_PRE00
  implicit none

  real(RP), intent(inout) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: rho,rhot
  real(RP) :: qdry,CVtot,Rtot,CVovCP
  integer :: i,j,k,iv3d

  do iv3d = 1, nv3d
    if (POSITIVE_DEFINITE_Q .and. iv3d == iv3d_q) then
      v3dg(:,:,:,iv3d) = max(v3dg(:,:,:,iv3d), 0.0d0)
    else if (POSITIVE_DEFINITE_QHYD .and. &
             (iv3d == iv3d_qc .or. iv3d == iv3d_qr .or. iv3d == iv3d_qi .or. iv3d == iv3d_qs .or. iv3d == iv3d_qg)) then
      v3dg(:,:,:,iv3d) = max(v3dg(:,:,:,iv3d), 0.0d0)
    end if
  end do

!$OMP PARALLEL DO PRIVATE(i,j,k,iv3d,qdry,CVtot,Rtot,CVovCP,rho,rhot) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev
       qdry  = 1.0d0
       CVtot = 0.0d0
       do iv3d = iv3d_q, nv3d ! loop over all moisture variables
         qdry  = qdry - v3dg(k,i,j,iv3d)
         CVtot = CVtot + v3dg(k,i,j,iv3d) * TRACER_CV(iv3d-iv3d_q+1)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * v3dg(k,i,j,iv3d_q)
       CVovCP = CVtot / ( CVtot + Rtot )

       rho = v3dg(k,i,j,iv3d_p) / (Rtot * v3dg(k,i,j,iv3d_t))
       rhot = PRE00 / Rtot * (v3dg(k,i,j,iv3d_p) / PRE00)**CVovCP

       v3dg(k,i,j,iv3d_rhot) = rhot
       v3dg(k,i,j,iv3d_rhow) = v3dg(k,i,j,iv3d_w) * rho !!!!!! inaccurate! do not consider staggered grid !!!!!!
       v3dg(k,i,j,iv3d_rhov) = v3dg(k,i,j,iv3d_v) * rho !!!!!!
       v3dg(k,i,j,iv3d_rhou) = v3dg(k,i,j,iv3d_u) * rho !!!!!!
       v3dg(k,i,j,iv3d_rho) = rho
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

  return
end subroutine state_trans_inv

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
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA
  use scale_comm_cartesC, only: &
      COMM_vars8, &
      COMM_wait
  implicit none

  real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP), intent(in) :: v2dg(nlon,nlat,nv2d)
  real(r_size), intent(in) :: topo(nlon,nlat)
  real(r_size), intent(out) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  real(r_size), intent(out) :: v2dgh(nlonh,nlath,nv2dd)
  real(RP) :: v3dgh_RP(nlevh,nlonh,nlath,nv3dd)
  real(RP) :: v2dgh_RP(nlonh,nlath,nv2dd)

  real(r_size) :: height(nlev,nlon,nlat)
  integer :: i, j, k, iv3d, iv2d

  ! Variables that can be directly copied
  !---------------------------------------------------------

  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_u) = v3dg(:,:,:,iv3d_u)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_v) = v3dg(:,:,:,iv3d_v)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_w) = v3dg(:,:,:,iv3d_w)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_t) = v3dg(:,:,:,iv3d_t)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_p) = v3dg(:,:,:,iv3d_p)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_q) = v3dg(:,:,:,iv3d_q)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qc) = v3dg(:,:,:,iv3d_qc)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qr) = v3dg(:,:,:,iv3d_qr)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qi) = v3dg(:,:,:,iv3d_qi)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qs) = v3dg(:,:,:,iv3d_qs)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qg) = v3dg(:,:,:,iv3d_qg)

  ! RH
  !---------------------------------------------------------

!  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_rh) = [[RH calculator]]

  ! Calculate height based the the topography and vertical coordinate
  !---------------------------------------------------------

  call scale_calc_z(topo, height)
  v3dgh_RP(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_hgt) = height
  v3dgh_RP(KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_hgt) = topo
  v3dgh_RP(1,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_hgt) = 0.0_RP




  ! Surface variables: use the 1st level as the surface (although it is not)
  !---------------------------------------------------------

  v2dgh_RP(:,:,iv2dd_topo) = v3dgh_RP(1+KHALO,:,:,iv3dd_hgt)                ! Use the first model level as topography (is this good?)
!  v2dgh_RP(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_topo) = topo(:,:) ! Use the real topography

  v2dgh_RP(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_ps) = v3dg(1,:,:,iv3d_p)
  v2dgh_RP(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_u10m) = v3dg(1,:,:,iv3d_u)
  v2dgh_RP(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_v10m) = v3dg(1,:,:,iv3d_v)
  v2dgh_RP(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_t2m) = v3dg(1,:,:,iv3d_t)
  v2dgh_RP(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_q2m) = v3dg(1,:,:,iv3d_q)

!  v2dgh_RP(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_rain) = [[No way]]


  ! Communicate the lateral halo areas
  !---------------------------------------------------------

  do iv3d = 1, nv3dd
    call COMM_vars8( v3dgh_RP(:,:,:,iv3d), iv3d )
  end do
  do iv3d = 1, nv3dd
    call COMM_wait ( v3dgh_RP(:,:,:,iv3d), iv3d )
  end do

  do iv2d = 1, nv2dd
    call COMM_vars8( v2dgh_RP(:,:,iv2d), iv2d )
  end do
  do iv2d = 1, nv2dd
    call COMM_wait ( v2dgh_RP(:,:,iv2d), iv2d )
  end do


  ! Pad the upper and lower halo areas
  !---------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(i,j,iv3d) SCHEDULE(STATIC) COLLAPSE(2)
  do iv3d = 1, nv3dd
    do j  = JS, JE
      do i  = IS, IE
        v3dgh_RP(   1:KS-1,i,j,iv3d) = v3dgh_RP(KS,i,j,iv3d)
        v3dgh_RP(KE+1:KA,  i,j,iv3d) = v3dgh_RP(KE,i,j,iv3d)
      end do
    end do
  end do
!$OMP END PARALLEL DO

  v3dgh = real(v3dgh_RP, kind=r_size)
  v2dgh = real(v2dgh_RP, kind=r_size)

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

!-------------------------------------------------------------------------------
! Calculate 3D height coordinate given the topography height (on original grids)
!-------------------------------------------------------------------------------
! [INPUT]
!   nij         : scattered grid numbers
!   topo(nij)   : topography height (on scattered grids)
! [OUTPUT]
!   z(nij,nlev) : 3D height coordinate (on scattered grids)
!-------------------------------------------------------------------------------
subroutine scale_calc_z(topo, z)
  use scale_atmos_grid_cartesC, only: &
     ATMOS_GRID_CARTESC_CZ, &
     ATMOS_GRID_CARTESC_FZ
  use scale_atmos_grid_cartesC_index, only: &
     KHALO, KS, KE
  implicit none

  real(r_size), intent(in) :: topo(nlon,nlat)
  real(r_size), intent(out) :: z(nlev,nlon,nlat)
  real(r_size) :: ztop
  integer :: i, j, k

  ztop = ATMOS_GRID_CARTESC_FZ(KE) - ATMOS_GRID_CARTESC_FZ(KS-1)
!$OMP PARALLEL DO PRIVATE(j,i,k) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev
        z(k, i, j) = (ztop - topo(i,j)) / ztop * ATMOS_GRID_CARTESC_CZ(k+KHALO) + topo(i,j)
      end do
    enddo
  enddo
!$OMP END PARALLEL DO

  return
end subroutine scale_calc_z

!-------------------------------------------------------------------------------
! Calculate 3D height coordinate given the topography height (on scattered grids)
!-------------------------------------------------------------------------------
! [INPUT]
!   nij         : scattered grid numbers
!   topo(nij)   : topography height (on scattered grids)
! [OUTPUT]
!   z(nij,nlev) : 3D height coordinate (on scattered grids)
!-------------------------------------------------------------------------------
subroutine scale_calc_z_grd(nij, topo, z)
  use scale_atmos_grid_cartesC, only: &
     ATMOS_GRID_CARTESC_CZ, &
     ATMOS_GRID_CARTESC_FZ
  use scale_atmos_grid_cartesC_index, only: &
     KHALO, KS, KE
  implicit none

  integer, intent(in) :: nij
  real(r_size), intent(in) :: topo(nij)
  real(r_size), intent(out) :: z(nij,nlev)
  real(r_size) :: ztop
  integer :: k, i

  ztop = ATMOS_GRID_CARTESC_FZ(KE) - ATMOS_GRID_CARTESC_FZ(KS-1)
!$OMP PARALLEL DO PRIVATE(i,k)
  do k = 1, nlev
    do i = 1, nij
      z(i,k) = (ztop - topo(i)) / ztop * ATMOS_GRID_CARTESC_CZ(k+KHALO) + topo(i)
    end do
  end do
!$OMP END PARALLEL DO

  return
end subroutine scale_calc_z_grd

!-------------------------------------------------------------------------------
! Calculate ensemble mean (on scattered grids)
!-------------------------------------------------------------------------------
! [INPUT]
!   mem                     : ensemble size
!   nens                    : ensemble demension of state variables
!   nij                     : scattered grid numbers
!   v3d(nij,nlev,nens,nv3d) : 3D ensemble state variables (on scattered grids)
!                             inputted by (:,:,1..mem,:)
!   v2d(nij,     nens,nv3d) : 2D ensemble state variables (on scattered grids)
!                             inputted by (:,  1..mem,:)
! [OUTPUT]
!   v3d(nij,nlev,nens,nv3d) : ensemble mean of 3D state variables (on scattered grids)
!                             outputted by (:,:,mem+1,:)
!   v2d(nij,     nens,nv3d) : ensemble mean of 2D state variables (on scattered grids)
!                             outputted by (:  ,mem+1,:)
!-------------------------------------------------------------------------------
subroutine ensmean_grd(mem, nens, nij, v3d, v2d)
  implicit none
  integer, intent(in) :: mem
  integer, intent(in) :: nens
  integer, intent(in) :: nij
  real(r_size), intent(inout) :: v3d(nij,nlev,nens,nv3d)
  real(r_size), intent(inout) :: v2d(nij,nens,nv2d)
  integer :: i, k, m, n, mmean

  mmean = mem + 1

!$OMP PARALLEL PRIVATE(i,k,m,n)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  do n = 1, nv3d
    do k = 1, nlev
      do i = 1, nij
        v3d(i,k,mmean,n) = v3d(i,k,1,n)
        do m = 2, mem
          v3d(i,k,mmean,n) = v3d(i,k,mmean,n) + v3d(i,k,m,n)
        end do
        v3d(i,k,mmean,n) = v3d(i,k,mmean,n) / real(mem, r_size)
      end do
    end do
  end do
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  do n = 1, nv2d
    do i = 1, nij
      v2d(i,mmean,n) = v2d(i,1,n)
      do m = 2, mem
        v2d(i,mmean,n) = v2d(i,mmean,n) + v2d(i,m,n)
      end do
      v2d(i,mmean,n) = v2d(i,mmean,n) / real(mem, r_size)
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
end subroutine ensmean_grd

!-------------------------------------------------------------------------------
! Calculate ensemble spread (on scattered grids)
! * the ensemble mean has to be already calculated
!-------------------------------------------------------------------------------
! [INPUT]
!   mem                     : ensemble size
!   nens                    : ensemble demension of state variables
!   nij                     : scattered grid numbers
!   v3d(nij,nlev,nens,nv3d) : 3D ensemble state variables with mean (on scattered grids)
!                             inputted by (:,:,1..mem+1,:)
!   v2d(nij,     nens,nv3d) : 2D ensemble state variables with mean (on scattered grids)
!                             inputted by (:,  1..mem+1,:)
! [OUTPUT]
!   v3ds(nij,nlev,nv3d)     : ensemble spread of 3D state variables (on scattered grids)
!   v2ds(nij,     nv3d)     : ensemble spread of 2D state variables (on scattered grids)
!-------------------------------------------------------------------------------
subroutine enssprd_grd(mem, nens, nij, v3d, v2d, v3ds, v2ds)
  implicit none
  integer, intent(in) :: mem
  integer, intent(in) :: nens
  integer, intent(in) :: nij
  real(r_size), intent(in) :: v3d(nij,nlev,nens,nv3d)
  real(r_size), intent(in) :: v2d(nij,nens,nv2d)
  real(r_size), intent(out) :: v3ds(nij,nlev,nv3d)
  real(r_size), intent(out) :: v2ds(nij,nv2d)
  integer :: i, k, m, n, mmean

  mmean = mem + 1

!$OMP PARALLEL PRIVATE(i,k,m,n)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  do n = 1, nv3d
    do k = 1, nlev
      do i = 1, nij
        v3ds(i,k,n) = (v3d(i,k,1,n) - v3d(i,k,mmean,n)) ** 2
        do m = 2, mem
          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n) - v3d(i,k,mmean,n)) ** 2
        end do
        v3ds(i,k,n) = sqrt(v3ds(i,k,n) / real(mem-1, r_size))
      end do
    end do
  end do
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  do n = 1, nv2d
    do i = 1, nij
      v2ds(i,n) = (v2d(i,1,n) - v2d(i,mmean,n)) ** 2
      do m = 2, mem
        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n) - v2d(i,mmean,n)) ** 2
      end do
      v2ds(i,n) = sqrt(v2ds(i,n) / real(mem-1, r_size))
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
end subroutine enssprd_grd

!-------------------------------------------------------------------------------
! Convert 1D rank of process to 2D rank
!-------------------------------------------------------------------------------
subroutine rank_1d_2d(rank, rank_i, rank_j)
  use scale_prc_cartesC, only: PRC_2Drank
  implicit none
  integer, intent(in) :: rank
  integer, intent(out) :: rank_i, rank_j

  rank_i = PRC_2Drank(rank,1)
  rank_j = PRC_2Drank(rank,2)

  return  
end subroutine rank_1d_2d

!-------------------------------------------------------------------------------
! Convert 2D rank of process to 1D rank
!-------------------------------------------------------------------------------
subroutine rank_2d_1d(rank_i, rank_j, rank)
  use scale_prc_cartesC, only: PRC_NUM_X
  implicit none
  integer, intent(in) :: rank_i, rank_j
  integer, intent(out) :: rank

  rank = rank_j * PRC_NUM_X + rank_i

  return  
end subroutine rank_2d_1d

!-------------------------------------------------------------------------------
! Convert <integer> global grid coordinates (i,j) to local given the 1D rank of process
!-------------------------------------------------------------------------------
subroutine ij_g2l(rank, ig, jg, il, jl)
  implicit none
  integer, intent(in) :: rank
  integer, intent(in) :: ig
  integer, intent(in) :: jg
  integer, intent(out) :: il
  integer, intent(out) :: jl
  integer :: rank_i, rank_j

  call rank_1d_2d(rank, rank_i, rank_j)
  il = ig - rank_i * nlon
  jl = jg - rank_j * nlat

  return  
end subroutine ij_g2l

!-------------------------------------------------------------------------------
! Convert <integer> local grid coordinates (i,j) to global given the 1D rank of process
!-------------------------------------------------------------------------------
subroutine ij_l2g(rank, il, jl, ig, jg)
  implicit none
  integer, intent(in) :: rank
  integer, intent(in) :: il
  integer, intent(in) :: jl
  integer, intent(out) :: ig
  integer, intent(out) :: jg
  integer :: rank_i, rank_j

  call rank_1d_2d(rank, rank_i, rank_j)
  ig = il + rank_i * nlon
  jg = jl + rank_j * nlat

  return  
end subroutine ij_l2g

!-------------------------------------------------------------------------------
! Convert <real> global grid coordinates (i,j) to local given the 1D rank of process
!-------------------------------------------------------------------------------
subroutine rij_g2l(rank, ig, jg, il, jl)
  implicit none
  integer, intent(in) :: rank
  real(r_size), intent(in) :: ig
  real(r_size), intent(in) :: jg
  real(r_size), intent(out) :: il
  real(r_size), intent(out) :: jl
  integer :: rank_i, rank_j

  call rank_1d_2d(rank, rank_i, rank_j)
  il = ig - real(rank_i * nlon, r_size)
  jl = jg - real(rank_j * nlat, r_size)

  return  
end subroutine rij_g2l

!-------------------------------------------------------------------------------
! Convert <real> local grid coordinates (i,j) to global given the 1D rank of process
!-------------------------------------------------------------------------------
subroutine rij_l2g(rank, il, jl, ig, jg)
  implicit none
  integer, intent(in) :: rank
  real(r_size), intent(in) :: il
  real(r_size), intent(in) :: jl
  real(r_size), intent(out) :: ig
  real(r_size), intent(out) :: jg
  integer :: rank_i, rank_j

  call rank_1d_2d(rank, rank_i, rank_j)
  ig = il + real(rank_i * nlon, r_size)
  jg = jl + real(rank_j * nlat, r_size)

  return  
end subroutine rij_l2g

!-------------------------------------------------------------------------------
! Given <real> global grid coordinates (i,j), return the 1D rank of process 
! * HALO grids are used
!-------------------------------------------------------------------------------
! [INPUT]
!   ig, jg : global grid coordinates
! [OUTPUT]
!   rank   : the 1D rank of process where the grid resides;
!            * return -1 if the grid is outside of the global domain
!-------------------------------------------------------------------------------
subroutine rij_rank(ig, jg, rank)
  use scale_prc_cartesC, only: &
      PRC_NUM_X, PRC_NUM_Y
#ifdef DEBUG
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, &
      IA, JA
  use scale_prc, only: &
      PRC_myrank
  use scale_atmos_grid_cartesC, only: &
      GRID_CX => ATMOS_GRID_CARTESC_CX, &
      GRID_CY => ATMOS_GRID_CARTESC_CY, &
      GRID_CXG => ATMOS_GRID_CARTESC_CXG, &
      GRID_CYG => ATMOS_GRID_CARTESC_CYG, &
      DX, &
      DY
#else
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
#endif
  implicit none
  real(r_size), intent(in) :: ig
  real(r_size), intent(in) :: jg
  integer, intent(out) :: rank
  integer :: rank_i, rank_j

  if (ig < real(1+IHALO,r_size) .or. ig > real(nlon*PRC_NUM_X+IHALO,r_size) .or. &
      jg < real(1+JHALO,r_size) .or. jg > real(nlat*PRC_NUM_Y+JHALO,r_size)) then
    rank = -1
    return
  end if

  rank_i = ceiling((ig-real(IHALO,r_size)-0.5d0) / real(nlon,r_size)) - 1
  rank_j = ceiling((jg-real(JHALO,r_size)-0.5d0) / real(nlat,r_size)) - 1
  call rank_2d_1d(rank_i, rank_j, rank)

#ifdef DEBUG
  if (PRC_myrank == rank) then
    if (ig < (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0 .or. &
        ig > (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0 .or. &
        jg < (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0 .or. &
        jg > (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0) then
      write (6,'(A)') '[Error] Process assignment fails!'
      write (6,'(3F10.2)') ig, (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0, (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0
      write (6,'(3F10.2)') jg, (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0, (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0
      stop
    end if
  end if
#endif

  return
end subroutine rij_rank

!-------------------------------------------------------------------------------
! Convert <real> global grid coordinates (i,j) to local where the grid resides
! * HALO grids are used
!-------------------------------------------------------------------------------
! [INPUT]
!   ig, jg : global grid coordinates
! [OUTPUT]
!   rank   : the 1D rank of process where the grid resides;
!            * return -1 if the grid is outside of the global domain
!   il, jl : local grid coordinates
!-------------------------------------------------------------------------------
subroutine rij_rank_g2l(ig, jg, rank, il, jl)
  use scale_prc_cartesC, only: &
      PRC_NUM_X, PRC_NUM_Y
#ifdef DEBUG
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO, &
      IA, JA
  use scale_prc, only: &
      PRC_myrank
  use scale_atmos_grid_cartesC, only: &
      GRID_CX => ATMOS_GRID_CARTESC_CX, &
      GRID_CY => ATMOS_GRID_CARTESC_CY, &
      GRID_CXG => ATMOS_GRID_CARTESC_CXG, &
      GRID_CYG => ATMOS_GRID_CARTESC_CYG, &
      DX, &
      DY
#else
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
#endif
  implicit none
  real(r_size), intent(in) :: ig
  real(r_size), intent(in) :: jg
  integer, intent(out) :: rank
  real(r_size), intent(out) :: il
  real(r_size), intent(out) :: jl
  integer :: rank_i, rank_j

  if (ig < real(1+IHALO,r_size) .or. ig > real(nlon*PRC_NUM_X+IHALO,r_size) .or. &
      jg < real(1+JHALO,r_size) .or. jg > real(nlat*PRC_NUM_Y+JHALO,r_size)) then
    il = -1.0d0
    jl = -1.0d0
    rank = -1
    return
  end if

  rank_i = ceiling((ig-real(IHALO,r_size)-0.5d0) / real(nlon,r_size)) - 1
  rank_j = ceiling((jg-real(JHALO,r_size)-0.5d0) / real(nlat,r_size)) - 1
  il = ig - real(rank_i * nlon, r_size)
  jl = jg - real(rank_j * nlat, r_size)
  call rank_2d_1d(rank_i, rank_j, rank)

#ifdef DEBUG
  if (PRC_myrank == rank) then
    if (ig < (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0 .or. &
        ig > (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0 .or. &
        jg < (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0 .or. &
        jg > (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0) then
      write (6,'(A)') '[Error] Process assignment fails!'
      write (6,'(3F10.2)') ig, (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0, (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0
      write (6,'(3F10.2)') jg, (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0, (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0
      stop
    end if
  end if
#endif

  return
end subroutine rij_rank_g2l

!-------------------------------------------------------------------------------
! Compute the current time labels for various files
!-------------------------------------------------------------------------------
subroutine timelabel_update(lcycle)
  use scale_time, only: &
    TIME_gettimelabel, &
    TIME_time2label, &
    TIME_NOWDATE, &
    TIME_NOWMS
  use scale_calendar, only: &
    CALENDAR_date2daysec, &
    CALENDAR_daysec2date
  implicit none
  real(DP), intent(in) :: lcycle
  integer :: absday, hist_date(6)
  real(DP) :: abssec, hist_subsec

  if (OBS_POSTFIX_TIMELABEL) then
    timelabel_obs = '_???????????????????.dat'
    call TIME_gettimelabel(timelabel_obs(2:20))
  else
    timelabel_obs = ''
  end if
  if (LOG_LEVEL >= 4) then
    write (6, '(2A)') 'Timelabel for observation files: ', timelabel_obs
  end if

  if (GUES_ANAL_POSTFIX_TIMELABEL) then
    timelabel_anal = '_???????????????????'
    call TIME_gettimelabel(timelabel_anal(2:20))
  else
    timelabel_anal = ''
  end if
  if (LOG_LEVEL >= 4) then
    write (6, '(2A)') 'Timelabel for guess/analysis files: ', timelabel_anal
  end if

  if (HISTORY_POSTFIX_TIMELABEL) then
    timelabel_hist = '_???????????????????'
    call CALENDAR_date2daysec( absday,       & ! [OUT]
                               abssec,       & ! [OUT]
                               TIME_NOWDATE, & ! [IN]
                               TIME_NOWMS,   & ! [IN]
                               0             ) ! [IN]
    abssec = abssec - lcycle
    call CALENDAR_daysec2date( hist_date,   & ! [OUT]
                               hist_subsec, & ! [OUT]
                               absday,      & ! [IN]
                               abssec,      & ! [IN]
                               0            ) ! [IN]
    call TIME_time2label( hist_date,           & ! [IN]
                          hist_subsec,         & ! [IN]
                          timelabel_hist(2:20) ) ! [OUT]
  else
    timelabel_hist = ''
  end if
  if (LOG_LEVEL >= 4) then
    write (6, '(2A)') 'Timelabel for history files: ', timelabel_hist
  end if

end subroutine timelabel_update

subroutine set_lonlat2d()
  use scale_atmos_grid_cartesC, only: &
      CX => ATMOS_GRID_CARTESC_CX, &
      CY => ATMOS_GRID_CARTESC_CY, &
      DX, DY
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  implicit none

  integer :: i, j
  real(r_size) :: ri, rj
  real(RP) :: lon2d_RP(nlon,nlat)
  real(RP) :: lat2d_RP(nlon,nlat)

  if (.not.(allocated(lon2d))) allocate(lon2d(nlong,nlatg))
  if (.not.(allocated(lat2d))) allocate(lat2d(nlong,nlatg))

!$OMP PARALLEL DO PRIVATE(i,j,ri,rj) COLLAPSE(2)
      do j = 1, nlat
        do i = 1, nlon
          ri = real(i + IHALO, r_size)
          rj = real(j + JHALO, r_size)
          call MAPPROJECTION_xy2lonlat(real((ri-1.0_r_size) * DX + CX(1),RP), &
                                       real((rj-1.0_r_size) * DY + CY(1),RP), &
                                       lon2d_RP(i,j), lat2d_RP(i,j))
          lon2d(i,j) = lon2d_RP(i,j) * rad2deg
          lat2d(i,j) = lat2d_RP(i,j) * rad2deg
        end do
      end do
!$OMP END PARALLEL DO


  return
end subroutine set_lonlat2d

subroutine jst2utc(jyear, jmonth, jday, jhour, jminute, jsecond, jtime_ms, utime)
  use scale_calendar, only: &
      CALENDAR_date2daysec, &
      CALENDAR_daysec2date, &
      CALENDAR_adjust_daysec
  implicit none

  integer, intent(in) :: jyear, jmonth, jday
  integer, intent(in) :: jhour, jminute, jsecond
  real(DP) :: jtime_ms
  integer, intent(out) :: utime(6)
  integer :: jtime(6)
  integer :: absday
  real(DP) :: abssec, utime_ms

  jtime(1) = jyear
  jtime(2) = jmonth
  jtime(3) = jday
  jtime(4) = jhour
  jtime(5) = jminute
  jtime(6) = jsecond

  call CALENDAR_date2daysec( absday,       & ! [OUT]
                             abssec,       & ! [OUT]
                             jtime,        & ! [IN]
                             jtime_ms,     & ! [IN]
                             0             ) ! [IN]

  abssec = abssec - real(3600*9, kind=DP)

  call CALENDAR_adjust_daysec( absday,   & ! [INOUT]
                               abssec )    ! [INOUT]

  call CALENDAR_daysec2date( utime,   & ! [OUT]
                             utime_ms, & ! [OUT]
                             absday,      & ! [IN]
                             abssec,      & ! [IN]
                             0            ) ! [IN]

  return
end subroutine jst2utc

subroutine advance_nowdate( date, dsec )
  use scale_calendar, only: &
      CALENDAR_date2daysec, &
      CALENDAR_daysec2date, &
      CALENDAR_adjust_daysec
  implicit none

  integer, intent(inout) :: date(6)
  real(DP) :: dsec
  integer :: absday
  real(DP) :: abssec
  real(DP) :: date_ms = 0.0_DP

  call CALENDAR_date2daysec( absday,       & ! [OUT]
                             abssec,       & ! [OUT]
                             date,         & ! [IN]
                             date_ms,       & ! [IN]
                             0             ) ! [IN]

  abssec = abssec + dsec

  call CALENDAR_adjust_daysec( absday,   & ! [INOUT]
                               abssec )    ! [INOUT]

  call CALENDAR_daysec2date( date,        & ! [OUT]
                             date_ms,     & ! [OUT]
                             absday,      & ! [IN]
                             abssec,      & ! [IN]
                             0            ) ! [IN]


  return
end subroutine advance_nowdate

!===============================================================================
END MODULE common_scale
