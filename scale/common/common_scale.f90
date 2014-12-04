MODULE common_scale
!=======================================================================
!
! [PURPOSE:] Common Information for SCALE
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   10/03/2012 Guo-Yuan Lien     modified for GFS model
!   07/24/2014 Guo-Yuan Lien     modified for SCALE model
!
!=======================================================================
!$USE OMP_LIB
  USE common
!  USE common_ncio

  use scale_stdio
!  use scale_stdio, only: H_MID

  use scale_precision, only: RP

  use scale_prof

!  use common_mpi

!  use common_mpi_scale, only: &
!    nitmax, &
!    proc2mem

!  use common_letkf, only: nbv

!  use scale_precision
!  use scale_stdio
!  use scale_prof
!  use scale_grid_index


!  use dc_log, only: &
!    loginit
!  use gtool_file, only: &
!!     fileread, &
!     filecloseall
!  use scale_grid_index, only: &
!    KHALO, IHALO, JHALO


  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: nlonsub=30
  INTEGER,PARAMETER :: nlatsub=30
  INTEGER,PARAMETER :: nlonns=2
  INTEGER,PARAMETER :: nlatns=2
  INTEGER,PARAMETER :: nlon=nlonsub*nlonns
  INTEGER,PARAMETER :: nlat=nlatsub*nlatns
  INTEGER,PARAMETER :: nlev=30

  integer,parameter :: nlonhalo=nlonsub+4
  integer,parameter :: nlathalo=nlatsub+4
  integer,parameter :: nlevhalo=nlev+4

  INTEGER,PARAMETER :: nv3d=11   ! 3D state variables (in SCALE init/restart files)
  INTEGER,PARAMETER :: nv3dd=13  ! 3D diagnostic variables (in SCALE history files)
  INTEGER,PARAMETER :: nv3dx=13  ! 3D diagnostic variables
  INTEGER,PARAMETER :: nv2d=0    ! 2D state variables (in SCALE init/restart files)
  INTEGER,PARAMETER :: nv2dd=7   ! 2D diagnostic variables (in SCALE history files)
  INTEGER,PARAMETER :: nv2dx=7   ! 2D diagnostic variables
!  INTEGER,PARAMETER :: iv3d_rho=1
!  INTEGER,PARAMETER :: iv3d_rhou=2
!  INTEGER,PARAMETER :: iv3d_rhov=3
!  INTEGER,PARAMETER :: iv3d_rhow=4
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_w=3
  INTEGER,PARAMETER :: iv3d_t=4
  INTEGER,PARAMETER :: iv3d_rhot=5 ! (~pressure)
  INTEGER,PARAMETER :: iv3d_q=6
  INTEGER,PARAMETER :: iv3d_qc=7
  INTEGER,PARAMETER :: iv3d_qr=8
  INTEGER,PARAMETER :: iv3d_qi=9
  INTEGER,PARAMETER :: iv3d_qs=10
  INTEGER,PARAMETER :: iv3d_qg=11
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
!  INTEGER,PARAMETER :: iv2dd_hgt=1
  INTEGER,PARAMETER :: iv2dd_ps=1
  INTEGER,PARAMETER :: iv2dd_tsfc=2
  INTEGER,PARAMETER :: iv2dd_rain=3
  INTEGER,PARAMETER :: iv2dd_u10m=4
  INTEGER,PARAMETER :: iv2dd_v10m=5
  INTEGER,PARAMETER :: iv2dd_t2m=6
  INTEGER,PARAMETER :: iv2dd_q2m=7
  INTEGER,PARAMETER :: nij0=nlonsub*nlatsub
  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,PARAMETER :: nlevalld=nlev*nv3dd+nv2dd
  INTEGER,PARAMETER :: nlevallx=nlev*nv3dx+nv2dx
  INTEGER,PARAMETER :: ngpv=nij0*nlevall
  INTEGER,PARAMETER :: ngpvd=nij0*nlevalld
  INTEGER,PARAMETER :: ngpvx=nij0*nlevallx
  REAL(r_size),SAVE :: lon(nlonsub,nlatsub)
  REAL(r_size),SAVE :: lat(nlonsub,nlatsub)
  REAL(r_size),SAVE :: lonu(nlonsub,nlatsub)
  REAL(r_size),SAVE :: latu(nlonsub,nlatsub)
  REAL(r_size),SAVE :: lonv(nlonsub,nlatsub)
  REAL(r_size),SAVE :: latv(nlonsub,nlatsub)
!  REAL(r_size),SAVE :: dx(nlatsub)
!  REAL(r_size),SAVE :: dy(nlatsub)
!  REAL(r_size),SAVE :: dy2(nlatsub)
  REAL(r_size),SAVE :: fcori(nlatsub)
!  REAL(r_size),SAVE :: wg(nlonsub,nlatsub)
  INTEGER,PARAMETER :: vname_max = 10
  CHARACTER(vname_max),SAVE :: v3d_name(nv3d)
  CHARACTER(vname_max),SAVE :: v3dd_name(nv3dx)
  CHARACTER(vname_max),SAVE :: v2d_name(nv2d)
  CHARACTER(vname_max),SAVE :: v2dd_name(nv2dx)

  character(len=H_MID), parameter :: MODELNAME = "SCALE-LES"

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_scale
  IMPLICIT NONE
!  REAL(r_sngl) :: slat(nlat), wlat(nlat)
!  REAL(r_size) :: totalwg, wgtmp, latm1, latm2
  INTEGER :: i,j

  WRITE(6,'(A)') 'Hello from set_common_scale'
  !
  ! Variable names (same as in the NetCDF file)
  !
  ! state variables (in 'restart' files, for LETKF)
!  v3d_name(iv3d_rho)  = 'DENS'
!  v3d_name(iv3d_rhou) = 'MOMX'
!  v3d_name(iv3d_rhov) = 'MOMY'
!  v3d_name(iv3d_rhow) = 'MOMZ'
  v3d_name(iv3d_rhot) = 'RHOT'
  v3d_name(iv3d_q)    = 'QV'
  v3d_name(iv3d_qc)   = 'QC'
  v3d_name(iv3d_qr)   = 'QR'
  v3d_name(iv3d_qi)   = 'QI'
  v3d_name(iv3d_qs)   = 'QS'
  v3d_name(iv3d_qg)   = 'QG'
  !
  ! diagnostic variables (in 'history' files, for observation operators)
  v3dd_name(iv3dd_u)    = 'U'
  v3dd_name(iv3dd_v)    = 'V'
  v3dd_name(iv3dd_w)    = 'W'
  v3dd_name(iv3dd_t)    = 'T'
  v3dd_name(iv3dd_p)    = 'PRES'
  v3dd_name(iv3dd_q)    = 'QV'
  v3dd_name(iv3dd_qc)   = 'QC'
  v3dd_name(iv3dd_qr)   = 'QR'
  v3dd_name(iv3dd_qi)   = 'QI'
  v3dd_name(iv3dd_qs)   = 'QS'
  v3dd_name(iv3dd_qg)   = 'QG'
  v3dd_name(iv3dd_rh)   = 'RH'
  v3dd_name(iv3dd_hgt)   = 'height'
  !
!  v2dd_name(iv2dd_hgt) = 'height'
  v2dd_name(iv2dd_ps) = 'SFC_PRES'
  v2dd_name(iv2dd_tsfc) = 'SFC_TEMP'
  v2dd_name(iv2dd_rain) = 'PREC'
  v2dd_name(iv2dd_u10m) = 'U10'
  v2dd_name(iv2dd_v10m) = 'V10'
  v2dd_name(iv2dd_t2m) = 'T2'
  v2dd_name(iv2dd_q2m) = 'Q2'
  !
  ! Lon, Lat
  !
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

!-----------------------------------------------------------------------
! Start using SCALE library
!-----------------------------------------------------------------------
subroutine set_scalelib(mem_np, nitmax, nprocs, proc2mem)
!subroutine set_scalelib(mem_np)
!subroutine set_scalelib(nitmax, proc2mem)
!subroutine set_scalelib
!  use common_nml, only: &
!    MEM_NP

  use scale_precision
!  use scale_stdio
!  use scale_prof

  use gtool_history, only: &
    historyinit
  use dc_log, only: &
    LogInit
  use scale_process, only: &
    PRC_setup,    &
    PRC_MPIstart, &
!      PRC_mpifinish, &
    PRC_master, &
    PRC_myrank, &
    PRC_myrank_world, &
    PRC_2Drank, &
    PRC_NUM_X, &
    PRC_NUM_Y
!    prc_nu, &
  use scale_const, only: &
    CONST_setup
  use scale_calendar, only: &
    CALENDAR_setup
  use scale_random, only: &
    RANDOM_setup
  use scale_time, only: &
    TIME_setup
  use scale_grid, only: &
    GRID_setup, &
    GRID_DOMAIN_CENTER_X, &
    GRID_DOMAIN_CENTER_Y

  use scale_grid_index

!  use scale_grid_nest, only: &
!    NEST_setup
!  use scale_land_grid_index, only: &
!    LAND_GRID_INDEX_setup
!  use scale_land_grid, only: &
!    LAND_GRID_setup
!  use scale_urban_grid_index, only: &
!    URBAN_GRID_INDEX_setup
!  use scale_urban_grid, only: &
!    URBAN_GRID_setup
!  use scale_tracer, only: &
!    TRACER_setup
  use scale_fileio, only: &
     FILEIO_setup
  use scale_comm, only: &
    COMM_setup
!  use scale_topography, only: &
!    TOPO_setup
!  use scale_landuse, only: &
!    LANDUSE_setup
!  use scale_grid_real, only: &
!    REAL_setup
  use scale_mapproj, only: &
    MPRJ_setup
  implicit none

!  integer,intent(in) :: nitmax ! maximum number of model files processed by a process
!  integer,intent(in) :: proc2mem(2,nitmax,nprocs)
!  integer,intent(in) :: proc2mem(:,:,:)

!  integer,intent(in) :: mem_np
  integer, intent(in) :: mem_np, nitmax, nprocs
  integer, intent(in) :: proc2mem(2,nitmax,nprocs)

!    character(len=H_MID) :: DATATYPE = 'DEFAULT' !< REAL4 or REAL8
  integer :: rankidx(2)

  !-----------------------------------------------------------------------------

  ! start SCALE MPI
  call PRC_MPIstart(mem_np, nitmax, nprocs, proc2mem)
!  call PRC_MPIstart(nbv, mem_np, nitmax, nprocs, proc2mem)

  ! setup process
  call PRC_setup

  ! setup Log
  call LogInit(IO_FID_CONF, IO_FID_LOG, IO_L)

  ! setup constants
  call CONST_setup

  ! setup time
  call TIME_setup( setup_TimeIntegration = .false. )

  call PROF_rapstart('Initialize')

  ! setup horizontal/vertical grid coordinates
  call GRID_INDEX_setup
  call GRID_setup

  ! check if the namelist seetings are consistent
  if (mem_np /= PRC_NUM_X * PRC_NUM_Y) then
    write(6,*) 'mem_np should be equal to PRC_NUM_X * PRC_NUM_Y.'
    stop
  else if (IMAX /= nlonsub) then
    write(6,*) 'IMAX should be equal to nlonsub.'
    stop
  else if (JMAX /= nlatsub) then
    write(6,*) 'JMAX should be equal to nlatsub.'
    stop
  else if (KMAX /= nlev) then
    write(6,*) 'KMAX should be equal to nlev.'
    stop
  end if

!  call LAND_GRID_INDEX_setup
!  call LAND_GRID_setup

!  call URBAN_GRID_INDEX_setup
!  call URBAN_GRID_setup

  ! setup file I/O
  call FILEIO_setup

  ! setup mpi communication
  call COMM_setup

  ! setup topography
!  call TOPO_setup
  ! setup land use category index/fraction
!  call LANDUSE_setup

  ! setup grid coordinates (real world)
!  call REAL_setup

  ! setup map projection
  call MPRJ_setup( GRID_DOMAIN_CENTER_X, GRID_DOMAIN_CENTER_Y )

  ! setup history file I/O
  rankidx(1) = PRC_2Drank(PRC_myrank, 1)
  rankidx(2) = PRC_2Drank(PRC_myrank, 2)
  call HistoryInit('','','',IMAX*JMAX*KMAX,PRC_master,PRC_myrank,rankidx)

  call PROF_rapend('Initialize')

  call PROF_rapstart('Main')

  return
end subroutine set_scalelib

!-----------------------------------------------------------------------
! Finish using SCALE library
!-----------------------------------------------------------------------
subroutine unset_scalelib
  use scale_process, only: &
    PRC_MPIfinish
  use gtool_file, only: &
    FileCloseAll
  implicit none

  call PROF_rapend('Main')

  call PROF_rapreport

  call FileCloseAll

  ! stop SCALE MPI
  call PRC_MPIfinish

  return
end subroutine unset_scalelib

!!-----------------------------------------------------------------------
!! File I/O
!!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE read_restart(filename,v3dg,v2dg)
  use gtool_file, only: FileRead
  use scale_process, only: PRC_myrank
  use common_mpi, only: myrank
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlonsub,nlatsub,nv2d)
  INTEGER :: iv3d,iv2d

  WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',PRC_myrank,'.nc'

  if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', 'DENS'
  call FileRead(v3dg(:,:,:,6), filename, 'DENS', 1, PRC_myrank)         ! temporarily save rho

  if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', 'MOMX'
  call FileRead(v3dg(:,:,:,iv3d_u), filename, 'MOMX', 1, PRC_myrank)
  v3dg(:,:,:,iv3d_u) = v3dg(:,:,:,iv3d_u) / v3dg(:,:,:,6)               ! -- get U

  if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', 'MOMY'
  call FileRead(v3dg(:,:,:,iv3d_v), filename, 'MOMY', 1, PRC_myrank)
  v3dg(:,:,:,iv3d_v) = v3dg(:,:,:,iv3d_v) / v3dg(:,:,:,6)               ! -- get V

  if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', 'MOMZ'
  call FileRead(v3dg(:,:,:,iv3d_w), filename, 'MOMZ', 1, PRC_myrank)
  v3dg(:,:,:,iv3d_w) = v3dg(:,:,:,iv3d_w) / v3dg(:,:,:,6)               ! -- get W

  if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', 'RHOT'
  call FileRead(v3dg(:,:,:,iv3d_rhot), filename, 'RHOT', 1, PRC_myrank) ! -- get rho * T
  v3dg(:,:,:,iv3d_t) = v3dg(:,:,:,iv3d_rhot) / v3dg(:,:,:,6)            ! -- get T

  do iv3d = 6, nv3d
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
    call FileRead(v3dg(:,:,:,iv3d), filename, trim(v3d_name(iv3d)), 1, PRC_myrank)
  end do

  do iv2d = 1, nv2d
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
    call FileRead(v2dg(:,:,iv2d), filename, trim(v2d_name(iv2d)), 1, PRC_myrank)
  end do

  RETURN
END SUBROUTINE read_restart
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE write_restart(filename,v3dg,v2dg)
  use netcdf, only: NF90_WRITE
  use scale_process, only: PRC_myrank
  use common_mpi, only: myrank
  use common_ncio
  implicit none

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlonsub,nlatsub,nv2d)
  REAL(RP) :: rhotmp(nlev,nlonsub,nlatsub)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: iv3d,iv2d,ncid

  WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',PRC_myrank,'.nc'
  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  call ncio_open(filename // filesuffix, NF90_WRITE, ncid)

  rhotmp = v3dg(:,:,:,iv3d_rhot) / v3dg(:,:,:,iv3d_t)

  if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Write 3D var: ', 'DENS'
  call ncio_write_3d_r8(ncid, 'DENS', nlev, nlonsub, nlatsub, 1, rhotmp) ! -- write rho

  if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Write 3D var: ', 'MOMX'
  call ncio_write_3d_r8(ncid, 'MOMX', nlev, nlonsub, nlatsub, 1, rhotmp * v3dg(:,:,:,iv3d_u)) ! -- write rho*u

  if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Write 3D var: ', 'MOMY'
  call ncio_write_3d_r8(ncid, 'MOMY', nlev, nlonsub, nlatsub, 1, rhotmp * v3dg(:,:,:,iv3d_v)) ! -- write rho*v

  if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Write 3D var: ', 'MOMZ'
  call ncio_write_3d_r8(ncid, 'MOMZ', nlev, nlonsub, nlatsub, 1, rhotmp * v3dg(:,:,:,iv3d_w)) ! -- write rho*w

  if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Write 3D var: ', 'RHOT'
  call ncio_write_3d_r8(ncid, 'RHOT', nlev, nlonsub, nlatsub, 1, v3dg(:,:,:,iv3d_rhot)) ! -- write rho*T

  DO iv3d = 6, nv3d
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
    call ncio_write_3d_r8(ncid, trim(v3d_name(iv3d)), nlev, nlonsub, nlatsub, 1, v3dg(:,:,:,iv3d))
  END DO

  DO iv2d = 1, nv2d
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
    call ncio_write_2d_r8(ncid, trim(v2d_name(iv2d)), nlonsub, nlatsub, 1, v2dg(:,:,iv2d))
  END DO

  call ncio_close(ncid)

  RETURN
END SUBROUTINE write_restart


!subroutine FileWrite_combine_3D(var,basename,varname,step,myrank)

!  use gtool_file, only: &
!    FileWrite, &
!    File
!  use scale_process, only: &
!    PRC_master, &
!    PRC_myrank, &
!    PRC_2Drank
!  implicit none

!  real(DP),         intent(in)           :: var(:,:,:)
!  character(LEN=*), intent(in)           :: basename
!  character(LEN=*), intent(in)           :: varname
!  integer,          intent(in)           :: step
!  integer,          intent(in)           :: myrank

!  real(RP),          intent(in)  :: var(:,:,:) !< value of the variable
!  character(len=*),  intent(in)  :: basename   !< basename of the file
!  character(len=*),  intent(in)  :: title      !< title    of the file
!  character(len=*),  intent(in)  :: varname    !< name        of the variable
!  character(len=*),  intent(in)  :: desc       !< description of the variable
!  character(len=*),  intent(in)  :: unit       !< unit        of the variable
!  character(len=*),  intent(in)  :: axistype   !< axis type (Z/X/Y)
!  character(len=*),  intent(in)  :: datatype   !< data type (REAL8/REAL4/default)
!  logical, optional, intent(in)  :: append     !< append existing (closed) file?

!!  integer          :: dtype
!!  character(len=2) :: dims(3)
!!  integer          :: dim1_max, dim1_S, dim1_E
!!  integer          :: dim2_max, dim2_S, dim2_E
!!  integer          :: dim3_max, dim3_S, dim3_E

!  real(RP), allocatable :: var3D(:,:,:)

!  integer :: rankidx(2)
!  logical :: append_sw
!  logical :: fileexisted
!  integer :: fid, vid
!  !---------------------------------------------------------------------------

!  rankidx(1) = PRC_2Drank(PRC_myrank,1)
!  rankidx(2) = PRC_2Drank(PRC_myrank,2)
!  if ( RP == 8 ) then
!    dtype = File_REAL8
!  elseif( RP == 4 ) then
!    dtype = File_REAL4
!  endif
!  call FileCreate( fid,               & ! [OUT]
!                   fileexisted,       & ! [OUT]
!                   basename,          & ! [IN]
!                   '', & !title,             & ! [IN]
!                   '', & !H_SOURCE,          & ! [IN]
!                   '', & !H_INSTITUTE,       & ! [IN]
!                   PRC_master,        & ! [IN]
!                   PRC_myrank,        & ! [IN]
!                   rankidx,           & ! [IN]
!                   .true. ) ! [IN]

!    if ( .NOT. fileexisted ) then ! only once
!       call FILEIO_set_axes( fid, dtype ) ! [IN]
!    endif

!    call FileAddVariable( vid, fid, varname, desc, '', dims, dtype ) ! [IN]

!    allocate( var3D(dim1_max,dim2_max,dim3_max) )

!    var3D(1:dim1_max,1:dim2_max,1:dim3_max) = var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E)
!    call FileWrite( vid, var3D(:,:,:), NOWSEC, NOWSEC ) ! [IN]

!    deallocate( var3D )
!  end subroutine FILEIO_write_3D



!!-- Read a grid file ---------------------------------------------------
!SUBROUTINE read_grd(filename,inv3d,inv2d,iv3dname,iv2dname,v3d,v2d)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: filename
!  INTEGER,INTENT(IN) :: inv3d,inv2d
!  CHARACTER(vname_max),INTENT(IN) :: iv3dname(inv3d)
!  CHARACTER(vname_max),INTENT(IN) :: iv2dname(inv2d)
!  REAL(r_size),INTENT(OUT) :: v3d(nlonsub,nlatsub,nlev,inv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nlonsub,nlatsub,inv2d)
!  REAL(r_dble) :: buf3d8(nlonsub,nlatsub,nlev)
!  REAL(r_dble) :: buf2d8(nlonsub,nlatsub)
!  INTEGER :: iunit,n

!  CALL ncio_open(filename, nf90_nowrite, iunit)
!  DO n=1,inv3d
!    IF(r_size == r_dble) THEN
!      CALL ncio_read_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,v3d(:,:,:,n))
!    ELSE
!      CALL ncio_read_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,buf3d8)
!      v3d(:,:,:,n) = REAL(buf3d8,r_size)
!    END IF
!  END DO
!  DO n=1,inv2d
!    IF(r_size == r_dble) THEN
!      CALL ncio_read_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,v2d(:,:,n))
!    ELSE
!      CALL ncio_read_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,buf2d8)
!      v2d(:,:,n) = REAL(buf2d8,r_size)
!    END IF
!  END DO
!  CALL ncio_close(iunit)

!  RETURN
!END SUBROUTINE read_grd

!SUBROUTINE read_grd4(filename,inv3d,inv2d,iv3dname,iv2dname,v3d,v2d)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: filename
!  INTEGER,INTENT(IN) :: inv3d,inv2d
!  CHARACTER(vname_max),INTENT(IN) :: iv3dname(inv3d)
!  CHARACTER(vname_max),INTENT(IN) :: iv2dname(inv2d)
!  REAL(r_sngl),INTENT(OUT) :: v3d(nlonsub,nlatsub,nlev,inv3d)
!  REAL(r_sngl),INTENT(OUT) :: v2d(nlonsub,nlatsub,inv2d)
!  REAL(r_dble) :: buf3d8(nlonsub,nlatsub,nlev)
!  REAL(r_dble) :: buf2d8(nlonsub,nlatsub)
!  INTEGER :: iunit,n

!  CALL ncio_open(filename, nf90_nowrite, iunit)
!  DO n=1,inv3d
!    CALL ncio_read_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,buf3d8)
!    v3d(:,:,:,n) = REAL(buf3d8,r_sngl)
!  END DO
!  DO n=1,inv2d
!    CALL ncio_read_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,buf2d8)
!    v2d(:,:,n) = REAL(buf2d8,r_sngl)
!  END DO
!  CALL ncio_close(iunit)

!  RETURN
!END SUBROUTINE read_grd4

!!-- Write a grid file -------------------------------------------------
!SUBROUTINE write_grd(filename,inv3d,inv2d,iv3dname,iv2dname,v3d,v2d)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: filename
!  INTEGER,INTENT(IN) :: inv3d,inv2d
!  CHARACTER(vname_max),INTENT(IN) :: iv3dname(inv3d)
!  CHARACTER(vname_max),INTENT(IN) :: iv2dname(inv2d)
!  REAL(r_size),INTENT(IN) :: v3d(nlonsub,nlatsub,nlev,inv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nlonsub,nlatsub,inv2d)
!  REAL(r_dble) :: buf3d8(nlonsub,nlatsub,nlev)
!  REAL(r_dble) :: buf2d8(nlonsub,nlatsub)
!  INTEGER :: iunit,n

!  CALL ncio_open(filename, nf90_write, iunit)
!  DO n=1,inv3d
!    IF(r_size == r_dble) THEN
!      CALL ncio_write_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,v3d(:,:,:,n))
!    ELSE
!      buf3d8 = REAL(v3d(:,:,:,n),r_dble)
!      CALL ncio_write_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,buf3d8)
!    END IF
!  END DO
!  DO n=1,inv2d
!    IF(r_size == r_dble) THEN
!      CALL ncio_write_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,v2d(:,:,n))
!    ELSE
!      buf2d8 = REAL(v2d(:,:,n),r_dble)
!      CALL ncio_write_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,buf2d8)
!    END IF
!  END DO
!  CALL ncio_close(iunit)

!  RETURN
!END SUBROUTINE write_grd

!SUBROUTINE write_grd4(filename,inv3d,inv2d,iv3dname,iv2dname,v3d,v2d)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: filename
!  INTEGER,INTENT(IN) :: inv3d,inv2d
!  CHARACTER(vname_max),INTENT(IN) :: iv3dname(inv3d)
!  CHARACTER(vname_max),INTENT(IN) :: iv2dname(inv2d)
!  REAL(r_sngl),INTENT(IN) :: v3d(nlonsub,nlatsub,nlev,inv3d)
!  REAL(r_sngl),INTENT(IN) :: v2d(nlonsub,nlatsub,inv2d)
!  REAL(r_dble) :: buf3d8(nlonsub,nlatsub,nlev)
!  REAL(r_dble) :: buf2d8(nlonsub,nlatsub)
!  INTEGER :: iunit,n

!  CALL ncio_open(filename, nf90_write, iunit)
!  DO n=1,inv3d
!    buf3d8 = REAL(v3d(:,:,:,n),r_dble)
!    CALL ncio_write_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,buf3d8)
!  END DO
!  DO n=1,inv2d
!    buf2d8 = REAL(v2d(:,:,n),r_dble)
!    CALL ncio_write_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,buf2d8)
!  END DO
!  CALL ncio_close(iunit)

!  RETURN
!END SUBROUTINE write_grd4
!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: i,k,m,n

  DO n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        DO m=2,member
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        END DO
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO

  DO n=1,nv2d
    DO i=1,nij
      v2dm(i,n) = v2d(i,1,n)
      DO m=2,member
        v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
      END DO
      v2dm(i,n) = v2dm(i,n) / REAL(member,r_size)
    END DO
  END DO

  RETURN
END SUBROUTINE ensmean_grd

END MODULE common_scale
