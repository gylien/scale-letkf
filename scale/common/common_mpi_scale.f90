module common_mpi_scale
!===============================================================================
!
! [PURPOSE:] MPI procedures
!
! [ATTENTION:]
!   DO NOT COMPILE WITH BOTH INLINE EXPANSION AND OMP OPTIONS TOGETHER
!    (Use ONE if you want, but DON'T USE BOTH AT THE SAME TIME)
!
! [HISTORY:]
!   01/23/2009 Takemasa Miyoshi  created
!   10/03/2012 Guo-Yuan Lien     modified for GFS model
!   12/30/2013 Guo-Yuan Lien     add get_nobs_mpi and read_obs2_mpi
!   08/14/2014 Guo-Yuan Lien     modified for SCALE model
!   .......... See git history for the following revisions
!
!===============================================================================
!$USE OMP_LIB
  use common
  use common_nml
  use common_mpi
  use common_scale

  use scale_precision, only: RP
  use scale_comm_cartesC, only: COMM_datatype
  use scale_file, only: FILE_AGGREGATE

  implicit none
  public

  integer,save :: nnodes
  integer,save :: nprocs_m

  integer,save :: nij1
  integer,save :: nij1max
  integer,allocatable,save :: nij1node(:)
  real(r_size),allocatable,save :: rig1(:),rjg1(:)
  real(r_size),allocatable,save :: topo1(:)
  real(r_size),allocatable,save :: hgt1(:,:)

  integer,save :: n_mem
  integer,save :: n_mempn
  integer,save :: nitmax ! maximum number of model files processed by a process

  integer,allocatable,save :: mempe_to_node(:,:)   ! No use in the LETKF code
  integer,allocatable,save :: mempe_to_rank(:,:)
  integer,allocatable,save :: rank_to_mem(:,:)
  integer,allocatable,save :: rank_to_pe(:)
  integer,allocatable,save :: rank_to_mempe(:,:,:) ! Deprecated except for the use in PRC_MPIsplit_letkf
  integer,allocatable,save :: ranke_to_mem(:,:)
  integer,allocatable,save :: myrank_to_mem(:)
  integer,save :: myrank_to_pe
  logical,save :: myrank_use = .false.
  logical,save :: myrank_use_da = .false.
  logical,save :: myrank_use_obs = .false.

  integer,save :: mydom = -1

  integer,save :: nens = -1
  integer,save :: nensobs = -1

  integer,save :: mmean = -99    ! use a value different from -1 to avoid equivalence to (my)rank_to_mem
  integer,save :: mmdet = -99    ! when no member is corresponded to a rank/iteration
  integer,save :: mmdetin = -99  ! 
  integer,save :: mmdetobs = -99 ! 

  integer,save :: mmean_rank_e = -1
  integer,save :: mmdet_rank_e = -1
  integer,save :: msprd_rank_e = -1

  integer,save :: MPI_COMM_u, nprocs_u, myrank_u
  integer,save :: MPI_COMM_a, nprocs_a, myrank_a
  integer,save :: MPI_COMM_da, nprocs_da, myrank_da
  integer,save :: MPI_COMM_d, nprocs_d, myrank_d
  integer,save :: MPI_COMM_e, nprocs_e, myrank_e
  integer,save :: MPI_COMM_ef, nprocs_ef, myrank_ef
  integer,save :: MPI_COMM_o, nprocs_o, myrank_o

  integer, parameter :: max_timer_levels = 5
  integer, parameter :: timer_name_width = 50
  real(r_dble), private, parameter :: timer_neglect = 1.0d-3
  real(r_dble), private, save :: timer_save(max_timer_levels) = -9.0d10

  integer, save :: MPI_RP

contains

!-------------------------------------------------------------------------------
! initialize_mpi_scale
!-------------------------------------------------------------------------------
subroutine initialize_mpi_scale
  use scale_prc, only: &
     PRC_MPIstart, &
     PRC_UNIVERSAL_setup, &
     PRC_UNIVERSAL_myrank
  implicit none
  integer :: universal_comm   ! dummy
  integer :: universal_nprocs ! dummy
  logical :: universal_master ! dummy
!  integer :: ierr

  call PRC_MPIstart( universal_comm ) ! [OUT]

!  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
!  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call PRC_UNIVERSAL_setup( universal_comm,   & ! [IN]
                            universal_nprocs, & ! [OUT]
                            universal_master  ) ! [OUT]
  nprocs = universal_nprocs
  myrank = PRC_UNIVERSAL_myrank

  if (myrank == 0) write(6,'(A,I6.6,A,I6.6)') 'Hello from MYRANK ', myrank, '/', nprocs-1
  if (r_size == r_dble) then
    MPI_r_size = MPI_DOUBLE_PRECISION
  else if (r_size == r_sngl) then
    MPI_r_size = MPI_REAL
  end if

  if (RP == r_dble) then
    MPI_RP = MPI_DOUBLE_PRECISION
  else if (RP == r_sngl) then
    MPI_RP = MPI_REAL
  end if

  return
end subroutine initialize_mpi_scale

!-------------------------------------------------------------------------------
! finalize_mpi_scale
!-------------------------------------------------------------------------------
subroutine finalize_mpi_scale
!  use scale_prc, only: PRC_MPIfinish
  implicit none
  integer :: ierr

!  call PRC_MPIfinish
  call MPI_Finalize(ierr)

  return
end subroutine finalize_mpi_scale

!-------------------------------------------------------------------------------
! set_common_mpi_scale
!-------------------------------------------------------------------------------
subroutine set_common_mpi_scale
  use scale_atmos_grid_cartesC, only: &
      ATMOS_GRID_CARTESC_CX, ATMOS_GRID_CARTESC_CY, &
      DX, DY
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  implicit none
  integer :: color, key
  integer :: ierr
  character(len=filelenmax) :: filename
  real(RP), allocatable :: height3dtmp(:,:,:)
  real(RP), allocatable :: lon2dtmp(:,:)
  real(RP), allocatable :: lat2dtmp(:,:)
  real(RP) :: lon2d_RP(nlon,nlat)
  real(RP) :: lat2d_RP(nlon,nlat)
  integer :: i, j
  real(r_size) :: ri, rj

  call mpi_timer('', 2)

  ! Communicator for 1-iteration ensemble member groups
  !-----------------------------------------------------------------------------

  color = myrank_to_pe
  key = myrank_to_mem(1) - 1

  if (myrank_use_da) then
    ! DA members
    call MPI_COMM_SPLIT(MPI_COMM_da, color, key, MPI_COMM_e, ierr)
  
    call MPI_COMM_SIZE(MPI_COMM_e, nprocs_e, ierr)
    call MPI_COMM_RANK(MPI_COMM_e, myrank_e, ierr)

  else
    ! No COMM_e for dacycle-forecast members
    return
!    ! exntended forecast members
!    call MPI_COMM_SPLIT(MPI_COMM_da, color, key, MPI_COMM_ef, ierr)
!  
!    call MPI_COMM_SIZE(MPI_COMM_ef, nprocs_ef, ierr)
!    call MPI_COMM_RANK(MPI_COMM_ef, myrank_ef, ierr)
!    
!    return
  endif

#ifdef DEBUG
  if (nprocs_e /= n_mem*n_mempn) then
    write (6, '(A)'), '[Error] XXXXXX wrong!!'
    stop
  end if
#endif

  call mpi_timer('set_common_mpi_scale:mpi_comm_split_e:', 2)

  ! Read/calculate model coordinates
  !-----------------------------------------------------------------------------

  if (VERIFY_COORD .and. (.not. DIRECT_TRANSFER)) then
    if (myrank_e == 0) then
!      allocate (height3d(nlev,nlon,nlat))
      allocate (lon2d(nlon,nlat))
      allocate (lat2d(nlon,nlat))
      allocate (height3dtmp(nlev,nlon,nlat))
      allocate (lon2dtmp(nlon,nlat))
      allocate (lat2dtmp(nlon,nlat))

      if (.not. allocated(topo2d)) then
        allocate (topo2d(nlon,nlat))
        if (FILE_AGGREGATE) then
#ifdef PNETCDF
          call read_topo_par(LETKF_TOPO_IN_BASENAME, topo2d, MPI_COMM_d)
#endif
        else
          call read_topo(LETKF_TOPO_IN_BASENAME, topo2d)
        end if
      end if
!      call scale_calc_z(topo2d, height3d)

!$OMP PARALLEL DO PRIVATE(i,j,ri,rj) COLLAPSE(2)
      do j = 1, nlat
        do i = 1, nlon
          ri = real(i + IHALO, r_size)
          rj = real(j + JHALO, r_size)
          call MAPPROJECTION_xy2lonlat(real((ri-1.0_r_size) * DX + ATMOS_GRID_CARTESC_CX(1),RP), &
                                       real((rj-1.0_r_size) * DY + ATMOS_GRID_CARTESC_CY(1),RP), lon2d_RP(i,j), lat2d_RP(i,j))
          lon2d(i,j) = lon2d_RP(i,j) * rad2deg
          lat2d(i,j) = lat2d_RP(i,j) * rad2deg
        end do
      end do
!$OMP END PARALLEL DO

      filename = trim(GUES_IN_BASENAME) // trim(timelabel_anal)
      call filename_replace_mem(filename, myrank_to_mem(1))
      call read_restart_coor(trim(filename), lon2dtmp, lat2dtmp, height3dtmp)

      if (maxval(abs(lon2dtmp - lon2d)) > 1.0d-6 .or. maxval(abs(lat2dtmp - lat2d)) > 1.0d-6) then
        write (6, '(A,F15.7,A,F15.7)') '[Error] Map projection settings are incorrect! -- maxdiff(lon) = ', &
                                       maxval(abs(lon2dtmp - lon2d)), ', maxdiff(lat) = ', maxval(abs(lat2dtmp - lat2d))
        stop
      end if
!      if (maxval(abs(height3dtmp - height3d)) > 1.0d-6) then
!        write (6, '(A,F15.7)') '[Error] 3D height calculation are incorrect, possibily due to inconsistent topography files! -- maxdiff(height) = ', &
!                               maxval(abs(height3dtmp - height3d))
!        stop
!      end if

      write (6, '(A)') 'VERIFY_COORD: Model coordinate calculation is good.'

      call mpi_timer('set_common_mpi_scale:verify_coord:', 2)
    end if
  end if

  return
end subroutine set_common_mpi_scale

!-------------------------------------------------------------------------------
! unset_common_mpi_scale
!-------------------------------------------------------------------------------
subroutine unset_common_mpi_scale
  implicit none
  integer:: ierr

  if (.not. myrank_use_da) return
  call MPI_COMM_FREE(MPI_COMM_e, ierr)

  return
end subroutine unset_common_mpi_scale

!-------------------------------------------------------------------------------
! set_common_mpi_grid
!-------------------------------------------------------------------------------
subroutine set_common_mpi_grid
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, &
    JHALO
  use scale_prc, only: &
    PRC_myrank

  implicit none
  REAL(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:)
  INTEGER :: i, j, n
  integer :: iproc, jproc
#ifdef DEBUG
  real(r_size) :: topo2dtmp(nlon,nlat)
#endif

  call mpi_timer('', 2)

  ! Compute nij1, nij1max, nij1node
  !-----------------------------------------------------------------------------

  i = mod(nlon*nlat, nprocs_e)
  nij1max = (nlon*nlat - i) / nprocs_e + 1
  if (myrank_e < i) then
    nij1 = nij1max
  else
    nij1 = nij1max - 1
  end if
!  write (6,'(A,I6.6,A,I7)') 'MYRANK ', myrank, ' number of grid points: nij1 =', nij1

  allocate (nij1node(nprocs_e))
  do n = 1, nprocs_e
    if (n-1 < i) then
      nij1node(n) = nij1max
    else
      nij1node(n) = nij1max - 1
    end if
  end do

  ALLOCATE(rig1(nij1))
  ALLOCATE(rjg1(nij1))
  ALLOCATE(topo1(nij1))

  ALLOCATE(hgt1(nij1,nlev))

  ALLOCATE(v3d(nij1,nlev,nv3d))
  ALLOCATE(v2d(nij1,nv2d))

  call mpi_timer('set_common_mpi_grid:nij1_cal:', 2)

!!!!!! ----- need to be replaced by more native communication !!!!!!

  if (myrank_e == mmean_rank_e) then
    v3dg = 0.0d0
    v2dg = 0.0d0

    call rank_1d_2d(PRC_myrank, iproc, jproc)
    do j = 1, nlat
      do i = 1, nlon
        v3dg(1,i,j,1) = real(i + iproc * nlon + IHALO, RP)
        v3dg(1,i,j,2) = real(j + jproc * nlat + JHALO, RP)
      end do
    end do

    call mpi_timer('set_common_mpi_grid:rij_cal:', 2)

    if (allocated(topo2d)) then
      write (6, '(1x,A,A15,A)') '*** Read 2D var: ', trim(topo2d_name), ' -- skipped because it was read previously'
#ifdef DEBUG
      if (.not. DIRECT_TRANSFER) then
        if (FILE_AGGREGATE) then
#ifdef PNETCDF
          call read_topo_par(LETKF_TOPO_IN_BASENAME, topo2dtmp, MPI_COMM_d)
#endif
        else
          call read_topo(LETKF_TOPO_IN_BASENAME, topo2dtmp)
        end if
        if (maxval(abs(topo2dtmp - topo2d)) > 1.0d-6) then
          write (6, '(A,F15.7)') '[Error] topo height in history files and restart files are inconsistent; maxdiff = ', maxval(abs(topo2dtmp - topo2d))
          stop
        end if
      end if
#endif
    else
      allocate (topo2d(nlon,nlat))
      if (DIRECT_TRANSFER) then
!        if (trim(LETKF_TOPO_IN_BASENAME) /= trim(TOPO_IN_BASENAME)) then
!          write (6, '(A)') '[Error] Direct transfer error: filenames mismatch.'
!          write (6, '(3A)') "        Filename in SCALE = '", trim(TOPO_IN_BASENAME), "'"
!          write (6, '(3A)') "        Filename in LETKF = '", trim(LETKF_TOPO_IN_BASENAME), "'"
!          stop
!        end if
        call read_topo(LETKF_TOPO_IN_BASENAME, topo2d)
        if (FILE_AGGREGATE) then
#ifdef PNETCDF
          call read_topo_par(LETKF_TOPO_IN_BASENAME, topo2d, MPI_COMM_d)
#endif
        else
          call read_topo(LETKF_TOPO_IN_BASENAME, topo2d)
        end if
      end if
    end if

    v3dg(1,:,:,3) = topo2d

    call mpi_timer('set_common_mpi_grid:read_topo:', 2)
  end if

  call mpi_timer('', 2, barrier=MPI_COMM_e)

  call scatter_grd_mpi(mmean_rank_e,v3dg,v2dg,v3d,v2d)

  rig1   = v3d(:,1,1)
  rjg1   = v3d(:,1,2)
  topo1  = v3d(:,1,3)

  call mpi_timer('set_common_mpi_grid:scatter:', 2)

  call scale_calc_z_grd(nij1, topo1, hgt1)

  call mpi_timer('set_common_mpi_grid:scale_calc_z:', 2)

  return
end subroutine set_common_mpi_grid

!-------------------------------------------------------------------------------
! set_mem_node_proc
!-------------------------------------------------------------------------------
SUBROUTINE set_mem_node_proc(mem)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: mem
  INTEGER :: tppn,tppnt,tmod
  INTEGER :: n,nn,m,q,qs,i,j,it,ip,ie

  call mpi_timer('', 2)

  if (mod(nprocs, PPN) /= 0) then
    write(6,'(A,I10)') '[Info] Total number of MPI processes      = ', nprocs
    write(6,'(A,I10)') '[Info] Number of processes per node (PPN) = ', PPN
    write(6,'(A)') '[Error] Total number of MPI processes should be an exact multiple of PPN.'
    stop
  end if
  nnodes = nprocs / PPN

  nprocs_m = sum(PRC_DOMAINS(1:NUM_DOMAIN))

  if (LOG_LEVEL >= 2 .and. myrank == 0) then
    write(6,'(A,I10)') '[Info] Total number of MPI processes                = ', nprocs
    write(6,'(A,I10)') '[Info] Number of nodes (NNODES)                     = ', nnodes
    write(6,'(A,I10)') '[Info] Number of processes per node (PPN)           = ', PPN
    write(6,'(A,I10)') '[Info] Number of processes per member (all domains) = ', nprocs_m
  end if

  if (MEM_NODES == 0) then
    MEM_NODES = (nprocs_m-1) / PPN + 1
  end if
  IF(MEM_NODES > 1) THEN
    n_mem = nnodes / MEM_NODES
    n_mempn = 1
  ELSE
    n_mem = nnodes
    n_mempn = PPN / nprocs_m
  END IF
  nitmax = (mem - 1) / (n_mem * n_mempn) + 1
  tppn = nprocs_m / MEM_NODES
  tmod = MOD(nprocs_m, MEM_NODES)

  ALLOCATE(mempe_to_node(nprocs_m,mem))
  ALLOCATE(mempe_to_rank(nprocs_m,mem))
  ALLOCATE(rank_to_mem(nitmax,nprocs))
  ALLOCATE(rank_to_pe(nprocs))
  ALLOCATE(rank_to_mempe(2,nitmax,nprocs))
  ALLOCATE(ranke_to_mem(nitmax,n_mem*n_mempn))
  ALLOCATE(myrank_to_mem(nitmax))

  rank_to_mem = -1
  rank_to_pe = -1
  rank_to_mempe = -1
  ranke_to_mem = -1
  m = 1
mem_loop: DO it = 1, nitmax
    ie = 1
    DO i = 0, n_mempn-1
      n = 0
      DO j = 0, n_mem-1
        IF(m > mem .and. it > 1) EXIT mem_loop
        qs = 0
        DO nn = 0, MEM_NODES-1
          IF(nn < tmod) THEN
            tppnt = tppn + 1
          ELSE
            tppnt = tppn
          END IF
          DO q = 0, tppnt-1
            ip = (n+nn)*PPN + i*nprocs_m + q
            if (m <= mem) then
              mempe_to_node(qs+1,m) = n+nn
              mempe_to_rank(qs+1,m) = ip
            end if
            rank_to_mem(it,ip+1) = m      ! These lines are outside of (m <= mem) condition
            if (it == 1) then             ! in order to cover over the entire first iteration
              rank_to_pe(ip+1) = qs       ! 
            end if                        ! 
            rank_to_mempe(1,it,ip+1) = m  ! 
            rank_to_mempe(2,it,ip+1) = qs ! 
            qs = qs + 1
          END DO
        END DO
        if (m <= mem) then
          ranke_to_mem(it,ie) = m
        end if
        ie = ie + 1
        m = m + 1
        n = n + MEM_NODES
      END DO
    END DO
  END DO mem_loop

  DO it = 1, nitmax
    myrank_to_mem(it) = rank_to_mem(it,myrank+1)
  END DO
  myrank_to_pe = rank_to_pe(myrank+1)

  if (myrank_to_mem(1) >= 1) then
    myrank_use = .true.
  end if

  ! settings related to mean/mdet (only valid when ENS_WITH_MEAN/ENS_WITH_MDET = .true.)
  !----------------------------------------------------------------
  if (ENS_WITH_MEAN) then
    mmean = MEMBER+1
    mmean_rank_e = mod(mmean-1, n_mem*n_mempn)
#ifdef DEBUG
    if (mmean_rank_e /= rank_to_mem(1,mempe_to_rank(1,mmean)+1)-1) then
      write (6, '(A)'), '[Error] XXXXXX wrong!!'
      stop
    end if
#endif
    msprd_rank_e = mmean_rank_e

    if (ENS_WITH_MDET) then
      nens = MEMBER+2
      nensobs = MEMBER+1

      mmdet = MEMBER+2
      mmdet_rank_e = mod(mmdet-1, n_mem*n_mempn)
#ifdef DEBUG
      if (mmdet_rank_e /= rank_to_mem(1,mempe_to_rank(1,mmdet)+1)-1) then
        write (6, '(A)'), '[Error] XXXXXX wrong!!'
        stop
      end if
#endif

      mmdetobs = MEMBER+1
      if (MDET_CYCLED) then
        mmdetin = mmdet
      else
        mmdetin = mmean
      end if
    else
      nens = MEMBER+1
      nensobs = MEMBER
    end if
  end if ! [ ENS_WITH_MEAN ]

  call mpi_timer('set_mem_node_proc:', 2)

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
! Start using SCALE library
!-------------------------------------------------------------------------------
subroutine set_scalelib(execname)
  use scale_io, only: &
    IO_setup, &
    IO_LOG_setup, &
    IO_FID_CONF, &
    IO_FID_LOG, &
    IO_L, &
    H_LONG
  use scale_prc, only: &
    PRC_mpi_alive, &
!    PRC_MPIstart, &
!    PRC_UNIVERSAL_setup, &
!    PRC_MPIsplit_letkf, &
    PRC_MPIsplit, &
    PRC_GLOBAL_setup, &
    PRC_LOCAL_setup, &
    PRC_UNIVERSAL_IsMaster, &
    PRC_nprocs, &
    PRC_myrank, &
    PRC_masterrank, &
    PRC_DOMAIN_nlim
  use scale_prc_cartesC, only: &
    PRC_2Drank, &
    PRC_CARTESC_setup
  use scale_const, only: &
    CONST_setup
  use scale_calendar, only: &
    CALENDAR_setup
  use scale_random, only: &
    RANDOM_setup
!  use scale_time, only: &
!    TIME_setup
  use scale_time, only: &
    TIME_DTSEC,       &
    TIME_STARTDAYSEC
  use scale_atmos_grid_cartesC, only: &
    ATMOS_GRID_CARTESC_setup, &
    ATMOS_GRID_CARTESC_DOMAIN_CENTER_X, &
    ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
    DX, &
    DY
  use scale_atmos_grid_cartesC_index
!  use scale_atmos_grid_cartesC_nest, only: &
!    NEST_setup
#ifdef PNETCDF
  use scale_land_grid_cartesC_index, only: &
    LAND_GRID_CARTESC_INDEX_setup
#endif
!  use scale_land_grid, only: &
!    LAND_GRID_setup
#ifdef PNETCDF
  use scale_urban_grid_cartesC_index, only: &
    URBAN_GRID_CARTESC_INDEX_setup
#endif
!  use scale_urban_grid, only: &
!    URBAN_GRID_setup
  use scale_file, only: &
    FILE_setup
  use scale_comm_cartesC, only: &
    COMM_setup
!  use scale_topography, only: &
!    TOPO_setup
!  use scale_landuse, only: &
!    LANDUSE_setup
!  use scale_atmos_grid_cartesC_real, only: &
!    REAL_setup
!  use scale_atmos_grid_cartesCtrans, only: &
!    GTRANS_setup
!  use scale_atmos_hydrostatic, only: &
!    ATMOS_HYDROSTATIC_setup
!  use scale_atmos_thermodyn, only: &
!    ATMOS_THERMODYN_setup
  use scale_atmos_hydrometeor, only: &
    ATMOS_HYDROMETEOR_setup
!  use mod_atmos_driver, only: &
!    ATMOS_driver_config
!  use scale_atmos_phy_mp, only: &
!    ATMOS_PHY_MP_config
!  use mod_atmos_admin, only: &
!    ATMOS_PHY_MP_TYPE, &
!    ATMOS_sw_phy_mp
!  use mod_user, only: &
!    USER_config
!  use mod_admin_time, only: &
!    ADMIN_TIME_setup
  use scale_mapprojection, only: &
    MAPPROJECTION_setup
  implicit none

  character(len=*), intent(in), optional :: execname

!  integer :: universal_comm
!  integer :: universal_nprocs
!  logical :: universal_master
  integer :: global_comm
  integer :: local_comm
  integer :: local_myrank
  logical :: local_ismaster
  integer :: intercomm_parent
  integer :: intercomm_child
  character(len=H_LONG) :: confname_domains(PRC_DOMAIN_nlim)
  character(len=H_LONG) :: confname_mydom

  integer :: color, key, idom, ierr
  integer :: rankidx(2)

  integer :: HIST_item_limit    ! dummy
  integer :: HIST_variant_limit ! dummy

  character(len=7) :: execname_ = ''

  if (present(execname)) execname_ = execname

  call mpi_timer('', 2, barrier=MPI_COMM_WORLD)

  ! Communicator for all processes used
  !-----------------------------------------------------------------------------

  if (myrank_use) then
    color = 0
    key   = (myrank_to_mem(1) - 1) * nprocs_m + myrank_to_pe
!    key   = myrank
  else
    color = MPI_UNDEFINED
    key   = MPI_UNDEFINED
  end if

  call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, MPI_COMM_u, ierr)

  if (.not. myrank_use) then
    write (6, '(A,I6.6,A)') 'MYRANK=', myrank, ': This process is not used!'
    return
  end if

  call MPI_COMM_SIZE(MPI_COMM_u, nprocs_u, ierr)
  call MPI_COMM_RANK(MPI_COMM_u, myrank_u, ierr)

  call mpi_timer('set_scalelib:mpi_comm_split_u:', 2)

  ! Communicator for all domains of single members
  !-----------------------------------------------------------------------------

  ! start SCALE MPI
!  call PRC_MPIstart( universal_comm ) ! [OUT]

  PRC_mpi_alive = .true.
!  universal_comm = MPI_COMM_u

!  call PRC_UNIVERSAL_setup( universal_comm,   & ! [IN]
!                            universal_nprocs, & ! [OUT]
!                            universal_master  ) ! [OUT]

!  if (myrank_to_mem(1) >= 1) then
    color = myrank_to_mem(1) - 1
    key   = myrank_to_pe
!  else
!    color = MPI_UNDEFINED
!    key   = MPI_UNDEFINED
!  endif

  call MPI_COMM_SPLIT(MPI_COMM_u, color, key, global_comm, ierr)

  call PRC_GLOBAL_setup( .false.,    & ! [IN]
                         global_comm ) ! [IN]

  call mpi_timer('set_scalelib:mpi_comm_split_d_global:', 2)

  ! Communicator for one domain
  !-----------------------------------------------------------------------------

  do idom = 1, NUM_DOMAIN
    confname_domains(idom) = trim(CONF_FILES)
    call filename_replace_dom(confname_domains(idom), idom)
  end do

  !--- split for nesting
  ! communicator split for nesting domains
  call PRC_MPIsplit( global_comm,      & ! [IN]
                     NUM_DOMAIN,       & ! [IN]
                     PRC_DOMAINS(:),   & ! [IN]
                     confname_domains(:), & ! [IN]
                     .false.,          & ! [IN]
                     .false.,          & ! [IN] flag bulk_split
                     .false.,          & ! [IN] no reordering
                     local_comm,       & ! [OUT]
                     intercomm_parent, & ! [OUT]
                     intercomm_child,  & ! [OUT]
                     confname_mydom    ) ! [OUT]

  MPI_COMM_d = local_comm

  do idom = 1, NUM_DOMAIN
    if (trim(confname_mydom) == trim(confname_domains(idom))) then
      mydom = idom
      exit
    end if
  end do

#ifdef DEBUG
  if (mydom <= 0) then
    write(6, '(A)'), '[Error] Cannot determine my domain ID.'
    stop
  end if
#endif

  !-----------------------------------------------------------------------------

  if (mydom >= 2) then ! In d01, keep using the original launcher config file; skip re-opening config files here
    call IO_setup( modelname, confname_mydom )

!    call read_nml_log
!    call read_nml_model
!    call read_nml_ensemble
!    call read_nml_process
  end if

  call PRC_LOCAL_setup( local_comm, local_myrank, local_ismaster )

!  call MPI_COMM_SIZE(MPI_COMM_d, nprocs_d, ierr)
  nprocs_d = PRC_nprocs
!  call MPI_COMM_RANK(MPI_COMM_d, myrank_d, ierr)
!  myrank_d = PRC_myrank
  myrank_d = local_myrank

  call mpi_timer('set_scalelib:mpi_comm_split_d_local:', 2)


  select case (execname_)
  case ('LETKF  ')
    call read_nml_obs_error
    call read_nml_obsope
    call read_nml_letkf
    call read_nml_letkf_obs
    call read_nml_letkf_var_local
    call read_nml_letkf_monitor
    call read_nml_letkf_radar
    call read_nml_letkf_h08
  case ('OBSOPE ', 'OBSMAKE')
    call read_nml_obs_error
    call read_nml_obsope
    call read_nml_letkf_radar
    call read_nml_letkf_h08
  case ('OBSSIM ')
    call read_nml_obssim
    call read_nml_letkf_radar
    call read_nml_letkf_h08
  end select

  ! Communicator for all processes for single domains
  !-----------------------------------------------------------------------------

!  if (mydom > 0) then
    color = mydom - 1
    key   = (myrank_to_mem(1) - 1) * nprocs_m + myrank_to_pe
!    key   = myrank
!  else
!    color = MPI_UNDEFINED
!    key   = MPI_UNDEFINED
!  end if

  call MPI_COMM_SPLIT(MPI_COMM_u, color, key, MPI_COMM_a, ierr)

  call MPI_COMM_SIZE(MPI_COMM_a, nprocs_a, ierr)
  call MPI_COMM_RANK(MPI_COMM_a, myrank_a, ierr)

  call mpi_timer('set_scalelib:mpi_comm_split_a:', 2)

  ! Setup scalelib LOG output (only for the universal master rank)
  !-----------------------------------------------------------------------------

  ! setup Log
  call IO_LOG_setup( local_myrank, PRC_UNIVERSAL_IsMaster )
!  call LogInit( IO_FID_CONF, IO_FID_LOG, IO_L )

  call mpi_timer('set_scalelib:log_setup_init:', 2)

  ! Other minimal scalelib setups for LETKF
  !-----------------------------------------------------------------------------

  ! setup process
  call PRC_CARTESC_setup

  ! setup PROF
!  call PROF_setup

  ! profiler start
!  call PROF_setprefx('INIT')
!  call PROF_rapstart('Initialize', 0)

  ! setup constants
  call CONST_setup

  ! setup calendar
!  call CALENDAR_setup

  ! setup random number
!  call RANDOM_setup

  ! setup time
!  call ADMIN_TIME_setup( setup_TimeIntegration = .true. )

  ! setup horizontal/vertical grid coordinates
  call ATMOS_GRID_CARTESC_INDEX_setup
  call ATMOS_GRID_CARTESC_setup
#ifdef PNETCDF
  call LAND_GRID_CARTESC_INDEX_setup
#endif
!  call LAND_GRID_setup
#ifdef PNETCDF
  call URBAN_GRID_CARTESC_INDEX_setup
#endif
!  call URBAN_GRID_setup

  ! setup tracer index
  call ATMOS_HYDROMETEOR_setup
!    call ATMOS_PHY_MP_config('TOMITA08') !!!!!!!!!!!!!!! tentative
!    if ( ATMOS_sw_phy_mp ) then
!       call ATMOS_PHY_MP_config( ATMOS_PHY_MP_TYPE )
!    end if
!  call ATMOS_driver_config
!  call USER_config

  ! setup file I/O
  call FILE_setup( PRC_myrank )

  ! setup mpi communication
  call COMM_setup

  ! setup topography
!  call TOPO_setup

  ! setup land use category index/fraction
!  call LANDUSE_setup

  ! setup grid coordinates (real world)
!  call REAL_setup
    ! setup map projection [[ in REAL_setup ]]
     call MAPPROJECTION_setup( ATMOS_GRID_CARTESC_DOMAIN_CENTER_X, ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y )

  ! setup grid transfer metrics (uses in ATMOS_dynamics)
!  call GTRANS_setup

  ! setup Z-ZS interpolation factor (uses in History)
!  call INTERP_setup

  ! setup restart
!  call ADMIN_restart_setup

  ! setup statistics
!  call STAT_setup

  ! setup history I/O
!  call HIST_setup
    ! setup history file I/O [[ in HIST_setup ]]
    rankidx(1) = PRC_2Drank(PRC_myrank, 1)
    rankidx(2) = PRC_2Drank(PRC_myrank, 2)

! tentative 11/26/2018 TH
!  call FILE_HISTORY_setup
!    call HistoryInit( HIST_item_limit,                  & ! [OUT]
!                      HIST_variant_limit,               & ! [OUT]
!                      IMAX, JMAX, KMAX,                 & ! [IN]
!                      PRC_masterrank,                   & ! [IN]
!                      PRC_myrank,                       & ! [IN]
!                      rankidx,                          & ! [IN]
!                      '',                               & ! [IN]
!                      '',                               & ! [IN]
!                      '',                               & ! [IN]
!                      0.0d0,                            & ! [IN]
!                      1.0d0,                            & ! [IN]
!                      default_basename='history',       & ! [IN]
!                      default_zcoord = 'model',         & ! [IN]
!                      default_tinterval = 1.0d0,        & ! [IN]
!                      namelist_fid=IO_FID_CONF          ) ! [IN]

  ! setup monitor I/O
!  call MONIT_setup

  ! setup nesting grid
!  call NEST_setup ( intercomm_parent, intercomm_child )

  ! setup common tools
!  call ATMOS_HYDROSTATIC_setup
!  call ATMOS_THERMODYN_setup
!  call ATMOS_SATURATION_setup
!  call BULKFLUX_setup
!  call ROUGHNESS_setup

  ! setup submodel administrator
!  call ATMOS_admin_setup
!  call OCEAN_admin_setup
!  call LAND_admin_setup
!  call URBAN_admin_setup
!  call CPL_admin_setup

  ! setup variable container
!  call ATMOS_vars_setup
!  call OCEAN_vars_setup
!  call LAND_vars_setup
!  call URBAN_vars_setup
!  call CPL_vars_setup

  call mpi_timer('set_scalelib:other_setup:', 2)

  return
end subroutine set_scalelib

!-------------------------------------------------------------------------------
! Finish using SCALE library
!-------------------------------------------------------------------------------
subroutine unset_scalelib
  use scale_file, only: &
     FILE_Close_All
  use scale_io, only: &
    IO_FID_CONF, &
    IO_FID_LOG, &
    IO_L, &
    IO_FID_STDOUT
  implicit none
  integer :: ierr

  if (myrank_use) then
!    call MONIT_finalize
    call FILE_Close_All

    ! Close logfile, configfile
    if ( IO_L ) then
      if( IO_FID_LOG /= IO_FID_STDOUT ) close(IO_FID_LOG)
    endif
    close(IO_FID_CONF)

    call MPI_COMM_FREE(MPI_COMM_d, ierr)
    call MPI_COMM_FREE(MPI_COMM_a, ierr)
    call MPI_COMM_FREE(MPI_COMM_u, ierr)
  end if

  return
end subroutine unset_scalelib

!-------------------------------------------------------------------------------
! Scatter gridded data to processes (nrank -> all)
!-------------------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(RP) :: bufs(nij1max,nlevall,nprocs_e)
  REAL(RP) :: bufr(nij1max,nlevall)
  INTEGER :: j,k,n,ierr,ns,nr

  ns = nij1max * nlevall
  nr = ns
  IF(myrank_e == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL grd_to_buf(nprocs_e,v3dg(k,:,:,n),bufs(:,j,:))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL grd_to_buf(nprocs_e,v2dg(:,:,n),bufs(:,j,:))
    END DO
  END IF

!  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_SCATTER(bufs,ns,COMM_datatype,&
                 & bufr,nr,COMM_datatype,nrank,MPI_COMM_e,ierr)

  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      v3d(:,k,n) = REAL(bufr(1:nij1,j),r_size)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
  END DO

!  CALL MPI_BARRIER(MPI_COMM_e,ierr)

  RETURN
END SUBROUTINE scatter_grd_mpi

!-------------------------------------------------------------------------------
! Gather gridded data (all -> nrank)
!-------------------------------------------------------------------------------
SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(RP) :: bufs(nij1max,nlevall)
  REAL(RP) :: bufr(nij1max,nlevall,nprocs_e)
  INTEGER :: j,k,n,ierr,ns,nr

  ns = nij1max * nlevall
  nr = ns
  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      bufs(1:nij1,j) = REAL(v3d(:,k,n),RP)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    bufs(1:nij1,j) = REAL(v2d(:,n),RP)
  END DO

!  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_GATHER(bufs,ns,COMM_datatype,&
                & bufr,nr,COMM_datatype,nrank,MPI_COMM_e,ierr)

  IF(myrank_e == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL buf_to_grd(nprocs_e,bufr(:,j,:),v3dg(k,:,:,n))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL buf_to_grd(nprocs_e,bufr(:,j,:),v2dg(:,:,n))
    END DO
  END IF

!  CALL MPI_BARRIER(MPI_COMM_e,ierr)

  RETURN
END SUBROUTINE gather_grd_mpi

!-------------------------------------------------------------------------------
! Read ensemble SCALE history files, one file per time (iter)
!-------------------------------------------------------------------------------
subroutine read_ens_history_iter(iter, step, v3dg, v2dg)
  implicit none
  integer, intent(in) :: iter
  integer, intent(in) :: step
  real(r_size), intent(out) :: v3dg(nlevh,nlonh,nlath,nv3dd)
  real(r_size), intent(out) :: v2dg(nlonh,nlath,nv2dd)
  character(filelenmax) :: filename
  integer :: im

  im = myrank_to_mem(iter)
  if (im >= 1 .and. im <= nens) then
    if (im <= MEMBER) then
      filename = trim(HISTORY_IN_BASENAME) // trim(timelabel_hist)
      call filename_replace_mem(filename, im)
    else if (im == mmean) then
      filename = trim(HISTORY_MEAN_IN_BASENAME) // trim(timelabel_hist)
    else if (im == mmdet) then
      filename = trim(HISTORY_MDET_IN_BASENAME) // trim(timelabel_hist)
    end if

    if (DIRECT_TRANSFER) then
!      !!!!!! Unable to check filename consistency here because FILE_HISTORY_DEFAULT_BASENAME and 
!      !!!!!! FILE_HISTORY_DEFAULT_POSTFIX_TIMELABEL are not public variables
!      if (trim(filename) /= trim(FILE_HISTORY_DEFAULT_BASENAME)//trim(timelabel_hist)) then
!        write (6, '(A)') '[Error] Direct transfer error: filenames mismatch.'
!        write (6, '(3A)') "        Output filename in SCALE = '", trim(FILE_HISTORY_DEFAULT_BASENAME)//trim(timelabel_hist), "'"
!        write (6, '(3A)') "        Input  filename in LETKF = '", trim(filename), "'"
!        stop
!      end if
      if (step /= 1) then
        write (6, '(A)') '[Error] Direct transfer of history data can only be done with 3D-LETKF (step = 1)'
        write (6, '(A,I5)') '        step =', step
        stop
      end if
      call read_history_direct(v3dg, v2dg)
    else
      if (FILE_AGGREGATE) then
#ifdef PNETCDF
        call read_history_par(trim(filename), step, v3dg, v2dg, MPI_COMM_d)
#endif
      else
        call read_history(trim(filename), step, v3dg, v2dg)
      endif
    endif ! DIRECT_TRANSFER
  end if

  return
end subroutine read_ens_history_iter

!-------------------------------------------------------------------------------
! Read ensemble first guess data and distribute to processes
!-------------------------------------------------------------------------------
subroutine read_ens_mpi(v3d, v2d)
  use mod_atmos_vars, only: &
    ATMOS_RESTART_OUT_BASENAME, &
    ATMOS_RESTART_OUT_POSTFIX_TIMELABEL
  implicit none
  real(r_size), intent(out) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(out) :: v2d(nij1,nens,nv2d)
  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  character(len=filelenmax) :: filename
  integer :: it, im, mstart, mend

  call mpi_timer('', 2)

  do it = 1, nitmax
    im = myrank_to_mem(it)

    ! Note: read all members + mdetin
    ! 
    if ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin) then
      if (im <= MEMBER) then
        filename = trim(GUES_IN_BASENAME) // trim(timelabel_anal)
        call filename_replace_mem(filename, im)
      else if (im == mmean) then
        filename = trim(GUES_MEAN_IN_BASENAME) // trim(timelabel_anal)
      else if (im == mmdet) then
        filename = trim(GUES_MDET_IN_BASENAME) // trim(timelabel_anal)
      end if

      if (DIRECT_TRANSFER) then
        if (ATMOS_RESTART_OUT_POSTFIX_TIMELABEL) then
          if ((myrank_a == 0) .and. trim(filename) /= trim(ATMOS_RESTART_OUT_BASENAME)//trim(timelabel_anal) .and. LOG_LEVEL >= 4) then
            write (6, '(A)') '[Error] Direct transfer error: filenames mismatch.'
            write (6, '(3A)') "        Output filename in SCALE = '", trim(ATMOS_RESTART_OUT_BASENAME)//trim(timelabel_anal), "'"
            write (6, '(3A)') "        Input  filename in LETKF = '", trim(filename), "'"
            stop
          end if
        else
          if ((myrank_a == 0) .and. trim(filename) /= trim(ATMOS_RESTART_OUT_BASENAME) .and. LOG_LEVEL >= 4) then
            write (6, '(A)') '[Error] Direct transfer error: filenames mismatch.'
            write (6, '(3A)') "        Output filename in SCALE = '", trim(ATMOS_RESTART_OUT_BASENAME), "'"
            write (6, '(3A)') "        Input  filename in LETKF = '", trim(filename), "'"
            stop
          end if
        endif ! ATMOS_RESTART_OUT_POSTFIX_TIMELABEL

        call read_restart_direct(v3dg, v2dg)

      else
!        write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',trim(filename),'.pe',myrank_d,'.nc'
        if (FILE_AGGREGATE) then
#ifdef PNETCDF
          call read_restart_par(trim(filename), v3dg, v2dg, MPI_COMM_d)
#endif
        else
          call read_restart(trim(filename), v3dg, v2dg)
        end if
      end if

      call mpi_timer('read_ens_mpi:read_restart:', 2)

      call state_trans(v3dg)

      call mpi_timer('read_ens_mpi:state_trans:', 2)
    else if (im <= nens) then ! This is to avoid the undefined value problem;
      v3dg = undef            ! it has no impact to the results
      v2dg = undef            ! 
    end if

    call mpi_timer('', 2, barrier=MPI_COMM_e)

    mstart = 1 + (it-1)*nprocs_e
    mend = min(it*nprocs_e, nens)
    if (mstart <= mend) then
      call scatter_grd_mpi_alltoall(mstart, mend, v3dg, v2dg, v3d, v2d)
    end if

    call mpi_timer('read_ens_mpi:scatter_grd_mpi_alltoall:', 2)
  end do ! [ it = 1, nitmax ]

  return
end subroutine read_ens_mpi

!-------------------------------------------------------------------------------
! Read ensemble additive inflation parameter and distribute to processes
!-------------------------------------------------------------------------------
subroutine read_ens_mpi_addiinfl(v3d, v2d)
  implicit none
  real(r_size), intent(out) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(out) :: v2d(nij1,nens,nv2d)
  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  character(len=filelenmax) :: filename
  integer :: it, im, mstart, mend

  do it = 1, nitmax
    im = myrank_to_mem(it)

    ! Note: read all members
    ! 
    if (im >= 1 .and. im <= MEMBER) then
      filename = INFL_ADD_IN_BASENAME
      call filename_replace_mem(filename, im)

!      write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',myrank_d,'.nc'
      if (FILE_AGGREGATE) then
#ifdef PNETCDF
        call read_restart_par(filename, v3dg, v2dg, MPI_COMM_d)
#endif
      else
      call read_restart(filename, v3dg, v2dg)
      end if
!      call state_trans(v3dg)
    end if

    mstart = 1 + (it-1)*nprocs_e
    mend = min(it*nprocs_e, MEMBER)
    if (mstart <= mend) then
      call scatter_grd_mpi_alltoall(mstart, mend, v3dg, v2dg, v3d, v2d)
    end if
  end do ! [ it = 1, nitmax ]

  return
end subroutine read_ens_mpi_addiinfl

!-------------------------------------------------------------------------------
! Write ensemble analysis data after collecting from processes
!-------------------------------------------------------------------------------
subroutine write_ens_mpi(v3d, v2d, mean3d, mean2d)
  use mod_atmos_vars, only: &
    ATMOS_RESTART_IN_BASENAME, &
    ATMOS_RESTART_IN_POSTFIX_TIMELABEL
  implicit none
  real(r_size), intent(in) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(in) :: v2d(nij1,nens,nv2d)
  real(RP), intent(out), optional :: mean3d(nlev,nlon,nlat,nv3d)
  real(RP), intent(out), optional :: mean2d(nlon,nlat,nv2d)

  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  character(len=filelenmax) :: filename
  integer :: it, im, mstart, mend


  call mpi_timer('', 2)

  do it = 1, nitmax
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    im = myrank_to_mem(it)

    mstart = 1 + (it-1)*nprocs_e
    mend = min(it*nprocs_e, nens)
    if (mstart <= mend) then
      call gather_grd_mpi_alltoall(mstart, mend, v3d, v2d, v3dg, v2dg)
    end if

    call mpi_timer('write_ens_mpi:gather_grd_mpi_alltoall:', 2)

    ! Note: write all members + mean + mdet
    ! 
    if ((im >= 1 .and. im <= MEMBER) .or. im == mmean .or. im == mmdet) then
      if (im <= MEMBER) then
        filename = trim(ANAL_OUT_BASENAME) // trim(timelabel_anal)
        call filename_replace_mem(filename, im)
      else if (im == mmean) then
        filename = trim(ANAL_MEAN_OUT_BASENAME) // trim(timelabel_anal)
        if (present(mean3d) .and. nv3d > 0) then
          mean3d = v3dg
        end if
        if (present(mean2d) .and. nv2d > 0) then
          mean2d = v2dg
        end if
      else if (im == mmdet) then
        filename = trim(ANAL_MDET_OUT_BASENAME) // trim(timelabel_anal)
      end if

      call state_trans_inv(v3dg)

      call mpi_timer('write_ens_mpi:state_trans_inv:', 2)

      if (DIRECT_TRANSFER) then
        if (ATMOS_RESTART_IN_POSTFIX_TIMELABEL) then
          if (trim(filename) /= trim(ATMOS_RESTART_IN_BASENAME)//trim(timelabel_anal)) then

            if (LOG_LEVEL >= 4 .and. myrank_da == 0) then
              write (6, '(A)') '[Error] Direct transfer error: filenames mismatch.'
              write (6, '(3A)') "        Output filename in LETKF = '", trim(filename), "'"
              write (6, '(3A)') "        Input  filename in SCALE = '", trim(ATMOS_RESTART_IN_BASENAME)//trim(timelabel_anal), "'"
            endif
            ! File names can be different if INDIR and OUTDIR are different
            !stop
          endif
        else
          if (trim(filename) /= trim(ATMOS_RESTART_IN_BASENAME)) then

            if (LOG_LEVEL >= 4 .and. myrank_da == 0) then
                write (6, '(A)') '[Error] Direct transfer error: filenames mismatch.'
                write (6, '(3A)') "        Output filename in LETKF = '", trim(filename), "'"
                write (6, '(3A)') "        Input  filename in SCALE = '", trim(ATMOS_RESTART_IN_BASENAME), "'"
            endif

            ! File names can be different if INDIR and OUTDIR are different
            !stop
          end if
        end if
 
        call write_restart_direct(v3dg, v2dg) 

      else
        if (FILE_AGGREGATE) then
#ifdef PNETCDF
          call write_restart_par(filename, v3dg, v2dg, MPI_COMM_d)
#endif
        else
          call write_restart(trim(filename), v3dg, v2dg)
        end if
      end if ! DIRECT_TRANSFER

      call mpi_timer('write_ens_mpi:write_restart:', 2)
    end if
  end do ! [ it = 1, nitmax ]

  return
end subroutine write_ens_mpi

!-------------------------------------------------------------------------------
! Scatter gridded data using MPI_ALLTOALL(V) (mstart~mend -> all)
!-------------------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi_alltoall(mstart,mend,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: mstart,mend
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(INOUT) :: v3d(nij1,nlev,nens,nv3d)
  REAL(r_size),INTENT(INOUT) :: v2d(nij1,nens,nv2d)
  REAL(RP) :: bufs(nij1max,nlevall,nprocs_e)
  REAL(RP) :: bufr(nij1max,nlevall,nprocs_e)
  INTEGER :: k,n,j,m,mcount,ierr
  INTEGER :: ns(nprocs_e),nst(nprocs_e),nr(nprocs_e),nrt(nprocs_e)

  mcount = mend - mstart + 1
#ifdef DEBUG
  IF(mcount > nprocs_e .OR. mcount <= 0) STOP
#endif

  IF(myrank_e < mcount) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL grd_to_buf(nprocs_e,v3dg(k,:,:,n),bufs(:,j,:))
      END DO
    END DO
    DO n=1,nv2d
      j = j+1
      CALL grd_to_buf(nprocs_e,v2dg(:,:,n),bufs(:,j,:))
    END DO
  END IF

!  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  IF(mcount == nprocs_e) THEN
    CALL MPI_ALLTOALL(bufs, nij1max*nlevall, COMM_datatype, &
                      bufr, nij1max*nlevall, COMM_datatype, MPI_COMM_e, ierr)
  ELSE
    CALL set_alltoallv_counts(mcount,nij1max*nlevall,nprocs_e,nr,nrt,ns,nst)
    CALL MPI_ALLTOALLV(bufs, ns, nst, COMM_datatype, &
                       bufr, nr, nrt, COMM_datatype, MPI_COMM_e, ierr)
  END IF

  DO m = mstart,mend
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        v3d(:,k,m,n) = REAL(bufr(1:nij1,j,m-mstart+1),r_size)
      END DO
    END DO
    DO n=1,nv2d
      j = j+1
      v2d(:,m,n) = REAL(bufr(1:nij1,j,m-mstart+1),r_size)
    END DO
  END DO

!  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  RETURN
END SUBROUTINE scatter_grd_mpi_alltoall

!-------------------------------------------------------------------------------
! Gather gridded data using MPI_ALLTOALL(V) (all -> mstart~mend)
!-------------------------------------------------------------------------------
SUBROUTINE gather_grd_mpi_alltoall(mstart,mend,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: mstart,mend
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nens,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nens,nv2d)
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(RP) :: bufs(nij1max,nlevall,nprocs_e)
  REAL(RP) :: bufr(nij1max,nlevall,nprocs_e)
  INTEGER :: k,n,j,m,mcount,ierr
  INTEGER :: ns(nprocs_e),nst(nprocs_e),nr(nprocs_e),nrt(nprocs_e)

  mcount = mend - mstart + 1
#ifdef DEBUG
  IF(mcount > nprocs_e .OR. mcount <= 0) STOP
#endif

  DO m = mstart,mend
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        bufs(1:nij1,j,m-mstart+1) = REAL(v3d(:,k,m,n),RP)
      END DO
    END DO
    DO n=1,nv2d
      j = j+1
      bufs(1:nij1,j,m-mstart+1) = REAL(v2d(:,m,n),RP)
    END DO
  END DO

!  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  IF(mcount == nprocs_e) THEN
    CALL MPI_ALLTOALL(bufs, nij1max*nlevall, COMM_datatype, &
                      bufr, nij1max*nlevall, COMM_datatype, MPI_COMM_e, ierr)
  ELSE
    CALL set_alltoallv_counts(mcount,nij1max*nlevall,nprocs_e,ns,nst,nr,nrt)
    CALL MPI_ALLTOALLV(bufs, ns, nst, COMM_datatype, &
                       bufr, nr, nrt, COMM_datatype, MPI_COMM_e, ierr)
  END IF

  IF(myrank_e < mcount) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL buf_to_grd(nprocs_e,bufr(:,j,:),v3dg(k,:,:,n))
      END DO
    END DO
    DO n=1,nv2d
      j = j+1
      CALL buf_to_grd(nprocs_e,bufr(:,j,:),v2dg(:,:,n))
    END DO
  END IF

!  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  RETURN
END SUBROUTINE gather_grd_mpi_alltoall

!-------------------------------------------------------------------------------
! Set the send/recieve counts of MPI_ALLTOALLV
!-------------------------------------------------------------------------------
SUBROUTINE set_alltoallv_counts(mcount,ngpblock,np,n_ens,nt_ens,n_mem,nt_mem)
  INTEGER,INTENT(IN) :: mcount,ngpblock
  INTEGER,INTENT(IN) :: np
  INTEGER,INTENT(OUT) :: n_ens(np),nt_ens(np),n_mem(np),nt_mem(np)
  INTEGER :: p

  n_ens = 0
  nt_ens = 0
  n_mem = 0
  nt_mem = 0
  DO p=1,mcount
    n_ens(p) = ngpblock
    IF(myrank_e+1 == p) THEN
      n_mem(:) = ngpblock
    END IF
  END DO
  DO p=2,np
    nt_ens(p) = nt_ens(p-1) + n_ens(p-1)
    nt_mem(p) = nt_mem(p-1) + n_mem(p-1)
  END DO

  RETURN
END SUBROUTINE set_alltoallv_counts

!-------------------------------------------------------------------------------
! gridded data -> buffer
!-------------------------------------------------------------------------------
SUBROUTINE grd_to_buf(np,grd,buf)
  INTEGER,INTENT(IN) :: np
  REAL(RP),INTENT(IN) :: grd(nlon,nlat)
  REAL(RP),INTENT(OUT) :: buf(nij1max,np)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,np
    DO i=1,nij1node(m)
      j = m-1 + np * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
#ifdef DEBUG
if (i < 1 .or. i > nij1max .or. m < 1 .or. m > np .or. ilon < 1 .or. ilon > nlon .or. ilat < 1 .or. ilat > nlat) then
  write(6, *), '[Error] ######', np, nij1max
  write(6, *), '[Error] ######', i, m, ilon, ilat
  stop
end if
#endif
      buf(i,m) = grd(ilon,ilat)
    END DO
  END DO

  DO m=1,np
    IF(nij1node(m) < nij1max) buf(nij1max,m) = undef
  END DO

  RETURN
END SUBROUTINE grd_to_buf

!-------------------------------------------------------------------------------
! buffer -> gridded data
!-------------------------------------------------------------------------------
SUBROUTINE buf_to_grd(np,buf,grd)
  INTEGER,INTENT(IN) :: np
  REAL(RP),INTENT(IN) :: buf(nij1max,np)
  REAL(RP),INTENT(OUT) :: grd(nlon,nlat)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,np
    DO i=1,nij1node(m)
      j = m-1 + np * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      grd(ilon,ilat) = buf(i,m)
    END DO
  END DO

  RETURN
END SUBROUTINE buf_to_grd

!-------------------------------------------------------------------------------
! Gather ensemble mean to {mmean_rank_e} and write SCALE restart files
!-------------------------------------------------------------------------------
subroutine write_ensmean(filename, v3d, v2d, calced, mean_out, mean3d, mean2d)
  implicit none
  character(len=*), intent(in) :: filename
  real(r_size), intent(inout) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(inout) :: v2d(nij1,nens,nv2d)
  logical, intent(in), optional :: calced
  logical, intent(in), optional :: mean_out
  real(RP), intent(out), optional :: mean3d(nlev,nlon,nlat,nv3d)
  real(RP), intent(out), optional :: mean2d(nlon,nlat,nv2d)

  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  logical :: calced_
  logical :: mean_out_

  call mpi_timer('', 2)

  calced_ = .false.
  if (present(calced)) then
    calced_ = calced
  end if
  mean_out_ = .true.
  if (present(mean_out)) then
    mean_out_ = mean_out
  end if

  if (.not. calced_) then
    call ensmean_grd(MEMBER, nens, nij1, v3d, v2d)

    call mpi_timer('write_ensmean:ensmean_grd:', 2)
  end if

  call mpi_timer('', 2, barrier=MPI_COMM_e)

  call gather_grd_mpi(mmean_rank_e, v3d(:,:,mmean,:), v2d(:,mmean,:), v3dg, v2dg)

  call mpi_timer('write_ensmean:gather_grd_mpi:', 2)

  if (myrank_e == mmean_rank_e) then
    if (present(mean3d) .and. nv3d > 0) then
      mean3d = v3dg
    end if
    if (present(mean2d) .and. nv2d > 0) then
      mean2d = v2dg
    end if

    call state_trans_inv(v3dg)

    call mpi_timer('write_ensmean:state_trans_inv:', 2)

    if (mean_out_) then
      if (FILE_AGGREGATE) then
#ifdef PNETCDF
        call write_restart_par(filename, v3dg, v2dg, MPI_COMM_d)
#endif
      else
        call write_restart(filename, v3dg, v2dg)
      end if
    end if

    call mpi_timer('write_ensmean:write_restart:', 2)
  end if

  return
end subroutine write_ensmean

!-------------------------------------------------------------------------------
! Gather ensemble spread to {msprd_rank_e} and write SCALE restart files
!-------------------------------------------------------------------------------
subroutine write_enssprd(filename, v3d, v2d)
  implicit none
  character(len=*), intent(in) :: filename
  real(r_size), intent(in) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(in) :: v2d(nij1,nens,nv2d)
  real(r_size) :: v3ds(nij1,nlev,nv3d)
  real(r_size) :: v2ds(nij1,nv2d)
  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)

  call mpi_timer('', 2)

  call enssprd_grd(MEMBER, nens, nij1, v3d, v2d, v3ds, v2ds)

  call mpi_timer('write_enssprd:enssprd_grd:', 2, barrier=MPI_COMM_e)

  call gather_grd_mpi(msprd_rank_e, v3ds, v2ds, v3dg, v2dg)

  call mpi_timer('write_enssprd:gather_grd_mpi:', 2)

  if (myrank_e == msprd_rank_e) then
!    call state_trans_inv(v3dg)              !! do not transform the spread output
    if (FILE_AGGREGATE) then
#ifdef PNETCDF
      call write_restart_par(filename, v3dg, v2dg, MPI_COMM_d)
#endif
    else
      call write_restart(filename, v3dg, v2dg) !!
    end if

    call mpi_timer('write_enssprd:write_restart:', 2)
  end if

  return
end subroutine write_enssprd

!-------------------------------------------------------------------------------
! MPI timer
!-------------------------------------------------------------------------------
subroutine mpi_timer(sect_name, level, barrier)
  implicit none
  character(len=*), intent(in) :: sect_name
  integer, intent(in) :: level
  integer, intent(in), optional :: barrier

  character(len=timer_name_width) :: sect_name_tmp
  character(len=14) :: sect_prefix_1
  character(len=12) :: sect_prefix_2
  real(r_dble) :: timer_before_barrier
  real(r_dble) :: timer_after_barrier
  integer :: i, ierr
  logical :: initialized

  if (LOG_LEVEL < level ) return

  timer_before_barrier = MPI_WTIME()
  timer_after_barrier = timer_before_barrier

  if (USE_MPI_BARRIER .and. present(barrier)) then
    if (barrier /= MPI_COMM_NULL) then
      call MPI_BARRIER(barrier, ierr)
      timer_after_barrier = MPI_WTIME()
    end if
  end if

  initialized = .true.
  if (timer_save(level) < 0.0d0) then
    initialized = .false.
    do i = level-1, 1, -1
      if (timer_save(i) >= 0.0d0) then
        timer_save(level) = timer_save(i)
        exit
      end if
    end do
  end if

  do i = max_timer_levels, level, -1
    if (timer_save(i) >= 0.0d0) then
      select case (i)
      case (1)
        sect_prefix_1 = '##### TIMER # '
        sect_prefix_2 = ''
      case (2)
        sect_prefix_1 = ' #### TIMER # '
        sect_prefix_2 = '...'
      case (3)
        sect_prefix_1 = '  ### TIMER # '
        sect_prefix_2 = '......'
      case (4)
        sect_prefix_1 = '   ## TIMER # '
        sect_prefix_2 = '.........'
      case (5)
        sect_prefix_1 = '    # TIMER # '
        sect_prefix_2 = '............'
      end select

      if (myrank == 0) then

        if (i == level .and. initialized .and. trim(sect_name) /= '') then
          sect_name_tmp = sect_name ! to left-align the text
          write (6,'(3A,2F14.6,A)') sect_prefix_1, trim(sect_prefix_2), sect_name_tmp, &
                                    timer_before_barrier - timer_save(i), &
                                    timer_after_barrier - timer_save(i)
        else if (timer_after_barrier - timer_save(i) >= timer_neglect) then
          if (i == level .and. initialized) then
            sect_name_tmp = ' (wait)'
          else
            sect_name_tmp = ' (unknown)'
          end if
          if (timer_before_barrier - timer_save(i) >= timer_neglect) then
            write (6,'(3A,2F14.6,A)') sect_prefix_1, trim(sect_prefix_2), sect_name_tmp, &
                                      timer_before_barrier - timer_save(i), &
                                      timer_after_barrier - timer_save(i)
          else
            write (6,'(3A,14x,F14.6,A)') sect_prefix_1, trim(sect_prefix_2), sect_name_tmp, &
                                         timer_after_barrier - timer_save(i)
          end if
        end if
      end if
    end if

    if (i == level) then
      timer_save(i) = timer_after_barrier
    else
      timer_save(i) = -9.0d10 ! reset the timer for all levels under this level
    end if

  end do

  return
end subroutine mpi_timer

!-------------------------------------------------------------------------------
! [Direct transfer] Send SCALE restart (analysis) data
!-------------------------------------------------------------------------------
subroutine send_emean_direct(v3dg,v2dg,fcst_cnt)
  implicit none

  real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP), intent(in) :: v2dg(nlon,nlat,nv2d)
  integer, intent(in) :: fcst_cnt
  integer :: iv3d, iv2d
  integer :: k, i, j

  integer :: rrank_a ! rank_a of receiving rank
  integer :: srank_a ! rank_a of sending rank ! Ensemble mean
  integer :: ierr, tag

  rrank_a = nens * nprocs_d + nprocs_d * (fcst_cnt - 1) + myrank_d ! dacycle-forecast member
  srank_a = MEMBER * nprocs_d + myrank_d

  if ( myrank_a /= srank_a ) return

  tag = myrank_d 

  if (nv3d > 0) then
    call MPI_Send(v3dg,nlev*nlon*nlat*nv3d,MPI_RP,rrank_a,tag,MPI_COMM_a,ierr) 
  endif

  if (nv2d > 0) then
    call MPI_Send(v2dg,nlon*nlat*nv2d,MPI_RP,rrank_a,tag+1,MPI_COMM_a,ierr) 
  endif

  return
end subroutine send_emean_direct

!-------------------------------------------------------------------------------
! [Direct transfer] Receive SCALE restart (analysis) data
!-------------------------------------------------------------------------------
subroutine receive_emean_direct()
  use mod_atmos_vars, only: &
    DENS, &
    MOMX, &
    MOMY, &
    MOMZ, &
    RHOT, &
    QTRC, QV, Qe, &
    ATMOS_vars_fillhalo
  use mod_atmos_phy_mp_driver, only: &
    ATMOS_PHY_MP_driver_qhyd2qtrc
  use mod_atmos_phy_mp_vars, only: &
    QS_MP, QE_MP
  use scale_atmos_hydrometeor, only: &
    I_HC, I_HR, I_HI, I_HS, I_HG
  use scale_atmos_grid_cartesC_index, only: &
    IS, IE, IA, JS, JE, JA, KS, KE, KA
  implicit none

  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  integer :: iv3d, iv2d

  integer :: srank_a ! rank_a of sending rank ! Ensemble mean
  integer :: ierr, tag
  integer, allocatable :: istat(:)

  srank_a = MEMBER * nprocs_d + myrank_d

  tag = myrank_d 

  allocate(istat(MPI_STATUS_SIZE))

  if (nv3d > 0) then
    call MPI_Recv(v3dg,nlev*nlon*nlat*nv3d,MPI_RP,srank_a,tag,MPI_COMM_a,istat,ierr)
    ! Ensemble mean data (mean3d & mean 2d) are not transformed into 
    ! SCALE prognostic variables (e.g., MOMX) in write_ens_mpi
    call state_trans_inv(v3dg) 
  endif

  if (nv2d > 0) then
    call MPI_Recv(v2dg,nlon*nlat*nv2d,MPI_RP,srank_a,tag+1,MPI_COMM_a,istat,ierr)
  endif
  
  deallocate(istat)

  do iv3d = 1, nv3d
    if (LOG_LEVEL >= 3) then
      write(6,'(1x,A,A15)') '*** Update 3D var [direct transfer]: ', trim(v3d_name(iv3d))
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

  do iv2d = 1, nv2d
    if (LOG_LEVEL >= 3) then
      write(6,'(1x,A,A15)') '*** Update 2D var [direct transfer]: ', trim(v2d_name(iv2d))
    end if
!    select case (iv2d)
!    case default
!      write (6, '(3A)') "[Error] Variable '", trim(v2d_name(iv2d)), "' is not
!      recognized."
!      stop
!    end select
  end do

  return
end subroutine receive_emean_direct
!-------------------------------------------------------------------------------
! Write the subdomain model data into a single GrADS file from DACYCLE (additional) forecasts
!-------------------------------------------------------------------------------
subroutine write_grd_dafcst_mpi(timelabel, ref3d, step)
  use mod_atmos_vars, only: &
    TEMP
!  use scale_atmos_hydrometeor, only: &
!    I_QV, I_HC, I_HR, I_HI, I_HS, I_HG
  use scale_atmos_grid_cartesC_index, only: &
    IS, IE, JS, JE, KS, KE, &
    KHALO
  use scale_io, only: &
    H_LONG

  implicit none
  character(15), intent(in) :: timelabel
  real(r_size), intent(in) :: ref3d(nlev,nlon,nlat)
  integer, intent(in) :: step

  character(len=H_LONG) :: filename
  real(r_sngl) :: bufs4(nlong,nlatg)
  real(r_sngl) :: bufr4(nlong,nlatg)
  integer :: iunit, iolen
  integer :: k, n, irec, ierr
  integer :: proc_i, proc_j
  integer :: ishift, jshift

!  real(r_sngl) :: v2d_ref(nlong,nlatg,nv3dd)

  call rank_1d_2d(myrank_d, proc_i, proc_j)
  ishift = proc_i * nlon
  jshift = proc_j * nlat

  if (myrank_d == 0) then
    filename = trim(DACYCLE_RUN_FCST_OUTNAME)//"/fcst_ref3d_"//trim(timelabel)//".grd"
    iunit = 55
    inquire (iolength=iolen) iolen
    open (iunit, file=trim(filename), form='unformatted', access='direct', &
          status='unknown', convert='native', recl=nlong*nlatg*iolen)
    irec = (step - 1)*nlev*2 ! 2 variable (nlev*2 record) output 
  end if

  do k = 1, nlev 
    bufs4(:,:) = 0.0
    bufs4(1+ishift:nlon+ishift, 1+jshift:nlat+jshift) = real(ref3d(k,1:nlon,1:nlat), r_sngl)
    call MPI_REDUCE(bufs4, bufr4, nlong*nlatg, MPI_REAL, MPI_SUM, 0, MPI_COMM_d, ierr)

    if (myrank_d == 0) then
      irec = irec + 1
      write (iunit, rec=irec) bufr4
    end if

  enddo

  ! debug
  do k = 1, nlev 
    bufs4(:,:) = 0.0
    bufs4(1+ishift:nlon+ishift, 1+jshift:nlat+jshift) = real(TEMP(KHALO+k,IS:IE,JS:JE), r_sngl)
    call MPI_REDUCE(bufs4, bufr4, nlong*nlatg, MPI_REAL, MPI_SUM, 0, MPI_COMM_d, ierr)

    if (myrank_d == 0) then
      irec = irec + 1
      write (iunit, rec=irec) bufr4
    end if

  enddo

  if (myrank_d == 0) then
    close (iunit)
  end if

  return
end subroutine write_grd_dafcst_mpi

!-------------------------------------------------------------------------------
! Write the analysis/guess ensemble mean into a single GrADS file 
!-------------------------------------------------------------------------------
subroutine write_grd_all_mpi(timelabel, v3dg, step)
  use scale_io, only: &
    H_LONG

  implicit none
  character(15), intent(in) :: timelabel
  real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  integer, intent(in) :: step ! 1:guess, 2:analysis

  character(4) :: head
  character(len=H_LONG) :: filename
  real(r_sngl) :: buf(nlev,nlong,nlatg)
  integer :: iunit, iolen
  integer :: k, n, irec, ierr
  integer :: proc_i, proc_j
  integer :: ishift, jshift

  integer :: iv3d

  call rank_1d_2d(myrank_d, proc_i, proc_j)
  ishift = proc_i * nlon
  jshift = proc_j * nlat

  if (myrank_d == 0) then

    if ( step == 1 ) head = "gues_"
    if ( step == 2 ) head = "anal_"

    filename = trim(OUT_GRADS_DA_ALL_PATH) // "/" //head // trim(timelabel) // ".grd"
    iunit = 55
    inquire (iolength=iolen) iolen
    open (iunit, file=trim(filename), form='unformatted', access='direct', &
          status='unknown', convert='native', recl=nlong*nlatg*iolen)
    irec = 0
  end if

  do iv3d = 1, nv3d
    buf(:,:,:) = 0.0
    buf(:,1+ishift:nlon+ishift, 1+jshift:nlat+jshift) = real( v3dg(:,:,:,iv3d), kind=r_sngl )
    call MPI_ALLREDUCE(MPI_IN_PLACE, buf, nlev*nlong*nlatg, MPI_RP, MPI_SUM, MPI_COMM_d, ierr)

    if (myrank_d == 0) then
      do k = 1, nlev
        irec = irec + 1
        write (iunit, rec=irec) buf(k,:,:)
      enddo
    end if

  enddo

  if (myrank_d == 0) then
    close (iunit)
  end if

  return
end subroutine write_grd_all_mpi

#ifdef PLOT_DCL
!-------------------------------------------------------------------------------
! Plot forecast/analysis 3D data by DCL
!-------------------------------------------------------------------------------
subroutine plot_dafcst_mpi(timelabel, ref3d, step)
  use mod_admin_time, only: &
    TIME_DTSEC_ATMOS_RESTART
  use scale_atmos_grid_cartesC, only: &
     CZ => ATMOS_GRID_CARTESC_CZ
!  use scale_atmos_hydrometeor, only: &
!    I_QV, I_HC, I_HR, I_HI, I_HS, I_HG
  use scale_atmos_grid_cartesC_index, only: &
    IS, IE, JS, JE, KS, KE, &
    KHALO
  use scale_io, only: &
    H_LONG

  implicit none
  character(15), intent(in) :: timelabel
  real(r_size), intent(in) :: ref3d(nlev,nlon,nlat)
  integer, intent(in), optional :: step

  character(len=H_LONG) :: filename
  integer :: nlev_plot
  real(r_sngl), allocatable :: bufr3d(:,:,:)
  integer :: iunit, iolen
  integer :: k, n, irec, ierr
  integer :: proc_i, proc_j
  integer :: ishift, jshift
  character(4) :: ftsec ! forecast time (second)
  character(5) :: cheight ! height (m)

  character(len=8) :: date
  character(len=10) :: time

  character(len=H_LONG) :: plotname

  logical :: fcst_
  character(len=4) :: header
  character(len=8) :: footer_fcst

  call mpi_timer('', 2)

  fcst_ = .false.
  header = "anal"
  footer_fcst = ""
  ftsec="anal"
  if (present(step)) then
    fcst_ = .true.
    header = "fcst"
    write(ftsec,'(I4.4)')  (step - 1) * int(TIME_DTSEC_ATMOS_RESTART)
    footer_fcst = "_FT" // ftsec // "s"
  end if

!  call MPI_BARRIER(MPI_COMM_d, ierr)
!  call date_and_time(date=date, time=time)
!  if (myrank_d == 0) then
!    write (6, '(3a,1x,a,1x,a,1x,a)') '[Info:plot] ',header,' start plotting: ', date, time, trim(timelabel) // trim(footer_fcst)
!  endif


  ! Count the number of plot levels
  nlev_plot = 0
  do k = plot_zlev_min, plot_zlev_max, plot_zlev_intv

    if (CZ(k+KHALO) > real(RADAR_ZMAX,kind=RP)) cycle ! Do not draw the stratosphere
    nlev_plot = nlev_plot + 1
  enddo

  allocate(bufr3d(nlev_plot,nlong,nlatg))
  bufr3d(:,:,:) = 0.0

  call rank_1d_2d(myrank_d, proc_i, proc_j)
  ishift = proc_i * nlon
  jshift = proc_j * nlat


  nlev_plot = 0
  do k = plot_zlev_min, plot_zlev_max, plot_zlev_intv

    if (CZ(k+KHALO) > real(RADAR_ZMAX,kind=RP)) cycle ! Do not draw the stratosphere
    nlev_plot = nlev_plot + 1

    bufr3d(nlev_plot,1+ishift:nlon+ishift, 1+jshift:nlat+jshift) = real(ref3d(k,1:nlon,1:nlat), r_sngl)
  enddo

  call MPI_ALLREDUCE(MPI_IN_PLACE, bufr3d, nlong*nlatg*nlev_plot, MPI_REAL, MPI_SUM, MPI_COMM_d, ierr)


  ! Gather required data for reflectivity computation

  nlev_plot = 0
  do k = plot_zlev_min, plot_zlev_max, plot_zlev_intv

    if (CZ(k+KHALO) > real(RADAR_ZMAX,kind=RP)) cycle ! Do not draw the stratosphere
    nlev_plot = nlev_plot + 1

    if ( ( myrank_d + 1 ) == k) then
      write(cheight,'(I5.5)')  int(CZ(k+KHALO)) ! tentative

      plotname = header // "_dbz_"//trim(timelabel) // trim(footer_fcst) // "_z" // cheight // "m"
      call plot_dbz_DCL (bufr3d(nlev_plot,1:nlong,1:nlatg),trim(plotname),cheight,ftsec)
    elseif ( ( myrank_d + 1 ) < k) then
      exit
    end if

  enddo

  deallocate(bufr3d)

  call MPI_BARRIER(MPI_COMM_d, ierr)
  call date_and_time(date=date, time=time)
  if (myrank_d == 0 ) then
    write (6, '(3a,1x,a,1x,a,1x,a)') '[Info:plot] ',header,' finish plotting: ', date, time, trim(timelabel) // trim(footer_fcst)
  endif

  return
end subroutine plot_dafcst_mpi
#endif

! Broadcast restart data among the dacycle (extended) forecast members
subroutine bcast_restart_efcst_mpi()
  use scale_atmos_grid_cartesC_index, only: &
     IA, JA, KA
  use scale_tracer, only: &
     QA
  use mod_atmos_admin, only: &
     ATMOS_do, &
     ATMOS_sw_phy_mp, &
     ATMOS_sw_phy_rd, &
     ATMOS_sw_phy_sf
  use mod_cpl_admin, only: &
     CPL_sw
  use mod_ocean_admin, only: &
     OCEAN_do, &
     OCEAN_ICE_TYPE
  use mod_land_admin, only: &
     LAND_do
  use mod_urban_admin, only: &
     URBAN_do
  use mod_atmos_vars, only: &
     DENS, MOMX, MOMY, MOMZ, &
     RHOT, QTRC, &
     ATMOS_vars_fillhalo
  use mod_atmos_phy_mp_vars, only: &
     ATMOS_PHY_MP_SFLX_rain, &
     ATMOS_PHY_MP_SFLX_snow, &
     ATMOS_PHY_MP_vars_fillhalo
  use mod_atmos_phy_rd_vars, only: &
     ATMOS_PHY_RD_SFLX_LW_up, &
     ATMOS_PHY_RD_SFLX_LW_dn, &
     ATMOS_PHY_RD_SFLX_SW_up, &
     ATMOS_PHY_RD_SFLX_SW_dn, &
     ATMOS_PHY_RD_TOAFLX_LW_up, &
     ATMOS_PHY_RD_TOAFLX_LW_dn, &
     ATMOS_PHY_RD_TOAFLX_SW_up, &
     ATMOS_PHY_RD_TOAFLX_SW_dn, &
     ATMOS_PHY_RD_SFLX_down,    &
     ATMOS_PHY_RD_vars_fillhalo
  use scale_cpl_sfc_index, only: &
     I_R_direct, I_R_diffuse, &
     I_R_IR, I_R_NIR, I_R_VIS
  use mod_atmos_phy_sf_vars, only: &
     ATMOS_PHY_SF_SFC_TEMP, &
     ATMOS_PHY_SF_SFC_albedo, &
     ATMOS_PHY_SF_SFC_Z0M, &
     ATMOS_PHY_SF_SFC_Z0H, &
     ATMOS_PHY_SF_SFC_Z0E, &
     ATMOS_PHY_SF_vars_fillhalo
  use scale_ocean_grid_cartesC_index, only: &
     OKMAX, OIA, OJA, OIS, OIE, OJS, OJE
  use mod_ocean_vars, only : &
     OCEAN_TEMP,     &
     OCEAN_OCN_Z0M,  &
     OCEAN_ICE_TEMP, &
     OCEAN_ICE_MASS, &
     OCEAN_SFC_TEMP, &
     OCEAN_SFC_albedo, &
     OCEAN_SFC_Z0M,    &
     OCEAN_SFC_Z0H,    &
     OCEAN_SFC_Z0E,    &
     OCEAN_ICE_FRAC,   &
     OCEAN_vars_total
  use scale_ocean_phy_ice_simple, only: &
     OCEAN_PHY_ICE_freezetemp, &
     OCEAN_PHY_ICE_fraction
  use mod_land_vars, only : &
     LAND_TEMP, &
     LAND_WATER, &
     LAND_SFC_TEMP, &
     LAND_SFC_albedo
  use scale_land_grid_cartesC_index, only: &
     LKMAX, LIA, LJA
  use mod_urban_vars, only : &
     URBAN_TR, &
     URBAN_TB, &
     URBAN_TG, &
     URBAN_TC, &
     URBAN_QC, &
     URBAN_UC, &
     URBAN_TRL, &
     URBAN_TBL, &
     URBAN_TGL, &
     URBAN_RAINR, &
     URBAN_RAINB, &
     URBAN_RAING, &
     URBAN_ROFF, &
     URBAN_SFC_TEMP, &
     URBAN_SFC_albedo, &
     URBAN_vars_total
  use scale_urban_grid_cartesC_index, only: &
     UIA, UJA, UKS, UKE
  implicit none

  integer :: ierr

  if ( ATMOS_do ) then
    call MPI_BCAST(DENS(1,1,1), KA*IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(MOMX(1,1,1), KA*IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(MOMY(1,1,1), KA*IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(MOMZ(1,1,1), KA*IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(RHOT(1,1,1), KA*IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(QTRC(1,1,1,1), KA*IA*JA*max(QA,1), MPI_RP, 0, MPI_COMM_ef, ierr)

    call ATMOS_vars_fillhalo

    if ( ATMOS_sw_phy_mp ) then
      call MPI_BCAST(ATMOS_PHY_MP_SFLX_rain(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_MP_SFLX_snow(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call ATMOS_PHY_MP_vars_fillhalo
    endif
  
    if ( ATMOS_sw_phy_rd ) then
      call MPI_BCAST(ATMOS_PHY_RD_SFLX_LW_up(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_RD_SFLX_LW_dn(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call MPI_BCAST(ATMOS_PHY_RD_SFLX_SW_up(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_RD_SFLX_SW_dn(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call MPI_BCAST(ATMOS_PHY_RD_TOAFLX_LW_up(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_RD_TOAFLX_LW_dn(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call MPI_BCAST(ATMOS_PHY_RD_TOAFLX_SW_up(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_RD_TOAFLX_SW_dn(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call MPI_BCAST(ATMOS_PHY_RD_SFLX_down(1,1,I_R_direct,I_R_IR), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_RD_SFLX_down(1,1,I_R_diffuse,I_R_IR), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call MPI_BCAST(ATMOS_PHY_RD_SFLX_down(1,1,I_R_direct,I_R_NIR), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_RD_SFLX_down(1,1,I_R_diffuse,I_R_NIR), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call MPI_BCAST(ATMOS_PHY_RD_SFLX_down(1,1,I_R_direct,I_R_VIS), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_RD_SFLX_down(1,1,I_R_diffuse,I_R_VIS), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call ATMOS_PHY_RD_vars_fillhalo
    endif
  
    if ( ATMOS_sw_phy_sf .and. (.not. CPL_sw) ) then
      call MPI_BCAST(ATMOS_PHY_SF_SFC_TEMP(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call MPI_BCAST(ATMOS_PHY_SF_SFC_albedo(1,1,I_R_direct,I_R_IR), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_SF_SFC_albedo(1,1,I_R_diffuse,I_R_IR), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call MPI_BCAST(ATMOS_PHY_SF_SFC_albedo(1,1,I_R_direct,I_R_NIR), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_SF_SFC_albedo(1,1,I_R_diffuse,I_R_NIR), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call MPI_BCAST(ATMOS_PHY_SF_SFC_albedo(1,1,I_R_direct,I_R_VIS), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_SF_SFC_albedo(1,1,I_R_diffuse,I_R_VIS), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call MPI_BCAST(ATMOS_PHY_SF_SFC_Z0M(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_SF_SFC_Z0H(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
      call MPI_BCAST(ATMOS_PHY_SF_SFC_Z0E(1,1), IA*JA, MPI_RP, 0, MPI_COMM_ef, ierr)
  
      call ATMOS_PHY_SF_vars_fillhalo
    endif
  
  !    if ( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_read
  !    if ( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_read
  !    if ( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_read
  !    if ( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_read
  !    if ( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_read

  endif ! ATMOS_do

  if ( OCEAN_do ) then
    call MPI_BCAST(OCEAN_TEMP(1,1,1), OKMAX*OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(OCEAN_OCN_Z0M(1,1), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(OCEAN_ICE_TEMP(1,1), OKMAX*OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(OCEAN_ICE_MASS(1,1), OKMAX*OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(OCEAN_SFC_TEMP(1,1), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(OCEAN_SFC_albedo(1,1,I_R_direct,I_R_IR), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(OCEAN_SFC_albedo(1,1,I_R_diffuse,I_R_IR), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(OCEAN_SFC_albedo(1,1,I_R_direct,I_R_NIR), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(OCEAN_SFC_albedo(1,1,I_R_diffuse,I_R_NIR), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(OCEAN_SFC_albedo(1,1,I_R_direct,I_R_VIS), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(OCEAN_SFC_albedo(1,1,I_R_diffuse,I_R_VIS), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(OCEAN_SFC_Z0M(1,1), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(OCEAN_SFC_Z0H(1,1), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(OCEAN_SFC_Z0E(1,1), OIA*OJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    if ( OCEAN_ICE_TYPE == 'NONE' ) then
       OCEAN_ICE_TEMP(:,:) = OCEAN_PHY_ICE_freezetemp
       OCEAN_ICE_FRAC(:,:) = 0.0_RP
    else
       OCEAN_ICE_TEMP(:,:) = min( OCEAN_ICE_TEMP(:,:), OCEAN_PHY_ICE_freezetemp )
       call OCEAN_PHY_ICE_fraction( OIA, OIS, OIE,       & ! [IN]
                                    OJA, OJS, OJE,       & ! [IN]
                                    OCEAN_ICE_MASS(:,:), & ! [IN]
                                    OCEAN_ICE_FRAC(:,:)  ) ! [OUT]
    endif

    call OCEAN_vars_total

  endif ! OCEAN_do

  if ( LAND_do ) then
    call MPI_BCAST(LAND_TEMP(1,1,1), LKMAX*LIA*LJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(LAND_WATER(1,1,1), LKMAX*LIA*LJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(LAND_SFC_TEMP(1,1), LIA*LJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(LAND_SFC_albedo(1,1,I_R_direct,I_R_IR), LIA*LJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(LAND_SFC_albedo(1,1,I_R_diffuse,I_R_IR), LIA*LJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(LAND_SFC_albedo(1,1,I_R_direct,I_R_NIR), LIA*LJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(LAND_SFC_albedo(1,1,I_R_diffuse,I_R_NIR), LIA*LJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(LAND_SFC_albedo(1,1,I_R_direct,I_R_VIS), LIA*LJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(LAND_SFC_albedo(1,1,I_R_diffuse,I_R_VIS), LIA*LJA, MPI_RP, 0, MPI_COMM_ef, ierr)

  endif ! LAND_do

  if ( URBAN_do ) then
    call MPI_BCAST(URBAN_TR(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_TB(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_TG(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_TC(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_QC(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_UC(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(URBAN_TRL(1,1,1), UIA*UJA*(UKE-UKS+1), MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_TBL(1,1,1), UIA*UJA*(UKE-UKS+1), MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_TGL(1,1,1), UIA*UJA*(UKE-UKS+1), MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(URBAN_RAINR(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_RAINB(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_RAING(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_ROFF(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(URBAN_SFC_TEMP(1,1), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(URBAN_SFC_albedo(1,1,I_R_direct,I_R_IR), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_SFC_albedo(1,1,I_R_diffuse,I_R_IR), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(URBAN_SFC_albedo(1,1,I_R_direct,I_R_NIR), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_SFC_albedo(1,1,I_R_diffuse,I_R_NIR), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call MPI_BCAST(URBAN_SFC_albedo(1,1,I_R_direct,I_R_VIS), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)
    call MPI_BCAST(URBAN_SFC_albedo(1,1,I_R_diffuse,I_R_VIS), UIA*UJA, MPI_RP, 0, MPI_COMM_ef, ierr)

    call URBAN_vars_total

  endif ! URBAN_do

  return
end subroutine bcast_restart_efcst_mpi

subroutine send_recv_emean_others(fcst_cnt)
  use scale_atmos_grid_cartesC_index, only: &
     IA, JA, KA
  use mod_atmos_admin, only: &
     ATMOS_sw_phy_mp, &
     ATMOS_sw_phy_rd, &
     ATMOS_sw_phy_sf
  use mod_atmos_phy_mp_vars, only: &
     ATMOS_PHY_MP_SFLX_rain, &
     ATMOS_PHY_MP_SFLX_snow, &
     ATMOS_PHY_MP_vars_fillhalo
  use mod_atmos_phy_rd_vars, only: &
     ATMOS_PHY_RD_SFLX_LW_up, &
     ATMOS_PHY_RD_SFLX_LW_dn, &
     ATMOS_PHY_RD_SFLX_SW_up, &
     ATMOS_PHY_RD_SFLX_SW_dn, &
     ATMOS_PHY_RD_TOAFLX_LW_up, &
     ATMOS_PHY_RD_TOAFLX_LW_dn, &
     ATMOS_PHY_RD_TOAFLX_SW_up, &
     ATMOS_PHY_RD_TOAFLX_SW_dn, &
     ATMOS_PHY_RD_SFLX_down,    &
     ATMOS_PHY_RD_vars_fillhalo
  use scale_cpl_sfc_index, only: &
     I_R_direct, I_R_diffuse, &
     I_R_IR, I_R_NIR, I_R_VIS
  use mod_atmos_phy_sf_vars, only: &
     ATMOS_PHY_SF_SFC_TEMP, &
     ATMOS_PHY_SF_SFC_albedo, &
     ATMOS_PHY_SF_SFC_Z0M, &
     ATMOS_PHY_SF_SFC_Z0H, &
     ATMOS_PHY_SF_SFC_Z0E, &
     ATMOS_PHY_SF_vars_fillhalo
  use mod_cpl_admin, only: &
     CPL_sw
  use mod_ocean_admin, only: &
     OCEAN_do, &
     OCEAN_ICE_TYPE
  use scale_ocean_grid_cartesC_index, only: &
     OKMAX, OIA, OJA, OIS, OIE, OJS, OJE
  use mod_ocean_vars, only : &
     OCEAN_TEMP,     &
     OCEAN_OCN_Z0M,  &
     OCEAN_ICE_TEMP, &
     OCEAN_ICE_MASS, &
     OCEAN_SFC_TEMP, &
     OCEAN_SFC_albedo, &
     OCEAN_SFC_Z0M,    &
     OCEAN_SFC_Z0H,    &
     OCEAN_SFC_Z0E,    &
     OCEAN_ICE_FRAC,   &
     OCEAN_vars_total
  use scale_ocean_phy_ice_simple, only: &
     OCEAN_PHY_ICE_freezetemp, &
     OCEAN_PHY_ICE_fraction
  use mod_land_admin, only: &
     LAND_do
  use scale_land_grid_cartesC_index, only: &
     LKMAX, LIA, LJA
  use mod_land_vars, only : &
     LAND_TEMP, &
     LAND_WATER, &
     LAND_SFC_TEMP, &
     LAND_SFC_albedo, &
     LAND_vars_total
  use mod_urban_admin, only: &
     URBAN_do
  use scale_urban_grid_cartesC_index, only: &
     UIA, UJA, UKS, UKE
  use mod_urban_vars, only : &
     URBAN_TR, &
     URBAN_TB, &
     URBAN_TG, &
     URBAN_TC, &
     URBAN_QC, &
     URBAN_UC, &
     URBAN_TRL, &
     URBAN_TBL, &
     URBAN_TGL, &
     URBAN_RAINR, &
     URBAN_RAINB, &
     URBAN_RAING, &
     URBAN_ROFF, &
     URBAN_SFC_TEMP, &
     URBAN_SFC_albedo, &
     URBAN_vars_total

  implicit none

  integer, intent(in) :: fcst_cnt

  integer :: srank_a ! send rank_a
  integer :: rrank_a ! receive rank_a

  integer :: atmos_nv2d
  integer :: ocean_nv2d
  integer :: land_nv2d
  integer :: urban_nv2d

  real(RP), allocatable :: atmos2d(:,:,:)
  real(RP), allocatable :: ocean2d(:,:,:)
  real(RP), allocatable :: land2d(:,:,:)
  real(RP), allocatable :: urban2d(:,:,:)

  integer :: ukmax_

  integer :: ierr, tag
  integer, allocatable :: istat(:)

  integer :: iv3d, iv2d

  srank_a = MEMBER * nprocs_d + myrank_d
  rrank_a = nens * nprocs_d + nprocs_d * (fcst_cnt - 1) + myrank_d ! dacycle-forecast member

  if ( myrank_a /= srank_a .and. myrank_a /= rrank_a ) return

  atmos_nv2d = 0
  ocean_nv2d = 0
  land_nv2d = 0
  urban_nv2d = 0

  allocate(istat(MPI_STATUS_SIZE))

  !!! Send/receive variables except for the LETKF control variables
  !!!  To do: Include ATMOS_DYN
  if ( ATMOS_sw_phy_mp ) then

     atmos_nv2d = 2  ! ATMOS_PHY_MP_SFLX_rain
                     ! ATMOS_PHY_MP_SFLX_snow

     allocate(atmos2d(atmos_nv2d,IA,JA))

     tag = myrank_d 
     if ( myrank_a == srank_a ) then
       atmos2d(1,1:IA,1:JA) = ATMOS_PHY_MP_SFLX_rain(1:IA,1:JA)
       atmos2d(2,1:IA,1:JA) = ATMOS_PHY_MP_SFLX_snow(1:IA,1:JA)

       call MPI_Send(atmos2d,atmos_nv2d*IA*JA,MPI_RP,rrank_a,tag,MPI_COMM_a,ierr) 

     elseif ( myrank_a == rrank_a ) then

       call MPI_Recv(atmos2d,atmos_nv2d*IA*JA,MPI_RP,srank_a,tag,MPI_COMM_a,istat,ierr)

       ATMOS_PHY_MP_SFLX_rain(1:IA,1:JA) = atmos2d(1,1:IA,1:JA)
       ATMOS_PHY_MP_SFLX_snow(1:IA,1:JA) = atmos2d(2,1:IA,1:JA)

       call ATMOS_PHY_MP_vars_fillhalo
     endif

     deallocate(atmos2d)

  endif

  if ( ATMOS_sw_phy_rd ) then
    atmos_nv2d = 14 ! ATMOS_PHY_RD_SFLX_LW_up
                    ! ATMOS_PHY_RD_SFLX_LW_dn
                    ! ATMOS_PHY_RD_SFLX_SW_up
                    ! ATMOS_PHY_RD_SFLX_SW_dn
                    ! ATMOS_PHY_RD_TOAFLX_LW_up
                    ! ATMOS_PHY_RD_TOAFLX_LW_dn
                    ! ATMOS_PHY_RD_TOAFLX_SW_up
                    ! ATMOS_PHY_RD_TOAFLX_SW_dn
                    ! ATMOS_PHY_RD_SFLX_down I_R_direct,I_R_IR
                    ! ATMOS_PHY_RD_SFLX_down I_R_diffuse,I_R_IR
                    ! ATMOS_PHY_RD_SFLX_down I_R_direct,I_R_NIR
                    ! ATMOS_PHY_RD_SFLX_down I_R_diffuse,I_R_NIR
                    ! ATMOS_PHY_RD_SFLX_down I_R_direct,I_R_VIS
                    ! ATMOS_PHY_RD_SFLX_down I_R_diffuse,I_R_VIS

     allocate(atmos2d(atmos_nv2d,IA,JA))

     tag = myrank_d 
     if ( myrank_a == srank_a ) then
       atmos2d(1, 1:IA,1:JA) = ATMOS_PHY_RD_SFLX_LW_up(1:IA,1:JA)
       atmos2d(2, 1:IA,1:JA) = ATMOS_PHY_RD_SFLX_LW_dn(1:IA,1:JA)
       atmos2d(3, 1:IA,1:JA) = ATMOS_PHY_RD_SFLX_SW_up(1:IA,1:JA)
       atmos2d(4, 1:IA,1:JA) = ATMOS_PHY_RD_SFLX_SW_dn(1:IA,1:JA)
       atmos2d(5, 1:IA,1:JA) = ATMOS_PHY_RD_TOAFLX_LW_up(1:IA,1:JA)
       atmos2d(6, 1:IA,1:JA) = ATMOS_PHY_RD_TOAFLX_LW_dn(1:IA,1:JA)
       atmos2d(7, 1:IA,1:JA) = ATMOS_PHY_RD_TOAFLX_SW_up(1:IA,1:JA)
       atmos2d(8, 1:IA,1:JA) = ATMOS_PHY_RD_TOAFLX_SW_dn(1:IA,1:JA)
       atmos2d(9, 1:IA,1:JA) = ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_direct,I_R_IR)
       atmos2d(10,1:IA,1:JA) = ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_diffuse,I_R_IR)
       atmos2d(11,1:IA,1:JA) = ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_direct,I_R_NIR)
       atmos2d(12,1:IA,1:JA) = ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_diffuse,I_R_NIR)
       atmos2d(13,1:IA,1:JA) = ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_direct,I_R_VIS)
       atmos2d(14,1:IA,1:JA) = ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_diffuse,I_R_VIS)

       call MPI_Send(atmos2d,atmos_nv2d*IA*JA,MPI_RP,rrank_a,tag,MPI_COMM_a,ierr) 

     elseif ( myrank_a == rrank_a ) then

       call MPI_Recv(atmos2d,atmos_nv2d*IA*JA,MPI_RP,srank_a,tag,MPI_COMM_a,istat,ierr)

       ATMOS_PHY_RD_SFLX_LW_up(1:IA,1:JA)   = atmos2d(1, 1:IA,1:JA)
       ATMOS_PHY_RD_SFLX_LW_dn(1:IA,1:JA)   = atmos2d(2, 1:IA,1:JA)
       ATMOS_PHY_RD_SFLX_SW_up(1:IA,1:JA)   = atmos2d(3, 1:IA,1:JA)
       ATMOS_PHY_RD_SFLX_SW_dn(1:IA,1:JA)   = atmos2d(4, 1:IA,1:JA) 
       ATMOS_PHY_RD_TOAFLX_LW_up(1:IA,1:JA) = atmos2d(5, 1:IA,1:JA) 
       ATMOS_PHY_RD_TOAFLX_LW_dn(1:IA,1:JA) = atmos2d(6, 1:IA,1:JA) 
       ATMOS_PHY_RD_TOAFLX_SW_up(1:IA,1:JA) = atmos2d(7, 1:IA,1:JA) 
       ATMOS_PHY_RD_TOAFLX_SW_dn(1:IA,1:JA) = atmos2d(8, 1:IA,1:JA) 
       ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_direct,I_R_IR)   = atmos2d(9, 1:IA,1:JA) 
       ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_diffuse,I_R_IR)  = atmos2d(10,1:IA,1:JA)
       ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_direct,I_R_NIR)  = atmos2d(11,1:IA,1:JA) 
       ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_diffuse,I_R_NIR) = atmos2d(12,1:IA,1:JA) 
       ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_direct,I_R_VIS)  = atmos2d(13,1:IA,1:JA) 
       ATMOS_PHY_RD_SFLX_down(1:IA,1:JA,I_R_diffuse,I_R_VIS) = atmos2d(14,1:IA,1:JA) 

       call ATMOS_PHY_RD_vars_fillhalo

     endif

     deallocate(atmos2d)

  endif

  if ( ATMOS_sw_phy_sf .and. (.not. CPL_sw) ) then

    atmos_nv2d = 10 ! ATMOS_PHY_SF_SFC_TEMP
                    ! ATMOS_PHY_SF_SFC_albedo I_R_direct,I_R_IR
                    ! ATMOS_PHY_SF_SFC_albedo I_R_diffuse,I_R_IR
                    ! ATMOS_PHY_SF_SFC_albedo I_R_direct,I_R_NIR
                    ! ATMOS_PHY_SF_SFC_albedo I_R_diffuse,I_R_NIR
                    ! ATMOS_PHY_SF_SFC_albedo I_R_direct,I_R_VIS
                    ! ATMOS_PHY_SF_SFC_albedo I_R_diffuse,I_R_VIS
                    ! ATMOS_PHY_SF_SFC_Z0M
                    ! ATMOS_PHY_SF_SFC_Z0H
                    ! ATMOS_PHY_SF_SFC_Z0E

     allocate(atmos2d(atmos_nv2d,IA,JA))

     tag = myrank_d 
     if ( myrank_a == srank_a ) then
       atmos2d(1, 1:IA,1:JA) = ATMOS_PHY_SF_SFC_TEMP(1:IA,1:JA)
       atmos2d(2, 1:IA,1:JA) = ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_direct,I_R_IR)
       atmos2d(3, 1:IA,1:JA) = ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_diffuse,I_R_IR)
       atmos2d(4, 1:IA,1:JA) = ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_direct,I_R_NIR)
       atmos2d(5, 1:IA,1:JA) = ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_diffuse,I_R_NIR)
       atmos2d(6, 1:IA,1:JA) = ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_direct,I_R_VIS)
       atmos2d(7, 1:IA,1:JA) = ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_diffuse,I_R_VIS)
       atmos2d(8, 1:IA,1:JA) = ATMOS_PHY_SF_SFC_Z0M(1:IA,1:JA)
       atmos2d(9, 1:IA,1:JA) = ATMOS_PHY_SF_SFC_Z0H(1:IA,1:JA)
       atmos2d(10,1:IA,1:JA) = ATMOS_PHY_SF_SFC_Z0E(1:IA,1:JA)

       call MPI_Send(atmos2d,atmos_nv2d*IA*JA,MPI_RP,rrank_a,tag,MPI_COMM_a,ierr) 

     elseif ( myrank_a == rrank_a ) then

       call MPI_Recv(atmos2d,atmos_nv2d*IA*JA,MPI_RP,srank_a,tag,MPI_COMM_a,istat,ierr)

       ATMOS_PHY_SF_SFC_TEMP(1:IA,1:JA) = atmos2d(1,1:IA,1:JA) 
       ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_direct,I_R_IR)   = atmos2d(2,1:IA,1:JA) 
       ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_diffuse,I_R_IR)  = atmos2d(3,1:IA,1:JA)
       ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_direct,I_R_NIR)  = atmos2d(4,1:IA,1:JA)
       ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_diffuse,I_R_NIR) = atmos2d(5,1:IA,1:JA)
       ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_direct,I_R_VIS)  = atmos2d(6,1:IA,1:JA) 
       ATMOS_PHY_SF_SFC_albedo(1:IA,1:JA,I_R_diffuse,I_R_VIS) = atmos2d(7,1:IA,1:JA) 
       ATMOS_PHY_SF_SFC_Z0M(1:IA,1:JA) = atmos2d(8,1:IA,1:JA) 
       ATMOS_PHY_SF_SFC_Z0H(1:IA,1:JA) = atmos2d(9,1:IA,1:JA) 
       ATMOS_PHY_SF_SFC_Z0E(1:IA,1:JA) = atmos2d(10,1:IA,1:JA) 

       call ATMOS_PHY_SF_vars_fillhalo

     endif

     deallocate(atmos2d)

  endif

  if ( OCEAN_do ) then
    ocean_nv2d = OKMAX ! OCEAN_TEMP (3D)

    ocean_nv2d = ocean_nv2d +  13 ! OCEAN_OCN_Z0M
                                  ! OCEAN_ICE_TEMP
                                  ! OCEAN_ICE_MASS
                                  ! OCEAN_SFC_TEMP 
                                  ! OCEAN_SFC_albedo I_R_direct,I_R_IR
                                  ! OCEAN_SFC_albedo I_R_diffuse,I_R_IR
                                  ! OCEAN_SFC_albedo I_R_direct,I_R_NIR
                                  ! OCEAN_SFC_albedo I_R_diffuse,I_R_NIR
                                  ! OCEAN_SFC_albedo I_R_direct,I_R_VIS
                                  ! OCEAN_SFC_albedo I_R_diffuse,I_R_VIS
                                  ! OCEAN_SFC_Z0M
                                  ! OCEAN_SFC_Z0H
                                  ! OCEAN_SFC_Z0E

     allocate(ocean2d(ocean_nv2d,OIA,OJA))

     tag = myrank_d 
     if ( myrank_a == srank_a ) then
       ocean2d(1:OKMAX,1:OIA,1:OJA) = OCEAN_TEMP(1:OKMAX,1:OIA,1:OJA)

       ocean2d(OKMAX+1, 1:OIA,1:OJA) = OCEAN_OCN_Z0M(1:OIA,1:OJA)
       ocean2d(OKMAX+2, 1:OIA,1:OJA) = OCEAN_ICE_TEMP(1:OIA,1:OJA)
       ocean2d(OKMAX+3, 1:OIA,1:OJA) = OCEAN_ICE_MASS(1:OIA,1:OJA)
       ocean2d(OKMAX+4, 1:OIA,1:OJA) = OCEAN_SFC_TEMP(1:OIA,1:OJA)
       ocean2d(OKMAX+5, 1:OIA,1:OJA) = OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_direct,I_R_IR)
       ocean2d(OKMAX+6, 1:OIA,1:OJA) = OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_diffuse,I_R_IR)
       ocean2d(OKMAX+7, 1:OIA,1:OJA) = OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_direct,I_R_NIR)
       ocean2d(OKMAX+8, 1:OIA,1:OJA) = OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_diffuse,I_R_NIR)
       ocean2d(OKMAX+9, 1:OIA,1:OJA) = OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_direct,I_R_VIS)
       ocean2d(OKMAX+10,1:OIA,1:OJA) = OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_diffuse,I_R_VIS)
       ocean2d(OKMAX+11,1:OIA,1:OJA) = OCEAN_SFC_Z0M(1:OIA,1:OJA)
       ocean2d(OKMAX+12,1:OIA,1:OJA) = OCEAN_SFC_Z0H(1:OIA,1:OJA)
       ocean2d(OKMAX+13,1:OIA,1:OJA) = OCEAN_SFC_Z0E(1:OIA,1:OJA)

       call MPI_Send(ocean2d,ocean_nv2d*OIA*OJA,MPI_RP,rrank_a,tag,MPI_COMM_a,ierr) 

     elseif ( myrank_a == rrank_a ) then

       call MPI_Recv(ocean2d,ocean_nv2d*OIA*OJA,MPI_RP,srank_a,tag,MPI_COMM_a,istat,ierr)

       OCEAN_TEMP(1:OKMAX,1:OIA,1:OJA) = ocean2d(1:OKMAX,1:OIA,1:OJA) 

       OCEAN_OCN_Z0M(1:OIA,1:OJA)  = ocean2d(OKMAX+1,1:OIA,1:OJA) 
       OCEAN_ICE_TEMP(1:OIA,1:OJA) = ocean2d(OKMAX+2,1:OIA,1:OJA)  
       OCEAN_ICE_MASS(1:OIA,1:OJA) = ocean2d(OKMAX+3,1:OIA,1:OJA)  
       OCEAN_SFC_TEMP(1:OIA,1:OJA) = ocean2d(OKMAX+4,1:OIA,1:OJA) 
       OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_direct,I_R_IR)   = ocean2d(OKMAX+5, 1:OIA,1:OJA)  
       OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_diffuse,I_R_IR)  = ocean2d(OKMAX+6, 1:OIA,1:OJA)  
       OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_direct,I_R_NIR)  = ocean2d(OKMAX+7, 1:OIA,1:OJA)  
       OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_diffuse,I_R_NIR) = ocean2d(OKMAX+8, 1:OIA,1:OJA)  
       OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_direct,I_R_VIS)  = ocean2d(OKMAX+9, 1:OIA,1:OJA) 
       OCEAN_SFC_albedo(1:OIA,1:OJA,I_R_diffuse,I_R_VIS) = ocean2d(OKMAX+10,1:OIA,1:OJA)  
       OCEAN_SFC_Z0M(1:OIA,1:OJA) = ocean2d(OKMAX+11,1:OIA,1:OJA) 
       OCEAN_SFC_Z0H(1:OIA,1:OJA) = ocean2d(OKMAX+12,1:OIA,1:OJA)  
       OCEAN_SFC_Z0E(1:OIA,1:OJA) = ocean2d(OKMAX+13,1:OIA,1:OJA)  

       if ( OCEAN_ICE_TYPE == 'NONE' ) then
          OCEAN_ICE_TEMP(:,:) = OCEAN_PHY_ICE_freezetemp
          OCEAN_ICE_FRAC(:,:) = 0.0_RP
       else
          OCEAN_ICE_TEMP(:,:) = min( OCEAN_ICE_TEMP(:,:), OCEAN_PHY_ICE_freezetemp )
          call OCEAN_PHY_ICE_fraction( OIA, OIS, OIE,       & ! [IN]
                                       OJA, OJS, OJE,       & ! [IN]
                                       OCEAN_ICE_MASS(:,:), & ! [IN]
                                       OCEAN_ICE_FRAC(:,:)  ) ! [OUT]
          call OCEAN_vars_total
 
       endif
     endif

     deallocate(ocean2d)

  endif

  if ( LAND_do ) then
    land_nv2d = LKMAX*2 ! LAND_TEMP, LAND_WATER (3D)

    land_nv2d = land_nv2d +  7  ! LAND_SFC_TEMP
                                ! LAND_SFC_albedo I_R_direct,I_R_IR
                                ! LAND_SFC_albedo I_R_diffuse,I_R_IR
                                ! LAND_SFC_albedo I_R_direct,I_R_NIR
                                ! LAND_SFC_albedo I_R_diffuse,I_R_NIR
                                ! LAND_SFC_albedo I_R_direct,I_R_VIS
                                ! LAND_SFC_albedo I_R_diffuse,I_R_VIS

     allocate(land2d(land_nv2d,LIA,LJA))

     tag = myrank_d 
     if ( myrank_a == srank_a ) then

       land2d(1:LKMAX,        1:LIA,1:LJA) = LAND_TEMP (1:LKMAX,1:LIA,1:LJA)
       land2d(1+LKMAX:LKMAX*2,1:LIA,1:LJA) = LAND_WATER(1:LKMAX,1:LIA,1:LJA)

       land2d(LKMAX*2+1,      1:LIA,1:LJA) = LAND_SFC_TEMP  (1:LIA,1:LJA)
       land2d(LKMAX*2+2,      1:LIA,1:LJA) = LAND_SFC_albedo(1:LIA,1:LJA,I_R_direct,I_R_IR)
       land2d(LKMAX*2+3,      1:LIA,1:LJA) = LAND_SFC_albedo(1:LIA,1:LJA,I_R_diffuse,I_R_IR)
       land2d(LKMAX*2+4,      1:LIA,1:LJA) = LAND_SFC_albedo(1:LIA,1:LJA,I_R_direct,I_R_NIR)
       land2d(LKMAX*2+5,      1:LIA,1:LJA) = LAND_SFC_albedo(1:LIA,1:LJA,I_R_diffuse,I_R_NIR)
       land2d(LKMAX*2+6,      1:LIA,1:LJA) = LAND_SFC_albedo(1:LIA,1:LJA,I_R_direct,I_R_VIS)
       land2d(LKMAX*2+7,      1:LIA,1:LJA) = LAND_SFC_albedo(1:LIA,1:LJA,I_R_diffuse,I_R_VIS)

       call MPI_Send(land2d,land_nv2d*LIA*LJA,MPI_RP,rrank_a,tag,MPI_COMM_a,ierr) 

     elseif ( myrank_a == rrank_a ) then

       call MPI_Recv(land2d,land_nv2d*LIA*LJA,MPI_RP,srank_a,tag,MPI_COMM_a,istat,ierr)

       LAND_TEMP (1:LKMAX,1:LIA,1:LJA) = land2d(1:LKMAX,        1:LIA,1:LJA)  
       LAND_WATER(1:LKMAX,1:LIA,1:LJA) = land2d(1+LKMAX:LKMAX*2,1:LIA,1:LJA)  

       LAND_SFC_TEMP(1:LIA,1:LJA)                       = land2d(LKMAX*2+1,1:LIA,1:LJA)  
       LAND_SFC_albedo(1:LIA,1:LJA,I_R_direct,I_R_IR)   = land2d(LKMAX*2+2,1:LIA,1:LJA)  
       LAND_SFC_albedo(1:LIA,1:LJA,I_R_diffuse,I_R_IR)  = land2d(LKMAX*2+3,1:LIA,1:LJA)  
       LAND_SFC_albedo(1:LIA,1:LJA,I_R_direct,I_R_NIR)  = land2d(LKMAX*2+4,1:LIA,1:LJA)  
       LAND_SFC_albedo(1:LIA,1:LJA,I_R_diffuse,I_R_NIR) = land2d(LKMAX*2+5,1:LIA,1:LJA)  
       LAND_SFC_albedo(1:LIA,1:LJA,I_R_direct,I_R_VIS)  = land2d(LKMAX*2+6,1:LIA,1:LJA)
       LAND_SFC_albedo(1:LIA,1:LJA,I_R_diffuse,I_R_VIS) = land2d(LKMAX*2+7,1:LIA,1:LJA)  

       call LAND_vars_total
 
     endif

     deallocate(land2d)

  endif ! LAND_do

  if ( URBAN_do ) then

    ukmax_ = UKE - UKS + 1

    urban_nv2d = ukmax_*3   ! URBAN_TRL, URBAN_TBL, URBAN_TGL

    urban_nv2d = urban_nv2d +  17 ! URBAN_TR
                                  ! URBAN_TB
                                  ! URBAN_TG
                                  ! URBAN_TC
                                  ! URBAN_QC
                                  ! URBAN_UC
                                  ! URBAN_RAINR
                                  ! URBAN_RAINB
                                  ! URBAN_RAING
                                  ! URBAN_ROFF
                                  ! URBAN_SFC_TEMP
                                  ! URBAN_SFC_albedo I_R_direct,I_R_IR
                                  ! URBAN_SFC_albedo I_R_diffuse,I_R_IR
                                  ! URBAN_SFC_albedo I_R_direct,I_R_NIR
                                  ! URBAN_SFC_albedo I_R_diffuse,I_R_NIR
                                  ! URBAN_SFC_albedo I_R_direct,I_R_VIS
                                  ! URBAN_SFC_albedo I_R_diffuse,I_R_VIS


     allocate(urban2d(urban_nv2d,UIA,UJA))

     tag = myrank_d 
     if ( myrank_a == srank_a ) then
       urban2d(1:ukmax_,           1:UIA,1:UJA) = URBAN_TRL(1:ukmax_,1:UIA,1:UJA)
       urban2d(1+ukmax_:ukmax_*2,  1:UIA,1:UJA) = URBAN_TBL(1:ukmax_,1:UIA,1:UJA)
       urban2d(1+ukmax_*2:ukmax_*3,1:UIA,1:UJA) = URBAN_TGL(1:ukmax_,1:UIA,1:UJA)

       urban2d(ukmax_*3+1, 1:UIA,1:UJA) = URBAN_TR(1:UIA,1:UJA)
       urban2d(ukmax_*3+2, 1:UIA,1:UJA) = URBAN_TB(1:UIA,1:UJA)
       urban2d(ukmax_*3+3, 1:UIA,1:UJA) = URBAN_TG(1:UIA,1:UJA)
       urban2d(ukmax_*3+4, 1:UIA,1:UJA) = URBAN_TC(1:UIA,1:UJA)
       urban2d(ukmax_*3+5, 1:UIA,1:UJA) = URBAN_QC(1:UIA,1:UJA)
       urban2d(ukmax_*3+6, 1:UIA,1:UJA) = URBAN_UC(1:UIA,1:UJA)
       urban2d(ukmax_*3+7, 1:UIA,1:UJA) = URBAN_RAINR(1:UIA,1:UJA)
       urban2d(ukmax_*3+8, 1:UIA,1:UJA) = URBAN_RAINB(1:UIA,1:UJA)
       urban2d(ukmax_*3+9, 1:UIA,1:UJA) = URBAN_RAING(1:UIA,1:UJA)
       urban2d(ukmax_*3+10,1:UIA,1:UJA) = URBAN_ROFF(1:UIA,1:UJA)
       urban2d(ukmax_*3+11,1:UIA,1:UJA) = URBAN_SFC_TEMP(1:UIA,1:UJA)
       urban2d(ukmax_*3+12,1:UIA,1:UJA) = URBAN_SFC_albedo(1:UIA,1:UJA,I_R_direct,I_R_IR)
       urban2d(ukmax_*3+13,1:UIA,1:UJA) = URBAN_SFC_albedo(1:UIA,1:UJA,I_R_diffuse,I_R_IR)
       urban2d(ukmax_*3+14,1:UIA,1:UJA) = URBAN_SFC_albedo(1:UIA,1:UJA,I_R_direct,I_R_NIR)
       urban2d(ukmax_*3+15,1:UIA,1:UJA) = URBAN_SFC_albedo(1:UIA,1:UJA,I_R_diffuse,I_R_NIR)
       urban2d(ukmax_*3+16,1:UIA,1:UJA) = URBAN_SFC_albedo(1:UIA,1:UJA,I_R_direct,I_R_VIS)
       urban2d(ukmax_*3+17,1:UIA,1:UJA) = URBAN_SFC_albedo(1:UIA,1:UJA,I_R_diffuse,I_R_VIS)

       call MPI_Send(urban2d,urban_nv2d*UIA*UJA,MPI_RP,rrank_a,tag,MPI_COMM_a,ierr) 

     elseif ( myrank_a == rrank_a ) then

       call MPI_Recv(urban2d,urban_nv2d*UIA*UJA,MPI_RP,srank_a,tag,MPI_COMM_a,istat,ierr)

       URBAN_TRL(1:ukmax_,1:UIA,1:UJA) = urban2d(1:ukmax_,           1:UIA,1:UJA)  
       URBAN_TBL(1:ukmax_,1:UIA,1:UJA) = urban2d(1+ukmax_:ukmax_*2,  1:UIA,1:UJA)  
       URBAN_TGL(1:ukmax_,1:UIA,1:UJA) = urban2d(1+ukmax_*2:ukmax_*3,1:UIA,1:UJA)  
      
       URBAN_TR(1:UIA,1:UJA)       = urban2d(ukmax_*3+1, 1:UIA,1:UJA)  
       URBAN_TB(1:UIA,1:UJA)       = urban2d(ukmax_*3+2, 1:UIA,1:UJA)  
       URBAN_TG(1:UIA,1:UJA)       = urban2d(ukmax_*3+3, 1:UIA,1:UJA)  
       URBAN_TC(1:UIA,1:UJA)       = urban2d(ukmax_*3+4, 1:UIA,1:UJA)  
       URBAN_QC(1:UIA,1:UJA)       = urban2d(ukmax_*3+5, 1:UIA,1:UJA)  
       URBAN_UC(1:UIA,1:UJA)       = urban2d(ukmax_*3+6, 1:UIA,1:UJA)  
       URBAN_RAINR(1:UIA,1:UJA)    = urban2d(ukmax_*3+7, 1:UIA,1:UJA)  
       URBAN_RAINB(1:UIA,1:UJA)    = urban2d(ukmax_*3+8, 1:UIA,1:UJA)  
       URBAN_RAING(1:UIA,1:UJA)    = urban2d(ukmax_*3+9, 1:UIA,1:UJA)  
       URBAN_ROFF(1:UIA,1:UJA)     = urban2d(ukmax_*3+10,1:UIA,1:UJA)  
       URBAN_SFC_TEMP(1:UIA,1:UJA) = urban2d(ukmax_*3+11,1:UIA,1:UJA)  

       URBAN_SFC_albedo(1:UIA,1:UJA,I_R_direct,I_R_IR)   = urban2d(ukmax_*3+12,1:UIA,1:UJA)  
       URBAN_SFC_albedo(1:UIA,1:UJA,I_R_diffuse,I_R_IR)  = urban2d(ukmax_*3+13,1:UIA,1:UJA)  
       URBAN_SFC_albedo(1:UIA,1:UJA,I_R_direct,I_R_NIR)  = urban2d(ukmax_*3+14,1:UIA,1:UJA)  
       URBAN_SFC_albedo(1:UIA,1:UJA,I_R_diffuse,I_R_NIR) = urban2d(ukmax_*3+15,1:UIA,1:UJA)  
       URBAN_SFC_albedo(1:UIA,1:UJA,I_R_direct,I_R_VIS)  = urban2d(ukmax_*3+16,1:UIA,1:UJA)  
       URBAN_SFC_albedo(1:UIA,1:UJA,I_R_diffuse,I_R_VIS) = urban2d(ukmax_*3+17,1:UIA,1:UJA)  

       call URBAN_vars_total

     endif
     deallocate(urban2d)

  endif ! URBAN_do

  deallocate(istat)

  return
end subroutine send_recv_emean_others

subroutine pawr_toshiba_hd_mpi(lon0, lat0, z0, missing, range_res, na, nr, ne, &
                            AZDIM, ELDIM, n_type, RDIM, az, el, rtdat, utime_obs)
  implicit none

  real(r_size), intent(inout) :: lon0, lat0, z0, missing
  integer, intent(inout) :: range_res, na, nr, ne
  integer, intent(in) :: AZDIM, ELDIM, n_type, RDIM
  real(r_sngl), intent(inout) :: az(AZDIM, ELDIM, n_type), el(AZDIM, ELDIM, n_type)
  real(r_sngl), intent(inout) :: rtdat(RDIM, AZDIM, ELDIM, n_type)
  integer, intent(inout) :: utime_obs(6)

  integer :: ierr

  if ( nprocs_o < 2 ) return

  call MPI_BCAST(lon0, 1, MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(lat0, 1, MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(z0, 1, MPI_r_size, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(missing, 1, MPI_r_size, 0, MPI_COMM_o, ierr)

  call MPI_BCAST(range_res, 1, MPI_INTEGER, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(na, 1, MPI_INTEGER, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(nr, 1, MPI_INTEGER, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(ne, 1, MPI_INTEGER, 0, MPI_COMM_o, ierr)

  call MPI_BCAST(az,    AZDIM*ELDIM*n_type, MPI_REAL, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(el,    AZDIM*ELDIM*n_type, MPI_REAL, 0, MPI_COMM_o, ierr)
  call MPI_BCAST(rtdat, RDIM*AZDIM*ELDIM*n_type, MPI_REAL, 0, MPI_COMM_o, ierr)

  call MPI_BCAST(utime_obs, 6, MPI_INTEGER, 0, MPI_COMM_o, ierr)

  return
end subroutine pawr_toshiba_hd_mpi

!subroutine pawr_toshiba_scattv_mpi(RDIM, AZDIM, ELDIM, n_type, ne_lmax, rtdat, az, el, rtdat_l, az_l, el_l)
subroutine pawr_toshiba_scattv_mpi(RDIM, AZDIM, ELDIM, n_type, ne_lmax, rtdat, rtdat_l)
  implicit none
  integer, intent(in) :: AZDIM, ELDIM, n_type, RDIM, ne_lmax
  real(r_sngl), intent(in) :: rtdat(RDIM, AZDIM, ELDIM, n_type)
!  real(r_sngl), intent(in) :: az(AZDIM, ELDIM, n_type)
!  real(r_sngl), intent(in) :: el(AZDIM, ELDIM, n_type)

  real(r_sngl), intent(out) :: rtdat_l(RDIM, AZDIM, ne_lmax, n_type)
!  real(r_sngl), intent(out) :: az_l(AZDIM, ne_lmax, n_type)
!  real(r_sngl), intent(out) :: el_l(AZDIM, ne_lmax, n_type)

  real(r_sngl), allocatable :: sbuf1(:), rbuf1(:)
!  real(r_sngl), allocatable :: sbuf2(:), rbuf2(:)
!  real(r_sngl), allocatable :: sbuf3(:), rbuf3(:)
  integer :: ierr

  integer :: it, ie, ia
  integer :: cnt

  allocate(sbuf1(ne_lmax*nprocs_o*RDIM*AZDIM))
  allocate(rbuf1(ne_lmax*RDIM*AZDIM))
!  allocate(sbuf2(ne_lmax*nprocs_o*RDIM*AZDIM))
!  allocate(rbuf2(ne_lmax*RDIM*AZDIM))
!  allocate(sbuf3(ne_lmax*nprocs_o*RDIM*AZDIM))
!  allocate(rbuf3(ne_lmax*RDIM*AZDIM))

  do it = 1, n_type

    if ( myrank_o == 0 ) then
      do ie = 1, ELDIM
        do ia = 1, AZDIM
          cnt = RDIM * AZDIM * (ie - 1) + RDIM *  (ia - 1)
          sbuf1(1+cnt:RDIM+cnt) = rtdat(1:RDIM,ia,ie,it)
!          sbuf2(1+cnt:RDIM+cnt) = az(1:RDIM,ia,ie)
!          sbuf3(1+cnt:RDIM+cnt) = el(1:RDIM,ia,ie)
        enddo
      enddo
    endif

    call MPI_Scatter(sbuf1, ne_lmax*RDIM*AZDIM, MPI_REAL, rbuf1, ne_lmax*RDIM*AZDIM, &
                     MPI_REAL, 0, MPI_COMM_o, ierr)
!    call MPI_Scatter(sbuf2, ne_lmax*RDIM*AZDIM, MPI_REAL, rbuf2, ne_lmax*RDIM*AZDIM, &
!                     MPI_REAL, 0, MPI_COMM_o, ierr)
!    call MPI_Scatter(sbuf3, ne_lmax*RDIM*AZDIM, MPI_REAL, rbuf3, ne_lmax*RDIM*AZDIM, &
!                     MPI_REAL, 0, MPI_COMM_o, ierr)

    do ie = 1, ne_lmax
      do ia = 1, AZDIM
        cnt = RDIM * AZDIM * (ie - 1) + RDIM *  (ia - 1)
        rtdat_l(1:RDIM,ia,ie,it) = rbuf1(1+cnt:RDIM+cnt)
!        az_l(1:RDIM,ia,ie) = rbuf2(1+cnt:RDIM+cnt)
!        el_l(1:RDIM,ia,ie) = rbuf3(1+cnt:RDIM+cnt)
      enddo
    enddo
  enddo ! nt

  deallocate(sbuf1,rbuf1)
!  deallocate(sbuf2,rbuf2)
!  deallocate(sbuf3,rbuf3)

  return
end subroutine pawr_toshiba_scattv_mpi

subroutine pawr_i8_allreduce(array,size)
  implicit none

  integer(8),intent(in) :: size
  integer(8),intent(inout) :: array(size)
  real(r_dble) :: rarray(size)
  
  integer :: ierr

  if (nprocs_o < 2) return
  rarray = real(array,kind=r_dble)
  call MPI_ALLREDUCE(MPI_IN_PLACE, rarray, size, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_o, ierr)

  array = rarray

  return
end subroutine pawr_i8_allreduce

!subroutine pawr_3dvar_allreduce(na,nr,ne,radlon)
!  implicit none
!
!  integer,intent(in) :: na, nr, ne
!  real(r_size),intent(inout) :: radlon(na, nr, ne)
!
!  real(r_size) :: rbuf3d(na, nr, ne)
!
!  integer :: ierr
!
!  if (nprocs_o < 2) return
!  rbuf3d = radlon
!  call MPI_ALLREDUCE(MPI_IN_PLACE, rbuf3d, na*nr*ne, MPI_r_size, MPI_SUM, MPI_COMM_o, ierr)
!  radlon = rbuf3d 
!
!  return
!end subroutine pawr_3dvar_allreduce

!SUBROUTINE get_nobs_mpi(obsfile,nrec,nn)
!SUBROUTINE read_obs2_mpi(obsfile,nn,nbv,elem,rlon,rlat,rlev,odat,oerr,otyp,tdif,hdxf,iqc)
!SUBROUTINE allreduce_obs_mpi(n,nbv,hdxf,iqc)

!===============================================================================
END MODULE common_mpi_scale
