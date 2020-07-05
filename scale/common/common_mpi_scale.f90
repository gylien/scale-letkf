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
  use common_obs_scale

  use scale_precision, only: RP
  use scale_comm_cartesC, only: COMM_datatype
#ifdef PNETCDF
  use scale_file, only: FILE_AGGREGATE
#endif

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
  integer,save :: MPI_COMM_d, nprocs_d, myrank_d
  integer,save :: MPI_COMM_e, nprocs_e, myrank_e

  integer, parameter :: max_timer_levels = 5
  integer, parameter :: timer_name_width = 50
  real(r_dble), private, parameter :: timer_neglect = 1.0d-3
  real(r_dble), private, save :: timer_save(max_timer_levels) = -9.0d10

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

!  write(6,'(A,I6.6,A,I6.6)') 'Hello from MYRANK ', myrank, '/', nprocs-1
  if (r_size == r_dble) then
    MPI_r_size = MPI_DOUBLE_PRECISION
  else if (r_size == r_sngl) then
    MPI_r_size = MPI_REAL
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
      CX => ATMOS_GRID_CARTESC_CX, &
      CY => ATMOS_GRID_CARTESC_CY, &
      DX, DY
  use scale_atmos_grid_cartesC_index, only: &
      IHALO, JHALO
  use scale_mapprojection, only: &
      MAPPROJECTION_xy2lonlat
  implicit none
  integer :: color, key
  integer :: ierr
  character(len=filelenmax) :: filename
  real(r_size), allocatable :: height3dtmp(:,:,:)
  real(r_size), allocatable :: lon2dtmp(:,:)
  real(r_size), allocatable :: lat2dtmp(:,:)
  real(RP), allocatable :: lon2dtmp_RP(:,:)
  real(RP), allocatable :: lat2dtmp_RP(:,:)
  real(RP), allocatable :: height3dtmp_RP(:,:,:)
  real(RP) :: lon_RP, lat_RP
  integer :: i, j
  real(r_size) :: ri, rj

  call mpi_timer('', 2)

  ! Communicator for 1-iteration ensemble member groups
  !-----------------------------------------------------------------------------

  color = myrank_to_pe
  key = myrank_to_mem(1) - 1

  call MPI_COMM_SPLIT(MPI_COMM_a, color, key, MPI_COMM_e, ierr)

  call MPI_COMM_SIZE(MPI_COMM_e, nprocs_e, ierr)
  call MPI_COMM_RANK(MPI_COMM_e, myrank_e, ierr)

#ifdef DEBUG
  if (nprocs_e /= n_mem*n_mempn) then
    write (6, '(A)'), '[Error] XXXXXX wrong!!'
    stop
  end if
#endif

  call mpi_timer('set_common_mpi_scale:mpi_comm_split_e:', 2)

  ! Read/calculate model coordinates
  !-----------------------------------------------------------------------------

  if (VERIFY_COORD) then
    if (myrank_e == 0) then
!      allocate (height3d(nlev,nlon,nlat))
      allocate (lon2d(nlon,nlat))
      allocate (lat2d(nlon,nlat))
      allocate (height3dtmp(nlev,nlon,nlat))
      allocate (lon2dtmp(nlon,nlat))
      allocate (lat2dtmp(nlon,nlat))
      allocate (lon2dtmp_RP(nlon,nlat))
      allocate (lat2dtmp_RP(nlon,nlat))
      allocate (height3dtmp_RP(nlev,nlon,nlat))

      if (.not. allocated(topo2d)) then
        allocate (topo2d(nlon,nlat))
        call read_topo(LETKF_TOPO_IN_BASENAME, topo2d)
      end if
!      call scale_calc_z(topo2d, height3d)

!$OMP PARALLEL DO PRIVATE(i,j,ri,rj) COLLAPSE(2)
      do j = 1, nlat
        do i = 1, nlon
          ri = real(i + IHALO, r_size)
          rj = real(j + JHALO, r_size)
          call MAPPROJECTION_xy2lonlat( real(ri-1.0_r_size, kind=RP)*DX + CX(1), &
                                        real(rj-1.0_r_size, kind=RP)*DY + CY(1), &
                                        lon_RP, lat_RP )
!                                        lon2d(i,j), lat2d(i,j))

          lon2d(i,j) = real(lon_RP, kind=r_size)*rad2deg
          lat2d(i,j) = real(lat_RP, kind=r_size)*rad2deg
        end do
      end do
!$OMP END PARALLEL DO

      filename = GUES_IN_BASENAME
      call filename_replace_mem(filename, myrank_to_mem(1))
      call read_restart_coor( filename, lon2dtmp_RP, lat2dtmp_RP, height3dtmp_RP )
     
      lon2dtmp = real(lon2dtmp_RP, kind=r_size)
      lat2dtmp = real(lat2dtmp_RP, kind=r_size)
      height3dtmp = real(height3dtmp_RP, kind=r_size)

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
  write (6,'(A,I6.6,A,I7)') 'MYRANK ', myrank, ' number of grid points: nij1 =', nij1

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
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call read_topo_par(LETKF_TOPO_IN_BASENAME, topo2dtmp, MPI_COMM_d)
      else
#endif
        call read_topo(LETKF_TOPO_IN_BASENAME, topo2dtmp)
#ifdef PNETCDF
      end if
#endif
      if (maxval(abs(topo2dtmp - topo2d)) > 1.0d-6) then
        write (6, '(A,F15.7)') '[Error] topo height in history files and restart files are inconsistent; maxdiff = ', maxval(abs(topo2dtmp - topo2d))
        stop
      end if
#endif
    else
      allocate (topo2d(nlon,nlat))
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call read_topo_par(LETKF_TOPO_IN_BASENAME, topo2d, MPI_COMM_d)
      else
#endif
        call read_topo(LETKF_TOPO_IN_BASENAME, topo2d)
#ifdef PNETCDF
      end if
#endif
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

  if (LOG_LEVEL >= 1) then
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

  ! settings related to mean (only valid when mem >= MEMBER+1)
  !----------------------------------------------------------------
  if (mem >= MEMBER+1) then
    nens = mem
    mmean = MEMBER+1

    mmean_rank_e = mod(mmean-1, n_mem*n_mempn)
#ifdef DEBUG
    if (mmean_rank_e /= rank_to_mem(1,mempe_to_rank(1,mmean)+1)-1) then
      write (6, '(A)'), '[Error] XXXXXX wrong!!'
      stop
    end if
#endif

    msprd_rank_e = mmean_rank_e

    if (DET_RUN) then
      nensobs = MEMBER+1
      mmdetobs = MEMBER+1
    else
      nensobs = MEMBER
    end if
  end if

  ! settings related to mdet (only valid when mem >= MEMBER+2)
  !----------------------------------------------------------------
  if (mem >= MEMBER+2 .and. DET_RUN) then
    mmdet = MEMBER+2
    if (DET_RUN_CYCLED) then
      mmdetin = mmdet
    else
      mmdetin = mmean
    end if

    mmdet_rank_e = mod(mmdet-1, n_mem*n_mempn)
#ifdef DEBUG
    if (mmdet_rank_e /= rank_to_mem(1,mempe_to_rank(1,mmdet)+1)-1) then
      write (6, '(A)'), '[Error] XXXXXX wrong!!'
      stop
    end if
#endif
  end if

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
  use scale_atmos_grid_cartesC, only: &
    ATMOS_GRID_CARTESC_setup, &
    ATMOS_GRID_CARTESC_DOMAIN_CENTER_X, &
    ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y
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
      filename = HISTORY_IN_BASENAME
      call filename_replace_mem(filename, im)
    else if (im == mmean) then
      filename = HISTORY_MEAN_IN_BASENAME
    else if (im == mmdet) then
      filename = HISTORY_MDET_IN_BASENAME
    end if

#ifdef PNETCDF
    if (FILE_AGGREGATE) then
      call read_history_par(trim(filename), step, v3dg, v2dg, MPI_COMM_d)
    else
#endif
      call read_history(trim(filename), step, v3dg, v2dg)
#ifdef PNETCDF
    end if
#endif
  end if

  return
end subroutine read_ens_history_iter

!-------------------------------------------------------------------------------
! Read ensemble first guess data and distribute to processes
!-------------------------------------------------------------------------------
subroutine read_ens_mpi(v3d, v2d)
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
        filename = GUES_IN_BASENAME
        call filename_replace_mem(filename, im)
      else if (im == mmean) then
        filename = GUES_MEAN_INOUT_BASENAME
      else if (im == mmdet) then
        filename = GUES_MDET_IN_BASENAME
      end if

!      write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',myrank_d,'.nc'
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call read_restart_par(filename, v3dg, v2dg, MPI_COMM_d)
      else
#endif
        call read_restart(filename, v3dg, v2dg)
#ifdef PNETCDF
      end if
#endif

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
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call read_restart_par(filename, v3dg, v2dg, MPI_COMM_d)
      else
#endif
      call read_restart(filename, v3dg, v2dg)
#ifdef PNETCDF
      end if
#endif
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
subroutine write_ens_mpi(v3d, v2d, monit_step)
  implicit none
  real(r_size), intent(in) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(in) :: v2d(nij1,nens,nv2d)
  integer, intent(in), optional :: monit_step
  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  character(len=filelenmax) :: filename
  integer :: it, im, mstart, mend
  integer :: monit_step_

  monit_step_ = 0
  if (present(monit_step)) then
    monit_step_ = monit_step
  end if

  do it = 1, nitmax
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    im = myrank_to_mem(it)

    mstart = 1 + (it-1)*nprocs_e
    mend = min(it*nprocs_e, nens)
    if (mstart <= mend) then
      call gather_grd_mpi_alltoall(mstart, mend, v3d, v2d, v3dg, v2dg)
    end if

    call mpi_timer('write_ens_mpi:gather_grd_mpi_alltoall:', 2)

    if (monit_step_ > 0 .and. mstart <= mmean .and. mmean <= mend) then
      call monit_obs_mpi(v3dg, v2dg, monit_step_)

      call mpi_timer('write_ens_mpi:monit_obs_mpi:', 2)
    end if

    ! Note: write all members + mean + mdet
    ! 
    if ((im >= 1 .and. im <= MEMBER) .or. im == mmean .or. im == mmdet) then
      if (im <= MEMBER) then
        filename = ANAL_OUT_BASENAME
        call filename_replace_mem(filename, im)
      else if (im == mmean) then
        filename = ANAL_MEAN_OUT_BASENAME
      else if (im == mmdet) then
        filename = ANAL_MDET_OUT_BASENAME
      end if

!      write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',myrank_d,'.nc'
      call state_trans_inv(v3dg)

      call mpi_timer('write_ens_mpi:state_trans_inv:', 2)

#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call write_restart_par(filename, v3dg, v2dg, MPI_COMM_d)
      else
#endif
        call write_restart(filename, v3dg, v2dg)
#ifdef PNETCDF
      end if
#endif

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
! MPI driver for monitoring observation departure statistics
!-------------------------------------------------------------------------------
subroutine monit_obs_mpi(v3dg, v2dg, monit_step)
  implicit none
  real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP), intent(in) :: v2dg(nlon,nlat,nv2d)
  integer, intent(in) :: monit_step

  integer :: nobs(nid_obs)
  integer :: nobs_g(nid_obs)
  real(r_size) :: bias(nid_obs)
  real(r_size) :: bias_g(nid_obs)
  real(r_size) :: rmse(nid_obs)
  real(r_size) :: rmse_g(nid_obs)
  logical :: monit_type(nid_obs)
  integer :: obsdep_g_nobs
  integer, allocatable :: obsdep_g_set(:)
  integer, allocatable :: obsdep_g_idx(:)
  integer, allocatable :: obsdep_g_qc(:)
  real(r_size), allocatable :: obsdep_g_omb(:)
  real(r_size), allocatable :: obsdep_g_oma(:)
  integer :: cnts
  integer :: cntr(nprocs_d)
  integer :: dspr(nprocs_d)
  integer :: i, ip, ierr

  call mpi_timer('', 2)

  ! NOTE: need to use 'mmean_rank_e' processes to run this calculation
  !       because only these processes have read topo files in 'topo2d'
  ! 
  if (myrank_e == mmean_rank_e) then
    call monit_obs(v3dg, v2dg, topo2d, nobs, bias, rmse, monit_type, .true., monit_step)

    call mpi_timer('monit_obs_mpi:monit_obs:', 2)

    do i = 1, nid_obs
      if (monit_type(i)) then
        nobs_g(i) = nobs(i)
        if (nobs(i) == 0) then
          bias_g(i) = 0.0d0
          rmse_g(i) = 0.0d0
        else
          bias_g(i) = bias(i) * real(nobs(i), r_size)
          rmse_g(i) = rmse(i) * rmse(i) * real(nobs(i), r_size)
        end if
      end if
    end do

    if (nprocs_d > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, nobs_g, nid_obs, MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, bias_g, nid_obs, MPI_r_size,  MPI_SUM, MPI_COMM_d, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, rmse_g, nid_obs, MPI_r_size,  MPI_SUM, MPI_COMM_d, ierr)
    end if

    do i = 1, nid_obs
      if (monit_type(i)) then
        if (nobs_g(i) == 0) then
          bias_g(i) = undef
          rmse_g(i) = undef
        else
          bias_g(i) = bias_g(i) / REAL(nobs_g(i),r_size)
          rmse_g(i) = sqrt(rmse_g(i) / REAL(nobs_g(i),r_size))
        end if
      else
        nobs_g(i) = -1
        bias_g(i) = undef
        rmse_g(i) = undef
      end if
    end do

    call mpi_timer('monit_obs_mpi:stat:mpi_allreduce(domain):', 2)

    if (OBSDEP_OUT .and. monit_step == 2) then
      cnts = obsdep_nobs
      cntr = 0
      cntr(myrank_d+1) = cnts
      call MPI_ALLREDUCE(MPI_IN_PLACE, cntr, nprocs_d, MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
      dspr = 0
      do ip = 1, nprocs_d-1
        dspr(ip+1) = dspr(ip) + cntr(ip)
      end do

      obsdep_g_nobs = dspr(nprocs_d) + cntr(nprocs_d)
      allocate (obsdep_g_set(obsdep_g_nobs))
      allocate (obsdep_g_idx(obsdep_g_nobs))
      allocate (obsdep_g_qc (obsdep_g_nobs))
      allocate (obsdep_g_omb(obsdep_g_nobs))
      allocate (obsdep_g_oma(obsdep_g_nobs))

      if (obsdep_g_nobs > 0) then
        call MPI_GATHERV(obsdep_set, cnts, MPI_INTEGER, obsdep_g_set, cntr, dspr, MPI_INTEGER, 0, MPI_COMM_d, ierr)
        call MPI_GATHERV(obsdep_idx, cnts, MPI_INTEGER, obsdep_g_idx, cntr, dspr, MPI_INTEGER, 0, MPI_COMM_d, ierr)
        call MPI_GATHERV(obsdep_qc,  cnts, MPI_INTEGER, obsdep_g_qc,  cntr, dspr, MPI_INTEGER, 0, MPI_COMM_d, ierr)
        call MPI_GATHERV(obsdep_omb, cnts, MPI_r_size,  obsdep_g_omb, cntr, dspr, MPI_r_size,  0, MPI_COMM_d, ierr)
        call MPI_GATHERV(obsdep_oma, cnts, MPI_r_size,  obsdep_g_oma, cntr, dspr, MPI_r_size,  0, MPI_COMM_d, ierr)
      end if

      if (myrank_d == 0) then
        write (6,'(A,I6.6,2A)') 'MYRANK ', myrank,' is writing an obsda file ', trim(OBSDEP_OUT_BASENAME)//'.dat'
        call write_obs_dep(trim(OBSDEP_OUT_BASENAME)//'.dat', &
                           obsdep_g_nobs, obsdep_g_set, obsdep_g_idx, obsdep_g_qc, obsdep_g_omb, obsdep_g_oma)
      end if
      deallocate (obsdep_g_set)
      deallocate (obsdep_g_idx)
      deallocate (obsdep_g_qc )
      deallocate (obsdep_g_omb)
      deallocate (obsdep_g_oma)

      call mpi_timer('monit_obs_mpi:obsdep:mpi_allreduce(domain):', 2)
    end if ! [ OBSDEP_OUT .and. monit_step == 2 ]

    if (monit_step == 2) then
      deallocate (obsdep_set)
      deallocate (obsdep_idx)
      deallocate (obsdep_qc )
      deallocate (obsdep_omb)
      deallocate (obsdep_oma)
    end if
  end if ! [ myrank_e == mmean_rank_e ]

  if (DEPARTURE_STAT_ALL_PROCESSES) then
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    call MPI_BCAST(nobs,       nid_obs, MPI_INTEGER, mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(bias,       nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(rmse,       nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(nobs_g,     nid_obs, MPI_INTEGER, mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(bias_g,     nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(rmse_g,     nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(monit_type, nid_obs, MPI_LOGICAL, mmean_rank_e, MPI_COMM_e, ierr)

    call mpi_timer('monit_obs_mpi:mpi_allreduce(ens):', 2)
  end if

  if (DEPARTURE_STAT_ALL_PROCESSES .or. myrank_e == mmean_rank_e) then
    if (monit_step == 1) then
      write(6,'(2A)') 'OBSERVATIONAL DEPARTURE STATISTICS [GUESS] (IN THIS SUBDOMAIN):'
    else if (monit_step == 2) then
      write(6,'(2A)') 'OBSERVATIONAL DEPARTURE STATISTICS [ANALYSIS] (IN THIS SUBDOMAIN):'
    end if
    call monit_print(nobs, bias, rmse, monit_type)

    if (monit_step == 1) then
      write(6,'(2A)') 'OBSERVATIONAL DEPARTURE STATISTICS [GUESS] (GLOBAL):'
    else if (monit_step == 2) then
      write(6,'(2A)') 'OBSERVATIONAL DEPARTURE STATISTICS [ANALYSIS] (GLOBAL):'
    end if
    call monit_print(nobs_g, bias_g, rmse_g, monit_type)

    call mpi_timer('monit_obs_mpi:monit_print:', 2)
  end if

  return
end subroutine monit_obs_mpi

!-------------------------------------------------------------------------------
! Gather ensemble mean to {mmean_rank_e} and write SCALE restart files
!-------------------------------------------------------------------------------
subroutine write_ensmean(filename, v3d, v2d, calced, monit_step)
  implicit none
  character(len=*), intent(in) :: filename
  real(r_size), intent(inout) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(inout) :: v2d(nij1,nens,nv2d)
  logical, intent(in), optional :: calced
  integer, intent(in), optional :: monit_step

  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  logical :: calced_
  integer :: monit_step_

  call mpi_timer('', 2)

  calced_ = .false.
  if (present(calced)) then
    calced_ = calced
  end if
  monit_step_ = 0
  if (present(monit_step)) then
    monit_step_ = monit_step
  end if

  if (.not. calced) then
    call ensmean_grd(MEMBER, nens, nij1, v3d, v2d)

    call mpi_timer('write_ensmean:ensmean_grd:', 2)
  end if

  call mpi_timer('', 2, barrier=MPI_COMM_e)

  call gather_grd_mpi(mmean_rank_e, v3d(:,:,mmean,:), v2d(:,mmean,:), v3dg, v2dg)

  call mpi_timer('write_ensmean:gather_grd_mpi:', 2)

  if (monit_step_ > 0) then
    call monit_obs_mpi(v3dg, v2dg, monit_step_)

    call mpi_timer('write_ensmean:monit_obs_mpi:', 2)
  end if

  if (myrank_e == mmean_rank_e) then
    call state_trans_inv(v3dg)

    call mpi_timer('write_ensmean:state_trans_inv:', 2)

#ifdef PNETCDF
    if (FILE_AGGREGATE) then
      call write_restart_par(filename, v3dg, v2dg, MPI_COMM_d)
    else
#endif
      call write_restart(filename, v3dg, v2dg)
#ifdef PNETCDF
    end if
#endif

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
#ifdef PNETCDF
    if (FILE_AGGREGATE) then
      call write_restart_par(filename, v3dg, v2dg, MPI_COMM_d)
    else
#endif
      call write_restart(filename, v3dg, v2dg) !!
#ifdef PNETCDF
    end if
#endif

    call mpi_timer('write_enssprd:write_restart:', 2)
  end if

  return
end subroutine write_enssprd

!-------------------------------------------------------------------------------
! Read all observation data from files of various formats
!-------------------------------------------------------------------------------
subroutine read_obs_all_mpi(obs)
  implicit none
  type(obs_info), intent(out) :: obs(OBS_IN_NUM)
  integer :: iof, ierr

  call mpi_timer('', 2)

  if (myrank_a == 0) then
    call read_obs_all(obs)

    call mpi_timer('read_obs_all_mpi:read_obs_all:', 2)
  end if

  call mpi_timer('', 2, barrier=MPI_COMM_a)

  do iof = 1, OBS_IN_NUM
    call MPI_BCAST(obs(iof)%nobs, 1, MPI_INTEGER, 0, MPI_COMM_a, ierr)
    if (myrank_a /= 0) then
      call obs_info_allocate(obs(iof), extended=.true.)
    end if

    call MPI_BCAST(obs(iof)%elm, obs(iof)%nobs, MPI_INTEGER, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%lon, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%lat, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%lev, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%dat, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%err, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%typ, obs(iof)%nobs, MPI_INTEGER, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%dif, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%meta, max_obs_info_meta, MPI_r_size, 0, MPI_COMM_a, ierr)
  end do ! [ iof = 1, OBS_IN_NUM ]

  call mpi_timer('read_obs_all_mpi:mpi_bcast:', 2)

  return
end subroutine read_obs_all_mpi

!-------------------------------------------------------------------------------
! Read all observation data from files of various formats
!-------------------------------------------------------------------------------
subroutine get_nobs_da_mpi(nobs)
  implicit none
  integer, intent(out) :: nobs
  character(filelenmax) :: obsdafile
  character(11) :: obsda_suffix = '.000000.dat'
  integer :: ierr

! read from all available data by every processes
!-----------------------------
!  if ((myrank_to_mem(1) >= 1 .and. myrank_to_mem(1) <= MEMBER) .or. &
!      myrank_to_mem(1) == mmdetin) then
!    if (myrank_to_mem(1) <= MEMBER) then
!      obsdafile = OBSDA_IN_BASENAME
!      call filename_replace_mem(obsdafile, myrank_to_mem(1))
!    else if (myrank_to_mem(1) == mmean) then
!      obsdafile = OBSDA_MEAN_IN_BASENAME
!    else if (myrank_to_mem(1) == mmdet) then
!      obsdafile = OBSDA_MDET_IN_BASENAME
!    end if
!    write (obsda_suffix(2:7), '(I6.6)') myrank_d
!#ifdef H08
!    call get_nobs(trim(obsdafile) // obsda_suffix, 6, nobs) ! H08
!#else
!    call get_nobs(trim(obsdafile) // obsda_suffix, 4, nobs)
!#endif
!  end if

! read by process 0 and broadcast
!-----------------------------
  if (myrank_e == 0) then
    obsdafile = OBSDA_IN_BASENAME
    call filename_replace_mem(obsdafile, 1)
    write (obsda_suffix(2:7), '(I6.6)') myrank_d
#ifdef H08
    call get_nobs(trim(obsdafile) // obsda_suffix, 6, nobs) ! H08
#else
    call get_nobs(trim(obsdafile) // obsda_suffix, 4, nobs)
#endif
  end if
  call MPI_BCAST(nobs, 1, MPI_INTEGER, 0, MPI_COMM_e, ierr)
!-----------------------------

  return
end subroutine get_nobs_da_mpi

!-------------------------------------------------------------------------------
! Partially reduce observations processed in the same processes in the iteration
!-------------------------------------------------------------------------------
#ifdef H08
subroutine obs_da_value_partial_reduce_iter(obsda, iter, nstart, nobs, ensval, qc, lev, val2)
#else
subroutine obs_da_value_partial_reduce_iter(obsda, iter, nstart, nobs, ensval, qc)
#endif
  implicit none
  type(obs_da_value), intent(inout) :: obsda
  integer, intent(in)      :: iter
  integer, intent(in)      :: nstart
  integer, intent(in)      :: nobs
  real(r_size), intent(in) :: ensval(nobs)
  integer, intent(in)      :: qc(nobs)
#ifdef H08
  real(r_size), intent(in) :: lev(nobs)
  real(r_size), intent(in) :: val2(nobs)
#endif
  integer :: nend
  integer :: im

  if (nobs <= 0) then
    return
  end if
  nend = nstart + nobs - 1
  im = myrank_to_mem(iter)
  if (.not. ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin)) then
    return
  end if

  ! variables with an ensemble dimension
  obsda%ensval(iter,nstart:nend) = ensval

  ! variables without an ensemble dimension
  obsda%qc(nstart:nend) = max(obsda%qc(nstart:nend), qc)
#ifdef H08
  if (im <= MEMBER) then ! only consider lev, val2 from members, not from the means
    obsda%lev(nstart:nend) = obsda%lev(nstart:nend) + lev
    obsda%val2(nstart:nend) = obsda%val2(nstart:nend) + val2
  end if
#endif

  return
end subroutine obs_da_value_partial_reduce_iter

!-------------------------------------------------------------------------------
! Allreduce observations to obtain complete ensemble values
!-------------------------------------------------------------------------------
subroutine obs_da_value_allreduce(obsda)
  implicit none
  type(obs_da_value), intent(inout) :: obsda
  real(r_size), allocatable :: ensval_bufs(:,:)
  real(r_size), allocatable :: ensval_bufr(:,:)
  integer :: cnts
  integer :: cntr(nprocs_e)
  integer :: dspr(nprocs_e)
  integer :: current_shape(2)
  integer :: ie, it, im, imb, ierr

  if (obsda%nobs <= 0) then
    return
  end if

  call mpi_timer('', 3)

  ! variables with an ensemble dimension
  cntr(:) = 0
  do ie = 1, nprocs_e
    do it = 1, nitmax
      im = ranke_to_mem(it, ie)
      if ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin) then
        cntr(ie) = cntr(ie) + 1
      end if
    end do
  end do
  allocate (ensval_bufs(obsda%nobs, cntr(myrank_e+1)))
  allocate (ensval_bufr(obsda%nobs, nensobs))

  do im = 1, cntr(myrank_e+1)
    ensval_bufs(:,im) = obsda%ensval(im,:)
  end do

  cntr(:) = cntr(:) * obsda%nobs
  cnts = cntr(myrank_e+1)
  dspr(1) = 0
  do ie = 2, nprocs_e
    dspr(ie) = dspr(ie-1) + cntr(ie-1)
  end do

  call mpi_timer('obs_da_value_allreduce:copy_bufs:', 3, barrier=MPI_COMM_e)

  call MPI_ALLGATHERV(ensval_bufs, cnts, MPI_r_size, ensval_bufr, cntr, dspr, MPI_r_size, MPI_COMM_e, ierr)

  call mpi_timer('obs_da_value_allreduce:mpi_allgatherv:', 3)

  current_shape = shape(obsda%ensval)
  if (current_shape(1) < nensobs) then
    deallocate (obsda%ensval)
    allocate (obsda%ensval(nensobs, obsda%nobs))
  end if

  imb = 0
  do ie = 1, nprocs_e
    do it = 1, nitmax
      im = ranke_to_mem(it, ie)
      if ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin) then
        imb = imb + 1
        if (im == mmdetin) then
          obsda%ensval(mmdetobs,:) = ensval_bufr(:,imb)
        else
          obsda%ensval(im,:) = ensval_bufr(:,imb)
        end if
      end if
    end do
  end do
  deallocate(ensval_bufs, ensval_bufr)

  call mpi_timer('obs_da_value_allreduce:copy_bufr:', 3, barrier=MPI_COMM_e)

  ! variables without an ensemble dimension
  if (nprocs_e > 1) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, obsda%qc(:), obsda%nobs, MPI_INTEGER, MPI_MAX, MPI_COMM_e, ierr)
  end if
#ifdef H08
  if (nprocs_e > 1) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, obsda%lev(:), obsda%nobs, MPI_r_size, MPI_SUM, MPI_COMM_e, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, obsda%val2(:), obsda%nobs, MPI_r_size, MPI_SUM, MPI_COMM_e, ierr)
  end if
  obsda%lev(:) = obsda%lev(:) / real(MEMBER, r_size)
  obsda%val2(:) = obsda%val2(:) / real(MEMBER, r_size)
#endif

  call mpi_timer('obs_da_value_allreduce:mpi_allreduce:', 3)

  return
end subroutine obs_da_value_allreduce

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

    if (i == level) then
      timer_save(i) = timer_after_barrier
    else
      timer_save(i) = -9.0d10 ! reset the timer for all levels under this level
    end if

  end do

  return
end subroutine mpi_timer

!SUBROUTINE get_nobs_mpi(obsfile,nrec,nn)
!SUBROUTINE read_obs2_mpi(obsfile,nn,nbv,elem,rlon,rlat,rlev,odat,oerr,otyp,tdif,hdxf,iqc)
!SUBROUTINE allreduce_obs_mpi(n,nbv,hdxf,iqc)

!===============================================================================
END MODULE common_mpi_scale
