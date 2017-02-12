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
  use scale_comm, only: COMM_datatype

  implicit none
  public

  integer,save :: nij1
  integer,save :: nij1max
  integer,allocatable,save :: nij1node(:)




!  integer,save          :: universal_comm                         ! universal communicator
!  integer,save          :: universal_nprocs                       ! number of procs in universal communicator
!  logical,save          :: universal_master                       ! master process  in universal communicator?
!  logical,save          :: universal_myrank                       ! 

!  integer,save          :: global_comm                            ! communicator for each member

!  integer,save          :: local_comm                             ! assigned local communicator
!  integer,save          :: intercomm_parent                       ! inter communicator with parent
!  integer,save          :: intercomm_child                        ! inter communicator with child

!  integer,save          :: NUM_DOMAIN                   = 1       ! number of domains
!  integer,save          :: PRC_DOMAINS(PRC_DOMAIN_nlim) = 0       ! number of total process in each domain
!  character(len=H_LONG),save :: CONF_FILES (PRC_DOMAIN_nlim) = "" ! name of configulation files
!  logical,save          :: ABORT_ALL_JOBS               = .false. ! abort all jobs or not?
!  logical,save          :: LOG_SPLIT                    = .false. ! log-output for mpi splitting?



!!!!  real(r_size),allocatable,save :: phi1(:)
!!!  real(r_size),allocatable,save :: lon1(:),lat1(:)
!!!  real(r_size),allocatable,save :: lonu1(:),latu1(:)
!!!  real(r_size),allocatable,save :: lonv1(:),latv1(:)
!!!  real(r_size),allocatable,save :: ri1(:),rj1(:)
!!!!  real(r_size),allocatable,save :: wg1(:)
  real(r_size),allocatable,save :: topo2d(:,:)
  real(r_size),allocatable,save :: rig1(:),rjg1(:)
  real(r_size),allocatable,save :: topo1(:)
  real(r_size),allocatable,save :: hgt1(:,:)




  integer,save :: nitmax ! maximum number of model files processed by a process
  integer,allocatable,save :: procs(:)
  integer,allocatable,save :: mem2node(:,:)
  integer,allocatable,save :: mem2proc(:,:)
  integer,allocatable,save :: proc2mem(:,:,:)
  integer,save :: n_mem
  integer,save :: n_mempn

  integer,save :: ens_mygroup = -1
  integer,save :: ens_myrank = -1
  logical,save :: myrank_use = .false.

  integer,save :: nens
  integer,save :: nensobs

  integer,save :: mmean
  integer,save :: mmdet
  integer,save :: mmdetin
  integer,save :: mmdetobs

  integer,save :: mmean_rank_e
  integer,save :: mmdet_rank_e
  integer,save :: msprd_rank_e



  integer,save :: MPI_COMM_e, nprocs_e, myrank_e
  integer,save :: MPI_COMM_d, nprocs_d, myrank_d
  integer,save :: MPI_COMM_a, nprocs_a, myrank_a



!  character(9) scale_filename = 'file.0000'

contains

!-------------------------------------------------------------------------------
! initialize_mpi_scale
!-------------------------------------------------------------------------------
subroutine initialize_mpi_scale
  use scale_process, only: &
     PRC_MPIstart, &
     PRC_UNIVERSAL_setup, &
     PRC_UNIVERSAL_myrank
  implicit none
  integer :: universal_comm   ! dummy
  integer :: universal_nprocs ! dummy
  logical :: universal_master ! dummy
  integer :: ierr

  call PRC_MPIstart( universal_comm ) ! [OUT]

!  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
!  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call PRC_UNIVERSAL_setup( universal_comm,   & ! [IN]
                            universal_nprocs, & ! [OUT]
                            universal_master  ) ! [OUT]
  nprocs = universal_nprocs
  myrank = PRC_UNIVERSAL_myrank

  write(6,'(A,I6.6,A,I6.6)') 'Hello from MYRANK ', myrank, '/', nprocs-1
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
!  use scale_process, only: PRC_MPIfinish
  implicit none
  integer :: ierr

!  call PRC_MPIfinish
  call MPI_Finalize(ierr)

  return
end subroutine finalize_mpi_scale

!-------------------------------------------------------------------------------
! set_common_mpi_scale
!-------------------------------------------------------------------------------
SUBROUTINE set_common_mpi_scale
  use scale_grid_index, only: &
    IHALO, &
    JHALO
!  use scale_process, only: &
!    PRC_myrank, &

  implicit none
  INTEGER :: i,n
  INTEGER :: ierr !,buf(4)
!  LOGICAL :: ex

  integer :: MPI_G_WORLD, MPI_G
!  integer :: iproc,jproc
  integer,allocatable :: ranks(:)
  integer,allocatable :: ranks_a(:)

  integer :: ip


  real(r_dble) :: rrtimer00, rrtimer
  rrtimer00 = MPI_WTIME()


  WRITE(6,'(A)') 'Hello from set_common_mpi_scale'

  nprocs_e = n_mem*n_mempn
  nprocs_a = nprocs_e*MEM_NP

  allocate (ranks(nprocs_e))
  allocate (ranks_a(nprocs_a))

  call MPI_Comm_group(MPI_COMM_WORLD,MPI_G_WORLD,ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_scale:mpi_comm_group_world:   ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  do ip = 1, nprocs
    if (proc2mem(1,1,ip) >= 1) then
      if (proc2mem(2,1,ip) == proc2mem(2,1,myrank+1)) then
        ranks(proc2mem(1,1,ip)) = ip-1
      end if
      ranks_a((proc2mem(1,1,ip)-1)*MEM_NP+proc2mem(2,1,ip)+1) = ip-1
    end if
  end do

!write(6,'(A,7I6)') '######===', myrank, ranks(:)

!!!!!! rewrite using MPI_COMM_SPLIT ??? !!!!!!

  call MPI_Group_incl(MPI_G_WORLD,nprocs_e,ranks,MPI_G,ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_scale:mpi_group_incl_e:       ', rrtimer-rrtimer00
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  rrtimer00 = MPI_WTIME()


  call MPI_Comm_create(MPI_COMM_WORLD,MPI_G,MPI_COMM_e,ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_scale:mpi_comm_create_e:      ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  call MPI_Comm_size(MPI_COMM_e,nprocs_e,ierr)
  call MPI_Comm_rank(MPI_COMM_e,myrank_e,ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_scale:mpi_comm_size_rank_e:   ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


!--

  call MPI_Group_incl(MPI_G_WORLD,nprocs_e*MEM_NP,ranks_a,MPI_G,ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_scale:mpi_group_incl_a:       ', rrtimer-rrtimer00
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  rrtimer00 = MPI_WTIME()


  call MPI_Comm_create(MPI_COMM_WORLD,MPI_G,MPI_COMM_a,ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_scale:mpi_comm_create_a:      ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  call MPI_Comm_size(MPI_COMM_a,nprocs_a,ierr)
  call MPI_Comm_rank(MPI_COMM_a,myrank_a,ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_scale:mpi_comm_size_rank_a:   ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


!--

  call MPI_Comm_size(MPI_COMM_d,nprocs_d,ierr)
  call MPI_Comm_rank(MPI_COMM_d,myrank_d,ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_scale:mpi_comm_size_rank_d:   ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


!write(6,'(A,9I6)') '######===', myrank, myrank_e, nprocs_e, ranks(:)

  deallocate(ranks)

  i = MOD(nlon*nlat,nprocs_e)
  nij1max = (nlon*nlat - i)/nprocs_e + 1
  IF(myrank_e < i) THEN
    nij1 = nij1max
  ELSE
    nij1 = nij1max - 1
  END IF
  WRITE(6,'(A,I6.6,A,I6)') 'MYRANK ',myrank,' number of grid points: nij1= ',nij1
  ALLOCATE(nij1node(nprocs_e))
  DO n=1,nprocs_e
    IF(n-1 < i) THEN
      nij1node(n) = nij1max
    ELSE
      nij1node(n) = nij1max - 1
    END IF
  END DO


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_scale:nij1_cal:               ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  RETURN
END SUBROUTINE set_common_mpi_scale

!-------------------------------------------------------------------------------
! unset_common_mpi_scale
!-------------------------------------------------------------------------------
SUBROUTINE unset_common_mpi_scale
  implicit none
  integer:: ierr

  call MPI_Comm_free(MPI_COMM_e,ierr)
  call MPI_Comm_free(MPI_COMM_a,ierr)
!  call unset_scalelib

  RETURN
END SUBROUTINE unset_common_mpi_scale

!-------------------------------------------------------------------------------
! set_common_mpi_grid
!-------------------------------------------------------------------------------
subroutine set_common_mpi_grid
  use scale_grid_index, only: &
    IHALO, &
    JHALO
  use scale_process, only: &
    PRC_myrank

  implicit none
  REAL(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:)
  INTEGER :: i,j
  integer :: iproc, jproc


  integer :: ierr
  real(r_dble) :: rrtimer00, rrtimer
  rrtimer00 = MPI_WTIME()


  ALLOCATE(topo2d(nlon,nlat))

  ALLOCATE(rig1(nij1))
  ALLOCATE(rjg1(nij1))
  ALLOCATE(topo1(nij1))

  ALLOCATE(hgt1(nij1,nlev))

  ALLOCATE(v3d(nij1,nlev,nv3d))
  ALLOCATE(v2d(nij1,nv2d))

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


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_grid:rij_cal:                 ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  ! Note that only 'mmean_rank_e' processes have topo data in 'topo2d', 
  ! that will be used later in the obs departure monitor calculation
    call read_topo(LETKF_TOPO_IN_BASENAME, topo2d)
    v3dg(1,:,:,3) = topo2d


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_grid:read_topo:               ', rrtimer-rrtimer00


  end if


  call MPI_BARRIER(MPI_COMM_e, ierr)
  rrtimer00 = MPI_WTIME()


  CALL scatter_grd_mpi(mmean_rank_e,v3dg,v2dg,v3d,v2d)

  rig1   = v3d(:,1,1)
  rjg1   = v3d(:,1,2)
  topo1  = v3d(:,1,3)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_grid:scatter:                 ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  call scale_calc_z(nij1, topo1, hgt1)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_common_mpi_grid:scale_calc_z:            ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  return
end subroutine set_common_mpi_grid

!-------------------------------------------------------------------------------
! set_mem_node_proc
!-------------------------------------------------------------------------------
SUBROUTINE set_mem_node_proc(mem)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: mem
  INTEGER :: tppn,tppnt,tmod
  INTEGER :: n,ns,nn,m,q,qs,i,j,it,ip


!  integer :: ierr
  real(r_dble) :: rrtimer00, rrtimer
  rrtimer00 = MPI_WTIME()


  ALLOCATE(procs(nprocs))
  ns = 0
  DO n = 1, NNODES
    procs(ns+1:ns+PPN) = n
    ns = ns + PPN
  END DO

  IF(MEM_NODES > 1) THEN
    n_mem = NNODES / MEM_NODES
    n_mempn = 1
  ELSE
    n_mem = NNODES
    n_mempn = PPN / MEM_NP
  END IF
  nitmax = (mem - 1) / (n_mem * n_mempn) + 1
  tppn = MEM_NP / MEM_NODES
  tmod = MOD(MEM_NP, MEM_NODES)

  ALLOCATE(mem2node(MEM_NP,mem))
  ALLOCATE(mem2proc(MEM_NP,mem))
  ALLOCATE(proc2mem(2,nitmax,nprocs))
  proc2mem = -1
  m = 1
mem_loop: DO it = 1, nitmax
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
            ip = (n+nn)*PPN + i*MEM_NP + q
            if (m <= mem) then
              mem2node(qs+1,m) = n+nn
              mem2proc(qs+1,m) = ip
            end if
            proc2mem(1,it,ip+1) = m    ! These lines are outside of (m <= mem) condition
            proc2mem(2,it,ip+1) = qs   ! in order to cover over the entire first iteration
            qs = qs + 1
          END DO
        END DO
        m = m + 1
        n = n + MEM_NODES
      END DO
    END DO
  END DO mem_loop

  ens_mygroup = proc2mem(1,1,myrank+1)
  ens_myrank = proc2mem(2,1,myrank+1)
  if (ens_mygroup >= 1) then
    myrank_use = .true.
  end if

  ! settings related to mean, mdet (only valid when mem = MEMBER+2)
  !----------------------------------------------------------------
  if (mem == MEMBER+2) then
    nens = MEMBER+2
    mmean = MEMBER+1
    mmdet = MEMBER+2
    if (DET_RUN_CYCLED) then
      mmdetin = mmdet
    else
      mmdetin = mmean
    end if

    mmean_rank_e = mod(mmean-1, n_mem*n_mempn)
    mmdet_rank_e = mod(mmdet-1, n_mem*n_mempn)
    msprd_rank_e = mmean_rank_e !!!!!! may be changed to mmdet_rank_e for (potential) imporved performance !!!!!!
#ifdef DEBUG
    if (mmean_rank_e /= proc2mem(1,1,mem2proc(1,mmean)+1)-1) then
      write (6, '(A)'), '[Error] XXXXXX wrong!!'
      stop
    end if
    if (mmdet_rank_e /= proc2mem(1,1,mem2proc(1,mmdet)+1)-1) then
      write (6, '(A)'), '[Error] XXXXXX wrong!!'
      stop
    end if
#endif

    nensobs = MEMBER+1
    mmdetobs = MEMBER+1
  end if


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_mem_node_proc:                           ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
! Start using SCALE library
!-------------------------------------------------------------------------------
subroutine set_scalelib

  use scale_stdio, only: &
    IO_LOG_setup, &
    IO_FID_CONF, &
    IO_FID_LOG, &
    IO_L, &
    H_LONG

  use gtool_history, only: &
    HistoryInit
  use dc_log, only: &
    LogInit
  use scale_process, only: &
    PRC_UNIVERSAL_setup, &
    PRC_MPIstart, &
    PRC_MPIsplit_letkf, &
    PRC_MPIsplit, &
    PRC_GLOBAL_setup, &
    PRC_LOCAL_setup, &
    PRC_masterrank, &
    PRC_myrank, &
    PRC_mpi_alive, &
    PRC_DOMAIN_nlim
  use scale_rm_process, only: &
    PRC_setup, &
    PRC_2Drank, &
    PRC_NUM_X, &
    PRC_NUM_Y 

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
  use scale_tracer, only: &
    TRACER_setup
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

!  use scale_gridtrans, only: &
!    GTRANS_setup

!  use scale_atmos_hydrostatic, only: &
!     ATMOS_HYDROSTATIC_setup
  use scale_atmos_thermodyn, only: &
     ATMOS_THERMODYN_setup

!  use mod_admin_time, only: &
!     ADMIN_TIME_setup


  use scale_mapproj, only: &
    MPRJ_setup
  implicit none



!    character(len=H_MID) :: DATATYPE = 'DEFAULT' !< REAL4 or REAL8
  integer :: rankidx(2)

  integer :: local_myrank
  logical :: local_ismaster

  CHARACTER(len=H_LONG) :: confname_dummy

!  integer :: universal_comm
!  integer :: universal_nprocs
!  logical :: universal_master
  integer :: global_comm
  integer :: local_comm
  integer :: intercomm_parent
  integer :: intercomm_child

  integer :: NUM_DOMAIN
  integer :: PRC_DOMAINS(PRC_DOMAIN_nlim)
  character(len=H_LONG) :: CONF_FILES (PRC_DOMAIN_nlim)


  integer :: ierr
  real(r_dble) :: rrtimer00, rrtimer
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  rrtimer00 = MPI_WTIME()


  !-----------------------------------------------------------------------------

  NUM_DOMAIN = 1
  PRC_DOMAINS = 0
  CONF_FILES = ""

  ! start SCALE MPI
!  call PRC_MPIstart( universal_comm ) ! [OUT]

  PRC_mpi_alive = .true.
!  universal_comm = MPI_COMM_WORLD

!  call PRC_UNIVERSAL_setup( universal_comm,   & ! [IN]
!                            universal_nprocs, & ! [OUT]
!                            universal_master  ) ! [OUT]

  ! split MPI communicator for LETKF
  call PRC_MPIsplit_letkf( MPI_COMM_WORLD,                   & ! [IN]
                           MEM_NP, nitmax, nprocs, proc2mem, & ! [IN]
                           global_comm                       ) ! [OUT]


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:prc_mpisplit_letkf:             ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  if (global_comm == MPI_COMM_NULL) then
!    write (6, '(A,I6.6,A)') 'MYRANK=',myrank,': This process is not used!'
    return
  end if

  call PRC_GLOBAL_setup( .false.,    & ! [IN]
                         global_comm ) ! [IN]


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:prc_global_setup:               ', rrtimer-rrtimer00
  call MPI_BARRIER(global_comm, ierr)
  rrtimer00 = MPI_WTIME()


  !--- split for nesting
  ! communicator split for nesting domains
  call PRC_MPIsplit( global_comm,      & ! [IN]
                     NUM_DOMAIN,       & ! [IN]
                     PRC_DOMAINS(:),   & ! [IN]
                     CONF_FILES(:),    & ! [IN]
                     .false.,          & ! [IN]
                     .false.,          & ! [IN] flag bulk_split
                     .false.,          & ! [IN] no reordering
                     local_comm,       & ! [OUT]
                     intercomm_parent, & ! [OUT]
                     intercomm_child,  & ! [OUT]
                     confname_dummy    ) ! [OUT]

  MPI_COMM_d = local_comm


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:prc_mpisplit_local:             ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  ! setup standard I/O
!  call IO_setup( MODELNAME, .true., cnf_fname )

  ! setup MPI
  call PRC_LOCAL_setup( local_comm, local_myrank, local_ismaster )


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:prc_local_setup:                ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  ! setup Log
  call IO_LOG_setup( local_myrank, local_ismaster )
  call LogInit( IO_FID_CONF, IO_FID_LOG, IO_L )


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:log_setup:                      ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  ! setup process
  call PRC_setup


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:prc_setup:                      ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  ! setup PROF
!  call PROF_setup

  ! setup constants
  call CONST_setup


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:const_setup:                    ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  ! setup calendar
!  call CALENDAR_setup

  ! setup random number
!  call RANDOM_setup

  ! setup time
!  call ADMIN_TIME_setup( setup_TimeIntegration = .true. )

!  call PROF_setprefx('INIT')
!  call PROF_rapstart('Initialize')

  ! setup horizontal/vertical grid coordinates
  call GRID_INDEX_setup
  call GRID_setup


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:grid_setup:                     ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


!  call LAND_GRID_INDEX_setup
!  call LAND_GRID_setup

!  call URBAN_GRID_INDEX_setup
!  call URBAN_GRID_setup

  ! setup tracer index
  call TRACER_setup


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:tracer_setup:                   ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  ! setup file I/O
  call FILEIO_setup


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:fileio_setup:                   ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  ! setup mpi communication
  call COMM_setup


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:comm_setup:                     ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  ! setup topography
!  call TOPO_setup
  ! setup land use category index/fraction
!  call LANDUSE_setup
  ! setup grid coordinates (real world)
!  call REAL_setup
    ! setup map projection [[ in REAL_setup ]]
    call MPRJ_setup( GRID_DOMAIN_CENTER_X, GRID_DOMAIN_CENTER_Y )


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:mprj_setup:                     ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


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

    call HistoryInit('', '', '', IMAX*JMAX*KMAX, PRC_masterrank, PRC_myrank, rankidx, &
                     0.0d0, 1.0d0, &
                     namelist_fid=IO_FID_CONF, default_basename='history')


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:history_init:                   ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  ! setup monitor I/O
!  call MONIT_setup

  ! setup nesting grid
!  call NEST_setup ( intercomm_parent, intercomm_child )

  ! setup common tools
!  call ATMOS_HYDROSTATIC_setup
  call ATMOS_THERMODYN_setup
!  call ATMOS_SATURATION_setup


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### set_scalelib:atmos_thermodyn_setup:          ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


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


  return
end subroutine set_scalelib

!-------------------------------------------------------------------------------
! Finish using SCALE library
!-------------------------------------------------------------------------------
subroutine unset_scalelib
  use gtool_file, only: &
    FileCloseAll
  use scale_stdio, only: &
    IO_FID_CONF, &
    IO_FID_LOG, &
    IO_L, &
    IO_FID_STDOUT
  implicit none

!  call MONIT_finalize

  call FileCloseAll

  ! Close logfile, configfile
  if ( IO_L ) then
    if( IO_FID_LOG /= IO_FID_STDOUT ) close(IO_FID_LOG)
  endif
  close(IO_FID_CONF)

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

  im = proc2mem(1,iter,myrank+1)
  if (im >= 1 .and. im <= nens) then
    if (im <= MEMBER) then
      call file_member_replace(im, HISTORY_IN_BASENAME, filename)
    else if (im == mmean) then
      filename = HISTORY_MEAN_IN_BASENAME
    else if (im == mmdet) then
      filename = HISTORY_MDET_IN_BASENAME
    end if

    call read_history(trim(filename), step, v3dg, v2dg)
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


  integer :: ierr
  real(r_dble) :: rrtimer00, rrtimer
!  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)

    ! Note: read all members + mdetin
    ! 
    if ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin) then
      if (im <= MEMBER) then
        call file_member_replace(im, GUES_IN_BASENAME, filename)
      else if (im == mmean) then
        filename = GUES_MEAN_INOUT_BASENAME
      else if (im == mmdet) then
        filename = GUES_MDET_IN_BASENAME
      end if

!      write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
      call read_restart(filename, v3dg, v2dg)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### read_ens_mpi:read_restart:                   ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


      call state_trans(v3dg)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### read_ens_mpi:state_trans:                    ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


    end if


  call MPI_BARRIER(MPI_COMM_e, ierr)
  rrtimer00 = MPI_WTIME()


    mstart = 1 + (it-1)*nprocs_e
    mend = min(it*nprocs_e, nens)
    if (mstart <= mend) then
      call scatter_grd_mpi_alltoall(mstart, mend, v3dg, v2dg, v3d, v2d)
    end if


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### read_ens_mpi:scatter_grd_mpi_alltoall:       ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


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
    im = proc2mem(1,it,myrank+1)

    ! Note: read all members
    ! 
    if (im >= 1 .and. im <= MEMBER) then
      call file_member_replace(im, INFL_ADD_IN_BASENAME, filename)

!      write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
      call read_restart(filename, v3dg, v2dg)
      call state_trans(v3dg)
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
subroutine write_ens_mpi(v3d, v2d, monit, caption)
  implicit none
  real(r_size), intent(in) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(in) :: v2d(nij1,nens,nv2d)
  logical, intent(in), optional :: monit
  character(len=*), intent(in), optional :: caption
  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  character(len=filelenmax) :: filename
  integer :: it, im, mstart, mend
  logical :: monit_


  integer :: ierr
  real(r_dble) :: rrtimer00, rrtimer
!  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


  monit_ = .false.
  if (present(monit) .and. present(caption)) then
    monit_ = monit
  end if

  do it = 1, nitmax


  call MPI_BARRIER(MPI_COMM_e, ierr)
  rrtimer00 = MPI_WTIME()


    im = proc2mem(1,it,myrank+1)

    mstart = 1 + (it-1)*nprocs_e
    mend = min(it*nprocs_e, nens)
    if (mstart <= mend) then
      call gather_grd_mpi_alltoall(mstart, mend, v3d, v2d, v3dg, v2dg)
    end if


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_ens_mpi:gather_grd_mpi_alltoall:       ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


    if (monit_ .and. mstart <= mmean .and. mmean <= mend) then
      call monit_obs_mpi(v3dg, v2dg, caption)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_ens_mpi:monit_obs_mpi:                 ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


    end if

    ! Note: write all members + mean + mdet
    ! 
    if (im >= 1 .and. im <= nens) then
      if (im <= MEMBER) then
        call file_member_replace(im, ANAL_OUT_BASENAME, filename)
      else if (im == mmean) then
        filename = ANAL_MEAN_OUT_BASENAME
      else if (im == mmdet) then
        filename = ANAL_MDET_OUT_BASENAME
      end if

!      write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
      call state_trans_inv(v3dg)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_ens_mpi:state_trans_inv:               ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


      call write_restart(filename, v3dg, v2dg)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_ens_mpi:write_restart:                 ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


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
subroutine monit_obs_mpi(v3dg, v2dg, caption)
  implicit none
  real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP), intent(in) :: v2dg(nlon,nlat,nv2d)
  character(len=*), intent(in) :: caption

  integer :: nobs(nid_obs)
  integer :: nobs_g(nid_obs)
  real(r_size) :: bias(nid_obs)
  real(r_size) :: bias_g(nid_obs)
  real(r_size) :: rmse(nid_obs)
  real(r_size) :: rmse_g(nid_obs)
  logical :: monit_type(nid_obs)
  integer :: i, ierr


  REAL(r_dble) :: rrtimer00, rrtimer
!  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


  ! NOTE: need to use 'mmean_rank_e' processes to run this calculation
  !       because only these processes have read topo files in 'topo2d'
  ! 
  if (myrank_e == mmean_rank_e) then
    call monit_obs(v3dg, v2dg, topo2d, nobs, bias, rmse, monit_type, .true.)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### monit_obs_mpi:monit_obs:                     ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


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

    call MPI_ALLREDUCE(MPI_IN_PLACE, nobs_g, nid_obs, MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, bias_g, nid_obs, MPI_r_size,  MPI_SUM, MPI_COMM_d, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, rmse_g, nid_obs, MPI_r_size,  MPI_SUM, MPI_COMM_d, ierr)

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


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### monit_obs_mpi:mpi_allreduce(domain):         ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  end if

  if (DEPARTURE_STAT_ALL_PROCESSES) then


  call MPI_BARRIER(MPI_COMM_e, ierr)
  rrtimer00 = MPI_WTIME()


    call MPI_BCAST(nobs,       nid_obs, MPI_INTEGER, mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(bias,       nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(rmse,       nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(nobs_g,     nid_obs, MPI_INTEGER, mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(bias_g,     nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(rmse_g,     nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(monit_type, nid_obs, MPI_LOGICAL, mmean_rank_e, MPI_COMM_e, ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### monit_obs_mpi:mpi_allreduce(ens):            ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  end if

  if (DEPARTURE_STAT_ALL_PROCESSES .or. myrank_e == mmean_rank_e) then
!    write(6,'(3A)') 'OBSERVATIONAL DEPARTURE STATISTICS (IN THIS SUBDOMAIN) [', trim(file_mean), ']:'
    write(6,'(2A)') trim(caption), ' (IN THIS SUBDOMAIN):'
    call monit_print(nobs, bias, rmse, monit_type)
!    write(6,'(3A)') 'OBSERVATIONAL DEPARTURE STATISTICS (GLOBAL) [', trim(file_mean), ']:'
    write(6,'(2A)') trim(caption), ' (GLOBAL):'
    call monit_print(nobs_g, bias_g, rmse_g, monit_type)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### monit_obs_mpi:monit_print:                   ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  end if

  return
end subroutine monit_obs_mpi

!-------------------------------------------------------------------------------
! Gather ensemble mean to {mmean_rank_e} and write SCALE restart files
!-------------------------------------------------------------------------------
subroutine write_ensmean(filename, v3d, v2d, calced, monit, caption)
  implicit none
  character(len=*), intent(in) :: filename
  real(r_size), intent(inout) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(inout) :: v2d(nij1,nens,nv2d)
  logical, intent(in), optional :: calced
  logical, intent(in), optional :: monit
  character(len=*), intent(in), optional :: caption

  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  logical :: calced_, monit_


  integer :: ierr
  real(r_dble) :: rrtimer00, rrtimer
!  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


  calced_ = .false.
  if (present(calced)) then
    calced_ = calced
  end if
  monit_ = .false.
  if (present(monit) .and. present(caption)) then
    monit_ = monit
  end if

  if (.not. calced) then
    call ensmean_grd(MEMBER, nij1, v3d, v2d)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_ensmean:ensmean_grd:                   ', rrtimer-rrtimer00
  rrtimer00 = MPI_WTIME()


  end if


  call MPI_BARRIER(MPI_COMM_e, ierr)
  rrtimer00 = MPI_WTIME()


  call gather_grd_mpi(mmean_rank_e, v3d(:,:,mmean,:), v2d(:,mmean,:), v3dg, v2dg)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_ensmean:gather_grd_mpi:                ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  if (monit_) then
    call monit_obs_mpi(v3dg, v2dg, caption)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_ensmean:monit_obs_mpi:                 ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  end if

  if (myrank_e == mmean_rank_e) then
    call state_trans_inv(v3dg)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_ensmean:state_trans_inv:               ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


    call write_restart(filename, v3dg, v2dg)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_ensmean:write_restart:                 ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


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


  integer :: ierr
  real(r_dble) :: rrtimer00, rrtimer
!  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


  call enssprd_grd(MEMBER, nij1, v3d, v2d, v3ds, v2ds)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_enssprd:enssprd_grd:                   ', rrtimer-rrtimer00
  call MPI_BARRIER(MPI_COMM_e, ierr)
  rrtimer00 = MPI_WTIME()


  call gather_grd_mpi(msprd_rank_e, v3ds, v2ds, v3dg, v2dg)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_enssprd:gather_grd_mpi:                ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  if (myrank_e == msprd_rank_e) then
!    call state_trans_inv(v3dg)              !! do not transform the spread output
    call write_restart(filename, v3dg, v2dg) !!


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### write_enssprd:write_restart:                 ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


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


  real(r_dble) :: rrtimer00, rrtimer
!  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


  if (myrank_a == 0) then
    call read_obs_all(obs)


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### read_obs_all_mpi:read_obs_all:               ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  end if


  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


  do iof = 1, OBS_IN_NUM
    call MPI_BCAST(obs(iof)%nobs, 1, MPI_INTEGER, 0, MPI_COMM_a, ierr)
    if (myrank_a /= 0) then
      call obs_info_allocate(obs(iof))
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


  rrtimer = MPI_WTIME()
  write (6,'(A,F15.7)') '###### read_obs_all_mpi:mpi_bcast:                  ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


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
!  if ((proc2mem(1,1,myrank+1) >= 1 .and. proc2mem(1,1,myrank+1) <= MEMBER) .or. &
!      proc2mem(1,1,myrank+1) == mmdetin) then
!    if (proc2mem(1,1,myrank+1) <= MEMBER) then
!      call file_member_replace(proc2mem(1,1,myrank+1), OBSDA_IN_BASENAME, obsdafile)
!    else if (proc2mem(1,1,myrank+1) == mmean) then
!      obsdafile = OBSDA_MEAN_IN_BASENAME
!    else if (proc2mem(1,1,myrank+1) == mmdet) then
!      obsdafile = OBSDA_MDET_IN_BASENAME
!    end if
!    write (obsda_suffix(2:7), '(I6.6)') proc2mem(2,1,myrank+1)
!#ifdef H08
!    call get_nobs(trim(obsdafile) // obsda_suffix, 8, nobs) ! H08
!#else
!    call get_nobs(trim(obsdafile) // obsda_suffix, 6, nobs)
!#endif
!  end if

! read by process 0 and broadcast
!-----------------------------
  if (myrank_e == 0) then
    write (obsda_suffix(2:7), '(I6.6)') proc2mem(2,1,myrank+1)
#ifdef H08
    call get_nobs(trim(obsdafile) // obsda_suffix, 8, nobs) ! H08
#else
    call get_nobs(trim(obsdafile) // obsda_suffix, 6, nobs)
#endif
  end if
  call MPI_BCAST(nobs, 1, MPI_INTEGER, 0, MPI_COMM_e, ierr)
!-----------------------------

  return
end subroutine get_nobs_da_mpi

!SUBROUTINE get_nobs_mpi(obsfile,nrec,nn)
!SUBROUTINE read_obs2_mpi(obsfile,nn,nbv,elem,rlon,rlat,rlev,odat,oerr,otyp,tdif,hdxf,iqc)
!SUBROUTINE allreduce_obs_mpi(n,nbv,hdxf,iqc)

!===============================================================================
END MODULE common_mpi_scale
