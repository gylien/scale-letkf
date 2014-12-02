module common_mpi_scale
!=======================================================================
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
!
!=======================================================================
!$USE OMP_LIB
  use common
  use common_mpi
  use common_letkf, only: nbv
  use common_scale
  use common_obs_scale

  use scale_precision
  use scale_stdio
!  use scale_prof
  use scale_grid_index

!  use common_scalelib

  implicit none
  public

!  integer,parameter :: mpibufsize=1000000
  integer,save :: nij1
  integer,save :: nij1max
  integer,allocatable,save :: nij1node(:)
!  real(r_size),allocatable,save :: phi1(:)
  real(r_size),allocatable,save :: lon1(:),lat1(:)
  real(r_size),allocatable,save :: lonu1(:),latu1(:)
  real(r_size),allocatable,save :: lonv1(:),latv1(:)
  real(r_size),allocatable,save :: ri1(:),rj1(:)
!  real(r_size),allocatable,save :: wg1(:)

  integer,save :: nitmax ! maximum number of model files processed by a process
  integer,allocatable,save :: procs(:)
  integer,allocatable,save :: mem2node(:,:)
  integer,allocatable,save :: mem2proc(:,:)
  integer,allocatable,save :: proc2mem(:,:,:)
  integer,save :: n_mem
  integer,save :: n_mempn

  integer,save :: scale_IO_group_n = -1
!  integer,save :: scale_IO_proc_n = -1
  logical,save :: valid_member = .false.

  integer,save :: MPI_COMM_e, nprocs_e, myrank_e



contains



!subroutine set_mpi_along_domains
!  use common_nml, only: &
!    MEM_NP

!  use gtool_history, only: &
!    historyinit
!  use scale_process, only: &
!    PRC_setup,    &
!    PRC_MPIstart, &
!!      PRC_mpifinish, &
!    PRC_master, &
!    PRC_myrank, &
!    PRC_2Drank, &
!    PRC_NUM_X, &
!    PRC_NUM_Y
!!    prc_nu, &
!  use scale_comm, only: &
!    COMM_setup
!  implicit none

!  integer :: rankidx(2)

!  !-----------------------------------------------------------------------------

!  ! setup mpi communication
!  call COMM_setup

!  ! setup history file I/O
!  rankidx(1) = PRC_2Drank(PRC_myrank, 1)
!  rankidx(2) = PRC_2Drank(PRC_myrank, 2)
!  call HistoryInit('','','',IMAX*JMAX*KMAX,PRC_master,PRC_myrank,rankidx)

!  ! check if the namelist seetings are consistent
!  if (MEM_NP /= PRC_NUM_X * PRC_NUM_Y) then
!    write(6,*) 'MEM_NP should be equal to PRC_NUM_X * PRC_NUM_Y.'
!    stop
!  end if

!  return
!end subroutine set_mpi_along_domains


subroutine rank_1d_2d(proc, pi, pj)
  use scale_process, only: PRC_2Drank
  implicit none
  integer, intent(in) :: proc
  integer, intent(out) :: pi, pj

  pi = PRC_2Drank(proc,1)
  pj = PRC_2Drank(proc,2)
end subroutine rank_1d_2d


subroutine rank_2d_1d(pi, pj, proc)
  use scale_process, only: PRC_NUM_X
  implicit none
  integer, intent(in) :: pi, pj
  integer, intent(out) :: proc

  proc = pj * PRC_NUM_X + pi
end subroutine rank_2d_1d

!-----------------------------------------------------------------------
!
! proc = -1: outside the global domain
!-----------------------------------------------------------------------
SUBROUTINE ijproc(rig,rjg,ri,rj,proc)
  use scale_grid_index, only: &
      IMAX,JMAX, &
      IHALO,JHALO, &
      IA,JA                  ! [for validation]
  use scale_process, only: &
      PRC_NUM_X,PRC_NUM_Y, &
      PRC_myrank             ! [for validation]
  use scale_grid, only: &    ! [for validation]
      GRID_CX, &             ! [for validation]
      GRID_CY, &             ! [for validation]
      GRID_CXG, &            ! [for validation]
      GRID_CYG, &            ! [for validation]
      DX, &                  ! [for validation]
      DY                     ! [for validation]
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rig
  REAL(r_size),INTENT(IN) :: rjg
  REAL(r_size),INTENT(OUT) :: ri
  REAL(r_size),INTENT(OUT) :: rj
  integer,INTENT(OUT) :: proc
  integer :: iproc, jproc

  if (rig < real(1+IHALO,r_size) .or. rig > real(IMAX*PRC_NUM_X+IHALO,r_size) .or. &
      rjg < real(1+JHALO,r_size) .or. rjg > real(JMAX*PRC_NUM_Y+JHALO,r_size)) then
    ri = -1.0d0
    rj = -1.0d0
    proc = -1
    return
  end if

  iproc = ceiling((rig-real(IHALO,r_size)-0.5d0) / real(IMAX,r_size)) - 1
  jproc = ceiling((rjg-real(JHALO,r_size)-0.5d0) / real(JMAX,r_size)) - 1
  ri = rig - iproc * IMAX
  rj = rjg - jproc * JMAX
  call rank_2d_1d(iproc,jproc,proc)

  if (PRC_myrank == proc) then                                                                                    ! [for validation]
    if (rig < (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0 .or. &                                                      ! [for validation]
        rig > (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0 .or. &                                                     ! [for validation]
        rjg < (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0 .or. &                                                      ! [for validation]
        rjg > (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0) then                                                      ! [for validation]
      write (6,'(A)') 'Error: Process assignment fails!'                                                          ! [for validation]
      write (6,'(3F10.2)') rig, (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0, (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0 ! [for validation]
      write (6,'(3F10.2)') rjg, (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0, (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0 ! [for validation]
      stop                                                                                                        ! [for validation]
    end if                                                                                                        ! [for validation]
  end if                                                                                                          ! [for validation]

  RETURN
END SUBROUTINE ijproc

!-----------------------------------------------------------------------
! set_common_mpi_scale
!-----------------------------------------------------------------------
SUBROUTINE set_common_mpi_scale(nbv,nnodes,ppn,mem_nodes,mem_np)
  INTEGER,INTENT(IN) :: nbv
  INTEGER,INTENT(IN) :: nnodes,ppn,mem_nodes,mem_np
  REAL(r_sngl) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:)
  INTEGER :: i,n
  INTEGER :: ierr,buf(4)
  LOGICAL :: ex


!  if (myrank == 0) then
!  print *, procs
!  print *, mem2node
!  print *, mem2proc
!  print *, proc2mem
!  end if
  integer :: MPI_G_WORLD, MPI_G
  integer :: n_mem,n_mempn
!  integer :: iproc,jproc
  integer,allocatable :: ranks(:)

  integer :: ip



  WRITE(6,'(A)') 'Hello from set_common_mpi_scale'


  CALL set_mem_node_proc(nbv+1,NNODES,PPN,MEM_NODES,MEM_NP)


!!!!!!------
  if (scale_IO_group_n >= 1) then

    call set_scalelib(MEM_NP, nitmax, nprocs, proc2mem)


    IF(MEM_NODES > 1) THEN
      n_mem = NNODES / MEM_NODES
      n_mempn = 1
    ELSE
      n_mem = NNODES
      n_mempn = PPN / MEM_NP
    END IF
    nprocs_e = n_mem*n_mempn

    allocate (ranks(nprocs_e))

    call MPI_Comm_group(MPI_COMM_WORLD,MPI_G_WORLD,ierr)

    do ip = 1, nprocs
      if (proc2mem(2,1,ip) == proc2mem(2,1,myrank+1)) then
        if (proc2mem(1,1,ip) >= 1) then
          ranks(proc2mem(1,1,ip)) = ip-1
        end if
      end if
    end do

!write(6,'(A,7I6)') '######===', myrank, ranks(:)

    call MPI_Group_incl(MPI_G_WORLD,nprocs_e,ranks,MPI_G,ierr)
    call MPI_Comm_create(MPI_COMM_WORLD,MPI_G,MPI_COMM_e,ierr)

    call MPI_Comm_size(MPI_COMM_e,nprocs_e,ierr)
    call MPI_Comm_rank(MPI_COMM_e,myrank_e,ierr)

!write(6,'(A,9I6)') '######===', myrank, myrank_e, nprocs_e, ranks(:)
!stop

    deallocate(ranks)
  end if
!!!!!!------




!  i = MOD(nlon*nlat,nprocs)
!  nij1max = (nlon*nlat - i)/nprocs + 1
!  IF(myrank < i) THEN
!    nij1 = nij1max
!  ELSE
!    nij1 = nij1max - 1
!  END IF
!  WRITE(6,'(A,I3.3,A,I6)') 'MYRANK ',myrank,' number of grid points: nij1= ',nij1
!  ALLOCATE(nij1node(nprocs))
!  DO n=1,nprocs
!    IF(n-1 < i) THEN
!      nij1node(n) = nij1max
!    ELSE
!      nij1node(n) = nij1max - 1
!    END IF
!  END DO



!!  ALLOCATE(phi1(nij1))
!  ALLOCATE(lon1(nij1))
!  ALLOCATE(lat1(nij1))
!  ALLOCATE(lonu1(nij1))
!  ALLOCATE(latu1(nij1))
!  ALLOCATE(lonv1(nij1))
!  ALLOCATE(latv1(nij1))
!  ALLOCATE(ri1(nij1))
!  ALLOCATE(rj1(nij1))
!!  ALLOCATE(wg1(nij1))

!  ALLOCATE(v3d(nij1,nlev,nv3d))
!  ALLOCATE(v2d(nij1,nv2d))
!  v3dg(:,:,1,1) = SNGL(lon)
!  v3dg(:,:,1,2) = SNGL(lat)
!  v3dg(:,:,1,3) = SNGL(lonu)
!  v3dg(:,:,1,4) = SNGL(latu)
!  v3dg(:,:,1,5) = SNGL(lonv)
!  v3dg(:,:,1,6) = SNGL(latv)
!!  v3dg(:,:,2,1) = SNGL(wg(:,:))
!  DO j=1,nlat
!    DO i=1,nlon
!      v3dg(i,j,1,7) = REAL(i,r_sngl)
!      v3dg(i,j,1,8) = REAL(j,r_sngl)
!    END DO
!  END DO
!!  v2dg(:,:,1) = SNGL(phi0)
!  CALL scatter_grd_mpi(0,v3dg,v2dg,v3d,v2d)
!  lon1  = v3d(:,1,1)
!  lat1  = v3d(:,1,2)
!  lonu1 = v3d(:,1,3)
!  latu1 = v3d(:,1,4)
!  lonv1 = v3d(:,1,5)
!  latv1 = v3d(:,1,6)
!  ri1   = v3d(:,1,7)
!  rj1   = v3d(:,1,8)
!  phi1  = v2d(:,1)
!!  wg1(:) = v3d(:,2,1)




  RETURN
END SUBROUTINE set_common_mpi_scale
!-----------------------------------------------------------------------
! set_common_mpi_scale
!-----------------------------------------------------------------------
SUBROUTINE unset_common_mpi_scale
  implicit none
  integer:: ierr

  if (scale_IO_group_n >= 1) then
    call MPI_Comm_free(MPI_COMM_e,ierr)
    call unset_scalelib
  end if

  RETURN
END SUBROUTINE unset_common_mpi_scale
!-----------------------------------------------------------------------
! set_mem2proc
!-----------------------------------------------------------------------
SUBROUTINE set_mem_node_proc(mem,nnodes,ppn,mem_nodes,mem_np)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: mem,nnodes,ppn,mem_nodes,mem_np
  INTEGER :: tppn,tppnt,tmod
  INTEGER :: n,ns,nn,m,q,qs,i,j,it,ip

  ALLOCATE(procs(nprocs))
  ns = 0
  DO n = 1, nnodes
    procs(ns+1:ns+ppn) = n
    ns = ns + ppn
  END DO

  IF(mem_nodes > 1) THEN
    n_mem = nnodes / mem_nodes
    n_mempn = 1
  ELSE
    n_mem = nnodes
    n_mempn = ppn / mem_np
  END IF
  nitmax = (mem - 1) / (n_mem * n_mempn) + 1
  tppn = mem_np / mem_nodes
  tmod = MOD(mem_np, mem_nodes)

  ALLOCATE(mem2node(mem_np,mem))
  ALLOCATE(mem2proc(mem_np,mem))
  ALLOCATE(proc2mem(2,nitmax,nprocs))
  proc2mem = -1
  m = 1
mem_loop: DO it = 1, nitmax
    DO i = 0, n_mempn-1
      n = 0
      DO j = 0, n_mem-1
        IF(m > mem .and. it > 1) EXIT mem_loop
        qs = 0
        DO nn = 0, mem_nodes-1
          IF(nn < tmod) THEN
            tppnt = tppn + 1
          ELSE
            tppnt = tppn
          END IF
          DO q = 0, tppnt-1
            ip = (n+nn)*ppn + i*mem_np + q
            if (m <= mem) then
              mem2node(qs+1,m) = n+nn
              mem2proc(qs+1,m) = ip
            end if
            proc2mem(1,it,ip+1) = m
            proc2mem(2,it,ip+1) = qs
            qs = qs + 1
          END DO
        END DO
        m = m + 1
        n = n + mem_nodes
      END DO
    END DO
  END DO mem_loop

  scale_IO_group_n = proc2mem(1,1,myrank+1)
!  scale_IO_proc_n = proc2mem(2,1,myrank+1)
  if (scale_IO_group_n >= 1 .and. scale_IO_group_n <= mem) then
    valid_member = .true.
  end if

  RETURN
END SUBROUTINE
!-----------------------------------------------------------------------
! Scatter gridded data to processes (nrank -> all)
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)

!  IF(mpibufsize > nij1max) THEN
    CALL scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
!  ELSE
!    CALL scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
!  END IF

  RETURN
END SUBROUTINE scatter_grd_mpi

!SUBROUTINE scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
!  INTEGER,INTENT(IN) :: nrank
!  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
!  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
!  REAL(r_sngl) :: tmp(nij1max,nprocs)
!  REAL(r_sngl) :: bufs(mpibufsize,nprocs)
!  REAL(r_sngl) :: bufr(mpibufsize)
!  INTEGER :: i,j,k,n,ierr,ns,nr
!  INTEGER :: iter,niter

!  ns = mpibufsize
!  nr = ns
!  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))

!  DO n=1,nv3d
!    DO k=1,nlev
!      IF(myrank == nrank) CALL grd_to_buf(v3dg(:,:,k,n),tmp)
!      DO iter=1,niter
!        IF(myrank == nrank) THEN
!          i = mpibufsize * (iter-1)
!          DO j=1,mpibufsize
!            i=i+1
!            IF(i > nij1max) EXIT
!            bufs(j,:) = tmp(i,:)
!          END DO
!        END IF
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
!                       & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!        i = mpibufsize * (iter-1)
!        DO j=1,mpibufsize
!          i=i+1
!          IF(i > nij1) EXIT
!          v3d(i,k,n) = REAL(bufr(j),r_size)
!        END DO
!      END DO
!    END DO
!  END DO

!  DO n=1,nv2d
!    IF(myrank == nrank) CALL grd_to_buf(v2dg(:,:,n),tmp)
!    DO iter=1,niter
!      IF(myrank == nrank) THEN
!        i = mpibufsize * (iter-1)
!        DO j=1,mpibufsize
!          i=i+1
!          IF(i > nij1max) EXIT
!          bufs(j,:) = tmp(i,:)
!        END DO
!      END IF
!      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
!                     & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!      i = mpibufsize * (iter-1)
!      DO j=1,mpibufsize
!        i=i+1
!        IF(i > nij1) EXIT
!        v2d(i,n) = REAL(bufr(j),r_size)
!      END DO
!    END DO
!  END DO

!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!  RETURN
!END SUBROUTINE scatter_grd_mpi_safe

SUBROUTINE scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
  REAL(r_sngl) :: bufr(nij1max,nlevall)
  INTEGER :: j,k,n,ierr,ns,nr

  ns = nij1max * nlevall
  nr = ns
  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

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

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE scatter_grd_mpi_fast
!-----------------------------------------------------------------------
! Gather gridded data (all -> nrank)
!-----------------------------------------------------------------------
SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)

!  IF(mpibufsize > nij1max) THEN
    CALL gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
!  ELSE
!    CALL gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
!  END IF

  RETURN
END SUBROUTINE gather_grd_mpi

!SUBROUTINE gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
!  INTEGER,INTENT(IN) :: nrank
!  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
!  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
!  REAL(r_sngl) :: tmp(nij1max,nprocs)
!  REAL(r_sngl) :: bufs(mpibufsize)
!  REAL(r_sngl) :: bufr(mpibufsize,nprocs)
!  INTEGER :: i,j,k,n,ierr,ns,nr
!  INTEGER :: iter,niter

!  ns = mpibufsize
!  nr = ns
!  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))

!  DO n=1,nv3d
!    DO k=1,nlev
!      DO iter=1,niter
!        i = mpibufsize * (iter-1)
!        DO j=1,mpibufsize
!          i=i+1
!          IF(i > nij1) EXIT
!          bufs(j) = REAL(v3d(i,k,n),r_sngl)
!        END DO
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        CALL MPI_GATHER(bufs,ns,MPI_REAL,&
!                      & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!        IF(myrank == nrank) THEN
!          i = mpibufsize * (iter-1)
!          DO j=1,mpibufsize
!            i=i+1
!            IF(i > nij1max) EXIT
!            tmp(i,:) = bufr(j,:)
!          END DO
!        END IF
!      END DO
!      IF(myrank == nrank) CALL buf_to_grd(tmp,v3dg(:,:,k,n))
!    END DO
!  END DO

!  DO n=1,nv2d
!    DO iter=1,niter
!      i = mpibufsize * (iter-1)
!      DO j=1,mpibufsize
!        i=i+1
!        IF(i > nij1) EXIT
!        bufs(j) = REAL(v2d(i,n),r_sngl)
!      END DO
!      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      CALL MPI_GATHER(bufs,ns,MPI_REAL,&
!                    & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!      IF(myrank == nrank) THEN
!        i = mpibufsize * (iter-1)
!        DO j=1,mpibufsize
!          i=i+1
!          IF(i > nij1max) EXIT
!          tmp(i,:) = bufr(j,:)
!        END DO
!      END IF
!    END DO
!    IF(myrank == nrank) CALL buf_to_grd(tmp,v2dg(:,:,n))
!  END DO

!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!  RETURN
!END SUBROUTINE gather_grd_mpi_safe

SUBROUTINE gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl) :: bufs(nij1max,nlevall)
  REAL(r_sngl) :: bufr(nij1max,nlevall,nprocs)
  INTEGER :: j,k,n,ierr,ns,nr

  ns = nij1max * nlevall
  nr = ns
  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      bufs(1:nij1,j) = REAL(v3d(:,k,n),r_sngl)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    bufs(1:nij1,j) = REAL(v2d(:,n),r_sngl)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,n))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,n))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE gather_grd_mpi_fast
!-----------------------------------------------------------------------
! Read ensemble data and distribute to processes
!-----------------------------------------------------------------------

!subroutine set_scale_mpi_comm
!  implicit none

!!    real(RP), allocatable :: U(:,:,:), MOMX(:,:,:)

!  integer :: k, i, j

!  integer               :: dim1_max, dim1_S, dim1_E
!  integer               :: dim2_max, dim2_S, dim2_E
!  integer               :: dim3_max, dim3_S, dim3_E
!  integer               :: dim4_max, dim4_S, dim4_E
!  real(RP), allocatable :: var3D(:,:,:)

!  integer :: iolen

!  character(len=100) :: basename
!  character(len=100) :: varname
!  integer :: step
!  basename = 'history'
!  varname = 'U'
!  step = 1

!  if (PRC_nu == 0) then
!    basename = trim(basename) // '.u000000'
!  else if (PRC_nu == 1) THEN
!    basename = trim(basename) // '.u000001'
!  end if

!  allocate( U(KA,IA,JA) )
!  allocate( MOMX(KA,IA,JA) )
!  U = 0.0d0




!    ! Read file

!!    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(varname)

!       dim1_max = IMAX !KMAX
!       dim2_max = JMAX !IMAX
!       dim3_max = KMAX !JMAX
!       dim1_S   = IS !KS
!       dim1_E   = IE !KE
!       dim2_S   = JS !IS
!       dim2_E   = JE !IE
!       dim3_S   = KS !JS
!       dim3_E   = KE !JE

!    allocate( var3D(dim1_max,dim2_max,dim3_max) )


!    call HIST_get(var3D, trim(basename), trim(varname), step=step)


!    forall (i=1:IMAX, j=1:JMAX, k=1:KMAX) U(k+KHALO,i+IHALO,j+JHALO) = var3D(i,j,k)

!    deallocate( var3D )


!    if (PRC_nu == 0) then
!      call FILEIO_read( MOMX(:,:,:),                          & ! [OUT]
!                        'init.u000000', 'MOMX', 'ZXY', step=1 ) ! [IN]
!    else if (PRC_nu == 1) then
!      call FILEIO_read( MOMX(:,:,:),                          & ! [OUT]
!                        'init.u000001', 'MOMX', 'ZXY', step=1 ) ! [IN]
!    end if


!    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
!    do j  = JS, JE
!    do i  = IS, IE
!       U(   1:KS-1,i,j) = U(KS,i,j)
!       U(KE+1:KA,  i,j) = U(KE,i,j)
!       MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
!       MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
!    enddo
!    enddo

!    call COMM_vars8( U   (:,:,:), 1 )
!    call COMM_vars8( MOMX(:,:,:), 2 )
!    call COMM_wait ( U   (:,:,:), 1 )
!    call COMM_wait ( MOMX(:,:,:), 2 )


SUBROUTINE read_ens_history_mpi(file,iter,step,v3dg,v2dg,ensmean)
!    use scale_fileio, only: &
!       FILEIO_write, &
!       FILEIO_read
    use scale_history, only: &
       HIST_get
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait

  IMPLICIT NONE

  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: iter
  INTEGER,INTENT(IN) :: step
  REAL(r_size),INTENT(OUT) :: v3dg(nlevhalo,nlonhalo,nlathalo,nv3dd)
  REAL(r_size),INTENT(OUT) :: v2dg(nlonhalo,nlathalo,nv2dd)
  LOGICAL,INTENT(INOUT),OPTIONAL :: ensmean
  INTEGER :: i,j,k,iv3d,iv2d,nbvr

  CHARACTER(9) :: filename='file.0000'

  real(RP), allocatable :: var3D(:,:,:)
  real(RP), allocatable :: var2D(:,:)

  if (.not. present(ensmean) .or. (.not. ensmean)) then
    nbvr = nbv
  else
    nbvr = nbv+1
  end if

  IF (valid_member) then
    allocate( var3D(IMAX,JMAX,KMAX) )
    allocate( var2D(IMAX,JMAX) )

    IF(proc2mem(1,iter,myrank+1) >= 1 .and. proc2mem(1,iter,myrank+1) <= nbvr) THEN
      WRITE(filename(1:4),'(A4)') file

      if (proc2mem(1,iter,myrank+1) == nbv+1) then
        WRITE(filename(6:9),'(A4)') 'mean'
      else
        WRITE(filename(6:9),'(I4.4)') proc2mem(1,iter,myrank+1)
      end if

      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',proc2mem(2,iter,myrank+1),'.nc'

      DO iv3d = 1, nv3dd
        if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(v3dd_name(iv3d))
        call HIST_get(var3D, filename, trim(v3dd_name(iv3d)), step)
        FORALL (i=1:IMAX, j=1:JMAX, k=1:KMAX) v3dg(k+KHALO,i+IHALO,j+JHALO,iv3d) = var3D(i,j,k)

!!!!!$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
        do j  = JS, JE
          do i  = IS, IE
            v3dg(   1:KS-1,i,j,iv3d) = v3dg(KS,i,j,iv3d)
            v3dg(KE+1:KA,  i,j,iv3d) = v3dg(KE,i,j,iv3d)
          enddo
        enddo
      END DO

      DO iv3d = 1, nv3dd
        call COMM_vars8( v3dg(:,:,:,iv3d), iv3d )
      END DO
      DO iv3d = 1, nv3dd
        call COMM_wait ( v3dg(:,:,:,iv3d), iv3d )
      END DO

      DO iv2d = 1, nv2dd
        if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 2D var: ', trim(v2dd_name(iv2d))
        call HIST_get(var2D, filename, trim(v2dd_name(iv2d)), step)
        v2dg(1+IHALO:IMAX+IHALO,1+JHALO:JMAX+JHALO,iv2d) = var2D(:,:)
      END DO

      DO iv2d = 1, nv2dd
        call COMM_vars8( v2dg(:,:,iv2d), iv2d )
      END DO
      DO iv2d = 1, nv2dd
        call COMM_wait ( v2dg(:,:,iv2d), iv2d )
      END DO
    END IF

    deallocate( var3D )
    deallocate( var2D )
  END IF

  RETURN
END SUBROUTINE read_ens_history_mpi
!!-----------------------------------------------------------------------
!! Write ensemble data after collecting data from processes
!!-----------------------------------------------------------------------
!SUBROUTINE write_ens_mpi(file,member,v3d,v2d)
!  CHARACTER(4),INTENT(IN) :: file
!  INTEGER,INTENT(IN) :: member
!  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
!  REAL(r_sngl) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl) :: v2dg(nlon,nlat,nv2d)
!  INTEGER :: l,n,ll,im
!  CHARACTER(11) :: filename='file000.grd'

!  ll = CEILING(REAL(member)/REAL(nprocs))
!  DO l=1,ll
!    DO n=0,nprocs-1
!      im = n+1 + (l-1)*nprocs
!      IF(im <= member) THEN
!        CALL gather_grd_mpi(n,v3d(:,:,im,:),v2d(:,im,:),v3dg,v2dg)
!      END IF
!    END DO

!    im = myrank+1 + (l-1)*nprocs
!    IF(im <= member) THEN
!      WRITE(filename(1:7),'(A4,I3.3)') file,im
!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
!      CALL write_grd4(filename,v3dg,v2dg,0)
!    END IF
!  END DO

!  RETURN
!END SUBROUTINE write_ens_mpi
!-----------------------------------------------------------------------
! gridded data -> buffer
!-----------------------------------------------------------------------
SUBROUTINE grd_to_buf(grd,buf)
  REAL(r_sngl),INTENT(IN) :: grd(nlon,nlat)
  REAL(r_sngl),INTENT(OUT) :: buf(nij1max,nprocs)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      buf(i,m) = grd(ilon,ilat)
    END DO
  END DO

  DO m=1,nprocs
    IF(nij1node(m) < nij1max) buf(nij1max,m) = undef
  END DO

  RETURN
END SUBROUTINE grd_to_buf
!-----------------------------------------------------------------------
! buffer -> gridded data
!-----------------------------------------------------------------------
SUBROUTINE buf_to_grd(buf,grd)
  REAL(r_sngl),INTENT(IN) :: buf(nij1max,nprocs)
  REAL(r_sngl),INTENT(OUT) :: grd(nlon,nlat)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      grd(ilon,ilat) = buf(i,m)
    END DO
  END DO

  RETURN
END SUBROUTINE buf_to_grd
!!-----------------------------------------------------------------------
!! STORING DATA (ensemble mean and spread)
!!-----------------------------------------------------------------------
!SUBROUTINE write_ensmspr_mpi(file,member,v3d,v2d)
!  CHARACTER(4),INTENT(IN) :: file
!  INTEGER,INTENT(IN) :: member
!  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
!  REAL(r_size) :: v3dm(nij1,nlev,nv3d)
!  REAL(r_size) :: v2dm(nij1,nv2d)
!  REAL(r_size) :: v3ds(nij1,nlev,nv3d)
!  REAL(r_size) :: v2ds(nij1,nv2d)
!  REAL(r_sngl) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl) :: v2dg(nlon,nlat,nv2d)
!  INTEGER :: i,k,m,n
!  CHARACTER(11) :: filename='file000.grd'

!  CALL ensmean_grd(member,nij1,v3d,v2d,v3dm,v2dm)

!  CALL gather_grd_mpi(0,v3dm,v2dm,v3dg,v2dg)
!  IF(myrank == 0) THEN
!!  IF(myrank == nprocs-1) THEN  ! The last processor
!    WRITE(filename(1:7),'(A4,A3)') file,'_me'
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
!    CALL write_grd4(filename,v3dg,v2dg,0)
!  END IF

!  DO n=1,nv3d
!!$OMP PARALLEL DO PRIVATE(i,k,m)
!    DO k=1,nlev
!      DO i=1,nij1
!        v3ds(i,k,n) = (v3d(i,k,1,n)-v3dm(i,k,n))**2
!        DO m=2,member
!          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n)-v3dm(i,k,n))**2
!        END DO
!        v3ds(i,k,n) = SQRT(v3ds(i,k,n) / REAL(member-1,r_size))
!      END DO
!    END DO
!!$OMP END PARALLEL DO
!  END DO

!  DO n=1,nv2d
!!$OMP PARALLEL DO PRIVATE(i,k,m)
!    DO i=1,nij1
!      v2ds(i,n) = (v2d(i,1,n)-v2dm(i,n))**2
!      DO m=2,member
!        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n)-v2dm(i,n))**2
!      END DO
!      v2ds(i,n) = SQRT(v2ds(i,n) / REAL(member-1,r_size))
!    END DO
!!$OMP END PARALLEL DO
!  END DO

!  CALL gather_grd_mpi(0,v3ds,v2ds,v3dg,v2dg)
!  IF(myrank == 0) THEN
!!  IF(myrank == MOD(nprocs*2-2,nprocs)) THEN  ! The second last processor
!    WRITE(filename(1:7),'(A4,A3)') file,'_sp'
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
!    CALL write_grd4(filename,v3dg,v2dg,0)
!  END IF

!  RETURN
!END SUBROUTINE write_ensmspr_mpi
!!-----------------------------------------------------------------------
!! Get number of observations from ensemble obs2 data,
!! assuming all members have the identical obs records
!!  -- 12/30/2013, Guo-Yuan Lien
!!-----------------------------------------------------------------------
!SUBROUTINE get_nobs_mpi(obsfile,nrec,nn)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: obsfile
!  INTEGER,INTENT(IN) :: nrec
!  INTEGER,INTENT(OUT) :: nn
!  CHARACTER(LEN=LEN(obsfile)) :: obsfile1
!  INTEGER :: ms1,ms2,ierr

!  IF(myrank == 0) THEN
!    ms1 = LEN(obsfile)-6
!    ms2 = LEN(obsfile)-4
!    obsfile1 = obsfile
!    WRITE(obsfile1(ms1:ms2),'(I3.3)') 1
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',obsfile1
!    CALL get_nobs(obsfile1,nrec,nn)
!  END IF
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(nn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!  RETURN
!END SUBROUTINE get_nobs_mpi
!!-----------------------------------------------------------------------
!! Read ensemble obs2 observation data and ALLREDUCE of hdxf and qc
!!  -- 12/30/2013, Guo-Yuan Lien (do not consider mpibufsize)
!!-----------------------------------------------------------------------
!SUBROUTINE read_obs2_mpi(obsfile,nn,nbv,elem,rlon,rlat,rlev,odat,oerr,otyp,tdif,hdxf,iqc)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: obsfile
!  INTEGER,INTENT(IN) :: nn
!  INTEGER,INTENT(IN) :: nbv
!  REAL(r_size),INTENT(OUT) :: elem(nn)
!  REAL(r_size),INTENT(OUT) :: rlon(nn)
!  REAL(r_size),INTENT(OUT) :: rlat(nn)
!  REAL(r_size),INTENT(OUT) :: rlev(nn)
!  REAL(r_size),INTENT(OUT) :: odat(nn)
!  REAL(r_size),INTENT(OUT) :: oerr(nn)
!  REAL(r_size),INTENT(OUT) :: otyp(nn)
!  REAL(r_size),INTENT(OUT) :: tdif(nn)
!  REAL(r_size),INTENT(OUT) :: hdxf(nn,nbv)
!  INTEGER,INTENT(OUT) :: iqc(nn,nbv)
!  CHARACTER(LEN=LEN(obsfile)) :: obsfile1
!  INTEGER :: l,im,n
!  INTEGER :: MPI_C,MPI_G,MPI_G_WORLD,ierr
!  INTEGER :: ms1,ms2
!  INTEGER,ALLOCATABLE :: useranks(:)

!  hdxf = 0.0d0
!  iqc = 0
!  ms1 = LEN(obsfile)-6
!  ms2 = LEN(obsfile)-4
!  obsfile1 = obsfile
!  l=0
!  DO
!    im = myrank+1 + nprocs * l
!    IF(im > nbv) EXIT
!    WRITE(obsfile1(ms1:ms2),'(I3.3)') im
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',obsfile1
!    CALL read_obs2(obsfile1,nn,elem,rlon,rlat,rlev,odat,oerr,otyp,tdif,hdxf(:,im),iqc(:,im))
!    l = l+1
!  END DO
!!
!! if the total number of processors is greater then the ensemble size,
!! broadcast the first 8 observation records(elm/lon/lat/.../dif)
!! from myrank=nbv-1 to the rest of processors that didn't read anything.
!!
!  IF(nprocs > nbv) THEN
!    ALLOCATE(useranks(nprocs-nbv+1))
!    do n = nbv, nprocs
!      useranks(n-nbv+1) = n-1
!    end do
!    call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_G_WORLD,ierr)
!    call MPI_GROUP_INCL(MPI_G_WORLD,nprocs-nbv+1,useranks,MPI_G,ierr)
!    call MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_G,MPI_C,ierr)

!    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    IF(myrank+1 >= nbv) THEN
!      CALL MPI_BCAST(elem,nn,MPI_r_size,0,MPI_C,ierr)
!      CALL MPI_BCAST(rlon,nn,MPI_r_size,0,MPI_C,ierr)
!      CALL MPI_BCAST(rlat,nn,MPI_r_size,0,MPI_C,ierr)
!      CALL MPI_BCAST(rlev,nn,MPI_r_size,0,MPI_C,ierr)
!      CALL MPI_BCAST(odat,nn,MPI_r_size,0,MPI_C,ierr)
!      CALL MPI_BCAST(oerr,nn,MPI_r_size,0,MPI_C,ierr)
!      CALL MPI_BCAST(otyp,nn,MPI_r_size,0,MPI_C,ierr)
!      CALL MPI_BCAST(tdif,nn,MPI_r_size,0,MPI_C,ierr)
!    END IF
!    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    DEALLOCATE(useranks)
!  END IF

!  CALL allreduce_obs_mpi(nn,nbv,hdxf,iqc)

!  RETURN
!END SUBROUTINE read_obs2_mpi
!!!!!!!!!-----------------------------------------------------------------------
!!!!!!!!!
!!!!!!!!!-----------------------------------------------------------------------
!!!!!!!!SUBROUTINE obs_info_allreduce(obs)
!!!!!!!!  IMPLICIT NONE
!!!!!!!!  TYPE(obs_info),INTENT(INOUT) :: obs

!!!!!!!!  ALLOCATE( obs%dat (obs%nobs) )

!!!!!!!!  CALL MPI_ALLREDUCE(ibufs,ibufr,obs%nobs,MPI_INTEGER,MPI_MAX,&
!!!!!!!!          & MPI_COMM_WORLD,ierr)


!!!!!!!!  RETURN
!!!!!!!!END SUBROUTINE obs_info_allreduce
!!-----------------------------------------------------------------------
!! MPI_ALLREDUCE of hdxf and qc
!!-----------------------------------------------------------------------
!SUBROUTINE allreduce_obs_mpi(n,nbv,hdxf,iqc)
!  INTEGER,INTENT(IN) :: n
!  INTEGER,INTENT(IN) :: nbv
!  REAL(r_size),INTENT(INOUT) :: hdxf(n,nbv)
!  INTEGER,INTENT(INOUT) :: iqc(n,nbv)
!  REAL(r_size) :: bufs(mpibufsize)
!  REAL(r_size) :: bufr(mpibufsize)
!  REAL(r_size),ALLOCATABLE :: tmp(:,:)
!  INTEGER :: ibufs(mpibufsize)
!  INTEGER :: ibufr(mpibufsize)
!  INTEGER,ALLOCATABLE :: itmp(:,:)
!  INTEGER :: i,j,k
!  INTEGER :: iter,niter
!  INTEGER :: ierr

!  niter = CEILING(REAL(n*nbv)/REAL(mpibufsize))
!  ALLOCATE(tmp(mpibufsize,niter))
!  ALLOCATE(itmp(mpibufsize,niter))
!  bufs=0.0d0
!  ibufs=0
!  i=1
!  iter=1
!  DO k=1,nbv
!    DO j=1,n
!      bufs(i) = hdxf(j,k)
!      ibufs(i) = iqc(j,k)
!      i=i+1
!      IF(i > mpibufsize) THEN
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        CALL MPI_ALLREDUCE(bufs,bufr,mpibufsize,MPI_r_size,MPI_SUM,&
!          & MPI_COMM_WORLD,ierr)
!        tmp(:,iter) = bufr
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        CALL MPI_ALLREDUCE(ibufs,ibufr,mpibufsize,MPI_INTEGER,MPI_MAX,&
!          & MPI_COMM_WORLD,ierr)
!        itmp(:,iter) = ibufr
!        i=1
!        iter=iter+1
!      END IF
!    END DO
!  END DO
!  IF(iter == niter) THEN
!    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    CALL MPI_ALLREDUCE(bufs,bufr,mpibufsize,MPI_r_size,MPI_SUM,&
!      & MPI_COMM_WORLD,ierr)
!    tmp(:,iter) = bufr
!    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    CALL MPI_ALLREDUCE(ibufs,ibufr,mpibufsize,MPI_INTEGER,MPI_MAX,&
!      & MPI_COMM_WORLD,ierr)
!    itmp(:,iter) = ibufr
!  END IF

!  i=1
!  iter=1
!  DO k=1,nbv
!    DO j=1,n
!      hdxf(j,k) = tmp(i,iter)
!      iqc(j,k) = itmp(i,iter)
!      i=i+1
!      IF(i > mpibufsize) THEN
!        i=1
!        iter=iter+1
!      END IF
!    END DO
!  END DO
!  DEALLOCATE(tmp,itmp)

!  RETURN
!END SUBROUTINE allreduce_obs_mpi

END MODULE common_mpi_scale
