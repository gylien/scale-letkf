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

  use scale_precision, only: RP
  use scale_stdio
!  use scale_prof
!  use scale_grid_index

!  use common_scalelib

  implicit none
  public

!  integer,parameter :: mpibufsize=1000000
  integer,save :: nij1
  integer,save :: nij1max
  integer,allocatable,save :: nij1node(:)
!!!!  real(r_size),allocatable,save :: phi1(:)
!!!  real(r_size),allocatable,save :: lon1(:),lat1(:)
!!!  real(r_size),allocatable,save :: lonu1(:),latu1(:)
!!!  real(r_size),allocatable,save :: lonv1(:),latv1(:)
!!!  real(r_size),allocatable,save :: ri1(:),rj1(:)
!!!!  real(r_size),allocatable,save :: wg1(:)
  real(r_size),allocatable,save :: rig1(:),rjg1(:)

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
  integer,save :: lastmem_rank_e

  integer,save :: MPI_COMM_e, nprocs_e, myrank_e
  integer,save :: MPI_COMM_a, nprocs_a, myrank_a

!  character(9) scale_filename = 'file.0000'

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


subroutine rank_1d_2d(proc, iproc, jproc)
  use scale_process, only: PRC_2Drank
  implicit none
  integer, intent(in) :: proc
  integer, intent(out) :: iproc, jproc

  iproc = PRC_2Drank(proc,1)
  jproc = PRC_2Drank(proc,2)

  return  
end subroutine rank_1d_2d


subroutine rank_2d_1d(iproc, jproc, proc)
  use scale_process, only: PRC_NUM_X
  implicit none
  integer, intent(in) :: iproc, jproc
  integer, intent(out) :: proc

  proc = jproc * PRC_NUM_X + iproc

  return  
end subroutine rank_2d_1d


subroutine ij_g2l(proc, ig, jg, il, jl)
  use scale_grid_index, only: &
      IMAX, JMAX
  implicit none
  integer, intent(in) :: proc
  integer, intent(in) :: ig
  integer, intent(in) :: jg
  integer, intent(out) :: il
  integer, intent(out) :: jl
  integer :: iproc, jproc

  call rank_1d_2d(proc, iproc, jproc)
  il = ig - iproc * IMAX
  jl = jg - jproc * JMAX

  return  
end subroutine ij_g2l


subroutine ij_l2g(proc, il, jl, ig, jg)
  use scale_grid_index, only: &
      IMAX, JMAX
  implicit none
  integer, intent(in) :: proc
  integer, intent(in) :: il
  integer, intent(in) :: jl
  integer, intent(out) :: ig
  integer, intent(out) :: jg
  integer :: iproc, jproc

  call rank_1d_2d(proc, iproc, jproc)
  ig = il + iproc * IMAX
  jg = jl + jproc * JMAX

  return  
end subroutine ij_l2g



!-----------------------------------------------------------------------
! using halo!
! proc = -1: outside the global domain
!-----------------------------------------------------------------------
SUBROUTINE rij_g2l_auto(proc,ig,jg,il,jl)
  use scale_grid_index, only: &
      IMAX,JMAX, &
      IHALO,JHALO
!      IA,JA                  ! [for validation]
  use scale_process, only: &
      PRC_NUM_X,PRC_NUM_Y
!      PRC_myrank             ! [for validation]
!  use scale_grid, only: &    ! [for validation]
!      GRID_CX, &             ! [for validation]
!      GRID_CY, &             ! [for validation]
!      GRID_CXG, &            ! [for validation]
!      GRID_CYG, &            ! [for validation]
!      DX, &                  ! [for validation]
!      DY                     ! [for validation]
  IMPLICIT NONE
  integer,INTENT(OUT) :: proc
  REAL(r_size),INTENT(IN) :: ig
  REAL(r_size),INTENT(IN) :: jg
  REAL(r_size),INTENT(OUT) :: il
  REAL(r_size),INTENT(OUT) :: jl
  integer :: iproc, jproc

  if (ig < real(1+IHALO,r_size) .or. ig > real(IMAX*PRC_NUM_X+IHALO,r_size) .or. &
      jg < real(1+JHALO,r_size) .or. jg > real(JMAX*PRC_NUM_Y+JHALO,r_size)) then
    il = -1.0d0
    jl = -1.0d0
    proc = -1
    return
  end if

  iproc = ceiling((ig-real(IHALO,r_size)-0.5d0) / real(IMAX,r_size)) - 1
  jproc = ceiling((jg-real(JHALO,r_size)-0.5d0) / real(JMAX,r_size)) - 1
  il = ig - iproc * IMAX
  jl = jg - jproc * JMAX
  call rank_2d_1d(iproc,jproc,proc)

!  if (PRC_myrank == proc) then                                                                                    ! [for validation]
!    if (rig < (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0 .or. &                                                      ! [for validation]
!        rig > (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0 .or. &                                                     ! [for validation]
!        rjg < (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0 .or. &                                                      ! [for validation]
!        rjg > (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0) then                                                      ! [for validation]
!      write (6,'(A)') 'Error: Process assignment fails!'                                                          ! [for validation]
!      write (6,'(3F10.2)') rig, (GRID_CX(1) - GRID_CXG(1)) / DX + 1.0d0, (GRID_CX(IA) - GRID_CXG(1)) / DX + 1.0d0 ! [for validation]
!      write (6,'(3F10.2)') rjg, (GRID_CY(1) - GRID_CYG(1)) / DY + 1.0d0, (GRID_CY(JA) - GRID_CYG(1)) / DY + 1.0d0 ! [for validation]
!      stop                                                                                                        ! [for validation]
!    end if                                                                                                        ! [for validation]
!  end if                                                                                                          ! [for validation]

  RETURN
END SUBROUTINE rij_g2l_auto

!-----------------------------------------------------------------------
! set_common_mpi_scale
!-----------------------------------------------------------------------
SUBROUTINE set_common_mpi_scale(mem,nnodes,ppn,mem_nodes,mem_np)
  use scale_grid_index, only: &
      IHALO,JHALO, &
      IMAX, JMAX
  use scale_process, only: &
    PRC_myrank

  implicit none
  INTEGER,INTENT(IN) :: mem
  INTEGER,INTENT(IN) :: nnodes,ppn,mem_nodes,mem_np
  REAL(RP) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
  REAL(RP) :: v2dg(nlonsub,nlatsub,nv2d)
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:)
  INTEGER :: i,j,n
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
  integer,allocatable :: ranks_a(:)

  integer :: ip

  integer :: iproc, jproc



  WRITE(6,'(A)') 'Hello from set_common_mpi_scale'


  CALL set_mem_node_proc(mem+1,NNODES,PPN,MEM_NODES,MEM_NP)


  if (scale_IO_group_n >= 1) then

!!!!!!------
    call set_scalelib(MEM_NP, nitmax, nprocs, proc2mem)


    IF(MEM_NODES > 1) THEN
      n_mem = NNODES / MEM_NODES
      n_mempn = 1
    ELSE
      n_mem = NNODES
      n_mempn = PPN / MEM_NP
    END IF
    nprocs_e = n_mem*n_mempn
    nprocs_a = nprocs_e*MEM_NP

    allocate (ranks(nprocs_e))
    allocate (ranks_a(nprocs_a))

    call MPI_Comm_group(MPI_COMM_WORLD,MPI_G_WORLD,ierr)

    do ip = 1, nprocs
      if (proc2mem(1,1,ip) >= 1) then
        if (proc2mem(2,1,ip) == proc2mem(2,1,myrank+1)) then
          ranks(proc2mem(1,1,ip)) = ip-1
        end if
        ranks_a((proc2mem(1,1,ip)-1)*MEM_NP+proc2mem(2,1,ip)+1) = ip-1
      end if
    end do

!write(6,'(A,7I6)') '######===', myrank, ranks(:)

    call MPI_Group_incl(MPI_G_WORLD,nprocs_e,ranks,MPI_G,ierr)
    call MPI_Comm_create(MPI_COMM_WORLD,MPI_G,MPI_COMM_e,ierr)

    call MPI_Comm_size(MPI_COMM_e,nprocs_e,ierr)
    call MPI_Comm_rank(MPI_COMM_e,myrank_e,ierr)

!--

    call MPI_Group_incl(MPI_G_WORLD,nprocs_e*MEM_NP,ranks_a,MPI_G,ierr)
    call MPI_Comm_create(MPI_COMM_WORLD,MPI_G,MPI_COMM_a,ierr)

    call MPI_Comm_size(MPI_COMM_a,nprocs_a,ierr)
    call MPI_Comm_rank(MPI_COMM_a,myrank_a,ierr)


!write(6,'(A,9I6)') '######===', myrank, myrank_e, nprocs_e, ranks(:)
!stop

    deallocate(ranks)
!!!!!!------




    i = MOD(nlonsub*nlatsub,nprocs_e)
    nij1max = (nlonsub*nlatsub - i)/nprocs_e + 1
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
    ALLOCATE(rig1(nij1))
    ALLOCATE(rjg1(nij1))

    ALLOCATE(v3d(nij1,nlev,nv3d))
    ALLOCATE(v2d(nij1,nv2d))


    call rank_1d_2d(PRC_myrank, iproc, jproc)
    do j = 1, nlatsub
      do i = 1, nlonsub
        v3dg(1,i,j,1) = i + iproc * IMAX + IHALO
        v3dg(1,i,j,2) = j + jproc * JMAX + JHALO
      end do
    end do

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

!!!!!! ----- need to be replaced by more native communication!!!!
    v3dg = 0.0
    v2dg = 0.0
!!!!!!
    CALL scatter_grd_mpi(0,v3dg,v2dg,v3d,v2d)

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

    rig1   = v3d(:,1,1)
    rjg1   = v3d(:,1,2)


  end if ! [ scale_IO_group_n >= 1 ]


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
    call MPI_Comm_free(MPI_COMM_a,ierr)
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

  lastmem_rank_e = mod(mem-1, n_mem*n_mempn)
!print *, lastmem_rank_e

  if (lastmem_rank_e /= proc2mem(1,1,mem2proc(1,mem)+1)-1) then
    print *, 'XXXXXX wrong!!'
    stop
  end if


  RETURN
END SUBROUTINE
!-----------------------------------------------------------------------
! Scatter gridded data to processes (nrank -> all)
!-----------------------------------------------------------------------
!SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
!  INTEGER,INTENT(IN) :: nrank
!  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
!  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)

!!  IF(mpibufsize > nij1max) THEN
!    CALL scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
!!  ELSE
!!    CALL scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
!!  END IF

!  RETURN
!END SUBROUTINE scatter_grd_mpi

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

!SUBROUTINE scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
!  INTEGER,INTENT(IN) :: nrank
!  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
!  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
!  REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
!  REAL(r_sngl) :: bufr(nij1max,nlevall)
!  INTEGER :: j,k,n,ierr,ns,nr

!  ns = nij1max * nlevall
!  nr = ns
!  IF(myrank == nrank) THEN
!    j=0
!    DO n=1,nv3d
!      DO k=1,nlev
!        j = j+1
!        CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
!      END DO
!    END DO

!    DO n=1,nv2d
!      j = j+1
!      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
!    END DO
!  END IF

!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
!                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

!  j=0
!  DO n=1,nv3d
!    DO k=1,nlev
!      j = j+1
!      v3d(:,k,n) = REAL(bufr(1:nij1,j),r_size)
!    END DO
!  END DO

!  DO n=1,nv2d
!    j = j+1
!    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
!  END DO

!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!  RETURN
!END SUBROUTINE scatter_grd_mpi_fast

SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlonsub,nlatsub,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(RP) :: bufs(nij1max,nlevall,nprocs_e)
  REAL(RP) :: bufr(nij1max,nlevall)
  INTEGER :: j,k,n,ierr,ns,nr

  integer :: MPI_RP

!!!!!!
  IF(RP == r_dble) THEN
    MPI_RP = MPI_DOUBLE_PRECISION
  ELSE IF(RP == r_sngl) THEN
    MPI_RP = MPI_REAL
  END IF
!!!!!!

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

  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_RP,&
                 & bufr,nr,MPI_RP,nrank,MPI_COMM_e,ierr)

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

  CALL MPI_BARRIER(MPI_COMM_e,ierr)

  RETURN
END SUBROUTINE scatter_grd_mpi
!-----------------------------------------------------------------------
! Gather gridded data (all -> nrank)
!-----------------------------------------------------------------------
SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlonsub,nlatsub,nv2d)
  REAL(RP) :: bufs(nij1max,nlevall)
  REAL(RP) :: bufr(nij1max,nlevall,nprocs_e)
  INTEGER :: j,k,n,ierr,ns,nr

  integer :: MPI_RP

!!!!!!
  IF(RP == r_dble) THEN
    MPI_RP = MPI_DOUBLE_PRECISION
  ELSE IF(RP == r_sngl) THEN
    MPI_RP = MPI_REAL
  END IF
!!!!!!

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

  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_GATHER(bufs,ns,MPI_RP,&
                & bufr,nr,MPI_RP,nrank,MPI_COMM_e,ierr)

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

  CALL MPI_BARRIER(MPI_COMM_e,ierr)

  RETURN
END SUBROUTINE gather_grd_mpi

!!!!!!!!!!subroutine set_scale_mpi_comm
!!!!!!!!!!  implicit none

!!!!!!!!!!!    real(RP), allocatable :: U(:,:,:), MOMX(:,:,:)

!!!!!!!!!!  integer :: k, i, j

!!!!!!!!!!  integer               :: dim1_max, dim1_S, dim1_E
!!!!!!!!!!  integer               :: dim2_max, dim2_S, dim2_E
!!!!!!!!!!  integer               :: dim3_max, dim3_S, dim3_E
!!!!!!!!!!  integer               :: dim4_max, dim4_S, dim4_E
!!!!!!!!!!  real(RP), allocatable :: var3D(:,:,:)

!!!!!!!!!!  integer :: iolen

!!!!!!!!!!  character(len=100) :: basename
!!!!!!!!!!  character(len=100) :: varname
!!!!!!!!!!  integer :: step
!!!!!!!!!!  basename = 'history'
!!!!!!!!!!  varname = 'U'
!!!!!!!!!!  step = 1

!!!!!!!!!!  if (PRC_nu == 0) then
!!!!!!!!!!    basename = trim(basename) // '.u000000'
!!!!!!!!!!  else if (PRC_nu == 1) THEN
!!!!!!!!!!    basename = trim(basename) // '.u000001'
!!!!!!!!!!  end if

!!!!!!!!!!  allocate( U(KA,IA,JA) )
!!!!!!!!!!  allocate( MOMX(KA,IA,JA) )
!!!!!!!!!!  U = 0.0d0




!!!!!!!!!!    ! Read file

!!!!!!!!!!!    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(varname)

!!!!!!!!!!       dim1_max = IMAX !KMAX
!!!!!!!!!!       dim2_max = JMAX !IMAX
!!!!!!!!!!       dim3_max = KMAX !JMAX
!!!!!!!!!!       dim1_S   = IS !KS
!!!!!!!!!!       dim1_E   = IE !KE
!!!!!!!!!!       dim2_S   = JS !IS
!!!!!!!!!!       dim2_E   = JE !IE
!!!!!!!!!!       dim3_S   = KS !JS
!!!!!!!!!!       dim3_E   = KE !JE

!!!!!!!!!!    allocate( var3D(dim1_max,dim2_max,dim3_max) )


!!!!!!!!!!    call HIST_get(var3D, trim(basename), trim(varname), step=step)


!!!!!!!!!!    forall (i=1:IMAX, j=1:JMAX, k=1:KMAX) U(k+KHALO,i+IHALO,j+JHALO) = var3D(i,j,k)

!!!!!!!!!!    deallocate( var3D )


!!!!!!!!!!    if (PRC_nu == 0) then
!!!!!!!!!!      call FILEIO_read( MOMX(:,:,:),                          & ! [OUT]
!!!!!!!!!!                        'init.u000000', 'MOMX', 'ZXY', step=1 ) ! [IN]
!!!!!!!!!!    else if (PRC_nu == 1) then
!!!!!!!!!!      call FILEIO_read( MOMX(:,:,:),                          & ! [OUT]
!!!!!!!!!!                        'init.u000001', 'MOMX', 'ZXY', step=1 ) ! [IN]
!!!!!!!!!!    end if


!!!!!!!!!!    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
!!!!!!!!!!    do j  = JS, JE
!!!!!!!!!!    do i  = IS, IE
!!!!!!!!!!       U(   1:KS-1,i,j) = U(KS,i,j)
!!!!!!!!!!       U(KE+1:KA,  i,j) = U(KE,i,j)
!!!!!!!!!!       MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
!!!!!!!!!!       MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
!!!!!!!!!!    enddo
!!!!!!!!!!    enddo

!!!!!!!!!!    call COMM_vars8( U   (:,:,:), 1 )
!!!!!!!!!!    call COMM_vars8( MOMX(:,:,:), 2 )
!!!!!!!!!!    call COMM_wait ( U   (:,:,:), 1 )
!!!!!!!!!!    call COMM_wait ( MOMX(:,:,:), 2 )


SUBROUTINE read_ens_history_mpi(file,iter,step,v3dg,v2dg,ensmean)
  use scale_grid_index, only: &
      IMAX, JMAX, KMAX, &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA
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

  nbvr = nbv
  if (present(ensmean)) then
    if (ensmean) then
      nbvr = nbv+1
    end if
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
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!SUBROUTINE read_ens_restart_mpi(file,iter,v3dg,v2dg,ensmean)
!  use gtool_file, only: &
!    FileRead
!  use scale_process, only: &
!    PRC_myrank

!  IMPLICIT NONE

!  CHARACTER(4),INTENT(IN) :: file
!  INTEGER,INTENT(IN) :: iter
!  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
!  REAL(RP),INTENT(OUT) :: v2dg(nlonsub,nlatsub,nv2d)
!  LOGICAL,INTENT(INOUT),OPTIONAL :: ensmean
!  INTEGER :: i,j,k,iv3d,iv2d,nbvr

!  CHARACTER(9) :: filename='file.0000'

!  if (.not. present(ensmean) .or. (.not. ensmean)) then
!    nbvr = nbv
!  else
!    nbvr = nbv+1
!  end if

!  IF (valid_member) then
!    IF(proc2mem(1,iter,myrank+1) >= 1 .and. proc2mem(1,iter,myrank+1) <= nbvr) THEN
!      WRITE(filename(1:4),'(A4)') file

!      if (proc2mem(1,iter,myrank+1) == nbv+1) then
!        WRITE(filename(6:9),'(A4)') 'mean'
!      else
!        WRITE(filename(6:9),'(I4.4)') proc2mem(1,iter,myrank+1)
!      end if

!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',proc2mem(2,iter,myrank+1),'.nc'

!      call read_restart(filename,v3dg,v2dg)
!!      DO iv3d = 1, nv3d
!!        if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
!!        select case (iv3d)
!!        case (iv3d_u)
!!          call FileRead(v3dg(:,:,:,iv3d), filename, 'MOMX', 1, PRC_myrank)
!!          v3dg(:,:,:,iv3d) = v3dg(:,:,:,iv3d) / v3dg(:,:,:,iv3d_rho)
!!        case (iv3d_v)
!!          call FileRead(v3dg(:,:,:,iv3d), filename, 'MOMY', 1, PRC_myrank)
!!          v3dg(:,:,:,iv3d) = v3dg(:,:,:,iv3d) / v3dg(:,:,:,iv3d_rho)
!!        case (iv3d_w)
!!          call FileRead(v3dg(:,:,:,iv3d), filename, 'MOMZ', 1, PRC_myrank)
!!          v3dg(:,:,:,iv3d) = v3dg(:,:,:,iv3d) / v3dg(:,:,:,iv3d_rho)
!!        case default
!!          call FileRead(v3dg(:,:,:,iv3d), filename, trim(v3d_name(iv3d)), 1, PRC_myrank)
!!        end select
!!      END DO

!!      DO iv2d = 1, nv2d
!!        if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
!!        call FileRead(v2dg(:,:,iv2d), filename, trim(v2d_name(iv2d)), 1, PRC_myrank)
!!      END DO
!    END IF
!  END IF

!  RETURN
!END SUBROUTINE read_ens_restart_mpi
!-----------------------------------------------------------------------
! Read ensemble data and distribute to processes
!-----------------------------------------------------------------------
subroutine read_ens_mpi(file,v3d,v2d)
  implicit none
  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nbv,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nbv,nv2d)
  REAL(RP) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
  REAL(RP) :: v2dg(nlonsub,nlatsub,nv2d)
  CHARACTER(9) :: filename='file.0000'
  integer :: it,im,mstart,mend

  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    if (im >= 1 .and. im <= nbv) then
      WRITE(filename(1:4),'(A4)') file
      WRITE(filename(6:9),'(I4.4)') im
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
      call read_restart(filename,v3dg,v2dg)
      call state_trans(v3dg)
    end if
    mstart = 1 + (it-1)*nprocs_e
    mend = MIN(it*nprocs_e, nbv)
    CALL scatter_grd_mpi_alltoall(mstart,mend,nbv,v3dg,v2dg,v3d,v2d)
  end do ! [ it = 1, nitmax ]

  return
end subroutine read_ens_mpi


!-----------------------------------------------------------------------
! Write ensemble data after collecting data from processes
!-----------------------------------------------------------------------
SUBROUTINE write_ens_mpi(file,v3d,v2d)
  implicit none
  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nbv,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nbv,nv2d)
  REAL(RP) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
  REAL(RP) :: v2dg(nlonsub,nlatsub,nv2d)
  CHARACTER(9) :: filename='file.0000'
  integer :: it,im,mstart,mend

  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    mstart = 1 + (it-1)*nprocs_e
    mend = MIN(it*nprocs_e, nbv)
    CALL gather_grd_mpi_alltoall(mstart,mend,nbv,v3d,v2d,v3dg,v2dg)
    if (im >= 1 .and. im <= nbv) then
      WRITE(filename(1:4),'(A4)') file
      WRITE(filename(6:9),'(I4.4)') im
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
      call state_trans_inv(v3dg)
      call write_restart(filename,v3dg,v2dg)
    end if
  end do ! [ it = 1, nitmax ]

  return
END SUBROUTINE write_ens_mpi

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!SUBROUTINE write_ens_restart_mpi(file,iter,v3dg,v2dg,ensmean)
!  use gtool_file, only: &
!    FileWrite
!  use scale_process, only: &
!    PRC_myrank

!  IMPLICIT NONE

!  CHARACTER(4),INTENT(IN) :: file
!  INTEGER,INTENT(IN) :: iter
!  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
!  REAL(RP),INTENT(OUT) :: v2dg(nlonsub,nlatsub,nv2d)
!  LOGICAL,INTENT(INOUT),OPTIONAL :: ensmean
!  INTEGER :: i,j,k,iv3d,iv2d,nbvr

!  CHARACTER(9) :: filename='file.0000'

!  if (.not. present(ensmean) .or. (.not. ensmean)) then
!    nbvr = nbv
!  else
!    nbvr = nbv+1
!  end if

!  IF (valid_member) then
!    IF(proc2mem(1,iter,myrank+1) >= 1 .and. proc2mem(1,iter,myrank+1) <= nbvr) THEN
!      WRITE(filename(1:4),'(A4)') file

!      if (proc2mem(1,iter,myrank+1) == nbv+1) then
!        WRITE(filename(6:9),'(A4)') 'mean'
!      else
!        WRITE(filename(6:9),'(I4.4)') proc2mem(1,iter,myrank+1)
!      end if

!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',proc2mem(2,iter,myrank+1),'.nc'

!      DO iv3d = 1, nv3d
!        if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
!        select case (iv3d)
!        case (iv3d_u)
!          call FileWrite(v3dg(:,:,:,iv3d), filename, 'MOMX', 1, PRC_myrank)
!          v3dg(:,:,:,iv3d) = v3dg(:,:,:,iv3d) / v3dg(:,:,:,iv3d_rho)
!        case (iv3d_v)
!          call FileWrite(v3dg(:,:,:,iv3d), filename, 'MOMY', 1, PRC_myrank)
!          v3dg(:,:,:,iv3d) = v3dg(:,:,:,iv3d) / v3dg(:,:,:,iv3d_rho)
!        case (iv3d_w)
!          call FileWrite(v3dg(:,:,:,iv3d), filename, 'MOMZ', 1, PRC_myrank)
!          v3dg(:,:,:,iv3d) = v3dg(:,:,:,iv3d) / v3dg(:,:,:,iv3d_rho)
!        case default
!          call FileWrite(v3dg(:,:,:,iv3d), filename, trim(v3d_name(iv3d)), 1, PRC_myrank)
!        end select
!      END DO

!      DO iv2d = 1, nv2d
!        if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
!        call FileWrite(v2dg(:,:,iv2d), filename, trim(v2d_name(iv2d)), 1, PRC_myrank)
!      END DO
!    END IF
!  END IF

!  RETURN
!END SUBROUTINE write_ens_restart_mpi



SUBROUTINE scatter_grd_mpi_alltoall(mstart,mend,mem,v3dg,v2dg,v3d,v2d)!,ngp,ngpmax,ngpnode)
  INTEGER,INTENT(IN) :: mstart,mend,mem !,ngp,ngpmax,ngpnode(nprocs_e)
!  INTEGER,INTENT(IN) :: np
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlonsub,nlatsub,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,mem,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,mem,nv2d)
  REAL(RP) :: bufs(nij1max,nlevall,nprocs_e)
  REAL(RP) :: bufr(nij1max,nlevall,nprocs_e)
!  REAL(r_sngl),ALLOCATABLE :: bufs3(:,:,:) , bufr3(:,:,:)
!  REAL(r_sngl),ALLOCATABLE :: bufs2(:,:)   , bufr2(:,:)
  INTEGER :: k,n,j,m,mcount,ierr
  INTEGER :: ns(nprocs_e),nst(nprocs_e),nr(nprocs_e),nrt(nprocs_e)

  integer :: MPI_RP

!!!!!!
  IF(RP == r_dble) THEN
    MPI_RP = MPI_DOUBLE_PRECISION
  ELSE IF(RP == r_sngl) THEN
    MPI_RP = MPI_REAL
  END IF
!!!!!!

  mcount = mend - mstart + 1
  IF(mcount > nprocs_e .OR. mcount <= 0) STOP

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

  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  IF(mcount == nprocs_e) THEN
    CALL MPI_ALLTOALL(bufs, nij1max*nlevall, MPI_RP, &
                      bufr, nij1max*nlevall, MPI_RP, MPI_COMM_e, ierr)
  ELSE
    CALL set_alltoallv_counts(mcount,nij1max*nlevall,nprocs_e,nr,nrt,ns,nst)
    CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_RP, &
                       bufr, nr, nrt, MPI_RP, MPI_COMM_e, ierr)
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

  CALL MPI_BARRIER(MPI_COMM_e,ierr)
!  DEALLOCATE(bufr,bufs)
  RETURN
END SUBROUTINE scatter_grd_mpi_alltoall

!-----------------------------------------------------------------------
! Gather gridded data using MPI_ALLTOALL(V) (all -> mstart~mend)
!-----------------------------------------------------------------------

SUBROUTINE gather_grd_mpi_alltoall(mstart,mend,mem,v3d,v2d,v3dg,v2dg) !,ngp,ngpmax,ngpnode)
  INTEGER,INTENT(IN) :: mstart,mend,mem !,ngpmax,ngpnode(nprocs_e)
!  INTEGER,INTENT(IN) :: np
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,mem,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,mem,nv2d)
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlonsub,nlatsub,nv2d)
  REAL(RP) :: bufs(nij1max,nlevall,nprocs_e)
  REAL(RP) :: bufr(nij1max,nlevall,nprocs_e)
!  REAL(r_sngl),ALLOCATABLE :: bufs3(:,:,:) , bufr3(:,:,:)
!  REAL(r_sngl),ALLOCATABLE :: bufs2(:,:)   , bufr2(:,:)
  INTEGER :: k,n,j,m,mcount,ierr
  INTEGER :: ns(nprocs_e),nst(nprocs_e),nr(nprocs_e),nrt(nprocs_e)

  integer :: MPI_RP

!!!!!!
  IF(RP == r_dble) THEN
    MPI_RP = MPI_DOUBLE_PRECISION
  ELSE IF(RP == r_sngl) THEN
    MPI_RP = MPI_REAL
  END IF
!!!!!!

  mcount = mend - mstart + 1
  IF(mcount > nprocs_e .OR. mcount <= 0) STOP

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

  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  IF(mcount == nprocs_e) THEN
    CALL MPI_ALLTOALL(bufs, nij1max*nlevall, MPI_RP, &
                      bufr, nij1max*nlevall, MPI_RP, MPI_COMM_e, ierr)
  ELSE
    CALL set_alltoallv_counts(mcount,nij1max*nlevall,nprocs_e,ns,nst,nr,nrt)
    CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_RP, &
                       bufr, nr, nrt, MPI_RP, MPI_COMM_e, ierr)
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

  CALL MPI_BARRIER(MPI_COMM_e,ierr)
!  DEALLOCATE(bufr,bufs)
  RETURN
END SUBROUTINE gather_grd_mpi_alltoall

!-----------------------------------------------------------------------
! Set the send/recieve counts of MPI_ALLTOALLV
!-----------------------------------------------------------------------
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


!-----------------------------------------------------------------------
! gridded data -> buffer
!-----------------------------------------------------------------------
SUBROUTINE grd_to_buf(np,grd,buf)
  INTEGER,INTENT(IN) :: np
  REAL(RP),INTENT(IN) :: grd(nlonsub,nlatsub)
  REAL(RP),INTENT(OUT) :: buf(nij1max,np)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,np
    DO i=1,nij1node(m)
      j = m-1 + np * (i-1)
      ilon = MOD(j,nlonsub) + 1
      ilat = (j-ilon+1) / nlonsub + 1
!if (i < 1 .or. i > nij1max .or. m < 1 .or. m > np .or. ilon < 1 .or. ilon > nlonsub .or. ilat < 1 .or. ilat > nlatsub) then
!print *, '######', np, nij1max
!print *, '########', i, m, ilon, ilat
!stop
!end if
      buf(i,m) = grd(ilon,ilat)
    END DO
  END DO

  DO m=1,np
    IF(nij1node(m) < nij1max) buf(nij1max,m) = undef
  END DO

  RETURN
END SUBROUTINE grd_to_buf
!-----------------------------------------------------------------------
! buffer -> gridded data
!-----------------------------------------------------------------------
SUBROUTINE buf_to_grd(np,buf,grd)
  INTEGER,INTENT(IN) :: np
  REAL(RP),INTENT(IN) :: buf(nij1max,np)
  REAL(RP),INTENT(OUT) :: grd(nlonsub,nlatsub)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,np
    DO i=1,nij1node(m)
      j = m-1 + np * (i-1)
      ilon = MOD(j,nlonsub) + 1
      ilat = (j-ilon+1) / nlonsub + 1
      grd(ilon,ilat) = buf(i,m)
    END DO
  END DO

  RETURN
END SUBROUTINE buf_to_grd
!-----------------------------------------------------------------------
! STORING DATA (ensemble mean and spread)
!-----------------------------------------------------------------------
SUBROUTINE write_ensmspr_mpi(file,v3d,v2d)
  use scale_process, only: PRC_myrank
  implicit none

  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nbv,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nbv,nv2d)
  REAL(r_size) :: v3dm(nij1,nlev,nv3d)
  REAL(r_size) :: v2dm(nij1,nv2d)
  REAL(r_size) :: v3ds(nij1,nlev,nv3d)
  REAL(r_size) :: v2ds(nij1,nv2d)
  REAL(RP) :: v3dg(nlev,nlonsub,nlatsub,nv3d)
  REAL(RP) :: v2dg(nlonsub,nlatsub,nv2d)
  INTEGER :: i,k,m,n
  CHARACTER(9) :: filename='file.0000'


  REAL(r_size) :: timer
  INTEGER :: ierr

  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  CALL CPU_TIME(timer)
  if (myrank == 0) print *, '######', timer


  CALL ensmean_grd(nbv,nij1,v3d,v2d,v3dm,v2dm)


  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  CALL CPU_TIME(timer)
  if (myrank == 0) print *, '######', timer


  CALL gather_grd_mpi(lastmem_rank_e,v3dm,v2dm,v3dg,v2dg)


  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  CALL CPU_TIME(timer)
  if (myrank == 0) print *, '######', timer


  IF(myrank_e == lastmem_rank_e) THEN
    WRITE(filename(1:4),'(A4)') file
    WRITE(filename(6:9),'(A4)') 'mean'
!    WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
    call state_trans_inv(v3dg)
    call write_restart(filename,v3dg,v2dg)
  END IF


  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  CALL CPU_TIME(timer)
  if (myrank == 0) print *, '######', timer


  DO n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij1
        v3ds(i,k,n) = (v3d(i,k,1,n)-v3dm(i,k,n))**2
        DO m=2,nbv
          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n)-v3dm(i,k,n))**2
        END DO
        v3ds(i,k,n) = SQRT(v3ds(i,k,n) / REAL(nbv-1,r_size))
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO


  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  CALL CPU_TIME(timer)
  if (myrank == 0) print *, '######', timer


  DO n=1,nv2d
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO i=1,nij1
      v2ds(i,n) = (v2d(i,1,n)-v2dm(i,n))**2
      DO m=2,nbv
        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n)-v2dm(i,n))**2
      END DO
      v2ds(i,n) = SQRT(v2ds(i,n) / REAL(nbv-1,r_size))
    END DO
!$OMP END PARALLEL DO
  END DO


  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  CALL CPU_TIME(timer)
  if (myrank == 0) print *, '######', timer


  CALL gather_grd_mpi(lastmem_rank_e,v3ds,v2ds,v3dg,v2dg)


  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  CALL CPU_TIME(timer)
  if (myrank == 0) print *, '######', timer


  IF(myrank_e == lastmem_rank_e) THEN
    WRITE(filename(1:4),'(A4)') file
    WRITE(filename(6:9),'(A4)') 'sprd'
!    WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
!    call state_trans_inv(v3dg)             !!
    call write_restart(filename,v3dg,v2dg)  !! not transformed to rho,rhou,rhov,rhow,rhot before writing.
  END IF


  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  CALL CPU_TIME(timer)
  if (myrank == 0) print *, '######', timer


  RETURN
END SUBROUTINE write_ensmspr_mpi
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
