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
!   01/08/2015 Guo-Yuan Lien     modified for SCALE model
!
!=======================================================================
!$USE OMP_LIB
  use common

  use common_nml

  use common_mpi
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
  real(r_size),allocatable,save :: topo1(:)
  real(r_size),allocatable,save :: hgt1(:,:)

  integer,save :: MPI_COMM_e, nprocs_e, myrank_e
  integer,save :: MPI_COMM_a, nprocs_a, myrank_a

!  character(9) scale_filename = 'file.0000'

contains
!-----------------------------------------------------------------------
! set_common_mpi_scale
!-----------------------------------------------------------------------
SUBROUTINE set_common_mpi_scale
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
  INTEGER :: i,j,n
  INTEGER :: ierr !,buf(4)
!  LOGICAL :: ex


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


!  CALL set_mem_node_proc(mem+1,NNODES,PPN,MEM_NODES,MEM_NP)



!!!!!!------
!    call set_scalelib(MEM_NP, nitmax, nprocs, proc2mem)


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
  ALLOCATE(topo1(nij1))

  ALLOCATE(hgt1(nij1,nlev))

  ALLOCATE(v3d(nij1,nlev,nv3d))
  ALLOCATE(v2d(nij1,nv2d))


!!!!!! ----- need to be replaced by more native communication!!!!
  v3dg = 0.0d0
  v2dg = 0.0d0

  call rank_1d_2d(PRC_myrank, iproc, jproc)
  do j = 1, nlat
    do i = 1, nlon
      v3dg(1,i,j,1) = real(i + iproc * nlon + IHALO, RP)
      v3dg(1,i,j,2) = real(j + jproc * nlat + JHALO, RP)
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

  if (myrank_e == lastmem_rank_e) then
!    call read_topo('topo', v2dg(:,:,1))  !!! nv2d = 0 in scale...
    call read_topo('topo', v3dg(1,:,:,3))

!print *, v3dg(1,:,:,3)

  end if

  CALL scatter_grd_mpi(lastmem_rank_e,v3dg,v2dg,v3d,v2d)
!!!!!! -----

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

!  topo1  = v2d(:,1)
  topo1  = v3d(:,1,3)

  call scale_calc_z(nij1, topo1, hgt1)

!!print *, hgt1(1,:)

!do n = 1, nij1
!if (hgt1(n,1) > 120.0d0) then
!print *, hgt1(n,:)
!end if
!end do

  RETURN
END SUBROUTINE set_common_mpi_scale
!-----------------------------------------------------------------------
! set_common_mpi_scale
!-----------------------------------------------------------------------
SUBROUTINE unset_common_mpi_scale
  implicit none
  integer:: ierr

  call MPI_Comm_free(MPI_COMM_e,ierr)
  call MPI_Comm_free(MPI_COMM_a,ierr)
!  call unset_scalelib

  RETURN
END SUBROUTINE unset_common_mpi_scale
!-----------------------------------------------------------------------
! Scatter gridded data to processes (nrank -> all)
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
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
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
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



SUBROUTINE read_ens_history_iter(file,iter,step,v3dg,v2dg,ensmean)
  IMPLICIT NONE

  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: iter
  INTEGER,INTENT(IN) :: step
  REAL(r_size),INTENT(OUT) :: v3dg(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(OUT) :: v2dg(nlonh,nlath,nv2dd)
  LOGICAL,INTENT(INOUT),OPTIONAL :: ensmean
  CHARACTER(9) :: filename='file.0000'

  integer :: mem

  mem = MEMBER
  if (present(ensmean)) then
    if (ensmean) then
      mem = MEMBER + 1
    end if
  end if

  IF (scale_IO_use) then

    IF(proc2mem(1,iter,myrank+1) >= 1 .and. proc2mem(1,iter,myrank+1) <= mem) THEN
      WRITE(filename(1:4),'(A4)') file

      if (proc2mem(1,iter,myrank+1) == MEMBER+1) then
        WRITE(filename(6:9),'(A4)') 'mean'
      else
        WRITE(filename(6:9),'(I4.4)') proc2mem(1,iter,myrank+1)
      end if

      call read_history(filename,step,v3dg,v2dg)
    END IF

  END IF

  RETURN
END SUBROUTINE read_ens_history_iter




!SUBROUTINE read_ens_history_mpi(file,iter,step,v3dg,v2dg,ensmean)
!  use scale_grid_index, only: &
!      IHALO, JHALO, KHALO, &
!      IS, IE, JS, JE, KS, KE, KA
!  use scale_history, only: &
!     HIST_get
!  use scale_comm, only: &
!     COMM_vars8, &
!     COMM_wait

!  IMPLICIT NONE

!  CHARACTER(4),INTENT(IN) :: file
!  INTEGER,INTENT(IN) :: iter
!  INTEGER,INTENT(IN) :: step
!  REAL(r_size),INTENT(OUT) :: v3dg(nlevh,nlonh,nlath,nv3dd)
!  REAL(r_size),INTENT(OUT) :: v2dg(nlonh,nlath,nv2dd)
!  LOGICAL,INTENT(INOUT),OPTIONAL :: ensmean
!  INTEGER :: i,j,k,iv3d,iv2d
!  CHARACTER(9) :: filename='file.0000'

!  real(RP), allocatable :: var3D(:,:,:)
!  real(RP), allocatable :: var2D(:,:)

!  integer :: mem

!  mem = MEMBER
!  if (present(ensmean)) then
!    if (ensmean) then
!      mem = MEMBER + 1
!    end if
!  end if

!  IF (scale_IO_use) then
!    allocate( var3D(nlon,nlat,nlev) )
!    allocate( var2D(nlon,nlat) )

!    IF(proc2mem(1,iter,myrank+1) >= 1 .and. proc2mem(1,iter,myrank+1) <= mem) THEN
!      WRITE(filename(1:4),'(A4)') file

!      if (proc2mem(1,iter,myrank+1) == MEMBER+1) then
!        WRITE(filename(6:9),'(A4)') 'mean'
!      else
!        WRITE(filename(6:9),'(I4.4)') proc2mem(1,iter,myrank+1)
!      end if

!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',proc2mem(2,iter,myrank+1),'.nc'

!      DO iv3d = 1, nv3dd
!        if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(v3dd_name(iv3d))
!        call HIST_get(var3D, filename, trim(v3dd_name(iv3d)), step)
!        FORALL (i=1:nlon, j=1:nlat, k=1:nlev) v3dg(k+KHALO,i+IHALO,j+JHALO,iv3d) = var3D(i,j,k)

!!!!!!$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
!        do j  = JS, JE
!          do i  = IS, IE
!            v3dg(   1:KS-1,i,j,iv3d) = v3dg(KS,i,j,iv3d)
!            v3dg(KE+1:KA,  i,j,iv3d) = v3dg(KE,i,j,iv3d)
!          enddo
!        enddo
!      END DO

!      DO iv3d = 1, nv3dd
!        call COMM_vars8( v3dg(:,:,:,iv3d), iv3d )
!      END DO
!      DO iv3d = 1, nv3dd
!        call COMM_wait ( v3dg(:,:,:,iv3d), iv3d )
!      END DO

!      DO iv2d = 1, nv2dd
!        if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 2D var: ', trim(v2dd_name(iv2d))
!        call HIST_get(var2D, filename, trim(v2dd_name(iv2d)), step)
!        v2dg(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2d) = var2D(:,:)
!      END DO

!      DO iv2d = 1, nv2dd
!        call COMM_vars8( v2dg(:,:,iv2d), iv2d )
!      END DO
!      DO iv2d = 1, nv2dd
!        call COMM_wait ( v2dg(:,:,iv2d), iv2d )
!      END DO
!    END IF

!    deallocate( var3D )
!    deallocate( var2D )
!  END IF

!  RETURN
!END SUBROUTINE read_ens_history_mpi
!-----------------------------------------------------------------------
! Read ensemble data and distribute to processes
!-----------------------------------------------------------------------
subroutine read_ens_mpi(file,v3d,v2d)
  implicit none
  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,MEMBER,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,MEMBER,nv2d)
  REAL(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP) :: v2dg(nlon,nlat,nv2d)
  CHARACTER(9) :: filename='file.0000'
  integer :: it,im,mstart,mend


  integer :: ierr
  REAL(r_dble) :: rrtimer00,rrtimer

  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer00 = MPI_WTIME()


  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    if (im >= 1 .and. im <= MEMBER) then
      WRITE(filename(1:4),'(A4)') file
      WRITE(filename(6:9),'(I4.4)') im
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
      call read_restart(filename,v3dg,v2dg)


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,F10.2)') '###### read_ens_mpi:read_restart:              ',rrtimer-rrtimer00
  rrtimer00=rrtimer


      call state_trans(v3dg)


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,F10.2)') '###### read_ens_mpi:state_trans:               ',rrtimer-rrtimer00
  rrtimer00=rrtimer


    end if
    mstart = 1 + (it-1)*nprocs_e
    mend = MIN(it*nprocs_e, MEMBER)
    CALL scatter_grd_mpi_alltoall(mstart,mend,v3dg,v2dg,v3d,v2d)
  end do ! [ it = 1, nitmax ]


  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,F10.2)') '###### read_ens_mpi:scatter_grd_mpi_alltoall:  ',rrtimer-rrtimer00
  rrtimer00=rrtimer


  return
end subroutine read_ens_mpi


!-----------------------------------------------------------------------
! Write ensemble data after collecting data from processes
!-----------------------------------------------------------------------
SUBROUTINE write_ens_mpi(file,v3d,v2d)
  implicit none
  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,MEMBER,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,MEMBER,nv2d)
  REAL(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP) :: v2dg(nlon,nlat,nv2d)
  CHARACTER(9) :: filename='file.0000'
  integer :: it,im,mstart,mend


  integer :: ierr
  REAL(r_dble) :: rrtimer00,rrtimer

  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer00 = MPI_WTIME()


  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    mstart = 1 + (it-1)*nprocs_e
    mend = MIN(it*nprocs_e, MEMBER)
    CALL gather_grd_mpi_alltoall(mstart,mend,v3d,v2d,v3dg,v2dg)


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,F10.2)') '###### write_ens_mpi:gather_grd_mpi_alltoall:  ',rrtimer-rrtimer00
  rrtimer00=rrtimer


    if (im >= 1 .and. im <= MEMBER) then
      WRITE(filename(1:4),'(A4)') file
      WRITE(filename(6:9),'(I4.4)') im
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
      call state_trans_inv(v3dg)


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,F10.2)') '###### write_ens_mpi:state_trans_inv:          ',rrtimer-rrtimer00
  rrtimer00=rrtimer


      call write_restart(filename,v3dg,v2dg)
    end if
  end do ! [ it = 1, nitmax ]


  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,F10.2)') '###### write_ens_mpi:write_restart:            ',rrtimer-rrtimer00
  rrtimer00=rrtimer


  return
END SUBROUTINE write_ens_mpi





SUBROUTINE scatter_grd_mpi_alltoall(mstart,mend,v3dg,v2dg,v3d,v2d)!,ngp,ngpmax,ngpnode)
  INTEGER,INTENT(IN) :: mstart,mend !,ngp,ngpmax,ngpnode(nprocs_e)
!  INTEGER,INTENT(IN) :: np
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,MEMBER,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,MEMBER,nv2d)
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

SUBROUTINE gather_grd_mpi_alltoall(mstart,mend,v3d,v2d,v3dg,v2dg) !,ngp,ngpmax,ngpnode)
  INTEGER,INTENT(IN) :: mstart,mend !,ngpmax,ngpnode(nprocs_e)
!  INTEGER,INTENT(IN) :: np
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,MEMBER,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,MEMBER,nv2d)
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
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
  REAL(RP),INTENT(IN) :: grd(nlon,nlat)
  REAL(RP),INTENT(OUT) :: buf(nij1max,np)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,np
    DO i=1,nij1node(m)
      j = m-1 + np * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
!if (i < 1 .or. i > nij1max .or. m < 1 .or. m > np .or. ilon < 1 .or. ilon > nlon .or. ilat < 1 .or. ilat > nlat) then
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
!-----------------------------------------------------------------------
! STORING DATA (ensemble mean and spread)
!-----------------------------------------------------------------------
SUBROUTINE write_ensmspr_mpi(file,v3d,v2d,obs,obsda)
  use scale_process, only: PRC_myrank



!  use scale_process, only: &
!      PRC_myrank
  use scale_grid_index, only: &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA
  use scale_comm, only: &
      COMM_vars8, &
      COMM_wait



  implicit none

  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,MEMBER,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,MEMBER,nv2d)
  REAL(r_size) :: v3dm(nij1,nlev,nv3d)
  REAL(r_size) :: v2dm(nij1,nv2d)
  REAL(r_size) :: v3ds(nij1,nlev,nv3d)
  REAL(r_size) :: v2ds(nij1,nv2d)
  REAL(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP) :: v2dg(nlon,nlat,nv2d)
  INTEGER :: i,k,m,n
  CHARACTER(9) :: filename='file.0000'


  type(obs_info),intent(in) :: obs(nobsfiles)
  type(obs_da_value),intent(in) :: obsda

  REAL(r_size) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size) :: v2dgh(nlonh,nlath,nv2dd)
  integer :: proc,iv3d,iv2d,j
  real(r_size) :: ri,rj,rk

  real(r_size),allocatable :: tmpelm(:)
  real(r_size),allocatable :: ohx(:)
  integer,allocatable :: oqc(:)

  REAL(r_size) :: timer
  INTEGER :: ierr

!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  CALL CPU_TIME(timer)
!  if (myrank == 0) print *, '######', timer


  CALL ensmean_grd(MEMBER,nij1,v3d,v2d,v3dm,v2dm)


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  CALL CPU_TIME(timer)
!  if (myrank == 0) print *, '######', timer


  CALL gather_grd_mpi(lastmem_rank_e,v3dm,v2dm,v3dg,v2dg)


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  CALL CPU_TIME(timer)
!  if (myrank == 0) print *, '######', timer


  IF(myrank_e == lastmem_rank_e) THEN
    WRITE(filename(1:4),'(A4)') file
    WRITE(filename(6:9),'(A4)') 'mean'
!    WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
    call state_trans_inv(v3dg)
    call write_restart(filename,v3dg,v2dg)

CALL MPI_BARRIER(MPI_COMM_a,ierr)
CALL CPU_TIME(timer)
if (myrank == 0) print *, '######', timer

    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_u) = v3dg(:,:,:,iv3d_u)
    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_v) = v3dg(:,:,:,iv3d_v)
    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_w) = v3dg(:,:,:,iv3d_w)
    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_t) = v3dg(:,:,:,iv3d_t)
    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_p) = v3dg(:,:,:,iv3d_p)
    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_q) = v3dg(:,:,:,iv3d_q)
!    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qc) = v3dg(:,:,:,iv3d_qc)
!    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qr) = v3dg(:,:,:,iv3d_qr)
!    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qi) = v3dg(:,:,:,iv3d_qi)
!    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qs) = v3dg(:,:,:,iv3d_qs)
!    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qg) = v3dg(:,:,:,iv3d_qg)
!    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_rh) =
!    v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_hgt) =

!    v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_topo) =
!    v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_ps) =
!    v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_rain) =
!    v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_u10m) =
!    v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_v10m) =
!    v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_t2m) =
!    v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_q2m) =


    DO iv3d = 1, nv3dd
!!!!!$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
      do j  = JS, JE
        do i  = IS, IE
          v3dgh(   1:KS-1,i,j,iv3d) = v3dgh(KS,i,j,iv3d)
          v3dgh(KE+1:KA,  i,j,iv3d) = v3dgh(KE,i,j,iv3d)
        enddo
      enddo
    END DO

    DO iv3d = 1, nv3dd
      call COMM_vars8( v3dgh(:,:,:,iv3d), iv3d )
    END DO
    DO iv3d = 1, nv3dd
      call COMM_wait ( v3dgh(:,:,:,iv3d), iv3d )
    END DO

    DO iv2d = 1, nv2dd
      call COMM_vars8( v2dgh(:,:,iv2d), iv2d )
    END DO
    DO iv2d = 1, nv2dd
      call COMM_wait ( v2dgh(:,:,iv2d), iv2d )
    END DO



    allocate (tmpelm(obsda%nobs))
    allocate (ohx(obsda%nobs))
    allocate (oqc(obsda%nobs))

    oqc = -1

    do n = 1, obsda%nobs


      tmpelm(n) = obs(obsda%set(n))%elm(obsda%idx(n))

      if (obsda%qc(n) /= iqc_good) write(6, *) '############', obsda%qc(n)

      call rij_g2l_auto(proc,obsda%ri(n),obsda%rj(n),ri,rj)

      if (PRC_myrank /= proc) then
        write(6, *) '############ Error!'
        stop
      end if

!      .. create the v2dg....... (for ps)

      if (obs(obsda%set(n))%dif(obsda%idx(n)) >= -3600.0 .and. &   ! ###### 3600.0 as a variable
          obs(obsda%set(n))%dif(obsda%idx(n)) <= 3600.0 .and. &    ! ######
          (obs(obsda%set(n))%elm(obsda%idx(n)) == id_u_obs .or. &
           obs(obsda%set(n))%elm(obsda%idx(n)) == id_v_obs .or. &
           obs(obsda%set(n))%elm(obsda%idx(n)) == id_t_obs .or. &
           obs(obsda%set(n))%elm(obsda%idx(n)) == id_tv_obs .or. &
           obs(obsda%set(n))%elm(obsda%idx(n)) == id_q_obs)) then
!           obs(obsda%set(n))%elm(obsda%idx(n)) == id_rh_obs .or. &
!           obs(obsda%set(n))%elm(obsda%idx(n)) == id_ps_obs .or. &

        call phys2ijk(v3dg(:,:,:,iv3dd_p),obs(obsda%set(n))%elm(obsda%idx(n)), &
                      ri,rj,obs(obsda%set(n))%lev(obsda%idx(n)),rk,oqc(n))

        if (oqc(n) == iqc_good) then
          if (obs(obsda%set(n))%elm(obsda%idx(n)) == id_u_obs) then
            ri = ri - 0.5
!            if (ri < 1.0000001) ri = 1.0000001  ! ###### should modity itpl_3d to prevent '1.0' problem....
          else if (obs(obsda%set(n))%elm(obsda%idx(n)) == id_v_obs) then
            rj = rj - 0.5
!            if (rj < 1.0000001) rj = 1.0000001  ! ######
          end if
          call Trans_XtoY(obs(obsda%set(n))%elm(obsda%idx(n)),ri,rj,rk,v3dg,v2dg,ohx(n),oqc(n))
        end if

      end if

    end do


    CALL monit_dep(obsda%nobs,tmpelm,ohx,oqc,0)

    deallocate (tmpelm)
    deallocate (ohx)
    deallocate (oqc)


CALL MPI_BARRIER(MPI_COMM_a,ierr)
CALL CPU_TIME(timer)
if (myrank == 0) print *, '######', timer


  END IF


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  CALL CPU_TIME(timer)
!  if (myrank == 0) print *, '######', timer


  DO n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij1
        v3ds(i,k,n) = (v3d(i,k,1,n)-v3dm(i,k,n))**2
        DO m=2,MEMBER
          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n)-v3dm(i,k,n))**2
        END DO
        v3ds(i,k,n) = SQRT(v3ds(i,k,n) / REAL(MEMBER-1,r_size))
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  CALL CPU_TIME(timer)
!  if (myrank == 0) print *, '######', timer


  DO n=1,nv2d
!$OMP PARALLEL DO PRIVATE(i,m)
    DO i=1,nij1
      v2ds(i,n) = (v2d(i,1,n)-v2dm(i,n))**2
      DO m=2,MEMBER
        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n)-v2dm(i,n))**2
      END DO
      v2ds(i,n) = SQRT(v2ds(i,n) / REAL(MEMBER-1,r_size))
    END DO
!$OMP END PARALLEL DO
  END DO


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  CALL CPU_TIME(timer)
!  if (myrank == 0) print *, '######', timer


  CALL gather_grd_mpi(lastmem_rank_e,v3ds,v2ds,v3dg,v2dg)


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  CALL CPU_TIME(timer)
!  if (myrank == 0) print *, '######', timer


  IF(myrank_e == lastmem_rank_e) THEN
    WRITE(filename(1:4),'(A4)') file
    WRITE(filename(6:9),'(A4)') 'sprd'
!    WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
!    call state_trans_inv(v3dg)             !!
    call write_restart(filename,v3dg,v2dg)  !! not transformed to rho,rhou,rhov,rhow,rhot before writing.
  END IF


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  CALL CPU_TIME(timer)
!  if (myrank == 0) print *, '######', timer


  RETURN
END SUBROUTINE write_ensmspr_mpi



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
  implicit none
  integer, intent(in) :: proc
  integer, intent(in) :: ig
  integer, intent(in) :: jg
  integer, intent(out) :: il
  integer, intent(out) :: jl
  integer :: iproc, jproc

  call rank_1d_2d(proc, iproc, jproc)
  il = ig - iproc * nlon
  jl = jg - jproc * nlat

  return  
end subroutine ij_g2l


subroutine ij_l2g(proc, il, jl, ig, jg)
  implicit none
  integer, intent(in) :: proc
  integer, intent(in) :: il
  integer, intent(in) :: jl
  integer, intent(out) :: ig
  integer, intent(out) :: jg
  integer :: iproc, jproc

  call rank_1d_2d(proc, iproc, jproc)
  ig = il + iproc * nlon
  jg = jl + jproc * nlat

  return  
end subroutine ij_l2g


!-----------------------------------------------------------------------
! using halo!
! proc = -1: outside the global domain
!-----------------------------------------------------------------------
SUBROUTINE rij_g2l_auto(proc,ig,jg,il,jl)
  use scale_grid_index, only: &
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

  if (ig < real(1+IHALO,r_size) .or. ig > real(nlon*PRC_NUM_X+IHALO,r_size) .or. &
      jg < real(1+JHALO,r_size) .or. jg > real(nlat*PRC_NUM_Y+JHALO,r_size)) then
    il = -1.0d0
    jl = -1.0d0
    proc = -1
    return
  end if

  iproc = ceiling((ig-real(IHALO,r_size)-0.5d0) / real(nlon,r_size)) - 1
  jproc = ceiling((jg-real(JHALO,r_size)-0.5d0) / real(nlat,r_size)) - 1
  il = ig - iproc * nlon
  jl = jg - jproc * nlat
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
