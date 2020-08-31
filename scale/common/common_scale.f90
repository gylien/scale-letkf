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
  use common_nml

  use scale_precision, only: RP
  use scale_prof

  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

  ! Parameter 'nv3d' is set in common_nml.f90 ; 3D state variables (in SCALE restart files)
  ! Parameter 'nv2d' is set in common_nml.f90 ; 2D state variables (in SCALE restart files)
!  INTEGER,PARAMETER :: nv3dd=17  ! 3D diagnostic variables (in SCALE history files)
  INTEGER,PARAMETER :: nv3dd=25 ! for obssim ! 3D diagnostic variables (in SCALE history files)
#ifdef H08
  INTEGER,PARAMETER :: nv2dd=9  ! H08  ! 2D diagnostic variables (in SCALE history files)
#else
  INTEGER,PARAMETER :: nv2dd=7  ! 2D diagnostic variables (in SCALE history files)
#endif
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
  INTEGER,PARAMETER :: iv3d_q=6
  INTEGER,PARAMETER :: iv3d_qc=7
  INTEGER,PARAMETER :: iv3d_qr=8
  INTEGER,PARAMETER :: iv3d_qi=9
  INTEGER,PARAMETER :: iv3d_qs=10
  INTEGER,PARAMETER :: iv3d_qg=11

!  ! Number concentration (SN14 only) 
  INTEGER,PARAMETER :: iv3d_nc=18
  INTEGER,PARAMETER :: iv3d_nr=19
  INTEGER,PARAMETER :: iv3d_ni=20
  INTEGER,PARAMETER :: iv3d_ns=21
  INTEGER,PARAMETER :: iv3d_ng=22

  ! Lightning variables
  INTEGER,PARAMETER :: iv3d_cc=12
  INTEGER,PARAMETER :: iv3d_cr=13
  INTEGER,PARAMETER :: iv3d_ci=14
  INTEGER,PARAMETER :: iv3d_cs=15
  INTEGER,PARAMETER :: iv3d_cg=16
  INTEGER,PARAMETER :: iv3d_oep=17 

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

  INTEGER,PARAMETER :: iv3dd_pfl=14
  INTEGER,PARAMETER :: iv3dd_nfl=15
  INTEGER,PARAMETER :: iv3dd_fp=16 ! flash point
  INTEGER,PARAMETER :: iv3dd_ltp=17 ! LT path

! Not used for DA, only for obssim
  INTEGER,PARAMETER :: iv3dd_cc=18
  INTEGER,PARAMETER :: iv3dd_cr=19
  INTEGER,PARAMETER :: iv3dd_ci=20
  INTEGER,PARAMETER :: iv3dd_cs=21
  INTEGER,PARAMETER :: iv3dd_cg=22

  INTEGER,PARAMETER :: iv3dd_ex=23
  INTEGER,PARAMETER :: iv3dd_ey=24
  INTEGER,PARAMETER :: iv3dd_ez=25

  INTEGER,PARAMETER :: iv2dd_topo=1
  INTEGER,PARAMETER :: iv2dd_ps=2
  INTEGER,PARAMETER :: iv2dd_rain=3
  INTEGER,PARAMETER :: iv2dd_u10m=4
  INTEGER,PARAMETER :: iv2dd_v10m=5
  INTEGER,PARAMETER :: iv2dd_t2m=6
  INTEGER,PARAMETER :: iv2dd_q2m=7
#ifdef H08
  INTEGER,PARAMETER :: iv2dd_lsmask=8 ! H08
  INTEGER,PARAMETER :: iv2dd_skint=9 ! H08
#endif
!  INTEGER,PARAMETER :: iv2dd_tsfc=8

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

  INTEGER,PARAMETER :: vname_max = 10
  CHARACTER(vname_max),SAVE :: v3d_name(nv3d)
  CHARACTER(vname_max),SAVE :: v3dd_name(nv3dd)
  CHARACTER(vname_max),SAVE :: v2d_name(nv2d)
  CHARACTER(vname_max),SAVE :: v2dd_name(nv2dd)


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

  real(r_size), allocatable, save :: count_CDF(:)
  real(r_size), allocatable, save :: flash_CDF(:)


CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
! mem:
!  >0: mem
!  0:  MEMBER
!  -1: MEMBER+1
!-----------------------------------------------------------------------
SUBROUTINE set_common_scale
  use scale_rm_process, only: &
    PRC_NUM_X, &
    PRC_NUM_Y
  use scale_grid_index, only: &
    IMAX, &
    JMAX, &
    KMAX, &
    IHALO, &
    JHALO, &
    KHALO

  IMPLICIT NONE
!  integer, intent(in) :: mem
!  REAL(r_sngl) :: slat(nlat), wlat(nlat)
!  REAL(r_size) :: totalwg, wgtmp, latm1, latm2
!  INTEGER :: i,j

  WRITE(6,'(A)') 'Hello from set_common_scale'

!  ! setup standard I/O
!  call IO_setup( MODELNAME, .false.)

  !
  ! Set up node and process distribution
  !
!  if (mem == 0) then
!    call set_mem_node_proc(MEMBER,NNODES,PPN,MEM_NODES,MEM_NP)
!  else if (mem == -1) then
!    call set_mem_node_proc(MEMBER+1,NNODES,PPN,MEM_NODES,MEM_NP)
!  else
!    call set_mem_node_proc(mem,NNODES,PPN,MEM_NODES,MEM_NP)
!  end if

!! print process distribution!!!

!  if (scale_IO_mygroup <= 0) then

!    write (6, '(A,I6.6,A)') 'MYRANK=',myrank,': This process is not used!'

!  else

!    call set_scalelib

    if (MEM_NP /= PRC_NUM_X * PRC_NUM_Y) then
      write(6,'(A,I10)') 'MEM_NP    = ', MEM_NP
      write(6,'(A,I10)') 'PRC_NUM_X = ', PRC_NUM_X
      write(6,'(A,I10)') 'PRC_NUM_Y = ', PRC_NUM_Y
      write(6,'(A)') 'MEM_NP should be equal to PRC_NUM_X * PRC_NUM_Y.'
      stop
    end if

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

    !
    ! Variable names (same as in the NetCDF file)
    !
    ! state variables (in 'restart' files, for LETKF)
    v3d_name(iv3d_rho)  = 'DENS'
    v3d_name(iv3d_rhou) = 'MOMX'
    v3d_name(iv3d_rhov) = 'MOMY'
    v3d_name(iv3d_rhow) = 'MOMZ'
    v3d_name(iv3d_rhot) = 'RHOT'
    v3d_name(iv3d_q)    = 'QV'
    v3d_name(iv3d_qc)   = 'QC'
    v3d_name(iv3d_qr)   = 'QR'
    v3d_name(iv3d_qi)   = 'QI'
    v3d_name(iv3d_qs)   = 'QS'
    v3d_name(iv3d_qg)   = 'QG'

!    ! SN14 only
    v3d_name(iv3d_nc)   = 'NC'
    v3d_name(iv3d_nr)   = 'NR'
    v3d_name(iv3d_ni)   = 'NI'
    v3d_name(iv3d_ns)   = 'NS'
    v3d_name(iv3d_ng)   = 'NG'

    v3d_name(iv3d_cc)   = 'CDNS_QC'
    v3d_name(iv3d_cr)   = 'CDNS_QR'
    v3d_name(iv3d_ci)   = 'CDNS_QI'
    v3d_name(iv3d_cs)   = 'CDNS_QS'
    v3d_name(iv3d_cg)   = 'CDNS_QG'
    v3d_name(iv3d_oep)   = 'Epot_old'

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

    v3dd_name(iv3dd_pfl)   = 'PosFLASH'
    v3dd_name(iv3dd_nfl)   = 'NegFLASH'
    v3dd_name(iv3dd_fp)   = 'FlashPoint'
    v3dd_name(iv3dd_ltp)   = 'LTpath'

!    ! For obssim
    v3dd_name(iv3dd_cc)   = 'QCRG_C'
    v3dd_name(iv3dd_cr)   = 'QCRG_R'
    v3dd_name(iv3dd_ci)   = 'QCRG_I'
    v3dd_name(iv3dd_cs)   = 'QCRG_S'
    v3dd_name(iv3dd_cg)   = 'QCRG_G'

    v3dd_name(iv3dd_ex)   = 'Ex'
    v3dd_name(iv3dd_ey)   = 'Ey'
    v3dd_name(iv3dd_ez)   = 'Ez'
    !
    v2dd_name(iv2dd_topo) = 'topo'
    v2dd_name(iv2dd_ps) = 'SFC_PRES'
    v2dd_name(iv2dd_rain) = 'PREC'
    v2dd_name(iv2dd_u10m) = 'U10'
    v2dd_name(iv2dd_v10m) = 'V10'
    v2dd_name(iv2dd_t2m) = 'T2'
    v2dd_name(iv2dd_q2m) = 'Q2'
#ifdef H08
    v2dd_name(iv2dd_lsmask) = 'lsmask' ! H08
    v2dd_name(iv2dd_skint) = 'SFC_TEMP' ! H08
#endif
!    v2dd_name(iv2dd_tsfc) = 'SFC_TEMP'
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


!  end if ! [ scale_IO_mygroup <= 0 ]


  RETURN
END SUBROUTINE set_common_scale

!SUBROUTINE unset_common_scale
!  IMPLICIT NONE

!  if (scale_IO_mygroup > 0) then
!    call unset_scalelib
!  end if

!  RETURN
!END SUBROUTINE unset_common_scale



subroutine set_common_conf(nprocs)
  use scale_stdio

  implicit none
  integer, intent(in) :: nprocs
  character(len=H_MID), parameter :: MODELNAME = "SCALE-LETKF"

  ! setup standard I/O
  call IO_setup( MODELNAME, .false.)

  call read_nml_ensemble
  call read_nml_letkf_prc

  if (nprocs /= NNODES * PPN) then
    write(6,'(A,I10)') 'Number of MPI processes = ', nprocs
    write(6,'(A,I10)') 'NNODES = ', NNODES
    write(6,'(A,I10)') 'PPN    = ', PPN
    write(6,'(A)') 'Number of MPI processes should be equal to NNODES * PPN.'
    stop
  end if

  return
end subroutine set_common_conf


!!-----------------------------------------------------------------------
!! set_mem2proc
!!-----------------------------------------------------------------------
!SUBROUTINE set_mem_node_proc(mem,nnodes,ppn,mem_nodes,mem_np)
!  IMPLICIT NONE
!  INTEGER,INTENT(IN) :: mem,nnodes,ppn,mem_nodes,mem_np
!  INTEGER :: tppn,tppnt,tmod
!  INTEGER :: n,ns,nn,m,q,qs,i,j,it,ip

!  ALLOCATE(procs(nprocs))
!  ns = 0
!  DO n = 1, nnodes
!    procs(ns+1:ns+ppn) = n
!    ns = ns + ppn
!  END DO

!  IF(mem_nodes > 1) THEN
!    n_mem = nnodes / mem_nodes
!    n_mempn = 1
!  ELSE
!    n_mem = nnodes
!    n_mempn = ppn / mem_np
!  END IF
!  nitmax = (mem - 1) / (n_mem * n_mempn) + 1
!  tppn = mem_np / mem_nodes
!  tmod = MOD(mem_np, mem_nodes)

!  ALLOCATE(mem2node(mem_np,mem))
!  ALLOCATE(mem2proc(mem_np,mem))
!  ALLOCATE(proc2mem(2,nitmax,nprocs))
!  proc2mem = -1
!  m = 1
!mem_loop: DO it = 1, nitmax
!    DO i = 0, n_mempn-1
!      n = 0
!      DO j = 0, n_mem-1
!        IF(m > mem .and. it > 1) EXIT mem_loop
!        qs = 0
!        DO nn = 0, mem_nodes-1
!          IF(nn < tmod) THEN
!            tppnt = tppn + 1
!          ELSE
!            tppnt = tppn
!          END IF
!          DO q = 0, tppnt-1
!            ip = (n+nn)*ppn + i*mem_np + q
!            if (m <= mem) then
!              mem2node(qs+1,m) = n+nn
!              mem2proc(qs+1,m) = ip
!            end if
!            proc2mem(1,it,ip+1) = m
!            proc2mem(2,it,ip+1) = qs
!            qs = qs + 1
!          END DO
!        END DO
!        m = m + 1
!        n = n + mem_nodes
!      END DO
!    END DO
!  END DO mem_loop

!  ens_mygroup = proc2mem(1,1,myrank+1)
!  ens_myrank = proc2mem(2,1,myrank+1)
!  if (ens_mygroup >= 1) then
!    myrank_use = .true.
!  end if

!  lastmem_rank_e = mod(mem-1, n_mem*n_mempn)
!!  if (lastmem_rank_e /= proc2mem(1,1,mem2proc(1,mem)+1)-1) then
!!    print *, 'XXXXXX wrong!!'
!!    stop
!!  end if

!  RETURN
!END SUBROUTINE



!!-----------------------------------------------------------------------
!! File I/O
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!!-----------------------------------------------------------------------
!SUBROUTINE read_restart(filename,v3dg,v2dg)
!  use gtool_file, only: FileRead, FileCloseAll
!  use scale_process, only: PRC_myrank
!  use common_mpi, only: myrank
!  IMPLICIT NONE

!  CHARACTER(*),INTENT(IN) :: filename
!  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
!  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
!  INTEGER :: iv3d,iv2d

!  WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',PRC_myrank,'.nc'

!  do iv3d = 1, nv3d
!    write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
!    call FileRead(v3dg(:,:,:,iv3d), filename, trim(v3d_name(iv3d)), 1, PRC_myrank)
!  end do

!  do iv2d = 1, nv2d
!    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
!    call FileRead(v2dg(:,:,iv2d), filename, trim(v2d_name(iv2d)), 1, PRC_myrank)
!  end do

!  call FileCloseAll

!  RETURN
!END SUBROUTINE read_restart
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE read_restart(filename,v3dg,v2dg)
  use netcdf, only: NF90_NOWRITE
  use scale_process, only: &
    PRC_myrank
  use scale_rm_process, only: &
    PRC_HAS_W,  &
    PRC_HAS_E,  &
    PRC_HAS_S,  &
    PRC_HAS_N
  use scale_grid_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX, &
    IMAXB, JMAXB
  use common_mpi, only: myrank
  use common_ncio
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: iv3d,iv2d,ncid
  integer :: is, ie, js, je
  real(RP) :: v3dgtmp(KMAX,IMAXB,JMAXB)
  real(RP) :: v2dgtmp(IMAXB,JMAXB)

  is = 1
  ie = IMAX
  js = 1
  je = JMAX
  if ( .not. PRC_HAS_W ) then
    is = is + IHALO
    ie = ie + IHALO
  end if
  if ( .not. PRC_HAS_S ) then
    js = js + JHALO
    je = je + JHALO
  end if

!  write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',PRC_myrank,'.nc'

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)

  do iv3d = 1, nv3d
    write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
    call ncio_read(ncid, trim(v3d_name(iv3d)), KMAX, IMAXB, JMAXB, 1, v3dgtmp)
    v3dg(:,:,:,iv3d) = v3dgtmp(:,is:ie,js:je)
  end do

  do iv2d = 1, nv2d
    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
    call ncio_read(ncid, trim(v2d_name(iv2d)), IMAXB, JMAXB, 1, v2dgtmp)
    v2dg(:,:,iv2d) = v2dgtmp(is:ie,js:je)
  end do

  call ncio_close(ncid)

  RETURN
END SUBROUTINE read_restart
!!-----------------------------------------------------------------------
!!
!!-----------------------------------------------------------------------
!SUBROUTINE write_restart(filename,v3dg,v2dg)
!  use netcdf, only: NF90_WRITE

!!  use gtool_file, only: FileOpen, FileClose, FileWrite
!!  use gtool_file_h


!  use scale_process, only: PRC_myrank
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
!    write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
!    call ncio_write(ncid, trim(v3d_name(iv3d)), nlev, nlon, nlat, 1, v3dg(:,:,:,iv3d))
!    call ncio_read (ncid, trim(v3d_name(iv3d)), nlev, nlon, nlat, 1, v3dgtmp)
!    call ncio_write(ncid, trim(v3d_name(iv3d)), nlev, nlon, nlat, 1, v3dgtmp)
!  END DO

!  DO iv2d = 1, nv2d
!    write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
!    call ncio_write(ncid, trim(v2d_name(iv2d)), nlon, nlat, 1, v2dg(:,:,iv2d))
!    call ncio_read (ncid, trim(v2d_name(iv2d)), nlon, nlat, 1, v2dgtmp)
!    call ncio_write(ncid, trim(v2d_name(iv2d)), nlon, nlat, 1, v2dgtmp)
!  END DO

!!  DO iv3d = 1, nv3d
!!    write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
!!    call ncio_check(nf90_inq_varid(ncid, trim(v3d_name(iv3d)), vid))

!!!print *, '######### vid:', myrank, vid

!!!    call FileWrite(vid, v3dg(:,:,:,iv3d), -1.0d0, -1.0d0)
!!    call file_write_data( vid, v3dg(:,:,:,iv3d), -1.0d0, -1.0d0, RP, & ! (in)
!!         error                                     ) ! (out)
!!  END DO

!!  DO iv2d = 1, nv2d
!!    write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
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
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE write_restart(filename,v3dg,v2dg)
  use netcdf, only: NF90_WRITE
  use scale_process, only: &
    PRC_myrank
  use scale_rm_process, only: &
    PRC_HAS_W,  &
    PRC_HAS_E,  &
    PRC_HAS_S,  &
    PRC_HAS_N
  use scale_grid_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX, &
    IMAXB, JMAXB
  use common_mpi, only: myrank
  use common_ncio
  implicit none

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: iv3d,iv2d,ncid
  integer :: is, ie, js, je
  real(RP) :: v3dgtmp(KMAX,IMAXB,JMAXB)
  real(RP) :: v2dgtmp(IMAXB,JMAXB)

  is = 1
  ie = IMAX
  js = 1
  je = JMAX
  if ( .not. PRC_HAS_W ) then
    is = is + IHALO
    ie = ie + IHALO
  end if
  if ( .not. PRC_HAS_S ) then
    js = js + JHALO
    je = je + JHALO
  end if

!  write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',PRC_myrank,'.nc'

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is writing a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_WRITE, ncid)

  do iv3d = 1, nv3d
!    write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
    call ncio_read (ncid, trim(v3d_name(iv3d)), KMAX, IMAXB, JMAXB, 1, v3dgtmp)
    v3dgtmp(:,is:ie,js:je) = v3dg(:,:,:,iv3d)
    call ncio_write(ncid, trim(v3d_name(iv3d)), KMAX, IMAXB, JMAXB, 1, v3dgtmp)
!    call ncio_read (ncid, trim(v3d_name(iv3d)), KMAX, IMAXB, JMAXB, 1, v3dgtmp)  !!! read and write again to work around the endian problem on the K computer
!    call ncio_write(ncid, trim(v3d_name(iv3d)), KMAX, IMAXB, JMAXB, 1, v3dgtmp)  !
  end do

  do iv2d = 1, nv2d
!    write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
    call ncio_read (ncid, trim(v2d_name(iv2d)), IMAXB, JMAXB, 1, v2dgtmp)
    v2dgtmp(is:ie,js:je) = v2dg(:,:,iv2d)
    call ncio_write(ncid, trim(v2d_name(iv2d)), IMAXB, JMAXB, 1, v2dgtmp)
!    call ncio_read (ncid, trim(v2d_name(iv2d)), IMAXB, JMAXB, 1, v2dgtmp)  !!! read and write again to work around the endian problem on the K computer
!    call ncio_write(ncid, trim(v2d_name(iv2d)), IMAXB, JMAXB, 1, v2dgtmp)  !
  end do

  call ncio_close(ncid)

  RETURN
END SUBROUTINE write_restart

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE read_topo(filename,topo)
  use netcdf, only: NF90_NOWRITE
  use scale_process, only: &
    PRC_myrank
  use scale_rm_process, only: &
    PRC_HAS_W,  &
    PRC_HAS_E,  &
    PRC_HAS_S,  &
    PRC_HAS_N
  use scale_grid_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, &
    IMAXB, JMAXB
  use common_mpi, only: myrank
  use common_ncio
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: topo(nlon,nlat)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: ncid
  integer :: is, ie, js, je
  real(RP) :: v2dgtmp(IMAXB,JMAXB)

  is = 1
  ie = IMAX
  js = 1
  je = JMAX
  if ( .not. PRC_HAS_W ) then
    is = is + IHALO
    ie = ie + IHALO
  end if
  if ( .not. PRC_HAS_S ) then
    js = js + JHALO
    je = je + JHALO
  end if

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)

  write(6,'(1x,A,A15)') '*** Read 2D var: ', 'TOPO'
  call ncio_read(ncid, 'TOPO', IMAXB, JMAXB, 1, v2dgtmp)
  topo = v2dgtmp(is:ie,js:je)

  call ncio_close(ncid)

  RETURN
END SUBROUTINE read_topo




!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine read_history(filename,step,v3dg,v2dg)
  use scale_process, only: &
      PRC_myrank
  use scale_grid_index, only: &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA
  use gtool_history, only: &
      HistoryGet
  use scale_comm, only: &
      COMM_vars8, &
      COMM_wait
  use common_mpi, only: myrank
  implicit none

  character(*),intent(in) :: filename
  integer,intent(in) :: step
  real(r_size),intent(out) :: v3dg(nlevh,nlonh,nlath,nv3dd)
  real(r_size),intent(out) :: v2dg(nlonh,nlath,nv2dd)
  integer :: i,j,k,iv3d,iv2d
  character(len=12) :: filesuffix = '.pe000000.nc'
  real(RP) :: var3d(nlon,nlat,nlev)
  real(RP) :: var2d(nlon,nlat)

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix

  ! 3D variables
  !-------------
  do iv3d = 1, nv3dd
    write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3dd_name(iv3d))
    call HistoryGet( var3D,                 & ! [OUT]
                     filename,              & ! [IN]
                     trim(v3dd_name(iv3d)), & ! [IN]
                     step                   ) ! [IN]
    forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dg(k+KHALO,i+IHALO,j+JHALO,iv3d) = var3D(i,j,k) ! use FORALL to change order of dimensions
  end do

  do iv3d = 1, nv3dd
    call COMM_vars8( v3dg(:,:,:,iv3d), iv3d )
  end do
  do iv3d = 1, nv3dd
    call COMM_wait ( v3dg(:,:,:,iv3d), iv3d )
  end do

  do iv3d = 1, nv3dd
!!!!!!!$OMP PARALLEL DO PRIVATE(i,j) OMP_SCHEDULE_ COLLAPSE(2)
    do j = JS, JE
      do i = IS, IE
        v3dg(   1:KS-1,i,j,iv3d) = v3dg(KS,i,j,iv3d)
        v3dg(KE+1:KA,  i,j,iv3d) = v3dg(KE,i,j,iv3d)
      end do
    end do
  end do

  ! 2D variables
  !-------------
  do iv2d = 1, nv2dd
    write(6,'(1x,A,A15)') '*** NOT Read 2D var: ', trim(v2dd_name(iv2d))
!    call HistoryGet( var2D,                 & ! [OUT]
!                     filename,              & ! [IN]
!                     trim(v2dd_name(iv2d)), & ! [IN]
!                     step                   ) ! [IN]
!    v2dg(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2d) = var2D(:,:)
  end do

!  do iv2d = 1, nv2dd
!    call COMM_vars8( v2dg(:,:,iv2d), iv2d )
!  end do
!  do iv2d = 1, nv2dd
!    call COMM_wait ( v2dg(:,:,iv2d), iv2d )
!  end do

  return
end subroutine read_history
!!-----------------------------------------------------------------------
!! (version handling data type conversion)
!!-----------------------------------------------------------------------
!subroutine read_history(filename,step,v3dg,v2dg)
!  use scale_process, only: &
!      PRC_myrank
!  use scale_grid_index, only: &
!      IHALO, JHALO, KHALO, &
!      IS, IE, JS, JE, KS, KE, KA
!  use gtool_history, only: &
!      HistoryGet
!  use scale_comm, only: &
!      COMM_vars8, &
!      COMM_wait
!  use common_mpi, only: myrank
!  implicit none

!  character(*),intent(in) :: filename
!  integer,intent(in) :: step
!  real(r_size),intent(out) :: v3dg(nlevh,nlonh,nlath,nv3dd)
!  real(r_size),intent(out) :: v2dg(nlonh,nlath,nv2dd)
!  integer :: i,j,k,iv3d,iv2d
!  character(len=12) :: filesuffix = '.pe000000.nc'
!  real(RP) :: var3d(nlon,nlat,nlev)
!  real(RP) :: var2d(nlon,nlat)
!  real(RP) :: v3dgtmp(nlevh,nlonh,nlath,nv3dd)
!  real(RP) :: v2dgtmp(nlonh,nlath,nv2dd)

!  write (filesuffix(4:9),'(I6.6)') PRC_myrank
!  write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix

!  ! 3D variables
!  !-------------
!  do iv3d = 1, nv3dd
!    write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3dd_name(iv3d))
!    call HistoryGet( var3D,                 & ! [OUT]
!                     filename,              & ! [IN]
!                     trim(v3dd_name(iv3d)), & ! [IN]
!                     step                   ) ! [IN]
!    forall (i=1:nlon, j=1:nlat, k=1:nlev) v3dgtmp(k+KHALO,i+IHALO,j+JHALO,iv3d) = var3D(i,j,k) ! use FORALL to change order of dimensions
!  end do

!  do iv3d = 1, nv3dd
!    call COMM_vars8( v3dgtmp(:,:,:,iv3d), iv3d )
!  end do
!  do iv3d = 1, nv3dd
!    call COMM_wait ( v3dgtmp(:,:,:,iv3d), iv3d )
!  end do
!  v3dg = real(v3dgtmp, r_size)

!  do iv3d = 1, nv3dd
!!!!!!$OMP PARALLEL DO PRIVATE(i,j) OMP_SCHEDULE_ COLLAPSE(2)
!    do j = JS, JE
!      do i = IS, IE
!        v3dg(   1:KS-1,i,j,iv3d) = v3dg(KS,i,j,iv3d)
!        v3dg(KE+1:KA,  i,j,iv3d) = v3dg(KE,i,j,iv3d)
!      end do
!    end do
!  end do

!  ! 2D variables
!  !-------------
!  do iv2d = 1, nv2dd
!    write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2dd_name(iv2d))
!    call HistoryGet( var2D,                 & ! [OUT]
!                     filename,              & ! [IN]
!                     trim(v2dd_name(iv2d)), & ! [IN]
!                     step                   ) ! [IN]
!    v2dgtmp(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2d) = var2D(:,:)
!  end do

!  do iv2d = 1, nv2dd
!    call COMM_vars8( v2dgtmp(:,:,iv2d), iv2d )
!  end do
!  do iv2d = 1, nv2dd
!    call COMM_wait ( v2dgtmp(:,:,iv2d), iv2d )
!  end do
!  v2dg = real(v2dgtmp, r_size)

!  return
!end subroutine read_history



!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine scale_calc_z(nij,topo,z)
  use scale_grid, only: &
     GRID_CZ, &
     GRID_FZ
  use scale_grid_index, only: &
     KHALO, KS, KE
  implicit none

  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: topo(nij)
  REAL(RP),INTENT(OUT) :: z(nij,nlev)
  real(r_size) :: ztop
  integer  :: k, i

  ztop = GRID_FZ(KE) - GRID_FZ(KS-1)
!$OMP PARALLEL DO PRIVATE(i,k)
  do k = 1, nlev
    do i = 1, nij
      z(i,k) = (ztop - topo(i)) / ztop * GRID_CZ(k+KHALO) + topo(i)
    end do
  end do
!$OMP END PARALLEL DO

  return
end subroutine scale_calc_z


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
SUBROUTINE ensmean_grd(mem,nij,v3d,v2d,v3dm,v2dm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: mem
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,mem,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,mem,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: i,k,m,n

  DO n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        DO m=2,mem
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        END DO
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(mem,r_size)
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO

  DO n=1,nv2d
!$OMP PARALLEL DO PRIVATE(i,m)
    DO i=1,nij
      v2dm(i,n) = v2d(i,1,n)
      DO m=2,mem
        v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
      END DO
      v2dm(i,n) = v2dm(i,n) / REAL(mem,r_size)
    END DO
!$OMP END PARALLEL DO
  END DO

  RETURN
END SUBROUTINE ensmean_grd




SUBROUTINE state_trans(v3dg)
!  use gtool_file, only: FileRead
!  use scale_process, only: PRC_myrank
!  use common_mpi, only: myrank
!  use scale_history, only: &
!     HIST_get

  use scale_atmos_thermodyn, only: AQ_CV
    use scale_const, only: &
       Rdry   => CONST_Rdry, &
       Rvap   => CONST_Rvap, &
       CVdry  => CONST_CVdry, &
       PRE00 => CONST_PRE00

  IMPLICIT NONE

  REAL(RP),INTENT(INOUT) :: v3dg(nlev,nlon,nlat,nv3d)
!  REAL(RP) :: pres(nlev,nlon,nlat)
!  REAL(RP) :: temp(nlev,nlon,nlat)
  REAL(RP) :: rho,pres,temp
  real(RP) :: qdry,CVtot,Rtot,CPovCV
  integer :: i,j,k,iv3d


!  REAL(RP) :: pres2(nlev,nlon,nlat)

!!!$OMP PARALLEL DO PRIVATE(i,j,k,iv3d,qdry,CVtot,Rtot,CPovCV,rho,pres,temp) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev
       qdry  = 1.0d0
       CVtot = 0.0d0
!       do iv3d = iv3d_q, nv3d ! loop over all moisture variables
       do iv3d = iv3d_q, iv3d_qg ! loop over all moisture variables
         qdry  = qdry - v3dg(k,i,j,iv3d)
         CVtot = CVtot + v3dg(k,i,j,iv3d) * AQ_CV(iv3d-iv3d_q+1)
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

       v3dg(k,i,j,iv3d_cc) = v3dg(k,i,j,iv3d_cc) !* rho  ! fC/kg => fC/m^3
       v3dg(k,i,j,iv3d_cr) = v3dg(k,i,j,iv3d_cr) !* rho  ! fC/kg => fC/m^3
       v3dg(k,i,j,iv3d_ci) = v3dg(k,i,j,iv3d_ci) !* rho  ! fC/kg => fC/m^3
       v3dg(k,i,j,iv3d_cs) = v3dg(k,i,j,iv3d_cs) !* rho  ! fC/kg => fC/m^3
       v3dg(k,i,j,iv3d_cg) = v3dg(k,i,j,iv3d_cg) !* rho  ! fC/kg => fC/m^3

       if ( N_LOG_TRANS ) then
         v3dg(k,i,j,iv3d_nc) = log( v3dg(k,i,j,iv3d_nc) + 1.0d0 )
         v3dg(k,i,j,iv3d_nr) = log( v3dg(k,i,j,iv3d_nr) + 1.0d0 )
         v3dg(k,i,j,iv3d_ni) = log( v3dg(k,i,j,iv3d_ni) + 1.0d0 )
         v3dg(k,i,j,iv3d_ns) = log( v3dg(k,i,j,iv3d_ns) + 1.0d0 )
         v3dg(k,i,j,iv3d_ng) = log( v3dg(k,i,j,iv3d_ng) + 1.0d0 )
       endif

      enddo
    enddo
  enddo
!!!$OMP END PARALLEL DO

!  v3dg(:,:,:,5) = v3dg(:,:,:,iv3d_rho) ! temporarily store density here
!  v3dg(:,:,:,iv3d_u) = v3dg(:,:,:,iv3d_rhou) / v3dg(:,:,:,5)
!  v3dg(:,:,:,iv3d_v) = v3dg(:,:,:,iv3d_rhov) / v3dg(:,:,:,5)
!  v3dg(:,:,:,iv3d_w) = v3dg(:,:,:,iv3d_rhow) / v3dg(:,:,:,5)
!  v3dg(:,:,:,iv3d_t) = temp
!  v3dg(:,:,:,iv3d_p) = pres

!!write(6,*) pres(:,10,10)
!write(6,*) v3dg(:,10,10,iv3d_t)
!!write(6,*) v3dg(6,:,:,iv3d_rhou)

!!write(6,*)

!if(myrank == 0) then

!  call HIST_get(pres2, 'hist.0001', 'T', 15)
!  write(6,*) pres2(10,10,:)

!!  call HIST_get(pres2, 'hist.0001', 'MOMX', 15)
!!  write(6,*) pres2(:,:,6)

!!  do i = 1, 15
!!  call HIST_get(pres2, 'hist.0001', 'MOMX', i)
!!  write(6,*) pres2(10,10,6)
!!  end do

!end if

!write(6,*)


  RETURN
END SUBROUTINE state_trans


SUBROUTINE state_trans_inv(v3dg)
  use scale_atmos_thermodyn, only: AQ_CP, AQ_CV
    use scale_const, only: &
       Rdry   => CONST_Rdry, &
       Rvap   => CONST_Rvap, &
       CVdry  => CONST_CVdry, &
       PRE00 => CONST_PRE00
  IMPLICIT NONE

  REAL(RP),INTENT(INOUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP) :: rho,rhot
  real(RP) :: qdry,CVtot,Rtot,CVovCP
  integer :: i,j,k,iv3d

  do iv3d = 1, nv3d
    if (POSITIVE_DEFINITE_Q .and. iv3d == iv3d_q) then
      v3dg(:,:,:,iv3d) = max(v3dg(:,:,:,iv3d), 0.0d0)
    else if (POSITIVE_DEFINITE_QHYD .and. (.not. POSITIVE_DEFINITE_QHYD_QCRG) .and. (iv3d == iv3d_qc .or. iv3d == iv3d_qr .or. iv3d == iv3d_qi .or. iv3d == iv3d_qs .or. iv3d == iv3d_qg)) then
      v3dg(:,:,:,iv3d) = max(v3dg(:,:,:,iv3d), 0.0d0)
    end if
  end do

  if (POSITIVE_DEFINITE_QHYD .and. POSITIVE_DEFINITE_QHYD_QCRG .and. .not. N_LOG_TRANS) then
!$OMP PARALLEL DO PRIVATE(i,j,k) COLLAPSE(2)
    do j = 1, nlat
      do i = 1, nlon
        do k = 1, nlev
          if (v3dg(k,i,j,iv3d_qc) <= 0.0d0 .or. v3dg(k,i,j,iv3d_nc) <= 0.0d0 ) then
            v3dg(k,i,j,iv3d_cc) = 0.0d0
            v3dg(k,i,j,iv3d_qc) = 0.0d0
            v3dg(k,i,j,iv3d_nc) = 0.0d0
          endif
          if (v3dg(k,i,j,iv3d_qr) <= 0.0d0 .or. v3dg(k,i,j,iv3d_nr) <= 0.0d0 ) then
            v3dg(k,i,j,iv3d_cr) = 0.0d0
            v3dg(k,i,j,iv3d_qr) = 0.0d0
            v3dg(k,i,j,iv3d_nr) = 0.0d0
          endif
          if (v3dg(k,i,j,iv3d_qi) <= 0.0d0 .or. v3dg(k,i,j,iv3d_ni) <= 0.0d0 ) then
            v3dg(k,i,j,iv3d_ci) = 0.0d0
            v3dg(k,i,j,iv3d_qi) = 0.0d0
            v3dg(k,i,j,iv3d_ni) = 0.0d0
          endif
          if (v3dg(k,i,j,iv3d_qs) <= 0.0d0 .or. v3dg(k,i,j,iv3d_ns) <= 0.0d0 ) then
            v3dg(k,i,j,iv3d_cs) = 0.0d0
            v3dg(k,i,j,iv3d_qs) = 0.0d0
            v3dg(k,i,j,iv3d_ns) = 0.0d0
          endif
          if (v3dg(k,i,j,iv3d_qg) <= 0.0d0 .or. v3dg(k,i,j,iv3d_ng) <= 0.0d0 ) then
            v3dg(k,i,j,iv3d_cg) = 0.0d0
            v3dg(k,i,j,iv3d_qg) = 0.0d0
            v3dg(k,i,j,iv3d_ng) = 0.0d0
          endif
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO
  endif

!$OMP PARALLEL DO PRIVATE(i,j,k,iv3d,qdry,CVtot,Rtot,CVovCP,rho,rhot) COLLAPSE(2)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev
       qdry  = 1.0d0
       CVtot = 0.0d0
!       do iv3d = iv3d_q, nv3d ! loop over all moisture variables
       do iv3d = iv3d_q, iv3d_qg ! loop over all moisture variables
         qdry  = qdry - v3dg(k,i,j,iv3d)
         CVtot = CVtot + v3dg(k,i,j,iv3d) * AQ_CV(iv3d-iv3d_q+1)
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

       v3dg(k,i,j,iv3d_cc) = v3dg(k,i,j,iv3d_cc) !/ rho  ! fC/m^3 => fC/kg
       v3dg(k,i,j,iv3d_cr) = v3dg(k,i,j,iv3d_cr) !/ rho  ! fC/m^3 => fC/kg
       v3dg(k,i,j,iv3d_ci) = v3dg(k,i,j,iv3d_ci) !/ rho  ! fC/m^3 => fC/kg
       v3dg(k,i,j,iv3d_cs) = v3dg(k,i,j,iv3d_cs) !/ rho  ! fC/m^3 => fC/kg
       v3dg(k,i,j,iv3d_cg) = v3dg(k,i,j,iv3d_cg) !/ rho  ! fC/m^3 => fC/kg

        if ( N_LOG_TRANS ) then
          v3dg(k,i,j,iv3d_nc) = exp( v3dg(k,i,j,iv3d_nc) ) - 1.0d0
          v3dg(k,i,j,iv3d_nr) = exp( v3dg(k,i,j,iv3d_nr) ) - 1.0d0
          v3dg(k,i,j,iv3d_ni) = exp( v3dg(k,i,j,iv3d_ni) ) - 1.0d0
          v3dg(k,i,j,iv3d_ns) = exp( v3dg(k,i,j,iv3d_ns) ) - 1.0d0
          v3dg(k,i,j,iv3d_ng) = exp( v3dg(k,i,j,iv3d_ng) ) - 1.0d0
        endif
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

  RETURN
END SUBROUTINE state_trans_inv


!-------------------------------------------------------------------------------
! Transform the state variables (from SCALE restart files) to 
! the variables in SCALE history files (with HALO), 
! so that they can be used for observation operator calculation
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

  real(r_size) :: ztop
  integer :: i, j, k, iv3d, iv2d

  ! Variables that can be directly copied
  !-----------------------------------------------------------------------------

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

  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_cc) = v3dg(:,:,:,iv3d_cc)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_cr) = v3dg(:,:,:,iv3d_cr)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_ci) = v3dg(:,:,:,iv3d_ci)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_cs) = v3dg(:,:,:,iv3d_cs)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_cg) = v3dg(:,:,:,iv3d_cg)

  ! RH
  !-----------------------------------------------------------------------------

!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_rh) = [[RH calculator]]

  ! Calculate height based the the topography and vertical coordinate
  !-----------------------------------------------------------------------------

  ztop = GRID_FZ(KE) - GRID_FZ(KS-1)
!!!$OMP PARALLEL DO PRIVATE(j,i,k)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev+KHALO ! T. Honda fullfilling the top halo levels
        v3dgh(k+KHALO, i+IHALO, j+JHALO, iv3dd_hgt) = (ztop - topo(i,j)) / ztop * GRID_CZ(k+KHALO) + topo(i,j)
      end do
    enddo
  enddo
!!!$OMP END PARALLEL DO

  ! Surface variables: use the 1st level as the surface (although it is not)
  !-----------------------------------------------------------------------------

  v2dgh(:,:,iv2dd_topo) = v3dgh(1+KHALO,:,:,iv3dd_hgt)                ! Use the first model level as topography (is this good?)
!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_topo) = topo(:,:) ! Use the real topography

  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_ps) = v3dg(1,:,:,iv3d_p)
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_u10m) = v3dg(1,:,:,iv3d_u)
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_v10m) = v3dg(1,:,:,iv3d_v)
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_t2m) = v3dg(1,:,:,iv3d_t)
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_q2m) = v3dg(1,:,:,iv3d_q)

!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_rain) = [[No way]]

#ifdef H08
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_skint) = v3dg(1,:,:,iv3d_t)

  ! Assume the point where terrain height is less than 10 m is the ocean. T.Honda (02/09/2016)
  !-----------------------------------------------------------------------------

!!!$OMP PARALLEL DO PRIVATE(j,i)
  do j = 1, nlat
    do i = 1, nlon
      v2dgh(i+IHALO,j+JHALO,iv2dd_lsmask) = min(max(topo(i,j) - 10.0d0, 0.0d0), 1.0d0)
    enddo
  enddo
!!!$OMP END PARALLEL DO
#endif

  ! Pad the upper and lower halo areas
  !-----------------------------------------------------------------------------

  do iv3d = 1, nv3dd
    do j  = JS, JE
      do i  = IS, IE
        v3dgh(   1:KS-1,i,j,iv3d) = v3dgh(KS,i,j,iv3d)
        v3dgh(KE+1:KA,  i,j,iv3d) = v3dgh(KE,i,j,iv3d)
      end do
    end do
  end do

  ! Communicate the lateral halo areas
  !-----------------------------------------------------------------------------

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

  return
end subroutine state_to_history
!-------------------------------------------------------------------------------




subroutine rank_1d_2d(proc, iproc, jproc)
  use scale_rm_process, only: PRC_2Drank
  implicit none
  integer, intent(in) :: proc
  integer, intent(out) :: iproc, jproc

  iproc = PRC_2Drank(proc,1)
  jproc = PRC_2Drank(proc,2)

  return  
end subroutine rank_1d_2d


subroutine rank_2d_1d(iproc, jproc, proc)
  use scale_rm_process, only: PRC_NUM_X
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



subroutine rij_g2l(proc, ig, jg, il, jl)
  implicit none
  integer, intent(in) :: proc
  real(r_size), intent(in) :: ig
  real(r_size), intent(in) :: jg
  real(r_size), intent(out) :: il
  real(r_size), intent(out) :: jl
  integer :: iproc, jproc

  call rank_1d_2d(proc, iproc, jproc)
  il = ig - real(iproc * nlon,r_size)
  jl = jg - real(jproc * nlat,r_size)

  return  
end subroutine rij_g2l


subroutine rij_l2g(proc, il, jl, ig, jg)
  implicit none
  integer, intent(in) :: proc
  real(r_size), intent(in) :: il
  real(r_size), intent(in) :: jl
  real(r_size), intent(out) :: ig
  real(r_size), intent(out) :: jg
  integer :: iproc, jproc

  call rank_1d_2d(proc, iproc, jproc)
  ig = il + real(iproc * nlon,r_size)
  jg = jl + real(jproc * nlat,r_size)

  return  
end subroutine rij_l2g


!-----------------------------------------------------------------------
! using halo!
! proc = -1: outside the global domain
!-----------------------------------------------------------------------
SUBROUTINE rij_g2l_auto(proc,ig,jg,il,jl)
  use scale_rm_process, only: &
      PRC_NUM_X, PRC_NUM_Y
#ifdef DEBUG
  use scale_grid_index, only: &
      IHALO, JHALO, &
      IA, JA
  use scale_process, only: &
      PRC_myrank
  use scale_grid, only: &
      GRID_CX, &
      GRID_CY, &
      GRID_CXG, &
      GRID_CYG, &
      DX, &
      DY
#else
  use scale_grid_index, only: &
      IHALO, JHALO
#endif
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
  il = ig - real(iproc * nlon,r_size)
  jl = jg - real(jproc * nlat,r_size)
  call rank_2d_1d(iproc,jproc,proc)

#ifdef DEBUG
  if (PRC_myrank == proc) then
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

  RETURN
END SUBROUTINE rij_g2l_auto

subroutine read_CDF( dimlen, count_CDF, flash_CDF )
  use netcdf
  use common_ncio
  implicit none

  integer :: ncid
  integer :: dimlen
  real(r_size), intent(out) :: count_CDF(dimlen)  
  real(r_size), intent(out) :: flash_CDF(dimlen)  

  call ncio_open( trim(CDF_FP_FILENAME), NF90_NOWRITE, ncid)

  call ncio_read_1d_r8(ncid, "count", dimlen, 1, count_CDF )
  call ncio_read_1d_r8(ncid, "flash", dimlen, 1, flash_CDF )

  call ncio_close(ncid)

  return
end subroutine read_CDF

subroutine get_CDF_len( dimlen )
  use netcdf
  use common_ncio
  implicit none

  integer :: ncid
  integer, intent(out) :: dimlen
  
  call ncio_open( trim(CDF_FP_FILENAME), NF90_NOWRITE, ncid)

  call ncio_read_dim(ncid, "flash", dimlen)

  call ncio_close(ncid)

  return
end subroutine get_CDF_len

subroutine get_P( y, p )
  ! F(y) in Lien et al. 2018
  implicit none

  real(r_size), intent(in) :: y
  real(r_size), intent(out) :: p
  real(r_size) :: p0, pc
  integer :: i, imax

  imax = size(flash_CDF)

  ! climatology
  pc = count_CDF(1) / sum( count_CDF ) * 0.5d0
  p0 = pc * 0.5d0

  if ( y < flash_CDF(2) ) then
    p = p0
    return
  else
    do i = 3, imax
      if ( flash_CDF(i-1) <= y .and. y < flash_CDF(i) ) then
        p = ( sum( count_CDF(1:i-1) ) + count_CDF(i)*0.5d0 ) / sum( count_CDF )
      endif
    enddo
  endif

  p = min( p, 0.999999999d0)

  return
end subroutine get_P

subroutine F_GT_inv( p, yt )
  implicit none

  real(r_size), intent(in) :: p
  real(r_size), intent(out) :: yt
  real(r_size) :: yt_tmp

  call vderfinv( 1, 2.0d0*p - 1.0d0, yt_tmp )
  yt = dsqrt(2.0d0) * yt_tmp 

  return
end subroutine F_GT_inv

END MODULE common_scale
