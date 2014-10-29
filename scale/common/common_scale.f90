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
  USE common_ncio
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: nlonsub=24
  INTEGER,PARAMETER :: nlatsub=24
  INTEGER,PARAMETER :: nlonns=2
  INTEGER,PARAMETER :: nlatns=2
  INTEGER,PARAMETER :: nlon=nlonsub*nlonns
  INTEGER,PARAMETER :: nlat=nlatsub*nlatns
  INTEGER,PARAMETER :: nlev=30
  INTEGER,PARAMETER :: nv3d=11   ! 3D state variables (in SCALE init/restart files)
  INTEGER,PARAMETER :: nv3dd=12  ! 3D diagnostic variables (in SCALE history files)
  INTEGER,PARAMETER :: nv3dx=12  ! 3D diagnostic variables
  INTEGER,PARAMETER :: nv2d=0    ! 2D state variables (in SCALE init/restart files)
  INTEGER,PARAMETER :: nv2dd=3   ! 2D diagnostic variables (in SCALE history files)
  INTEGER,PARAMETER :: nv2dx=3   ! 2D diagnostic variables
  INTEGER,PARAMETER :: iv3d_rho=1
  INTEGER,PARAMETER :: iv3d_rhou=2
  INTEGER,PARAMETER :: iv3d_rhov=3
  INTEGER,PARAMETER :: iv3d_rhow=4
  INTEGER,PARAMETER :: iv3d_rhot=5
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
  INTEGER,PARAMETER :: iv2dd_hgt=1
  INTEGER,PARAMETER :: iv2dd_tsfc=2
  INTEGER,PARAMETER :: iv2dd_rain=3
  INTEGER,PARAMETER :: nij0=nlon*nlat
  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,PARAMETER :: nlevalld=nlev*nv3dd+nv2dd
  INTEGER,PARAMETER :: nlevallx=nlev*nv3dx+nv2dx
  INTEGER,PARAMETER :: ngpv=nij0*nlevall
  INTEGER,PARAMETER :: ngpvd=nij0*nlevalld
  INTEGER,PARAMETER :: ngpvx=nij0*nlevallx
  REAL(r_size),SAVE :: lon(nlon,nlat)
  REAL(r_size),SAVE :: lat(nlon,nlat)
  REAL(r_size),SAVE :: lonu(nlon,nlat)
  REAL(r_size),SAVE :: latu(nlon,nlat)
  REAL(r_size),SAVE :: lonv(nlon,nlat)
  REAL(r_size),SAVE :: latv(nlon,nlat)
!  REAL(r_size),SAVE :: dx(nlat)
!  REAL(r_size),SAVE :: dy(nlat)
!  REAL(r_size),SAVE :: dy2(nlat)
  REAL(r_size),SAVE :: fcori(nlat)
!  REAL(r_size),SAVE :: wg(nlon,nlat)
  INTEGER,PARAMETER :: vname_max = 10
  CHARACTER(vname_max),SAVE :: v3d_name(nv3d)
  CHARACTER(vname_max),SAVE :: v3dd_name(nv3dx)
  CHARACTER(vname_max),SAVE :: v2d_name(nv2d)
  CHARACTER(vname_max),SAVE :: v2dd_name(nv2dx)

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
  ! state variables (for LETKF)
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
  !
  ! diagnostic variables (for observation operators)
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
  !
  v2dd_name(iv2dd_tsfc) = 'SFC_TEMP'
  v2dd_name(iv2dd_rain) = 'PREC'
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
! File I/O
!-----------------------------------------------------------------------
!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_grd(filename,inv3d,inv2d,iv3dname,iv2dname,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  INTEGER,INTENT(IN) :: inv3d,inv2d
  CHARACTER(vname_max),INTENT(IN) :: iv3dname(inv3d)
  CHARACTER(vname_max),INTENT(IN) :: iv2dname(inv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nlonsub,nlatsub,nlev,inv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlonsub,nlatsub,inv2d)
  REAL(r_dble) :: buf3d8(nlonsub,nlatsub,nlev)
  REAL(r_dble) :: buf2d8(nlonsub,nlatsub)
  INTEGER :: iunit,n

  CALL ncio_open(filename, nf90_nowrite, iunit)
  DO n=1,inv3d
    IF(r_size == r_dble) THEN
      CALL ncio_read_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,v3d(:,:,:,n))
    ELSE
      CALL ncio_read_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,buf3d8)
      v3d(:,:,:,n) = REAL(buf3d8,r_size)
    END IF
  END DO
  DO n=1,inv2d
    IF(r_size == r_dble) THEN
      CALL ncio_read_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,v2d(:,:,n))
    ELSE
      CALL ncio_read_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,buf2d8)
      v2d(:,:,n) = REAL(buf2d8,r_size)
    END IF
  END DO
  CALL ncio_close(iunit)

  RETURN
END SUBROUTINE read_grd

SUBROUTINE read_grd4(filename,inv3d,inv2d,iv3dname,iv2dname,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  INTEGER,INTENT(IN) :: inv3d,inv2d
  CHARACTER(vname_max),INTENT(IN) :: iv3dname(inv3d)
  CHARACTER(vname_max),INTENT(IN) :: iv2dname(inv2d)
  REAL(r_sngl),INTENT(OUT) :: v3d(nlonsub,nlatsub,nlev,inv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlonsub,nlatsub,inv2d)
  REAL(r_dble) :: buf3d8(nlonsub,nlatsub,nlev)
  REAL(r_dble) :: buf2d8(nlonsub,nlatsub)
  INTEGER :: iunit,n

  CALL ncio_open(filename, nf90_nowrite, iunit)
  DO n=1,inv3d
    CALL ncio_read_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,buf3d8)
    v3d(:,:,:,n) = REAL(buf3d8,r_sngl)
  END DO
  DO n=1,inv2d
    CALL ncio_read_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,buf2d8)
    v2d(:,:,n) = REAL(buf2d8,r_sngl)
  END DO
  CALL ncio_close(iunit)

  RETURN
END SUBROUTINE read_grd4

!-- Write a grid file -------------------------------------------------
SUBROUTINE write_grd(filename,inv3d,inv2d,iv3dname,iv2dname,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  INTEGER,INTENT(IN) :: inv3d,inv2d
  CHARACTER(vname_max),INTENT(IN) :: iv3dname(inv3d)
  CHARACTER(vname_max),INTENT(IN) :: iv2dname(inv2d)
  REAL(r_size),INTENT(IN) :: v3d(nlonsub,nlatsub,nlev,inv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlonsub,nlatsub,inv2d)
  REAL(r_dble) :: buf3d8(nlonsub,nlatsub,nlev)
  REAL(r_dble) :: buf2d8(nlonsub,nlatsub)
  INTEGER :: iunit,n

  CALL ncio_open(filename, nf90_write, iunit)
  DO n=1,inv3d
    IF(r_size == r_dble) THEN
      CALL ncio_write_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,v3d(:,:,:,n))
    ELSE
      buf3d8 = REAL(v3d(:,:,:,n),r_dble)
      CALL ncio_write_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,buf3d8)
    END IF
  END DO
  DO n=1,inv2d
    IF(r_size == r_dble) THEN
      CALL ncio_write_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,v2d(:,:,n))
    ELSE
      buf2d8 = REAL(v2d(:,:,n),r_dble)
      CALL ncio_write_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,buf2d8)
    END IF
  END DO
  CALL ncio_close(iunit)

  RETURN
END SUBROUTINE write_grd

SUBROUTINE write_grd4(filename,inv3d,inv2d,iv3dname,iv2dname,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  INTEGER,INTENT(IN) :: inv3d,inv2d
  CHARACTER(vname_max),INTENT(IN) :: iv3dname(inv3d)
  CHARACTER(vname_max),INTENT(IN) :: iv2dname(inv2d)
  REAL(r_sngl),INTENT(IN) :: v3d(nlonsub,nlatsub,nlev,inv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlonsub,nlatsub,inv2d)
  REAL(r_dble) :: buf3d8(nlonsub,nlatsub,nlev)
  REAL(r_dble) :: buf2d8(nlonsub,nlatsub)
  INTEGER :: iunit,n

  CALL ncio_open(filename, nf90_write, iunit)
  DO n=1,inv3d
    buf3d8 = REAL(v3d(:,:,:,n),r_dble)
    CALL ncio_write_3d_r8(iunit,trim(iv3dname(n)),nlonsub,nlatsub,nlev,1,buf3d8)
  END DO
  DO n=1,inv2d
    buf2d8 = REAL(v2d(:,:,n),r_dble)
    CALL ncio_write_2d_r8(iunit,trim(iv2dname(n)),nlonsub,nlatsub,1,buf2d8)
  END DO
  CALL ncio_close(iunit)

  RETURN
END SUBROUTINE write_grd4
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
!SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
!  IMPLICIT NONE
!  INTEGER,INTENT(IN) :: member
!  INTEGER,INTENT(IN) :: nij
!  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
!  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
!  INTEGER :: i,k,m,n

!  DO n=1,nv3d
!!$OMP PARALLEL DO PRIVATE(i,k,m)
!    DO k=1,nlev
!      DO i=1,nij
!        v3dm(i,k,n) = v3d(i,k,1,n)
!        DO m=2,member
!          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
!        END DO
!        v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
!      END DO
!    END DO
!!$OMP END PARALLEL DO
!  END DO

!  DO n=1,nv2d
!    DO i=1,nij
!      v2dm(i,n) = v2d(i,1,n)
!      DO m=2,member
!        v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
!      END DO
!      v2dm(i,n) = v2dm(i,n) / REAL(member,r_size)
!    END DO
!  END DO

!  RETURN
!END SUBROUTINE ensmean_grd

END MODULE common_scale
