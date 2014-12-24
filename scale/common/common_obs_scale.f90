MODULE common_obs_scale
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   10/04/2012 Guo-Yuan Lien     modified for GFS model
!   07/25/2014 Guo-Yuan Lien     modified for SCALE model
!
!=======================================================================
!
! [LETKF observation format]
!   (In file, all stored in single precision float)
!
!  column  description
!     (1)  variable type (1 to nid_obs; see 'id_*_obs' parameters)
!     (2)  longitude (degree)
!     (3)  latitude (degree)
!     (4)  level/height
!            u,v,t,tv,q,rh: level (hPa)
!            ps: station elevation (m)
!     (5)  observation value
!            wind (m/s)
!            temperature (K)
!            specific humidity (kg/kg)
!            relative humidity (%)
!            surface pressure (hPa)
!     (6)  observation error
!            unit same as observation value
!     (7)  observation platform type (1 to nobtype+1; see 'obtypelist' array)
!
!  --- columns below only exist in obs2 (after the observation operator processing)
!     (8)  observation time relative to analysis time (hour)
!     (9)  h(x) observation in model background
!            unit same as observation value except
!            surface pressure (Pa)
!    (10)  quality control mark (1=pass; others=do not pass)
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_scale

  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: nid_obs=8
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_tv_obs=3074
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
!
! surface observations codes > 9999
!
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_rain_obs=19999
  INTEGER,PARAMETER :: id_tclon_obs=99991  ! not used
  INTEGER,PARAMETER :: id_tclat_obs=99992  ! not used
  INTEGER,PARAMETER :: id_tcmip_obs=99993  ! not used

  INTEGER,PARAMETER :: elem_uid(nid_obs)= &
     (/id_u_obs, id_v_obs, id_t_obs, id_tv_obs, id_q_obs, id_rh_obs, &
       id_ps_obs, id_rain_obs/)
!       id_tclon_obs, id_tclat_obs, id_tcmip_obs/)

  INTEGER,PARAMETER :: nobtype = 21
  CHARACTER(6),PARAMETER :: obtypelist(nobtype)= &
     (/'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR', &
       'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG', &
       'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND', &
       'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW', &
       'TMPAPR'/)
       
  CHARACTER(3),PARAMETER :: obelmlist(nid_obs)= &
     (/'  U', '  V', '  T', ' Tv', '  Q', ' RH', ' PS', 'PRC'/)
!     (/'  U', '  V', '  T', ' Tv', '  Q', ' RH', ' PS', 'PRC', 'TCX', 'TCY', 'TCP'/)

  TYPE obs_info
    INTEGER :: nobs = 0
    INTEGER,ALLOCATABLE :: elm(:)
    REAL(r_size),ALLOCATABLE :: lon(:)
    REAL(r_size),ALLOCATABLE :: lat(:)
    REAL(r_size),ALLOCATABLE :: lev(:)
    REAL(r_size),ALLOCATABLE :: dat(:)
    REAL(r_size),ALLOCATABLE :: err(:)
    INTEGER,ALLOCATABLE :: typ(:)
    REAL(r_size),ALLOCATABLE :: dif(:)
  END TYPE obs_info

  TYPE obs_da_value
    INTEGER :: nobs = 0
    INTEGER,ALLOCATABLE :: idx(:)
    REAL(r_size),ALLOCATABLE :: val(:)
    REAL(r_size),ALLOCATABLE :: ensval(:,:)
    INTEGER,ALLOCATABLE :: qc(:)
    REAL(r_size),ALLOCATABLE :: ri(:)
    REAL(r_size),ALLOCATABLE :: rj(:)
  END TYPE obs_da_value




!  TYPE obs_csort
!    INTEGER :: nobs = 0
!    LOGICAL :: sorted = .false.
!    INTEGER,ALLOCATABLE :: idx(:)
!    INTEGER,ALLOCATABLE :: nobsgrd(:,:)
!  END TYPE obs_csort




  INTEGER,PARAMETER :: iqc_good=0
  INTEGER,PARAMETER :: iqc_gross_err=5
  INTEGER,PARAMETER :: iqc_ps_ter=10
  INTEGER,PARAMETER :: iqc_out_vhi=20
  INTEGER,PARAMETER :: iqc_out_vlo=21
  INTEGER,PARAMETER :: iqc_out_h=22
  INTEGER,PARAMETER :: iqc_otype=90
  INTEGER,PARAMETER :: iqc_time=91

CONTAINS
!-----------------------------------------------------------------------
! Transformation from model variables to an observation
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY(elm,ri,rj,rk,v3d,v2d,yobs,qc)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_size),INTENT(IN) :: v3d(nlevhalo,nlonhalo,nlathalo,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonhalo,nlathalo,nv2dd)
  REAL(r_size),INTENT(OUT) :: yobs
  INTEGER,INTENT(OUT) :: qc
  REAL(r_size) :: tg,qg
  REAL(r_size) :: qq
!  REAL(r_size) :: dummy(3)
!  INTEGER :: i,j,k
!  INTEGER :: is,ie,js,je,ks,ke
!  ie = CEILING( ri )
!  is = ie-1
!  je = CEILING( rj )
!  js = je-1
!  ke = CEILING( rk )
!  ks = ke-1

  qc = 0
  SELECT CASE (elm)
  CASE(id_u_obs)  ! U
    CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri,rj,yobs)
  CASE(id_v_obs)  ! V
    CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj,yobs)
  CASE(id_t_obs)  ! T
    CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,yobs)
  CASE(id_tv_obs)  ! Tv
    CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,yobs)
    CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,qq)
    yobs = yobs * (1.0d0 + fvirt * qq)
  CASE(id_q_obs)  ! Q
    CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,yobs)
  CASE(id_ps_obs) ! PS   !##################################

!    CALL itpl_2d(v2d(:,:,iv2dd_t2m),ri,rj,tg)
!    CALL itpl_2d(v2d(:,:,iv2dd_q2m),ri,rj,qg)
!    CALL itpl_2d(v2d(:,:,iv2dd_ps),ri,rj,yobs)
!    CALL prsadj(yobs,rk,tg,qg)

    CALL itpl_2d(v2d(:,:,iv2dd_ps),ri,rj,yobs)

!  CASE(id_rain_obs) ! RAIN                        ############# (not finished)
!    CALL itpl_2d(v2d(:,:,iv2dd_rain),ri,rj,yobs) !#############
  CASE(id_rh_obs) ! RH
    CALL itpl_3d(v3d(:,:,:,iv3dd_rh),rk,ri,rj,yobs)
!  CASE(id_tclon_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(1)
!  CASE(id_tclat_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(2)
!  CASE(id_tcmip_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(3)
  CASE DEFAULT
    yobs = undef
    qc = iqc_otype
  END SELECT

  RETURN
END SUBROUTINE Trans_XtoY
!-----------------------------------------------------------------------
! TC center search
!  [AUTHORS:] T. Miyoshi and M. Kunii
!-----------------------------------------------------------------------
!SUBROUTINE tctrk(ps,t2,ri,rj,trk)
!  IMPLICIT NONE
!  INTEGER,PARAMETER :: isearch = 10 !search radius [grid points]
!  REAL(r_size),INTENT(IN) :: ps(nlon,nlat)
!  REAL(r_size),INTENT(IN) :: t2(nlon,nlat)
!  REAL(r_size),INTENT(IN) :: ri,rj
!  REAL(r_size),INTENT(OUT) :: trk(3) !1:lon, 2:lat, 3:minp
!  REAL(r_size) :: wk(100,3)
!  REAL(r_size) :: slp(nlon,nlat)
!  REAL(r_size) :: p1,p2,p3,p4,p5,a,c,d,xx,yy
!  INTEGER :: i,j,i0,i1,j0,j1,n

!  i0 = MAX(1,FLOOR(ri)-isearch)
!  i1 = MIN(nlon,CEILING(ri)+isearch)
!  j0 = MAX(1,FLOOR(rj)-isearch)
!  j1 = MIN(nlat,CEILING(rj)+isearch)
!  trk = undef

!  DO j=j0,j1
!    DO i=i0,i1
!      slp(i,j) = ps(i,j) * (1.0d0 - 0.0065d0 * phi0(i,j) / &
!        & (t2(i,j) + 0.0065d0 * phi0(i,j))) ** -5.257d0
!    END DO
!  END DO

!  n=0
!  DO j=j0+1,j1-1
!    DO i=i0+1,i1-1
!      IF(slp(i,j) > slp(i  ,j-1)) CYCLE
!      IF(slp(i,j) > slp(i  ,j+1)) CYCLE
!      IF(slp(i,j) > slp(i-1,j  )) CYCLE
!      IF(slp(i,j) > slp(i+1,j  )) CYCLE
!      IF(slp(i,j) > slp(i-1,j-1)) CYCLE
!      IF(slp(i,j) > slp(i+1,j-1)) CYCLE
!      IF(slp(i,j) > slp(i-1,j+1)) CYCLE
!      IF(slp(i,j) > slp(i+1,j+1)) CYCLE
!      p1 = slp(i,j)
!      p2 = slp(i-1,j)
!      p3 = slp(i+1,j)
!      p4 = slp(i,j-1)
!      p5 = slp(i,j+1)
!      c = (p3-p2)*0.5d0
!      d = (p5-p4)*0.5d0
!      a = (p2+p3+p4+p5)*0.25d0 - p1
!      IF(a == 0.0d0) CYCLE
!      xx = -0.5d0 * c / a
!      yy = -0.5d0 * d / a
!      n = n+1
!      wk(n,3) = p1 - a*(xx*xx + yy*yy)
!      wk(n,2) = lat(i,j) * (1.0d0 - yy) + lat(i,j+1) * yy
!      wk(n,1) = lon(i,j) * (1.0d0 - xx) + lon(i+1,j) * xx
!    END DO
!  END DO

!  j=1
!  IF(n > 1) THEN
!    a = wk(1,3)
!    DO i=2,n
!      IF(wk(i,3) < a) THEN
!        a = wk(i,3)
!        j = i
!      END IF
!    END DO
!  END IF
!  trk = wk(j,:)

!END SUBROUTINE tctrk
!-----------------------------------------------------------------------
! Compute relative humidity (RH)
!-----------------------------------------------------------------------
SUBROUTINE calc_rh(t,q,p,rh)
  IMPLICIT NONE
  REAL(r_size),PARAMETER :: t0=273.15d0
  REAL(r_size),PARAMETER :: e0c=6.11d0
  REAL(r_size),PARAMETER :: al=17.3d0
  REAL(r_size),PARAMETER :: bl=237.3d0
  REAL(r_size),PARAMETER :: e0i=6.1121d0
  REAL(r_size),PARAMETER :: ai=22.587d0
  REAL(r_size),PARAMETER :: bi=273.86d0
  REAL(r_size),INTENT(IN) :: t,q,p
  REAL(r_size),INTENT(OUT) :: rh
  REAL(r_size) :: e,es,tc

  e = q * p * 0.01d0 / (0.378d0 * q + 0.622d0)

  tc = t-t0
  IF(tc >= 0.0d0) THEN
    es = e0c * exp(al*tc/(bl+tc))
  ELSE IF(tc <= -15.d0) THEN
    es = e0i * exp(ai*tc/(bi+tc))
  ELSE
    es = e0c * exp(al*tc/(bl+tc)) * (15.0d0+tc)/15.0d0 &
       + e0i * exp(ai*tc/(bi+tc)) * (-tc) / 15.0d0
  END IF

  rh = e/es

  RETURN
END SUBROUTINE calc_rh
!-----------------------------------------------------------------------
! Pressure adjustment for a different height level
!-----------------------------------------------------------------------
SUBROUTINE prsadj(p,dz,t,q)
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: p
  REAL(r_size),INTENT(IN) :: dz ! height difference (target - original) [m]
  REAL(r_size),INTENT(IN) :: t  ! temperature [K] at original level
  REAL(r_size),INTENT(IN) :: q  ! humidity [kg/kg] at original level
  REAL(r_size),PARAMETER :: gamma=5.0d-3 ! lapse rate [K/m]
  REAL(r_size) :: tv

  IF(dz /= 0) THEN
    tv = t * (1.0d0 + 0.608d0 * q)
    p = p * ((-gamma*dz+tv)/tv)**(gg/(gamma*rd)) !tv is at original level
!    p = p * (tv/(tv+gamma*dz))**(gg/(gamma*rd)) !tv is at target level
  END IF

  RETURN
END SUBROUTINE prsadj
!-----------------------------------------------------------------------
! Coordinate conversion
!
! rk = 0.0d0  : surface observation
! rk = -1.0d0 : too high
! rk = -2.0d0 : too low
! rk = -3.0d0 : horizontally outside (should not happen)
!-----------------------------------------------------------------------
SUBROUTINE phys2ijk(p_full,elem,ri,rj,rlev,rk)
  use scale_grid_index, only: &
      KHALO
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: p_full(nlevhalo,nlonhalo,nlathalo)
  INTEGER,INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels
  REAL(r_size),INTENT(OUT) :: rk
  REAL(r_size) :: ak
  REAL(r_size) :: lnps(nlonhalo,nlathalo)
  REAL(r_size) :: plev(nlevhalo)
  REAL(r_size) :: ptop
  INTEGER :: i,j,k, ii, jj, ks
!
! rlev -> rk
!
  if (ri < 1.0d0 .or. ri > nlonhalo .or. rj < 1.0d0 .or. rj > nlathalo) then
    write (6,'(A)') 'warning: observation is outside of the horizontal domain'
    rk = -3.0d0
    return
  end if
  !
  IF(elem > 9999) THEN ! surface observation
    rk = 0.0d0
  ELSE
    !
    ! horizontal interpolation
    !
    i = CEILING(ri)
    j = CEILING(rj)
    !
    ! Find the lowest valid level
    !
    ks = 1+KHALO
    do jj = j-1, j
      do ii = i-1, i
        DO k=1+KHALO,nlev+KHALO
          if (p_full(k,ii,jj) >= 0.0d0) exit
        END DO
        if (k > ks) ks = k
      end do
    end do
!    DO k=1,nlevhalo
    DO k=1+KHALO,nlev+KHALO
!      IF(i <= nlon+IHALO) THEN
        lnps(i-1:i,j-1:j) = LOG(p_full(k,i-1:i,j-1:j))
!      ELSE
!        lnps(i-1,j-1:j) = LOG(p_full(k,i-1,j-1:j))
!        lnps(1,j-1:j) = LOG(p_full(k,1,j-1:j))
!      END IF
      CALL itpl_2d(lnps,ri,rj,plev(k))
    END DO
    !
    ! Log pressure
    !
    rk = LOG(rlev)
    !
    ! determine if rk is within bound.
    !
    IF(rk < plev(nlev+KHALO)) THEN
      call itpl_2d(p_full(nlev+KHALO,:,:),ri,rj,ptop)
      write(6,'(A,F8.1,A,F8.1)') 'warning: observation is too high: ptop=', ptop, ', lev=', rlev
      rk = -1.0d0
      RETURN
    END IF
    IF(rk > plev(ks)) THEN
      call itpl_2d(p_full(ks,:,:),ri,rj,ptop)
      write(6,'(A,F8.1,A,F8.1)') 'warning: observation is too low: ptop=', ptop, ', lev=', rlev
      rk = -2.0d0
      RETURN
    END IF
    !
    ! find rk
    !
    DO k=ks+1,nlev+KHALO
      IF(plev(k) < rk) EXIT ! assuming descending order of plev
    END DO
    ak = (rk - plev(k-1)) / (plev(k) - plev(k-1))
    rk = REAL(k-1,r_size) + ak
  END IF

  RETURN
END SUBROUTINE phys2ijk
!-----------------------------------------------------------------------
! Coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE phys2ij(rlon,rlat,rig,rjg)
  use scale_grid, only: &
      GRID_CXG, &
      GRID_CYG, &
      DX, &
      DY
  use scale_mapproj, only: &
      MPRJ_lonlat2xy
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(OUT) :: rig
  REAL(r_size),INTENT(OUT) :: rjg
!
! rlon,rlat -> ri,rj
!
  call MPRJ_lonlat2xy(rlon*pi/180.0d0,rlat*pi/180.0d0,rig,rjg)
  rig = (rig - GRID_CXG(1)) / DX + 1.0d0
  rjg = (rjg - GRID_CYG(1)) / DY + 1.0d0

  RETURN
END SUBROUTINE phys2ij
!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlonhalo,nlathalo)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

!  IF(i <= nlon+IHALO) THEN
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj
!  ELSE
!    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
!       & + var(1  ,j-1) *    ai  * (1-aj) &
!       & + var(i-1,j  ) * (1-ai) *    aj  &
!       & + var(1  ,j  ) *    ai  *    aj
!  END IF

  RETURN
END SUBROUTINE itpl_2d

SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlevhalo,nlonhalo,nlathalo)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rk
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj,ak
  INTEGER :: i,j,k

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)
  k = CEILING(rk)
  ak = rk - REAL(k-1,r_size)

!  IF(i <= nlon+IHALO) THEN
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak
!  ELSE
!    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
!       & + var(1  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
!       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
!       & + var(1  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
!       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
!       & + var(1  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
!       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
!       & + var(1  ,j  ,k  ) *    ai  *    aj  *    ak
!  END IF

  RETURN
END SUBROUTINE itpl_3d
!-----------------------------------------------------------------------
! Monitor departure
!  ofmt: output format
!    0: U,V,T(Tv),Q,RH,PS (default)
!    1: U,V,T(Tv),Q,RH,PS,RAIN
!-----------------------------------------------------------------------
SUBROUTINE monit_dep(nn,elm,dep,qc,ofmt)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elm(nn)
  REAL(r_size),INTENT(IN) :: dep(nn)
  INTEGER,INTENT(IN) :: qc(nn)
  INTEGER,INTENT(IN),OPTIONAL :: ofmt
  INTEGER :: ofmt1
  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_q,rmse_rh,rmse_ps,rmse_rain
  REAL(r_size) :: bias_u,bias_v,bias_t,bias_q,bias_rh,bias_ps,bias_rain
  INTEGER :: n,iu,iv,it,iq,irh,ips,irain

  ofmt1 = 0
  IF(PRESENT(ofmt)) ofmt1 = ofmt

  rmse_u = 0.0d0
  rmse_v = 0.0d0
  rmse_t = 0.0d0
  rmse_q = 0.0d0
  rmse_ps = 0.0d0
  rmse_rh = 0.0d0
  rmse_rain = 0.0d0
  bias_u = 0.0d0
  bias_v = 0.0d0
  bias_t = 0.0d0
  bias_q = 0.0d0
  bias_ps = 0.0d0
  bias_rh = 0.0d0
  bias_rain = 0.0d0
  iu = 0
  iv = 0
  it = 0
  iq = 0
  ips = 0
  irh = 0
  irain = 0
  DO n=1,nn
    IF(qc(n) <= 0) CYCLE
    SELECT CASE(NINT(elm(n)))
    CASE(id_u_obs)
      rmse_u = rmse_u + dep(n)**2
      bias_u = bias_u + dep(n)
      iu = iu + 1
    CASE(id_v_obs)
      rmse_v = rmse_v + dep(n)**2
      bias_v = bias_v + dep(n)
      iv = iv + 1
    CASE(id_t_obs,id_tv_obs) ! compute T, Tv together
      rmse_t = rmse_t + dep(n)**2
      bias_t = bias_t + dep(n)
      it = it + 1
    CASE(id_q_obs)
      rmse_q = rmse_q + dep(n)**2
      bias_q = bias_q + dep(n)
      iq = iq + 1
    CASE(id_rh_obs)
      rmse_rh = rmse_rh + dep(n)**2
      bias_rh = bias_rh + dep(n)
      irh = irh + 1
    CASE(id_ps_obs)
      rmse_ps = rmse_ps + dep(n)**2
      bias_ps = bias_ps + dep(n)
      ips = ips + 1
    CASE(id_rain_obs)
      rmse_rain = rmse_rain + dep(n)**2
      bias_rain = bias_rain + dep(n)
      irain = irain + 1
    END SELECT
  END DO
  IF(iu == 0) THEN
    rmse_u = undef
    bias_u = undef
  ELSE
    rmse_u = SQRT(rmse_u / REAL(iu,r_size))
    bias_u = bias_u / REAL(iu,r_size)
  END IF
  IF(iv == 0) THEN
    rmse_v = undef
    bias_v = undef
  ELSE
    rmse_v = SQRT(rmse_v / REAL(iv,r_size))
    bias_v = bias_v / REAL(iv,r_size)
  END IF
  IF(it == 0) THEN
    rmse_t = undef
    bias_t = undef
  ELSE
    rmse_t = SQRT(rmse_t / REAL(it,r_size))
    bias_t = bias_t / REAL(it,r_size)
  END IF
  IF(iq == 0) THEN
    rmse_q = undef
    bias_q = undef
  ELSE
    rmse_q = SQRT(rmse_q / REAL(iq,r_size))
    bias_q = bias_q / REAL(iq,r_size)
  END IF
  IF(irh == 0) THEN
    rmse_rh = undef
    bias_rh = undef
  ELSE
    rmse_rh = SQRT(rmse_rh / REAL(irh,r_size))
    bias_rh = bias_rh / REAL(irh,r_size)
  END IF
  IF(ips == 0) THEN
    rmse_ps = undef
    bias_ps = undef
  ELSE
    rmse_ps = SQRT(rmse_ps / REAL(ips,r_size))
    bias_ps = bias_ps / REAL(ips,r_size)
  END IF
  IF(irain == 0) THEN
    rmse_rain = undef
    bias_rain = undef
  ELSE
    rmse_rain = SQRT(rmse_rain / REAL(irain,r_size))
    bias_rain = bias_rain / REAL(irain,r_size)
  END IF

  IF(ofmt1 == 0) THEN
    WRITE(6,'(A)') '=============================================================================='
    WRITE(6,'(6x,6A12)') 'U','V','T(Tv)','Q','RH','PS'
    WRITE(6,'(A)') '------------------------------------------------------------------------------'
    WRITE(6,'(A6,6ES12.3)') 'BIAS  ',bias_u,bias_v,bias_t,bias_q,bias_rh,bias_ps
    WRITE(6,'(A6,6ES12.3)') 'RMSE  ',rmse_u,rmse_v,rmse_t,rmse_q,rmse_rh,rmse_ps
    WRITE(6,'(A6,6I12)') 'NUMBER',iu,iv,it,iq,irh,ips
    WRITE(6,'(A)') '=============================================================================='
  ELSE IF(ofmt1 == 1) THEN
    WRITE(6,'(A)') '=========================================================================================='
    WRITE(6,'(6x,7A12)') 'U','V','T(Tv)','Q','RH','PS','RAIN'
    WRITE(6,'(A)') '------------------------------------------------------------------------------------------'
    WRITE(6,'(A6,6ES12.3)') 'BIAS  ',bias_u,bias_v,bias_t,bias_q,bias_rh,bias_ps
    WRITE(6,'(A6,6ES12.3)') 'RMSE  ',rmse_u,rmse_v,rmse_t,rmse_q,rmse_rh,rmse_ps
    WRITE(6,'(A6,7I12)') 'NUMBER',iu,iv,it,iq,irh,ips,irain
    WRITE(6,'(A)') '=========================================================================================='
  END IF

  RETURN
END SUBROUTINE monit_dep
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_info_allocate(obs)
  IMPLICIT NONE
  TYPE(obs_info),INTENT(INOUT) :: obs

  call obs_info_deallocate(obs)

  ALLOCATE( obs%elm (obs%nobs) )
  ALLOCATE( obs%lon (obs%nobs) )
  ALLOCATE( obs%lat (obs%nobs) )
  ALLOCATE( obs%lev (obs%nobs) )
  ALLOCATE( obs%dat (obs%nobs) )
  ALLOCATE( obs%err (obs%nobs) )
  ALLOCATE( obs%typ (obs%nobs) )
  ALLOCATE( obs%dif (obs%nobs) )

  obs%elm = 0
  obs%lon = 0.0d0
  obs%lat = 0.0d0
  obs%lev = 0.0d0
  obs%dat = 0.0d0
  obs%err = 0.0d0
  obs%typ = 0
  obs%dif = 0.0d0

  RETURN
END SUBROUTINE obs_info_allocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_info_deallocate(obs)
  IMPLICIT NONE
  TYPE(obs_info),INTENT(INOUT) :: obs

  IF(ALLOCATED(obs%elm)) DEALLOCATE(obs%elm)
  IF(ALLOCATED(obs%lon)) DEALLOCATE(obs%lon)
  IF(ALLOCATED(obs%lat)) DEALLOCATE(obs%lat)
  IF(ALLOCATED(obs%lev)) DEALLOCATE(obs%lev)
  IF(ALLOCATED(obs%dat)) DEALLOCATE(obs%dat)
  IF(ALLOCATED(obs%err)) DEALLOCATE(obs%err)
  IF(ALLOCATED(obs%typ)) DEALLOCATE(obs%typ)
  IF(ALLOCATED(obs%dif)) DEALLOCATE(obs%dif)

  RETURN
END SUBROUTINE obs_info_deallocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_da_value_allocate(obs,member)
  IMPLICIT NONE
  TYPE(obs_da_value),INTENT(INOUT) :: obs
  INTEGER,INTENT(IN) :: member

  call obs_da_value_deallocate(obs)

  ALLOCATE( obs%idx    (obs%nobs) )
  ALLOCATE( obs%val    (obs%nobs) )
  ALLOCATE( obs%qc     (obs%nobs) )
  ALLOCATE( obs%ri     (obs%nobs) )
  ALLOCATE( obs%rj     (obs%nobs) )

  obs%idx = 0
  obs%val = 0.0d0
  obs%qc = 0
  obs%ri = 0.0d0
  obs%rj = 0.0d0

  if (member > 0) then
    ALLOCATE( obs%ensval (member,obs%nobs) )
    obs%ensval = 0.0d0
  end if

  RETURN
END SUBROUTINE obs_da_value_allocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_da_value_deallocate(obs)
  IMPLICIT NONE
  TYPE(obs_da_value),INTENT(INOUT) :: obs

  IF(ALLOCATED(obs%idx    )) DEALLOCATE(obs%idx    )
  IF(ALLOCATED(obs%val    )) DEALLOCATE(obs%val    )
  IF(ALLOCATED(obs%ensval )) DEALLOCATE(obs%ensval )
  IF(ALLOCATED(obs%qc     )) DEALLOCATE(obs%qc     )
  IF(ALLOCATED(obs%ri     )) DEALLOCATE(obs%ri     )
  IF(ALLOCATED(obs%rj     )) DEALLOCATE(obs%rj     )

  RETURN
END SUBROUTINE obs_da_value_deallocate
!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nrec,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nrec
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl),ALLOCATABLE :: wk(:)
  INTEGER :: ios
  INTEGER :: iu,iv,it,iq,irh,ips,itc
  INTEGER :: iunit
  LOGICAL :: ex

  ALLOCATE(wk(nrec))
  nn = 0
  iu = 0
  iv = 0
  it = 0
  iq = 0
  irh = 0
  ips = 0
  itc = 0
  iunit=91
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    DO
      READ(iunit,IOSTAT=ios) wk
      IF(ios /= 0) EXIT
!      SELECT CASE(NINT(wk(1)))
!      CASE(id_u_obs)
!        iu = iu + 1
!      CASE(id_v_obs)
!        iv = iv + 1
!      CASE(id_t_obs,id_tv_obs)
!        it = it + 1
!      CASE(id_q_obs)
!        iq = iq + 1
!      CASE(id_rh_obs)
!        irh = irh + 1
!      CASE(id_ps_obs)
!        ips = ips + 1
!      CASE(id_tclon_obs)
!        itc = itc + 1
!      END SELECT
      nn = nn + 1
    END DO
!    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
!    WRITE(6,'(A12,I10)') '          U:',iu
!    WRITE(6,'(A12,I10)') '          V:',iv
!    WRITE(6,'(A12,I10)') '      T(Tv):',it
!    WRITE(6,'(A12,I10)') '          Q:',iq
!    WRITE(6,'(A12,I10)') '         RH:',irh
!    WRITE(6,'(A12,I10)') '         Ps:',ips
!    WRITE(6,'(A12,I10)') '   TC TRACK:',itc
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF
  DEALLOCATE(wk)

  RETURN
END SUBROUTINE get_nobs

SUBROUTINE read_obs(cfile,obs)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(INOUT) :: obs
  REAL(r_sngl) :: wk(8)
  INTEGER :: n,iunit

!  call obs_info_allocate(obs)

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,obs%nobs
    READ(iunit) wk
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_v_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_t_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_tv_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_q_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_ps_obs)
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = wk(6) * 100.0 ! hPa -> Pa
    CASE(id_rh_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = wk(5) * 0.01 ! percent input
      wk(6) = wk(6) * 0.01 ! percent input
    CASE(id_tcmip_obs)
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = wk(6) * 100.0 ! hPa -> Pa
    END SELECT
    obs%elm(n) = NINT(wk(1))
    obs%lon(n) = REAL(wk(2),r_size)
    obs%lat(n) = REAL(wk(3),r_size)
    obs%lev(n) = REAL(wk(4),r_size)
    obs%dat(n) = REAL(wk(5),r_size)
    obs%err(n) = REAL(wk(6),r_size)
    obs%typ(n) = NINT(wk(7))
    obs%dif(n) = REAL(wk(8),r_size)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs

SUBROUTINE write_obs(cfile,obs,append)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(IN) :: obs
  INTEGER,INTENT(IN),OPTIONAL :: append
  INTEGER :: appendr
  REAL(r_sngl) :: wk(8)
  INTEGER :: n,iunit

  iunit=92
  appendr = 0
  IF(present(append)) appendr = append

  IF(appendr == 1) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  END IF
  DO n=1,obs%nobs
    wk(1) = REAL(obs%elm(n),r_sngl)
    wk(2) = REAL(obs%lon(n),r_sngl)
    wk(3) = REAL(obs%lat(n),r_sngl)
    wk(4) = REAL(obs%lev(n),r_sngl)
    wk(5) = REAL(obs%dat(n),r_sngl)
    wk(6) = REAL(obs%err(n),r_sngl)
    wk(7) = REAL(obs%typ(n),r_sngl)
    wk(8) = REAL(obs%dif(n),r_sngl)
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    CASE(id_v_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    CASE(id_t_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    CASE(id_tv_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    CASE(id_q_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
    CASE(id_ps_obs)
      wk(5) = wk(5) * 0.01 ! Pa -> hPa
      wk(6) = wk(6) * 0.01 ! Pa -> hPa
    CASE(id_rh_obs)
      wk(4) = wk(4) * 0.01 ! Pa -> hPa
      wk(5) = wk(5) * 100.0 ! percent output
      wk(6) = wk(6) * 100.0 ! percent output
    CASE(id_tcmip_obs)
      wk(5) = wk(5) * 0.01 ! Pa -> hPa
      wk(6) = wk(6) * 0.01 ! Pa -> hPa
    END SELECT
    WRITE(iunit) wk
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs

!
! check = .false.: no check, overwrite anyway
!         .true.:  check, stop if inconsistency occurs, use maximum qc
!
SUBROUTINE read_obs_da(cfile,obs,im,check)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_da_value),INTENT(INOUT) :: obs
  INTEGER,INTENT(IN) :: im
  LOGICAL,INTENT(IN) :: check
  REAL(r_sngl) :: wk(5)
  INTEGER :: n,iunit

!  call obs_da_value_allocate(obs)

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,obs%nobs
    READ(iunit) wk
    if (check .and. obs%idx(n) /= NINT(wk(1))) then
      write (6,'(A)') 'error: obs_da_value%idx are inconsistent among the ensemble'
      stop
    end if
    obs%idx(n) = NINT(wk(1))  !!!!!! will overflow......
    if (im == 0) then
      obs%val(n) = REAL(wk(2),r_size)
    else
      obs%ensval(im,n) = REAL(wk(2),r_size)
    end if
    if ((.not. check) .or. (check .and. obs%qc(n) < NINT(wk(3)))) then ! choose the maximum qc value if check = .true.
      obs%qc(n) = NINT(wk(3))
    end if
    if (check .and. obs%ri(n) /= REAL(wk(4),r_size)) then
      write (6,'(A)') 'error: obs_da_value%ri are inconsistent among the ensemble'
      stop
    end if
    obs%ri(n) = REAL(wk(4),r_size)
    if (check .and. obs%rj(n) /= REAL(wk(5),r_size)) then
      write (6,'(A)') 'error: obs_da_value%rj are inconsistent among the ensemble'
      stop
    end if
    obs%rj(n) = REAL(wk(5),r_size)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_da

SUBROUTINE write_obs_da(cfile,obs,im,append)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_da_value),INTENT(IN) :: obs
  INTEGER,INTENT(IN) :: im
  INTEGER,INTENT(IN),OPTIONAL :: append
  INTEGER :: appendr
  REAL(r_sngl) :: wk(5)
  INTEGER :: n,iunit

  iunit=92
  appendr = 0
  IF(present(append)) appendr = append
  IF(appendr == 1) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append',STATUS='replace')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',STATUS='replace')
  END IF
  DO n=1,obs%nobs
    wk(1) = REAL(obs%idx(n),r_sngl)  !!!!!! will overflow......
    if (im == 0) then
      wk(2) = REAL(obs%val(n),r_sngl)
    else
      wk(2) = REAL(obs%ensval(im,n),r_sngl)
    end if
    wk(3) = REAL(obs%qc(n),r_sngl)
    wk(4) = REAL(obs%ri(n),r_sngl)
    wk(5) = REAL(obs%rj(n),r_sngl)
    WRITE(iunit) wk
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_da

!SUBROUTINE read_obs2(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,otyp,tdif,ohx,oqc)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: cfile
!  INTEGER,INTENT(IN) :: nn
!  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
!  REAL(r_size),INTENT(OUT) :: rlon(nn)
!  REAL(r_size),INTENT(OUT) :: rlat(nn)
!  REAL(r_size),INTENT(OUT) :: rlev(nn)
!  REAL(r_size),INTENT(OUT) :: odat(nn)
!  REAL(r_size),INTENT(OUT) :: oerr(nn)
!  REAL(r_size),INTENT(OUT) :: otyp(nn)
!  REAL(r_size),INTENT(OUT) :: tdif(nn)
!  REAL(r_size),INTENT(OUT) :: ohx(nn)
!  INTEGER,INTENT(OUT) :: oqc(nn)
!  REAL(r_sngl) :: wk(10)
!  INTEGER :: n,iunit

!  iunit=91
!  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
!  DO n=1,nn
!    READ(iunit) wk
!    SELECT CASE(NINT(wk(1)))
!    CASE(id_u_obs)
!      wk(4) = wk(4) * 100.0 ! hPa -> Pa
!    CASE(id_v_obs)
!      wk(4) = wk(4) * 100.0 ! hPa -> Pa
!    CASE(id_t_obs)
!      wk(4) = wk(4) * 100.0 ! hPa -> Pa
!    CASE(id_tv_obs)
!      wk(4) = wk(4) * 100.0 ! hPa -> Pa
!    CASE(id_q_obs)
!      wk(4) = wk(4) * 100.0 ! hPa -> Pa
!    CASE(id_ps_obs)
!      wk(5) = wk(5) * 100.0 ! hPa -> Pa
!      wk(6) = wk(6) * 100.0 ! hPa -> Pa
!    CASE(id_rh_obs)
!      wk(4) = wk(4) * 100.0 ! hPa -> Pa
!      wk(5) = wk(5) * 0.01 ! percent input
!      wk(6) = wk(6) * 0.01 ! percent input
!    CASE(id_tcmip_obs)
!      wk(5) = wk(5) * 100.0 ! hPa -> Pa
!      wk(6) = wk(6) * 100.0 ! hPa -> Pa
!    END SELECT
!    elem(n) = REAL(wk(1),r_size)
!    rlon(n) = REAL(wk(2),r_size)
!    rlat(n) = REAL(wk(3),r_size)
!    rlev(n) = REAL(wk(4),r_size)
!    odat(n) = REAL(wk(5),r_size)
!    oerr(n) = REAL(wk(6),r_size)
!    otyp(n) = REAL(wk(7),r_size)
!    tdif(n) = REAL(wk(8),r_size)
!    ohx(n) = REAL(wk(9),r_size)
!    oqc(n) = NINT(wk(10))
!  END DO
!  CLOSE(iunit)

!  RETURN
!END SUBROUTINE read_obs2

!SUBROUTINE write_obs2(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,otyp,tdif,ohx,oqc,append)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: cfile
!  INTEGER,INTENT(IN) :: nn
!  REAL(r_size),INTENT(IN) :: elem(nn) ! element number
!  REAL(r_size),INTENT(IN) :: rlon(nn)
!  REAL(r_size),INTENT(IN) :: rlat(nn)
!  REAL(r_size),INTENT(IN) :: rlev(nn)
!  REAL(r_size),INTENT(IN) :: odat(nn)
!  REAL(r_size),INTENT(IN) :: oerr(nn)
!  REAL(r_size),INTENT(IN) :: otyp(nn)
!  REAL(r_size),INTENT(IN) :: tdif(nn)
!  REAL(r_size),INTENT(IN) :: ohx(nn)
!  INTEGER,INTENT(IN) :: oqc(nn)
!  INTEGER,INTENT(IN) :: append
!  REAL(r_sngl) :: wk(10)
!  INTEGER :: n,iunit

!  iunit=92
!  IF(append == 0) THEN
!    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
!  ELSE
!    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append')
!  END IF
!  DO n=1,nn
!    wk(1) = REAL(elem(n),r_sngl)
!    wk(2) = REAL(rlon(n),r_sngl)
!    wk(3) = REAL(rlat(n),r_sngl)
!    wk(4) = REAL(rlev(n),r_sngl)
!    wk(5) = REAL(odat(n),r_sngl)
!    wk(6) = REAL(oerr(n),r_sngl)
!    wk(7) = REAL(otyp(n),r_sngl)
!    wk(8) = REAL(tdif(n),r_sngl)
!    wk(9) = REAL(ohx(n),r_sngl)
!    wk(10) = REAL(oqc(n),r_sngl)
!    SELECT CASE(NINT(wk(1)))
!    CASE(id_u_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!    CASE(id_v_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!    CASE(id_t_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!    CASE(id_tv_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!    CASE(id_q_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!    CASE(id_ps_obs)
!      wk(5) = wk(5) * 0.01 ! Pa -> hPa
!      wk(6) = wk(6) * 0.01 ! Pa -> hPa
!    CASE(id_rh_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!      wk(5) = wk(5) * 100.0 ! percent output
!      wk(6) = wk(6) * 100.0 ! percent output
!    CASE(id_tcmip_obs)
!      wk(5) = wk(5) * 0.01 ! Pa -> hPa
!      wk(6) = wk(6) * 0.01 ! Pa -> hPa
!    END SELECT
!    WRITE(iunit) wk
!  END DO
!  CLOSE(iunit)

!  RETURN
!END SUBROUTINE write_obs2

FUNCTION uid_obs(id_obs)
  IMPLICIT NONE
  INTEGER :: id_obs
  INTEGER :: uid_obs

  SELECT CASE(id_obs)
  CASE(id_u_obs)
    uid_obs = 1
  CASE(id_v_obs)
    uid_obs = 2
  CASE(id_t_obs)
    uid_obs = 3
  CASE(id_tv_obs)
    uid_obs = 4
  CASE(id_q_obs)
    uid_obs = 5
  CASE(id_rh_obs)
    uid_obs = 6
  CASE(id_ps_obs)
    uid_obs = 7
  CASE(id_rain_obs)
    uid_obs = 8
!  CASE(id_tclon_obs)
!    uid_obs = 9
!  CASE(id_tclat_obs)
!    uid_obs = 10
!  CASE(id_tcmip_obs)
!    uid_obs = 11
  CASE DEFAULT
    uid_obs = -1 ! error
  END SELECT
END FUNCTION uid_obs

END MODULE common_obs_scale
