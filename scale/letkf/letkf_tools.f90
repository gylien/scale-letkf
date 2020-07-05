MODULE letkf_tools
!=======================================================================
!
! [PURPOSE:] Module for LETKF with GFS
!
! [HISTORY:]
!   01/26/2009   Takemasa Miyoshi  created
!   10/04/2012   Guo-Yuan Lien     modified for GFS model
!   07/01/2013   Daisuke Hotta     ported EFSO code from Y.Ohta's code
!   01/01/2014   Guo-Yuan Lien     merged to GFS-LETKF main development
!   October 2014 Guo-Yuan Lien     modified for SCALE model
!   ............ See git history for the following revisions
!
!=======================================================================
!$USE OMP_LIB
  USE common
  use common_nml
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
  USE common_letkf

  USE letkf_obs
!  USE efso_tools

  use scale_precision, only: RP
#ifdef PNETCDF
  use scale_file, only: FILE_AGGREGATE
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: das_letkf !, das_efso

  real(r_size),save :: var_local(nv3d+nv2d,nid_obs_varlocal)

  integer,save :: ctype_merge(nid_obs,nobtype)
  integer,allocatable,save :: n_merge(:)
  integer,allocatable,save :: ic_merge(:,:)
  integer,save :: n_merge_max
  logical,save :: radar_only

  integer,parameter :: n_search_incr = 8

CONTAINS
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE das_letkf(gues3d,gues2d,anal3d,anal2d)
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  use common_rand
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: gues3d(nij1,nlev,nens,nv3d) ! background ensemble
  REAL(r_size),INTENT(INOUT) :: gues2d(nij1,nens,nv2d)      !  output: destroyed
  REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nens,nv3d)   ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal2d(nij1,nens,nv2d)

!  REAL(r_size) :: mean3d(nij1,nlev,nv3d)
!  REAL(r_size) :: mean2d(nij1,nv2d)
  REAL(r_size) :: work3d(nij1,nlev,nv3d)
  REAL(r_size) :: work2d(nij1,nv2d)
  REAL(r_size),ALLOCATABLE :: work3da(:,:,:)     !GYL
  REAL(r_size),ALLOCATABLE :: work2da(:,:)       !GYL
  REAL(r_size),ALLOCATABLE :: work3dn(:,:,:,:)   !GYL
  REAL(r_size),ALLOCATABLE :: work2dn(:,:,:)     !GYL
  REAL(RP),ALLOCATABLE :: work3dg(:,:,:,:)
  REAL(RP),ALLOCATABLE :: work2dg(:,:,:)

  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size),ALLOCATABLE :: depd(:)            !GYL

  integer :: var_local_n2nc_max
  integer :: var_local_n2nc(nv3d+nv2d)
  integer :: var_local_n2n(nv3d+nv2d)
  logical :: var_local_n2n_found
  integer :: n2n, n2nc

  real(r_size) :: parm
  real(r_size),allocatable :: trans(:,:,:)
  real(r_size),allocatable :: transm(:,:)
  real(r_size),allocatable :: transmd(:,:)
  real(r_size),allocatable :: pa(:,:,:)
  real(r_size) :: transrlx(MEMBER,MEMBER)
  logical :: trans_done(nv3d+nv2d)

  INTEGER :: ij,ilev,n,m,i,k,nobsl
  INTEGER :: nobsl_t(nid_obs,nobtype)            !GYL
  REAL(r_size) :: cutd_t(nid_obs,nobtype)        !GYL
  REAL(r_size) :: beta                           !GYL
  REAL(r_size) :: tmpinfl                        !GYL
  REAL(r_size) :: q_mean,q_sprd                  !GYL
  REAL(r_size) :: q_anal(MEMBER)                 !GYL

  INTEGER :: mshuf,ierr                          !GYL
  INTEGER :: ishuf(MEMBER)                       !GYL
  real(r_size), allocatable :: addinfl_weight(:) !GYL
  real(r_size) :: rdx,rdy,rdxy,ref_min_dist      !GYL
  integer :: ic,ic2,iob                          !GYL

  integer,allocatable :: search_q0(:,:,:)

  integer :: OMP_GET_NUM_THREADS, omp_chunk

  character(len=timer_name_width) :: timer_str

  call mpi_timer('', 2)

  WRITE(6,'(A)') 'Hello from das_letkf'
  WRITE(6,'(A,F15.2)') '  INFL_MUL = ',INFL_MUL

  WRITE(6,'(A,I8)') 'Target observation numbers (global) : NOBS=',nobstotalg
  WRITE(6,'(A,I8)') 'Target observation numbers processed in this subdomian : NOBS=',nobstotal
!!  !
!!  ! In case of no obs
!!  !
!!  IF(nobstotal == 0) THEN
!!    WRITE(6,'(A)') 'No observation assimilated'
!!    anal3d = gues3d
!!    anal2d = gues2d
!!    RETURN
!!  END IF
  !
  ! Variable localization
  !
  var_local(:,1) = VAR_LOCAL_UV(:)
  var_local(:,2) = VAR_LOCAL_T(:)
  var_local(:,3) = VAR_LOCAL_Q(:)
  var_local(:,4) = VAR_LOCAL_PS(:)
  var_local(:,5) = VAR_LOCAL_RAIN(:)
  var_local(:,6) = VAR_LOCAL_TC(:)
  var_local(:,7) = VAR_LOCAL_RADAR_REF(:)
  var_local(:,8) = VAR_LOCAL_RADAR_VR(:)
  var_local(:,9) = VAR_LOCAL_H08(:) ! H08
  var_local_n2nc_max = 1
  var_local_n2nc(1) = 1
  var_local_n2n(1) = 1
  do n = 2, nv3d+nv2d
    var_local_n2n_found = .false.
    do i = 1, var_local_n2nc_max
      if (maxval(abs(var_local(var_local_n2nc(i),:) - var_local(n,:))) < tiny(var_local)) then
        var_local_n2nc(n) = var_local_n2nc(i)
        var_local_n2n(n) = var_local_n2n(var_local_n2nc(n))
        var_local_n2n_found = .true.
        exit
      end if
    end do
    if (.not. var_local_n2n_found) then
      var_local_n2nc_max = var_local_n2nc_max + 1
      var_local_n2nc(n) = var_local_n2nc_max
      var_local_n2n(n) = n
    end if
  end do
  if (LOG_LEVEL >= 2) then
    write (6, '(A,I3)') '[Info] var_local_n2nc_max=', var_local_n2nc_max
    do n = 1, nv3d+nv2d
      write (6, '(A,I3,A,I3,A,I3,A,I3)') '[Info] var_local_n2n(', n, ')=', var_local_n2n(n), '; var_local_n2nc(', n, ')=', var_local_n2nc(n)
    end do
  end if
  !
  ! Observation number limit (*to be moved to namelist*)
  !
  ctype_merge(:,:) = 0
  ctype_merge(uid_obs(id_radar_ref_obs),22) = 1
  ctype_merge(uid_obs(id_radar_ref_zero_obs),22) = 1

  allocate (n_merge(nctype))
  allocate (ic_merge(nid_obs*nobtype,nctype))
  n_merge(:) = 1
  do ic = 1, nctype
    if (n_merge(ic) > 0) then
      ic_merge(1,ic) = ic
      if (ctype_merge(elm_u_ctype(ic),typ_ctype(ic)) > 0) then
        do ic2 = ic+1, nctype
          if (ctype_merge(elm_u_ctype(ic2),typ_ctype(ic2)) == ctype_merge(elm_u_ctype(ic),typ_ctype(ic))) then
            n_merge(ic) = n_merge(ic) + 1
            ic_merge(n_merge(ic),ic) = ic2
            n_merge(ic2) = 0
            if (LOG_LEVEL >= 2) then
              write(6, '(9A)') '[Info] Observation number limit: Consider obs types (', obtypelist(typ_ctype(ic)), ', ', obelmlist(elm_u_ctype(ic)), &
                               ') and (', obtypelist(typ_ctype(ic2)), ', ', obelmlist(elm_u_ctype(ic2)), ') together'
            end if
          end if
        end do
      end if ! [ ctype_merge(elm_u_ctype(ic),typ_ctype(ic)) > 0 ]
    end if ! [ n_merge(ic) > 0 ]
  end do ! [ ic = 1, nctype ]
  n_merge_max = maxval(n_merge)

  allocate (search_q0(nctype,nv3d+1,nij1))
  search_q0(:,:,:) = 1
  !
  radar_only = .true.
  do ic = 1, nctype
    if (obtypelist(typ_ctype(ic)) /= 'PHARAD') then
      radar_only = .false.
      exit
    end if
  end do
  !
  ! FCST PERTURBATIONS
  !
!  .... this has been done by write_ensmean in letkf.f90
!  CALL ensmean_grd(MEMBER,nens,nij1,gues3d,gues2d,mean3d,mean2d)
!$OMP PARALLEL PRIVATE(n,m,k,i)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  DO n=1,nv3d
    DO m=1,MEMBER
      DO k=1,nlev
        DO i=1,nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - gues3d(i,k,mmean,n)
        END DO
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  DO n=1,nv2d
    DO m=1,MEMBER
      DO i=1,nij1
        gues2d(i,m,n) = gues2d(i,m,n) - gues2d(i,mmean,n)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  call mpi_timer('das_letkf:fcst_perturbation:', 2)

  !
  ! multiplicative inflation
  !
  IF(INFL_MUL > 0.0d0) THEN  ! fixed multiplicative inflation parameter
    work3d = INFL_MUL
    work2d = INFL_MUL
  ELSE  ! 3D parameter values are read-in
    allocate (work3dg(nlon,nlat,nlev,nv3d))
    allocate (work2dg(nlon,nlat,nv2d))
    IF(myrank_e == mmean_rank_e) THEN
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',INFL_MUL_IN_BASENAME,'.pe',myrank_d,'.nc'
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call read_restart_par(INFL_MUL_IN_BASENAME,work3dg,work2dg,MPI_COMM_d)
      else
#endif
        call read_restart(INFL_MUL_IN_BASENAME,work3dg,work2dg)
#ifdef PNETCDF
      end if
#endif

      call mpi_timer('das_letkf:adaptive_infl_read_restart:', 2)
    END IF

    call mpi_timer('', 2, barrier=MPI_COMM_e)

    CALL scatter_grd_mpi(mmean_rank_e,work3dg,work2dg,work3d,work2d)

    call mpi_timer('das_letkf:adaptive_infl_scatter:', 2)
  END IF
  IF(INFL_MUL_MIN > 0.0d0) THEN
    work3d = max(work3d, INFL_MUL_MIN)
    work2d = max(work2d, INFL_MUL_MIN)
  END IF
  !
  ! RTPS relaxation: inflation output
  !
  IF(RELAX_SPREAD_OUT) THEN
    allocate (work3da(nij1,nlev,nv3d))
    allocate (work2da(nij1,nv2d))
    work3da = 1.0d0
    work2da = 1.0d0
  END IF
  !
  ! NOBS output
  !
  IF(NOBS_OUT) THEN
    allocate (work3dn(nobtype+6,nij1,nlev,nv3d))
    allocate (work2dn(nobtype,nij1,nv2d))
    work3dn = 0.0d0
    work2dn = 0.0d0
  END IF

  call mpi_timer('das_letkf:allocation_shared_vars:', 2)

!$OMP PARALLEL PRIVATE(ilev,ij,n,m,k,hdxf,rdiag,rloc,dep,depd,nobsl,nobsl_t,cutd_t,parm,beta,n2n,n2nc,trans,transm,transmd,transrlx,pa,trans_done,tmpinfl,q_mean,q_sprd,q_anal,timer_str)
  omp_chunk = min(4, max(1, (nij1-1) / OMP_GET_NUM_THREADS() + 1))

  allocate (hdxf (nobstotal,MEMBER))
  allocate (rdiag(nobstotal))
  allocate (rloc (nobstotal))
  allocate (dep  (nobstotal))
  if (DET_RUN) then
    allocate (depd (nobstotal))
  end if
  allocate (trans  (MEMBER,MEMBER,var_local_n2nc_max))
  allocate (transm (MEMBER,       var_local_n2nc_max))
  allocate (transmd(MEMBER,       var_local_n2nc_max))
  allocate (pa     (MEMBER,MEMBER,var_local_n2nc_max))


!$OMP MASTER
  call mpi_timer('das_letkf:allocation_private_vars:', 2)
  call mpi_timer('', 3)

  write (6, '(A,I3)') 'OpenMP chunk for dynamic schedule =', omp_chunk
!$OMP END MASTER
  !
  ! MAIN ASSIMILATION LOOP
  !
  DO ilev=1,nlev
    if (LOG_LEVEL >= 3) then
      call mpi_timer('', 4)
    end if

!$OMP DO SCHEDULE(DYNAMIC,omp_chunk)
    DO ij=1,nij1

      trans_done(:) = .false.                                                          !GYL

      ! weight parameter based on grid locations (not for covariance inflation purpose)
      ! if the weight is zero, no analysis update is needed
      call relax_beta(rig1(ij),rjg1(ij),hgt1(ij,ilev),beta)

      if (LOG_LEVEL >= 3) then
        write (timer_str, '(A25,I4,A4,I4,A2)') 'das_letkf:relax_beta(lev=', ilev, ',ij=', ij, '):'
        call mpi_timer(trim(timer_str), 4)
      end if

      if (beta == 0.0d0) then
        do n = 1, nv3d
          do m = 1, MEMBER
            anal3d(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n) + gues3d(ij,ilev,m,n)
          end do
          if (DET_RUN) then
            anal3d(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n)
          end if
        end do
        if (ilev == 1) then
          do n = 1, nv2d
            do m = 1, MEMBER
              anal2d(ij,m,n) = gues2d(ij,mmean,n) + gues2d(ij,m,n)
            end do
            if (DET_RUN) then
              anal2d(ij,mmdet,n) = gues2d(ij,mmdet,n)
            end if
          end do
        end if

        if (LOG_LEVEL >= 3) then
          write (timer_str, '(A30,I4,A4,I4,A7)') 'das_letkf:letkf_core_skip(lev=', ilev, ',ij=', ij, ',allv):'
          call mpi_timer(trim(timer_str), 4)
        end if

        cycle
      end if

      if (LOG_LEVEL >= 3) then
        call mpi_timer('', 5)
      end if

      ! update 3D variables
      DO n=1,nv3d

        n2nc = var_local_n2nc(n)
        n2n = var_local_n2n(n)

        if (gues3d(ij,ilev,mmean,iv3d_p) < Q_UPDATE_TOP .and. n >= iv3d_q .and. n <= iv3d_qg) then !GYL - Upper bound of Q update levels
          do m = 1, MEMBER                                                             !GYL
            anal3d(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n) + gues3d(ij,ilev,m,n)        !GYL
          end do                                                                       !GYL
          if (DET_RUN) then                                                            !GYL
            anal3d(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n)                          !GYL
          end if                                                                       !GYL

          if (LOG_LEVEL >= 3) then
            write (timer_str, '(A30,I4,A4,I4,A3,I2,A2)') 'das_letkf:letkf_core_skip(lev=', ilev, ',ij=', ij, ',n=', n, '):'
            call mpi_timer(trim(timer_str), 5)
          end if

          cycle                                                                        !GYL
        end if                                                                         !GYL

        if (RELAX_TO_INFLATED_PRIOR) then
          parm = work3d(ij,ilev,n)
        else
          parm = 1.0d0
        end if

        ! calculate mean and perturbation weights
        if (trans_done(n2nc)) then
          ! if weights already computed for other variables can be re-used(no variable localization), do not need to compute again
          if (INFL_MUL_ADAPTIVE) then
            work3d(ij,ilev,n) = work3d(ij,ilev,n2n)
          end if
          if (NOBS_OUT) then
            work3dn(:,ij,ilev,n) = work3dn(:,ij,ilev,n2n)
          end if

          if (LOG_LEVEL >= 3) then
            write (timer_str, '(A30,I4,A4,I4,A3,I2,A2)') 'das_letkf:letkf_core_copy(lev=', ilev, ',ij=', ij, ',n=', n, '):'
            call mpi_timer(trim(timer_str), 5)
          end if
        ELSE
          ! compute weights with localized observations
          if (DET_RUN) then                                                            !GYL
            CALL obs_local(rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),n, & !GYL
                           hdxf,rdiag,rloc,dep,nobsl,depd=depd,nobsl_t=nobsl_t,cutd_t=cutd_t,srch_q0=search_q0(:,n,ij)) !GYL
          else                                                                         !GYL
            CALL obs_local(rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),n, & !GYL
                           hdxf,rdiag,rloc,dep,nobsl,nobsl_t=nobsl_t,cutd_t=cutd_t,srch_q0=search_q0(:,n,ij)) !GYL
          end if                                                                       !GYL
          IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                         !GYL
            if (DET_RUN) then                                                          !GYL
              CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & !GYL
                              trans(:,:,n2nc),transm=transm(:,n2nc),pao=pa(:,:,n2nc), & !GYL
                              rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE, &       !GYL
                              depd=depd,transmd=transmd(:,n2nc))                       !GYL
            else                                                                       !GYL
              CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & !GYL
                              trans(:,:,n2nc),transm=transm(:,n2nc),pao=pa(:,:,n2nc), & !GYL
                              rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)         !GYL
            end if                                                                     !GYL
          ELSE                                                                         !GYL
            if (DET_RUN) then                                                          !GYL
              CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & !GYL
                              trans(:,:,n2nc),transm=transm(:,n2nc),           &       !GYL
                              rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE, &       !GYL
                              depd=depd,transmd=transmd(:,n2nc))                       !GYL
            else                                                                       !GYL
              CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & !GYL
                              trans(:,:,n2nc),transm=transm(:,n2nc),           &       !GYL
                              rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)         !GYL
            end if                                                                     !GYL
          END IF                                                                       !GYL
          trans_done(n2nc) = .true.                                                    !GYL
          IF(NOBS_OUT) THEN                                                            !GYL
            work3dn(:,ij,ilev,n) = real(sum(nobsl_t, dim=1),r_size)                    !GYL !!! NOBS: sum over all variables for each report type
            work3dn(nobtype+1,ij,ilev,n) = real(nobsl_t(9,22),r_size)                  !GYL !!! NOBS: ref
            work3dn(nobtype+2,ij,ilev,n) = real(nobsl_t(10,22),r_size)                 !GYL !!! NOBS: re0
            work3dn(nobtype+3,ij,ilev,n) = real(nobsl_t(11,22),r_size)                 !GYL !!! NOBS: vr
            work3dn(nobtype+4,ij,ilev,n) = real(cutd_t(9,22),r_size)                   !GYL !!! CUTOFF_DIST: ref
            work3dn(nobtype+5,ij,ilev,n) = real(cutd_t(10,22),r_size)                  !GYL !!! CUTOFF_DIST: re0
            work3dn(nobtype+6,ij,ilev,n) = real(cutd_t(11,22),r_size)                  !GYL !!! CUTOFF_DIST: vr
          END IF                                                                       !GYL

          if (LOG_LEVEL >= 3) then
            write (timer_str, '(A30,I4,A4,I4,A3,I2,A2)') 'das_letkf:letkf_core_calc(lev=', ilev, ',ij=', ij, ',n=', n, '):'
            call mpi_timer(trim(timer_str), 5)
          end if
        END IF

        ! relaxation via LETKF weight
        IF(RELAX_ALPHA /= 0.0d0) THEN                                                  !GYL - RTPP method (Zhang et al. 2004)
          CALL weight_RTPP(trans(:,:,n2nc),parm,transrlx)                              !GYL
        ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                      !GYL - RTPS method (Whitaker and Hamill 2012)
          IF(RELAX_SPREAD_OUT) THEN                                                    !GYL
            CALL weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues3d(ij,ilev,:,n), &       !GYL
                             parm,transrlx,work3da(ij,ilev,n))                         !GYL
          ELSE                                                                         !GYL
            CALL weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues3d(ij,ilev,:,n), &       !GYL
                             parm,transrlx,tmpinfl)                                    !GYL
          END IF                                                                       !GYL
        ELSE                                                                           !GYL
          transrlx = trans(:,:,n2nc)                                                   !GYL - No relaxation
        END IF                                                                         !GYL

        ! total weight matrix
        DO m=1,MEMBER                                                                  !GYL
          DO k=1,MEMBER                                                                !GYL
            transrlx(k,m) = (transrlx(k,m) + transm(k,n2nc)) * beta                    !GYL
          END DO                                                                       !GYL
          transrlx(m,m) = transrlx(m,m) + (1.0d0-beta)                                 !GYL
        END DO                                                                         !GYL

        ! analysis update of members
        DO m=1,MEMBER
          anal3d(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n)                                !GYL
          DO k=1,MEMBER
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &                                !GYL
                                + gues3d(ij,ilev,k,n) * transrlx(k,m)                  !GYL
          END DO  
        END DO

        ! analysis update of deterministic run
        if (DET_RUN) then                                                              !GYL
          anal3d(ij,ilev,mmdet,n) = 0.0d0                                              !GYL
          DO k=1,MEMBER                                                                !GYL
            anal3d(ij,ilev,mmdet,n) = anal3d(ij,ilev,mmdet,n) &                        !GYL
                                    + gues3d(ij,ilev,k,n) * transmd(k,n2nc)            !GYL
          END DO                                                                       !GYL
          anal3d(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n) &                          !GYL
                                  + anal3d(ij,ilev,mmdet,n) * beta                     !GYL
        end if                                                                         !GYL

        ! limit q spread
        IF(Q_SPRD_MAX > 0.0d0 .and. n == iv3d_q) THEN                                  !GYL
          q_mean = SUM(anal3d(ij,ilev,1:MEMBER,n)) / REAL(MEMBER,r_size)               !GYL
          q_sprd = 0.0d0                                                               !GYL
          DO m=1,MEMBER                                                                !GYL
            q_anal(m) = anal3d(ij,ilev,m,n) - q_mean                                   !GYL
            q_sprd = q_sprd + q_anal(m)**2                                             !GYL
          END DO                                                                       !GYL
          q_sprd = SQRT(q_sprd / REAL(MEMBER-1,r_size)) / q_mean                       !GYL
          IF(q_sprd > Q_SPRD_MAX) THEN                                                 !GYL
            DO m=1,MEMBER                                                              !GYL
              anal3d(ij,ilev,m,n) = q_mean + q_anal(m) * Q_SPRD_MAX / q_sprd           !GYL
            END DO                                                                     !GYL
          END IF                                                                       !GYL
        END IF                                                                         !GYL

        if (LOG_LEVEL >= 3) then
          write (timer_str, '(A30,I4,A4,I4,A3,I2,A2)') 'das_letkf:letkf_core_anal(lev=', ilev, ',ij=', ij, ',n=', n, '):'
          call mpi_timer(trim(timer_str), 5)
        end if

      END DO ! [ n=1,nv3d ]

      if (LOG_LEVEL >= 3) then
        write (timer_str, '(A25,I4,A4,I4,A2)') 'das_letkf:letkf_core(lev=', ilev, ',ij=', ij, '):'
        call mpi_timer(trim(timer_str), 4)
      end if

      ! update 2D variables at ilev = 1
      IF(ilev == 1) THEN 

        DO n=1,nv2d

          n2nc = var_local_n2nc(nv3d+n)
          n2n = var_local_n2n(nv3d+n)

          if (RELAX_TO_INFLATED_PRIOR) then
            parm = work2d(ij,n)
          else
            parm = 1.0d0
          end if

          ! calculate mean and perturbation weights
          if (trans_done(n2nc)) then
            ! if weights already computed for other variables can be re-used(no variable localization), do not need to compute again
            IF(n2n <= nv3d) then
              if (INFL_MUL_ADAPTIVE) then
                work2d(ij,n) = work3d(ij,ilev,n2n)
              end if
              if (NOBS_OUT) then
                work2dn(:,ij,n) = work3dn(:,ij,ilev,n2n)
              end if
            else
              if (INFL_MUL_ADAPTIVE) then
                work2d(ij,n) = work2d(ij,n2n-nv3d)
              end if
              if (NOBS_OUT) then
                work2dn(:,ij,n) = work2dn(:,ij,n2n-nv3d)
              end if
            end if

            if (LOG_LEVEL >= 3) then
              write (timer_str, '(A32,I4,A3,I2,A2)') 'das_letkf:letkf_core_copy(2d,ij=', ij, ',n=', n, '):'
              call mpi_timer(trim(timer_str), 5)
            end if
          ELSE
            ! compute weights with localized observations
            if (DET_RUN) then                                                          !GYL
              CALL obs_local(rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),nv3d+n,hdxf,rdiag,rloc,dep,nobsl,depd=depd,nobsl_t=nobsl_t,cutd_t=cutd_t,srch_q0=search_q0(:,nv3d+1,ij))
            else                                                                       !GYL
              CALL obs_local(rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),nv3d+n,hdxf,rdiag,rloc,dep,nobsl,nobsl_t=nobsl_t,cutd_t=cutd_t,srch_q0=search_q0(:,nv3d+1,ij))
            end if                                                                     !GYL
            IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                       !GYL
              if (DET_RUN) then                                                        !GYL
                CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & !GYL
                                trans(:,:,n2nc),transm=transm(:,n2nc),pao=pa(:,:,n2nc), & !GYL
                                rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE, &     !GYL
                                depd=depd,transmd=transmd(:,n2nc))                     !GYL
              else                                                                     !GYL
                CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & !GYL
                                trans(:,:,n2nc),transm=transm(:,n2nc),pao=pa(:,:,n2nc), & !GYL
                                rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)       !GYL
              end if                                                                   !GYL
            ELSE                                                                       !GYL
              if (DET_RUN) then                                                        !GYL
                CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & !GYL
                                trans(:,:,n2nc),transm=transm(:,n2nc),           &     !GYL
                                rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE, &     !GYL
                                depd=depd,transmd=transmd(:,n2nc))                     !GYL
              else                                                                     !GYL
                CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & !GYL
                                trans(:,:,n2nc),transm=transm(:,n2nc),           &     !GYL
                                rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)       !GYL
              end if                                                                   !GYL
            END IF                                                                     !GYL
            trans_done(n2nc) = .true.                                                  !GYL
            IF(NOBS_OUT) THEN                                                          !GYL
              work2dn(:,ij,n) = real(sum(nobsl_t,dim=1),r_size)                        !GYL !!! NOBS: sum over all variables for each report type
            END IF                                                                     !GYL

            if (LOG_LEVEL >= 3) then
              write (timer_str, '(A32,I4,A3,I2,A2)') 'das_letkf:letkf_core_calc(2d,ij=', ij, ',n=', n, '):'
              call mpi_timer(trim(timer_str), 5)
            end if
          END IF

          ! relaxation via LETKF weight
          IF(RELAX_ALPHA /= 0.0d0) THEN                                              !GYL - RTPP method (Zhang et al. 2004)
            CALL weight_RTPP(trans(:,:,n2nc),parm,transrlx)                          !GYL
          ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                  !GYL - RTPS method (Whitaker and Hamill 2012)
            IF(RELAX_SPREAD_OUT) THEN                                                !GYL
              CALL weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues2d(ij,:,n), &        !GYL
                               parm,transrlx,work2da(ij,n))                          !GYL
            ELSE                                                                     !GYL
              CALL weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues2d(ij,:,n), &        !GYL
                               parm,transrlx,tmpinfl)                                !GYL
            END IF                                                                   !GYL
          ELSE                                                                       !GYL
            transrlx = trans(:,:,n2nc)                                               !GYL - No relaxation
          END IF                                                                     !GYL

          ! total weight matrix
          DO m=1,MEMBER                                                              !GYL
            DO k=1,MEMBER                                                            !GYL
              transrlx(k,m) = (transrlx(k,m) + transm(k,n2nc)) * beta                !GYL
            END DO                                                                   !GYL
            transrlx(m,m) = transrlx(m,m) + (1.0d0-beta)                             !GYL
          END DO                                                                     !GYL

          ! analysis update of members
          DO m=1,MEMBER
            anal2d(ij,m,n) = gues2d(ij,mmean,n)                                      !GYL
            DO k=1,MEMBER
              anal2d(ij,m,n) = anal2d(ij,m,n) &                                      !GYL
                             + gues2d(ij,k,n) * transrlx(k,m)                        !GYL
            END DO
          END DO

          ! analysis update of deterministic run
          if (DET_RUN) then                                                          !GYL
            anal2d(ij,mmdet,n) = 0.0d0                                               !GYL
            DO k=1,MEMBER                                                            !GYL
              anal2d(ij,mmdet,n) = anal2d(ij,mmdet,n) &                              !GYL
                                 + gues2d(ij,k,n) * transmd(k,n2nc)                  !GYL
            END DO                                                                   !GYL
            anal2d(ij,mmdet,n) = gues2d(ij,mmdet,n) &                                !GYL
                               + anal2d(ij,mmdet,n) * beta                           !GYL
          end if                                                                     !GYL

          if (LOG_LEVEL >= 3) then
            write (timer_str, '(A32,I4,A3,I2,A2)') 'das_letkf:letkf_core_anal(2d,ij=', ij, ',n=', n, '):'
            call mpi_timer(trim(timer_str), 5)
          end if

        END DO ! [ n=1,nv2d ]

        if (LOG_LEVEL >= 3) then
          write (timer_str, '(A27,I4,A2)') 'das_letkf:letkf_core(2d,ij=', ij, '):'
          call mpi_timer(trim(timer_str), 4)
        end if

      END IF ! [ ilev == 1 ]

    END DO ! [ ij=1,nij1 ]
!$OMP END DO

!$OMP MASTER
    if (LOG_LEVEL >= 3) then
      if (maxval(MAX_NOBS_PER_GRID(:)) > 0 .and. MAX_NOBS_PER_GRID_CRITERION == 1) then
        do ic = 1, nctype
          write (6, '(A,I4,A)') 'ic=', ic, ': search_q0='
          write (6, *) search_q0(ic,1,:)
        end do
        if (ilev == 1) then
          write (6, '(A)') 'For 2d variables:'
          do ic = 1, nctype
            write (6, '(A,I4,A)') 'ic=', ic, ': search_q0='
            write (6, *) search_q0(ic,nv3d+1,:)
          end do
        end if
      end if
    end if

    write (timer_str, '(A25,I4,A2)') 'das_letkf:letkf_core(lev=', ilev, '):'
    call mpi_timer(trim(timer_str), 3)
!$OMP END MASTER

  END DO ! [ ilev=1,nlev ]

  deallocate (hdxf,rdiag,rloc,dep)
  if (DET_RUN) then
    deallocate (depd)
  end if
  deallocate (trans,transm,transmd,pa)
!$OMP END PARALLEL

  call mpi_timer('das_letkf:letkf_core:', 2)

  deallocate (n_merge,ic_merge)
  deallocate (search_q0)
  !
  ! Compute analyses of observations (Y^a)
  !
!!  IF(obsanal_output) THEN
!!    call das_letkf_obs(work3dg,work2dg)
!!  END IF
  !
  ! Write updated inflation parameters
  !
  IF(INFL_MUL_ADAPTIVE) THEN
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    if (.not. allocated(work3dg)) allocate (work3dg(nlon,nlat,nlev,nv3d))
    if (.not. allocated(work2dg)) allocate (work2dg(nlon,nlat,nv2d))
    CALL gather_grd_mpi(mmean_rank_e,work3d,work2d,work3dg,work2dg)

    call mpi_timer('das_letkf:adaptive_infl_gather:', 2)

    IF(myrank_e == mmean_rank_e) THEN
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',INFL_MUL_OUT_BASENAME,'.pe',myrank_d,'.nc'
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call write_restart_par(INFL_MUL_OUT_BASENAME,work3dg,work2dg,MPI_COMM_d)
      else
#endif
        call write_restart(INFL_MUL_OUT_BASENAME,work3dg,work2dg)
#ifdef PNETCDF
      end if
#endif

      call mpi_timer('das_letkf:adaptive_infl_write_restart:', 2)
    END IF
  END IF
  !
  ! Write inflation parameter (in analysis) corresponding to the RTPS method
  !
  IF(RELAX_SPREAD_OUT) THEN
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    if (.not. allocated(work3dg)) allocate (work3dg(nlon,nlat,nlev,nv3d))
    if (.not. allocated(work2dg)) allocate (work2dg(nlon,nlat,nv2d))
    CALL gather_grd_mpi(mmean_rank_e,work3da,work2da,work3dg,work2dg)

    call mpi_timer('das_letkf:relax_spread_out_gather:', 2)

    IF(myrank_e == mmean_rank_e) THEN
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',RELAX_SPREAD_OUT_BASENAME,'.pe',myrank_d,'.nc'
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call write_restart_par(RELAX_SPREAD_OUT_BASENAME,work3dg,work2dg,MPI_COMM_d)
      else
#endif
        call write_restart(RELAX_SPREAD_OUT_BASENAME,work3dg,work2dg)
#ifdef PNETCDF
      end if
#endif

      call mpi_timer('das_letkf:relax_spread_out_write_restart:', 2)
    END IF
    DEALLOCATE(work3da,work2da)
  END IF
  !
  ! Write observation numbers
  !
  IF(NOBS_OUT) THEN
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    if (.not. allocated(work3dg)) allocate (work3dg(nlon,nlat,nlev,nv3d))
    if (.not. allocated(work2dg)) allocate (work2dg(nlon,nlat,nv2d))
    work3d(:,:,1) = work3dn(1,:,:,iv3d_t)  !!! Assuming variable localization is not used so that obs numbers used are the same over variables,
    work3d(:,:,2) = work3dn(3,:,:,iv3d_t)  !!! use "variable dimenstion" to save obs numbers of different observation types
    work3d(:,:,3) = work3dn(4,:,:,iv3d_t)
    work3d(:,:,4) = work3dn(8,:,:,iv3d_t)
    work3d(:,:,5) = work3dn(22,:,:,iv3d_t)
    work3d(:,:,6) = work3dn(nobtype+1,:,:,iv3d_t)
    work3d(:,:,7) = work3dn(nobtype+2,:,:,iv3d_t)
    work3d(:,:,8) = work3dn(nobtype+3,:,:,iv3d_t)
    work3d(:,:,9) = work3dn(nobtype+4,:,:,iv3d_t)
    work3d(:,:,10) = work3dn(nobtype+5,:,:,iv3d_t)
    work3d(:,:,11) = work3dn(nobtype+6,:,:,iv3d_t)
    CALL gather_grd_mpi(mmean_rank_e,work3d,work2d,work3dg,work2dg)

    call mpi_timer('das_letkf:nobs_out_gather:', 2)

    IF(myrank_e == mmean_rank_e) THEN
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',NOBS_OUT_BASENAME,'.pe',myrank_d,'.nc'
#ifdef PNETCDF
      if (FILE_AGGREGATE) then
        call write_restart_par(NOBS_OUT_BASENAME,work3dg,work2dg,MPI_COMM_d)
      else
#endif
        call write_restart(NOBS_OUT_BASENAME,work3dg,work2dg)
#ifdef PNETCDF
      end if
#endif

      call mpi_timer('das_letkf:nobs_out_write_restart:', 2)
    END IF
    DEALLOCATE(work3dn,work2dn)
  END IF
  IF (allocated(work3dg)) deallocate (work3dg)
  IF (allocated(work2dg)) deallocate (work2dg)
  !
  ! Additive inflation
  !
  IF(INFL_ADD > 0.0d0) THEN
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    if (INFL_ADD_Q_RATIO) then
      work3d(:,:,:) = gues3d(:,:,mmean,:)
    else
      work3d(:,:,:) = 1.0d0
    end if

    allocate (addinfl_weight(nij1))
    if (INFL_ADD_REF_ONLY) then
      addinfl_weight(:) = 0.0d0
      ic = ctype_elmtyp(uid_obs(id_radar_ref_obs), 22)
      if (ic > 0) then
        do ij = 1, nij1
          ref_min_dist = 1.0d33
          !!!!!! save this (ref_min_dist) information when doing DA
          do iob = obsgrd(ic)%ac_ext(0, 1) + 1, obsgrd(ic)%ac_ext(obsgrd(ic)%ngrdext_i, obsgrd(ic)%ngrdext_j)
            rdx = (rig1(ij) - obs(obsda_sort%set(iob))%ri(obsda_sort%idx(iob))) * DX
            rdy = (rjg1(ij) - obs(obsda_sort%set(iob))%rj(obsda_sort%idx(iob))) * DY
            rdxy = rdx*rdx + rdy*rdy
            if (rdxy < ref_min_dist) then
              ref_min_dist = rdxy
            end if
          end do

          ref_min_dist = ref_min_dist / (hori_loc_ctype(ic) * hori_loc_ctype(ic))
          if (ref_min_dist <= dist_zero_fac_square) then
            addinfl_weight(ij) = EXP(-0.5d0 * ref_min_dist)
          end if
        end do
      end if
    else
      addinfl_weight(:) = 1.0d0
    end if

    call mpi_timer('das_letkf:additive_infl_addinfl_weight:', 2)

    CALL read_ens_mpi_addiinfl(gues3d,gues2d)

    call mpi_timer('das_letkf:additive_infl_read_ens_mpi:', 2)

    CALL ensmean_grd(MEMBER,nens,nij1,gues3d,gues2d)

    call mpi_timer('das_letkf:additive_infl_ensmean_grd:', 2)

    write (6, '(A)') '===== Additive covariance inflation ====='
    write (6, '(A,F10.4)') '  parameter:', INFL_ADD
    if (INFL_ADD_SHUFFLE) then
      if (myrank_a == 0) then
        call Knuth_Shuffle(MEMBER, ishuf)
      end if
      call MPI_BCAST(ishuf, MEMBER, MPI_INTEGER, 0, MPI_COMM_a, ierr)
      write (6, '(A)') '  suffle members: on'
      write (6, *) ' suffle sequence: ', ishuf
    end if
    write (6,'(A)') '========================================='

!$OMP PARALLEL PRIVATE(n,m,k,i,mshuf)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2) 
    DO n=1,nv3d
      DO m=1,MEMBER
        DO k=1,nlev
          DO i=1,nij1
            gues3d(i,k,m,n) = gues3d(i,k,m,n) - gues3d(i,k,mmean,n)
          END DO
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    DO n=1,nv2d
      DO m=1,MEMBER
        DO i=1,nij1
          gues2d(i,m,n) = gues2d(i,m,n) - gues2d(i,mmean,n)
        END DO
      END DO
    END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC) COLLAPSE(2) 
    DO n=1,nv3d
      DO m=1,MEMBER
        if (INFL_ADD_SHUFFLE) then
          mshuf = ishuf(m)
        else
          mshuf = m
        end if
        if (n == iv3d_q .or. n == iv3d_qc .or. n == iv3d_qr .or. n == iv3d_qi .or. n == iv3d_qs .or. n == iv3d_qg) then
          DO k=1,nlev
            DO i=1,nij1
              anal3d(i,k,m,n) = anal3d(i,k,m,n) &
                & + gues3d(i,k,mshuf,n) * INFL_ADD * addinfl_weight(i) * work3d(i,k,n)
            END DO
          END DO
        else
          DO k=1,nlev
            DO i=1,nij1
              anal3d(i,k,m,n) = anal3d(i,k,m,n) &
                & + gues3d(i,k,mshuf,n) * INFL_ADD * addinfl_weight(i)
            END DO
          END DO
        end if
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    DO n=1,nv2d
      DO m=1,MEMBER
        if (INFL_ADD_SHUFFLE) then
          mshuf = ishuf(m)
        else
          mshuf = m
        end if
        DO i=1,nij1
          anal2d(i,m,n) = anal2d(i,m,n) + gues2d(i,mshuf,n) * INFL_ADD * addinfl_weight(i)
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    deallocate (addinfl_weight)

    call mpi_timer('das_letkf:additive_infl_cal:', 2)
  END IF ! [ INFL_ADD > 0.0d0 ]

  RETURN
END SUBROUTINE das_letkf
!!-----------------------------------------------------------------------
!! Data assimilation for observations: Compute analyses of observations (Y^a)
!! * currently only support multiplicative and adaptive inflation
!!  -- 01/01/2014, Guo-Yuan Lien, 
!!-----------------------------------------------------------------------
!SUBROUTINE das_letkf_obs(v3dinfl,v2dinfl)
!  IMPLICIT NONE
!  REAL(r_sngl),INTENT(IN) :: v3dinfl(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(IN) :: v2dinfl(nlon,nlat,nv2d)
!  REAL(r_size),ALLOCATABLE :: v3dinflx(:,:,:,:)
!  REAL(r_size),ALLOCATABLE :: v2dinflx(:,:,:)
!  REAL(r_size),ALLOCATABLE :: v3dtmp(:,:,:,:)
!  REAL(r_size),ALLOCATABLE :: v2dtmp(:,:,:)
!  REAL(r_size),ALLOCATABLE :: tmpps(:)
!  REAL(r_size),ALLOCATABLE :: tmptv(:,:)
!  REAL(r_size),ALLOCATABLE :: tmpp(:,:)
!  REAL(r_size),ALLOCATABLE :: obsanal(:,:)
!  REAL(r_size),ALLOCATABLE :: obsanalmean(:)
!  REAL(r_size) :: hdxf(nobstotal,MEMBER)
!  REAL(r_size) :: rdiag(nobstotal)
!  REAL(r_size) :: rloc(nobstotal)
!  REAL(r_size) :: dep(nobstotal)
!  REAL(r_size) :: ohx(nobs)
!  REAL(r_size) :: parm
!  REAL(r_size) :: trans(MEMBER,MEMBER)
!  REAL(r_size) :: ri,rj,rk
!  REAL(r_size) :: rlev,p_update_q
!  REAL(r_size) :: q_sprd
!  REAL(r_size) :: q_anal(MEMBER)
!  INTEGER :: n,nn,m,k,nobsl,ierr,iret
!  INTEGER :: inflelem,irank,nobsp,nobspmax
!  CHARACTER(14) :: obsanalfile='obsanalNNN.dat'

!  WRITE(6,'(A)') 'Hello from das_letkf_obs: Compute [Y^a]'
!  !
!  ! If adaptive inflation is used, prepare a global array of inflation parameter
!  !
!  IF(COV_INFL_MUL <= 0.0d0) THEN
!    ALLOCATE(v3dinflx(nlon,nlat,nlev,nv3dx))
!    ALLOCATE(v2dinflx(nlon,nlat,nv2dx))
!    IF(myrank == 0) THEN
!      ALLOCATE(v3dtmp(nlon,nlat,nlev,nv3d))
!      ALLOCATE(v2dtmp(nlon,nlat,nv2d))
!      ALLOCATE(tmpps(nlon*nlat))
!      ALLOCATE(tmptv(nlon*nlat,nlev))
!      ALLOCATE(tmpp(nlon*nlat,nlev))
!      CALL read_grd('gues_me.grd',v3dtmp,v2dtmp,0)  ! read ensemble mean into a temporary array
!      CALL read_grdx('gues001.grd',v3dinflx,v2dinflx) ! only the orography is used, P will be recalulated
!      v3dinflx(:,:,:,iv3d_u) = v3dinfl(:,:,:,iv3d_u)
!      v3dinflx(:,:,:,iv3d_v) = v3dinfl(:,:,:,iv3d_v)
!      v3dinflx(:,:,:,iv3d_t) = v3dinfl(:,:,:,iv3d_t)
!      v3dinflx(:,:,:,iv3d_q) = v3dinfl(:,:,:,iv3d_q)
!      v3dinflx(:,:,:,iv3d_qc) = v3dinfl(:,:,:,iv3d_qc)
!!      v2dinflx(:,:,iv2d_ps) = v2dinfl(:,:,iv2d_ps)
!      v2dinflx(:,:,iv2d_ps) = v3dinfl(:,:,1,iv3d_u)
!      tmpps = reshape(v2dtmp(:,:,iv2d_ps),(/nlon*nlat/))
!      tmptv = reshape(v3dtmp(:,:,:,iv3d_t) * (1.0d0 + fvirt * v3dtmp(:,:,:,iv3d_q)),(/nlon*nlat,nlev/))
!      call sigio_modprd(nlon*nlat,nlon*nlat,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl, &
!                        gfs_vcoord,iret,tmpps,tmptv,pm=tmpp)
!      v3dinflx(:,:,:,iv3d_p) = reshape(tmpp,(/nlon,nlat,nlev/))
!      DEALLOCATE(v3dtmp,v2dtmp,tmpps,tmptv,tmpp)
!    END IF
!    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    call MPI_BCAST(v3dinflx,nlon*nlat*nlev*nv3dx,MPI_r_size,0,MPI_COMM_WORLD,ierr)
!    call MPI_BCAST(v2dinflx,nlon*nlat*nv2dx,MPI_r_size,0,MPI_COMM_WORLD,ierr)
!  END IF
!  !
!  ! Define the partition of observations for parallel computation
!  !
!  nn = MOD(nobs,nprocs)
!  nobspmax = (nobs - nn)/nprocs + 1
!  IF(myrank < nn) THEN
!    nobsp = nobspmax
!  ELSE
!    nobsp = nobspmax-1
!  END IF
!  WRITE(6,'(A,I3.3,A,I8)') 'MYRANK ',myrank,' process obs number=', nobsp
!  !
!  ! Main LETKF loop
!  !
!  ALLOCATE(obsanal(nobs,MEMBER))
!  ALLOCATE(obsanalmean(nobs))
!  obsanal = 0.0d0
!  obsanalmean = 0.0d0
!  nn = myrank+1
!  DO
!    IF(nn > nobs) EXIT
!!    WRITE(6,'(A,I8)') 'nn = ',nn
!    !
!    ! The observation variable type is different from the grid variable type.
!    ! To compute the analyses of observations as regular grids,
!    ! what grid variable will the observation variable be regarded as?
!    !
!    ! Also determine the pressure level will the observation variable be regarded as?
!    !
!    SELECT CASE(NINT(obselm(nn)))
!    CASE(id_u_obs)
!      n = iv3d_u          ! for variable localization, what grid variable to be regarded as? 
!      inflelem = id_u_obs ! for inflation parameter,   what grid variable to be regarded as?
!      rlev = obslev(nn)
!    CASE(id_v_obs)
!      n = iv3d_v
!      inflelem = id_v_obs
!      rlev = obslev(nn)
!    CASE(id_t_obs,id_tv_obs)
!      n = iv3d_t
!      inflelem = id_t_obs
!      rlev = obslev(nn)
!    CASE(id_q_obs,id_rh_obs)
!      n = iv3d_q
!      inflelem = id_q_obs
!      rlev = obslev(nn)
!    CASE(id_ps_obs)
!      n = nv3d+iv2d_ps
!      inflelem = id_ps_obs
!      rlev = obsdat(nn)   ! for ps variable, use the observed pressure value
!    CASE(id_rain_obs)
!      n = 0
!      inflelem = id_q_obs
!      rlev = base_obsv_rain ! for precipitation, assigh the level 'base_obsv_rain'
!    CASE DEFAULT
!      n = 0
!      IF(NINT(obselm(nn)) > 9999) THEN
!        inflelem = id_ps_obs
!        CALL itpl_2d(v3dinflx(:,:,1,iv3d_p),ri,rj,rlev)
!      ELSE
!        inflelem = id_u_obs
!        rlev = obslev(nn)
!      END IF
!    END SELECT
!    !
!    ! Determine the inflation parameter
!    !
!    IF(COV_INFL_MUL > 0.0d0) THEN
!      parm = COV_INFL_MUL
!    ELSE
!      CALL phys2ijk(v3dinflx(:,:,:,iv3d_p),real(inflelem,r_size),obslon(nn),obslat(nn),rlev,ri,rj,rk)
!      IF(CEILING(rk) > nlev) THEN
!        rk = REAL(nlev,r_size)
!      END IF
!      IF(CEILING(rk) < 2 .AND. inflelem /= id_ps_obs) THEN
!        IF(inflelem > 9999) THEN
!          rk = 0.0d0
!        ELSE
!          rk = 1.00001d0
!        END IF
!      END IF
!      IF(inflelem == id_ps_obs) THEN
!        CALL itpl_2d(v2dinflx(:,:,iv2d_orog),ri,rj,rk)
!        rk = obslev(nn) - rk
!      END IF
!      CALL Trans_XtoY(real(inflelem,r_size),ri,rj,rk,v3dinflx,v2dinflx,parm)
!    END IF
!    !
!    ! LETKF computation
!    !
!    CALL obs_local(obslon(nn),obslat(nn),rlev,n,hdxf,rdiag,rloc,dep,nobsl)
!    CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans,RELAX_ALPHA)

!    IF(n == iv3d_q .OR. n == iv3d_qc) THEN
!      CALL itpl_2d(v3dinflx(:,:,LEV_UPDATE_Q,iv3d_p),ri,rj,p_update_q)
!    END IF
!    IF((n == iv3d_q .OR. n == iv3d_qc) .AND. obslev(nn) < p_update_q) THEN
!      obsanal(nn,:) = obsdat(nn) - obsdep(nn) + obshdxf(nn,:)
!      obsanalmean(nn) = obsdat(nn) - obsdep(nn)
!    ELSE
!      DO m=1,MEMBER
!        obsanal(nn,m) = obsdat(nn) - obsdep(nn)
!        DO k=1,MEMBER
!          obsanal(nn,m) = obsanal(nn,m) + obshdxf(nn,k) * trans(k,m)
!        END DO
!        obsanalmean(nn) = obsanalmean(nn) + obsanal(nn,m)
!      END DO
!      obsanalmean(nn) = obsanalmean(nn) / real(MEMBER,r_size)
!    END IF
!    IF(n == iv3d_q .AND. obslev(nn) >= p_update_q) THEN
!      q_sprd = 0.0d0
!      DO m=1,MEMBER
!        q_anal(m) = obsanal(nn,m) - obsanalmean(nn)
!        q_sprd = q_sprd + q_anal(m)**2
!      END DO
!      q_sprd = SQRT(q_sprd / REAL(MEMBER-1,r_size)) / obsanalmean(nn)
!      IF(q_sprd > Q_SPRD_MAX) THEN
!        DO m=1,MEMBER
!          obsanal(nn,m) = obsanalmean(nn) + q_anal(m) * Q_SPRD_MAX / q_sprd
!        END DO
!      END IF
!    END IF

!    nn = nn + nprocs
!  END DO
!  !
!  ! MPI_REDUCE and output obsanalfiles
!  !
!  ! mean
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  CALL MPI_REDUCE(obsanalmean,ohx,nobs,MPI_r_size,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!  IF(myrank == 0) THEN
!    WRITE(obsanalfile(8:10),'(A3)') '_me'
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsanalfile
!    CALL write_obs2(obsanalfile,nobs,obselm,obslon,obslat,obslev, &
!                    obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!  END IF
!  ! members
!  irank = 0
!  DO m=1,MEMBER
!    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    CALL MPI_REDUCE(obsanal(:,m),ohx,nobs,MPI_r_size,MPI_SUM,irank,MPI_COMM_WORLD,ierr)
!    IF(myrank == irank) THEN
!      WRITE(obsanalfile(8:10),'(I3.3)') m
!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsanalfile
!      CALL write_obs2(obsanalfile,nobs,obselm,obslon,obslat,obslev, &
!                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!    END IF
!    irank = irank + 1
!    IF(irank >= nprocs) irank = 0
!  END DO

!  DEALLOCATE(obsanal)
!  IF(COV_INFL_MUL <= 0.0d0) THEN
!    DEALLOCATE(v3dinflx,v2dinflx)
!  END IF
!  RETURN
!END SUBROUTINE das_letkf_obs
!!-----------------------------------------------------------------------
!! Subroutine for observation sensitivity computation
!! Ported from Y.Ohta's SPEEDY-LETKF system by D.Hotta, 07/01/2013
!! [ref: Eq.(6,7,9), Ota et al. 2013]
!!-----------------------------------------------------------------------
!! [INPUT]
!!  gues3d,gues2d: xmean^g_0
!!  fcst3d,fcst2d: C^(1/2)*X^f_t                    [(J/kg)^(1/2)]
!!  fcer3d,fcer2d: C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)  [(J/kg)^(1/2)]
!! (save variables)
!!  obshdxf:
!! [OUTPUT]
!!-----------------------------------------------------------------------
!SUBROUTINE das_efso(gues3d,gues2d,fcst3d,fcst2d,fcer3d,fcer2d)
!  IMPLICIT NONE
!  REAL(r_size),INTENT(IN) :: gues3d(nij1,nlev,nv3d)     !
!  REAL(r_size),INTENT(IN) :: gues2d(nij1,nv2d)          !
!  REAL(r_size),INTENT(IN) :: fcst3d(nij1,nlev,MEMBER,nv3d) ! forecast ensemble
!  REAL(r_size),INTENT(IN) :: fcst2d(nij1,MEMBER,nv2d)      !
!  REAL(r_size),INTENT(IN) :: fcer3d(nij1,nlev,nv3d) ! forecast error
!  REAL(r_size),INTENT(IN) :: fcer2d(nij1,nv2d)      !
!  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
!  REAL(r_size),ALLOCATABLE :: hdxa_rinv(:,:)
!  REAL(r_size),ALLOCATABLE :: rdiag(:)
!  REAL(r_size),ALLOCATABLE :: rloc(:)
!  REAL(r_size),ALLOCATABLE :: dep(:)
!  REAL(r_size),ALLOCATABLE :: tmptv(:,:)
!  REAL(r_size),ALLOCATABLE :: pfull(:,:)
!  REAL(r_size),ALLOCATABLE :: djdy(:,:)
!  REAL(r_size),ALLOCATABLE :: recbuf(:,:)
!  REAL(r_size) :: work1(nterm,MEMBER)
!  INTEGER,ALLOCATABLE :: oindex(:)
!  INTEGER :: ij,k,ilev,m,nob,nobsl,ierr,iret,iterm

!  WRITE(6,'(A)') 'Hello from das_obsense'
!  nobstotal = nobs !+ ntvs
!  WRITE(6,'(A,I8)') 'Target observation numbers : NOBS=',nobs!,', NTVS=',ntvs
!  !
!  ! In case of no obs
!  !
!  IF(nobstotal == 0) THEN
!    WRITE(6,'(A)') 'No observation assimilated'
!    RETURN
!  END IF
!  ALLOCATE(djdy(nterm,nobstotal))
!  djdy = 0.0_r_size
!  !
!  ! p_full for background ensemble mean
!  !
!  ALLOCATE( tmptv(nij1,nlev) )
!  ALLOCATE( pfull(nij1,nlev) )
!  tmptv = gues3d(:,:,iv3d_t) * (1.0d0 + fvirt * gues3d(:,:,iv3d_q))
!  call sigio_modprd(nij1,nij1,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl, &
!                    gfs_vcoord,iret,gues2d(:,iv2d_ps),tmptv,pm=pfull)
!  DEALLOCATE(tmptv)
!  !
!  ! MAIN ASSIMILATION LOOP
!  !
!!$OMP PARALLEL PRIVATE(ij,ilev,k,hdxf,rdiag,rloc,dep,nobsl,oindex, &
!!$                     work1,m,nob)
!  ALLOCATE( hdxf(1:nobstotal,1:MEMBER),rdiag(1:nobstotal),rloc(1:nobstotal), &
!       & dep(1:nobstotal) )
!  ALLOCATE(oindex(1:nobstotal))
!!--- For ILEV = 1 - NLEV
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO ilev=1,nlev
!    WRITE(6,'(A,I3)') 'ilev = ',ilev
!    DO ij=1,nij1
!      IF(ABS(locadv_rate) > TINY(locadv_rate)) THEN
!        CALL obs_local(lon2(ij,ilev),lat2(ij,ilev),pfull(ij,ilev),0,hdxf,rdiag,rloc,dep,nobsl,oindex)
!      ELSE
!        CALL obs_local(lon1(ij),lat1(ij),pfull(ij,ilev),0,hdxf,rdiag,rloc,dep,nobsl,oindex)
!      END IF
!      IF( nobsl /= 0 ) THEN
!        ! Forecast error
!        work1 = 0.0_r_size
!        DO k=1,nv3d
!          SELECT CASE(k)
!          CASE(iv3d_u,iv3d_v)
!            iterm = 1
!          CASE(iv3d_t)
!            iterm = 2
!          CASE(iv3d_q)
!            iterm = 3
!          CASE DEFAULT
!            iterm = 0
!          END SELECT
!          IF(iterm > 0) THEN
!            DO m=1,MEMBER
!              work1(iterm,m) = work1(iterm,m) + fcst3d(ij,ilev,m,k) * fcer3d(ij,ilev,k)
!            END DO
!          END IF
!        END DO
!        IF(ilev == 1) THEN
!          DO k=1,nv2d
!            IF(k == iv2d_ps) THEN
!              DO m=1,MEMBER
!                work1(2,m) = work1(2,m) + fcst2d(ij,m,k) * fcer2d(ij,k)
!              END DO
!            END IF
!          END DO
!        END IF
!        !!! work1: [1/2(K-1)](X^f_t)^T*C*(e^f_t+e^g_t)  [J/kg]
!        ! Hdxa Rinv
!        ALLOCATE(hdxa_rinv(nobsl,MEMBER))
!        DO m=1,MEMBER
!          DO nob=1,nobsl
!            hdxa_rinv(nob,m) = hdxf(nob,m) / rdiag(nob) * rloc(nob)
!          END DO
!        END DO
!        !!! hdxa_rinv: rho*R^(-1)*Y^a_0 = rho*R^(-1)*(H X^a_0)
!        ! dJ/dy
!        DO nob=1,nobsl
!          DO m=1,MEMBER
!            djdy(:,oindex(nob)) = djdy(:,oindex(nob)) + work1(:,m) * hdxa_rinv(nob,m)
!          END DO
!        END DO
!        !!! djdy: [1/2(K-1)]rho*R^(-1)*Y^a_0*(X^f_t)^T*C*(e^f_t+e^g_t)
!        DEALLOCATE(hdxa_rinv)
!      END IF
!    END DO
!  END DO
!!$OMP END DO
!  DEALLOCATE(hdxf,rdiag,rloc,dep,oindex)
!!$OMP END PARALLEL
!  !
!  ! Calculate observation sensitivity
!  !
!!$OMP PARALLEL PRIVATE(nob)
!!$OMP DO
!  DO nob=1,nobstotal
!    obsense(:,nob) = djdy(:,nob) * obsdep(nob)
!  END DO
!  !!! obsense: delta e^{f-g}_t = [1/2(K-1)][y_o-H(xmean^b_0)]^T*rho*R^(-1)*Y^a_0*(X^f_t)^T*C*(e^f_t+e^g_t)
!!$OMP END DO
!!$OMP END PARALLEL
!  ! Gather observation sensitivity informations to the root
!  ALLOCATE(recbuf(nterm,nobstotal))
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  CALL MPI_REDUCE(obsense(:,1:nobstotal),recbuf,nterm*nobstotal,MPI_r_size,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!  IF(myrank == 0) obsense(:,1:nobstotal) = recbuf(:,:)
!  DEALLOCATE(recbuf)
!  DEALLOCATE(djdy)
!  DEALLOCATE(pfull)
!  RETURN
!END SUBROUTINE das_efso

!-------------------------------------------------------------------------------
! Find local observations to be used for a targeted grid
!-------------------------------------------------------------------------------
! [INPUT]
!   ri      : horizontal i-grid cooridnate of the targeted grid
!   rj      : horizontal j-grid cooridnate of the targeted grid
!   rlev    : vertical pressure of the targeted grid
!   rz      : vertical height   of the targeted grid
!   nvar    : variable index of the targeted grid
!   srch_q0 : (optional) first guess of the multiplier of incremental search distances
! [OUT]
!   hdxf    : fcstast ensemble perturbations in the observation space
!   rdiag   : localization-weighted observation error variances
!   rloc    : localization weights
!   dep     : observation departure
!   nobsl   : number of valid observations (in hdxf, rdiag, rloc, dep)
!   depd    : (optional) observation departure for the deterministic run
!   nobsl_t : (optional) number of assimilated observations wrt. observation variables/types
!   cutd_t  : (optional) cutoff distance of assimilated observations wrt. observation variables/types
!   srch_q0 : (optional) revised first guess of the multiplier of incremental search distances for the next call
!-------------------------------------------------------------------------------
subroutine obs_local(ri, rj, rlev, rz, nvar, hdxf, rdiag, rloc, dep, nobsl, depd, nobsl_t, cutd_t, srch_q0)
  use common_sort
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  implicit none

  real(r_size), intent(in) :: ri, rj, rlev, rz
  integer, intent(in) :: nvar
  real(r_size), intent(out) :: hdxf(nobstotal,MEMBER)
  real(r_size), intent(out) :: rdiag(nobstotal)
  real(r_size), intent(out) :: rloc(nobstotal)
  real(r_size), intent(out) :: dep(nobstotal)
  integer, intent(out) :: nobsl
  real(r_size), intent(out), optional :: depd(nobstotal)
  integer, intent(out), optional :: nobsl_t(nid_obs,nobtype)
  real(r_size), intent(out), optional :: cutd_t(nid_obs,nobtype)
  integer, intent(inout), optional :: srch_q0(nctype)

  integer, allocatable :: nobs_use(:)
  integer, allocatable :: nobs_use2(:)
  real(r_size), allocatable :: dist_tmp(:)
  real(r_size), allocatable :: rloc_tmp(:)
  real(r_size), allocatable :: rdiag_tmp(:)

  real(r_size) :: nrloc, nrdiag
  real(r_size) :: ndist_dummy
  integer :: iob, ityp, ielm
  integer :: imin, imax, jmin, jmax
  integer :: ic, ic2, icm
  integer :: n, nn, nn_prev
  integer :: nobsl_prev, nobsl_incr
  integer :: nobsl_max_master
  integer :: ielm_u_master, ityp_master

  integer :: q
  logical :: loop
  integer :: nn_steps(n_merge_max+1)
  logical :: reach_cutoff
  integer :: imin_cutoff(n_merge_max), imax_cutoff(n_merge_max)
  integer :: jmin_cutoff(n_merge_max), jmax_cutoff(n_merge_max)
  real(r_size) :: search_incr(n_merge_max)
  real(r_size) :: search_incr_i(n_merge_max), search_incr_j(n_merge_max)
  real(r_size) :: dist_cutoff_fac, dist_cutoff_fac_square

  real(r_size) :: cutd
  integer :: nobsl_cm(n_merge_max)

  !-----------------------------------------------------------------------------
  ! Initialize
  !-----------------------------------------------------------------------------

  nobsl = 0
  if (present(nobsl_t)) then
    nobsl_t(:,:) = 0
  end if
  if (present(cutd_t)) then
    cutd_t(:,:) = 0.0d0
    if (MAX_NOBS_PER_GRID_CRITERION == 1) then
      do ic = 1, nctype
        cutd_t(elm_u_ctype(ic),typ_ctype(ic)) = hori_loc_ctype(ic) * dist_zero_fac
      end do
    end if
  end if

#ifdef DEBUG
  if (present(depd)) then
    if (.not. DET_RUN) then
      write (6, '(A)') "[Error] If 'depd' optional input is given, 'DET_RUN' needs to be enabled."
      stop 99
    end if
  end if
#endif

  if (nobstotal == 0) then
    return
  end if

  if (maxval(MAX_NOBS_PER_GRID(:)) > 0) then
    allocate (nobs_use (nobstotal))
    allocate (nobs_use2(nobstotal))
    allocate (rloc_tmp (nobstotal))
    allocate (rdiag_tmp(nobstotal))
    if (MAX_NOBS_PER_GRID_CRITERION == 1) then
      allocate (dist_tmp (nobstotal))
    end if
!    dist_tmp(:) = -1.0d6
    rloc_tmp(:) = -1.0d6
!    rdiag_tmp(:) = -1.0d6
  else
    allocate (nobs_use(maxnobs_per_ctype))
  end if

  !-----------------------------------------------------------------------------
  ! For each observation type,
  ! do rough data search by a rectangle using the sorting mesh, and then
  ! do precise data search by normalized 3D distance and variable localization.
  !-----------------------------------------------------------------------------

  do ic = 1, nctype
    if (n_merge(ic) == 0) then
      if (present(cutd_t)) then
        cutd_t(elm_u_ctype(ic),typ_ctype(ic)) = 0.0d0
      end if
      cycle
    end if

    nobsl_max_master = MAX_NOBS_PER_GRID(typ_ctype(ic)) ! Use the number limit setting of the "master" obs type for all group of obs types
    ielm_u_master = elm_u_ctype(ic)                     ! Count observation numbers    at the "master" obs type for all group of obs types
    ityp_master = typ_ctype(ic)                         !

    if (nobsl_max_master <= 0) then
    !---------------------------------------------------------------------------
    ! When obs number limit is not enabled,
    ! directly prepare (hdxf, dep, depd, rdiag, rloc) output.
    !---------------------------------------------------------------------------

      nobsl_prev = nobsl

      do icm = 1, n_merge(ic)
        ic2 = ic_merge(icm,ic)
        ielm = elm_ctype(ic2)
        ityp = typ_ctype(ic2)

        if (obsgrd(ic2)%tot_ext > 0) then
          nn = 0
          call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
          call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)

          do n = 1, nn
            iob = nobs_use(n)

            call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ic2, ndist_dummy, nrloc, nrdiag)
            if (nrloc == 0.0d0) cycle

            nobsl = nobsl + 1
            hdxf(nobsl,:) = obsda_sort%ensval(1:MEMBER,iob)
            rdiag(nobsl) = nrdiag
            rloc(nobsl) = nrloc
            dep(nobsl) = obsda_sort%val(iob)
            if (present(depd)) then
              depd(nobsl) = obsda_sort%ensval(mmdetobs,iob)
            end if
          end do ! [ n = 1, nn ]
        end if ! [ obsgrd(ic2)%tot_ext > 0 ]

        if (present(nobsl_t)) then
          nobsl_t(elm_u_ctype(ic2),ityp) = nobsl - nobsl_prev
        end if
      end do ! [ do icm = 1, n_merge(ic) ]

    !---------------------------------------------------------------------------
    else if (MAX_NOBS_PER_GRID_CRITERION == 1) then
    !---------------------------------------------------------------------------
    ! When obs number limit is enabled and the sorting criterion is simply distance,
    ! try the incremental observation location search
    ! (within the localization cutoff area) before conduct selection.
    !---------------------------------------------------------------------------

      nn = 0
      do icm = 1, n_merge(ic)
        ic2 = ic_merge(icm,ic)

        if (obsgrd(ic2)%tot_ext > 0) then
          call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
          call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)
        end if ! [ obsgrd(ic2)%tot_ext > 0 ]
      end do ! [ do icm = 1, n_merge(ic) ]

      if (LOG_LEVEL >= 3) then
        write (6, '(A,14x,I8)') '--- ALL      : ', nn
      end if

      if (nn == 0) cycle

      ! Determine the search_incr based on the master obs ctype:
      ! (zero-weight distance) / (# of search increment), but not smaller than the sorting mesh size
      search_incr(1) = hori_loc_ctype(ic) * dist_zero_fac / real(n_search_incr, r_size)  ! (unit: m)
      search_incr(1) = max(search_incr(1), obsgrd(ic)%grdspc_i, obsgrd(ic)%grdspc_j)     ! (unit: m)

      do icm = 1, n_merge(ic)
        ic2 = ic_merge(icm,ic)
        ! Determine the search_incr of the other obs ctypes based on its horizontal localization scale
        ! relative to that of the master obs ctype
        if (icm > 1) then
          search_incr(icm) = search_incr(1) / hori_loc_ctype(ic) * hori_loc_ctype(ic2)
        end if
        search_incr_i(icm) = search_incr(icm) / DX  ! (unit: x-grid)
        search_incr_j(icm) = search_incr(icm) / DY  ! (unit: y-grid)
        call obs_local_range(ic2, ri, rj, imin_cutoff(icm), imax_cutoff(icm), jmin_cutoff(icm), jmax_cutoff(icm))
      end do

      nobsl_incr = 0
      if (present(srch_q0)) then
        q = srch_q0(ic) - 1
      else
        q = 0
      end if
      loop = .true.

      do while (loop)
        q = q + 1
        nn = 0
        reach_cutoff = .true.

        do icm = 1, n_merge(ic)
          ic2 = ic_merge(icm,ic)
          nn_steps(icm) = nn

          if (obsgrd(ic2)%tot_ext > 0) then
            call ij_obsgrd_ext(ic2, ri-search_incr_i(icm)*q, rj-search_incr_j(icm)*q, imin, jmin)
            call ij_obsgrd_ext(ic2, ri+search_incr_i(icm)*q, rj+search_incr_j(icm)*q, imax, jmax)

            if (imin <= imin_cutoff(icm) .and. imax >= imax_cutoff(icm) .and. &
                jmin <= jmin_cutoff(icm) .and. jmax >= jmax_cutoff(icm)) then
              imin = imin_cutoff(icm)
              imax = imax_cutoff(icm)
              jmin = jmin_cutoff(icm)
              jmax = jmax_cutoff(icm)
            else
              reach_cutoff = .false.
            end if

            call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)
          end if ! [ obsgrd(ic2)%tot_ext > 0 ]
        end do ! [ do icm = 1, n_merge(ic) ]

        nn_steps(n_merge(ic)+1) = nn

        if (LOG_LEVEL >= 3) then
          write (6, '(A,I4,A,F12.3,L2,I8,8x,10F12.3)') '--- TRY #', q, ': ', search_incr(1)*q, reach_cutoff, nn, search_incr(2:n_merge(ic))*q
        end if

        if ((.not. reach_cutoff) .and. nn < nobsl_max_master) cycle
        if (reach_cutoff) then
          loop = .false.
          if (nn == 0) exit
        end if

        nobsl_incr = 0

        ! Determine the cutoff fraction in this incremental search,
        ! which should be the same for all obs ctypes based on its definition
        dist_cutoff_fac = search_incr(1) * q / hori_loc_ctype(ic)
        dist_cutoff_fac_square = dist_cutoff_fac * dist_cutoff_fac

        do icm = 1, n_merge(ic)
          ic2 = ic_merge(icm,ic)

          if (nn_steps(icm+1) > nn_steps(icm)) then
            do n = nn_steps(icm)+1, nn_steps(icm+1)
              iob = nobs_use(n)

              if (rloc_tmp(iob) == 0.0d0) cycle
                
              if (rloc_tmp(iob) < 0.0d0) then
                call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ic2, dist_tmp(iob), rloc_tmp(iob), rdiag_tmp(iob))
                if (rloc_tmp(iob) == 0.0d0) cycle
              end if

              if (.not. reach_cutoff) then
                if (dist_tmp(iob) > dist_cutoff_fac_square) cycle
              end if

              nobsl_incr = nobsl_incr + 1
              nobs_use2(nobsl_incr) = iob
            end do
          end if ! [ nn_steps(icm+1) > nn_steps(icm) ]
        end do ! [ do icm = 1, n_merge(ic) ]

        if (LOG_LEVEL >= 3) then
          write (6, '(A,I4,A,F12.3,L2,2I8,10F12.3)') '--- TRY #', q, ': ', search_incr(1)*q, reach_cutoff, nn, nobsl_incr, search_incr(2:n_merge(ic))*q
        end if

        if (nobsl_incr >= nobsl_max_master) loop = .false.
      end do ! [ loop ]

      if (present(srch_q0)) then
        if (q == srch_q0(ic) .and. nobsl_incr > nobsl_max_master * 3) then ! when (nobsl_incr >= nobsl_max_master) too soon, decrease srch_q0
          srch_q0(ic) = q - 1
        else if (q > srch_q0(ic)) then ! when (nobsl_incr >= nobsl_max_master) too late, increase srch_q0
          srch_q0(ic) = q
        end if
      end if

      if (nobsl_incr == 0) cycle

      if (nobsl_incr > nobsl_max_master) then
        call QUICKSELECT_arg(dist_tmp, nobs_use2, 1, nobsl_incr, nobsl_max_master)
        nobsl_incr = nobsl_max_master
      end if

      if (LOG_LEVEL >= 3) then
        nobsl_cm(:) = 0
      end if

      do n = 1, nobsl_incr
        nobsl = nobsl + 1
        iob = nobs_use2(n)
        hdxf(nobsl,:) = obsda_sort%ensval(1:MEMBER,iob)
        rdiag(nobsl) = rdiag_tmp(iob)
        rloc(nobsl) = rloc_tmp(iob)
        dep(nobsl) = obsda_sort%val(iob)
        if (present(depd)) then
          depd(nobsl) = obsda_sort%ensval(mmdetobs,iob)
        end if

        if (LOG_LEVEL >= 3) then
          ic2 = ctype_elmtyp(uid_obs(obs(obsda_sort%set(iob))%elm(obsda_sort%idx(iob))), obs(obsda_sort%set(iob))%typ(obsda_sort%idx(iob)))
          do icm = 1, n_merge(ic)
            if (ic2 == ic_merge(icm,ic)) then
              nobsl_cm(icm) = nobsl_cm(icm) + 1
            end if
          end do
        end if
      end do

      if (LOG_LEVEL >= 3) then
        if (nobsl_incr == nobsl_max_master) then
          cutd = hori_loc_ctype(ic) * sqrt(dist_tmp(nobs_use2(nobsl_incr)))
        else
          cutd = hori_loc_ctype(ic) * dist_zero_fac
        end if
        write (6, '(A,F12.3,L2,8x,I8,10I8)') '--- FNL      : ', cutd, reach_cutoff, nobsl_incr, nobsl_cm(1:n_merge(ic))
      end if

      if (present(nobsl_t)) then
        nobsl_t(ielm_u_master,ityp_master) = nobsl_incr
      end if
      if (present(cutd_t)) then
        if (nobsl_incr == nobsl_max_master) then
          cutd_t(ielm_u_master,ityp_master) = hori_loc_ctype(ic) * sqrt(dist_tmp(nobs_use2(nobsl_incr)))
        end if
      end if

    !---------------------------------------------------------------------------
    else
    !---------------------------------------------------------------------------
    ! When obs number limit is enabled and the sorting criterion is NOT distance,
    ! conduct selection to all observations within the localization cutoff area.
    !---------------------------------------------------------------------------

      nn = 0
      nobsl_incr = 0

      do icm = 1, n_merge(ic)
        ic2 = ic_merge(icm,ic)

        if (obsgrd(ic2)%tot_ext > 0) then
          nn_prev = nn
          call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
          call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)

          do n = nn_prev+1, nn
            iob = nobs_use(n)

            call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ic2, ndist_dummy, rloc_tmp(iob), rdiag_tmp(iob))
            if (rloc_tmp(iob) == 0.0d0) cycle

            nobsl_incr = nobsl_incr + 1
            nobs_use2(nobsl_incr) = iob
          end do
        end if ! [ obsgrd(ic2)%tot_ext > 0 ]
      end do ! [ do icm = 1, n_merge(ic) ]

      if (nobsl_incr == 0) cycle

      if (nobsl_incr > nobsl_max_master) then
        if (MAX_NOBS_PER_GRID_CRITERION == 2) then
          call QUICKSELECT_desc_arg(rloc_tmp, nobs_use2, 1, nobsl_incr, nobsl_max_master)
        else if (MAX_NOBS_PER_GRID_CRITERION == 3) then
          call QUICKSELECT_arg(rdiag_tmp, nobs_use2, 1, nobsl_incr, nobsl_max_master)
        else
          write (6, '(A,I6)') "[Error] Unsupported 'MAX_NOBS_PER_GRID_CRITERION':", MAX_NOBS_PER_GRID_CRITERION
          stop 99
        end if
        nobsl_incr = nobsl_max_master
      end if

      do n = 1, nobsl_incr
        nobsl = nobsl + 1
        iob = nobs_use2(n)
        hdxf(nobsl,:) = obsda_sort%ensval(1:MEMBER,iob)
        rdiag(nobsl) = rdiag_tmp(iob)
        rloc(nobsl) = rloc_tmp(iob)
        dep(nobsl) = obsda_sort%val(iob)
        if (present(depd)) then
          depd(nobsl) = obsda_sort%ensval(mmdetobs,iob)
        end if
      end do

      if (present(nobsl_t)) then
        nobsl_t(ielm_u_master,ityp_master) = nobsl_incr
      end if
      if (present(cutd_t)) then
        if (nobsl_incr == nobsl_max_master) then
          if (MAX_NOBS_PER_GRID_CRITERION == 2) then
            cutd_t(ielm_u_master,ityp_master) = rloc_tmp(nobs_use2(nobsl_incr))
          else if (MAX_NOBS_PER_GRID_CRITERION == 3) then
            cutd_t(ielm_u_master,ityp_master) = rdiag_tmp(nobs_use2(nobsl_incr))
          end if
        end if
      end if

    !---------------------------------------------------------------------------
    end if

  end do ! [ ic = 1, nctype ]

  !-----------------------------------------------------------------------------
  ! Finalize
  !-----------------------------------------------------------------------------

#ifdef DEBUG
  if (nobsl > nobstotal) then
    write (6,'(A,I5,A,I5)') '[Error] nobsl=', nobsl, ' > nobstotal=', nobstotal
    write (6,*) 'ri,rj,lev,rz=', ri, rj, rlev, rz
    stop 99
  end if
#endif

  deallocate (nobs_use)
  if (maxval(MAX_NOBS_PER_GRID(:)) > 0) then
    deallocate (nobs_use2)
    deallocate (rloc_tmp)
    deallocate (rdiag_tmp)
    if (MAX_NOBS_PER_GRID_CRITERION == 1) then
      deallocate (dist_tmp)
    end if
  end if

  return
end subroutine obs_local

!-------------------------------------------------------------------------------
! Calculate the range of the rectangle that covers the (horizontal) localization
! cut-off length in the extended subdomain, given the observation type
!-------------------------------------------------------------------------------
subroutine obs_local_range(ctype, ri, rj, imin, imax, jmin, jmax)
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  implicit none
  integer, intent(in) :: ctype
  real(r_size), intent(in) :: ri, rj
  integer, intent(out) :: imin, imax, jmin, jmax

  real(r_size) :: dist_zero_i, dist_zero_j

  dist_zero_i = hori_loc_ctype(ctype) * dist_zero_fac / DX
  dist_zero_j = hori_loc_ctype(ctype) * dist_zero_fac / DY
  call ij_obsgrd_ext(ctype, ri - dist_zero_i, rj - dist_zero_j, imin, jmin)
  call ij_obsgrd_ext(ctype, ri + dist_zero_i, rj + dist_zero_j, imax, jmax)
#ifdef DEBUG
  if (imin < 1 .or. imax > obsgrd(ctype)%ngrdext_i .or. &
      jmin < 1 .or. jmax > obsgrd(ctype)%ngrdext_j) then
    write (6, '(A)') '[Error] The extended subdomain is not wide enough.'
    stop 99
  end if
#endif

  return
end subroutine obs_local_range

!-------------------------------------------------------------------------------
! Subroutine for main calculation of obs_local
!-------------------------------------------------------------------------------
subroutine obs_local_cal(ri, rj, rlev, rz, nvar, iob, ic, ndist, nrloc, nrdiag)
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  implicit none
  real(r_size), intent(in) :: ri, rj, rlev, rz ! coordinate of the targeted model grid
  integer, intent(in) :: nvar         ! index of targeted model variable
  integer, intent(in) :: iob          ! index of observation in obsda_sort
  integer, intent(in) :: ic           ! observation combined type
  real(r_size), intent(out) :: ndist  ! normalized 3D distance SQUARE       (in case of rejected obs: -1.)
  real(r_size), intent(out) :: nrloc  ! localization weight                 (in case of rejected obs:  0.)
  real(r_size), intent(out) :: nrdiag ! weighted observation error variance (in case of rejected obs: -1.)

  integer :: obelm           ! observation variable type
  integer :: obtyp           ! observation report type
  integer :: obset
  integer :: obidx
  real(r_size) :: rdx, rdy
  real(r_size) :: nd_h, nd_v ! normalized horizontal/vertical distances

  nrloc = 0.0d0
  nrdiag = -1.0d0
  ndist = -1.0d0

  obelm = elm_ctype(ic)
  obtyp = typ_ctype(ic)
  obset = obsda_sort%set(iob)
  obidx = obsda_sort%idx(iob)
#ifdef DEBUG
  if (obelm /= obs(obset)%elm(obidx)) then
    write (6, '(A)') '[Error] inconsistent observation variable type !!!'
    stop 99
  end if
  if (obtyp /= obs(obset)%typ(obidx)) then
    write (6, '(A)') '[Error] inconsistent observation report type !!!'
    stop 99
  end if
#endif
  !
  ! Calculate variable localization
  !
  if (nvar > 0) then  ! use variable localization only when nvar > 0
#ifdef DEBUG
    if (uid_obs_varlocal(obelm) <= 0) then
      write (6,'(A)') '[Error] unsupport observation type in variable localization.'
      stop 1
    end if
#endif
    nrloc = var_local(nvar,uid_obs_varlocal(obelm))

    !--- reject obs by variable localization
    if (nrloc < tiny(var_local)) then
      nrloc = 0.0d0
      return
    end if
  end if
  !
  ! Calculate normalized vertical distances
  !
  if (vert_loc_ctype(ic) == 0.0d0) then
    nd_v = 0.0d0                                                            ! no vertical localization
  else if (obelm == id_ps_obs) then
    nd_v = ABS(LOG(obs(obset)%dat(obidx)) - LOG(rlev)) / vert_loc_ctype(ic) ! for ps, use observed ps value for the base of vertical localization
  else if (obelm == id_rain_obs) then
    nd_v = ABS(LOG(VERT_LOCAL_RAIN_BASE) - LOG(rlev)) / vert_loc_ctype(ic)  ! for rain, use VERT_LOCAL_RAIN_BASE for the base of vertical localization
  else if (obtyp == 22) then ! obtypelist(obtyp) == 'PHARAD'
    nd_v = ABS(obs(obset)%lev(obidx) - rz) / vert_loc_ctype(ic)             ! for PHARAD, use z-coordinate for vertical localization
#ifdef H08
  else if (obtyp == 23) then ! obtypelist(obtyp) == 'H08IRB'                ! H08
    nd_v = ABS(LOG(obsda_sort%lev(iob)) - LOG(rlev)) / vert_loc_ctype(ic)   ! H08 for H08IRB, use obsda2%lev(iob) for the base of vertical localization
#endif
  else
    nd_v = ABS(LOG(obs(obset)%lev(obidx)) - LOG(rlev)) / vert_loc_ctype(ic)
  end if

  !--- reject obs by normalized vertical distance
  !    (do this first because there is large possibility to reject obs by the vertical distrance)
  if (nd_v > dist_zero_fac) then
    nrloc = 0.0d0
    return
  end if
  !
  ! Calculate normalized horizontal distances
  !
  rdx = (ri - obs(obset)%ri(obidx)) * DX
  rdy = (rj - obs(obset)%rj(obidx)) * DY
  nd_h = sqrt(rdx*rdx + rdy*rdy) / hori_loc_ctype(ic)

  !--- reject obs by normalized horizontal distance
  if (nd_h > dist_zero_fac) then
    nrloc = 0.0d0
    return
  end if
  !
  ! Calculate (normalized 3D distances)^2
  !
  ndist = nd_h * nd_h + nd_v * nd_v

  !--- reject obs by normalized 3D distance
  if (ndist > dist_zero_fac_square) then
    nrloc = 0.0d0
    ndist = -1.0d0
    return
  end if
  !
  ! Calculate observational localization
  !
  nrloc = nrloc * EXP(-0.5d0 * ndist)
  !
  ! Calculate (observation variance / localization)
  !
  nrdiag = obs(obset)%err(obidx) * obs(obset)%err(obidx) / nrloc

  return
end subroutine obs_local_cal

!-------------------------------------------------------------------------------
! Relaxation parameter based on grid locations (not for covariance inflation purpose)
!-------------------------------------------------------------------------------
subroutine relax_beta(ri, rj, rz, beta)
  use scale_atmos_grid_cartesC, only: &
    DX, DY
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO
  implicit none
  real(r_size), intent(in) :: ri, rj, rz
  real(r_size), intent(out) :: beta
  real(r_size) :: dist_bdy

  beta = 1.0d0
  !
  ! Upper bound of updates when RADAR_ZMAX is set and only radar observations are assimilated
  !
  if (radar_only .and. rz > RADAR_ZMAX + max(VERT_LOCAL(22), VERT_LOCAL_RADAR_VR) * dist_zero_fac) then
    beta = 0.0d0
    return
  end if
  !
  ! Boundary buffer
  !
  if (BOUNDARY_BUFFER_WIDTH > 0.0d0) then
    dist_bdy = min(min(ri-IHALO, nlong+IHALO+1-ri) * DX, &
                   min(rj-JHALO, nlatg+JHALO+1-rj) * DY) / BOUNDARY_BUFFER_WIDTH
#ifdef DEBUG
    if (dist_bdy < 0.0d0) then
      write (6, '(A,4F10.3)') '[Error] Wrong dist_bdy:', &
            ri-IHALO, nlong+IHALO+1-ri, rj-JHALO, nlatg+JHALO+1-rj
      stop 99
    end if
#endif
    if (dist_bdy < 1.0d0) then
      beta = max(dist_bdy, 0.0d0)
    end if
  end if

  return
end subroutine relax_beta

!-------------------------------------------------------------------------------
! Relaxation via LETKF weight - RTPP method
!-------------------------------------------------------------------------------
subroutine weight_RTPP(w, infl, wrlx)
  implicit none
  real(r_size), intent(in) :: w(MEMBER,MEMBER)
  real(r_size), intent(in) :: infl
  real(r_size), intent(out) :: wrlx(MEMBER,MEMBER)
  integer :: m

  wrlx = (1.0d0 - RELAX_ALPHA) * w
  do m = 1, MEMBER
    wrlx(m,m) = wrlx(m,m) + RELAX_ALPHA * sqrt(infl)
  end do

  return
end subroutine weight_RTPP

!-------------------------------------------------------------------------------
! Relaxation via LETKF weight - RTPS method
!-------------------------------------------------------------------------------
subroutine weight_RTPS(w, pa, xb, infl, wrlx, infl_out)
  implicit none
  real(r_size), intent(in) :: w(MEMBER,MEMBER)
  real(r_size), intent(in) :: pa(MEMBER,MEMBER)
  real(r_size), intent(in) :: xb(MEMBER)
  real(r_size), intent(in) :: infl
  real(r_size), intent(out) :: wrlx(MEMBER,MEMBER)
  real(r_size), intent(out) :: infl_out
  real(r_size) :: var_g, var_a
  integer :: m, k

  var_g = 0.0d0
  var_a = 0.0d0
  do m = 1, MEMBER
    var_g = var_g + xb(m) * xb(m)
    do k = 1, MEMBER
      var_a = var_a + xb(k) * pa(k,m) * xb(m)
    end do
  end do
  if (var_g > 0.0d0 .and. var_a > 0.0d0) then
    infl_out = RELAX_ALPHA_SPREAD * sqrt(var_g * infl / (var_a * real(MEMBER-1,r_size))) &  ! Whitaker and Hamill 2012
             - RELAX_ALPHA_SPREAD + 1.0d0                                                   !
!    infl_out = sqrt(RELAX_ALPHA_SPREAD * (var_g * infl / (var_a * real(MEMBER-1,r_size))) & ! Hamrud et al. 2015 (slightly modified)
!                  - RELAX_ALPHA_SPREAD + 1.0d0)                                             !
    wrlx = w * infl_out
  else
    wrlx = w
    infl_out = 1.0d0
  end if

  return
end subroutine weight_RTPS

END MODULE letkf_tools
