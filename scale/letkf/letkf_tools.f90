MODULE letkf_tools
!=======================================================================
!
! [PURPOSE:] Module for LETKF with GFS
!
! [HISTORY:]
!   01/26/2009 Takemasa Miyoshi  created
!   10/04/2012 Guo-Yuan Lien     modified for GFS model
!   07/01/2013 Daisuke Hotta     ported EFSO code from Y.Ohta's code
!   01/01/2014 Guo-Yuan Lien     merged to GFS-LETKF main development
!
!=======================================================================
  USE common
  use common_nml
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
  USE common_letkf

  USE letkf_obs
!  USE efso_tools

  use scale_precision, only: RP

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: das_letkf !, das_efso

  real(r_size),save :: var_local(nv3d+nv2d,nid_obs_varlocal)
  integer,save :: var_local_n2n(nv3d+nv2d)

  integer,save :: ctype_merge(nid_obs,nobtype)

CONTAINS
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE das_letkf(gues3d,gues2d,anal3d,anal2d)
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: gues3d(nij1,nlev,MEMBER,nv3d) ! background ensemble
  REAL(r_size),INTENT(INOUT) :: gues2d(nij1,MEMBER,nv2d)      !  output: destroyed
  REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,MEMBER,nv3d)   ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal2d(nij1,MEMBER,nv2d)

  REAL(r_size) :: mean3d(nij1,nlev,nv3d)
  REAL(r_size) :: mean2d(nij1,nv2d)
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

  REAL(r_size) :: parm
  REAL(r_size) :: trans(MEMBER,MEMBER,nv3d+nv2d)
  REAL(r_size) :: transm(MEMBER,nv3d+nv2d)       !GYL
  REAL(r_size) :: transrlx(MEMBER,MEMBER)        !GYL
  REAL(r_size) :: pa(MEMBER,MEMBER,nv3d+nv2d)    !GYL

  INTEGER :: ij,ilev,n,m,i,k,nobsl
  INTEGER :: nobsl_t(nid_obs,nobtype)            !GYL
  REAL(r_size) :: beta                           !GYL
  REAL(r_size) :: tmpinfl                        !GYL
  REAL(r_size) :: q_mean,q_sprd                  !GYL
  REAL(r_size) :: q_anal(MEMBER)                 !GYL

  WRITE(6,'(A)') 'Hello from das_letkf'
  WRITE(6,'(A,F15.2)') '  INFL_MUL = ',INFL_MUL

  WRITE(6,'(A,I8)') 'Target observation numbers (global) : NOBS=',nobstotalg
  WRITE(6,'(A,I8)') 'Target observation numbers processed in this subdomian : NOBS=',nobstotal
  !
  ! In case of no obs
  !
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
  var_local_n2n(1) = 1
  DO n=2,nv3d+nv2d
    DO i=1,n
      var_local_n2n(n) = i
      IF(MAXVAL(ABS(var_local(i,:)-var_local(n,:))) < TINY(var_local)) EXIT
    END DO
  END DO
  !
  ! Observation number limit (*to be moved to namelist*)
  !
  ctype_merge(:,:) = 0
  ctype_merge(uid_obs(id_radar_ref_obs),22) = 1
  ctype_merge(uid_obs(id_radar_ref_zero_obs),22) = 1
  !
  ! FCST PERTURBATIONS
  !
  CALL ensmean_grd(MEMBER,nij1,gues3d,gues2d,mean3d,mean2d)
  DO n=1,nv3d
    DO m=1,MEMBER
!$OMP PARALLEL DO PRIVATE(i,k)
      DO k=1,nlev
        DO i=1,nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - mean3d(i,k,n)
        END DO
      END DO
!$OMP END PARALLEL DO
    END DO
  END DO
  DO n=1,nv2d
    DO m=1,MEMBER
!$OMP PARALLEL DO PRIVATE(i)
      DO i=1,nij1
        gues2d(i,m,n) = gues2d(i,m,n) - mean2d(i,n)
      END DO
!$OMP END PARALLEL DO
    END DO
  END DO
  !
  ! multiplicative inflation
  !
  IF(INFL_MUL > 0.0d0) THEN  ! fixed multiplicative inflation parameter
    work3d = INFL_MUL
    work2d = INFL_MUL
  ELSE  ! 3D parameter values are read-in
    allocate (work3dg(nlon,nlat,nlev,nv3d))
    allocate (work2dg(nlon,nlat,nv2d))
    IF(myrank_e == lastmem_rank_e) THEN
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',INFL_MUL_IN_BASENAME,'.pe',proc2mem(2,1,myrank+1),'.nc'
      call read_restart(INFL_MUL_IN_BASENAME,work3dg,work2dg)
    END IF
    CALL scatter_grd_mpi(lastmem_rank_e,work3dg,work2dg,work3d,work2d)
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
    allocate (work3dn(nobtype,nij1,nlev,nv3d))
    allocate (work2dn(nobtype,nij1,nv2d))
    work3dn = 0.0d0
    work2dn = 0.0d0
  END IF
  !
  ! MAIN ASSIMILATION LOOP
  !
  ALLOCATE(hdxf (nobstotal,MEMBER))
  ALLOCATE(rdiag(nobstotal))
  ALLOCATE(rloc (nobstotal))
  ALLOCATE(dep  (nobstotal))

  DO ilev=1,nlev
    WRITE(6,'(A,I3,F18.3)') 'ilev = ',ilev, MPI_WTIME()

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ij,n,m,k,hdxf,rdiag,rloc,dep,nobsl,nobsl_t,parm,beta,trans,transm,transrlx,pa,tmpinfl,q_mean,q_sprd,q_anal)
    DO ij=1,nij1
!      WRITE(6,'(A,I3,A,I8,F18.3)') 'ilev = ',ilev, ', ij = ',ij, MPI_WTIME()

      ! update 3D variables
      DO n=1,nv3d

        ! calculate mean and perturbation weights
        IF(var_local_n2n(n) < n) THEN
          ! if weights already computed for other variables can be re-used(no variable localization), copy from there 
          trans(:,:,n) = trans(:,:,var_local_n2n(n))
          transm(:,n) = transm(:,var_local_n2n(n))                                     !GYL
          IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                         !GYL
            pa(:,:,n) = pa(:,:,var_local_n2n(n))                                       !GYL
          END IF                                                                       !GYL
          work3d(ij,ilev,n) = work3d(ij,ilev,var_local_n2n(n))
          IF(NOBS_OUT) THEN                                                            !GYL
            work3dn(:,ij,ilev,n) = work3dn(:,ij,ilev,var_local_n2n(n))                 !GYL
          END IF                                                                       !GYL
        ELSE
          ! compute weights with localized observations
          CALL obs_local(rig1(ij),rjg1(ij),mean3d(ij,ilev,iv3d_p),hgt1(ij,ilev),n,hdxf,rdiag,rloc,dep,nobsl,nobsl_t=nobsl_t)
          IF(RELAX_TO_INFLATED_PRIOR) THEN                                             !GYL
            parm = work3d(ij,ilev,n)                                                   !GYL
          ELSE                                                                         !GYL
            parm = 1.0d0                                                               !GYL
          END IF                                                                       !GYL
          IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                         !GYL
            CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & !GYL
                            trans(:,:,n),transm=transm(:,n),pao=pa(:,:,n),   &         !GYL
                            rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)           !GYL
          ELSE                                                                         !GYL
            CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & !GYL
                            trans(:,:,n),transm=transm(:,n),                 &         !GYL
                            rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)           !GYL
          END IF                                                                       !GYL
          IF(NOBS_OUT) THEN                                                            !GYL
            work3dn(:,ij,ilev,n) = real(sum(nobsl_t, dim=1),r_size)                    !GYL
            work3dn(21,ij,ilev,n) = real(nobsl_t(9,22),r_size)                         !GYL !!! addtionally save ref nobs in a special place
          END IF                                                                       !GYL
        END IF

        ! weight parameter based on grid locations (not for cov inflation purpose)     !GYL
        CALL relax_beta(rig1(ij),rjg1(ij),mean3d(ij,ilev,iv3d_p),n,beta)               !GYL

        IF(beta == 0.0d0) THEN                                                         !GYL
          ! no analysis update needed
          anal3d(ij,ilev,:,n) = mean3d(ij,ilev,n) + gues3d(ij,ilev,:,n)                !GYL
        ELSE                                                                           !GYL
          ! relaxation via LETKF weight
          IF(RELAX_ALPHA /= 0.0d0) THEN                                                !GYL - RTPP method (Zhang et al. 2004)
            CALL weight_RTPP(trans(:,:,n),parm,transrlx)                               !GYL
          ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                    !GYL - RTPS method (Whitaker and Hamill 2012)
            IF(RELAX_SPREAD_OUT) THEN                                                  !GYL
              CALL weight_RTPS(trans(:,:,n),pa(:,:,n),gues3d(ij,ilev,:,n), &           !GYL
                               parm,transrlx,work3da(ij,ilev,n))                       !GYL
            ELSE                                                                       !GYL
              CALL weight_RTPS(trans(:,:,n),pa(:,:,n),gues3d(ij,ilev,:,n), &           !GYL
                               parm,transrlx,tmpinfl)                                  !GYL
            END IF                                                                     !GYL
          ELSE                                                                         !GYL
            transrlx = trans(:,:,n)                                                    !GYL - No relaxation
          END IF                                                                       !GYL

          ! total weight matrix
          DO m=1,MEMBER                                                                !GYL
            DO k=1,MEMBER                                                              !GYL
              transrlx(k,m) = (transrlx(k,m) + transm(k,n)) * beta                     !GYL
            END DO                                                                     !GYL
            transrlx(m,m) = transrlx(m,m) + (1.0d0-beta)                               !GYL
          END DO                                                                       !GYL

          ! analysis update
          DO m=1,MEMBER
            anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)
            DO k=1,MEMBER
              anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &                              !GYL
                                  + gues3d(ij,ilev,k,n) * transrlx(k,m)                !GYL
            END DO  
          END DO
        END IF ! [ beta == 0.0d0 ]                                                     !GYL

        ! limit q spread
        IF(Q_SPRD_MAX > 0.0d0 .and. n == iv3d_q) THEN                                  !GYL
          q_mean = SUM(anal3d(ij,ilev,:,n)) / REAL(MEMBER,r_size)                      !GYL
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

      END DO ! [ n=1,nv3d ]

      ! update 2D variables at ilev = 1
      IF(ilev == 1) THEN 
        DO n=1,nv2d

          ! calculate mean and perturbation weights
          IF(var_local_n2n(nv3d+n) < nv3d+n) THEN
            ! if weights already computed for other variables can be re-used(no variable localization), copy from there 
            trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
            transm(:,nv3d+n) = transm(:,var_local_n2n(nv3d+n))                         !GYL
            IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                       !GYL
              pa(:,:,nv3d+n) = pa(:,:,var_local_n2n(nv3d+n))                           !GYL
            END IF                                                                     !GYL
            IF(var_local_n2n(nv3d+n) <= nv3d) THEN                                     !GYL - correct the bug of the 2d variable update
              work2d(ij,n) = work3d(ij,ilev,var_local_n2n(nv3d+n))                     !GYL
              IF(NOBS_OUT) THEN                                                        !GYL
                work2dn(:,ij,n) = work3dn(:,ij,ilev,var_local_n2n(nv3d+n))             !GYL
              END IF                                                                   !GYL
            ELSE                                                                       !GYL
              work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n)-nv3d)                     !GYL
              IF(NOBS_OUT) THEN                                                        !GYL
                work2dn(:,ij,n) = work2dn(:,ij,var_local_n2n(nv3d+n)-nv3d)             !GYL
              END IF                                                                   !GYL
            END IF                                                                     !GYL
          ELSE
            ! compute weights with localized observations
            CALL obs_local(rig1(ij),rjg1(ij),mean3d(ij,ilev,iv3d_p),hgt1(ij,ilev),nv3d+n,hdxf,rdiag,rloc,dep,nobsl,nobsl_t=nobsl_t)
            IF(RELAX_TO_INFLATED_PRIOR) THEN                                           !GYL
              parm = work2d(ij,n)                                                      !GYL
            ELSE                                                                       !GYL
              parm = 1.0d0                                                             !GYL
            END IF                                                                     !GYL
            IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                       !GYL
              CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & !GYL
                              trans(:,:,nv3d+n),transm=transm(:,nv3d+n),pao=pa(:,:,nv3d+n), & !GYL
                              rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)         !GYL
            ELSE                                                                       !GYL
              CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & !GYL
                              trans(:,:,nv3d+n),transm=transm(:,nv3d+n),       &       !GYL
                              rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)         !GYL
            END IF                                                                     !GYL
            IF(NOBS_OUT) THEN                                                          !GYL
              work2dn(:,ij,n) = real(sum(nobsl_t,dim=1),r_size)                        !GYL
              work2dn(21,ij,n) = real(nobsl_t(9,22),r_size)                            !GYL !!! addtionally save ref nobs in a special place
            END IF                                                                     !GYL
          END IF

          ! weight parameter based on grid locations (not for cov inflation purpose)   !GYL
          CALL relax_beta(rig1(ij),rjg1(ij),mean3d(ij,ilev,iv3d_p),nv3d+n,beta)        !GYL

          IF(beta == 0.0d0) THEN                                                       !GYL
            ! no analysis update needed
            anal2d(ij,:,n) = mean2d(ij,n) + gues2d(ij,:,n)                             !GYL
          ELSE                                                                         !GYL
            ! relaxation via LETKF weight
            IF(RELAX_ALPHA /= 0.0d0) THEN                                              !GYL - RTPP method (Zhang et al. 2004)
              CALL weight_RTPP(trans(:,:,nv3d+n),parm,transrlx)                        !GYL
            ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                  !GYL - RTPS method (Whitaker and Hamill 2012)
              IF(RELAX_SPREAD_OUT) THEN                                                !GYL
                CALL weight_RTPS(trans(:,:,nv3d+n),pa(:,:,nv3d+n),gues2d(ij,:,n), &    !GYL
                                 parm,transrlx,work2da(ij,n))                          !GYL
              ELSE                                                                     !GYL
                CALL weight_RTPS(trans(:,:,nv3d+n),pa(:,:,nv3d+n),gues2d(ij,:,n), &    !GYL
                                 parm,transrlx,tmpinfl)                                !GYL
              END IF                                                                   !GYL
            ELSE                                                                       !GYL
              transrlx = trans(:,:,nv3d+n)                                             !GYL - No relaxation
            END IF                                                                     !GYL

            ! total weight matrix
            DO m=1,MEMBER                                                              !GYL
              DO k=1,MEMBER                                                            !GYL
                transrlx(k,m) = (transrlx(k,m) + transm(k,nv3d+n)) * beta              !GYL
              END DO                                                                   !GYL
              transrlx(m,m) = transrlx(m,m) + (1.0d0-beta)                             !GYL
            END DO                                                                     !GYL

            ! analysis update
            DO m=1,MEMBER
              anal2d(ij,m,n) = mean2d(ij,n)
              DO k=1,MEMBER
                anal2d(ij,m,n) = anal2d(ij,m,n) &                                      !GYL - sum trans and transm here
                               + gues2d(ij,k,n) * transrlx(k,m)                        !GYL
              END DO
            END DO
          END IF ! [ beta == 0.0d0 ]                                                   !GYL

        END DO ! [ n=1,nv2d ]
      END IF ! [ ilev == 1 ]

    END DO ! [ ij=1,nij1 ]
!$OMP END PARALLEL DO

  END DO ! [ ilev=1,nlev ]

  DEALLOCATE(hdxf,rdiag,rloc,dep)
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
    if (.not. allocated(work3dg)) allocate (work3dg(nlon,nlat,nlev,nv3d))
    if (.not. allocated(work2dg)) allocate (work2dg(nlon,nlat,nv2d))
    CALL gather_grd_mpi(lastmem_rank_e,work3d,work2d,work3dg,work2dg)
    IF(myrank_e == lastmem_rank_e) THEN
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',INFL_MUL_OUT_BASENAME,'.pe',proc2mem(2,1,myrank+1),'.nc'
      call write_restart(INFL_MUL_OUT_BASENAME,work3dg,work2dg)
    END IF
  END IF
  !
  ! Write inflation parameter (in analysis) corresponding to the RTPS method
  !
  IF(RELAX_SPREAD_OUT) THEN
    if (.not. allocated(work3dg)) allocate (work3dg(nlon,nlat,nlev,nv3d))
    if (.not. allocated(work2dg)) allocate (work2dg(nlon,nlat,nv2d))
    CALL gather_grd_mpi(lastmem_rank_e,work3da,work2da,work3dg,work2dg)
    IF(myrank_e == lastmem_rank_e) THEN
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',RELAX_SPREAD_OUT_BASENAME,'.pe',proc2mem(2,1,myrank+1),'.nc'
      call write_restart(RELAX_SPREAD_OUT_BASENAME,work3dg,work2dg)
    END IF
    DEALLOCATE(work3da,work2da)
  END IF
  !
  ! Write observation numbers
  !
  IF(NOBS_OUT) THEN
    if (.not. allocated(work3dg)) allocate (work3dg(nlon,nlat,nlev,nv3d))
    if (.not. allocated(work2dg)) allocate (work2dg(nlon,nlat,nv2d))
    work3d(:,:,1) = work3dn(1,:,:,iv3d_t)  !!! Assuming variable localization is not used so that obs numbers used are the same over variables,
    work3d(:,:,2) = work3dn(3,:,:,iv3d_t)  !!! use "variable dimenstion" to save obs numbers of different observation types
    work3d(:,:,3) = work3dn(4,:,:,iv3d_t)
    work3d(:,:,4) = work3dn(8,:,:,iv3d_t)
    work3d(:,:,5) = work3dn(21,:,:,iv3d_t)
    work3d(:,:,6) = work3dn(22,:,:,iv3d_t)
    CALL gather_grd_mpi(lastmem_rank_e,work3d,work2d,work3dg,work2dg)
    IF(myrank_e == lastmem_rank_e) THEN
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',NOBS_OUT_BASENAME,'.pe',proc2mem(2,1,myrank+1),'.nc'
      call write_restart(NOBS_OUT_BASENAME,work3dg,work2dg)
    END IF
    DEALLOCATE(work3dn,work2dn)
  END IF
  IF (allocated(work3dg)) deallocate (work3dg)
  IF (allocated(work2dg)) deallocate (work2dg)
  !
  ! Additive inflation
  !
  IF(INFL_ADD > 0.0d0) THEN
    CALL read_ens_mpi(INFL_ADD_IN_BASENAME,gues3d,gues2d)
    CALL ensmean_grd(MEMBER,nij1,gues3d,gues2d,work3d,work2d)
    DO n=1,nv3d
      DO m=1,MEMBER
!$OMP PARALLEL DO PRIVATE(i,k)
        DO k=1,nlev
          DO i=1,nij1
            gues3d(i,k,m,n) = gues3d(i,k,m,n) - work3d(i,k,n)
          END DO
        END DO
!$OMP END PARALLEL DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,MEMBER
!$OMP PARALLEL DO PRIVATE(i)
        DO i=1,nij1
          gues2d(i,m,n) = gues2d(i,m,n) - work2d(i,n)
        END DO
!$OMP END PARALLEL DO
      END DO
    END DO

    WRITE(6,'(A)') '===== Additive covariance inflation ====='
    WRITE(6,'(A,F10.4)') '  parameter:',INFL_ADD
    WRITE(6,'(A)') '========================================='
!    parm = 0.7d0
!    DO ilev=1,nlev
!      parm_infl_damp(ilev) = 1.0d0 + parm &
!        & + parm * REAL(1-ilev,r_size)/REAL(nlev_dampinfl,r_size)
!      parm_infl_damp(ilev) = MAX(parm_infl_damp(ilev),1.0d0)
!    END DO
    DO n=1,nv3d
      DO m=1,MEMBER
!$OMP PARALLEL DO PRIVATE(ij,ilev)
        DO ilev=1,nlev
          DO ij=1,nij1
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
              & + gues3d(ij,ilev,m,n) * INFL_ADD
          END DO
        END DO
!$OMP END PARALLEL DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,MEMBER
!$OMP PARALLEL DO PRIVATE(ij)
        DO ij=1,nij1
          anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,m,n) * INFL_ADD
        END DO
!$OMP END PARALLEL DO
      END DO
    END DO
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
! [OUT]
!   hdxf    : fcstast ensemble perturbations in the observation space
!   rdiag   : localization-weighted observation error variances
!   rloc    : localization weights
!   dep     : observation departure
!   nobsl   : number of valid observations (in hdxf, rdiag, rloc, dep)
!   nobsl_t : (optional) number of valid observations wrt. observation variables/types
!-------------------------------------------------------------------------------
subroutine obs_local(ri, rj, rlev, rz, nvar, hdxf, rdiag, rloc, dep, nobsl, nobsl_t)
  use common_sort
  use scale_grid, only: &
    DX, DY
  use scale_rm_process, only: &
    PRC_NUM_X, &
    PRC_NUM_Y
  implicit none
  real(r_size), intent(in) :: ri, rj, rlev, rz
  integer, intent(in) :: nvar
  real(r_size), intent(out) :: hdxf(nobstotal,member)
  real(r_size), intent(out) :: rdiag(nobstotal)
  real(r_size), intent(out) :: rloc(nobstotal)
  real(r_size), intent(out) :: dep(nobstotal)
  integer, intent(out) :: nobsl
  integer, intent(out), optional :: nobsl_t(nid_obs,nobtype)

  integer, allocatable :: nobs_use(:)

  integer, allocatable :: iob_tmp(:)
  real(r_size), allocatable :: dist_tmp(:)
  real(r_size), allocatable :: rloc_tmp(:)
  real(r_size), allocatable :: rdiag_tmp(:)

!  real(r_size), allocatable :: dist_t(:,:)
!  real(r_size), allocatable :: rloc_t(:,:)
!  real(r_size), allocatable :: rdiag_t(:,:)
!  integer, allocatable :: isort_t(:,:)
!  integer, allocatable :: iob_t(:,:)

  real(r_size) :: ndist, nrloc, nrdiag
  integer :: iob, ityp, ielm, ielm_u
  integer :: imin, imax, jmin, jmax
  integer :: ic, ic2, icm
  integer :: n, nn, nn_prev
  integer :: nobsl_prev, nobsl_incr
  integer :: nobsl_max_master
  integer :: ielm_u_master, ityp_master

  logical :: ctype_skip(nctype)
  integer :: ic_merge(nid_obs*nobtype)
  integer :: n_merge

  integer :: q
  logical :: loop
  integer :: nn_steps(nctype+1)
  logical :: reach_cutoff
  integer :: imin_cutoff, imax_cutoff, jmin_cutoff, jmax_cutoff
  real(r_size) :: search_incr, search_incr_i, search_incr_j

!  logical :: condition

  !-----------------------------------------------------------------------------

  if (nobstotal == 0) then
    nobsl = 0
    if (present(nobsl_t)) then
      nobsl_t(:,:) = 0
    end if  
  end if

  !-----------------------------------------------------------------------------
  ! Initialize
  !-----------------------------------------------------------------------------


  if (maxval(MAX_NOBS_PER_GRID(:)) > 0) then

    allocate (nobs_use (nobstotal))
    allocate (iob_tmp  (nobstotal))
    allocate (dist_tmp (nobstotal))
    allocate (rloc_tmp (nobstotal))
    allocate (rdiag_tmp(nobstotal))


!    allocate (isort_t(maxval(MAX_NOBS_PER_GRID(:)), nid_obs))
!    allocate (iob_t  (maxval(MAX_NOBS_PER_GRID(:)), nid_obs))
!    allocate (rloc_t (maxval(MAX_NOBS_PER_GRID(:)), nid_obs))
!    allocate (rdiag_t(maxval(MAX_NOBS_PER_GRID(:)), nid_obs))
!    if (MAX_NOBS_PER_GRID_CRITERION == 1) then
!      allocate (dist_t(maxval(MAX_NOBS_PER_GRID(:)), nid_obs))
!    end if

  else

    allocate (nobs_use(maxnobs_per_ctype))

  end if

  nobsl = 0

  !-----------------------------------------------------------------------------
  ! For each observation type,
  ! do rough data search by a rectangle using the sorting mesh, and then
  ! do precise data search by normalized 3D distance and variable localization
  !-----------------------------------------------------------------------------

  ctype_skip(:) = .false.

  do ic = 1, nctype
    if (ctype_skip(ic)) cycle

    n_merge = 1
    ic_merge(1) = ic
    if (ctype_merge(elm_u_ctype(ic),typ_ctype(ic)) > 0) then
      do ic2 = ic+1, nctype
        if (ctype_merge(elm_u_ctype(ic2),typ_ctype(ic2)) == ctype_merge(elm_u_ctype(ic),typ_ctype(ic))) then
          n_merge = n_merge + 1
          ic_merge(n_merge) = ic2
          ctype_skip(ic2) = .true.
!          write(6, '(9A)') '[Info] Observation number limit: Consider obs types (', obtypelist(typ_ctype(ic)), ', ', obelmlist(elm_u_ctype(ic)), &
!                           ') and (', obtypelist(typ_ctype(ic2)), ', ', obelmlist(elm_u_ctype(ic2)), ') together'
        end if
      end do
    end if

    nobsl_max_master = MAX_NOBS_PER_GRID(typ_ctype(ic)) ! Use the number limit setting of the "master" obs type for all group of obs types
    ielm_u_master = elm_u_ctype(ic)                     ! Count observation numbers    at the "master" obs type for all group of obs types
    ityp_master = typ_ctype(ic)                         !

    if (nobsl_max_master <= 0) then
    !-------------------------------------------------------------------------
    ! When obs number limit is not enabled,
    ! Directly prepare (hdxf, dep, rdiag, rloc) output.
    !-------------------------------------------------------------------------

      nobsl_prev = nobsl

      do icm = 1, n_merge
        ic2 = ic_merge(icm)

        if (obsgrd(ic2)%tot_ext > 0) then
          ielm = elm_ctype(ic2)
          ityp = typ_ctype(ic2)

          nn = 0
          call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
          call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)

          do n = 1, nn
            iob = nobs_use(n)

            call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ielm, ityp, ndist, nrloc, nrdiag)
            if (nrloc == 0.0d0) cycle

            nobsl = nobsl + 1
            hdxf(nobsl,:) = obsda2%ensval(:,iob)
            dep(nobsl) = obsda2%val(iob)
            rdiag(nobsl) = nrdiag
            rloc(nobsl) = nrloc
          end do ! [ n = 1, nn ]
        end if ! [ obsgrd(ic2)%tot_ext > 0 ]
      end do ! [ do icm = 1, n_merge ]

      if (present(nobsl_t)) then
        nobsl_t(ielm_u_master,ityp_master) = nobsl - nobsl_prev
      end if

    !-------------------------------------------------------------------------
    else if (MAX_NOBS_PER_GRID_CRITERION == 1) then
    !-------------------------------------------------------------------------
    ! XXX
    !-------------------------------------------------------------------------

      nn = 0
      do icm = 1, n_merge
        ic2 = ic_merge(icm)

        if (obsgrd(ic2)%tot_ext > 0) then
          call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
          call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)
        end if ! [ obsgrd(ic2)%tot_ext > 0 ]
      end do ! [ do icm = 1, n_merge ]


write (6, '(A,14x,I8)') '--- ALL      : ', nn


      if (nn == 0) cycle




      search_incr = max(obsgrd(ic)%grdspc_i, obsgrd(ic)%grdspc_j)
      search_incr_i = search_incr / DX
      search_incr_j = search_incr / DY

      nobsl_incr = 0


      q = 0
      loop = .true.
      do while (loop)
        q = q + 1


        nn = 0


        reach_cutoff = .true.

        do icm = 1, n_merge
          ic2 = ic_merge(icm)

          nn_steps(icm) = nn

          if (obsgrd(ic2)%tot_ext > 0) then
            call ij_obsgrd_ext(ic2, ri-search_incr_i*q, rj-search_incr_j*q, imin, jmin)
            call ij_obsgrd_ext(ic2, ri+search_incr_i*q, rj+search_incr_j*q, imax, jmax)

            call obs_local_range(ic2, ri, rj, imin_cutoff, imax_cutoff, jmin_cutoff, jmax_cutoff)
            if (imin < imin_cutoff .or. imax > imax_cutoff .or. &
                jmin < jmin_cutoff .or. jmax > jmax_cutoff) then
              imin = imin_cutoff
              imax = imax_cutoff
              jmin = jmin_cutoff
              jmax = jmax_cutoff
            else
              reach_cutoff = .false.
            end if

            call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)
          end if ! [ obsgrd(ic2)%tot_ext > 0 ]
        end do ! [ do icm = 1, n_merge ]

        nn_steps(n_merge+1) = nn

write (6, '(A,I4,A,F12.3,L2,I8)') '--- Try #', q, ': ', search_incr*q, reach_cutoff, nn

        if ((.not. reach_cutoff) .and. nn < nobsl_max_master) cycle

        if (reach_cutoff) then
          loop = .false.
          if (nn == 0) exit
        end if

        nobsl_incr = 0

        do icm = 1, n_merge
          ic2 = ic_merge(icm)

          if (nn_steps(icm+1) > nn_steps(icm)) then
            ielm = elm_ctype(ic2)
            ityp = typ_ctype(ic2)

            do n = nn_steps(icm)+1, nn_steps(icm+1)
              iob = nobs_use(n)

              if (reach_cutoff) then
                call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ielm, ityp, ndist, nrloc, nrdiag)
              else
                call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ielm, ityp, ndist, nrloc, nrdiag, cutoff=search_incr*q)
              end if
              if (nrloc == 0.0d0) cycle

              nobsl_incr = nobsl_incr + 1
              iob_tmp(nobsl_incr) = iob
              dist_tmp(nobsl_incr) = ndist
              rloc_tmp(nobsl_incr) = nrloc
              rdiag_tmp(nobsl_incr) = nrdiag
            end do
          end if ! [ nn_steps(icm+1) > nn_steps(icm) ]
        end do ! [ do icm = 1, n_merge ]

write (6, '(A,I4,A,F12.3,L2,2I8)') '--- Try #', q, ': ', search_incr*q, reach_cutoff, nn, nobsl_incr

        if (nobsl_incr >= nobsl_max_master) loop = .false.

      end do ! [ loop ]

      if (nobsl_incr == 0) cycle

      if (nobsl_incr > nobsl_max_master) then ! there is a small chance that nobsl_incr = nobsl_max_master, which does not need to perform selection
        call QUICKSELECT(dist_tmp, 1, nobsl_incr, nobsl_max_master, B=rloc_tmp, C=rdiag_tmp, I=iob_tmp)
        nobsl_incr = nobsl_max_master
      end if

      rloc(nobsl+1:nobsl+nobsl_incr) = rloc_tmp(1:nobsl_incr)
      rdiag(nobsl+1:nobsl+nobsl_incr) = rdiag_tmp(1:nobsl_incr)
      do n = 1, nobsl_incr
        nobsl = nobsl + 1
        iob = iob_tmp(n)
        hdxf(nobsl,:) = obsda2%ensval(:,iob)
        dep(nobsl) = obsda2%val(iob)
      end do

      if (present(nobsl_t)) then
        nobsl_t(ielm_u_master,ityp_master) = nobsl_incr
      end if

    !-------------------------------------------------------------------------
    else
    !-------------------------------------------------------------------------
    ! XXX
    !-------------------------------------------------------------------------

      nn = 0
      nobsl_incr = 0
      do icm = 1, n_merge
        ic2 = ic_merge(icm)

        if (obsgrd(ic2)%tot_ext > 0) then
          ielm = elm_ctype(ic2)
          ityp = typ_ctype(ic2)

          nn_prev = nn
          call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
          call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)

          do n = nn_prev+1, nn
            iob = nobs_use(n)

            call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ielm, ityp, ndist, nrloc, nrdiag)
            if (nrloc == 0.0d0) cycle

            nobsl_incr = nobsl_incr + 1
            iob_tmp(nobsl_incr) = iob
            rloc_tmp(nobsl_incr) = nrloc
            rdiag_tmp(nobsl_incr) = nrdiag
          end do
        end if ! [ obsgrd(ic2)%tot_ext > 0 ]
      end do ! [ do icm = 1, n_merge ]

      if (nobsl_incr == 0) cycle

      if (nobsl_incr > nobsl_max_master) then
        call QUICKSELECT(rdiag_tmp, 1, nobsl_incr, nobsl_max_master, B=rloc_tmp, I=iob_tmp)
!!!!!! only valid for MAX_NOBS_PER_GRID_CRITERION = 3. MAX_NOBS_PER_GRID_CRITERION = 2 needs to be considered.
        nobsl_incr = nobsl_max_master
      end if

      rloc(nobsl+1:nobsl+nobsl_incr) = rloc_tmp(1:nobsl_incr)
      rdiag(nobsl+1:nobsl+nobsl_incr) = rdiag_tmp(1:nobsl_incr)
      do n = 1, nobsl_incr
        nobsl = nobsl + 1
        iob = iob_tmp(n)
        hdxf(nobsl,:) = obsda2%ensval(:,iob)
        dep(nobsl) = obsda2%val(iob)
      end do

      if (present(nobsl_t)) then
        nobsl_t(ielm_u_master,ityp_master) = nobsl_incr
      end if

    !-------------------------------------------------------------------------
    end if

  end do ! [ ic = 1, nctype ]



!!!!!!!!      else
!!!!!!!!      !-------------------------------------------------------------------------
!!!!!!!!      ! When obs number limit is enabled,
!!!!!!!!      ! Save only the observations within the number limit in temporary arrays,
!!!!!!!!      ! and prepare (hdxf, dep, rdiag, rloc) output later.
!!!!!!!!      !-------------------------------------------------------------------------

!!!!!!        isort_t(:,:) = 0
!!!!!!        nobsl_t_(:) = 0

!!!!!!        call obs_local_range(ityp, ri, rj, imin, imax, jmin, jmax)

!!!!!!        search_incr_i = obsgrd(ityp)%grdspc_i / DX
!!!!!!        search_incr_j = obsgrd(ityp)%grdspc_j / DY

!!!!!!        q = 0
!!!!!!        do
!!!!!!          call ij_obsgrd_ext(ityp, ri-search_incr_i*q, rj-search_incr_j*q, imin_tmp, jmin_tmp)
!!!!!!          call ij_obsgrd_ext(ityp, ri+search_incr_i*q, rj+search_incr_j*q, imax_tmp, jmax_tmp)

!!!!!!          if (imin_tmp < imin .or. imax_tmp > imax .or. &
!!!!!!              jmin_tmp < jmin .or. jmax_tmp > jmax) then
!!!!!!            XXXXXX
!!!!!!          end if
!!!!!!          nn = 0
!!!!!!          call obs_choose_ext(ityp, imin, imax, jmin, jmax, nn, nobs_use)

!!!!!!!        call QUICKSELECT(MAX_NOBS_PER_GRID(ityp), obsda2%val(1:nn), nn, n)

!!!!!!        ! Sort the first [MAX_NOBS_PER_GRID(ityp)] observations in terms of:
!!!!!!        !  MAX_NOBS_PER_GRID_CRITERION = 1: smallest ndist
!!!!!!        !  MAX_NOBS_PER_GRID_CRITERION = 2: largest nrloc
!!!!!!        !  MAX_NOBS_PER_GRID_CRITERION = 3: smallest nrdiag

!!!!!!        do n = 1, nn  ! loop over observations within the search rectangle

!!!!!!          iob = nobs_use(n)
!!!!!!          ielm_u = uid_obs(obs(obsda2%set(iob))%elm(obsda2%idx(iob)))

!!!!!!          call obs_local_cal(ri, rj, rlev, rz, nvar, iob, ndist, nrloc, nrdiag)
!!!!!!          if (nrdiag < 0.0d0) cycle ! observation rejected XXX

!!!!!!          !---------------------------------------------------------------------
!!!!!!          ! Case 0: If the number limit has been reached and the priority of
!!!!!!          !         this obs is lower than all of the obs in the current set of
!!!!!!          !         choice, skip right away.
!!!!!!          !---------------------------------------------------------------------
!!!!!!          if (nobsl_t_(ielm_u) >= MAX_NOBS_PER_GRID(ityp)) then
!!!!!!            condition = .false.
!!!!!!            if (MAX_NOBS_PER_GRID_CRITERION == 1) then
!!!!!!              if (ndist >= dist_t(isort_t(MAX_NOBS_PER_GRID(ityp),ielm_u),ielm_u)) condition = .true.
!!!!!!            else if (MAX_NOBS_PER_GRID_CRITERION == 2) then
!!!!!!              if (nrloc <= rloc_t(isort_t(MAX_NOBS_PER_GRID(ityp),ielm_u),ielm_u)) condition = .true.
!!!!!!            else if (MAX_NOBS_PER_GRID_CRITERION == 3) then
!!!!!!              if (nrdiag >= rdiag_t(isort_t(MAX_NOBS_PER_GRID(ityp),ielm_u),ielm_u)) condition = .true.
!!!!!!            end if
!!!!!!            if (condition) then
!!!!!!              cycle
!!!!!!            end if
!!!!!!          end if

!!!!!!          do s = 1, MAX_NOBS_PER_GRID(ityp)

!!!!!!            !-------------------------------------------------------------------
!!!!!!            ! Case 1: This obs is of the last priority,
!!!!!!            !         but the number limit has NOT been reached,
!!!!!!            !         save this obs in the spare space of the temporary arrays.
!!!!!!            !-------------------------------------------------------------------
!!!!!!            if (isort_t(s,ielm_u) == 0) then
!!!!!!              nobsl_t_(ielm_u) = nobsl_t_(ielm_u) + 1
!!!!!!              isort_t(s,ielm_u) = nobsl_t_(ielm_u)

!!!!!!              iob_t(isort_t(s,ielm_u),ielm_u) = iob  ! iob_t: indices to retrieve ensval(:,:) and val(:) later
!!!!!!                                                     !        [do not create the potentially very large ensval_t(:,:,:) array to save ensval(:)]
!!!!!!              rloc_t(isort_t(s,ielm_u),ielm_u) = nrloc
!!!!!!              rdiag_t(isort_t(s,ielm_u),ielm_u) = nrdiag
!!!!!!              if (MAX_NOBS_PER_GRID_CRITERION == 1) then
!!!!!!                dist_t(isort_t(s,ielm_u),ielm_u) = ndist
!!!!!!              end if
!!!!!!              exit  ! case matched, exit the loop
!!!!!!            !-------------------------------------------------------------------
!!!!!!            else
!!!!!!            !-------------------------------------------------------------------
!!!!!!              condition = .false.
!!!!!!              if (MAX_NOBS_PER_GRID_CRITERION == 1) then
!!!!!!                if (ndist < dist_t(isort_t(s,ielm_u),ielm_u)) condition = .true.
!!!!!!              else if (MAX_NOBS_PER_GRID_CRITERION == 2) then
!!!!!!                if (nrloc > rloc_t(isort_t(s,ielm_u),ielm_u)) condition = .true.
!!!!!!              else if (MAX_NOBS_PER_GRID_CRITERION == 3) then
!!!!!!                if (nrdiag < rdiag_t(isort_t(s,ielm_u),ielm_u)) condition = .true.
!!!!!!              end if
!!!!!!              if (condition) then
!!!!!!              !-----------------------------------------------------------------
!!!!!!                if (nobsl_t_(ielm_u) < MAX_NOBS_PER_GRID(ityp)) then
!!!!!!                !---------------------------------------------------------------
!!!!!!                ! Case 2: This obs has the priority higher than some obs in the current set of choice
!!!!!!                !         and the number limit has NOT been reached,
!!!!!!                !         save this obs in the spare space of the temporary arrays,
!!!!!!                !         and shift the sorting index array accordingly.
!!!!!!                !---------------------------------------------------------------
!!!!!!                  nobsl_t_(ielm_u) = nobsl_t_(ielm_u) + 1
!!!!!!                  do ss = nobsl_t_(ielm_u), s+1, -1
!!!!!!                    isort_t(ss,ielm_u) = isort_t(ss-1,ielm_u)
!!!!!!                  end do
!!!!!!                  isort_t(s,ielm_u) = nobsl_t_(ielm_u)

!!!!!!                  iob_t(isort_t(s,ielm_u),ielm_u) = iob
!!!!!!                  rloc_t(isort_t(s,ielm_u),ielm_u) = nrloc
!!!!!!                  rdiag_t(isort_t(s,ielm_u),ielm_u) = nrdiag
!!!!!!                  if (MAX_NOBS_PER_GRID_CRITERION == 1) then
!!!!!!                    dist_t(isort_t(s,ielm_u),ielm_u) = ndist
!!!!!!                  end if
!!!!!!                !---------------------------------------------------------------
!!!!!!                else
!!!!!!                !---------------------------------------------------------------
!!!!!!                ! Case 3: This obs has the priority higher than some obs in the current set of choice
!!!!!!                !         and the number limit has been reached,
!!!!!!                !         save this obs by overwriting the temporary arrays at where the obs of the last priority is,
!!!!!!                !         and shift the sorting index array accordingly.
!!!!!!                !---------------------------------------------------------------
!!!!!!                  tmpisort = isort_t(MAX_NOBS_PER_GRID(ityp),ielm_u)
!!!!!!                  do ss = MAX_NOBS_PER_GRID(ityp), s+1, -1
!!!!!!                    isort_t(ss,ielm_u) = isort_t(ss-1,ielm_u)
!!!!!!                  end do
!!!!!!                  isort_t(s,ielm_u) = tmpisort

!!!!!!                  iob_t(isort_t(s,ielm_u),ielm_u) = iob
!!!!!!                  rloc_t(isort_t(s,ielm_u),ielm_u) = nrloc
!!!!!!                  rdiag_t(isort_t(s,ielm_u),ielm_u) = nrdiag
!!!!!!                  if (MAX_NOBS_PER_GRID_CRITERION == 1) then
!!!!!!                    dist_t(isort_t(s,ielm_u),ielm_u) = ndist
!!!!!!                  end if
!!!!!!                !---------------------------------------------------------------
!!!!!!                end if
!!!!!!                exit  ! case matched, exit the loop

!!!!!!              !-----------------------------------------------------------------
!!!!!!              end if
!!!!!!            !-------------------------------------------------------------------
!!!!!!            ! Otherwise, skip this obs because of its too low priority.
!!!!!!            !-------------------------------------------------------------------
!!!!!!            end if
!!!!!!          end do ! [ s = 1, MAX_NOBS_PER_GRID(ityp) ]
!!!!!!  
!!!!!!        end do ! [ n = 1, nn ]

!!!!!!        !
!!!!!!        ! prepare (hdxf, dep, rdiag, rloc) output from the sort result
!!!!!!        !

!!!!!!        do ielm_u = 1, nid_obs
!!!!!!          do s = 1, nobsl_t_(ielm_u)
!!!!!!            nobsl = nobsl + 1
!!!!!!            iob = iob_t(s,ielm_u)

!!!!!!            hdxf(nobsl,:) = obsda2%ensval(:,iob)
!!!!!!            dep(nobsl) = obsda2%val(iob)
!!!!!!            rloc(nobsl) = rloc_t(s,ielm_u)
!!!!!!            rdiag(nobsl) = rdiag_t(s,ielm_u)
!!!!!!          end do ! [ s = 1, nobsl_t_(ielm_u) ]
!!!!!!        end do ! [ ielm_u = 1, nid_obs ]

!!!!!!        if (present(nobsl_t)) nobsl_t(:,ityp) = nobsl_t_(:)

!!!!!!        !-----------------------------------------------------------------------
!!!!!!      end if ! [ MAX_NOBS_PER_GRID(ityp) <= 0 ]

  !-----------------------------------------------------------------------------
  ! Finalize
  !-----------------------------------------------------------------------------

  if (nobsl > nobstotal) then
    write (6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=', nobsl, ' > NOBSTOTAL=', nobstotal
    write (6,*) 'RI,RJ,LEV,NOBSL,NOBSTOTAL=', ri, rj, rlev, rz, nobsl, nobstotal
    stop 99
  end if

  deallocate (nobs_use)

  if (maxval(MAX_NOBS_PER_GRID(:)) > 0) then

    deallocate (iob_tmp)
    deallocate (dist_tmp)
    deallocate (rloc_tmp)
    deallocate (rdiag_tmp)

!    deallocate (isort_t)
!    deallocate (iob_t)
!    deallocate (rloc_t)
!    deallocate (rdiag_t)
!    if (MAX_NOBS_PER_GRID_CRITERION == 1) then
!      deallocate (dist_t)
!    end if
  end if

  return
end subroutine obs_local

!-------------------------------------------------------------------------------
! Calculate the range of the rectangle that covers the (horizontal) localization
! cut-off length in the extended subdomain, given the observation type
!-------------------------------------------------------------------------------
subroutine obs_local_range(ctype, ri, rj, imin, imax, jmin, jmax)
  use scale_grid, only: &
    DX, DY
  implicit none
  integer, intent(in) :: ctype
  real(r_size), intent(in) :: ri, rj
  integer, intent(out) :: imin, imax, jmin, jmax

  real(r_size) :: hori_loc
  real(r_size) :: dist_zero_i, dist_zero_j

  hori_loc = HORI_LOCAL(typ_ctype(ctype))
  if (elm_ctype(ctype) == id_radar_ref_zero_obs) then
    hori_loc = HORI_LOCAL_RADAR_OBSNOREF
  end if

  dist_zero_i = hori_loc * dist_zero_fac / DX
  dist_zero_j = hori_loc * dist_zero_fac / DY
  call ij_obsgrd_ext(ctype, ri - dist_zero_i, rj - dist_zero_j, imin, jmin)
  call ij_obsgrd_ext(ctype, ri + dist_zero_i, rj + dist_zero_j, imax, jmax)
!  imin = max(1, imin)
!  imax = min(obsgrd(ctype)%ngrdext_i, imax)
!  jmin = max(1, jmin)
!  jmax = min(obsgrd(ctype)%ngrdext_j, jmax)
  if (imin < 1 .or. imax > obsgrd(ctype)%ngrdext_i .or. &
      jmin < 1 .or. jmax > obsgrd(ctype)%ngrdext_j) then
    write (6, '(A)') '[Error] The extended subdomain is not wide enough.'
    stop 99
  end if

  return
end subroutine obs_local_range

!-------------------------------------------------------------------------------
! Subroutine for main calculation of obs_local
!-------------------------------------------------------------------------------
subroutine obs_local_cal(ri, rj, rlev, rz, nvar, iob, obelm, obtyp, ndist, nrloc, nrdiag, cutoff)
  use scale_grid, only: &
    DX, DY
  implicit none
  real(r_size), intent(in) :: ri, rj, rlev, rz ! coordinate of the targeted model grid
  integer, intent(in) :: nvar         ! index of targeted model variable
  integer, intent(in) :: iob          ! index of observation in obsda2
  integer, intent(in) :: obelm        ! observation variable type
  integer, intent(in) :: obtyp        ! observation report type
  real(r_size), intent(out) :: ndist  ! normalized 3D distance SQUARE
  real(r_size), intent(out) :: nrloc  ! localization weight
  real(r_size), intent(out) :: nrdiag ! weighted observation error variance
  real(r_size), intent(in), optional :: cutoff

  real(r_size) :: dist_zero_fac_cutoff
  integer :: obset, obidx
  real(r_size) :: rdx, rdy
  real(r_size) :: nd_h, nd_v ! normalized horizontal/vertical distances

  nrloc = 0.0d0
  nrdiag = -1.0d0
  ndist = -1.0d0

  obset = obsda2%set(iob)
  obidx = obsda2%idx(iob)
  if (obelm /= obs(obset)%elm(obidx)) then
    write (6, '(A)') '[Error] inconsistent observation variable type !!!'
    stop 99
  end if
  if (obtyp /= obs(obset)%typ(obidx)) then
    write (6, '(A)') '[Error] inconsistent observation report type !!!'
    stop 99
  end if
  !
  ! Calculate variable localization
  !
  if (nvar > 0) then  ! use variable localization only when nvar > 0
    if (uid_obs_varlocal(obelm) <= 0) then
      write (6,'(A)') '[Error] unsupport observation type in variable localization.'
      stop 1
    end if
    nrloc = var_local(nvar,uid_obs_varlocal(obelm))

    !--- reject obs by variable localization
    if (nrloc < tiny(var_local)) then
      nrloc = 0.0d0
      return
    end if
  end if
  !
  ! Pre-process the cutoff distance
  !
  if (present(cutoff)) then
    if (obelm == id_radar_ref_zero_obs) then
      dist_zero_fac_cutoff = cutoff / HORI_LOCAL_RADAR_OBSNOREF
    else
      dist_zero_fac_cutoff = cutoff / HORI_LOCAL(obtyp)
    end if
  else
    dist_zero_fac_cutoff = dist_zero_fac
  end if
  !
  ! Calculate normalized vertical distances
  !
  if (VERT_LOCAL(obtyp) == 0.0d0) then
    nd_v = 0.0d0                                                           ! no vertical localization
  else if (obelm == id_ps_obs) then
    nd_v = ABS(LOG(obs(obset)%dat(obidx)) - LOG(rlev)) / VERT_LOCAL(obtyp) ! for ps, use observed ps value for the base of vertical localization
  else if (obelm == id_rain_obs) then
    nd_v = ABS(LOG(VERT_LOCAL_RAIN_BASE) - LOG(rlev)) / VERT_LOCAL(obtyp)  ! for rain, use VERT_LOCAL_RAIN_BASE for the base of vertical localization
  else if (obtyp == 22) then ! obtypelist(obtyp) == 'PHARAD'
    nd_v = ABS(obs(obset)%lev(obidx) - rz) / VERT_LOCAL(obtyp)             ! for PHARAD, use z-coordinate for vertical localization
#ifdef H08
  else if (obtyp == 23) then ! obtypelist(obtyp) == 'H08IRB'               ! H08
    nd_v = ABS(LOG(obsda2%lev(iob)) - LOG(rlev)) / VERT_LOCAL(obtyp)       ! H08 for H08IRB, use obsda2%lev(iob) for the base of vertical localization
#endif
  else
    nd_v = ABS(LOG(obs(obset)%lev(obidx)) - LOG(rlev)) / VERT_LOCAL(obtyp)
  end if

  !--- reject obs by normalized vertical distance
  !    (do this first because there is large possibility to reject obs by the vertical distrance)
  if (nd_v > dist_zero_fac_cutoff) then
    nrloc = 0.0d0
    return
  end if
  !
  ! Calculate normalized horizontal distances
  !
  rdx = (ri - obsda2%ri(iob)) * DX
  rdy = (rj - obsda2%rj(iob)) * DY
  if (obelm == id_radar_ref_zero_obs) then
    nd_h = sqrt(rdx*rdx + rdy*rdy) / HORI_LOCAL_RADAR_OBSNOREF
  else
    nd_h = sqrt(rdx*rdx + rdy*rdy) / HORI_LOCAL(obtyp)
  end if

  !--- reject obs by normalized horizontal distance
  if (nd_h > dist_zero_fac_cutoff) then
    nrloc = 0.0d0
    return
  end if
  !
  ! Calculate (normalized 3D distances)^2
  !
  ndist = nd_h * nd_h + nd_v * nd_v

  !--- reject obs by normalized 3D distance
  if (present(cutoff)) then
    if (ndist > dist_zero_fac_cutoff * dist_zero_fac_cutoff) then
      nrloc = 0.0d0
      return
    end if
  else
    if (ndist > dist_zero_fac_square) then
      nrloc = 0.0d0
      return
    end if
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
subroutine relax_beta(ri, rj, rlev, nvar, beta)
  use scale_grid, only: &
    DX, DY
  use scale_grid_index, only: &
    IHALO, JHALO
  implicit none
  real(r_size), intent(in) :: ri, rj, rlev
  integer, intent(in) :: nvar
  real(r_size), intent(out) :: beta
  real(r_size) :: dist_bdy

  beta = 1.0d0
  !
  ! Upper-limit of Q update levels
  !
  if (rlev < Q_UPDATE_TOP) then
    if (nvar >= iv3d_q .and. nvar <= iv3d_qg) then
      beta = 0.0d0
      return
    end if
  end if
  !
  ! Boundary buffer
  !
  if (BOUNDARY_BUFFER_WIDTH > 0.0d0) then
    dist_bdy = min(min(ri-IHALO, nlong+IHALO+1-ri) * DX, &
                   min(rj-JHALO, nlatg+JHALO+1-rj) * DY) / BOUNDARY_BUFFER_WIDTH
!    if (dist_bdy < 0.0d0) then
!      write (6, '(A,4F10.3)') '[Error] Wrong dist_bdy:', &
!            ri-IHALO, nlong+IHALO+1-ri, rj-JHALO, nlatg+JHALO+1-rj
!      stop 1
!    end if
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
