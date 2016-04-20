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
!  USE efso_nml
!  USE efso_tools

  use scale_precision, only: RP

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: das_letkf !, das_efso

  real(r_size),save :: var_local(nv3d+nv2d,nid_obs_varlocal)
  integer,save :: var_local_n2n(nv3d+nv2d)

CONTAINS
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE das_letkf(gues3d,gues2d,anal3d,anal2d)
  IMPLICIT NONE
  CHARACTER(8) :: inflfile='infl_mul'
  CHARACTER(20) :: inflfile_0='infl_mul.pe000000.nc'
  REAL(r_size),INTENT(INOUT) :: gues3d(nij1,nlev,MEMBER,nv3d) ! background ensemble
  REAL(r_size),INTENT(INOUT) :: gues2d(nij1,MEMBER,nv2d)      !  output: destroyed
  REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,MEMBER,nv3d) ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal2d(nij1,MEMBER,nv2d)
  REAL(r_size),ALLOCATABLE :: mean3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: mean2d(:,:)
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size),ALLOCATABLE :: work3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: work2d(:,:)
  REAL(r_size),ALLOCATABLE :: work3da(:,:,:)     !GYL
  REAL(r_size),ALLOCATABLE :: work2da(:,:)       !GYL

  REAL(r_size),ALLOCATABLE :: work3dl(:,:,:)
  REAL(r_size),ALLOCATABLE :: work2dl(:,:)

  REAL(RP),ALLOCATABLE :: work3dg(:,:,:,:)
  REAL(RP),ALLOCATABLE :: work2dg(:,:,:)
  REAL(r_size) :: parm
  REAL(r_size) :: trans(MEMBER,MEMBER,nv3d+nv2d)
  REAL(r_size) :: transm(MEMBER,nv3d+nv2d)       !GYL
  REAL(r_size) :: transrlx(MEMBER,MEMBER)        !GYL
  REAL(r_size) :: pa(MEMBER,MEMBER,nv3d+nv2d)    !GYL
  REAL(r_size) :: q_mean,q_sprd                  !GYL
  REAL(r_size) :: q_anal(MEMBER)                 !GYL
!  LOGICAL :: ex
!  INTEGER :: ij,ilev,n,m,i,j,k,nobsl,ierr,iret
  INTEGER :: ij,ilev,n,m,i,k,nobsl

  WRITE(6,'(A)') 'Hello from das_letkf'
  WRITE(6,'(A,F15.2)') '  COV_INFL_MUL = ',COV_INFL_MUL

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
  var_local(:,1) = VAR_LOCAL_UV(1:nv3d+nv2d)
  var_local(:,2) = VAR_LOCAL_T(1:nv3d+nv2d)
  var_local(:,3) = VAR_LOCAL_Q(1:nv3d+nv2d)
  var_local(:,4) = VAR_LOCAL_PS(1:nv3d+nv2d)
  var_local(:,5) = VAR_LOCAL_RAIN(1:nv3d+nv2d)
  var_local(:,6) = VAR_LOCAL_TC(1:nv3d+nv2d)
  var_local(:,7) = VAR_LOCAL_RADAR_REF(1:nv3d+nv2d)
  var_local(:,8) = VAR_LOCAL_RADAR_VR(1:nv3d+nv2d)
  var_local(:,9) = VAR_LOCAL_H08(1:nv3d+nv2d) ! H08
  var_local_n2n(1) = 1
  DO n=2,nv3d+nv2d
    DO i=1,n
      var_local_n2n(n) = i
      IF(MAXVAL(ABS(var_local(i,:)-var_local(n,:))) < TINY(var_local)) EXIT
    END DO
  END DO
  !
  ! FCST PERTURBATIONS
  !
  ALLOCATE(mean3d(nij1,nlev,nv3d))
  ALLOCATE(mean2d(nij1,nv2d))
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
  IF(COV_INFL_MUL > 0.0d0) THEN ! fixed multiplicative inflation parameter
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    work3d = COV_INFL_MUL
    work2d = COV_INFL_MUL
  END IF
  IF(COV_INFL_MUL <= 0.0d0) THEN ! 3D parameter values are read-in
    ALLOCATE( work3dg(nlon,nlat,nlev,nv3d) )
    ALLOCATE( work2dg(nlon,nlat,nv2d) )
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    IF(myrank_e == lastmem_rank_e) THEN

      IF(ADAPTIVE_INFL_INIT) THEN
        work3dg = -1.0d0 * COV_INFL_MUL
        work2dg = -1.0d0 * COV_INFL_MUL
      ELSE
!        WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',inflfile,'.pe',proc2mem(2,1,myrank+1),'.nc'
        call read_restart(inflfile,work3dg,work2dg)
!        call state_trans(work3dg)
      END IF

!      INQUIRE(FILE=inflfile_0,EXIST=ex)
!      IF(ex) THEN
!!        WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',inflfile,'.pe',proc2mem(2,1,myrank+1),'.nc'
!        call read_restart(inflfile,work3dg,work2dg)
!!        call state_trans(work3dg)
!      ELSE
!        WRITE(6,'(2A)') '!!WARNING: no such file exist: ',inflfile
!        work3dg = -1.0d0 * COV_INFL_MUL
!        work2dg = -1.0d0 * COV_INFL_MUL
!      END IF
    END IF
    CALL scatter_grd_mpi(lastmem_rank_e,work3dg,work2dg,work3d,work2d)
  END IF
  !
  ! RTPS relaxation
  !
  IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN
    ALLOCATE( work3da(nij1,nlev,nv3d) )
    ALLOCATE( work2da(nij1,nv2d) )
    work3da = 1.0d0
    work2da = 1.0d0
  END IF

  ALLOCATE( work3dl(nij1,nlev,nv3d) )
  ALLOCATE( work2dl(nij1,nv2d) )
  work3dl = 1.0d0
  work2dl = 1.0d0

  !
  ! MAIN ASSIMILATION LOOP
  !
  ALLOCATE( hdxf(1:nobstotal,1:MEMBER),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
  DO ilev=1,nlev
    WRITE(6,'(A,I3,F18.3)') 'ilev = ',ilev, MPI_WTIME()
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ij,n,hdxf,rdiag,rloc,dep,nobsl,parm,trans,transm,transrlx,pa,m,k,q_mean,q_sprd,q_anal)
    DO ij=1,nij1
!WRITE(6,'(A,I3,A,I8,F18.3)') 'ilev = ',ilev, ', ij = ',ij, MPI_WTIME()
      DO n=1,nv3d
        IF(var_local_n2n(n) < n) THEN
          trans(:,:,n) = trans(:,:,var_local_n2n(n))
          transm(:,n) = transm(:,var_local_n2n(n))                                     !GYL
          IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                         !GYL
            pa(:,:,n) = pa(:,:,var_local_n2n(n))                                       !GYL
          END IF                                                                       !GYL
          work3d(ij,ilev,n) = work3d(ij,ilev,var_local_n2n(n))

          work3dl(ij,ilev,n) = work3dl(ij,ilev,var_local_n2n(n))

        ELSE
          CALL obs_local(rig1(ij),rjg1(ij),mean3d(ij,ilev,iv3d_p),hgt1(ij,ilev),n,hdxf,rdiag,rloc,dep,nobsl)
          parm = work3d(ij,ilev,n)
          IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                         !GYL
            CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &         !GYL
                            trans(:,:,n),transm=transm(:,n),pao=pa(:,:,n),minfl=MIN_INFL_MUL) !GYL
          ELSE                                                                         !GYL
            CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &         !GYL
                            trans(:,:,n),transm=transm(:,n),minfl=MIN_INFL_MUL)        !GYL
          END IF                                                                       !GYL
          work3d(ij,ilev,n) = parm

          work3dl(ij,ilev,n) = real(nobsl,r_size)

        END IF
        IF((n == iv3d_q .OR. n == iv3d_qc .OR. n == iv3d_qr .OR. n == iv3d_qi .OR. n == iv3d_qs .OR. n == iv3d_qg) &
           .AND. ilev > LEV_UPDATE_Q) THEN !GYL, do not update upper-level q,qc
          anal3d(ij,ilev,:,n) = mean3d(ij,ilev,n) + gues3d(ij,ilev,:,n)                !GYL
        ELSE                                                                           !GYL
          IF(RELAX_ALPHA /= 0.0d0) THEN                                                !GYL - RTPP method (Zhang et al. 2005)
            CALL weight_RTPP(trans(:,:,n),transrlx)                                    !GYL
          ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                    !GYL - RTPS method (Whitaker and Hamill 2012)
            CALL weight_RTPS(trans(:,:,n),pa(:,:,n),gues3d(ij,ilev,:,n),transrlx,work3da(ij,ilev,n)) !GYL
          ELSE                                                                         !GYL
            transrlx = trans(:,:,n)                                                    !GYL
          END IF                                                                       !GYL
          DO m=1,MEMBER
            anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)
            DO k=1,MEMBER
              anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &                              !GYL - sum trans and transm here
                & + gues3d(ij,ilev,k,n) * (transrlx(k,m) + transm(k,n))                !GYL
            END DO
          END DO
        END IF                                                                         !GYL
        IF(n == iv3d_q .AND. ilev <= LEV_UPDATE_Q) THEN                                !GYL - limit the lower-level q spread
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
      IF(ilev == 1) THEN !update 2d variable at ilev=1
        DO n=1,nv2d
          IF(var_local_n2n(nv3d+n) < nv3d+n) THEN
            trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
            transm(:,nv3d+n) = transm(:,var_local_n2n(nv3d+n))                         !GYL
            IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                       !GYL
              pa(:,:,nv3d+n) = pa(:,:,var_local_n2n(nv3d+n))                           !GYL
            END IF                                                                     !GYL
            IF(var_local_n2n(nv3d+n) <= nv3d) THEN                                     !GYL - correct the bug of the 2d variable update
              work2d(ij,n) = work3d(ij,ilev,var_local_n2n(nv3d+n))                     !GYL

              work2dl(ij,n) = work3dl(ij,ilev,var_local_n2n(nv3d+n))                     !GYL

            ELSE                                                                       !GYL
              work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n)-nv3d)                     !GYL

              work2dl(ij,n) = work2dl(ij,var_local_n2n(nv3d+n)-nv3d)                     !GYL

            END IF                                                                     !GYL
          ELSE
            CALL obs_local(rig1(ij),rjg1(ij),mean3d(ij,ilev,iv3d_p),hgt1(ij,ilev),nv3d+n,hdxf,rdiag,rloc,dep,nobsl)
            parm = work2d(ij,n)
            IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                       !GYL
              CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &       !GYL
                              trans(:,:,nv3d+n),transm=transm(:,nv3d+n),pao=pa(:,:,nv3d+n),minfl=MIN_INFL_MUL) !GYL
            ELSE                                                                       !GYL
              CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &       !GYL
                              trans(:,:,nv3d+n),transm=transm(:,nv3d+n),minfl=MIN_INFL_MUL) !GYL
            END IF                                                                     !GYL
            work2d(ij,n) = parm

            work2dl(ij,n) = real(nobsl,r_size)

          END IF
          IF(RELAX_ALPHA /= 0.0d0) THEN                                                !GYL - RTPP method (Zhang et al. 2005)
            CALL weight_RTPP(trans(:,:,nv3d+n),transrlx)                               !GYL
          ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                    !GYL - RTPS method (Whitaker and Hamill 2012)
            CALL weight_RTPS(trans(:,:,nv3d+n),pa(:,:,nv3d+n),gues2d(ij,:,n),transrlx,work2da(ij,n)) !GYL
          ELSE                                                                         !GYL
            transrlx = trans(:,:,nv3d+n)                                               !GYL
          END IF                                                                       !GYL
          DO m=1,MEMBER
            anal2d(ij,m,n) = mean2d(ij,n)
            DO k=1,MEMBER
              anal2d(ij,m,n) = anal2d(ij,m,n) &                                        !GYL - sum trans and transm here
                & + gues2d(ij,k,n) * (transrlx(k,m) + transm(k,nv3d+n))                !GYL
            END DO
          END DO
        END DO ! [ n=1,nv2d ]
      END IF ! [ ilev == 1 ]
    END DO ! [ ij=1,nij1 ]
!$OMP END PARALLEL DO
  END DO ! [ ilev=1,nlev ]
  DEALLOCATE(hdxf,rdiag,rloc,dep)
!  !
!  ! Compute analyses of observations (Y^a)
!  !
!  IF(obsanal_output) THEN
!    call das_letkf_obs(work3dg,work2dg)
!  END IF
  !
  ! Write updated inflation parameters
  !
  IF(COV_INFL_MUL < 0.0d0) THEN
    CALL gather_grd_mpi(lastmem_rank_e,work3d,work2d,work3dg,work2dg)
    IF(myrank_e == lastmem_rank_e) THEN
!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',inflfile,'.pe',proc2mem(2,1,myrank+1),'.nc'
!      call state_trans_inv(work3dg)
      call write_restart(inflfile,work3dg,work2dg)
    END IF
    DEALLOCATE(work3d,work2d)
  END IF
  !
  ! Write inflation parameter (in analysis) corresponding to the RTPS method
  !
  IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN
    if (.not. allocated(work3dg)) allocate(work3dg(nlon,nlat,nlev,nv3d))
    if (.not. allocated(work2dg)) allocate(work2dg(nlon,nlat,nv2d))
    CALL gather_grd_mpi(lastmem_rank_e,work3da,work2da,work3dg,work2dg)
    IF(myrank_e == lastmem_rank_e) THEN
      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',INFL_OUT_BASENAME,'.pe',proc2mem(2,1,myrank+1),'.nc'
!      call state_trans_inv(work3dg)
      call write_restart(INFL_OUT_BASENAME,work3dg,work2dg)
    END IF
    DEALLOCATE(work3da,work2da)
  END IF

  if (.not. allocated(work3dg)) allocate(work3dg(nlon,nlat,nlev,nv3d))
  if (.not. allocated(work2dg)) allocate(work2dg(nlon,nlat,nv2d))
  CALL gather_grd_mpi(lastmem_rank_e,work3dl,work2dl,work3dg,work2dg)
  IF(myrank_e == lastmem_rank_e) THEN
    WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',NOBS_OUT_BASENAME,'.pe',proc2mem(2,1,myrank+1),'.nc'
!    call state_trans_inv(work3dg)
    call write_restart(NOBS_OUT_BASENAME,work3dg,work2dg)
  END IF
  DEALLOCATE(work3dl,work2dl)

  IF (allocated(work3dg)) deallocate(work3dg)
  IF (allocated(work2dg)) deallocate(work2dg)
  !
  ! Additive inflation
  !
  IF(SP_INFL_ADD > 0.0d0) THEN
    CALL read_ens_mpi('addi',gues3d,gues2d)
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
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

    DEALLOCATE(work3d,work2d)
    WRITE(6,'(A)') '===== Additive covariance inflation ====='
    WRITE(6,'(A,F10.4)') '  parameter:',SP_INFL_ADD
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
              & + gues3d(ij,ilev,m,n) * SP_INFL_ADD
          END DO
        END DO
!$OMP END PARALLEL DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,MEMBER
!$OMP PARALLEL DO PRIVATE(ij)
        DO ij=1,nij1
          anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,m,n) * SP_INFL_ADD
        END DO
!$OMP END PARALLEL DO
      END DO
    END DO
  END IF

  DEALLOCATE(mean3d,mean2d)
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
!    CALL letkf_core(MEMBER,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans,MIN_INFL_MUL,RELAX_ALPHA)

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
!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
! -- modified, using (rlon,rlat,rlev) instead of (ij,ilev), Guo-Yuan Lien
! -- optional oindex output, followed by D.Hotta
!-----------------------------------------------------------------------
SUBROUTINE obs_local(ri,rj,rlev,rz,nvar,hdxf,rdiag,rloc,dep,nobsl)
  use scale_grid, only: &
    DX, DY
  use scale_grid_index, only: &
    IHALO, JHALO
  use scale_les_process, only: &
    PRC_NUM_X, &
    PRC_NUM_Y

  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: ri,rj,rlev,rz
  INTEGER,INTENT(IN) :: nvar
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,MEMBER)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: nd_h,nd_v ! normalized horizontal/vertical distances
  REAL(r_size) :: ndist     ! normalized 3D distance SQUARE
  INTEGER,ALLOCATABLE:: nobs_use(:)
  integer :: ip,imin1,imax1,jmin1,jmax1,imin2,imax2,jmin2,jmax2
  integer :: iproc,jproc
  integer :: iidx, ityp
  integer :: ielm, ielm_u, ielm_varlocal
  integer :: n,nn,iob
  integer :: s,ss
  real(r_size) :: rdx,rdy


  real(r_size), allocatable :: ndist_t(:,:,:)
  integer, allocatable :: ip_t(:,:,:)
  integer, allocatable :: iob_t(:,:,:)
  integer, allocatable :: isort_t(:,:,:)
  integer :: nobsl_t(nid_obs,nobtype)


  real(r_size) :: radar_ri(OBS_IN_NUM)      !!!!!!
  real(r_size) :: radar_rj(OBS_IN_NUM)      !!!!!!
  real(r_size) :: dist_radar_edge, dist_bdy !!!!!!
  integer :: iof                            !!!!!!



  if (RADAR_EDGE_TAPER_WIDTH > 0.0d0 .and. RADAR_RANGE > 0.0d0) then                !!!!!!
    do iof = 1, OBS_IN_NUM                                                          !!!!!!
      if (OBS_IN_FORMAT(iof) == 2) then                                             !!!!!!
        call phys2ij(obs(iof)%meta(1),obs(iof)%meta(2),radar_ri(iof),radar_rj(iof)) !!!!!!
      end if                                                                        !!!!!!
    end do                                                                          !!!!!!
  end if                                                                            !!!!!!



!
! INITIALIZE
!
  nobsl_t(:,:) = 0

  IF( maxval(nobsgrd(nlon,nlat,:)) > 0 ) THEN
    ALLOCATE(nobs_use(maxval(nobsgrd(nlon,nlat,:))))
  END IF

  if (MAX_NOBS_PER_GRID > 0) then
    allocate(ndist_t(MAX_NOBS_PER_GRID, nid_obs, nobtype))
    allocate(ip_t(MAX_NOBS_PER_GRID, nid_obs, nobtype))
    allocate(iob_t(MAX_NOBS_PER_GRID, nid_obs, nobtype))
    allocate(isort_t(MAX_NOBS_PER_GRID, nid_obs, nobtype))
    isort_t(:,:,:) = 0
  end if

!
! data search
!
  imin1 = max(1, floor(ri - dlon_zero))
  imax1 = min(PRC_NUM_X*nlon, ceiling(ri + dlon_zero))
  jmin1 = max(1, floor(rj - dlat_zero))
  jmax1 = min(PRC_NUM_Y*nlat, ceiling(rj + dlat_zero))

  nobsl = 0
  nobsl_t(:,:) = 0

  do ip = 0, MEM_NP-1

    if (obsda2(ip)%nobs > 0) then

      call rank_1d_2d(ip, iproc, jproc)

      imin2 = max(1, imin1 - iproc*nlon)
      imax2 = min(nlon, imax1 - iproc*nlon)
      jmin2 = max(1, jmin1 - jproc*nlat)
      jmax2 = min(nlat, jmax1 - jproc*nlat)

      nn = 0
      call obs_choose(imin2,imax2,jmin2,jmax2,ip,nn,nobs_use)

!write(6,'(A,6I8)') '$$$==', imin2,imax2,jmin2,jmax2,ip,nn

      do n = 1, nn

        iidx = obsda2(ip)%idx(nobs_use(n))
        ielm = obs(obsda2(ip)%set(nobs_use(n)))%elm(iidx)
        ityp = obs(obsda2(ip)%set(nobs_use(n)))%typ(iidx)
        ielm_u = uid_obs(ielm)

        !
        ! localization
        !
!print *, '@@@', nobs_use(n), iidx, ielm

        rdx = (ri - obsda2(ip)%ri(nobs_use(n))) * DX
        rdy = (rj - obsda2(ip)%rj(nobs_use(n))) * DY
        nd_h = sqrt(rdx*rdx + rdy*rdy)

        ! calculate normalized horizontal/vertical distances
        select case (ielm)
        case (id_ps_obs)
          nd_h = nd_h / SIGMA_OBS
          nd_v = ABS(LOG(obs(obsda2(ip)%set(nobs_use(n)))%dat(iidx)) - LOG(rlev)) / SIGMA_OBSV
        case (id_rain_obs)
          nd_h = nd_h / SIGMA_OBS_RAIN
          nd_v = ABS(LOG(BASE_OBSV_RAIN) - LOG(rlev)) / SIGMA_OBSV_RAIN
        case (id_radar_ref_obs, id_radar_vr_obs, id_radar_prh_obs)
          if (ielm == id_radar_ref_obs .and. obs(obsda2(ip)%set(nobs_use(n)))%dat(iidx) <= RADAR_REF_THRES_DBZ+1.0d-6) then
            nd_h = nd_h / SIGMA_OBS_RADAR_OBSNOREF
          else
            nd_h = nd_h / SIGMA_OBS_RADAR
          end if
          nd_v = ABS(obs(obsda2(ip)%set(nobs_use(n)))%lev(iidx) - rz) / SIGMA_OBSZ_RADAR
        case (id_tclon_obs, id_tclat_obs, id_tcmip_obs)
          nd_h = nd_h / SIGMA_OBS_TC
          nd_v = ABS(LOG(obs(obsda2(ip)%set(nobs_use(n)))%lev(iidx)) - LOG(rlev)) / SIGMA_OBSV_TC
!          nd_v = 0.0d0
        case (id_H08IR_obs)                                                                                               ! H08       
          nd_h = nd_h / SIGMA_OBS_H08                                                                                     ! H08                       
          nd_v = ABS(LOG(obs(obsda2(ip)%set(nobs_use(n)))%lev(iidx)) - LOG(rlev)) / SIGMA_OBSV_H08 ! H08 ! bug fixed (02/09/2016) 
        case default
          nd_h = nd_h / SIGMA_OBS
          nd_v = ABS(LOG(obs(obsda2(ip)%set(nobs_use(n)))%lev(iidx)) - LOG(rlev)) / SIGMA_OBSV
        end select

        ndist = nd_h * nd_h + nd_v * nd_v

        if (ndist > dist_zero_fac * dist_zero_fac) cycle
        !
        ! variable localization
        !
        if (nvar > 0) then ! use variable localization only when nvar > 0
          ielm_varlocal = uid_obs_varlocal(ielm)
          if (ielm_varlocal <= 0) then
            write (6,'(A)') 'xxx Warning!!! unsupport observation type in variable localization!'
            ielm_varlocal = 1
          end if
          if (var_local(nvar,ielm_varlocal) < tiny(var_local)) cycle
        end if


        if (MAX_NOBS_PER_GRID <= 0) then

          nobsl = nobsl + 1
          nobsl_t(ielm_u,ityp) = nobsl_t(ielm_u,ityp) + 1
          hdxf(nobsl,:) = obsda2(ip)%ensval(:,nobs_use(n))
          dep(nobsl) = obsda2(ip)%val(nobs_use(n))
          rdiag(nobsl) = obs(obsda2(ip)%set(nobs_use(n)))%err(iidx) * obs(obsda2(ip)%set(nobs_use(n)))%err(iidx)

          ! Observational localization
          rloc(nobsl) = EXP(-0.5d0 * ndist)

          ! variable localization
          if (nvar > 0) then ! use variable localization only when nvar > 0
            rloc(nobsl) = rloc(nobsl) * var_local(nvar,ielm_varlocal)
          end if

          ! radar edge taper (localization)
          if (RADAR_EDGE_TAPER_WIDTH > 0.0d0 .and. RADAR_RANGE > 0.0d0 .and. &                 !!!!!!
              OBS_IN_FORMAT(obsda2(ip)%set(nobs_use(n))) == 2) then                            !!!!!!
            rdx = (ri - radar_ri(obsda2(ip)%set(nobs_use(n)))) * DX                            !!!!!!
            rdy = (rj - radar_rj(obsda2(ip)%set(nobs_use(n)))) * DY                            !!!!!!
            dist_radar_edge = (RADAR_RANGE - sqrt(rdx*rdx + rdy*rdy)) / RADAR_EDGE_TAPER_WIDTH !!!!!!
            if (dist_radar_edge < 0.0d0) then                                                  !!!!!!
              dist_radar_edge = 0.0d0                                                          !!!!!!
!              write (6,'(A,F12.2)') '[Warning] dist_radar_edge < 0. dist_radar_edge =', RADAR_RANGE - sqrt(rdx*rdx + rdy*rdy)
            end if                                                                             !!!!!!
            if (dist_radar_edge < 1.0d0) then                                                  !!!!!!
              rloc(nobsl) = rloc(nobsl) * dist_radar_edge                                      !!!!!!
            end if                                                                             !!!!!!
          end if                                                                               !!!!!!

        else

          do s = 1, MAX_NOBS_PER_GRID
            if (isort_t(s,ielm_u,ityp) == 0) then
              nobsl_t(ielm_u,ityp) = nobsl_t(ielm_u,ityp) + 1
              isort_t(s,ielm_u,ityp) = nobsl_t(ielm_u,ityp)
              ndist_t(isort_t(s,ielm_u,ityp),ielm_u,ityp) = ndist
              ip_t(isort_t(s,ielm_u,ityp),ielm_u,ityp) = ip
              iob_t(isort_t(s,ielm_u,ityp),ielm_u,ityp) = nobs_use(n)
            else if (ndist < ndist_t(isort_t(s,ielm_u,ityp),ielm_u,ityp)) then
              if (nobsl_t(ielm_u,ityp) < MAX_NOBS_PER_GRID) then
                nobsl_t(ielm_u,ityp) = nobsl_t(ielm_u,ityp) + 1
                do ss = nobsl_t(ielm_u,ityp), s+1, -1
                  isort_t(ss,ielm_u,ityp) = isort_t(ss-1,ielm_u,ityp)
                end do
                isort_t(s,ielm_u,ityp) = nobsl_t(ielm_u,ityp)
                ndist_t(isort_t(s,ielm_u,ityp),ielm_u,ityp) = ndist
                ip_t(isort_t(s,ielm_u,ityp),ielm_u,ityp) = ip
                iob_t(isort_t(s,ielm_u,ityp),ielm_u,ityp) = nobs_use(n)
              else
                isort_t(s,ielm_u,ityp) = isort_t(MAX_NOBS_PER_GRID,ielm_u,ityp)
                do ss = MAX_NOBS_PER_GRID, s+1, -1
                  isort_t(ss,ielm_u,ityp) = isort_t(ss-1,ielm_u,ityp)
                end do
                ndist_t(isort_t(s,ielm_u,ityp),ielm_u,ityp) = ndist
                ip_t(isort_t(s,ielm_u,ityp),ielm_u,ityp) = ip
                iob_t(isort_t(s,ielm_u,ityp),ielm_u,ityp) = nobs_use(n)
              end if
            end if
          end do

        end if

      end do ! [ n = 1, nn ]

    end if ! [ obsda2(ip)%nobs > 0 ]

  end do ! [ ip = 0, MEM_NP-1 ]


  if (MAX_NOBS_PER_GRID > 0) then
    do ityp = 1, nobtype
      do ielm_u = 1, nid_obs
        do s = 1, nobsl_t(ielm_u,ityp)

          nobsl = nobsl + 1
          ip = ip_t(s,ielm_u,ityp)
          iob = iob_t(s,ielm_u,ityp)
          iidx = obsda2(ip)%idx(iob)

          hdxf(nobsl,:) = obsda2(ip)%ensval(:,iob)
          dep(nobsl) = obsda2(ip)%val(iob)
          rdiag(nobsl) = obs(obsda2(ip)%set(iob))%err(iidx) * obs(obsda2(ip)%set(iob))%err(iidx)

          ! Observational localization
          rloc(nobsl) = EXP(-0.5d0 * ndist_t(s,ielm_u,ityp))

          ! variable localization
          if (nvar > 0) then ! use variable localization only when nvar > 0
            ielm_varlocal = uid_obs_varlocal(elem_uid(ielm_u))
            if (ielm_varlocal <= 0) then
              write (6,'(A)') 'xxx Warning!!! unsupport observation type in variable localization!'
              ielm_varlocal = 1
            end if
            rloc(nobsl) = rloc(nobsl) * var_local(nvar,ielm_varlocal)
          end if

          ! radar edge taper (localization)
          if (RADAR_EDGE_TAPER_WIDTH > 0.0d0 .and. RADAR_RANGE > 0.0d0 .and. &                 !!!!!!
              OBS_IN_FORMAT(obsda2(ip)%set(iob)) == 2) then                                    !!!!!!
            rdx = (ri - radar_ri(obsda2(ip)%set(iob))) * DX                                    !!!!!!
            rdy = (rj - radar_rj(obsda2(ip)%set(iob))) * DY                                    !!!!!!
            dist_radar_edge = (RADAR_RANGE - sqrt(rdx*rdx + rdy*rdy)) / RADAR_EDGE_TAPER_WIDTH !!!!!!
            if (dist_radar_edge < 0.0d0) then                                                  !!!!!!
              dist_radar_edge = 0.0d0                                                          !!!!!!
!              write (6,'(A,F12.2)') '[Warning] dist_radar_edge < 0. dist_radar_edge =', RADAR_RANGE - sqrt(rdx*rdx + rdy*rdy)
            end if                                                                             !!!!!!
            if (dist_radar_edge < 1.0d0) then                                                  !!!!!!
              rloc(nobsl) = rloc(nobsl) * dist_radar_edge                                      !!!!!!
            end if                                                                             !!!!!!
          end if                                                                               !!!!!!

        end do ! [ s = 1, nobsl_t(ielm_u,ityp) ]
      end do ! [ ielm_u = 1, nid_obs ]
    end do ! [ ityp = 1, nobtype ]
  end if ! [ MAX_NOBS_PER_GRID > 0 ]



  write(6, '(A,I8)') '******', nobsl
  write(6, '(360I3)') nobsl_t(:,:)



!
  IF( nobsl > nobstotal ) THEN
    WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
    WRITE(6,*) 'RI,RJ,LEV,NOBSL,NOBSTOTAL=', ri,rj,rlev,nobsl,nobstotal
    STOP 99
  END IF
!
  IF( allocated(nobs_use) ) DEALLOCATE(nobs_use)

  if (MAX_NOBS_PER_GRID > 0) then
    deallocate(ndist_t)
    deallocate(ip_t)
    deallocate(iob_t)
    deallocate(isort_t)
  end if
!


  ! boundary taper (localization)
  if (BOUNDARY_TAPER_WIDTH > 0.0d0) then                                              !!!!!!
    dist_bdy = min(min(ri - IHALO, nlong+IHALO+1 - ri) * DX, &                        !!!!!!
                   min(rj - JHALO, nlatg+JHALO+1 - rj) * DY) / BOUNDARY_TAPER_WIDTH   !!!!!!
    if (dist_bdy < 0.0d0) then                                                        !!!!!!
      write (6, '(A)') '[Error] ###### Wrong dist_bdy ######'                         !!!!!!
!write (6,*) '$$$$$$####', ri - IHALO, nlong+IHALO+1 - ri, rj - JHALO, nlatg+JHALO+1 - rj
      stop                                                                            !!!!!!
    else if (dist_bdy < 1.0d0) then                                                   !!!!!!
      rloc = rloc * dist_bdy                                                          !!!!!!
    end if                                                                            !!!!!!
  end if                                                                              !!!!!!


  RETURN
END SUBROUTINE obs_local

!-----------------------------------------------------------------------
! Relaxation via LETKF weight - RTPP method
!-----------------------------------------------------------------------
subroutine weight_RTPP(w, wrlx)
  implicit none
  real(r_size), intent(in) :: w(MEMBER,MEMBER)
  real(r_size), intent(out) :: wrlx(MEMBER,MEMBER)
  integer :: m

  wrlx = (1.0d0 - RELAX_ALPHA) * w
  do m = 1, MEMBER
    wrlx(m,m) = wrlx(m,m) + RELAX_ALPHA
  end do

  return
end subroutine weight_RTPP
!-----------------------------------------------------------------------
! Relaxation via LETKF weight - RTPS method
!-----------------------------------------------------------------------
subroutine weight_RTPS(w, pa, xb, wrlx, infl)
  implicit none
  real(r_size), intent(in) :: w(MEMBER,MEMBER)
  real(r_size), intent(in) :: pa(MEMBER,MEMBER)
  real(r_size), intent(in) :: xb(MEMBER)
  real(r_size), intent(out) :: wrlx(MEMBER,MEMBER)
  real(r_size), intent(out) :: infl
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
    infl = RELAX_ALPHA_SPREAD * sqrt(var_g / (var_a * real(MEMBER-1,r_size))) - RELAX_ALPHA_SPREAD + 1.0d0   ! Whitaker and Hamill 2012
!    infl = sqrt(RELAX_ALPHA_SPREAD * (var_g / (var_a * real(MEMBER-1,r_size))) - RELAX_ALPHA_SPREAD + 1.0d0) ! Hamrud et al. 2015 (slightly modified)
    wrlx = w * infl
  else
    wrlx = w
    infl = 1.0d0
  end if

  return
end subroutine weight_RTPS

!SUBROUTINE obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
!  INTEGER,INTENT(INOUT) :: nn, nobs_use(nobs)
!  INTEGER :: j,n,ib,ie,ip

!  DO j=jmin,jmax
!    IF(imin > 1) THEN
!      ib = nobsgrd(imin-1,j)+1
!    ELSE
!      IF(j > 1) THEN
!        ib = nobsgrd(nlon,j-1)+1
!      ELSE
!        ib = 1
!      END IF
!    END IF
!    ie = nobsgrd(imax,j)
!    n = ie - ib + 1
!    IF(n == 0) CYCLE
!    DO ip=ib,ie
!      IF(nn > nobs) THEN
!        WRITE(6,*) 'FATALERROR, NN > NOBS', NN, NOBS
!      END IF
!      nobs_use(nn) = ip
!      nn = nn + 1
!    END DO
!  END DO

!  RETURN
!END SUBROUTINE obs_local_sub

END MODULE letkf_tools
