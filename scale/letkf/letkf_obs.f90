MODULE letkf_obs
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   10/04/2012 Guo-Yuan Lien     modified for GFS model
!   08/30/2013 Guo-Yuan Lien     separating obs operator, following changes by Takemasa MIYOSHI
!   01/01/2014 Guo-Yuan Lien     add EFSO, following changes by Daisuke HOTTA
!
!=======================================================================
!$USE OMP_LIB
  USE common
  use common_nml
  USE common_mpi
  USE common_scale
  USE common_obs_scale
  USE common_mpi_scale
  USE common_letkf

  IMPLICIT NONE
  PUBLIC

  real(r_size),parameter :: dist_zero_fac = 3.651483717        ! SQRT(10.0d0/3.0d0) * 2.0d0
  real(r_size),parameter :: dist_zero_fac_square = 13.33333333 ! dist_zero_fac * dist_zero_fac

  type(obs_info),allocatable,save :: obs(:)
  type(obs_da_value),save :: obsda
  type(obs_da_value),save :: obsda2  ! sorted

  type obs_grid_type
    integer :: ngrd_i
    integer :: ngrd_j
    real(r_size) :: grdspc_i
    real(r_size) :: grdspc_j
    integer :: ngrdsch_i
    integer :: ngrdsch_j
    integer :: ngrdext_i
    integer :: ngrdext_j
    integer, allocatable :: n(:,:,:)
    integer, allocatable :: ac(:,:,:)
    integer, allocatable :: tot(:)
    integer, allocatable :: n_ext(:,:)
    integer, allocatable :: ac_ext(:,:)
    integer :: tot_ext
    integer :: tot_sub(2)             ! only for diagnostic print; 1: before QC; 2: after QC
    integer :: tot_g(2)               ! only for diagnostic print
    integer, allocatable :: next(:,:) ! temporary array
  end type obs_grid_type

  integer,save :: nobs_ext

  type(obs_grid_type),save :: obsgrd(nobtype)

  integer,save :: nobstotalg
  integer,save :: nobstotal
!  integer,save :: maxnobs_sub_per_type
  integer,save :: maxnobs_per_type
CONTAINS
!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------
SUBROUTINE set_letkf_obs
  use scale_grid, only: &
    DX, &
    DY
  use scale_grid_index, only: &
    IHALO,JHALO
!  use scale_process, only: &
!    MPI_COMM_d => LOCAL_COMM_WORLD, &
!    PRC_myrank
  use scale_rm_process, only: &
    PRC_NUM_X, &
    PRC_NUM_Y


  IMPLICIT NONE
  INTEGER :: n,i,j,ierr,im,iof,iidx

  integer :: mem_ref

  integer :: it,ip
  integer :: itype,ielm_u
  real(r_size) :: max_hori_local
  real(r_size) :: target_grdspc

#ifdef H08
  REAL(r_size):: ch_num ! H08
!  REAL(r_size),allocatable :: hx_sprd(:) ! H08
#endif
  integer :: iproc,jproc
  integer :: iproc2,jproc2

  integer :: nobs_sub(2),nobs_g(2)            ! 1: before QC; 2: after QC
  integer :: nobstypevar_g(2,nid_obs,nobtype) !

  character(len=3) :: use_obs_print
  character(4) :: nstr
  character(len=8) :: line1_str(nid_obs) = '========'
  character(len=8) :: line2_str(nid_obs) = '--------'
  


!---
  integer :: ns
  integer :: nr(MEM_NP)
  integer :: nrt(MEM_NP)
  integer, allocatable :: obsidx(:)

  integer :: nk
  integer :: nn_ext, nn_sub
  integer :: ns_ext, ne_ext, ns_bufr, ne_bufr
  integer :: ishift, jshift

  type(obs_da_value) :: obsbufs, obsbufr
  integer :: ip2, imin1,imax1,jmin1,jmax1,imin2,imax2,jmin2,jmax2
!---




  real(r_size),allocatable :: tmpelm(:)
  INTEGER :: monit_nobs(nid_obs)
  REAL(r_size) :: bias(nid_obs)
  REAL(r_size) :: rmse(nid_obs)

  character(filelenmax) :: obsdafile
  character(11) :: obsda_suffix = '.000000.dat'

  type(obs_da_value) :: obsda_ext


  WRITE(6,'(A)') 'Hello from set_letkf_obs'



!-------------------------------------------------------------------------------
! Read externally processed observations
!-------------------------------------------------------------------------------

  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    if (im >= 1 .and. im <= MEMBER) then
      if (it == 1) then
        WRITE(6,'(A,I10)') 'Internally processed observations: ', obsda%nobs - nobs_ext
        WRITE(6,'(A,I10)') 'Externally processed observations: ', nobs_ext
        WRITE(6,'(A,I10)') 'Total                observations: ', obsda%nobs
      end if

      if (OBSDA_IN .and. nobs_ext > 0) then
        obsda_ext%nobs = nobs_ext
        call obs_da_value_allocate(obsda_ext,0)
        write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is reading externally processed observations for member ', &
              im, ', subdomain id #', proc2mem(2,it,myrank+1)
        call file_member_replace(im, OBSDA_IN_BASENAME, obsdafile)
        write (obsda_suffix(2:7),'(I6.6)') proc2mem(2,it,myrank+1)
        call read_obs_da(trim(obsdafile)//obsda_suffix,obsda_ext,0)

        if (OBSDA_OUT) then
          write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is appending observations for member ', &
                im, ', subdomain id #', proc2mem(2,it,myrank+1)
          call file_member_replace(im, OBSDA_OUT_BASENAME, obsdafile)
!          write (obsda_suffix(2:7),'(I6.6)') proc2mem(2,it,myrank+1)
          call write_obs_da(trim(obsdafile)//obsda_suffix,obsda_ext,0,append=.true.)
        end if

        ! variables without an ensemble dimension
        if (it == 1) then
          obsda%set(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%set
          obsda%idx(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%idx
          obsda%ri(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%ri
          obsda%rj(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%rj
          obsda%qc(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%qc
#ifdef H08
          obsda%lev(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%lev
          obsda%val2(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%val2
#endif
        else
          if (maxval(abs(obsda%set(obsda%nobs-nobs_ext+1:obsda%nobs) - obsda_ext%set)) > 0) then
            write (6,'(A)') 'error: obsda%set are inconsistent among the ensemble'
            stop
          end if
          if (maxval(abs(obsda%idx(obsda%nobs-nobs_ext+1:obsda%nobs) - obsda_ext%idx)) > 0) then
            write (6,'(A)') 'error: obsda%idx are inconsistent among the ensemble'
            stop
          end if
          if (maxval(abs(obsda%ri(obsda%nobs-nobs_ext+1:obsda%nobs) - obsda_ext%ri)) > 1.e-6) then
            write (6,'(A)') 'error: obsda%ri are inconsistent among the ensemble'
            stop
          end if
          if (maxval(abs(obsda%rj(obsda%nobs-nobs_ext+1:obsda%nobs) - obsda_ext%rj)) > 1.e-6) then
            write (6,'(A)') 'error: obsda%rj are inconsistent among the ensemble'
            stop
          end if
          obsda%qc(obsda%nobs-nobs_ext+1:obsda%nobs) = max(obsda%qc(obsda%nobs-nobs_ext+1:obsda%nobs), obsda_ext%qc)
#ifdef H08
          obsda%lev(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda%lev(obsda%nobs-nobs_ext+1:obsda%nobs) + obsda_ext%lev
          obsda%val2(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda%val2(obsda%nobs-nobs_ext+1:obsda%nobs) + obsda_ext%val2
#endif
        end if

        ! variables with an ensemble dimension
        obsda%ensval(im,obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%val
      end if ! [ OBSDA_IN .and. nobs_ext > 0 ]
    end if ! [ im >= 1 .and. im <= MEMBER ]
  end do ! [ it = 1, nitmax ]

  ! AllREDUCE observations:
  !   if the number of processors is greater then the ensemble size,
  !   broadcast the observation indices and real grid numbers
  !   from myrank_e=MEMBER-1 to the rest of processors that didn't read anything
  !   !!!!!! now broadcast to all --> need to be revised
  !-----------------------------------------------------------------------------

  if (nprocs_e > MEMBER) then

!      ALLOCATE(ranks(nprocs_e-MEMBER+1))
!      do n = MEMBER, nprocs_e
!        ranks(n-MEMBER+1) = n-1
!      end do
!      call MPI_Comm_group(MPI_COMM_e,MPI_G_e,ierr)
!      call MPI_GROUP_INCL(MPI_G_e,nprocs_e-MEMBER+1,ranks,MPI_G_obstmp,ierr)
!      call MPI_COMM_CREATE(MPI_COMM_e,MPI_G_obstmp,MPI_COMM_obstmp,ierr)

!      IF(myrank_e+1 >= MEMBER) THEN
!        CALL MPI_BARRIER(MPI_COMM_obstmp,ierr)
!        call MPI_BCAST(obsda%nobs, 1, MPI_INTEGER, 0, MPI_COMM_obstmp, ierr)
!!    print *, myrank, obsda%nobs

!        if (myrank_e+1 > MEMBER) then
!          CALL obs_da_value_allocate(obsda,MEMBER)
!        end if

!        CALL MPI_BARRIER(MPI_COMM_obstmp,ierr)
!        call MPI_BCAST(obsda%set, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_obstmp, ierr)
!        call MPI_BCAST(obsda%idx, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_obstmp, ierr)
!        call MPI_BCAST(obsda%ri, obsda%nobs, MPI_r_size, 0, MPI_COMM_obstmp, ierr)
!        call MPI_BCAST(obsda%rj, obsda%nobs, MPI_r_size, 0, MPI_COMM_obstmp, ierr)
!!        CALL MPI_BARRIER(MPI_COMM_obstmp,ierr)
!      end if

!      deallocate(ranks)
!      call MPI_Comm_free(MPI_COMM_obstmp,ierr)


      CALL MPI_BARRIER(MPI_COMM_e,ierr)
      call MPI_BCAST(obsda%nobs, 1, MPI_INTEGER, 0, MPI_COMM_e, ierr)
!    print *, myrank, obsda%nobs

      if (myrank_e+1 > MEMBER) then
        CALL obs_da_value_allocate(obsda,MEMBER)
      end if

      CALL MPI_BARRIER(MPI_COMM_e,ierr)
      call MPI_BCAST(obsda%set, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_e, ierr)
      call MPI_BCAST(obsda%idx, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_e, ierr)
      call MPI_BCAST(obsda%ri, obsda%nobs, MPI_r_size, 0, MPI_COMM_e, ierr)
      call MPI_BCAST(obsda%rj, obsda%nobs, MPI_r_size, 0, MPI_COMM_e, ierr)

  end if ! [ nprocs_e > MEMBER ]


! obsda%val not used; averaged from obsda%ensval later

  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,obsda%ensval,obsda%nobs*MEMBER,MPI_r_size,MPI_SUM,MPI_COMM_e,ierr)

  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,obsda%qc,obsda%nobs,MPI_INTEGER,MPI_MAX,MPI_COMM_e,ierr)

#ifdef H08
!-- H08
! calculate the ensemble mean of obsda%lev
!
  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,obsda%lev,obsda%nobs,MPI_r_size,MPI_SUM,MPI_COMM_e,ierr)

  obsda%lev = obsda%lev / REAL(MEMBER,r_size)

! calculate the ensemble mean of obsda%val2 (clear sky BT)
!
  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,obsda%val2,obsda%nobs,MPI_r_size,MPI_SUM,MPI_COMM_e,ierr)

  obsda%val2 = obsda%val2 / REAL(MEMBER,r_size)

!-- H08
#endif

!  call MPI_Comm_free(MPI_COMM_e,ierr)

!-------------------------------------------------------------------------------
! Process observations and quality control (QC)
!-------------------------------------------------------------------------------

  ! Pre-process radar data
  !-----------------------------------------------------------------------------

  do iof = 1, OBS_IN_NUM
    do n = 1, obs(iof)%nobs
      if (obs(iof)%elm(n) == id_radar_ref_obs) then
        if (obs(iof)%dat(n) >= 0.0d0 .and. obs(iof)%dat(n) < 1.0d10) then
          if (obs(iof)%dat(n) < MIN_RADAR_REF) then
            obs(iof)%dat(n) = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT
          else
            obs(iof)%dat(n) = 10.0d0 * log10(obs(iof)%dat(n))
          end if
        else
          obs(iof)%dat(n) = undef
        end if
        if (USE_OBSERR_RADAR_REF) then
          obs(iof)%err(n) = OBSERR_RADAR_REF
        end if
      end if

      if (USE_OBSERR_RADAR_VR .AND. obs(iof)%elm(n) == id_radar_vr_obs) then
        obs(iof)%err(n) = OBSERR_RADAR_VR
      end if
    end do ! [ n = 1, obs(iof)%nobs ]
  end do ! [ iof = 1, OBS_IN_NUM ]

  ! Compute perturbation and departure
  !  -- gross error check
  !  -- QC based on background (radar reflectivity)
  !  -- process Himawari-8 data
  !-----------------------------------------------------------------------------

  allocate(tmpelm(obsda%nobs))

#ifdef H08
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i,iof,iidx,mem_ref,ch_num)
#else
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i,iof,iidx,mem_ref)
#endif
  do n = 1, obsda%nobs
    IF(obsda%qc(n) > 0) CYCLE

    iof = obsda%set(n)
    iidx = obsda%idx(n)

    tmpelm(n) = obs(iof)%elm(iidx)



!!!###### RADAR assimilation ######
    if (obs(iof)%elm(iidx) == id_radar_ref_obs) then
      if (.not. USE_RADAR_REF) then
        obsda%qc(n) = iqc_otype
        cycle
      end if

      if (obs(iof)%dat(iidx) == undef) then
        obsda%qc(n) = iqc_obs_bad
        cycle
      end if

      !!! obsda%ensval: already converted to dBZ
      mem_ref = 0
      do i = 1, MEMBER
        if (obsda%ensval(i,n) > RADAR_REF_THRES_DBZ+1.0d-6) then
          mem_ref = mem_ref + 1
        end if
      end do
      if (obs(iof)%dat(iidx) > RADAR_REF_THRES_DBZ+1.0d-6) then
        if (mem_ref < MIN_RADAR_REF_MEMBER_OBSREF) then
          obsda%qc(n) = iqc_ref_mem
!          write (6,'(A)') '* Reflectivity does not fit assimilation criterion'
!          write (6,'(A,F6.2,A,F6.2,A,I6,A,F7.3)') &
!                '*  (lon,lat)=(',obs(iof)%lon(iidx),',',obs(iof)%lat(iidx),'), mem_ref=', &
!                mem_ref,', ref_obs=', obs(iof)%dat(iidx)
          cycle
        end if
      else
        if (mem_ref < MIN_RADAR_REF_MEMBER) then
          obsda%qc(n) = iqc_ref_mem
!          write (6,'(A)') '* Reflectivity does not fit assimilation criterion'
!          write (6,'(A,F6.2,A,F6.2,A,I6,A,F7.3)') &
!                '*  (lon,lat)=(',obs(iof)%lon(iidx),',',obs(iof)%lat(iidx),'), mem_ref=', &
!                mem_ref,', ref_obs=', obs(iof)%dat(iidx)
          cycle
        end if
      end if
    end if

    if (obs(iof)%elm(iidx) == id_radar_vr_obs) then
      if (.not. USE_RADAR_VR) then
        obsda%qc(n) = iqc_otype
        cycle
      end if
    end if

!    if (obs(iof)%elm(iidx) == id_radar_prh_obs) then
!      if (.not. USE_RADAR_PSEUDO_RH) then
!        obsda%qc(n) = iqc_otype
!        cycle
!      end if
!    end if
!!!###### end RADAR assimilation ######



#ifdef H08
!!!###### Himawari-8 assimilation ###### ! H08
    if (obs(iof)%elm(iidx) == id_H08IR_obs) then
      if (obs(iof)%dat(iidx) == undef) then
        obsda%qc(n) = iqc_obs_bad
        cycle
      end if

! -- reject Himawari-8 obs sensitivie above H08_LIMIT_LEV (Pa) ! H08 --
      if (obsda%lev(n) < H08_LIMIT_LEV) then
        obsda%qc(n) = iqc_obs_bad
        cycle
      endif

!
! -- Counting how many members have cloud.
! -- Cloudy members should have negative values.
!
      mem_ref = 0
      do i = 1, MEMBER
        if (obsda%ensval(i,n) < 0.0d0) then
          mem_ref = mem_ref + 1
          obsda%ensval(i,n) = obsda%ensval(i,n) * (-1.0d0)
        end if
      end do

!
! -- reject Band #11(ch=5) & #12(ch=6) of Himawari-8 obs ! H08
! -- because these channels are sensitive to chemical tracers
! NOTE!!
!    channel num of Himawari-8 obs is stored in obs%lev (T.Honda 11/04/2015)
!      if ((int(obs(iof)%elm(iidx)) == 11) .or. &
!          (int(obs(iof)%lev(iidx)) == 12)) then
!        obsda%qc(n) = iqc_obs_bad
!        cycle
!      endif
    endif
!!!###### end Himawari-8 assimilation ###### ! H08
#endif



    obsda%val(n) = obsda%ensval(1,n)
    DO i=2,MEMBER
      obsda%val(n) = obsda%val(n) + obsda%ensval(i,n)
    END DO
    obsda%val(n) = obsda%val(n) / REAL(MEMBER,r_size)
#ifdef H08
! Compute CA (cloud effect average, Okamoto et al. 2014QJRMS)
! CA is stored in obsda%val2

    obsda%val2(n) = (abs(obsda%val(n) - obsda%val2(n)) & ! CM
                   + abs(obs(iof)%dat(iidx) - obsda%val2(n)) &! CO
                   &) * 0.5d0
#endif
    DO i=1,MEMBER
      obsda%ensval(i,n) = obsda%ensval(i,n) - obsda%val(n) ! Hdx
    END DO
    obsda%val(n) = obs(iof)%dat(iidx) - obsda%val(n) ! y-Hx

!   compute sprd in obs space ! H08

!    hx_sprd(n) = 0.0d0 !H08
!    DO i=1,MEMBER
!      hx_sprd(n) = hx_sprd(n) + obsda%ensval(i,n) * obsda%ensval(i,n)
!    ENDDO
!    hx_sprd(n) = dsqrt(hx_sprd(n) / REAL(MEMBER,r_size))

    select case (obs(iof)%elm(iidx)) !gross error
    case (id_rain_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RAIN * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_radar_ref_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_REF * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_radar_vr_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_VR * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_radar_prh_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_PRH * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_H08IR_obs)
      ! Adaptive QC depending on the sky condition in the background.
      ! !!Not finished yet!!
      ! 
      ! In config.nml.obsope,
      !  H08_CLDSKY_THRS  < 0.0: turn off ! all members are diagnosed as cloudy.
      !  H08_CLDSKY_THRS  > 0.0: turn on
      !
      IF(mem_ref < H08_MIN_CLD_MEMBER)THEN ! Clear sky
        IF(ABS(obsda%val(n)) > 1.0d0 * obs(iof)%err(iidx)) THEN
          obsda%qc(n) = iqc_gross_err
        END IF
      ELSE ! Cloudy sky
        IF(ABS(obsda%val(n)) > GROSS_ERROR_H08 * obs(iof)%err(iidx)) THEN
          obsda%qc(n) = iqc_gross_err
        END IF
      END IF

      IF(obs(iof)%dat(iidx) < H08_BT_MIN)THEN
        obsda%qc(n) = iqc_gross_err
      ENDIF

!      IF(ABS(obsda%val(n)) > GROSS_ERROR_H08 * obs(iof)%err(iidx)) THEN
!        obsda%qc(n) = iqc_gross_err
!      END IF
    case (id_tclon_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_TCX * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_tclat_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_TCY * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_tcmip_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_TCP * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case default
      IF(ABS(obsda%val(n)) > GROSS_ERROR * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    end select

    IF(obs(iof)%elm(iidx) == id_H08IR_obs)THEN

#ifdef H08
!
! Derived H08 obs height (based on the weighting function output from RTTOV fwd
! model) is substituted into obs%lev.
! Band num. is substituded into obsda%lev. This will be used in monit_obs.
!
      ch_num = obs(iof)%lev(iidx)

      IF(DEPARTURE_STAT_H08)THEN
!
! For obs err correlation statistics based on Desroziers et al. (2005, QJRMS).
!
        write(6, '(a,2I6,2F8.2,4F12.4,2I6,F10.4)')"H08-O-B", &
             obs(iof)%elm(iidx), &
             nint(ch_num), & ! obsda%lev includes band num.
             obs(iof)%lon(iidx), &
             obs(iof)%lat(iidx), &
             obsda%val(n),& ! O-B
             obsda%lev(n), & ! sensitive height
             obs(iof)%dat(iidx), &
             obs(iof)%err(iidx), &
             obsda%qc(n),        &
             mem_ref,  &  ! # of cloudy member
             obsda%val2(n)
      ELSE
        write(6, '(2I6,2F8.2,4F12.4,I3)') &
             obs(iof)%elm(iidx), & ! id
             nint(ch_num), & ! band num
             obs(iof)%lon(iidx), &
             obs(iof)%lat(n), &
             obsda%lev(iidx), & ! sensitive height
             obs(iof)%dat(iidx), &
             obs(iof)%err(iidx), &
             obsda%val(n), &
             obsda%qc(n)
      ENDIF !  [.not. DEPARTURE_STAT_H08]
#endif
    ELSE
      write (6, '(2I6,2F8.2,4F12.4,I3)') obs(iof)%elm(iidx), &
                                         obs(iof)%typ(iidx), &
                                         obs(iof)%lon(iidx), &
                                         obs(iof)%lat(iidx), &
                                         obs(iof)%lev(iidx), &
                                         obs(iof)%dat(iidx), &
                                         obs(iof)%err(iidx), &
                                         obsda%val(n), &
                                         obsda%qc(n)
    ENDIF


  END DO ! [ n = 1, obsda%nobs ]
!$OMP END PARALLEL DO

  ! Temporal observation localization !!!!!! not implemented yet !!!!!!
  !-----------------------------------------------------------------------------
!!
!!  DO n=1,nobs
!!    tmperr(n) = tmperr(n) * exp(0.25d0 * (tmpdif(n) / SIGMA_OBST)**2)
!!  END DO
!!
!! SELECT OBS IN THE NODE
!!

  ! Print departure statistics
  !-------------------------------------------------------------------------------

  write (6, *)  
  write (6,'(A,I6,A)') 'OBSERVATIONAL DEPARTURE STATISTICS (IN THIS SUBDOMAIN #', myrank_d, '):'

  call monit_dep(obsda%nobs, tmpelm, obsda%val, obsda%qc, monit_nobs, bias, rmse)
  call monit_print(monit_nobs, bias, rmse)
  deallocate(tmpelm)

!-------------------------------------------------------------------------------
! "Bucket sort" of observations of each type (with different sorting meshes)
!-------------------------------------------------------------------------------

  ! Determine mesh size for bucket sort
  !-----------------------------------------------------------------------------

  do itype = 1, nobtype
    max_hori_local = HORI_LOCAL(itype)
    if (itype == 22) then  !PHARAD
      max_hori_local = max(max_hori_local, HORI_LOCAL_RADAR_OBSNOREF)
    end if

    if (OBS_SORT_GRID_SPACING(itype) > 0) then
      target_grdspc = OBS_SORT_GRID_SPACING(itype)
    else if (MAX_NOBS_PER_GRID > 0) then
!    if (MAX_NOBS_PER_GRID(itype) > 0) then
      target_grdspc = 0.1d0 * sqrt(real(MAX_NOBS_PER_GRID, r_size)) * OBS_MIN_SPACING(itype)        ! need to be tuned
!      target_grdspc = 0.1d0 * sqrt(real(MAX_NOBS_PER_GRID(itype), r_size)) * OBS_MIN_SPACING(itype) ! need to be tuned
    else
      target_grdspc = max_hori_local * dist_zero_fac / 6.0d0                ! need to be tuned
    end if
    obsgrd(itype)%ngrd_i = min(ceiling(DX * real(nlon,r_size) / target_grdspc), nlon)
    obsgrd(itype)%ngrd_j = min(ceiling(DY * real(nlat,r_size) / target_grdspc), nlat)
    obsgrd(itype)%grdspc_i = DX * real(nlon,r_size) / obsgrd(itype)%ngrd_i
    obsgrd(itype)%grdspc_j = DY * real(nlat,r_size) / obsgrd(itype)%ngrd_j
    obsgrd(itype)%ngrdsch_i = ceiling(max_hori_local * dist_zero_fac / obsgrd(itype)%grdspc_i)
    obsgrd(itype)%ngrdsch_j = ceiling(max_hori_local * dist_zero_fac / obsgrd(itype)%grdspc_j)
    obsgrd(itype)%ngrdext_i = obsgrd(itype)%ngrd_i + obsgrd(itype)%ngrdsch_i * 2
    obsgrd(itype)%ngrdext_j = obsgrd(itype)%ngrd_j + obsgrd(itype)%ngrdsch_j * 2

    allocate (obsgrd(itype)%n (  obsgrd(itype)%ngrd_i, obsgrd(itype)%ngrd_j, 0:MEM_NP-1))
    allocate (obsgrd(itype)%ac(0:obsgrd(itype)%ngrd_i, obsgrd(itype)%ngrd_j, 0:MEM_NP-1))
    allocate (obsgrd(itype)%tot(0:MEM_NP-1))
    allocate (obsgrd(itype)%n_ext (  obsgrd(itype)%ngrdext_i, obsgrd(itype)%ngrdext_j))
    allocate (obsgrd(itype)%ac_ext(0:obsgrd(itype)%ngrdext_i, obsgrd(itype)%ngrdext_j))
    obsgrd(itype)%n (:,:,:) = 0
    obsgrd(itype)%ac(:,:,:) = 0
    obsgrd(itype)%tot(:) = 0
    obsgrd(itype)%n_ext (:,:) = 0
    obsgrd(itype)%ac_ext(:,:) = 0
    obsgrd(itype)%tot_ext = 0
    obsgrd(itype)%tot_sub(:) = 0
    obsgrd(itype)%tot_g(:) = 0

    allocate (obsgrd(itype)%next(obsgrd(itype)%ngrd_i, obsgrd(itype)%ngrd_j))
  end do

  ! Print observation usage settings
  !-----------------------------------------------------------------------------

  write (6, *)
  write (6, '(A)') 'OBSERVATION USAGE SETTINGS:'
  write (6, '(A)') '=============================================================================='
  write (6, '(A)') 'TYPE    USE HORI_LOC   VERT_LOC TIME_LOC MAX_NOBS MIN_SPAC SORT_MESH_X _MESH_Y'
  write (6, '(A)') '                (km) (lnP or m)      (s)              (km)        (km)    (km)'
  write (6, '(A)') '------------------------------------------------------------------------------'
  do itype = 1, nobtype
    if (USE_OBS(itype)) then
      use_obs_print = 'Yes'
    else
      use_obs_print = 'No'
    end if
    select case (itype)
    case (22) ! vertical localization in Z
      write (6, '(A6,1x,A4,F9.2,F7.2,A4,F9.2,I9,F9.2,F12.2,F8.2)') obtypelist(itype), use_obs_print, HORI_LOCAL(itype)/1000.0d0, &
!                VERT_LOCAL(itype)/1000.0d0, '[km]', TIME_LOCAL(itype)/1000.0d0, MAX_NOBS_PER_GRID(itype), &
                VERT_LOCAL(itype)/1000.0d0, '[km]', TIME_LOCAL(itype)/1000.0d0, MAX_NOBS_PER_GRID, &
                OBS_MIN_SPACING(itype)/1000.0d0, obsgrd(itype)%grdspc_i/1000.0d0, obsgrd(itype)%grdspc_j/1000.0d0
    case default ! vertical localization in ln(p)
      write (6, '(A6,1x,A4,F9.2,F11.3,F9.2,I9,F9.2,F12.2,F8.2)') obtypelist(itype), use_obs_print, HORI_LOCAL(itype)/1000.0d0, &
!                VERT_LOCAL(itype), TIME_LOCAL(itype)/1000.0d0, MAX_NOBS_PER_GRID(itype), &
                VERT_LOCAL(itype), TIME_LOCAL(itype)/1000.0d0, MAX_NOBS_PER_GRID, &
                OBS_MIN_SPACING(itype)/1000.0d0, obsgrd(itype)%grdspc_i/1000.0d0, obsgrd(itype)%grdspc_j/1000.0d0
    end select
  end do
  write (6, '(A)') '=============================================================================='

  ! First scan: count the observation numbers in each mesh (in each subdomian)
  !-----------------------------------------------------------------------------

  nobstypevar_g(:,:,:) = 0
  do n = 1, obsda%nobs
    ielm_u = uid_obs(obs(obsda%set(n))%elm(obsda%idx(n)))
    itype = obs(obsda%set(n))%typ(obsda%idx(n))

    if (obsda%qc(n) == iqc_good) then
      call ij_obsgrd(itype, obsda%ri(n), obsda%rj(n), i, j)
      if (i < 1) i = 1                                       ! Assume the process assignment was correct,
      if (i > obsgrd(itype)%ngrd_i) i = obsgrd(itype)%ngrd_i ! so this correction is only to remedy the round-off problem.
      if (j < 1) j = 1                                       !
      if (j > obsgrd(itype)%ngrd_j) j = obsgrd(itype)%ngrd_j !

      obsgrd(itype)%n(i,j,myrank_d) = obsgrd(itype)%n(i,j,myrank_d) + 1
      obsgrd(itype)%tot_sub(2) = obsgrd(itype)%tot_sub(2) + 1             ! only used for diagnostic print (obs number after qc)
      if (ielm_u > 0) then                                                !
        nobstypevar_g(2,ielm_u,itype) = nobstypevar_g(2,ielm_u,itype) + 1 !
      end if                                                              !
    end if

    obsgrd(itype)%tot_sub(1) = obsgrd(itype)%tot_sub(1) + 1               ! only used for diagnostic print (obs number prior to qc)
    if (ielm_u > 0) then                                                  !
      nobstypevar_g(1,ielm_u,itype) = nobstypevar_g(1,ielm_u,itype) + 1   !
    end if                                                                !
  end do

  ! Compute the accumulated numbers in each mesh
  !-----------------------------------------------------------------------------

  do itype = 1, nobtype
    if (itype > 1) then
      obsgrd(itype)%ac(0,1,myrank_d) = obsgrd(itype-1)%ac(obsgrd(itype-1)%ngrd_i,obsgrd(itype-1)%ngrd_j,myrank_d)
    end if
    do j = 1, obsgrd(itype)%ngrd_j
      if (j > 1) then
        obsgrd(itype)%ac(0,j,myrank_d) = obsgrd(itype)%ac(obsgrd(itype)%ngrd_i,j-1,myrank_d)
      end if
      do i = 1, obsgrd(itype)%ngrd_i
        obsgrd(itype)%ac(i,j,myrank_d) = obsgrd(itype)%ac(i-1,j,myrank_d) + obsgrd(itype)%n(i,j,myrank_d)
!        obsgrd(itype)%next(i,j) = obsgrd(itype)%ac(i-1,j,myrank_d)
      end do
    end do
    obsgrd(itype)%next(1:obsgrd(itype)%ngrd_i,:) = obsgrd(itype)%ac(0:obsgrd(itype)%ngrd_i-1,:,myrank_d)
  end do

  ! Second scan: save the indices of bucket-sorted observations in obsda%keys(:)
  !-----------------------------------------------------------------------------

  do n = 1, obsda%nobs
    if (obsda%qc(n) == iqc_good) then
      itype = obs(obsda%set(n))%typ(obsda%idx(n))
      call ij_obsgrd(itype, obsda%ri(n), obsda%rj(n), i, j)
      if (i < 1) i = 1                                       ! Assume the process assignment was correct,
      if (i > obsgrd(itype)%ngrd_i) i = obsgrd(itype)%ngrd_i ! so this correction is only to remedy the round-off problem.
      if (j < 1) j = 1                                       !
      if (j > obsgrd(itype)%ngrd_j) j = obsgrd(itype)%ngrd_j !

      obsgrd(itype)%next(i,j) = obsgrd(itype)%next(i,j) + 1
      obsda%key(obsgrd(itype)%next(i,j)) = n
    end if
  end do


!#ifdef H08
! -- H08
!  if((obs(3)%nobs >= 1) .and. OBS_IN_NUM >= 3)then
!    CALL MPI_BARRIER(MPI_COMM_d,ierr)
!    CALL MPI_ALLREDUCE(MPI_IN_PLACE,obs(3)%lev,obs(3)%nobs,MPI_r_size,MPI_MAX,MPI_COMM_d,ierr)
!  endif
! -- H08
!#endif


  ! ALLREDUCE observation number information from subdomains, and compute total numbers
  !-----------------------------------------------------------------------------

  nobs_sub(:) = 0
  nobs_g(:) = 0
  do itype = 1, nobtype
    call MPI_ALLREDUCE(MPI_IN_PLACE, obsgrd(itype)%n, obsgrd(itype)%ngrd_i*obsgrd(itype)%ngrd_j*MEM_NP, &
                       MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, obsgrd(itype)%ac(0:obsgrd(itype)%ngrd_i,:,:), (obsgrd(itype)%ngrd_i+1)*obsgrd(itype)%ngrd_j*MEM_NP, &
                       MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
    call MPI_ALLREDUCE(obsgrd(itype)%tot_sub, obsgrd(itype)%tot_g, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)

    if (itype == 1) then
      obsgrd(itype)%tot(:) = obsgrd(itype)%ac(obsgrd(itype)%ngrd_i,obsgrd(itype)%ngrd_j,:)
    else
      obsgrd(itype)%tot(:) = obsgrd(itype)%ac(obsgrd(itype)%ngrd_i,obsgrd(itype)%ngrd_j,:) &
                           - obsgrd(itype-1)%ac(obsgrd(itype-1)%ngrd_i,obsgrd(itype-1)%ngrd_j,:)
    end if
    if (obsgrd(itype)%tot(myrank_d) /= obsgrd(itype)%tot_sub(2)) then
      write (6, '(A)') '[Error] Observation counts are inconsistent !!!'
      stop 99
    end if

    nobs_sub(:) = nobs_sub(:) + obsgrd(itype)%tot_sub(:)
    nobs_g(:) = nobs_g(:) + obsgrd(itype)%tot_g(:)

    deallocate (obsgrd(itype)%next)
  end do

  call MPI_ALLREDUCE(MPI_IN_PLACE, nobstypevar_g, 2*nid_obs*nobtype, MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)

  if (obsgrd(nobtype)%ac(obsgrd(nobtype)%ngrd_i,obsgrd(nobtype)%ngrd_j,myrank_d) /= nobs_sub(2)) then
    write (6, '(A)') '[Error] Observation counts are inconsistent !!!'
    stop 99
  end if
  if (sum(obsgrd(nobtype)%ac(obsgrd(nobtype)%ngrd_i,obsgrd(nobtype)%ngrd_j,:)) /= nobs_g(2)) then
    write (6, '(A)') '[Error] Observation counts are inconsistent !!!'
    stop 99
  end if
  nobstotalg = nobs_g(2) ! total obs number in the global domain (all types)

  ! Print observation counts for each types
  !-----------------------------------------------------------------------------

  write(nstr, '(I4)') nid_obs

  write (6, *)
  write (6, '(A)') 'OBSERVATION COUNTS BEFORE QC (GLOABL):'
  write (6, '(A7,'//nstr//'A8,A)') '=======', line1_str(:), '=========='
  write (6, '(A6,1x,'//nstr//'A8,A10)') 'TYPE  ', obelmlist(:), '     TOTAL'
  write (6, '(A7,'//nstr//'A8,A)') '-------', line2_str(:), '----------'
  do itype = 1, nobtype
    write (6, '(A6,1x,'//nstr//'I8,I10)') obtypelist(itype), nobstypevar_g(1,:,itype), obsgrd(itype)%tot_g(1)
  end do
  write (6, '(A7,'//nstr//'A8,A)') '-------', line2_str(:), '----------'
  write (6, '(A6,1x,'//nstr//'I8,I10)') 'TOTAL ', sum(nobstypevar_g(1,:,:), dim=2), nobs_g(1)
  write (6, '(A7,'//nstr//'A8,A)') '=======', line1_str(:), '=========='

  write (6, *)
  write (6, '(A)') 'OBSERVATION COUNTS AFTER QC (GLOABL):'
  write (6, '(A7,'//nstr//'A8,A)') '=======', line1_str(:), '=========='
  write (6, '(A6,1x,'//nstr//'A8,A10)') 'TYPE  ', obelmlist(:), '     TOTAL'
  write (6, '(A7,'//nstr//'A8,A)') '-------', line2_str(:), '----------'
  do itype = 1, nobtype
    write (6, '(A6,1x,'//nstr//'I8,I10)') obtypelist(itype), nobstypevar_g(2,:,itype), obsgrd(itype)%tot_g(2)
  end do
  write (6, '(A7,'//nstr//'A8,A)') '-------', line2_str(:), '----------'
  write (6, '(A6,1x,'//nstr//'I8,I10)') 'TOTAL ', sum(nobstypevar_g(2,:,:), dim=2), nobs_g(2)
  write (6, '(A7,'//nstr//'A8,A)') '=======', line1_str(:), '=========='

  ! Calculate observation numbers in the extended (localization) subdomain,
  ! in preparation for communicating obsetvations in the extended subdomain
  !-----------------------------------------------------------------------------

  do itype = 1, nobtype
    call rank_1d_2d(myrank_d, iproc, jproc)
    imin1 = iproc*obsgrd(itype)%ngrd_i+1 - obsgrd(itype)%ngrdsch_i
    imax1 = (iproc+1)*obsgrd(itype)%ngrd_i + obsgrd(itype)%ngrdsch_i
    jmin1 = jproc*obsgrd(itype)%ngrd_j+1 - obsgrd(itype)%ngrdsch_j
    jmax1 = (jproc+1)*obsgrd(itype)%ngrd_j + obsgrd(itype)%ngrdsch_j

    do ip2 = 0, MEM_NP-1
      call rank_1d_2d(ip2, iproc2, jproc2)
      imin2 = max(1, imin1 - iproc2*obsgrd(itype)%ngrd_i)
      imax2 = min(obsgrd(itype)%ngrd_i, imax1 - iproc2*obsgrd(itype)%ngrd_i)
      jmin2 = max(1, jmin1 - jproc2*obsgrd(itype)%ngrd_j)
      jmax2 = min(obsgrd(itype)%ngrd_j, jmax1 - jproc2*obsgrd(itype)%ngrd_j)

      ishift = (iproc2 - iproc) * obsgrd(itype)%ngrd_i + obsgrd(itype)%ngrdsch_i
      jshift = (jproc2 - jproc) * obsgrd(itype)%ngrd_j + obsgrd(itype)%ngrdsch_j
      obsgrd(itype)%n_ext(imin2+ishift:imax2+ishift, jmin2+jshift:jmax2+jshift) = obsgrd(itype)%n(imin2:imax2, jmin2:jmax2, ip2)
    end do

    if (itype > 1) then
      obsgrd(itype)%ac_ext(0,1) = obsgrd(itype-1)%ac_ext(obsgrd(itype-1)%ngrdext_i,obsgrd(itype-1)%ngrdext_j)
    end if
    do j = 1, obsgrd(itype)%ngrdext_j
      if (j > 1) then
        obsgrd(itype)%ac_ext(0,j) = obsgrd(itype)%ac_ext(obsgrd(itype)%ngrdext_i,j-1)
      end if
      do i = 1, obsgrd(itype)%ngrdext_i
        obsgrd(itype)%ac_ext(i,j) = obsgrd(itype)%ac_ext(i-1,j) + obsgrd(itype)%n_ext(i,j)
      end do
    end do

    if (itype == 1) then
      obsgrd(itype)%tot_ext = obsgrd(itype)%ac_ext(obsgrd(itype)%ngrdext_i,obsgrd(itype)%ngrdext_j)
    else
      obsgrd(itype)%tot_ext = obsgrd(itype)%ac_ext(obsgrd(itype)%ngrdext_i,obsgrd(itype)%ngrdext_j) &
                            - obsgrd(itype-1)%ac_ext(obsgrd(itype-1)%ngrdext_i,obsgrd(itype-1)%ngrdext_j)
    end if
  end do ! [ itype = 1, nobtype ]

  nobstotal = obsgrd(nobtype)%ac_ext(obsgrd(nobtype)%ngrdext_i,obsgrd(nobtype)%ngrdext_j) ! total obs number in the extended subdomain (all types)

!  maxnobs_sub_per_type = maxval(obsgrd(1)%tot(:))
  maxnobs_per_type = obsgrd(1)%tot_ext
  do itype = 2, nobtype
!    maxnobs_sub_per_type = max(maxnobs_sub_per_type, maxval(obsgrd(itype)%tot(:)))
    maxnobs_per_type = max(maxnobs_per_type, obsgrd(itype)%tot_ext)
  end do

  ! Construct sorted obsda2: 
  !-----------------------------------------------------------------------------

  obsda2%nobs = nobstotal
  call obs_da_value_allocate(obsda2, MEMBER)

  ! 1) Copy the observation data in own subdomain to obsda2 with sorted order
  !-----------------------------------------------------------------------------
  nk = 0

  do itype = 1, nobtype
    ishift = obsgrd(itype)%ngrdsch_i
    jshift = obsgrd(itype)%ngrdsch_j

    do j = 1, obsgrd(itype)%ngrd_j
      do n = 1, obsgrd(itype)%ac_ext(obsgrd(itype)%ngrd_i+ishift,j+jshift) &
              - obsgrd(itype)%ac_ext(ishift,j+jshift)
        nn_ext = n + obsgrd(itype)%ac_ext(ishift,j+jshift)
        nn_sub = n + obsgrd(itype)%ac(0,j,myrank_d)
        nk = nk + 1

        obsda2%set(nn_ext) = obsda%set(obsda%key(nn_sub))
        obsda2%idx(nn_ext) = obsda%idx(obsda%key(nn_sub))
        obsda2%key(nk) = nn_ext  ! save the keys of observations within the subdomain (excluding the localization buffer area)
        obsda2%val(nn_ext) = obsda%val(obsda%key(nn_sub))
        obsda2%ensval(:,nn_ext) = obsda%ensval(:,obsda%key(nn_sub))
        obsda2%qc(nn_ext) = obsda%qc(obsda%key(nn_sub))
        obsda2%ri(nn_ext) = obsda%ri(obsda%key(nn_sub))
        obsda2%rj(nn_ext) = obsda%rj(obsda%key(nn_sub))
#ifdef H08
        obsda2%lev(nn_ext) = obsda%lev(obsda%key(nn_sub)) ! H08
#endif
      end do
    end do
  end do ! [ itype = 1, nobtype ]

  if (nk /= nobs_sub(2)) then
    write (6, '(A)') '[Error] Error with number of observations in the subdomain !!!'
    stop 99
  end if

  obsda2%nobs_in_key = nk

  ! 2) Communicate observations within the extended (localization) subdomains
  !-----------------------------------------------------------------------------
  obsbufs%nobs = nobs_sub(2)
  call obs_da_value_allocate(obsbufs, MEMBER)
  allocate (obsidx(nobs_sub(2)))

  do ip = 0, MEM_NP-1

    ! a) Make send buffer with sorted order
    !---------------------------------------------------------------------------
    ns = 0
    nr = 0
    nrt = 0

    do itype = 1, nobtype
      call rank_1d_2d(ip, iproc, jproc)
      imin1 = iproc*obsgrd(itype)%ngrd_i+1 - obsgrd(itype)%ngrdsch_i
      imax1 = (iproc+1)*obsgrd(itype)%ngrd_i + obsgrd(itype)%ngrdsch_i
      jmin1 = jproc*obsgrd(itype)%ngrd_j+1 - obsgrd(itype)%ngrdsch_j
      jmax1 = (jproc+1)*obsgrd(itype)%ngrd_j + obsgrd(itype)%ngrdsch_j

      do ip2 = 0, MEM_NP-1
        if (ip2 /= ip) then
          call rank_1d_2d(ip2, iproc2, jproc2)
          imin2 = max(1, imin1 - iproc2*obsgrd(itype)%ngrd_i)
          imax2 = min(obsgrd(itype)%ngrd_i, imax1 - iproc2*obsgrd(itype)%ngrd_i)
          jmin2 = max(1, jmin1 - jproc2*obsgrd(itype)%ngrd_j)
          jmax2 = min(obsgrd(itype)%ngrd_j, jmax1 - jproc2*obsgrd(itype)%ngrd_j)

          if (myrank_d == ip2) then
            call obs_choose(itype,ip2,imin2,imax2,jmin2,jmax2,nr(ip2+1),obsidx)
          else
            call obs_choose(itype,ip2,imin2,imax2,jmin2,jmax2,nr(ip2+1))
          end if
        end if ! [ ip2 /= ip ]
      end do ! [ ip2 = 0, MEM_NP-1 ]
    end do ! [ itype = 1, nobtype ]

    do ip2 = 1, MEM_NP-1
      nrt(ip2+1) = nrt(ip2) + nr(ip2)
    end do ! [ ip2 = 1, MEM_NP-1 ]

    ns = nr(myrank_d+1)  ! When myrank_d == ip, this should be 0.
    do n = 1, ns
      obsbufs%set(n) = obsda%set(obsda%key(obsidx(n)))
      obsbufs%idx(n) = obsda%idx(obsda%key(obsidx(n)))
      obsbufs%val(n) = obsda%val(obsda%key(obsidx(n)))
      obsbufs%ensval(:,n) = obsda%ensval(:,obsda%key(obsidx(n)))
      obsbufs%qc(n) = obsda%qc(obsda%key(obsidx(n)))
      obsbufs%ri(n) = obsda%ri(obsda%key(obsidx(n)))
      obsbufs%rj(n) = obsda%rj(obsda%key(obsidx(n)))
#ifdef H08
      obsbufs%lev(n) = obsda%lev(obsda%key(obsidx(n))) ! H08
#endif
    end do

    ! b) GATHERV observation data
    !---------------------------------------------------------------------------
    if (nrt(MEM_NP) + nr(MEM_NP) <= 0) cycle

    obsbufr%nobs = nrt(MEM_NP) + nr(MEM_NP)
    call obs_da_value_allocate(obsbufr, MEMBER)

    call MPI_GATHERV(obsbufs%set, ns, MPI_INTEGER, obsbufr%set, nr, nrt, MPI_INTEGER, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%idx, ns, MPI_INTEGER, obsbufr%idx, nr, nrt, MPI_INTEGER, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%val, ns, MPI_r_size, obsbufr%val, nr, nrt, MPI_r_size, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%ensval, ns*MEMBER, MPI_r_size, obsbufr%ensval, nr*MEMBER, nrt*MEMBER, MPI_r_size, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%qc, ns, MPI_INTEGER, obsbufr%qc, nr, nrt, MPI_INTEGER, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%ri, ns, MPI_r_size, obsbufr%ri, nr, nrt, MPI_r_size, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%rj, ns, MPI_r_size, obsbufr%rj, nr, nrt, MPI_r_size, ip, MPI_COMM_d, ierr)
#ifdef H08
    call MPI_GATHERV(obsbufs%lev, ns, MPI_r_size, obsbufr%lev, nr, nrt, MPI_r_size, ip, MPI_COMM_d, ierr) ! H08
#endif

    ! c) In the domain receiving data, copy the receive buffer to obsda2
    !---------------------------------------------------------------------------
    if (myrank_d == ip) then

      ne_bufr = 0

      do ip2 = 0, MEM_NP-1
        if (ip2 /= ip) then

          if (ne_bufr /= nrt(ip2+1)) then
            write (6, '(A)') '[Error] Error in copying receive buffer !!!'
            write (6, *) ip, ip2, ne_bufr, nrt(ip2+1)
            write (6, *) nrt(:)
            stop 99
          end if 

          do itype = 1, nobtype
            call rank_1d_2d(ip, iproc, jproc)
            imin1 = iproc*obsgrd(itype)%ngrd_i+1 - obsgrd(itype)%ngrdsch_i
            imax1 = (iproc+1)*obsgrd(itype)%ngrd_i + obsgrd(itype)%ngrdsch_i
            jmin1 = jproc*obsgrd(itype)%ngrd_j+1 - obsgrd(itype)%ngrdsch_j
            jmax1 = (jproc+1)*obsgrd(itype)%ngrd_j + obsgrd(itype)%ngrdsch_j

            call rank_1d_2d(ip2, iproc2, jproc2)
            imin2 = max(1, imin1 - iproc2*obsgrd(itype)%ngrd_i)
            imax2 = min(obsgrd(itype)%ngrd_i, imax1 - iproc2*obsgrd(itype)%ngrd_i)
            jmin2 = max(1, jmin1 - jproc2*obsgrd(itype)%ngrd_j)
            jmax2 = min(obsgrd(itype)%ngrd_j, jmax1 - jproc2*obsgrd(itype)%ngrd_j)

            ishift = (iproc2 - iproc) * obsgrd(itype)%ngrd_i + obsgrd(itype)%ngrdsch_i
            jshift = (jproc2 - jproc) * obsgrd(itype)%ngrd_j + obsgrd(itype)%ngrdsch_j

            do j = jmin2, jmax2
              ns_ext = obsgrd(itype)%ac_ext(imin2+ishift-1,j+jshift) + 1
              ne_ext = obsgrd(itype)%ac_ext(imax2+ishift  ,j+jshift)

              if (ne_ext - ns_ext + 1 /= obsgrd(itype)%ac(imax2,j,ip2) - obsgrd(itype)%ac(imin2-1,j,ip2)) then
                write (6, '(A)') '[Error] observation grid indices have errors !!!'
                write (6, *) itype, ip, ip2, j, imin1, imax1, jmin1, jmax1, imin2, imax2, jmin2, jmax2, ishift, jshift, &
                             obsgrd(itype)%ngrd_i, obsgrd(itype)%ngrd_j, obsgrd(itype)%ngrdsch_i, obsgrd(itype)%ngrdsch_j, ns_ext, ne_ext, ns_bufr, ne_bufr
                stop 99
              end if

              if (ns_ext > ne_ext) cycle
              ns_bufr = ne_bufr + 1
              ne_bufr = ns_bufr + ne_ext - ns_ext

              obsda2%set(ns_ext:ne_ext) = obsbufr%set(ns_bufr:ne_bufr)
              obsda2%idx(ns_ext:ne_ext) = obsbufr%idx(ns_bufr:ne_bufr)
              obsda2%val(ns_ext:ne_ext) = obsbufr%val(ns_bufr:ne_bufr)
              obsda2%ensval(:,ns_ext:ne_ext) = obsbufr%ensval(:,ns_bufr:ne_bufr)
              obsda2%qc(ns_ext:ne_ext) = obsbufr%qc(ns_bufr:ne_bufr)
              obsda2%ri(ns_ext:ne_ext) = obsbufr%ri(ns_bufr:ne_bufr)
              obsda2%rj(ns_ext:ne_ext) = obsbufr%rj(ns_bufr:ne_bufr)
#ifdef H08
              obsda2%lev(ns_ext:ne_ext) = obsbufr%lev(ns_bufr:ne_bufr) ! H08
#endif
            end do
          end do ! [ itype = 1, nobtype ]

        end if ! [ ip2 /= ip ]
      end do ! [ ip2 = 0, MEM_NP-1 ]

    end if ! [ myrank_d == ip ]

  end do ! [ ip = 0, MEM_NP-1 ]

  do itype = 1, nobtype
    deallocate (obsgrd(itype)%n)
    deallocate (obsgrd(itype)%ac)
    deallocate (obsgrd(itype)%n_ext)
  end do
  call obs_da_value_deallocate(obsda)
  call obs_da_value_deallocate(obsbufs)
  call obs_da_value_deallocate(obsbufr)
  deallocate (obsidx)

  ! Print observation counts
  !-----------------------------------------------------------------------------

  write (6, *)
  write (6, '(A,I6,A)') 'OBSERVATION COUNTS (GLOABL AND IN THIS SUBDOMAIN #', myrank_d, '):'
  write (6, '(A)') '================================================================='
  write (6, '(A)') 'TYPE        GLOBAL  SUBDOMAIN     GLOBAL  SUBDOMAIN EXT_SUBDOMAIN'
  write (6, '(A)') '         before QC  before QC   after QC   after QC      after QC'
  write (6, '(A)') '-----------------------------------------------------------------'
  do itype = 1, nobtype
    write (6, '(A6,1x,4I11,I14)') obtypelist(itype), obsgrd(itype)%tot_g(1), obsgrd(itype)%tot_sub(1), &
              obsgrd(itype)%tot_g(2), obsgrd(itype)%tot_sub(2), obsgrd(itype)%tot_ext
  end do
  write (6, '(A)') '-----------------------------------------------------------------'
    write (6, '(A6,1x,4I11,11x,A3)') 'TOTAL ', nobs_g(1), nobs_sub(1), nobs_g(2), nobs_sub(2), 'N/A'
  write (6, '(A)') '================================================================='

  RETURN
END SUBROUTINE set_letkf_obs

!-----------------------------------------------------------------------
! Convert grid (i,j) values to obsgrid (ogi, ogj) sorting mesh
!-----------------------------------------------------------------------
subroutine ij_obsgrd(obtype, ri, rj, ogi, ogj)
  use scale_grid_index, only: &
    IHALO,JHALO
!  use scale_process, only: &
!    PRC_myrank
  implicit none
  integer, intent(in) :: obtype
  real(r_size), intent(in) :: ri, rj
  integer, intent(out) :: ogi, ogj
  real(r_size) :: ril, rjl

  call rij_g2l(myrank_d, ri, rj, ril, rjl)
  ogi = ceiling((ril - real(IHALO,r_size) - 0.5) * real(obsgrd(obtype)%ngrd_i,r_size) / real(nlon,r_size))
  ogj = ceiling((rjl - real(JHALO,r_size) - 0.5) * real(obsgrd(obtype)%ngrd_i,r_size) / real(nlat,r_size))

  return
end subroutine ij_obsgrd

!-----------------------------------------------------------------------
! Convert grid (i,j) values to obsgrid (ogi, ogj) sorting mesh
! in the extended subdomain
!-----------------------------------------------------------------------
subroutine ij_obsgrd_ext(obtype, ri, rj, ogi, ogj)
  use scale_grid_index, only: &
    IHALO,JHALO
!  use scale_process, only: &
!    PRC_myrank
  implicit none
  integer, intent(in) :: obtype
  real(r_size), intent(in) :: ri, rj
  integer, intent(out) :: ogi, ogj
  real(r_size) :: ril, rjl

  call rij_g2l(myrank_d, ri, rj, ril, rjl)
  ogi = ceiling((ril - real(IHALO,r_size) - 0.5) * real(obsgrd(obtype)%ngrd_i,r_size) / real(nlon,r_size)) &
      + obsgrd(obtype)%ngrdsch_i
  ogj = ceiling((rjl - real(JHALO,r_size) - 0.5) * real(obsgrd(obtype)%ngrd_i,r_size) / real(nlat,r_size)) &
      + obsgrd(obtype)%ngrdsch_j

  return
end subroutine ij_obsgrd_ext

!-----------------------------------------------------------------------
! Choose observations in a rectangle using the bucket sort results
!-----------------------------------------------------------------------
subroutine obs_choose(obtype, proc, imin, imax, jmin, jmax, nn, nobs_use)
  implicit none
  integer, intent(in) :: obtype
  integer, intent(in) :: proc
  integer, intent(in) :: imin, imax, jmin, jmax
  integer, intent(inout) :: nn
  integer, intent(inout), optional :: nobs_use(:)
  integer :: n, j

  if (imin > imax .or. jmin > jmax) return
  if (obsgrd(obtype)%tot(proc) == 0) return

  do j = jmin, jmax
    if (present(nobs_use)) then
      do n = obsgrd(obtype)%ac(imin-1,j,proc)+1, obsgrd(obtype)%ac(imax,j,proc)
        nn = nn + 1
        nobs_use(nn) = n
      end do
    else
      nn = nn + obsgrd(obtype)%ac(imax,j,proc) - obsgrd(obtype)%ac(imin-1,j,proc)
    end if
  end do

  return
end subroutine obs_choose

!-----------------------------------------------------------------------
! Choose observations in a rectangle using the bucket sort results
! in the extended subdomain
!-----------------------------------------------------------------------
subroutine obs_choose_ext(obtype, imin, imax, jmin, jmax, nn, nobs_use)
  implicit none
  integer, intent(in) :: obtype
  integer, intent(in) :: imin, imax, jmin, jmax
  integer, intent(inout) :: nn
  integer, intent(inout), optional :: nobs_use(:)
  integer :: n, j

  if (imin > imax .or. jmin > jmax) return
  if (obsgrd(obtype)%tot_ext == 0) return

  do j = jmin, jmax
    if (present(nobs_use)) then
      do n = obsgrd(obtype)%ac_ext(imin-1,j)+1, obsgrd(obtype)%ac_ext(imax,j)
        nn = nn + 1
        nobs_use(nn) = n
      end do
    else
      nn = nn + obsgrd(obtype)%ac_ext(imax,j) - obsgrd(obtype)%ac_ext(imin-1,j)
    end if
  end do

  return
end subroutine obs_choose_ext


!!-----------------------------------------------------------------------
!! Monitor observation diagnostics after the letkf update
!!  Moved from main program code to this subroutine, 2013/12/26 Guo-Yuan Lien
!!-----------------------------------------------------------------------
!SUBROUTINE monit_obs
!  IMPLICIT NONE
!  REAL(r_size),ALLOCATABLE :: ohx(:)
!  REAL(r_size),ALLOCATABLE :: dep(:)
!  INTEGER,ALLOCATABLE :: oqc(:)
!  INTEGER :: l,im
!  CHARACTER(10) :: ombfile='omb.dat'
!  CHARACTER(10) :: omafile='oma.dat'
!  CHARACTER(14) :: obsguesfile='obsguesNNN.dat'
!  CHARACTER(14) :: obsanalfile='obsanalNNN.dat'
!
!  IF(omb_output .AND. myrank == 0) THEN
!    ALLOCATE(ohx(nobs),oqc(nobs),dep(nobs))
!    CALL monit_output('gues',0,ohx,oqc)
!    dep = obsdat - ohx
!    WRITE(6,'(A)') 'OBSERVATIONAL DEPARTURE STATISTICS [gues mean]:'
!    CALL monit_dep(nobs,obselm,dep,oqc)
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',ombfile
!    CALL write_obs2(ombfile,nobs,obselm,obslon,obslat,obslev, &
!                    obsdat,obserr,obstyp,obsdif,dep,obsqc,0)
!  END IF
!  IF(oma_output .AND. myrank == 0) THEN
!    IF(.NOT. ALLOCATED(ohx)) ALLOCATE(ohx(nobs),oqc(nobs),dep(nobs))
!    CALL monit_output('anal',0,ohx,oqc)
!    dep = obsdat - ohx
!    WRITE(6,'(A)') 'OBSERVATIONAL DEPARTURE STATISTICS [anal mean]:'
!    CALL monit_dep(nobs,obselm,dep,oqc)
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',omafile
!    CALL write_obs2(omafile,nobs,obselm,obslon,obslat,obslev, &
!                    obsdat,obserr,obstyp,obsdif,dep,obsqc,0)
!  END IF

!  IF(obsgues_output) THEN
!    IF(.NOT. ALLOCATED(ohx)) ALLOCATE(ohx(nobs),oqc(nobs))
!    IF(myrank == 0) THEN
!      ohx = obsdat - obsdep
!      oqc = 1
!      WRITE(obsguesfile(8:10),'(A3)') '_me'
!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsguesfile
!      CALL write_obs2(obsguesfile,nobs,obselm,obslon,obslat,obslev, &
!                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!    END IF
!    l=0
!    DO
!      im = myrank+1 + nprocs * l
!      IF(im > MEMBER) EXIT
!      ohx(:) = obsdat(:) - obsdep(:) + obshdxf(:,im)
!      WRITE(obsguesfile(8:10),'(I3.3)') im
!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsguesfile
!      CALL write_obs2(obsguesfile,nobs,obselm,obslon,obslat,obslev, &
!                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!      l = l+1
!    END DO
!  END IF

!! This is not an accurate estimate of obsanal.
!! To obtain an accurate estimate of obsanal, use code in [letkf_tools:das_letkf]
!!------
!!  IF(obsanal_output) THEN
!!    IF(.NOT. ALLOCATED(ohx)) ALLOCATE(ohx(nobs),oqc(nobs))
!!    CALL monit_output('anal',0,ohx,oqc)
!!    WRITE(obsanalfile(8:10),'(A3)') '_me'
!!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsanalfile
!!    CALL write_obs2(obsanalfile,nobs,obselm,obslon,obslat,obslev, &
!!                    obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!!    l=0
!!    DO
!!      im = myrank+1 + nprocs * l
!!      IF(im > MEMBER) EXIT
!!      CALL monit_output('anal',im,ohx,oqc)
!!      WRITE(obsanalfile(8:10),'(I3.3)') im
!!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsanalfile
!!      CALL write_obs2(obsanalfile,nobs,obselm,obslon,obslat,obslev, &
!!                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!!      l = l+1
!!    END DO
!!  END IF

!  IF(ALLOCATED(ohx)) DEALLOCATE(ohx,oqc)
!  IF(ALLOCATED(dep)) DEALLOCATE(dep)

!  RETURN
!END SUBROUTINE monit_obs
!-----------------------------------------------------------------------
! Monitor h(xb) or h(xa) from a LETKF output file
! Adopted from 'monit_mean' subroutine, 2013/12/24 Guo-Yuan Lien
!-----------------------------------------------------------------------
!  file: 'gues' or 'anal'
!  im:   member # (integer); 0 for ensmean (always called from myrank == 0)
!  ohx:  h(xb) or h(xa)
!  oqc:  quality of h(xb) or h(xa)
!-----------------------------------------------------------------------
!SUBROUTINE monit_output(file,im,ohx,oqc)
!  IMPLICIT NONE
!  CHARACTER(4),INTENT(IN) :: file
!  INTEGER,INTENT(IN) :: im
!  REAL(r_size),INTENT(OUT) :: ohx(nobs)
!  INTEGER,INTENT(OUT) :: oqc(nobs)
!  REAL(r_size) :: v3d(nlon,nlat,nlev,nv3dx)
!  REAL(r_size) :: v2d(nlon,nlat,nv2dx)
!  REAL(r_size) :: v3dtmp(nlon,nlat,nlev,nv3d)
!  REAL(r_size) :: v2dtmp(nlon,nlat,nv2d)
!  REAL(r_size) :: elem
!  REAL(r_size) :: bias_u,bias_v,bias_t,bias_ps,bias_q,bias_rh,bias_rain
!  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_ps,rmse_q,rmse_rh,rmse_rain
!  REAL(r_size) :: hdxf,dep,ri,rj,rk
!  INTEGER :: n,iu,iv,it,iq,ips,irh,irain
!  CHARACTER(11) :: filename='filexxx.grd'
!  INTEGER :: iret
!  REAL(r_size) :: tmpps(nlon*nlat)
!  REAL(r_size) :: tmptv(nlon*nlat,nlev)
!  REAL(r_size) :: tmpp(nlon*nlat,nlev)



!  CHARACTER(4),INTENT(IN) :: file
!  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,MEMBER,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nij1,MEMBER,nv2d)
!  REAL(RP) :: v3dg(nlev,nlon,nlat,nv3d)
!  REAL(RP) :: v2dg(nlon,nlat,nv2d)
!  CHARACTER(9) :: filename='file.0000'
!  integer :: it,im,mstart,mend




!  IF(im == 0) THEN
!    WRITE(filename(1:7),'(A4,A3)') file,'_me'
!  ELSE
!    WRITE(filename(1:7),'(A4,I3.3)') file,im
!  END IF
!  CALL read_grd(filename,v3dtmp,v2dtmp,0) ! read ensemble mean into a temporary array
!  IF(im == 0) THEN
!    WRITE(filename(1:7),'(A4,I3.3)') 'gues',1 ! assume called from myrank == 0
!  ELSE
!    WRITE(filename(1:7),'(A4,I3.3)') 'gues',im
!  END IF
!  CALL read_grdx(filename,v3d,v2d) ! only the orography is used, P will be recalulated






!!-----------------------------------------------------------------------




!      WRITE(filename(1:4),'(A4)') file
!      WRITE(filename(6:9),'(I4.4)') im
!!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
!      call read_restart(filename,v3dg,v2dg)






!!  v3d(:,:,:,iv3d_u) = v3dtmp(:,:,:,iv3d_u)
!!  v3d(:,:,:,iv3d_v) = v3dtmp(:,:,:,iv3d_v)
!!  v3d(:,:,:,iv3d_t) = v3dtmp(:,:,:,iv3d_t)
!!  v3d(:,:,:,iv3d_q) = v3dtmp(:,:,:,iv3d_q)
!!  v3d(:,:,:,iv3d_qc) = v3dtmp(:,:,:,iv3d_qc)
!!  v2d(:,:,iv2d_ps) = v2dtmp(:,:,iv2d_ps)
!!  tmpps = reshape(v2d(:,:,iv2d_ps),(/nlon*nlat/))
!!  tmptv = reshape(v3d(:,:,:,iv3d_t) * (1.0d0 + fvirt * v3d(:,:,:,iv3d_q)),(/nlon*nlat,nlev/))
!!  call sigio_modprd(nlon*nlat,nlon*nlat,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl, &
!!                    gfs_vcoord,iret,tmpps,tmptv,pm=tmpp)
!!  v3d(:,:,:,iv3d_p) = reshape(tmpp,(/nlon,nlat,nlev/))

!!  oqc = 1
!!  DO n=1,nobs
!!    CALL phys2ijk(v3d(:,:,:,iv3d_p),obselm(n),obslon(n),obslat(n),obslev(n),ri,rj,rk)
!!    !
!!    ! For monitoring, don't skip any observation below or above model vertical extent.
!!    ! Just put bad QC but still get estimate.
!!    !
!!    IF(CEILING(rk) > nlev) THEN
!!      rk = REAL(nlev,r_size)
!!      oqc(n) = 0
!!    END IF
!!    IF(CEILING(rk) < 2 .AND. NINT(obselm(n)) /= id_ps_obs) THEN
!!      IF(NINT(obselm(n)) > 9999) THEN
!!        rk = 0.0d0
!!      ELSE IF(NINT(obselm(n)) == id_u_obs .OR. NINT(obselm(n)) == id_v_obs) THEN
!!        rk = 1.00001d0
!!      ELSE
!!        rk = 1.00001d0
!!        oqc(n) = 0
!!      END IF
!!    END IF
!!    IF(NINT(obselm(n)) == id_ps_obs) THEN
!!      CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,rk)
!!      rk = obslev(n) - rk
!!    END IF
!!    IF(NINT(obselm(n)) == id_rain_obs) THEN ! No way to get the accumulated precipitation value
!!      ohx(n) = obsdat(n)
!!      oqc(n) = 0
!!    ELSE
!!      CALL Trans_XtoY(obselm(n),ri,rj,rk,v3d,v2d,ohx(n))
!!    END IF
!!  END DO

!  RETURN
!END SUBROUTINE monit_output
!!-----------------------------------------------------------------------
!! Read observation diagnostics for EFSO
!!  Adapted from Y.Ohta's EFSO code for SPEEDY-LETKF   2013/07/17 D.Hotta
!!  Modified, renamed from 'read_monit_obs' to 'set_efso_obs', 2013/12/26 Guo-Yuan Lien
!!-----------------------------------------------------------------------
!SUBROUTINE set_efso_obs
!  IMPLICIT NONE
!  REAL(r_size),ALLOCATABLE :: tmpdep(:)
!  INTEGER,ALLOCATABLE :: tmpqc0(:,:)
!  INTEGER,ALLOCATABLE :: tmpqc(:)
!  INTEGER :: nj(0:nlat-1)
!  INTEGER :: njs(1:nlat-1)
!  INTEGER :: l,im,i,j,n,nn,ierr
!  CHARACTER(14) :: obsguesfile='obsguesNNN.dat'
!  CHARACTER(14) :: obsanalfile='obsanalNNN.dat'

!  WRITE(6,'(A)') 'Hello from set_efso_obs'

!  dist_zero = SIGMA_OBS * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zero_rain = SIGMA_OBS_RAIN * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zerov = SIGMA_OBSV * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zerov_rain = SIGMA_OBSV_RAIN * SQRT(10.0d0/3.0d0) * 2.0d0

!  CALL get_nobs_mpi(obsanalfile,10,nobs)
!  WRITE(6,'(I10,A)') nobs,' TOTAL OBSERVATIONS INPUT'
!  IF(nobs == 0) RETURN
!!
!! INITIALIZE GLOBAL VARIABLES
!!
!  ALLOCATE( obselm(nobs) )
!  ALLOCATE( obslon(nobs) )
!  ALLOCATE( obslat(nobs) )
!  ALLOCATE( obslev(nobs) )
!  ALLOCATE( obsdat(nobs) )
!  ALLOCATE( obserr(nobs) )
!  ALLOCATE( obstyp(nobs) )
!  ALLOCATE( obsdif(nobs) )
!  ALLOCATE( obsdep(nobs) )
!  ALLOCATE( obshdxf(nobs,MEMBER) )
!  ALLOCATE( obsqc(nobs) )
!  ALLOCATE( tmpdep(nobs) )
!  ALLOCATE( tmpqc0(nobs,MEMBER) )
!!
!! reading background observation data and compute departure
!!
!  IF(myrank == 0) THEN
!    WRITE(obsguesfile(8:10),'(A3)') '_me'
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',obsguesfile
!    CALL read_obs2(obsguesfile,nobs,obselm,obslon,obslat,obslev, &
!                   obsdat,obserr,obstyp,obsdif,obsdep,obsqc)
!    obsdep = obsdat - obsdep
!  END IF
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(obsdep,nobs,MPI_r_size,0,MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(obsqc,nobs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!
!! reading ensemble analysis observations
!!
!  CALL read_obs2_mpi(obsanalfile,nobs,MEMBER,obselm,obslon,obslat,obslev, &
!                     obsdat,obserr,obstyp,obsdif,obshdxf,tmpqc0)

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i)
!  DO n=1,nobs
!    tmpdep(n) = obshdxf(n,1)
!    DO i=2,MEMBER
!      tmpdep(n) = tmpdep(n) + obshdxf(n,i)
!    END DO
!    tmpdep(n) = tmpdep(n) / REAL(MEMBER,r_size) ! mean
!    DO i=1,MEMBER
!      obshdxf(n,i) = obshdxf(n,i) - tmpdep(n)
!    END DO
!  END DO
!!$OMP END PARALLEL DO
!  DEALLOCATE(tmpdep,tmpqc0)
!!
!! Create observation box
!!
!  nobsgrd = 0
!  nj = 0
!!$OMP PARALLEL PRIVATE(i,j,n,nn)
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    DO n=1,nobs
!      IF(obslat(n) < lat(j) .OR. lat(j+1) <= obslat(n)) CYCLE
!      nj(j) = nj(j) + 1
!    END DO
!  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    njs(j) = SUM(nj(0:j-1))
!  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    nn = 0
!    DO n=1,nobs
!      IF(obslat(n) < lat(j) .OR. lat(j+1) <= obslat(n)) CYCLE
!      nn = nn + 1
!    END DO
!  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    IF(nj(j) == 0) THEN
!      nobsgrd(:,j) = njs(j)
!      CYCLE
!    END IF
!    nn = 0
!    DO i=1,nlon
!      DO n=njs(j)+1,njs(j)+nj(j)
!        IF(i < nlon) THEN
!          IF(obslon(n) < lon(i) .OR. lon(i+1) <= obslon(n)) CYCLE
!        ELSE
!          IF(obslon(n) < lon(nlon) .OR. 360.0d0 <= obslon(n)) CYCLE
!        END IF
!        nn = nn + 1
!      END DO
!      nobsgrd(i,j) = njs(j) + nn
!    END DO
!    IF(nn /= nj(j)) THEN
!!$OMP CRITICAL
!      WRITE(6,'(A,2I10)') 'OBS DATA SORT ERROR: ',nn,nj(j)
!      WRITE(6,'(F6.2,A,F6.2)') lat(j),'< LAT <',lat(j+1)
!      WRITE(6,'(F6.2,A,F6.2)') MINVAL(obslat(njs(j)+1:njs(j)+nj(j))),'< OBSLAT <',MAXVAL(obslat(njs(j)+1:njs(j)+nj(j)))
!!$OMP END CRITICAL
!    END IF
!  END DO
!!$OMP END DO
!!$OMP END PARALLEL

!  RETURN
!END SUBROUTINE set_efso_obs

END MODULE letkf_obs
