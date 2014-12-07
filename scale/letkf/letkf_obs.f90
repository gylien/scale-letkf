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
  USE common_mpi
  USE common_scale
  USE common_obs_scale
  USE common_mpi_scale
  USE common_letkf

  use common_nml

!  USE common_precip


!  use common_scalelib


  IMPLICIT NONE
  PUBLIC

  INTEGER,SAVE :: nobs
  LOGICAL,PARAMETER :: omb_output=.TRUE.
  LOGICAL,PARAMETER :: oma_output=.TRUE.
  LOGICAL,PARAMETER :: obsgues_output=.FALSE.
  LOGICAL,PARAMETER :: obsanal_output=.FALSE.
  REAL(r_size),PARAMETER :: sigma_obs=10.0d3
!  REAL(r_size),PARAMETER :: sigma_obs_rain=350.0d3   ! GYL
  REAL(r_size),PARAMETER :: sigma_obsv=0.4d0
!  REAL(r_size),PARAMETER :: sigma_obsv_rain=0.4d0    ! GYL
!  REAL(r_size),PARAMETER :: base_obsv_rain=85000.0d0 ! GYL
  REAL(r_size),PARAMETER :: sigma_obst=3.0d0
  REAL(r_size),SAVE :: dist_zero
!  REAL(r_size),SAVE :: dist_zero_rain
  REAL(r_size),SAVE :: dist_zerov
!  REAL(r_size),SAVE :: dist_zerov_rain
!  REAL(r_size),ALLOCATABLE,SAVE :: obselm(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obslon(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obslat(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obslev(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obserr(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obstyp(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obsdif(:)
!!  REAL(r_size),ALLOCATABLE,SAVE :: obsi(:)
!!  REAL(r_size),ALLOCATABLE,SAVE :: obsj(:)
!!  REAL(r_size),ALLOCATABLE,SAVE :: obsk(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)
!  INTEGER,ALLOCATABLE,SAVE :: obsqc(:) ! GYL: QC values in output diag files (could be any value >= 1)

!  INTEGER,SAVE :: nobsgrd(0:nlon,1:nlat) ! global count
  INTEGER,allocatable,SAVE :: nobsgrd(:,:,:)

  type(obs_info),save :: obs
  type(obs_da_value),save :: obsda
  type(obs_da_value),allocatable,save :: obsda2(:)  ! sorted


!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: nslots=11 ! number of time slots for 4D-LETKF
  INTEGER,PARAMETER :: nbslot=6 ! basetime slot
  REAL(r_size),PARAMETER :: slotint=60.0d0 ! time interval between slots in second

  CHARACTER(7) :: obsfile='obs.dat' !IN
  CHARACTER(21) :: obsdafile='obsda.0000.000000.dat' !IN

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
  use scale_process, only: &
    MPI_COMM_u, &
    PRC_myrank, &
    PRC_NUM_X, &
    PRC_NUM_Y

  IMPLICIT NONE
  REAL(r_size),PARAMETER :: gross_error=5.0d0
  REAL(r_size) :: dlon1,dlon2,dlon,dlat
  REAL(r_size),ALLOCATABLE :: wk2d(:,:)
  INTEGER,ALLOCATABLE :: iwk2d(:,:)
  REAL(r_size),ALLOCATABLE :: tmpelm(:)
  REAL(r_size),ALLOCATABLE :: tmplon(:)
  REAL(r_size),ALLOCATABLE :: tmplat(:)
  REAL(r_size),ALLOCATABLE :: tmplev(:)
  REAL(r_size),ALLOCATABLE :: tmpdat(:)
  REAL(r_size),ALLOCATABLE :: tmperr(:)
  REAL(r_size),ALLOCATABLE :: tmptyp(:)
  REAL(r_size),ALLOCATABLE :: tmpdif(:)
!  REAL(r_size),ALLOCATABLE :: tmpi(:)
!  REAL(r_size),ALLOCATABLE :: tmpj(:)
!  REAL(r_size),ALLOCATABLE :: tmpk(:)
  REAL(r_size),ALLOCATABLE :: tmpdep(:)
  REAL(r_size),ALLOCATABLE :: tmphdxf(:,:)
  INTEGER,ALLOCATABLE :: tmpqc0(:,:)
  INTEGER,ALLOCATABLE :: tmpqc(:)
  REAL(r_size),ALLOCATABLE :: tmp2elm(:)
  REAL(r_size),ALLOCATABLE :: tmp2lon(:)
  REAL(r_size),ALLOCATABLE :: tmp2lat(:)
  REAL(r_size),ALLOCATABLE :: tmp2lev(:)
  REAL(r_size),ALLOCATABLE :: tmp2dat(:)
  REAL(r_size),ALLOCATABLE :: tmp2err(:)
  REAL(r_size),ALLOCATABLE :: tmp2typ(:)
  REAL(r_size),ALLOCATABLE :: tmp2dif(:)
!  REAL(r_size),ALLOCATABLE :: tmp2i(:)
!  REAL(r_size),ALLOCATABLE :: tmp2j(:)
!  REAL(r_size),ALLOCATABLE :: tmp2k(:)
  REAL(r_size),ALLOCATABLE :: tmp2dep(:)
  REAL(r_size),ALLOCATABLE :: tmp2hdxf(:,:)
  INTEGER,ALLOCATABLE :: tmp2qc(:)
  INTEGER :: n,i,j,ierr,nn,l,im
  INTEGER :: nj(0:nlat-1)
  INTEGER :: njs(1:nlat-1)
  CHARACTER(10) :: obsfile='obsNNN.dat'
!  CHARACTER(8) :: cdffile_m='cdfm.grd'         ! GYL, PRECIP assimilation
!  CHARACTER(8) :: cdffile_o='cdfo.grd'         ! GYL
!  CHARACTER(10) :: maskfile='ppmask.grd'       ! GYL
!  REAL(r_size) :: ppcdf_m(nlon,nlat,0:ncdf)    ! GYL
!  REAL(r_size) :: ppcdf_o(nlon,nlat,0:ncdf)    ! GYL
!  REAL(r_size) :: ppzero_m(nlon,nlat)          ! GYL
!  REAL(r_size) :: ppzero_o(nlon,nlat)          ! GYL
!  REAL(r_size) :: ppmask(nlon,nlat)            ! GYL
!  INTEGER :: pp_mem, il, bg_lev, ob_lev        ! GYL
!  INTEGER :: pp_ntotal(pp_bg_nlev,pp_ob_nlev)  ! GYL
!  INTEGER :: zero_mem                          ! GYL
!  REAL(r_size) :: ym, sigma                    ! GYL
!  REAL(r_size) :: ri, rj                       ! GYL
!  REAL(r_size) :: tmpdat_o, obserr_p, obserr_n ! GYL
!  INTEGER :: ii, jj                            ! GYL

  integer :: it,ip
  logical :: check
!  REAL(r_size),allocatable :: bufr(:)
  REAL(r_size),allocatable :: bufr(:,:)
  INTEGER,allocatable :: bufri(:)
  INTEGER,allocatable :: bufri2(:,:,:)


  integer :: iproc,jproc
  integer,allocatable :: ranks(:)

  integer,allocatable :: nnext(:,:)

  integer :: MPI_G_e, MPI_G_obstmp, MPI_COMM_obstmp





!---
  integer :: ns
  integer, allocatable :: nr(:), nrt(:), obsidx(:)
  integer, allocatable :: obsbufs(:), obsbufr(:)
  integer :: ip2, imin1,imax1,jmin1,jmax1,imin2,imax2,jmin2,jmax2
!---





  WRITE(6,'(A)') 'Hello from set_letkf_obs'



  dist_zero = sigma_obs * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zero_rain = sigma_obs_rain * SQRT(10.0d0/3.0d0) * 2.0d0
  dist_zerov = sigma_obsv * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zerov_rain = sigma_obsv_rain * SQRT(10.0d0/3.0d0) * 2.0d0






!  if (valid_member) then
!    call unset_scalelib
!  end if

!!!!!!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!!!!!





  
  if (scale_IO_group_n >= 1) then





!--------------------


    check = .false.
    do it = 1, nitmax
      im = proc2mem(1,it,myrank+1)
      if (im >= 1 .and. im <= nbv) then
        write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is processing member ', &
              im, ', subdomain id #', proc2mem(2,it,myrank+1)

  !      nproc = 0

  !      obsda%nobs = nproc

  !      write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' finishes processing member ', &
  !            im, ', subdomain id #', proc2mem(2,it,myrank+1)
  !      write (6,'(A,I8,A)') ' -- ', nproc, ' observations found'

        write (obsdafile(7:10),'(I4.4)') im
        write (obsdafile(12:17),'(I6.6)') proc2mem(2,it,myrank+1)


        CALL get_nobs(obsdafile,5,obsda%nobs)
        WRITE(6,'(A,I9,A)') 'TOTAL: ', obsda%nobs, ' OBSERVATIONS'

        CALL obs_da_value_allocate(obsda,nbv)


        call read_obs_da(obsdafile,obsda,im,check)
        check = .true.

!if(myrank==2) obsda%qc(4) = 17
!if(myrank==6) obsda%qc(4) = 9


      end if
    end do ! [ it = 1, nitmax ]

!
! if the number of processors is greater then the ensemble size,
! broadcast the observation indices and real grid numbers
! from myrank_e=nbv-1 to the rest of processors that didn't read anything.
!
!####################### now broadcast to all, need to be corrected.
    if (nprocs_e > nbv) then

!      ALLOCATE(ranks(nprocs_e-nbv+1))
!      do n = nbv, nprocs_e
!        ranks(n-nbv+1) = n-1
!      end do
!      call MPI_Comm_group(MPI_COMM_e,MPI_G_e,ierr)
!      call MPI_GROUP_INCL(MPI_G_e,nprocs_e-nbv+1,ranks,MPI_G_obstmp,ierr)
!      call MPI_COMM_CREATE(MPI_COMM_e,MPI_G_obstmp,MPI_COMM_obstmp,ierr)

!      IF(myrank_e+1 >= nbv) THEN
!        CALL MPI_BARRIER(MPI_COMM_obstmp,ierr)
!        call MPI_BCAST(obsda%nobs, 1, MPI_INTEGER, 0, MPI_COMM_obstmp, ierr)
!!    print *, myrank, obsda%nobs

!        if (myrank_e+1 > nbv) then
!          CALL obs_da_value_allocate(obsda,nbv)
!        end if

!        CALL MPI_BARRIER(MPI_COMM_obstmp,ierr)
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

        if (myrank_e+1 > nbv) then
          CALL obs_da_value_allocate(obsda,nbv)
        end if

        CALL MPI_BARRIER(MPI_COMM_e,ierr)
        call MPI_BCAST(obsda%idx, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_e, ierr)
        call MPI_BCAST(obsda%ri, obsda%nobs, MPI_r_size, 0, MPI_COMM_e, ierr)
        call MPI_BCAST(obsda%rj, obsda%nobs, MPI_r_size, 0, MPI_COMM_e, ierr)
!        CALL MPI_BARRIER(MPI_COMM_e,ierr)



    end if


    allocate (bufr(nbv,obsda%nobs))
    bufr = 0.0d0
    CALL MPI_BARRIER(MPI_COMM_e,ierr)
    CALL MPI_ALLREDUCE(obsda%ensval,bufr,obsda%nobs*nbv,MPI_r_size,MPI_SUM,MPI_COMM_e,ierr)
    obsda%ensval = bufr
    deallocate(bufr)

    allocate (bufri(obsda%nobs))
    bufri = 0
    CALL MPI_BARRIER(MPI_COMM_e,ierr)
    CALL MPI_ALLREDUCE(obsda%qc,bufri,obsda%nobs,MPI_INTEGER,MPI_MAX,MPI_COMM_e,ierr)
    obsda%qc = bufri
    deallocate(bufri)

!  write (6,*) obsda%nobs
!  write (6,*) obsda%idx



!    call MPI_Comm_free(MPI_COMM_e,ierr)



!if(myrank==10) then
!    print *, '######======'
!do n = 1, obsda%nobs
!if (maxval(abs(obsda%ensval(:,n))) > 1.0d6) then
!    print *, n, obsda%qc(n)
!    print *, obsda%ensval(:,n)
!    print *, ' '
!end if
!end do
!end if
!stop

!!                                                                               ! GYL, PRECIP assimilation
!! reading precipitation transformation definition and mask                      ! GYL
!!                                                                               ! GYL
!  if (opt_pptrans >= 2) then                                                    ! GYL
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',cdffile_m          ! GYL
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',cdffile_o          ! GYL
!    call read_ppcdf(cdffile_m, cdffile_o, ppcdf_m, ppcdf_o, ppzero_m, ppzero_o) ! GYL
!  end if                                                                        ! GYL
!  WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',maskfile             ! GYL
!  call read_ppmask(maskfile, ppmask)                                            ! GYL
!  pp_ntotal = 0                                                                 ! GYL




!!!###### PRECIP assimilation ######
!!    if (tmpelm(n) == id_rain_obs) then

!!      CALL phys2ij(tmplon(n),tmplat(n),ri,rj)
!!      ii = CEILING(ri-0.5) ! nearest point
!!      jj = CEILING(rj-0.5) ! nearest point
!!      if (ii < 1)    ii = ii + nlon
!!      if (ii > nlon) ii = ii - nlon
!!      if (jj < 1)    jj = 1
!!      if (jj > nlat) jj = nlat

!!      if (ppmask(ii,jj) < mask_thres) then
!!        tmpqc(n) = 0
!!!        write (6,'(A)') '* Precipitation not used because of the mask file'
!!!        write (6,'(A,F6.2,A,F6.2,A)') &
!!!              '*  (lon,lat)=(',tmplon(n),',',tmplat(n),')'
!!!        cycle
!!      end if

!!      pp_mem = 0
!!      do i = 1, nbv
!!        if (tmphdxf(n,i) >= ppzero_thres) then
!!          pp_mem = pp_mem + 1
!!        end if
!!      end do

!!      bg_lev = pp_bg_nlev
!!      do il = 1, pp_bg_nlev-1
!!        if (pp_mem < pp_bg_levs(il)) then
!!          bg_lev = il
!!          exit
!!        end if
!!      end do
!!      ob_lev = pp_ob_nlev
!!      do il = 1, pp_ob_nlev-1
!!        if (tmpdat(n) < pp_ob_levs(il)) then
!!          ob_lev = il
!!          exit
!!        end if
!!      end do
!!      pp_ntotal(bg_lev,ob_lev) = pp_ntotal(bg_lev,ob_lev) + 1
!!      tmpqc(n) = 100 + pp_mem  !! For precip, qc = 100 + number of members with precip

!!      if (.not. pp_criterion(bg_lev,ob_lev)) then
!!        tmpqc(n) = 0
!!!        write (6,'(A)') '* Precipitation does not fit assimilation criterion'
!!!        write (6,'(A,F6.2,A,F6.2,A,I3,A,I2,A,F7.3,A,I2,A)') &
!!!              '*  (lon,lat)=(',tmplon(n),',',tmplat(n),'), pp_mem=', &
!!!              pp_mem, '(', bg_lev, '), pp_obs=', tmpdat(n), '(', ob_lev, ')'
!!        cycle
!!      end if

!!      if (opt_pptrans >= 1) then
!!        tmpdat_o = tmpdat(n)
!!        tmplev(n) = tmpdat_o  !! For precip, lev = original observed value if transformation is used
!!        if (opt_pptrans == 1) then ! log transformation
!!          do i = 1, nbv
!!            tmphdxf(n,i) = pptrans_log(tmphdxf(n,i))
!!          end do
!!          tmpdat(n) = pptrans_log(tmpdat_o)
!!        else if (opt_pptrans == 2) then ! Gaussian transformation with median zero rain
!!          do i = 1, nbv
!!            tmphdxf(n,i) = pptrans_normal(tmphdxf(n,i), ppcdf_m(ii,jj,:), ppzero_m(ii,jj))
!!          end do
!!          tmpdat(n) = pptrans_normal(tmpdat_o, ppcdf_o(ii,jj,:), ppzero_o(ii,jj))
!!        else if (opt_pptrans == 3) then ! Gaussian transformation with modified median zero rain
!!          call pptrans_normal_mdzero_def(tmphdxf(n,:), ppcdf_m(ii,jj,:), ppzero_m(ii,jj), zero_mem, ym, sigma)
!!          tmpdat(n) = pptrans_normal_mdzero(tmpdat_o, ppcdf_o(ii,jj,:), ppzero_o(ii,jj), ppzero_m(ii,jj), zero_mem, ym, sigma)
!!        end if

!!        if (opt_ppobserr == 1) then ! transformed obserr from obs data file
!!          if (opt_pptrans == 1) then ! log transformation
!!            tmperr(n) = tmperr(n) / (tmpdat_o + log_trans_tiny)
!!            if (tmperr(n) < min_ppobserr) tmperr(n) = min_ppobserr
!!          else if (opt_pptrans == 2) then ! Gaussian transformation with median zero rain
!!            obserr_p = pptrans_normal(tmpdat_o+tmperr(n), ppcdf_o(ii,jj,:), ppzero_o(ii,jj)) - tmpdat(n)
!!            if (obserr_p < min_ppobserr) obserr_p = min_ppobserr
!!            obserr_n = tmpdat(n) - pptrans_normal(tmpdat_o-tmperr(n), ppcdf_o(ii,jj,:), ppzero_o(ii,jj))
!!            if (obserr_n < min_ppobserr) obserr_n = min_ppobserr
!!            tmperr(n) = 0.5d0 * (obserr_p + obserr_n)
!!          else if (opt_pptrans == 3) then ! Gaussian transformation with modified median zero rain
!!            obserr_p = pptrans_normal_mdzero(tmpdat_o+tmperr(n), ppcdf_o(ii,jj,:), ppzero_o(ii,jj), ppzero_m(ii,jj), zero_mem, ym, sigma) - tmpdat(n)
!!            if (obserr_p < min_ppobserr) obserr_p = min_ppobserr
!!            obserr_n = tmpdat(n) - pptrans_normal_mdzero(tmpdat_o-tmperr(n), ppcdf_o(ii,jj,:), ppzero_o(ii,jj), ppzero_m(ii,jj), zero_mem, ym, sigma)
!!            if (obserr_n < min_ppobserr) obserr_n = min_ppobserr
!!            tmperr(n) = 0.5d0 * (obserr_p + obserr_n)
!!          end if
!!        else if (opt_ppobserr == 2) then ! constant obserr
!!          tmperr(n) = const_ppobserr
!!        end if
!!      end if ! [ opt_pptrans == 1 ]

!!    end if ! [ tmpelm(n) == id_rain_obs ]
!!!###### end PRECIP assimilation ######

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i)
    do n = 1, obsda%nobs
      IF(obsda%qc(n) > 0) CYCLE
      obsda%val(n) = obsda%ensval(1,n)
      DO i=2,nbv
        obsda%val(n) = obsda%val(n) + obsda%ensval(i,n)
      END DO
      obsda%val(n) = obsda%val(n) / REAL(nbv,r_size)
      DO i=1,nbv
        obsda%ensval(i,n) = obsda%ensval(i,n) - obsda%val(n) ! Hdx
      END DO
      obsda%val(n) = obs%dat(obsda%idx(n)) - obsda%val(n) ! y-Hx
!if (myrank == 0) print *, obsda%idx(n), obs%elm(obsda%idx(n)), obs%dat(obsda%idx(n)), obsda%val(n), obsda%ensval(:,n)
      IF(ABS(obsda%val(n)) > gross_error * obs%err(obsda%idx(n))) THEN !gross error
        obsda%qc(n) = iqc_gross_err
      END IF
    END DO
!!$OMP END PARALLEL DO





!  WRITE(6,'(A)') 'OBSERVATIONAL DEPARTURE STATISTICS:'
!  CALL monit_dep(nobs,tmpelm,tmpdep,tmpqc,1)

!!
!! temporal observation localization
!!
!!  DO n=1,nobs
!!    tmperr(n) = tmperr(n) * exp(0.25d0 * (tmpdif(n) / sigma_obst)**2)
!!  END DO
!!
!! SELECT OBS IN THE NODE
!!

!  nn = 0
!  DO n=1,nobs
!    IF(tmpqc(n) <= 0) CYCLE
!    nn = nn+1
!    tmpelm(nn) = tmpelm(n)
!    tmplon(nn) = tmplon(n)
!    tmplat(nn) = tmplat(n)
!    tmplev(nn) = tmplev(n)
!    tmpdat(nn) = tmpdat(n)
!    tmperr(nn) = tmperr(n)
!    tmptyp(nn) = tmptyp(n)
!    tmpdif(nn) = tmpdif(n)
!!    tmpi(nn) = tmpi(n)
!!    tmpj(nn) = tmpj(n)
!!    tmpk(nn) = tmpk(n)
!    tmpdep(nn) = tmpdep(n)
!    tmphdxf(nn,:) = tmphdxf(n,:)
!    tmpqc(nn) = tmpqc(n)
!  END DO
!  nobs = nn
!  WRITE(6,'(I10,A,I3.3)') nobs,' OBSERVATIONS TO BE ASSIMILATED IN MYRANK ',myrank




!    call set_scalelib


! Sorting

    allocate ( obsda2(0:MEM_NP-1) )

    obsda2(PRC_myrank)%nobs = 0
    allocate ( nobsgrd(0:nlonsub,1:nlatsub,0:MEM_NP-1) )
    nobsgrd = 0
    do n = 1, obsda%nobs
      if (obsda%qc(n) == 0) then
        obsda2(PRC_myrank)%nobs = obsda2(PRC_myrank)%nobs + 1
        call rank_1d_2d(PRC_myrank, iproc, jproc)
        i = ceiling(obsda%ri(n)-0.5) - IHALO - iproc * IMAX
        j = ceiling(obsda%rj(n)-0.5) - JHALO - jproc * JMAX
        nobsgrd(i,j,PRC_myrank) = nobsgrd(i,j,PRC_myrank) + 1
      end if
    end do

    write (6,'(I10,A)') obsda2(PRC_myrank)%nobs,' OBSERVATIONS TO BE ASSIMILATED'

! regional count
    do j = 1, nlatsub
      if (j > 1) then
        nobsgrd(0,j,PRC_myrank) = nobsgrd(nlonsub,j-1,PRC_myrank)
      end if
      do i = 1, nlonsub
        nobsgrd(i,j,PRC_myrank) = nobsgrd(i-1,j,PRC_myrank) + nobsgrd(i,j,PRC_myrank)
      end do
    end do

!do i = 0, MEM_NP-1
!do j = 1, nlatsub
!write (6,'(31I4)') nobsgrd(:,j,i)
!end do
!write (6,*)
!end do

!write (6, *) 'XXXXXX'




!do i = 0, MEM_NP-1
!do j = 1, nlatsub
!write (6,'(31I4)') nobsgrd(:,j,i)
!end do
!write (6,*)
!end do

!! global count
!    do j = 1, nlat
!      if (j > 1) then
!        nobsgrd(0,j) = nobsgrd(nlon,j-1)
!      end if
!      do i = 1, nlon
!        nobsgrd(i,j) = nobsgrd(i-1,j) + nobsgrd(i,j)
!      end do
!    end do

    allocate ( nnext (nlonsub,nlatsub) )
    nnext(1:nlonsub,:) = nobsgrd(0:nlonsub-1,:,PRC_myrank) + 1
    call obs_da_value_allocate(obsda2(PRC_myrank),nbv)
    do n = 1, obsda%nobs
      if (obsda%qc(n) == 0) then
        call rank_1d_2d(PRC_myrank, iproc, jproc)
        i = ceiling(obsda%ri(n)-0.5) - IHALO - iproc * IMAX
        j = ceiling(obsda%rj(n)-0.5) - JHALO - jproc * JMAX

        obsda2(PRC_myrank)%idx(nnext(i,j)) = obsda%idx(n)
        obsda2(PRC_myrank)%val(nnext(i,j)) = obsda%val(n)
        obsda2(PRC_myrank)%ensval(:,nnext(i,j)) = obsda%ensval(:,n)
        obsda2(PRC_myrank)%qc(nnext(i,j)) = obsda%qc(n)
        obsda2(PRC_myrank)%ri(nnext(i,j)) = obsda%ri(n)
        obsda2(PRC_myrank)%rj(nnext(i,j)) = obsda%rj(n)

        nnext(i,j) = nnext(i,j) + 1
      end if
    end do
    deallocate (nnext)

    call obs_da_value_deallocate(obsda)






! communication
    allocate ( bufri2 (0:nlonsub,1:nlatsub,0:MEM_NP-1) )
    call MPI_ALLREDUCE(nobsgrd,bufri2,(nlonsub+1)*nlatsub*MEM_NP,MPI_INTEGER,MPI_SUM,MPI_COMM_u,ierr)
    nobsgrd(0:nlonsub,1:nlatsub,0:MEM_NP-1) = bufri2(0:nlonsub,1:nlatsub,0:MEM_NP-1)
    deallocate ( bufri2 )




    allocate(nr(MEM_NP),nrt(MEM_NP))

!if (PRC_myrank == 0 .and. myrank_e == 0) then
!print *, nobsgrd(nlonsub,nlatsub,:)
!print *, maxval(nobsgrd(nlonsub,nlatsub,:))
!end if

    allocate ( obsbufs(maxval(nobsgrd(nlonsub,nlatsub,:))) )
    allocate ( obsidx(maxval(nobsgrd(nlonsub,nlatsub,:))) )

    do ip = 0, MEM_NP-1
!    do ip = 0, 0

      call rank_1d_2d(ip, iproc, jproc)

!if (PRC_myrank == 0 .and. myrank_e == 0) then
!print *, PRC_NUM_X,IMAX,iproc,jproc,sigma_obs,DX,DY,dist_zero
!end if

      imin1 = max(1, iproc*IMAX+1 - ceiling(dist_zero/DX))
      imax1 = min(PRC_NUM_X*IMAX, (iproc+1)*IMAX + ceiling(dist_zero/DX))
      jmin1 = max(1, jproc*JMAX+1 - ceiling(dist_zero/DY))
      jmax1 = min(PRC_NUM_Y*JMAX, (jproc+1)*JMAX + ceiling(dist_zero/DY))
!      i = ceiling(obsda%ri(n)-0.5) - IHALO - iproc * IMAX
!      j = ceiling(obsda%rj(n)-0.5) - JHALO - jproc * JMAX

!if (myrank_e == 0) then
!print *, PRC_myrank,imin1,imax1,jmin1,jmax1
!end if

      nr = 0
      nrt = 0

      do ip2 = 0, MEM_NP-1

        if (ip2 /= ip) then
          call rank_1d_2d(ip2, iproc, jproc)
!          i = ceiling(obsda%ri(n)-0.5) - IHALO - iproc * IMAX
!          j = ceiling(obsda%rj(n)-0.5) - JHALO - jproc * JMAX

          imin2 = max(1, imin1 - iproc*IMAX)
          imax2 = min(IMAX, imax1 - iproc*IMAX)
          jmin2 = max(1, jmin1 - jproc*JMAX)
          jmax2 = min(JMAX, jmax1 - jproc*JMAX)

!if (myrank_e == 0) then
!print *, PRC_myrank,ip2,imin2,imax2,jmin2,jmax2
!end if


          nr(ip2+1) = 0
          call obs_choose(imin2,imax2,jmin2,jmax2,ip2,nr(ip2+1),obsidx)

          if (ip2 == PRC_myrank) then

!              allocate ( obsidx(obsda2(PRC_myrank)%nobs) )
!              call obs_choose(imin2,imax2,jmin2,jmax2,PRC_myrank,ns,obsidx)

            ns = nr(ip2+1)
            do n = 1, ns
              obsbufs(n) = obsda2(PRC_myrank)%idx(obsidx(n))
            end do

          end if


        end if

        if (ip2 > 0) then
          nrt(ip2+1) = nrt(ip2) + nr(ip2)
        end if

      end do

      allocate ( obsbufr(nrt(MEM_NP)+nr(MEM_NP)) )

      call MPI_GATHERV(obsbufs, ns, MPI_INTEGER, obsbufr, nr, nrt, MPI_INTEGER, ip, MPI_COMM_u, ierr)


      if (PRC_myrank == ip) then

!        write(6,*) ns
        write(6,*) nr(:)
        write(6,*) nrt(:)
        write(6,*)
        write(6,*) nrt(MEM_NP)+nr(MEM_NP)
        write(6,*)
        write(6,*) obsbufr(:)

      end if

      deallocate (obsbufr)

    end do ! [ ip = 0, MEM_MP ]








!    call unset_scalelib


  end if ! [ scale_IO_group_n >= 1 ]




!!!
!!! SORT
!!!
!!  ALLOCATE( tmp2elm(nobs) )
!!  ALLOCATE( tmp2lon(nobs) )
!!  ALLOCATE( tmp2lat(nobs) )
!!  ALLOCATE( tmp2lev(nobs) )
!!  ALLOCATE( tmp2dat(nobs) )
!!  ALLOCATE( tmp2err(nobs) )
!!  ALLOCATE( tmp2typ(nobs) )
!!  ALLOCATE( tmp2dif(nobs) )
!!!  ALLOCATE( tmp2i(nobs) )
!!!  ALLOCATE( tmp2j(nobs) )
!!!  ALLOCATE( tmp2k(nobs) )
!!  ALLOCATE( tmp2dep(nobs) )
!!  ALLOCATE( tmp2hdxf(nobs,nbv) )
!!  ALLOCATE( tmp2qc(nobs) )
!!  ALLOCATE( obselm(nobs) )
!!  ALLOCATE( obslon(nobs) )
!!  ALLOCATE( obslat(nobs) )
!!  ALLOCATE( obslev(nobs) )
!!  ALLOCATE( obsdat(nobs) )
!!  ALLOCATE( obserr(nobs) )
!!  ALLOCATE( obstyp(nobs) )
!!  ALLOCATE( obsdif(nobs) )
!!!  ALLOCATE( obsi(nobs) )
!!!  ALLOCATE( obsj(nobs) )
!!!  ALLOCATE( obsk(nobs) )
!!  ALLOCATE( obsdep(nobs) )
!!  ALLOCATE( obshdxf(nobs,nbv) )
!!  ALLOCATE( obsqc(nobs) )
!!  nobsgrd = 0
!!  nj = 0
!!!$OMP PARALLEL PRIVATE(i,j,n,nn)
!!!$OMP DO SCHEDULE(DYNAMIC)
!!  DO j=1,nlat-1
!!    DO n=1,nobs
!!      IF(tmplat(n) < lat(j) .OR. lat(j+1) <= tmplat(n)) CYCLE
!!      nj(j) = nj(j) + 1
!!    END DO
!!  END DO
!!!$OMP END DO
!!!$OMP DO SCHEDULE(DYNAMIC)
!!  DO j=1,nlat-1
!!    njs(j) = SUM(nj(0:j-1))
!!  END DO
!!!$OMP END DO
!!!$OMP DO SCHEDULE(DYNAMIC)
!!  DO j=1,nlat-1
!!    nn = 0
!!    DO n=1,nobs
!!      IF(tmplat(n) < lat(j) .OR. lat(j+1) <= tmplat(n)) CYCLE
!!      nn = nn + 1
!!      tmp2elm(njs(j)+nn) = tmpelm(n)
!!      tmp2lon(njs(j)+nn) = tmplon(n)
!!      tmp2lat(njs(j)+nn) = tmplat(n)
!!      tmp2lev(njs(j)+nn) = tmplev(n)
!!      tmp2dat(njs(j)+nn) = tmpdat(n)
!!      tmp2err(njs(j)+nn) = tmperr(n)
!!      tmp2typ(njs(j)+nn) = tmptyp(n)
!!      tmp2dif(njs(j)+nn) = tmpdif(n)
!!!      tmp2i(njs(j)+nn) = tmpi(n)
!!!      tmp2j(njs(j)+nn) = tmpj(n)
!!!      tmp2k(njs(j)+nn) = tmpk(n)
!!      tmp2dep(njs(j)+nn) = tmpdep(n)
!!      tmp2hdxf(njs(j)+nn,:) = tmphdxf(n,:)
!!      tmp2qc(njs(j)+nn) = tmpqc(n)
!!    END DO
!!  END DO
!!!$OMP END DO
!!!$OMP DO SCHEDULE(DYNAMIC)
!!  DO j=1,nlat-1
!!    IF(nj(j) == 0) THEN
!!      nobsgrd(:,j) = njs(j)
!!      CYCLE
!!    END IF
!!    nn = 0
!!    DO i=1,nlon
!!      DO n=njs(j)+1,njs(j)+nj(j)
!!        IF(i < nlon) THEN
!!          IF(tmp2lon(n) < lon(i) .OR. lon(i+1) <= tmp2lon(n)) CYCLE
!!        ELSE
!!          IF(tmp2lon(n) < lon(nlon) .OR. 360.0d0 <= tmp2lon(n)) CYCLE
!!        END IF
!!        nn = nn + 1
!!        obselm(njs(j)+nn) = tmp2elm(n)
!!        obslon(njs(j)+nn) = tmp2lon(n)
!!        obslat(njs(j)+nn) = tmp2lat(n)
!!        obslev(njs(j)+nn) = tmp2lev(n)
!!        obsdat(njs(j)+nn) = tmp2dat(n)
!!        obserr(njs(j)+nn) = tmp2err(n)
!!        obstyp(njs(j)+nn) = tmp2typ(n)
!!        obsdif(njs(j)+nn) = tmp2dif(n)
!!!        obsi(njs(j)+nn) = tmp2i(n)
!!!        obsj(njs(j)+nn) = tmp2j(n)
!!!        obsk(njs(j)+nn) = tmp2k(n)
!!        obsdep(njs(j)+nn) = tmp2dep(n)
!!        obshdxf(njs(j)+nn,:) = tmp2hdxf(n,:)
!!        obsqc(njs(j)+nn) = tmp2qc(n)
!!      END DO
!!      nobsgrd(i,j) = njs(j) + nn
!!    END DO
!!    IF(nn /= nj(j)) THEN
!!!$OMP CRITICAL
!!      WRITE(6,'(A,2I)') 'OBS DATA SORT ERROR: ',nn,nj(j)
!!      WRITE(6,'(F6.2,A,F6.2)') lat(j),'< LAT <',lat(j+1)
!!      WRITE(6,'(F6.2,A,F6.2)') MINVAL(tmp2lat(njs(j)+1:njs(j)+nj(j))),'< OBSLAT <',MAXVAL(tmp2lat(njs(j)+1:njs(j)+nj(j)))
!!!$OMP END CRITICAL
!!    END IF
!!  END DO
!!!$OMP END DO
!!!$OMP END PARALLEL
!!  DEALLOCATE( tmp2elm )
!!  DEALLOCATE( tmp2lon )
!!  DEALLOCATE( tmp2lat )
!!  DEALLOCATE( tmp2lev )
!!  DEALLOCATE( tmp2dat )
!!  DEALLOCATE( tmp2err )
!!  DEALLOCATE( tmp2typ )
!!  DEALLOCATE( tmp2dif )
!!!  DEALLOCATE( tmp2i )
!!!  DEALLOCATE( tmp2j )
!!!  DEALLOCATE( tmp2k )
!!  DEALLOCATE( tmp2dep )
!!  DEALLOCATE( tmp2hdxf )
!!  DEALLOCATE( tmp2qc )
!!  DEALLOCATE( tmpelm )
!!  DEALLOCATE( tmplon )
!!  DEALLOCATE( tmplat )
!!  DEALLOCATE( tmplev )
!!  DEALLOCATE( tmpdat )
!!  DEALLOCATE( tmperr )
!!  DEALLOCATE( tmptyp )
!!  DEALLOCATE( tmpdif )
!!!  DEALLOCATE( tmpi )
!!!  DEALLOCATE( tmpj )
!!!  DEALLOCATE( tmpk )
!!  DEALLOCATE( tmpdep )
!!  DEALLOCATE( tmphdxf )
!!  DEALLOCATE( tmpqc )

  RETURN
END SUBROUTINE set_letkf_obs




SUBROUTINE obs_choose(imin,imax,jmin,jmax,proc,nn,nobs_use)

!  use scale_process, only: &
!    PRC_myrank


  implicit none

  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax,proc
  INTEGER,INTENT(INOUT) :: nn
  INTEGER,INTENT(INOUT),OPTIONAL :: nobs_use(:)
  INTEGER :: j,ip

  DO j = jmin, jmax

!if (myrank_e == 0 .and. PRC_myrank == 0) then
!print *, nobsgrd(imin-1,j,proc)+1, nobsgrd(imax,j,proc)
!end if


    DO ip = nobsgrd(imin-1,j,proc)+1, nobsgrd(imax,j,proc)
!      IF(nn > nobs) THEN
!        WRITE(6,*) 'FATALERROR, NN > NOBS', NN, NOBS
!      END IF
      nn = nn + 1
      if (present(nobs_use)) then
        nobs_use(nn) = ip
      end if
    END DO
  END DO

  RETURN
END SUBROUTINE obs_choose


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
!      IF(im > nbv) EXIT
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
!!      IF(im > nbv) EXIT
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
!!-----------------------------------------------------------------------
!! Monitor h(xb) or h(xa) from a LETKF output file
!! Adopted from 'monit_mean' subroutine, 2013/12/24 Guo-Yuan Lien
!!-----------------------------------------------------------------------
!!  file: 'gues' or 'anal'
!!  im:   member # (integer); 0 for ensmean (always called from myrank == 0)
!!  ohx:  h(xb) or h(xa)
!!  oqc:  quality of h(xb) or h(xa)
!!-----------------------------------------------------------------------
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

!  v3d(:,:,:,iv3d_u) = v3dtmp(:,:,:,iv3d_u)
!  v3d(:,:,:,iv3d_v) = v3dtmp(:,:,:,iv3d_v)
!  v3d(:,:,:,iv3d_t) = v3dtmp(:,:,:,iv3d_t)
!  v3d(:,:,:,iv3d_q) = v3dtmp(:,:,:,iv3d_q)
!  v3d(:,:,:,iv3d_qc) = v3dtmp(:,:,:,iv3d_qc)
!  v2d(:,:,iv2d_ps) = v2dtmp(:,:,iv2d_ps)
!  tmpps = reshape(v2d(:,:,iv2d_ps),(/nlon*nlat/))
!  tmptv = reshape(v3d(:,:,:,iv3d_t) * (1.0d0 + fvirt * v3d(:,:,:,iv3d_q)),(/nlon*nlat,nlev/))
!  call sigio_modprd(nlon*nlat,nlon*nlat,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl, &
!                    gfs_vcoord,iret,tmpps,tmptv,pm=tmpp)
!  v3d(:,:,:,iv3d_p) = reshape(tmpp,(/nlon,nlat,nlev/))

!  oqc = 1
!  DO n=1,nobs
!    CALL phys2ijk(v3d(:,:,:,iv3d_p),obselm(n),obslon(n),obslat(n),obslev(n),ri,rj,rk)
!    !
!    ! For monitoring, don't skip any observation below or above model vertical extent.
!    ! Just put bad QC but still get estimate.
!    !
!    IF(CEILING(rk) > nlev) THEN
!      rk = REAL(nlev,r_size)
!      oqc(n) = 0
!    END IF
!    IF(CEILING(rk) < 2 .AND. NINT(obselm(n)) /= id_ps_obs) THEN
!      IF(NINT(obselm(n)) > 9999) THEN
!        rk = 0.0d0
!      ELSE IF(NINT(obselm(n)) == id_u_obs .OR. NINT(obselm(n)) == id_v_obs) THEN
!        rk = 1.00001d0
!      ELSE
!        rk = 1.00001d0
!        oqc(n) = 0
!      END IF
!    END IF
!    IF(NINT(obselm(n)) == id_ps_obs) THEN
!      CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,rk)
!      rk = obslev(n) - rk
!    END IF
!    IF(NINT(obselm(n)) == id_rain_obs) THEN ! No way to get the accumulated precipitation value
!      ohx(n) = obsdat(n)
!      oqc(n) = 0
!    ELSE
!      CALL Trans_XtoY(obselm(n),ri,rj,rk,v3d,v2d,ohx(n))
!    END IF
!  END DO

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

!  dist_zero = sigma_obs * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zero_rain = sigma_obs_rain * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zerov = sigma_obsv * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zerov_rain = sigma_obsv_rain * SQRT(10.0d0/3.0d0) * 2.0d0

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
!  ALLOCATE( obshdxf(nobs,nbv) )
!  ALLOCATE( obsqc(nobs) )
!  ALLOCATE( tmpdep(nobs) )
!  ALLOCATE( tmpqc0(nobs,nbv) )
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
!  CALL read_obs2_mpi(obsanalfile,nobs,nbv,obselm,obslon,obslat,obslev, &
!                     obsdat,obserr,obstyp,obsdif,obshdxf,tmpqc0)

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i)
!  DO n=1,nobs
!    tmpdep(n) = obshdxf(n,1)
!    DO i=2,nbv
!      tmpdep(n) = tmpdep(n) + obshdxf(n,i)
!    END DO
!    tmpdep(n) = tmpdep(n) / REAL(nbv,r_size) ! mean
!    DO i=1,nbv
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
