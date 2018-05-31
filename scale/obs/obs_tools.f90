module obs_tools
!=======================================================================
!
! [PURPOSE:]
!
!=======================================================================
!$USE OMP_LIB
  use common
  use common_nml
  use common_scale
  use common_mpi_scale
  use common_obs_scale
  use radar_obs

  implicit none
  public

  integer, private :: obsdep_nobs                     ! obsdep information
  integer, allocatable, private :: obsdep_set(:)      ! 
  integer, allocatable, private :: obsdep_idx(:)      ! 
  integer, allocatable, private :: obsdep_qc(:)       ! 
  real(r_size), allocatable, private :: obsdep_omb(:) ! 
  real(r_size), allocatable, private :: obsdep_oma(:) ! 

contains

!-------------------------------------------------------------------------------
! MPI driver for monitoring observation departure statistics
!-------------------------------------------------------------------------------
subroutine monit_obs_mpi(v3dg, v2dg, monit_step)
  implicit none
  real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP), intent(in) :: v2dg(nlon,nlat,nv2d)
  integer, intent(in) :: monit_step

  integer :: nobs(nid_obs)
  integer :: nobs_g(nid_obs)
  real(r_size) :: bias(nid_obs)
  real(r_size) :: bias_g(nid_obs)
  real(r_size) :: rmse(nid_obs)
  real(r_size) :: rmse_g(nid_obs)
  logical :: monit_type(nid_obs)
  integer :: obsdep_g_nobs
  integer, allocatable :: obsdep_g_set(:)
  integer, allocatable :: obsdep_g_idx(:)
  integer, allocatable :: obsdep_g_qc(:)
  real(r_size), allocatable :: obsdep_g_omb(:)
  real(r_size), allocatable :: obsdep_g_oma(:)
  integer :: cnts
  integer :: cntr(nprocs_d)
  integer :: dspr(nprocs_d)
  integer :: i, ip, ierr

  character(len=14) :: timelabel ! YYYYMMDDHHNNSS

  call date2tlab_LETKF(adate, timelabel)

  call mpi_timer('', 2)

  ! NOTE: need to use 'mmean_rank_e' processes to run this calculation
  !       because only these processes have read topo files in 'topo2d'
  ! 
  if (myrank_e == mmean_rank_e) then
    call monit_obs(v3dg, v2dg, topo2d, nobs, bias, rmse, monit_type, .true., monit_step)

    call mpi_timer('monit_obs_mpi:monit_obs:', 2)

    do i = 1, nid_obs
      if (monit_type(i)) then
        nobs_g(i) = nobs(i)
        if (nobs(i) == 0) then
          bias_g(i) = 0.0d0
          rmse_g(i) = 0.0d0
        else
          bias_g(i) = bias(i) * real(nobs(i), r_size)
          rmse_g(i) = rmse(i) * rmse(i) * real(nobs(i), r_size)
        end if
      end if
    end do

    if (nprocs_d > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, nobs_g, nid_obs, MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, bias_g, nid_obs, MPI_r_size,  MPI_SUM, MPI_COMM_d, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, rmse_g, nid_obs, MPI_r_size,  MPI_SUM, MPI_COMM_d, ierr)
    end if

    do i = 1, nid_obs
      if (monit_type(i)) then
        if (nobs_g(i) == 0) then
          bias_g(i) = undef
          rmse_g(i) = undef
        else
          bias_g(i) = bias_g(i) / REAL(nobs_g(i),r_size)
          rmse_g(i) = sqrt(rmse_g(i) / REAL(nobs_g(i),r_size))
        end if
      else
        nobs_g(i) = -1
        bias_g(i) = undef
        rmse_g(i) = undef
      end if
    end do

    call mpi_timer('monit_obs_mpi:stat:mpi_allreduce(domain):', 2)

    if (OBSDEP_OUT .and. monit_step == 2) then
      cnts = obsdep_nobs
      cntr = 0
      cntr(myrank_d+1) = cnts
      call MPI_ALLREDUCE(MPI_IN_PLACE, cntr, nprocs_d, MPI_INTEGER, MPI_SUM, MPI_COMM_d, ierr)
      dspr = 0
      do ip = 1, nprocs_d-1
        dspr(ip+1) = dspr(ip) + cntr(ip)
      end do

      obsdep_g_nobs = dspr(nprocs_d) + cntr(nprocs_d)
      allocate (obsdep_g_set(obsdep_g_nobs))
      allocate (obsdep_g_idx(obsdep_g_nobs))
      allocate (obsdep_g_qc (obsdep_g_nobs))
      allocate (obsdep_g_omb(obsdep_g_nobs))
      allocate (obsdep_g_oma(obsdep_g_nobs))

      if (obsdep_g_nobs > 0) then
        call MPI_GATHERV(obsdep_set, cnts, MPI_INTEGER, obsdep_g_set, cntr, dspr, MPI_INTEGER, 0, MPI_COMM_d, ierr)
        call MPI_GATHERV(obsdep_idx, cnts, MPI_INTEGER, obsdep_g_idx, cntr, dspr, MPI_INTEGER, 0, MPI_COMM_d, ierr)
        call MPI_GATHERV(obsdep_qc,  cnts, MPI_INTEGER, obsdep_g_qc,  cntr, dspr, MPI_INTEGER, 0, MPI_COMM_d, ierr)
        call MPI_GATHERV(obsdep_omb, cnts, MPI_r_size,  obsdep_g_omb, cntr, dspr, MPI_r_size,  0, MPI_COMM_d, ierr)
        call MPI_GATHERV(obsdep_oma, cnts, MPI_r_size,  obsdep_g_oma, cntr, dspr, MPI_r_size,  0, MPI_COMM_d, ierr)
      end if

      if (myrank_d == 0) then
        write (6,'(A,I6.6,2A)') 'MYRANK ', myrank,' is writing an obsda file ', trim(OBSDEP_OUT_BASENAME)//'.dat'
        call write_obs_dep(trim(OBSDEP_OUT_BASENAME)//'.dat', &
                           obsdep_g_nobs, obsdep_g_set, obsdep_g_idx, obsdep_g_qc, obsdep_g_omb, obsdep_g_oma)
      end if
      deallocate (obsdep_g_set)
      deallocate (obsdep_g_idx)
      deallocate (obsdep_g_qc )
      deallocate (obsdep_g_omb)
      deallocate (obsdep_g_oma)

      call mpi_timer('monit_obs_mpi:obsdep:mpi_allreduce(domain):', 2)
    end if ! [ OBSDEP_OUT .and. monit_step == 2 ]

    if (monit_step == 2) then
      deallocate (obsdep_set)
      deallocate (obsdep_idx)
      deallocate (obsdep_qc )
      deallocate (obsdep_omb)
      deallocate (obsdep_oma)
    end if
  end if ! [ myrank_e == mmean_rank_e ]

  if (DEPARTURE_STAT_ALL_PROCESSES) then
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    call MPI_BCAST(nobs,       nid_obs, MPI_INTEGER, mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(bias,       nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(rmse,       nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(nobs_g,     nid_obs, MPI_INTEGER, mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(bias_g,     nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(rmse_g,     nid_obs, MPI_r_size,  mmean_rank_e, MPI_COMM_e, ierr)
    call MPI_BCAST(monit_type, nid_obs, MPI_LOGICAL, mmean_rank_e, MPI_COMM_e, ierr)

    call mpi_timer('monit_obs_mpi:mpi_allreduce(ens):', 2)
  end if

  if (DEPARTURE_STAT_ALL_PROCESSES .or. myrank_e == mmean_rank_e) then
    if (monit_step == 1) then
      write(6,'(2A)') 'OBSERVATIONAL DEPARTURE STATISTICS [GUESS] (IN THIS SUBDOMAIN): ',timelabel
    else if (monit_step == 2) then
      write(6,'(2A)') 'OBSERVATIONAL DEPARTURE STATISTICS [ANALYSIS] (IN THIS SUBDOMAIN): ',timelabel
    end if
    call monit_print(nobs, bias, rmse, monit_type)

    if (monit_step == 1) then
      write(6,'(2A)') 'OBSERVATIONAL DEPARTURE STATISTICS [GUESS] (GLOBAL): ',timelabel
    else if (monit_step == 2) then
      write(6,'(2A)') 'OBSERVATIONAL DEPARTURE STATISTICS [ANALYSIS] (GLOBAL): ',timelabel
    end if
    call monit_print(nobs_g, bias_g, rmse_g, monit_type)

    call mpi_timer('monit_obs_mpi:monit_print:', 2)
  end if

  return
end subroutine monit_obs_mpi

!-----------------------------------------------------------------------
! Monitor observation departure by giving the v3dg,v2dg data
!-----------------------------------------------------------------------
subroutine monit_obs(v3dg,v2dg,topo,nobs,bias,rmse,monit_type,use_key,step)
  use scale_process, only: &
    PRC_myrank

  implicit none

  REAL(RP),intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),intent(in) :: v2dg(nlon,nlat,nv2d)
  real(r_size),intent(in) :: topo(nlon,nlat)
  INTEGER,INTENT(OUT) :: nobs(nid_obs)
  REAL(r_size),INTENT(OUT) :: bias(nid_obs)
  REAL(r_size),INTENT(OUT) :: rmse(nid_obs)
  LOGICAL,INTENT(OUT) :: monit_type(nid_obs)
  logical,intent(in) :: use_key
  integer,intent(in) :: step

  REAL(r_size) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size) :: v2dgh(nlonh,nlath,nv2dd)
  integer :: nnobs
  integer :: n,nn
  integer :: iset,iidx
  real(r_size) :: ril,rjl,rk,rkz

  real(r_size),allocatable :: oelm(:)
  real(r_size),allocatable :: ohx(:)
  integer,allocatable :: oqc(:)

  integer :: OMP_GET_NUM_THREADS, omp_chunk

!  REAL(r_size) :: timer
!  INTEGER :: ierr

#ifdef H08
! -- for Himawari-8 obs --
  INTEGER :: nprof_H08 ! num of H08 obs
  REAL(r_size),ALLOCATABLE :: ril_H08(:),rjl_H08(:)
  REAL(r_size),ALLOCATABLE :: lon_H08(:),lat_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_ril_H08(:),tmp_rjl_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_lon_H08(:),tmp_lat_H08(:)
  INTEGER,ALLOCATABLE :: n2prof(:) ! obs num 2 prof num

  REAL(r_size),ALLOCATABLE :: yobs_H08(:),plev_obs_H08(:)
  REAL(r_size),ALLOCATABLE :: yobs_H08_clr(:)
  REAL(r_size),ALLOCATABLE :: CA(:) ! (Okamoto et al., 2014QJRMS)
  INTEGER :: ns
  INTEGER,ALLOCATABLE :: qc_H08(:)
#endif

! -- for TC vital assimilation --
!  INTEGER :: obs_idx_TCX, obs_idx_TCY, obs_idx_TCP ! obs index
!  INTEGER :: bTC_proc ! the process where the background TC is located.
! bTC: background TC in each subdomain
! bTC(1,:) : tcx (m), bTC(2,:): tcy (m), bTC(3,:): mslp (Pa)
!  REAL(r_size),ALLOCATABLE :: bTC(:,:)
!  REAL(r_size),ALLOCATABLE :: bufr(:,:)
!  REAL(r_size) :: bTC_mslp

!CALL MPI_BARRIER(MPI_COMM_a,ierr)
!CALL CPU_TIME(timer)
!if (myrank == 0) print *, '######', timer

  call state_to_history(v3dg, v2dg, topo, v3dgh, v2dgh)

  if (use_key) then
    nnobs = obsda_sort%nobs_in_key
  else
    nnobs = obsda_sort%nobs
  end if

  allocate (oelm(nnobs))
  allocate (ohx(nnobs))
  allocate (oqc(nnobs))

#ifdef DEBUG
  if (step < 0 .or. step > 2) then
    write (6, *) '[Error] monit_obs: step should be 0, 1, or 2.'
    stop
  end if
#endif
  if (step == 1) then
    obsdep_nobs = nnobs
    allocate (obsdep_set(obsdep_nobs))
    allocate (obsdep_idx(obsdep_nobs))
    allocate (obsdep_qc (obsdep_nobs))
    allocate (obsdep_omb(obsdep_nobs))
    allocate (obsdep_oma(obsdep_nobs))
  end if

  oqc = -1

!  obs_idx_TCX = -1
!  obs_idx_TCY = -1
!  obs_idx_TCP = -1

!$OMP PARALLEL PRIVATE(n,nn,iset,iidx,ril,rjl,rk,rkz)
  omp_chunk = min(4, max(1, (nnobs-1) / OMP_GET_NUM_THREADS() + 1))
!$OMP DO SCHEDULE(DYNAMIC,omp_chunk)
  do n = 1, nnobs

    if (use_key) then
      nn = obsda_sort%key(n)
    else
      nn = n
    end if

    iset = obsda_sort%set(nn)
    iidx = obsda_sort%idx(nn)

    if (step == 1) then
      obsdep_set(n) = iset
      obsdep_idx(n) = iidx
#ifdef DEBUG
    else if (step == 2) then
      if (obsdep_set(n) /= iset) then
        write (6, *) "[Error] 'set' for y_b and y_a are inconsistent!"
        stop
      end if
      if (obsdep_idx(n) /= iidx) then
        write (6, *) "[Error] 'idx' for y_b and y_a are inconsistent!"
        stop
      end if
#endif
    end if

#ifdef DEBUG
    if (obsda_sort%qc(nn) /= iqc_good) then
      write (6, *) "[Error] The QC value of this observation provided for monitoring is not good: ", &
                   obsda_sort%qc(nn)
      stop
    end if
#endif

    oelm(n) = obs(iset)%elm(iidx)

!    select case (int(oelm(n)))
!    case (id_tclon_obs)
!      obs_idx_TCX = n
!      cycle
!    case (id_tclat_obs)
!      obs_idx_TCY = n
!      cycle
!    case (id_tcmip_obs)
!      obs_idx_TCP = n
!      cycle
!    end select

!!!!!!#ifdef H08
!!!!!!    if(int(oelm(n)) == id_H08IR_obs)cycle
!!!!!!#endif

    call rij_g2l(PRC_myrank, obs(iset)%ri(iidx), obs(iset)%rj(iidx), ril, rjl)
#ifdef DEBUG
    if (PRC_myrank /= obs(iset)%rank(iidx) .or. obs(iset)%rank(iidx) == -1) then
      write (6, *) "[Error] This observation provided for monitoring does not reside in my rank: ", &
                   PRC_myrank, obs(iset)%rank(iidx), obs(iset)%ri(iidx), obs(iset)%rj(iidx), ril, rjl
      stop
    end if
#endif

    if (DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
        abs(obs(iset)%dif(iidx)) <= DEPARTURE_STAT_T_RANGE) then

      oqc(n) = iqc_otype

      select case (OBS_IN_FORMAT(iset))
      !=========================================================================
      case (obsfmt_prepbufr)
      !-------------------------------------------------------------------------
        call phys2ijk(v3dgh(:,:,:,iv3dd_p),obs(iset)%elm(iidx), &
                      ril,rjl,obs(iset)%lev(iidx),rk,oqc(n))
        if (oqc(n) == iqc_good) then
          call Trans_XtoY(obs(iset)%elm(iidx),ril,rjl,rk, &
                          obs(iset)%lon(iidx),obs(iset)%lat(iidx), &
                          v3dgh,v2dgh,ohx(n),oqc(n),stggrd=1)
        end if
      !=========================================================================
      case (obsfmt_radar, obsfmt_pawr_toshiba)
      !-------------------------------------------------------------------------
        if (DEPARTURE_STAT_RADAR) then
          call phys2ijkz(v3dgh(:,:,:,iv3dd_hgt),ril,rjl,obs(iset)%lev(iidx),rkz,oqc(n))
          if (oqc(n) == iqc_good) then
            call Trans_XtoY_radar(obs(iset)%elm(iidx),obs(iset)%meta(1), &
                                  obs(iset)%meta(2),obs(iset)%meta(3),ril,rjl,rkz, &
                                  obs(iset)%lon(iidx),obs(iset)%lat(iidx), &
                                  obs(iset)%lev(iidx),v3dgh,v2dgh,ohx(n),oqc(n),stggrd=1)
            if (oqc(n) == iqc_ref_low) oqc(n) = iqc_good ! when process the observation operator, we don't care if reflectivity is too small
          end if
        end if
#ifdef H08
      !=========================================================================
!      case (obsfmt_h08)
      !-------------------------------------------------------------------------

#endif
      !=========================================================================
      end select

      if (oqc(n) == iqc_good) then
        ohx(n) = obs(iset)%dat(iidx) - ohx(n)
      else
        ohx(n) = undef
      end if

      if (step == 1) then
        obsdep_qc(n) = oqc(n)
        obsdep_omb(n) = ohx(n)
      else if (step == 2) then
        if (obsdep_qc(n) == iqc_good) then ! Use the QC value of y_a only if the QC of y_b is good
          obsdep_qc(n) = oqc(n)            !
        end if                             !
        obsdep_oma(n) = ohx(n)
      end if

      if (LOG_LEVEL >= 3) then
        write (6, '(2I6,2F8.2,4F12.4,I3)') &
              obs(iset)%elm(iidx), &
              obs(iset)%typ(iidx), &
              obs(iset)%lon(iidx), &
              obs(iset)%lat(iidx), &
              obs(iset)%lev(iidx), &
              obs(iset)%dat(iidx), &
              obs(iset)%err(iidx), &
              ohx(n), &
              oqc(n)
      end if

    end if ! [ DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
           !   abs(obs(iset)%dif(iidx)) <= DEPARTURE_STAT_T_RANGE ]

  end do ! [ n = 1, nnobs ]
!$OMP END DO
!$OMP END PARALLEL


#ifdef H08
!
! -- Count the number of the Himawari-8 obs "location" (=num of prof).
!    Then, Trans_XtoY_H08 will be called without openMP.
!

  if (DEPARTURE_STAT_H08) then !-- [DEPARTURE_STAT_H08]

    ALLOCATE(tmp_ril_H08(nnobs))
    ALLOCATE(tmp_rjl_H08(nnobs))
    ALLOCATE(tmp_lon_H08(nnobs))
    ALLOCATE(tmp_lat_H08(nnobs))
    ALLOCATE(n2prof(nnobs))

    n2prof = 0
    nprof_H08 = 0
    do n = 1, nnobs
      if (use_key) then
        nn = obsda_sort%key(n)
      else
        nn = n
      end if
      iset = obsda_sort%set(nn)
      iidx = obsda_sort%idx(nn)

      if (step == 1) then
        obsdep_set(n) = iset
        obsdep_idx(n) = iidx
      end if

      oelm(n) = obs(iset)%elm(iidx)
      if(oelm(n) /= id_H08IR_obs)cycle

      call rij_g2l(PRC_myrank, obs(iset)%ri(iidx), obs(iset)%rj(iidx), ril, rjl)
#ifdef DEBUG
      if (PRC_myrank /= obs(iset)%rank(iidx) .or. obs(iset)%rank(iidx) == -1) then
        write (6, *) "[Error] This observation provided for monitoring does not reside in my rank: ", &
                     PRC_myrank, obs(iset)%rank(iidx), obs(iset)%ri(iidx), obs(iset)%rj(iidx), ril, rjl
        stop
      end if
#endif

      if(nprof_H08 > 1)then
        if((tmp_ril_H08(nprof_H08)==ril) .and. (tmp_rjl_H08(nprof_H08)==rjl))then
          n2prof(n) = nprof_H08
          cycle
        else
          nprof_H08 = nprof_H08 + 1
          tmp_ril_H08(nprof_H08) = ril
          tmp_rjl_H08(nprof_H08) = rjl
          tmp_lon_H08(nprof_H08) = obs(iset)%lon(iidx)
          tmp_lat_H08(nprof_H08) = obs(iset)%lat(iidx)
          n2prof(n) = nprof_H08
        endif
      else ! nprof_H08 <= 1
        nprof_H08 = nprof_H08 + 1
        tmp_ril_H08(nprof_H08) = ril
        tmp_rjl_H08(nprof_H08) = rjl
        tmp_lon_H08(nprof_H08) = obs(iset)%lon(iidx)
        tmp_lat_H08(nprof_H08) = obs(iset)%lat(iidx)
        n2prof(n) = nprof_H08
      endif
    end do ! [ n = 1, nnobs ]

    IF(nprof_H08 >=1)THEN ! [nprof_H08 >=1]
      ALLOCATE(ril_H08(nprof_H08))
      ALLOCATE(rjl_H08(nprof_H08))
      ALLOCATE(lon_H08(nprof_H08))
      ALLOCATE(lat_H08(nprof_H08))

      ril_H08 = tmp_ril_H08(1:nprof_H08)
      rjl_H08 = tmp_rjl_H08(1:nprof_H08)
      lon_H08 = tmp_lon_H08(1:nprof_H08)
      lat_H08 = tmp_lat_H08(1:nprof_H08)

      DEALLOCATE(tmp_ril_H08, tmp_rjl_H08)
      DEALLOCATE(tmp_lon_H08, tmp_lat_H08)

      ALLOCATE(yobs_H08(nprof_H08*nch))
      ALLOCATE(yobs_H08_clr(nprof_H08*nch))
      ALLOCATE(CA(nprof_H08*nch))
      ALLOCATE(plev_obs_H08(nprof_H08*nch))
      ALLOCATE(qc_H08(nprof_H08*nch))


      CALL Trans_XtoY_H08(nprof_H08,ril_H08,rjl_H08,&
                          lon_H08,lat_H08,v3dgh,v2dgh,&
                          yobs_H08,plev_obs_H08,&
                          qc_H08,stggrd=1,yobs_H08_clr=yobs_H08_clr)

      ! yobs here should be positive!!
      yobs_H08 = abs(yobs_H08)

      write (6, '(A)')"MEAN-HIMAWARI-8-STATISTICS"

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,5) PRIVATE(n,nn,iset,iidx,ns)
      do n = 1, nnobs
        if (use_key) then
          nn = obsda_sort%key(n)
        else
          nn = n
        end if
        iset = obsda_sort%set(nn)
        iidx = obsda_sort%idx(nn)

        oelm(n) = obs(iset)%elm(iidx)
        if(oelm(n) /= id_H08IR_obs)cycle
  
        ns = (n2prof(n) - 1) * nch + nint(obsda_sort%lev(nn) - 6.0) 

        if (DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. & 
          abs(obs(iset)%dif(iidx)) <= DEPARTURE_STAT_T_RANGE) then
!          oqc(n) = iqc_otype

          oqc(n) = qc_H08(ns)
          ohx(n) = yobs_H08(ns)
!          if(plev_obs_H08(ns) < H08_LIMIT_LEV) oqc(n) = iqc_obs_bad

          CA(n) =  (abs(yobs_H08(ns) - yobs_H08_clr(ns)) & ! CM
                   +  abs(obs(iset)%dat(iidx) - yobs_H08_clr(ns)) & ! CO
                   &) * 0.5d0 

          if (step == 1) then
            obsdep_qc(n) = oqc(n)
            obsdep_omb(n) = ohx(n)
          else if (step == 2) then
            if (obsdep_qc(n) == iqc_good) then ! Use the QC value of y_a only if the QC of y_b is good
              obsdep_qc(n) = oqc(n)            !
            end if                             !
            obsdep_oma(n) = ohx(n)
          end if
                   
          if(oqc(n) == iqc_good) then
            ohx(n) = obs(iset)%dat(iidx) - ohx(n) 
            write (6, '(A,2I6,2F8.2,5F11.4,I6,F10.4)')"H08-O-A-B",&
                  obs(iset)%elm(iidx), &
                  nint(obsda_sort%lev(nn)), & ! obsda_sort%lev includes the band num.
                  obs(iset)%lon(iidx), &
                  obs(iset)%lat(iidx), &
                  ohx(n), &! O-A
                  obsda_sort%val(nn), &! O-B
                  plev_obs_H08(ns), &
                  obs(iset)%dat(iidx), &
                  obs(iset)%err(iidx), &
                  oqc(n),&
                  CA(n) 
          endif

        endif ! [DEPARTURE_STAT_T_RANGE]
      end do ! [ n = 1, nnobs ]
!$OMP END PARALLEL DO

      DEALLOCATE(yobs_H08, yobs_H08_clr, plev_obs_H08, qc_H08, CA)
    ENDIF ! [nprof_H08 >=1]

    IF(ALLOCATED(n2prof)) DEALLOCATE(n2prof)
    IF(ALLOCATED(tmp_ril_H08)) DEALLOCATE(tmp_ril_H08)
    IF(ALLOCATED(tmp_rjl_H08)) DEALLOCATE(tmp_rjl_H08)
    IF(ALLOCATED(tmp_lon_H08)) DEALLOCATE(tmp_lon_H08)
    IF(ALLOCATED(tmp_lat_H08)) DEALLOCATE(tmp_lat_H08)
  endif !-- [DEPARTURE_STAT_H08]

#endif

! ###  -- TC vital assimilation -- ###
!  if (obs_idx_TCX > 0 .and. obs_idx_TCY > 0 .and. obs_idx_TCP > 0 .and.&
!    obs(obsda_sort%set(obs_idx_TCX))%dif(obsda_sort%idx(obs_idx_TCX)) == &
!    obs(obsda_sort%set(obs_idx_TCY))%dif(obsda_sort%idx(obs_idx_TCY)) .and. &
!    obs(obsda_sort%set(obs_idx_TCY))%dif(obsda_sort%idx(obs_idx_TCY)) == &
!    obs(obsda_sort%set(obs_idx_TCP))%dif(obsda_sort%idx(obs_idx_TCP)) .and. & 
!    (DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
!    abs(obs(obsda_sort%set(obs_idx_TCX))%dif(obsda_sort%idx(obs_idx_TCX))) <= DEPARTURE_STAT_T_RANGE))then
!
!    allocate(bTC(3,0:nprocs_d-1))
!    allocate(bufr(3,0:nproces_d-1))
!
!    bTC = 9.99d33
!    bufr = 9.99d33
!
!!!!!!    call search_tc_subdom(obsda_sort%ri(obs_idx_TCX),obsda_sort%rj(obs_idx_TCX),v2dg,bTC(1,PRC_myrank),bTC(2,PRC_myrank),bTC(3,PRC_myrank))
!
!    CALL MPI_BARRIER(MPI_COMM_d,ierr)
!    CALL MPI_ALLREDUCE(bTC,bufr,3*nprocs_d,MPI_r_size,MPI_MIN,MPI_COMM_d,ierr)
!    bTC = bufr
!
!    deallocate(bufr)
!
!
!    ! Assume MSLP of background TC is lower than 1100 (hPa). 
!    bTC_mslp = 1100.0d2
!    do n = 0, nprocs_d - 1
!      write(6,'(3e20.5)')bTC(1,n),bTC(2,n),bTC(3,n) ! debug
!      if (bTC(3,n) < bTC_mslp ) then
!        bTC_mslp = bTC(3,n)
!        bTC_proc = n
!      endif
!    enddo ! [ n = 0, nprocs_d - 1]
!
!    do n = 1, 3
!      if(n==1) i = obs_idx_TCX
!      if(n==2) i = obs_idx_TCY
!      if(n==3) i = obs_idx_TCP
!
!      ohx(i) = bTC(n,bTC_proc)
!      oqc(i) = iqc_otype
!      if(bTC_MSLP < 1100.0d2) oqc(i) = iqc_good
!
!      if (oqc(i) == iqc_good) then
!        ohx(i) = obs(obsda_sort%set(i))%dat(obsda_sort%idx(i)) - ohx(i)
!      end if
!    enddo
!
!  endif ! [DEPARTURE_STAT_T_RANGE]


  call monit_dep(nnobs,oelm,ohx,oqc,nobs,bias,rmse)

  monit_type = .false.
  monit_type(uid_obs(id_u_obs)) = .true.
  monit_type(uid_obs(id_v_obs)) = .true.
  monit_type(uid_obs(id_t_obs)) = .true.
  monit_type(uid_obs(id_tv_obs)) = .true.
  monit_type(uid_obs(id_q_obs)) = .true.
!  monit_type(uid_obs(id_rh_obs)) = .true.
  monit_type(uid_obs(id_ps_obs)) = .true.
  if (DEPARTURE_STAT_RADAR) then
    monit_type(uid_obs(id_radar_ref_obs)) = .true.
    monit_type(uid_obs(id_radar_ref_zero_obs)) = .true.
    monit_type(uid_obs(id_radar_vr_obs)) = .true.
!    monit_type(uid_obs(id_radar_prh_obs)) = .true.
  end if
  if (DEPARTURE_STAT_H08) then
    monit_type(uid_obs(id_H08IR_obs)) = .true.
  end if

  deallocate (oelm)
  deallocate (ohx)
  deallocate (oqc)

  return
end subroutine monit_obs

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
subroutine read_obs_all(obs)
  implicit none

  type(obs_info), intent(out) :: obs(OBS_IN_NUM)
  character(14) :: timelabel
  integer :: iof
  logical :: ex

  character(filelenmax) :: obsfile

  call date2tlab_LETKF(adate,timelabel)

  do iof = 1, OBS_IN_NUM
    if (OBS_IN_FORMAT(iof) /= obsfmt_pawr_toshiba) then
      inquire (file=trim(OBS_IN_NAME(iof))//trim(timelabel_obs), exist=ex)
      if (.not. ex) then
        write(6,*) '[Warning] FILE ',trim(OBS_IN_NAME(iof))//trim(timelabel_obs),' NOT FOUND'


        obs(iof)%nobs = 0
        call obs_info_allocate(obs(iof), extended=.true.) !!! check why this is necessary !!!


        cycle
      end if
    end if

    select case (OBS_IN_FORMAT(iof))
    case (obsfmt_prepbufr)
      call get_nobs(trim(OBS_IN_NAME(iof))//trim(timelabel_obs),8,obs(iof)%nobs)
    case (obsfmt_radar)
      call get_nobs_radar(trim(OBS_IN_NAME(iof))//trim(timelabel_obs), obs(iof)%nobs, obs(iof)%meta(1), obs(iof)%meta(2), obs(iof)%meta(3))
    case (obsfmt_pawr_toshiba)
      if (.not. OBS_USE_JITDT) then
        call read_obs_radar_toshiba(trim(OBS_IN_NAME(iof))//trim(timelabel_obs), obs(iof))
      end if
    case (obsfmt_h08)
      call get_nobs_H08(trim(OBS_IN_NAME(iof))//trim(timelabel_obs),obs(iof)%nobs) ! H08
    case default
      write(6,*) '[Error] Unsupported observation file format!'
      stop
    end select

    write(6,'(5A,I9,A)') 'OBS FILE [', trim(OBS_IN_NAME(iof))//trim(timelabel_obs), '] (FORMAT ', &
                         trim(OBS_IN_FORMAT(iof)), '): TOTAL ', &
                         obs(iof)%nobs, ' OBSERVATIONS'

    select case (OBS_IN_FORMAT(iof))
    case (obsfmt_prepbufr)
      call obs_info_allocate(obs(iof), extended=.true.)
      call read_obs(trim(OBS_IN_NAME(iof))//trim(timelabel_obs),obs(iof))
    case (obsfmt_radar)
      call obs_info_allocate(obs(iof), extended=.true.)
      call read_obs_radar(trim(OBS_IN_NAME(iof))//trim(timelabel_obs),obs(iof))
    !case (obsfmt_pawr_toshiba)
    !  data already read by 'read_obs_radar_toshiba' above
    case (obsfmt_h08)
      call obs_info_allocate(obs(iof), extended=.true.)
      call read_obs_H08(trim(OBS_IN_NAME(iof))//trim(timelabel_obs),obs(iof)) ! H08
    end select
  end do ! [ iof = 1, OBS_IN_NUM ]

  return
end subroutine read_obs_all

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
subroutine write_obs_all(obs, missing, file_suffix)
  implicit none

  type(obs_info), intent(in) :: obs(OBS_IN_NUM)
  logical, intent(in), optional :: missing
  character(len=*), intent(in), optional :: file_suffix
  logical :: missing_
  integer :: iof, strlen1, strlen2
  character(len=200) :: filestr

  missing_ = .true.
  IF(present(missing)) missing_ = missing

  do iof = 1, OBS_IN_NUM
    if (present(file_suffix)) then
      strlen1 = len(trim(OBS_IN_NAME(iof))//trim(timelabel_obs))
      strlen2 = len(trim(file_suffix))
      write (filestr(1:strlen1),'(A)') trim(OBS_IN_NAME(iof))//trim(timelabel_obs)
      write (filestr(strlen1+1:strlen1+strlen2),'(A)') trim(file_suffix)
    else
      filestr = trim(OBS_IN_NAME(iof))//trim(timelabel_obs)
    end if
    select case (OBS_IN_FORMAT(iof))
    case (obsfmt_prepbufr)
      call write_obs(trim(filestr),obs(iof),missing=missing_)
    case (obsfmt_radar)
      call write_obs_radar(trim(filestr),obs(iof),missing=missing_)
    case (obsfmt_h08)
      call write_obs_H08(trim(filestr),obs(iof),missing=missing_) ! H08
    end select
  end do ! [ iof = 1, OBS_IN_NUM ]

  return
end subroutine write_obs_all

!-------------------------------------------------------------------------------
! Read all observation data from files of various formats
!-------------------------------------------------------------------------------
subroutine read_obs_all_mpi(obs)
  implicit none
  type(obs_info), intent(out) :: obs(OBS_IN_NUM)
  integer :: iof, ierr

  call mpi_timer('', 2)

  if (myrank_a == 0) then
    call read_obs_all(obs)

    call mpi_timer('read_obs_all_mpi:read_obs_all:', 2)
  end if

  call mpi_timer('', 2, barrier=MPI_COMM_a)

  do iof = 1, OBS_IN_NUM
    if (OBS_IN_FORMAT(iof) == obsfmt_pawr_toshiba .and. OBS_USE_JITDT) then
      cycle
    end if

    call MPI_BCAST(obs(iof)%nobs, 1, MPI_INTEGER, 0, MPI_COMM_a, ierr)
    if (myrank_a /= 0) then
      call obs_info_allocate(obs(iof), extended=.true.)
    end if

    call MPI_BCAST(obs(iof)%elm, obs(iof)%nobs, MPI_INTEGER, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%lon, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%lat, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%lev, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%dat, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%err, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%typ, obs(iof)%nobs, MPI_INTEGER, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%dif, obs(iof)%nobs, MPI_r_size, 0, MPI_COMM_a, ierr)
    call MPI_BCAST(obs(iof)%meta, max_obs_info_meta, MPI_r_size, 0, MPI_COMM_a, ierr)
  end do ! [ iof = 1, OBS_IN_NUM ]

  call mpi_timer('read_obs_all_mpi:mpi_bcast:', 2)

  do iof = 1, OBS_IN_NUM
    select case (OBS_IN_FORMAT(iof))
    case (obsfmt_pawr_toshiba)
      if (OBS_USE_JITDT) then
        call read_obs_radar_toshiba(trim(OBS_IN_NAME(iof)), obs(iof))
      end if
    end select
  end do ! [ iof = 1, OBS_IN_NUM ]

  if (OBS_USE_JITDT) then
    call mpi_timer('read_obs_all_mpi:read_obs_jitdt:', 2)
  end if

  return
end subroutine read_obs_all_mpi

!-------------------------------------------------------------------------------
! Read all observation data from files of various formats
!-------------------------------------------------------------------------------
subroutine get_nobs_da_mpi(nobs)
  implicit none
  integer, intent(out) :: nobs
  character(filelenmax) :: obsdafile
  character(11) :: obsda_suffix = '.000000.dat'
  integer :: ierr

! read from all available data by every processes
!-----------------------------
!  if ((myrank_to_mem(1) >= 1 .and. myrank_to_mem(1) <= MEMBER) .or. &
!      myrank_to_mem(1) == mmdetin) then
!    if (myrank_to_mem(1) <= MEMBER) then
!      obsdafile = OBSDA_IN_BASENAME
!      call filename_replace_mem(obsdafile, myrank_to_mem(1))
!    else if (myrank_to_mem(1) == mmean) then
!      obsdafile = OBSDA_MEAN_IN_BASENAME
!    else if (myrank_to_mem(1) == mmdet) then
!      obsdafile = OBSDA_MDET_IN_BASENAME
!    end if
!    write (obsda_suffix(2:7), '(I6.6)') myrank_d
!#ifdef H08
!    call get_nobs(trim(obsdafile) // obsda_suffix, 6, nobs) ! H08
!#else
!    call get_nobs(trim(obsdafile) // obsda_suffix, 4, nobs)
!#endif
!  end if

! read by process 0 and broadcast
!-----------------------------
  if (myrank_e == 0) then
    obsdafile = OBSDA_IN_BASENAME
    call filename_replace_mem(obsdafile, 1)
    write (obsda_suffix(2:7), '(I6.6)') myrank_d
#ifdef H08
    call get_nobs(trim(obsdafile) // obsda_suffix, 6, nobs) ! H08
#else
    call get_nobs(trim(obsdafile) // obsda_suffix, 4, nobs)
#endif
  end if
  call MPI_BCAST(nobs, 1, MPI_INTEGER, 0, MPI_COMM_e, ierr)
!-----------------------------

  return
end subroutine get_nobs_da_mpi

!-------------------------------------------------------------------------------
! Partially reduce observations processed in the same processes in the iteration
!-------------------------------------------------------------------------------
#ifdef H08
subroutine obs_da_value_partial_reduce_iter(obsda, iter, nstart, nobs, ensval, qc, lev, val2)
#else
subroutine obs_da_value_partial_reduce_iter(obsda, iter, nstart, nobs, ensval, qc)
#endif
  implicit none
  type(obs_da_value), intent(inout) :: obsda
  integer, intent(in)      :: iter
  integer, intent(in)      :: nstart
  integer, intent(in)      :: nobs
  real(r_size), intent(in) :: ensval(nobs)
  integer, intent(in)      :: qc(nobs)
#ifdef H08
  real(r_size), intent(in) :: lev(nobs)
  real(r_size), intent(in) :: val2(nobs)
#endif
  integer :: nend
  integer :: im

  if (nobs <= 0) then
    return
  end if
  nend = nstart + nobs - 1
  im = myrank_to_mem(iter)
  if (.not. ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin)) then
    return
  end if

  ! variables with an ensemble dimension
  obsda%ensval(iter,nstart:nend) = ensval

  ! variables without an ensemble dimension
  obsda%qc(nstart:nend) = max(obsda%qc(nstart:nend), qc)
#ifdef H08
  if (im <= MEMBER) then ! only consider lev, val2 from members, not from the means
    obsda%lev(nstart:nend) = obsda%lev(nstart:nend) + lev
    obsda%val2(nstart:nend) = obsda%val2(nstart:nend) + val2
  end if
#endif

  return
end subroutine obs_da_value_partial_reduce_iter

!-------------------------------------------------------------------------------
! Allreduce observations to obtain complete ensemble values
!-------------------------------------------------------------------------------
subroutine obs_da_value_allreduce(obsda)
  implicit none
  type(obs_da_value), intent(inout) :: obsda
  real(r_size), allocatable :: ensval_bufs(:,:)
  real(r_size), allocatable :: ensval_bufr(:,:)
  integer :: cnts
  integer :: cntr(nprocs_e)
  integer :: dspr(nprocs_e)
  integer :: current_shape(2)
  integer :: ie, it, im, imb, ierr

  if (obsda%nobs <= 0) then
    return
  end if

  call mpi_timer('', 3)

  ! variables with an ensemble dimension
  cntr(:) = 0
  do ie = 1, nprocs_e
    do it = 1, nitmax
      im = ranke_to_mem(it, ie)
      if ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin) then
        cntr(ie) = cntr(ie) + 1
      end if
    end do
  end do
  allocate (ensval_bufs(obsda%nobs, cntr(myrank_e+1)))
  allocate (ensval_bufr(obsda%nobs, nensobs))

  do im = 1, cntr(myrank_e+1)
    ensval_bufs(:,im) = obsda%ensval(im,:)
  end do

  cntr(:) = cntr(:) * obsda%nobs
  cnts = cntr(myrank_e+1)
  dspr(1) = 0
  do ie = 2, nprocs_e
    dspr(ie) = dspr(ie-1) + cntr(ie-1)
  end do

  call mpi_timer('obs_da_value_allreduce:copy_bufs:', 3, barrier=MPI_COMM_e)

  call MPI_ALLGATHERV(ensval_bufs, cnts, MPI_r_size, ensval_bufr, cntr, dspr, MPI_r_size, MPI_COMM_e, ierr)

  call mpi_timer('obs_da_value_allreduce:mpi_allgatherv:', 3)

  current_shape = shape(obsda%ensval)
  if (current_shape(1) < nensobs) then
    deallocate (obsda%ensval)
    allocate (obsda%ensval(nensobs, obsda%nobs))
  end if

  imb = 0
  do ie = 1, nprocs_e
    do it = 1, nitmax
      im = ranke_to_mem(it, ie)
      if ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin) then
        imb = imb + 1
        if (im == mmdetin) then
          obsda%ensval(mmdetobs,:) = ensval_bufr(:,imb)
        else
          obsda%ensval(im,:) = ensval_bufr(:,imb)
        end if
      end if
    end do
  end do
  deallocate(ensval_bufs, ensval_bufr)

  call mpi_timer('obs_da_value_allreduce:copy_bufr:', 3, barrier=MPI_COMM_e)

  ! variables without an ensemble dimension
  if (nprocs_e > 1) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, obsda%qc(:), obsda%nobs, MPI_INTEGER, MPI_MAX, MPI_COMM_e, ierr)
  end if
#ifdef H08
  if (nprocs_e > 1) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, obsda%lev(:), obsda%nobs, MPI_r_size, MPI_SUM, MPI_COMM_e, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, obsda%val2(:), obsda%nobs, MPI_r_size, MPI_SUM, MPI_COMM_e, ierr)
  end if
  obsda%lev(:) = obsda%lev(:) / real(MEMBER, r_size)
  obsda%val2(:) = obsda%val2(:) / real(MEMBER, r_size)
#endif

  call mpi_timer('obs_da_value_allreduce:mpi_allreduce:', 3)

  return
end subroutine obs_da_value_allreduce

!=======================================================================
end module obs_tools
