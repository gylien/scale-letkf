MODULE obsope_tools
!=======================================================================
!
! [PURPOSE:] Observation operator tools
!
! [HISTORY:]
!   November 2014  Guo-Yuan Lien  created
!   .............  See git history for the following revisions
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
  USE common_obs_scale
  use common_nml
!  use scale_process, only: &
!    PRC_myrank
!    MPI_COMM_d => LOCAL_COMM_WORLD
  use scale_grid_index, only: &
    KHALO, IHALO, JHALO

  IMPLICIT NONE
  PUBLIC

CONTAINS

!-----------------------------------------------------------------------
! Observation operator calculation
!-----------------------------------------------------------------------
SUBROUTINE obsope_cal(obsda, obsda_return, nobs_extern)
  IMPLICIT NONE

  type(obs_da_value), intent(out) :: obsda
  logical, intent(in) :: obsda_return
  integer, optional, intent(in) :: nobs_extern

  integer :: it, im, iof, islot, ierr
  integer :: n, nn, nn_0, nsub, nmod, n1, n2

  integer :: nobs     ! observation number processed in this subroutine
  integer :: nobs_all
  integer :: nobs_max_per_file
  integer :: nobs_max_per_file_sub
  integer :: slot_nobsg

  integer :: ip, ibufs
  integer, allocatable :: cntr(:), dspr(:)
  integer, allocatable :: cnts(:), dsps(:)
  integer, allocatable :: bsn(:,:), bsna(:,:), bsnext(:,:)
  integer :: islot_time_out, islot_domain_out

  integer, allocatable :: obrank(:)
  real(r_size), allocatable :: obri(:)
  real(r_size), allocatable :: obrj(:)
  integer, allocatable :: obrank_bufs(:)
  real(r_size), allocatable :: ri_bufs(:)
  real(r_size), allocatable :: rj_bufs(:)

  integer, allocatable :: obset_bufs(:)
  integer, allocatable :: obidx_bufs(:)
  real(r_size), allocatable :: ri_bufs2(:)
  real(r_size), allocatable :: rj_bufs2(:)

  integer :: slot_id(SLOT_START:SLOT_END)
  real(r_size) :: slot_lb(SLOT_START:SLOT_END)
  real(r_size) :: slot_ub(SLOT_START:SLOT_END)

  real(r_size), allocatable :: v3dg(:,:,:,:)
  real(r_size), allocatable :: v2dg(:,:,:)

  integer, allocatable :: qc_p(:)
#ifdef H08
  real(r_size), allocatable :: lev_p(:)
  real(r_size), allocatable :: val2_p(:)
#endif

  real(r_size) :: ril, rjl, rk

  character(filelenmax) :: obsdafile
  character(len=11) :: obsda_suffix = '.000000.dat'
  character(len=4) :: nstr
  character(len=timer_name_width) :: timer_str

#ifdef H08
! -- for Himawari-8 obs --
  integer :: nallprof ! Maximum number of Him8 profiles required for RTTOV
  real(r_size), allocatable :: ri_H08(:),rj_H08(:)
  real(r_size), allocatable :: lon_H08(:),lat_H08(:)
  integer, allocatable :: nnB07(:) ! index of Him8 band 7
  integer :: nprof ! num of Him8 profile
  real(r_size), allocatable :: yobs_H08(:,:),plev_obs_H08(:,:)
  real(r_size), allocatable :: yobs_H08_clr(:,:)
  integer, allocatable :: qc_H08(:,:)
  integer :: ch

#endif

! -- for TC vital assimilation --
! bTC: background TC in each subdomain
! bTC(1,:) : tcx (m), bTC(2,:): tcy (m), bTC(3,:): mslp (Pa)
  real(r_size),allocatable :: bTC(:,:)
  real(r_size) :: bTC_mslp

! Multiple TCs are not considered (04/14/2017)
  real(r_size) :: TC_rij(2) = -1.0d0
  integer :: bTC_rank_d ! the process where the background TC is located.
  integer :: obs_nn_TCX, obs_nn_TCY, obs_nn_TCP ! TCX, TCY, TCP

!-------------------------------------------------------------------------------

  call mpi_timer('', 2)

!-------------------------------------------------------------------------------
! First scan of all observation data: Compute their horizontal location and time
!-------------------------------------------------------------------------------

  nobs_all = 0
  nobs_max_per_file = 0
  do iof = 1, OBS_IN_NUM
    if (OBSDA_RUN(iof)) then
      nobs_all = nobs_all + obs(iof)%nobs
      if (obs(iof)%nobs > nobs_max_per_file) then
        nobs_max_per_file = obs(iof)%nobs
      end if
    end if
  end do

  allocate (obrank(nobs_all))
  allocate (obri(nobs_all))
  allocate (obrj(nobs_all))

  nobs_max_per_file_sub = (nobs_max_per_file - 1) / nprocs_a + 1
  allocate (obrank_bufs(nobs_max_per_file_sub))
  allocate (ri_bufs(nobs_max_per_file_sub))
  allocate (rj_bufs(nobs_max_per_file_sub))

  allocate (cntr(nprocs_a))
  allocate (dspr(nprocs_a))

  ! Use all processes to compute the basic obsevration information
  ! (locations in model grids and the subdomains they belong to)
  !-----------------------------------------------------------------------------

  nn_0 = 0
  do iof = 1, OBS_IN_NUM
    if (OBSDA_RUN(iof) .and. obs(iof)%nobs > 0) then
      nsub = obs(iof)%nobs / nprocs_a
      nmod = mod(obs(iof)%nobs, nprocs_a)
      do ip = 1, nmod
        cntr(ip) = nsub + 1
      end do
      do ip = nmod+1, nprocs_a
        cntr(ip) = nsub
      end do
      dspr(1) = nn_0
      do ip = 2, nprocs_a
        dspr(ip) = dspr(ip-1) + cntr(ip-1)
      end do

      n1 = dspr(myrank_a+1) - nn_0 + 1
      n2 = dspr(myrank_a+1) - nn_0 + cntr(myrank_a+1)

      obrank_bufs(:) = -1

      ibufs = 0
      do n = n1, n2
        ibufs = ibufs + 1

        call phys2ij(obs(iof)%lon(n), obs(iof)%lat(n), ri_bufs(ibufs), rj_bufs(ibufs))
        call rij_g2l_auto(obrank_bufs(ibufs), ri_bufs(ibufs), rj_bufs(ibufs), ril, rjl) ! rij, rjl discarded here; re-computed later
#ifdef TCV
        if (obs(iof)%elm(n) == id_tcmip_obs)then
          TC_rij(1) = ri_bufs(ibufs)
          TC_rij(2) = rj_bufs(ibufs)
        endif
#endif
      end do ! [ n = n1, n2 ]

      call mpi_timer('obsope_cal:first_scan_cal:', 2, barrier=MPI_COMM_a)

      call MPI_GATHERV(obrank_bufs, cntr(myrank_a+1), MPI_INTEGER, obrank, cntr, dspr, MPI_INTEGER, 0, MPI_COMM_a, ierr)
      call MPI_GATHERV(ri_bufs,     cntr(myrank_a+1), MPI_r_size , obri,   cntr, dspr, MPI_r_size,  0, MPI_COMM_a, ierr)
      call MPI_GATHERV(rj_bufs,     cntr(myrank_a+1), MPI_r_size , obrj,   cntr, dspr, MPI_r_size,  0, MPI_COMM_a, ierr)

      call mpi_timer('obsope_cal:first_scan_reduce:', 2)

      nn_0 = nn_0 + obs(iof)%nobs
    end if ! [ OBSDA_RUN(iof) .and. obs(iof)%nobs > 0 ]
  end do ! [ do iof = 1, OBS_IN_NUM ]

  deallocate (cntr, dspr)
  deallocate (obrank_bufs, ri_bufs, rj_bufs)

#ifdef TCV
  call MPI_ALLREDUCE(MPI_IN_PLACE, TC_rij, 2, MPI_r_size, MPI_MAX, MPI_COMM_a, ierr)  ! 
  write(6,'(a,2f9.1)')"DEBUG",TC_rij(1),TC_rij(2)
#endif

  ! Bucket sort of observation wrt. time slots and subdomains using the process rank 0
  !-----------------------------------------------------------------------------

  islot_time_out = SLOT_END + 1   ! slot = SLOT_END+1 for observation not in the assimilation time window
  islot_domain_out = SLOT_END + 2 ! slot = SLOT_END+2 for observation outside of the model domain

  allocate (bsn   (SLOT_START  :SLOT_END+2, 0:nprocs_d-1))
  allocate (bsna  (SLOT_START-1:SLOT_END+2, 0:nprocs_d-1))

  if (myrank_e == 0) then
    allocate ( obset_bufs(nobs_all) )
    allocate ( obidx_bufs(nobs_all) )
    allocate ( ri_bufs2(nobs_all) )
    allocate ( rj_bufs2(nobs_all) )
  end if

  if (myrank_a == 0) then
    allocate (bsnext(SLOT_START  :SLOT_END+2, 0:nprocs_d-1))
    bsn(:,:) = 0
    bsna(:,:) = 0
    bsnext(:,:) = 0

    nn = 0
    do iof = 1, OBS_IN_NUM
      if (OBSDA_RUN(iof) .and. obs(iof)%nobs > 0) then
        do n = 1, obs(iof)%nobs
          nn = nn + 1
          if (obrank(nn) == -1) then
            ! process the observations outside of the model domain in process rank 0
            bsn(islot_domain_out, 0) = bsn(islot_domain_out, 0) + 1
          else
            islot = ceiling(obs(iof)%dif(n) / SLOT_TINTERVAL - 0.5d0) + SLOT_BASE
            if (islot < SLOT_START .or. islot > SLOT_END) then
              islot = islot_time_out
            end if
            bsn(islot, obrank(nn)) = bsn(islot, obrank(nn)) + 1
          end if
        end do ! [ n = 1, obs(iof)%nobs ]
      end if ! [ OBSDA_RUN(iof) .and. obs(iof)%nobs > 0 ]
    end do ! [ do iof = 1, OBS_IN_NUM ]

    do ip = 0, nprocs_d-1
      if (ip > 0) then
        bsna(SLOT_START-1, ip) = bsna(SLOT_END+2, ip-1)
      end if
      do islot = SLOT_START, SLOT_END+2
        bsna(islot, ip) = bsna(islot-1, ip) + bsn(islot, ip)
      end do
      bsnext(SLOT_START:SLOT_END+2, ip) = bsna(SLOT_START-1:SLOT_END+1, ip)
    end do

    nn = 0
    do iof = 1, OBS_IN_NUM
      if (OBSDA_RUN(iof) .and. obs(iof)%nobs > 0) then
        do n = 1, obs(iof)%nobs
          nn = nn + 1
          if (obrank(nn) == -1) then
            ! process the observations outside of the model domain in process rank 0
            bsnext(islot_domain_out, 0) = bsnext(islot_domain_out, 0) + 1
            obset_bufs(bsnext(islot_domain_out, 0)) = iof
            obidx_bufs(bsnext(islot_domain_out, 0)) = n
            ri_bufs2(bsnext(islot_domain_out, 0)) = obri(nn)
            rj_bufs2(bsnext(islot_domain_out, 0)) = obrj(nn)
          else
            islot = ceiling(obs(iof)%dif(n) / SLOT_TINTERVAL - 0.5d0) + SLOT_BASE
            if (islot < SLOT_START .or. islot > SLOT_END) then
              islot = islot_time_out
            end if
            bsnext(islot, obrank(nn)) = bsnext(islot, obrank(nn)) + 1
            obset_bufs(bsnext(islot, obrank(nn))) = iof
            obidx_bufs(bsnext(islot, obrank(nn))) = n
            ri_bufs2(bsnext(islot, obrank(nn))) = obri(nn)
            rj_bufs2(bsnext(islot, obrank(nn))) = obrj(nn)
          end if
        end do ! [ n = 1, obs(iof)%nobs ]
      end if ! [ OBSDA_RUN(iof) .and. obs(iof)%nobs > 0 ]
    end do ! [ do iof = 1, OBS_IN_NUM ]

    deallocate (bsnext)

    call mpi_timer('obsope_cal:bucket_sort:', 2)
  end if ! [ myrank_a == 0 ]

  deallocate ( obrank, obri, obrj )

  ! Broadcast the bucket-sort observation numbers to all processes and print
  !-----------------------------------------------------------------------------

  call mpi_timer('', 2, barrier=MPI_COMM_a)

  call MPI_BCAST(bsn,  (SLOT_END-SLOT_START+3)*nprocs_d, MPI_INTEGER, 0, MPI_COMM_a, ierr)
  call MPI_BCAST(bsna, (SLOT_END-SLOT_START+4)*nprocs_d, MPI_INTEGER, 0, MPI_COMM_a, ierr)

  call mpi_timer('obsope_cal:sort_info_bcast:', 2)

  do islot = SLOT_START, SLOT_END
    slot_id(islot) = islot - SLOT_START + 1
    slot_lb(islot) = (real(islot - SLOT_BASE, r_size) - 0.5d0) * SLOT_TINTERVAL
    slot_ub(islot) = (real(islot - SLOT_BASE, r_size) + 0.5d0) * SLOT_TINTERVAL
  end do

  write (nstr, '(I4)') SLOT_END - SLOT_START + 1
  write (6, *)
  write (6, '(A,I6,A)') 'OBSERVATION COUNTS BEFORE QC (FROM OBSOPE):'
  write (6, '(A,'//nstr//"('=========='),A)") '====================', '===================='
  write (6, '(A,'//nstr//'I10.4)') '            SLOT #  ', slot_id(:)
  write (6, '(A,'//nstr//'F10.1)') '            FROM (s)', slot_lb(:)
  write (6, '(A,'//nstr//'F10.1,A)') 'SUBDOMAIN #   TO (s)', slot_ub(:), '  OUT_TIME     TOTAL'
  write (6, '(A,'//nstr//"('----------'),A)") '--------------------', '--------------------'
  do ip = 0, nprocs_d-1
    write (6, '(I11.6,9x,'//nstr//'I10,2I10)') ip, bsn(SLOT_START:SLOT_END, ip), bsn(islot_time_out, ip), bsna(SLOT_END+1, ip) - bsna(SLOT_START-1, ip)
  end do
  write (6, '(A,'//nstr//"('----------'),A)") '--------------------', '--------------------'
  write (6, '(A,'//nstr//'(10x),10x,I10)') ' OUT_DOMAIN         ', bsn(islot_domain_out, 0)
  write (6, '(A,'//nstr//"('----------'),A)") '--------------------', '--------------------'
  write (6, '(A,'//nstr//'I10,2I10)') '      TOTAL         ', sum(bsn(SLOT_START:SLOT_END, :), dim=2), sum(bsn(islot_time_out, :)), bsna(SLOT_END+2, nprocs_d-1)
  write (6, '(A,'//nstr//"('=========='),A)") '====================', '===================='

  ! Scatter the basic obsevration information to processes group {myrank_e = 0},
  ! each of which only gets the data in its own subdomain
  !-----------------------------------------------------------------------------

  nobs = bsna(SLOT_END+2, myrank_d) - bsna(SLOT_START-1, myrank_d)
  if (present(nobs_extern)) then
    obsda%nobs = nobs + nobs_extern
  else
    obsda%nobs = nobs
  end if

  if (obsda_return) then
    call obs_da_value_allocate(obsda, nensobs)
  else
    call obs_da_value_allocate(obsda, 0)
  end if

  if (myrank_e == 0) then
    allocate (cnts(nprocs_d))
    allocate (dsps(nprocs_d))
    do ip = 0, nprocs_d-1
      dsps(ip+1) = bsna(SLOT_START-1, ip)
      cnts(ip+1) = bsna(SLOT_END+2, ip) - dsps(ip+1)
    end do

    call MPI_SCATTERV(obset_bufs, cnts, dsps, MPI_INTEGER, obsda%set, cnts(myrank_d+1), MPI_INTEGER, 0, MPI_COMM_d, ierr)
    call MPI_SCATTERV(obidx_bufs, cnts, dsps, MPI_INTEGER, obsda%idx, cnts(myrank_d+1), MPI_INTEGER, 0, MPI_COMM_d, ierr)
    call MPI_SCATTERV(ri_bufs2,   cnts, dsps, MPI_r_size,  obsda%ri,  cnts(myrank_d+1), MPI_r_size,  0, MPI_COMM_d, ierr)
    call MPI_SCATTERV(rj_bufs2,   cnts, dsps, MPI_r_size,  obsda%rj,  cnts(myrank_d+1), MPI_r_size,  0, MPI_COMM_d, ierr)

    call mpi_timer('obsope_cal:mpi_scatterv:', 2)

    deallocate (cnts, dsps)
    deallocate (obset_bufs, obidx_bufs, ri_bufs2, rj_bufs2)
  end if ! [ myrank_e == 0 ]

  ! Broadcast the basic obsevration information
  ! from processes group {myrank_e = 0} to all processes
  !-----------------------------------------------------------------------------

  call mpi_timer('', 2, barrier=MPI_COMM_e)

  call MPI_BCAST(obsda%set, nobs, MPI_INTEGER, 0, MPI_COMM_e, ierr)
  call MPI_BCAST(obsda%idx, nobs, MPI_INTEGER, 0, MPI_COMM_e, ierr)
  call MPI_BCAST(obsda%ri,  nobs, MPI_r_size,  0, MPI_COMM_e, ierr)
  call MPI_BCAST(obsda%rj,  nobs, MPI_r_size,  0, MPI_COMM_e, ierr)

  call mpi_timer('obsope_cal:mpi_broadcast:', 2)

!-------------------------------------------------------------------------------
! Second scan of observation data in own subdomain: Compute H(x), QC, ... etc.
!-------------------------------------------------------------------------------

  if (obsda_return .and. nitmax > 1 .and. nobs > 0) then
    allocate (qc_p(nobs))
#ifdef H08
    allocate (lev_p(nobs))
    allocate (val2_p(nobs))
#endif
  end if

  allocate ( v3dg (nlevh,nlonh,nlath,nv3dd) )
  allocate ( v2dg (nlonh,nlath,nv2dd) )

  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    if ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin) then

      write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is processing member ', &
            im, ', subdomain id #', proc2mem(2,it,myrank+1)

      if (nobs > 0) then
        obsda%qc(1:nobs) = iqc_undef
#ifdef H08
        obsda%lev(1:nobs) = 0.0d0
#endif
      end if

      ! Observations not in the assimilation time window
      ! 
      n1 = bsna(islot_time_out-1, myrank_d) - bsna(SLOT_START-1, myrank_d) + 1
      n2 = bsna(islot_time_out,   myrank_d) - bsna(SLOT_START-1, myrank_d)
      if (n1 <= n2) then
        obsda%qc(n1:n2) = iqc_time
      end if

      ! Observations outside of the model domain
      ! 
      n1 = bsna(islot_domain_out-1, myrank_d) - bsna(SLOT_START-1, myrank_d) + 1
      n2 = bsna(islot_domain_out,   myrank_d) - bsna(SLOT_START-1, myrank_d)
      if (n1 <= n2) then
        obsda%qc(n1:n2) = iqc_out_h
      end if

      ! Valid observations: loop over time slots
      ! 
      do islot = SLOT_START, SLOT_END
        write (6, '(A,I3,A,F9.1,A,F9.1,A)') 'Slot #', islot-SLOT_START+1, ': time window (', slot_lb(islot), ',', slot_ub(islot), '] sec'

        n1 = bsna(islot-1, myrank_d) - bsna(SLOT_START-1, myrank_d) + 1
        n2 = bsna(islot,   myrank_d) - bsna(SLOT_START-1, myrank_d)
        slot_nobsg = sum(bsn(islot, :))

        if (slot_nobsg <= 0) then
          write (6, '(A)') ' -- no observations found in this time slot... do not need to read model data'
          cycle
        end if

        write (6, '(A,I10)') ' -- # obs in the slot = ', slot_nobsg
        write (6, '(A,I6,A,I6,A,I10)') ' -- # obs in the slot and processed by rank ', myrank, ' (subdomain #', myrank_d, ') = ', bsn(islot, myrank_d)

        call read_ens_history_iter(it, islot, v3dg, v2dg)

        write (timer_str, '(A30,I4,A7,I4,A2)') 'obsope_cal:read_ens_history(t=', it, ', slot=', islot, '):'
        call mpi_timer(trim(timer_str), 2)

#ifdef TCV
        obs_nn_TCX = -1
        obs_nn_TCY = -1
        obs_nn_TCP = -1
#endif

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(nn,n,iof,ril,rjl,rk)
            do nn = n1, n2
              iof = obsda%set(nn)
              n = obsda%idx(nn)


#ifdef H08
              !! Him8 obs will be processed separatly !!
              if (obs(iof)%elm(n) == id_H08IR_obs) cycle 
#endif


#ifdef TCV
              !! TC vital !!
              select case (obs(iof)%elm(n))
              case (id_tclon_obs)
                obs_nn_TCX = nn
                cycle
              case (id_tclat_obs)
                obs_nn_TCY = nn
                cycle
              case (id_tcmip_obs)
                obs_nn_TCP = nn
                cycle
              end select
#endif

              call rij_g2l(myrank_d, obsda%ri(nn), obsda%rj(nn), ril, rjl)

              if (.not. USE_OBS(obs(iof)%typ(n))) then
                obsda%qc(nn) = iqc_otype
                cycle
              end if

              if (obs(iof)%elm(n) == id_radar_ref_obs .or. obs(iof)%elm(n) == id_radar_ref_zero_obs .or. obs(iof)%elm(n) == id_radar_vr_obs) then
                if (obs(iof)%lev(n) > RADAR_ZMAX) then
                  obsda%qc(nn) = iqc_radar_vhi
#ifdef DEBUG
                  write(6,'(A,F8.1,A,I5)') 'warning: radar observation is too high: lev=', obs(iof)%lev(n), ', elem=', obs(iof)%elm(n)
#endif
                else
                  call phys2ijkz(v3dg(:,:,:,iv3dd_hgt), ril, rjl, obs(iof)%lev(n), rk, obsda%qc(nn))
                end if
              else
                call phys2ijk(v3dg(:,:,:,iv3dd_p), obs(iof)%elm(n), ril, rjl, obs(iof)%lev(n), rk, obsda%qc(nn))
              end if

              if (obsda%qc(nn) == iqc_good) then
                select case (OBS_IN_FORMAT(iof))
                case (1)
                  call Trans_XtoY(obs(iof)%elm(n), ril, rjl, rk, &
                                  obs(iof)%lon(n), obs(iof)%lat(n), v3dg, v2dg, obsda%val(nn), obsda%qc(nn))
                case (2)
                  call Trans_XtoY_radar(obs(iof)%elm(n), obs(iof)%meta(1), obs(iof)%meta(2), obs(iof)%meta(3), ril, rjl, rk, &
                                        obs(iof)%lon(n), obs(iof)%lat(n), obs(iof)%lev(n), v3dg, v2dg, obsda%val(nn), obsda%qc(nn))
                  if (obsda%qc(nn) == iqc_ref_low) obsda%qc(nn) = iqc_good ! when process the observation operator, we don't care if reflectivity is too small

                  !!!!!! may not need to do this at this stage !!!!!!
                  !if (obs(iof)%elm(n) == id_radar_ref_obs) then
                  !  obsda%val(nn) = 10.0d0 * log10(obsda%val(nn))
                  !end if
                  !!!!!!
                end select
              end if

            end do ! [ nn = n1, n2 ]
!$OMP END PARALLEL DO

#ifdef H08

            ! Him8 observations: count the number of profiles required for RTTOV
            ! 

            nallprof = int((n2 - n1 + 1) / nch)
            allocate(nnB07(nallprof))
            allocate(ri_H08(nallprof))
            allocate(rj_H08(nallprof))
            allocate(lon_H08(nallprof))
            allocate(lat_H08(nallprof))

            nprof = 0
            do nn = n1, n2
              iof = obsda%set(nn)
              n = obsda%idx(nn)

              if (obs(iof)%elm(n) /= id_H08IR_obs) cycle

              if (nint(obs(iof)%lev(n)) == 7) then
                nprof = nprof + 1
                nnB07(nprof) = nn
                call rij_g2l(myrank_d, obsda%ri(nn), obsda%rj(nn), ril, rjl)

                ri_H08(nprof) = ril
                rj_H08(nprof) = rjl
                lon_H08(nprof) = obs(iof)%lon(n)
                lat_H08(nprof) = obs(iof)%lat(n)
              endif

            end do ! [ nn = n1, n2 ]

            ! Him8 observation: apply radiative transfer model (RTTOV) 
            !             Note: OpenMP will be used within SCALE_H08_fwd
            !                    

            if(nprof >=1) then
              allocate(yobs_H08(nch,nprof))
              allocate(yobs_H08_clr(nch,nprof))
              allocate(plev_obs_H08(nch,nprof))
              allocate(qc_H08(nch,nprof))

              call Trans_XtoY_H08(nprof,ri_H08(1:nprof),rj_H08(1:nprof),&
                                  lon_H08(1:nprof),lat_H08(1:nprof),v3dg,v2dg,&
                                  yobs_H08,plev_obs_H08,&
                                  qc_H08,yobs_H08_clr=yobs_H08_clr)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(nn,ch,iof,n)
              do nn = 1, nprof
                do ch = 1, nch
                  obsda%val(nnB07(nn)+ch-1) = yobs_H08(ch,nn) 
                  obsda%qc(nnB07(nn)+ch-1) = qc_H08(ch,nn) 
                  obsda%lev(nnB07(nn)+ch-1) = plev_obs_H08(ch,nn) 
                  !obsda%val2(nnB07(nn)+ch-1) = yobs_H08_clr(ch,nn) 
                  ! Get CA (Okamoto et al. 2014)

                  iof = obsda%set(nnB07(nn)+ch-1)
                  n = obsda%idx(nnB07(nn)+ch-1)

                  obsda%val2(nnB07(nn)+ch-1) = (abs(yobs_H08(ch,nn) - yobs_H08_clr(ch,nn)) &
                                                + abs(obs(iof)%dat(n) - yobs_H08_clr(ch,nn))) * 0.5d0
                  !
                  ! Simple bias correction depending on the sky condition
                  ! (diagnosed by CA)
                  if(H08_BIAS_SIMPLE)then
                    if(obsda%val2(nnB07(nn)+ch-1) > H08_CA_THRES)then
                      obsda%val(nnB07(nn)+ch-1) = obsda%val(nnB07(nn)+ch-1) + H08_BIAS_CLOUD(ch)
                    else
                      obsda%val(nnB07(nn)+ch-1) = obsda%val(nnB07(nn)+ch-1) + H08_BIAS_CLEAR(ch)
                    endif
                  endif
                enddo
              enddo
!$OMP END PARALLEL DO

              deallocate(ri_H08, rj_H08)
              deallocate(lon_H08, lat_H08)
              deallocate(yobs_H08, plev_obs_H08)
              deallocate(yobs_H08_clr)
              deallocate(qc_H08)

            endif ! [nprof >= 1]


#endif
!    -- End of Him8 DA part --

!   -- TC vital DA -- 
#ifdef TCV
        if(minval(TC_rij) > 0.0d0)then
          !
          ! TC vital obs should have 3 data (i.e., lon, lat, and MSLP)
          ! bTC(1,:) : tcx (m), bTC(2,:): tcy (m), bTC(3,:): mslp

          ! (1) Get "local (each subdomain's)" TC information 
          !
          !     *A local TC simply corresponds to the minimum SLP in each
          !     subdomain. 
          !     *If the diagnosed TC position in a subdomain is too far
          !     from the vital position, the subdomain is regarded as a no TC
          !     subdomain.

          allocate(bTC(3,0:MEM_NP-1))

          bTC = 9.99d33
          call search_tc_subdom(TC_rij(1),TC_rij(2),v2dg,bTC(1,myrank_d),bTC(2,myrank_d),bTC(3,myrank_d))
          ! Note: obs(iof)%dat(obs_idx_TCX) is not longitude (deg) but X (m).
          !       Units of the original TC vital position are converted in
          !       subroutine read_obs in common_obs_scale.f90.
          !

          ! (2) Compare local TCs in each subdomain and assign the strongest one as a background TC in each member
          !
          if (nprocs_d > 1) then
             CALL MPI_ALLREDUCE(MPI_IN_PLACE,bTC,3*MEM_NP,MPI_r_size,MPI_MIN,MPI_COMM_d,ierr)
          end if

          bTC_mslp = 9.99d33
          do n = 0, MEM_NP - 1
            if (bTC(3,n) < bTC_mslp ) then
              bTC_mslp = bTC(3,n)
              bTC_rank_d = n
            endif
          enddo ! [ n = 0, MEM_NP - 1]

          ! (3) Substitute Hx in a subdomain that covers an observed TC
          !
          if (obs_nn_TCP > 0) then
            do n = 1, 3
              if (n == 1) nn = obs_nn_TCX
              if (n == 2) nn = obs_nn_TCY
              if (n == 3) nn = obs_nn_TCP
              obsda%val(nn) = bTC(n,btc_rank_d)
              obsda%qc(nn) = iqc_good
            enddo ! [ n = 1, 3 ]
          endif ! [ obs_nn_TCP > 0]
          deallocate(bTC)

        endif ! [minval(TC_rij) > 0.0d0]
#endif
!   -- End of TC vital DA -- 

        write (timer_str, '(A30,I4,A7,I4,A2)') 'obsope_cal:obsope_step_2   (t=', it, ', slot=', islot, '):'
        call mpi_timer(trim(timer_str), 2)

      end do ! [ islot = SLOT_START, SLOT_END ]

      call mpi_timer('', 2)

      ! Write obsda data to files if OBSDA_OUT = .true.
      ! 
      if (OBSDA_OUT) then
!        write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is writing observations for member ', &
!              im, ', subdomain id #', proc2mem(2,it,myrank+1)
        if (im <= MEMBER) then
          call file_member_replace(im, OBSDA_OUT_BASENAME, obsdafile)
        else if (im == mmean) then
          obsdafile = OBSDA_MEAN_OUT_BASENAME
        else if (im == mmdet) then
          obsdafile = OBSDA_MDET_OUT_BASENAME
        end if
        write (obsda_suffix(2:7),'(I6.6)') proc2mem(2,it,myrank+1)
        write (6,'(A,I6.6,2A)') 'MYRANK ', myrank,' is writing an obsda file ', trim(obsdafile)//obsda_suffix
        call write_obs_da(trim(obsdafile)//obsda_suffix,obsda,0)

        write (timer_str, '(A30,I4,A2)') 'obsope_cal:write_obs_da    (t=', it, '):'
        call mpi_timer(trim(timer_str), 2)
      end if

      ! Prepare variables that will need to be communicated if obsda_return = .true.
      ! 
      if (obsda_return .and. nobs > 0) then
        ! variables with an ensemble dimension
        if (im == mmdetin) then
          obsda%ensval(mmdetobs, 1:nobs) = obsda%val(1:nobs)
        else
          obsda%ensval(im, 1:nobs) = obsda%val(1:nobs)
        end if

        ! variables without an ensemble dimension: merge them with the previous iterations
        if (nitmax > 1 .and. nobs > 0) then
          if (it == 1) then
            qc_p(:) = obsda%qc(1:nobs)
#ifdef H08
            lev_p(:) = obsda%lev(1:nobs)
            val2_p(:) = obsda%val2(1:nobs)
#endif
          else
            qc_p(:) = max(qc_p(:), obsda%qc(1:nobs))
#ifdef H08
            if (im <= MEMBER) then ! only consider lev from members, not from the means
              lev_p(:) = lev_p(:) + obsda%lev(1:nobs)
              val2_p(:) = val2_p(:) + obsda%val2(1:nobs)
            end if
#endif
          end if
        end if
      end if ! [ obsda_return .and. nobs > 0 ]

    end if ! [ (im >= 1 .and. im <= MEMBER) .or. im == mmdetin ]
  end do ! [ it = 1, nitmax ]

  deallocate ( v3dg, v2dg )
  deallocate ( bsn, bsna )

!-------------------------------------------------------------------------------
! If obsda_return = .true., allreduce the obsda data among ensemble members
!-------------------------------------------------------------------------------

  call mpi_timer('', 2, barrier=MPI_COMM_e)

  if (obsda_return .and. nobs > 0) then
    ! variables with an ensemble dimension
    if (nprocs_e > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, obsda%ensval(:,1:nobs), nensobs*nobs, MPI_r_size, MPI_SUM, MPI_COMM_e, ierr)
    end if

    ! variables without an ensemble dimension
    if (nitmax > 1) then
      obsda%qc(1:nobs) = qc_p(:)
      deallocate (qc_p)
#ifdef H08
      obsda%lev(1:nobs) = lev_p(:)
      deallocate (lev_p)
      obsda%val2(1:nobs) = val2_p(:)
      deallocate (val2_p)
#endif
    end if

    if (nprocs_e > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, obsda%qc(1:nobs), nobs, MPI_INTEGER, MPI_MAX, MPI_COMM_e, ierr)  ! maximum value of qc
    end if
#ifdef H08
    if (nprocs_e > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, obsda%lev(1:nobs), nobs, MPI_r_size, MPI_SUM, MPI_COMM_e, ierr)  ! ensemble mean of obsda%lev
      call MPI_ALLREDUCE(MPI_IN_PLACE, obsda%val2(1:nobs), nobs, MPI_r_size, MPI_SUM, MPI_COMM_e, ierr)  ! ensemble mean of obsda%val2
    end if
    obsda%lev(1:nobs) = obsda%lev(1:nobs) / REAL(MEMBER, r_size)                                                    !
    obsda%val2(1:nobs) = obsda%val2(1:nobs) / REAL(MEMBER, r_size)                                                    !
#endif

    call mpi_timer('obsope_cal:mpi_allreduce:', 2)
  end if ! [ obsda_return .and. nobs > 0 ]

  if (.not. obsda_return) then
    call obs_da_value_deallocate(obsda)
  end if

  return
end subroutine obsope_cal

!-----------------------------------------------------------------------
! Observation generator calculation
!-----------------------------------------------------------------------
SUBROUTINE obsmake_cal(obs)
  IMPLICIT NONE

  TYPE(obs_info),INTENT(INOUT) :: obs(OBS_IN_NUM)
  REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)

  integer :: islot,proc
  integer :: n,nslot,nobs,nobs_slot,ierr,iqc,iof
  integer :: nobsmax,nobsall
  real(r_size) :: rig,rjg,ri,rj,rk
  real(r_size) :: slot_lb,slot_ub
  real(r_size),allocatable :: bufr(:)
  real(r_size),allocatable :: error(:)

  CHARACTER(10) :: obsoutfile = 'obsout.dat'
  INTEGER :: ns 
#ifdef H08
! obsmake for H08 is not available !! (03/17/2016) T.Honda
! -- for Himawari-8 obs --
  INTEGER :: nallprof ! H08: Num of all profiles (entire domain) required by RTTOV
  INTEGER :: nprof_H08 ! num of H08 obs
  REAL(r_size),ALLOCATABLE :: ri_H08(:),rj_H08(:)
  REAL(r_size),ALLOCATABLE :: lon_H08(:),lat_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_ri_H08(:),tmp_rj_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_lon_H08(:),tmp_lat_H08(:)

  REAL(r_size),ALLOCATABLE :: yobs_H08(:),plev_obs_H08(:)
  INTEGER,ALLOCATABLE :: qc_H08(:)
  INTEGER,ALLOCATABLE :: idx_H08(:) ! index array
  INTEGER :: ich
#endif

!-----------------------------------------------------------------------

  write (6,'(A,I6.6,A,I6.6)') 'MYRANK ', myrank, ' is processing subdomain id #', myrank_d

  allocate ( v3dg (nlevh,nlonh,nlath,nv3dd) )
  allocate ( v2dg (nlonh,nlath,nv2dd) )

  do iof = 1, OBS_IN_NUM
    obs(iof)%dat = 0.0d0
  end do

  nobs = 0
  do islot = SLOT_START, SLOT_END
    slot_lb = (real(islot-SLOT_BASE,r_size) - 0.5d0) * SLOT_TINTERVAL
    slot_ub = (real(islot-SLOT_BASE,r_size) + 0.5d0) * SLOT_TINTERVAL
    write (6,'(A,I3,A,F9.1,A,F9.1,A)') 'Slot #', islot-SLOT_START+1, ': time window (', slot_lb, ',', slot_ub, '] sec'

    call read_ens_history_iter(1,islot,v3dg,v2dg)

    do iof = 1, OBS_IN_NUM
      IF(OBS_IN_FORMAT(iof) /= 3)THEN ! except H08 obs
        nslot = 0
        nobs_slot = 0
        do n = 1, obs(iof)%nobs

          if (obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub) then
            nslot = nslot + 1

            call phys2ij(obs(iof)%lon(n),obs(iof)%lat(n),rig,rjg)
            call rij_g2l_auto(proc,rig,rjg,ri,rj)

  !          if (myrank_d == 0) then
  !            print *, proc, rig, rjg, ri, rj
  !          end if

            if (proc < 0 .and. myrank_d == 0) then ! if outside of the global domain, processed by myrank_d == 0
              obs(iof)%dat(n) = undef
            end if

            if (myrank_d == proc) then
              nobs = nobs + 1
              nobs_slot = nobs_slot + 1

  !IF(NINT(elem(n)) == id_ps_obs) THEN
  !  CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,dz)
  !  rk = rlev(n) - dz
  !  IF(ABS(rk) > threshold_dz) THEN ! pressure adjustment threshold
  !    ! WRITE(6,'(A)') '* PS obs vertical adjustment beyond threshold'
  !    ! WRITE(6,'(A,F10.2,A,F6.2,A,F6.2,A)') '* dz=',rk,&
  !    ! & ', (lon,lat)=(',elon(n),',',elat(n),')'
  !    CYCLE
  !  END IF
  !END IF

              if (obs(iof)%elm(n) == id_radar_ref_obs .or. obs(iof)%elm(n) == id_radar_ref_zero_obs .or. obs(iof)%elm(n) == id_radar_vr_obs) then
                call phys2ijkz(v3dg(:,:,:,iv3dd_hgt),ri,rj,obs(iof)%lev(n),rk,iqc)
              else
                call phys2ijk(v3dg(:,:,:,iv3dd_p),obs(iof)%elm(n),ri,rj,obs(iof)%lev(n),rk,iqc)
              end if

              if (iqc /= iqc_good) then
                obs(iof)%dat(n) = undef
              else
                select case (OBS_IN_FORMAT(iof))
                case (1)
                  call Trans_XtoY(obs(iof)%elm(n),ri,rj,rk, &
                                  obs(iof)%lon(n),obs(iof)%lat(n),v3dg,v2dg,obs(iof)%dat(n),iqc)
                case (2)
                  call Trans_XtoY_radar(obs(iof)%elm(n),obs(iof)%meta(1),obs(iof)%meta(2),obs(iof)%meta(3),ri,rj,rk, &
                                        obs(iof)%lon(n),obs(iof)%lat(n),obs(iof)%lev(n),v3dg,v2dg,obs(iof)%dat(n),iqc)
                end select

 !!! For radar observation, when reflectivity value is too low, do not generate ref/vr observations
 !!! No consideration of the terrain blocking effects.....

                if (iqc /= iqc_good) then
                  obs(iof)%dat(n) = undef
                end if
              end if

            end if ! [ myrank_d == proc ]

          end if ! [ obs%dif(n) > slot_lb .and. obs%dif(n) <= slot_ub ]

        end do ! [ n = 1, obs%nobs ]

#ifdef H08
! -- H08 part --
      ELSEIF(OBS_IN_FORMAT(iof) == 3)THEN ! H08
        nslot = 0
        nobs_slot = 0
        nprof_H08 = 0

        nallprof = obs(iof)%nobs/nch

        ALLOCATE(tmp_ri_H08(nallprof))
        ALLOCATE(tmp_rj_H08(nallprof))
        ALLOCATE(tmp_lon_H08(nallprof))
        ALLOCATE(tmp_lat_H08(nallprof))
        ALLOCATE(idx_H08(nallprof))

        do n = 1, nallprof
          ns = (n - 1) * nch + 1
          if (obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub) then
            nslot = nslot + 1

            call phys2ij(obs(iof)%lon(ns),obs(iof)%lat(ns),rig,rjg)
            call rij_g2l_auto(proc,rig,rjg,ri,rj)


            if (proc < 0 .and. myrank_d == 0) then ! if outside of the global domain, processed by myrank_d == 0
              obs(iof)%dat(ns:ns+nch-1) = undef
            end if

            if (myrank_d == proc) then
              nprof_H08 = nprof_H08 + 1 ! num of prof in myrank node
              idx_H08(nprof_H08) = ns ! idx of prof in myrank node
              tmp_ri_H08(nprof_H08) = ri
              tmp_rj_H08(nprof_H08) = rj
              tmp_lon_H08(nprof_H08) = obs(iof)%lon(ns)
              tmp_lat_H08(nprof_H08) = obs(iof)%lat(ns)

              nobs = nobs + nch
              nobs_slot = nobs_slot + nch

            end if ! [ myrank_d == proc ]

          end if ! [ obs%dif(n) > slot_lb .and. obs%dif(n) <= slot_ub ]

        end do ! [ n = 1, nallprof ]

        IF(nprof_H08 >=1)THEN
          ALLOCATE(ri_H08(nprof_H08))
          ALLOCATE(rj_H08(nprof_H08))
          ALLOCATE(lon_H08(nprof_H08))
          ALLOCATE(lat_H08(nprof_H08))

          ri_H08 = tmp_ri_H08(1:nprof_H08)
          rj_H08 = tmp_rj_H08(1:nprof_H08)
          lon_H08 = tmp_lon_H08(1:nprof_H08)
          lat_H08 = tmp_lat_H08(1:nprof_H08)

          ALLOCATE(yobs_H08(nprof_H08*nch))
          ALLOCATE(plev_obs_H08(nprof_H08*nch))
          ALLOCATE(qc_H08(nprof_H08*nch))

          CALL Trans_XtoY_H08(nprof_H08,ri_H08,rj_H08,&
                              lon_H08,lat_H08,v3dg,v2dg,&
                              yobs_H08,plev_obs_H08,&
                              qc_H08)

          DO n = 1, nprof_H08
            ns = idx_H08(n)

            obs(iof)%lon(ns:ns+nch-1)=lon_H08(n:n+nch-1)
            obs(iof)%lat(ns:ns+nch-1)=lat_H08(n:n+nch-1)

            DO ich = 1, nch-1
              IF(qc_H08(n+ich-1) == iqc_good)THEN
                obs(iof)%dat(ns+ich-1)=undef
              ELSE
                obs(iof)%dat(ns+ich-1)=yobs_H08(n+ich-1)
              ENDIF
            ENDDO
          ENDDO

        ENDIF

        DEALLOCATE(tmp_ri_H08,tmp_rj_H08)
        DEALLOCATE(tmp_lon_H08,tmp_lat_H08)


! -- end of H08 part --
#endif
      ENDIF

    end do ! [ iof = 1, OBS_IN_NUM ]

    write (6,'(3A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot = ', nslot
    write (6,'(3A,I6,A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot and processed by rank ', myrank, ' = ', nobs_slot

  end do ! [ islot = SLOT_START, SLOT_END ]

  deallocate ( v3dg, v2dg )

  if (myrank_d == 0) then
    nobsmax = 0
    nobsall = 0
    do iof = 1, OBS_IN_NUM
      if (obs(iof)%nobs > nobsmax) nobsmax = obs(iof)%nobs
      nobsall = nobsall + obs(iof)%nobs
    end do

    allocate ( bufr(nobsmax) )
    allocate ( error(nobsall) )

    call com_randn(nobsall, error) ! generate all random numbers at the same time
    ns = 0
  end if

  do iof = 1, OBS_IN_NUM

    call MPI_REDUCE(obs(iof)%dat,bufr(1:obs(iof)%nobs),obs(iof)%nobs,MPI_r_size,MPI_SUM,0,MPI_COMM_d,ierr)

    if (myrank_d == 0) then
      obs(iof)%dat = bufr(1:obs(iof)%nobs)

      do n = 1, obs(iof)%nobs
        select case(obs(iof)%elm(n))
        case(id_u_obs)
          obs(iof)%err(n) = OBSERR_U
        case(id_v_obs)
          obs(iof)%err(n) = OBSERR_V
        case(id_t_obs,id_tv_obs)
          obs(iof)%err(n) = OBSERR_T
        case(id_q_obs)
          obs(iof)%err(n) = OBSERR_Q
        case(id_rh_obs)
          obs(iof)%err(n) = OBSERR_RH
        case(id_ps_obs)
          obs(iof)%err(n) = OBSERR_PS
        case(id_radar_ref_obs,id_radar_ref_zero_obs)
          obs(iof)%err(n) = OBSERR_RADAR_REF
        case(id_radar_vr_obs)
          obs(iof)%err(n) = OBSERR_RADAR_VR
!
! -- Not available (02/09/2015)
!        case(id_H08IR_obs) ! H08
!          obs(iof)%err(n) = OBSERR_H08(ch) !H08
!        case default
          write(6,'(A)') 'warning: skip assigning observation error (unsupported observation type)' 
        end select

        if (obs(iof)%dat(n) /= undef .and. obs(iof)%err(n) /= undef) then
          obs(iof)%dat(n) = obs(iof)%dat(n) + obs(iof)%err(n) * error(ns+n)
        end if

!print *, '######', obs%elm(n), obs%dat(n)
      end do ! [ n = 1, obs(iof)%nobs ]

      ns = ns + obs(iof)%nobs
    end if ! [ myrank_d == 0 ]

  end do ! [ iof = 1, OBS_IN_NUM ]

  if (myrank_d == 0) then
    deallocate ( bufr )
    deallocate ( error )

    call write_obs_all(obs, missing=.false., file_suffix='.out') ! only at the head node
  end if

end subroutine obsmake_cal

!-------------------------------------------------------------------------------
! Model-to-observation simulator calculation
!-------------------------------------------------------------------------------
subroutine obssim_cal(v3dgh, v2dgh, v3dgsim, v2dgsim, stggrd)
  use scale_grid, only: &
      GRID_CX, GRID_CY, &
      DX, DY
  use scale_grid_index, only: &
      IHALO, JHALO, KHALO
  use scale_mapproj, only: &
      MPRJ_xy2lonlat

  implicit none

  real(r_size), intent(in) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  real(r_size), intent(in) :: v2dgh(nlonh,nlath,nv2dd)
  real(r_size), intent(out) :: v3dgsim(nlev,nlon,nlat,OBSSIM_NUM_3D_VARS)
  real(r_size), intent(out) :: v2dgsim(nlon,nlat,OBSSIM_NUM_2D_VARS)
  integer, intent(in), optional :: stggrd

  integer :: i, j, k, iv3dsim, iv2dsim
  real(r_size) :: ri, rj, rk
  real(r_size) :: lon, lat, lev
  real(r_size) :: tmpobs
  integer :: tmpqc

#ifdef H08
! -- for Himawari-8 obs --
  real(r_size) :: ri_H08(nlon*nlat),rj_H08(nlon*nlat)
  real(r_size) :: lon_H08(nlon*nlat),lat_H08(nlon*nlat)
  integer :: np ! num of Him8 profile
  real(r_size) :: yobs_H08(nch,nlon*nlat),plev_obs_H08(nch,nlon*nlat)
  integer :: qc_H08(nch,nlon*nlat)
  integer :: ch

  np = 0

#endif

!-------------------------------------------------------------------------------

  write (6,'(A,I6.6,A,I6.6)') 'MYRANK ', myrank, ' is processing subdomain id #', myrank_d


  do j = 1, nlat
    rj = real(j + JHALO, r_size)

    do i = 1, nlon
      ri = real(i + IHALO, r_size)
      call MPRJ_xy2lonlat((ri-1.0d0) * DX + GRID_CX(1), (rj-1.0d0) * DY + GRID_CY(1), lon, lat)
      lon = lon * rad2deg
      lat = lat * rad2deg

#ifdef H08
      np = np + 1

      ri_H08(np) = ri
      rj_H08(np) = rj
      lon_H08(np) = lon
      lat_H08(np) = lat
#endif

      do k = 1, nlev
        rk = real(k + KHALO, r_size)

        do iv3dsim = 1, OBSSIM_NUM_3D_VARS
          select case (OBSSIM_3D_VARS_LIST(iv3dsim))
          case (id_radar_ref_obs, id_radar_ref_zero_obs, id_radar_vr_obs, id_radar_prh_obs)
            lev = v3dgh(k+KHALO, i+IHALO, j+JHALO, iv3dd_hgt)
            call Trans_XtoY_radar(OBSSIM_3D_VARS_LIST(iv3dsim), OBSSIM_RADAR_LON, OBSSIM_RADAR_LAT, OBSSIM_RADAR_Z, ri, rj, rk, &
                                  lon, lat, lev, v3dgh, v2dgh, tmpobs, tmpqc, stggrd)
            if (tmpqc == iqc_ref_low) tmpqc = iqc_good ! when process the observation operator, we don't care if reflectivity is too small
          case default
            call Trans_XtoY(OBSSIM_3D_VARS_LIST(iv3dsim), ri, rj, rk, &
                            lon, lat, v3dgh, v2dgh, tmpobs, tmpqc, stggrd)
          end select

          if (tmpqc == 0) then
            v3dgsim(k,i,j,iv3dsim) = real(tmpobs, r_sngl)
          else
            v3dgsim(k,i,j,iv3dsim) = real(undef, r_sngl)
          end if
        end do ! [ iv3dsim = 1, OBSSIM_NUM_3D_VARS ]

        ! 2D observations calculated when k = 1
        if (k == 1) then
          do iv2dsim = 1, OBSSIM_NUM_2D_VARS
            select case (OBSSIM_2D_VARS_LIST(iv2dsim))
            case (id_H08IR_obs)   
              cycle
!            case (id_tclon_obs, id_tclat_obs, id_tcmip_obs)
!              call ...
            case default
              call Trans_XtoY(OBSSIM_2D_VARS_LIST(iv2dsim), ri, rj, rk, &
                              lon, lat, v3dgh, v2dgh, tmpobs, tmpqc, stggrd)
            end select

            if (tmpqc == 0) then
              v2dgsim(i,j,iv2dsim) = real(tmpobs, r_sngl)
            else
              v2dgsim(i,j,iv2dsim) = real(undef, r_sngl)
            end if
          end do ! [ iv2dsim = 1, OBSSIM_NUM_2D_VARS ]
        end if ! [ k == 1 ]

      end do ! [ k = 1, nlev ]

    end do ! [ i = 1, nlon ]

  end do ! [ j = 1, nlat ]


#ifdef H08
  call Trans_XtoY_H08(np,ri_H08(1:np),rj_H08(1:np),&
                      lon_H08(1:np),lat_H08(1:np),v3dgh,v2dgh,&
                      yobs_H08,plev_obs_H08,&
                      qc_H08)

  np = 0
  do j = 1, nlat

    do i = 1, nlon

      np = np + 1
      ch = 0

      do iv2dsim = 1, OBSSIM_NUM_2D_VARS
        if(OBSSIM_2D_VARS_LIST(iv2dsim) /= id_H08IR_obs .and. ch <= nch)cycle

        ch = ch + 1
        v2dgsim(i,j,iv2dsim) = real(yobs_H08(ch,np), r_sngl)
      enddo ! [iv2dsim = 1, OBSSIM_NUM_2D_VARS]

    end do ! [ i = 1, nlon ]

  end do ! [ j = 1, nlat ]


#endif

!-------------------------------------------------------------------------------

end subroutine obssim_cal

!!!!!! it is not good to open/close a file many times for different steps !!!!!!
!-------------------------------------------------------------------------------
! Write the subdomain model data into a single GrADS file
!-------------------------------------------------------------------------------
subroutine write_grd_mpi(filename, nv3dgrd, nv2dgrd, step, v3d, v2d)
  implicit none
  character(*), intent(in) :: filename
  integer, intent(in) :: nv3dgrd
  integer, intent(in) :: nv2dgrd
  integer, intent(in) :: step
  real(r_size), intent(in) :: v3d(nlev,nlon,nlat,nv3dgrd)
  real(r_size), intent(in) :: v2d(nlon,nlat,nv2dgrd)

  real(r_sngl) :: bufs4(nlong,nlatg)
  real(r_sngl) :: bufr4(nlong,nlatg)
  integer :: iunit, iolen
  integer :: k, n, irec, ierr
  integer :: proc_i, proc_j
  integer :: ishift, jshift

  call rank_1d_2d(myrank_d, proc_i, proc_j)
  ishift = proc_i * nlon
  jshift = proc_j * nlat

  if (myrank_d == 0) then
    iunit = 55
    inquire (iolength=iolen) iolen
    open (iunit, file=trim(filename), form='unformatted', access='direct', &
          status='unknown', convert='native', recl=nlong*nlatg*iolen)
    irec = (nlev * nv3dgrd + nv2dgrd) * (step-1)
  end if

  do n = 1, nv3dgrd
    do k = 1, nlev
      bufs4(:,:) = 0.0
      bufs4(1+ishift:nlon+ishift, 1+jshift:nlat+jshift) = real(v3d(k,:,:,n), r_sngl)
      call MPI_REDUCE(bufs4, bufr4, nlong*nlatg, MPI_REAL, MPI_SUM, 0, MPI_COMM_d, ierr)
      if (myrank_d == 0) then
        irec = irec + 1
        write (iunit, rec=irec) bufr4
      end if
    end do
  end do

  do n = 1, nv2dgrd
    bufs4(:,:) = 0.0
    bufs4(1+ishift:nlon+ishift, 1+jshift:nlat+jshift) = real(v2d(:,:,n), r_sngl)
    call MPI_REDUCE(bufs4, bufr4, nlong*nlatg, MPI_REAL, MPI_SUM, 0, MPI_COMM_d, ierr)
    if (myrank_d == 0) then
      irec = irec + 1
      write (iunit, rec=irec) bufr4
    end if
  end do

  if (myrank_d == 0) then
    close (iunit)
  end if

  return
end subroutine write_grd_mpi

!=======================================================================

END MODULE obsope_tools
