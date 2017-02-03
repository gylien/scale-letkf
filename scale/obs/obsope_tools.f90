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

!  use common_scalelib

  use common_nml

!  use scale_process, only: &
!       PRC_myrank
!       MPI_COMM_d => LOCAL_COMM_WORLD

  use scale_grid_index, only: &
    KHALO, IHALO, JHALO
#ifdef H08
  use scale_grid, only: &
      DX, DY,           &
      BUFFER_DX,        &
      BUFFER_DY
#endif


  IMPLICIT NONE
  PUBLIC

CONTAINS

!-----------------------------------------------------------------------
! Observation operator calculation
!-----------------------------------------------------------------------
SUBROUTINE obsope_cal(obsda_return, nobs_extern)
  IMPLICIT NONE

  type(obs_da_value),optional,intent(inout) :: obsda_return
  integer,optional,intent(in) :: nobs_extern
  type(obs_da_value) :: obsda
  REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)

  integer :: it,islot,proc,im,iof
  integer :: n,nn,nslot,nobs,nobs_0,nobs_slot,nobs_alldomain


  integer :: nobs_max_per_file, nobs_max_per_file_sub
  integer :: nn_0, nsub, nmod, n1, n2
  integer :: ip, ibufs
  integer, allocatable :: cntr(:), dspr(:)
  integer, allocatable :: cnts(:), dsps(:)
  integer, allocatable :: bsn(:,:), bsna(:,:), bsnext(:,:)


!  real(r_size) :: rig,rjg,ri,rj,rk
  real(r_size) :: rig,rjg,rk
  real(r_size),allocatable :: ri(:),rj(:)


!  real(r_size),allocatable :: obri(:),obrj(:)
!  integer,allocatable :: obslot(:)
  integer,allocatable :: obrank(:)
  real(r_size),allocatable :: obri(:)
  real(r_size),allocatable :: obrj(:)
  integer,allocatable :: obrank_bufs(:)
  real(r_size),allocatable :: ri_bufs(:)
  real(r_size),allocatable :: rj_bufs(:)

  integer,allocatable :: obset_bufs(:)
  integer,allocatable :: obidx_bufs(:)
  real(r_size),allocatable :: ri_bufs2(:)
  real(r_size),allocatable :: rj_bufs2(:)


  real(r_size) :: ril, rjl
  real(r_size) :: slot_lb,slot_ub

#ifdef H08
! -- for Himawari-8 obs --
  INTEGER :: nallprof ! H08: Num of all profiles (entire domain) required by RTTOV
  INTEGER :: ns ! H08 obs count
  INTEGER :: nprof_H08 ! num of H08 obs
  REAL(r_size),ALLOCATABLE :: ri_H08(:),rj_H08(:)
  REAL(r_size),ALLOCATABLE :: lon_H08(:),lat_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_ri_H08(:),tmp_rj_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_lon_H08(:),tmp_lat_H08(:)

  REAL(r_size),ALLOCATABLE :: yobs_H08(:),plev_obs_H08(:)
  REAL(r_size),ALLOCATABLE :: yobs_H08_clr(:)
  INTEGER :: ch
  INTEGER,ALLOCATABLE :: qc_H08(:)

! -- Rejecting obs over the buffer regions. --
!
! bris: "ri" at the wetern end of the domain excluding buffer regions
! brie: "ri" at the eastern end of the domain excluding buffer regions
! bris: "rj" at the southern end of the domain excluding buffer regions
! bris: "rj" at the northern end of the domain excluding buffer regions
!
! e.g.,   ri:    ...bris...........brie...
!             buffer |  NOT buffer  | buffer
!
!
  REAL(r_size) :: bris, brie
  REAL(r_size) :: brjs, brje
#endif

! -- for TC vital assimilation --
  INTEGER :: obs_set_TCX, obs_set_TCY, obs_set_TCP ! obs set
  INTEGER :: obs_idx_TCX, obs_idx_TCY, obs_idx_TCP ! obs index
  INTEGER :: bTC_proc ! the process where the background TC is located.
! bTC: background TC in each subdomain
! bTC(1,:) : tcx (m), bTC(2,:): tcy (m), bTC(3,:): mslp (Pa)
  REAL(r_size),ALLOCATABLE :: bTC(:,:)
  REAL(r_size) :: bTC_mslp

  character(filelenmax) :: obsdafile
  character(11) :: obsda_suffix = '.000000.dat'

!-----------------------------------------------------------------------


  integer :: ierr
  real(r_dble) :: rrtimer00, rrtimer
!  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


#ifdef H08
!  call phys2ij(MSLP_TC_LON,MSLP_TC_LAT,MSLP_TC_rig,MSLP_TC_rjg)
  bris = real(BUFFER_DX/DX,r_size) + real(IHALO,r_size) 
  brjs = real(BUFFER_DY/DY,r_size) + real(JHALO,r_size)
  brie = (real(nlong+2*IHALO,r_size) - bris)
  brje = (real(nlatg+2*JHALO,r_size) - brjs)
#endif

  nobs_alldomain = 0
  nobs_max_per_file = 0
  do iof = 1, OBS_IN_NUM
    if (OBSDA_RUN(iof)) then
      nobs_alldomain = nobs_alldomain + obs(iof)%nobs
      if (obs(iof)%nobs > nobs_max_per_file) then
        nobs_max_per_file = obs(iof)%nobs
      end if
    end if
  end do
!!!  obsda%nobs = nobs_alldomain
!!!  call obs_da_value_allocate(obsda,0)
!!!  allocate ( ri(nobs_alldomain) )
!!!  allocate ( rj(nobs_alldomain) )



  allocate ( v3dg (nlevh,nlonh,nlath,nv3dd) )
  allocate ( v2dg (nlonh,nlath,nv2dd) )






  obs_set_TCX = -1
  obs_set_TCY = -1
  obs_set_TCP = -1
  obs_idx_TCX = -1
  obs_idx_TCY = -1
  obs_idx_TCP = -1

  allocate (cntr(nprocs_a))
  allocate (dspr(nprocs_a))

  allocate ( obrank(nobs_alldomain) )
  allocate ( obri(nobs_alldomain) )
  allocate ( obrj(nobs_alldomain) )

  nobs_max_per_file_sub = (nobs_max_per_file - 1) / nprocs_a + 1
  allocate ( obrank_bufs(nobs_max_per_file_sub) )
  allocate ( ri_bufs(nobs_max_per_file_sub) )
  allocate ( rj_bufs(nobs_max_per_file_sub) )


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
        select case (obs(iof)%elm(n))
        case (id_tclon_obs)
          obs_set_TCX = iof
          obs_idx_TCX = n
          cycle
        case (id_tclat_obs)
          obs_set_TCY = iof
          obs_idx_TCY = n
          cycle
        case (id_tcmip_obs)
          obs_set_TCP = iof
          obs_idx_TCP = n
          cycle
        end select

        islot = ceiling(obs(iof)%dif(n) / SLOT_TINTERVAL - 0.5d0) + SLOT_BASE

        if (islot >= SLOT_START .and. islot <= SLOT_END) then
          call phys2ij(obs(iof)%lon(n), obs(iof)%lat(n), ri_bufs(ibufs), rj_bufs(ibufs))
          call rij_g2l_auto(obrank_bufs(ibufs), ri_bufs(ibufs), rj_bufs(ibufs), ril, rjl) ! rij, rjl discarded here; re-computed later
        end if
      end do ! [ n = n1, n2 ]


  rrtimer = MPI_WTIME()
  write (6,'(A,4x,F15.7)') '###### obsope_cal:first_scan_cal:               ', rrtimer-rrtimer00
  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


      call MPI_GATHERV(obrank_bufs, cntr(myrank_a+1), MPI_INTEGER, obrank, cntr, dspr, MPI_INTEGER, 0, MPI_COMM_a, ierr)
      call MPI_GATHERV(ri_bufs,     cntr(myrank_a+1), MPI_r_size , obri,   cntr, dspr, MPI_r_size,  0, MPI_COMM_a, ierr)
      call MPI_GATHERV(rj_bufs,     cntr(myrank_a+1), MPI_r_size , obrj,   cntr, dspr, MPI_r_size,  0, MPI_COMM_a, ierr)



  rrtimer = MPI_WTIME()
  write (6,'(A,4x,F15.7)') '###### obsope_cal:first_scan_reduce:            ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


      nn_0 = nn_0 + obs(iof)%nobs
    end if ! [ OBSDA_RUN(iof) .and. obs(iof)%nobs > 0 ]
  end do ! [ do iof = 1, OBS_IN_NUM ]

  deallocate (cntr, dspr)
  deallocate (obrank_bufs, ri_bufs, rj_bufs)

  allocate (bsna  (SLOT_START-1:SLOT_END, 0:nprocs_d-1))

  if (myrank_e == 0) then
    allocate ( obset_bufs(nobs_alldomain) )
    allocate ( obidx_bufs(nobs_alldomain) )
    allocate ( ri_bufs2(nobs_alldomain) )
    allocate ( rj_bufs2(nobs_alldomain) )
  end if

  if (myrank_a == 0) then
    allocate (bsn   (SLOT_START  :SLOT_END, 0:nprocs_d-1))
    allocate (bsnext(SLOT_START  :SLOT_END, 0:nprocs_d-1))
    bsn(:,:) = 0
    bsna(:,:) = 0
    bsnext(:,:) = 0

    nn = 0
    do iof = 1, OBS_IN_NUM
      if (OBSDA_RUN(iof) .and. obs(iof)%nobs > 0) then
        do n = 1, obs(iof)%nobs
          nn = nn + 1
          islot = ceiling(obs(iof)%dif(n) / SLOT_TINTERVAL - 0.5d0) + SLOT_BASE
          if (islot >= SLOT_START .and. islot <= SLOT_END) then
!            nslot = nslot + 1
            if (obrank(nn) /= -1) then
!              nobs = nobs + 1
!              nobs_slot = nobs_slot + 1

              bsn(islot, obrank(nn)) = bsn(islot, obrank(nn)) + 1
            end if
          end if
        end do ! [ n = 1, obs(iof)%nobs ]
      end if ! [ OBSDA_RUN(iof) .and. obs(iof)%nobs > 0 ]
    end do ! [ do iof = 1, OBS_IN_NUM ]

    do ip = 0, nprocs_d-1
      if (ip > 0) then
        bsna(SLOT_START-1, ip) = bsna(SLOT_END, ip-1)
      end if
      do islot = SLOT_START, SLOT_END
        bsna(islot, ip) = bsna(islot-1, ip) + bsn(islot, ip)
      end do
      bsnext(SLOT_START:SLOT_END, ip) = bsna(SLOT_START-1:SLOT_END-1, ip)
    end do

    nn = 0
    do iof = 1, OBS_IN_NUM
      if (OBSDA_RUN(iof) .and. obs(iof)%nobs > 0) then
        do n = 1, obs(iof)%nobs
          nn = nn + 1
          islot = ceiling(obs(iof)%dif(n) / SLOT_TINTERVAL - 0.5d0) + SLOT_BASE
          if (islot >= SLOT_START .and. islot <= SLOT_END) then
!            nslot = nslot + 1
            if (obrank(nn) /= -1) then
!              nobs = nobs + 1
!              nobs_slot = nobs_slot + 1

              bsnext(islot, obrank(nn)) = bsnext(islot, obrank(nn)) + 1
              obset_bufs(bsnext(islot, obrank(nn))) = iof
              obidx_bufs(bsnext(islot, obrank(nn))) = n
              ri_bufs2(bsnext(islot, obrank(nn))) = obri(nn)
              rj_bufs2(bsnext(islot, obrank(nn))) = obrj(nn)
            end if
          end if
        end do ! [ n = 1, obs(iof)%nobs ]
      end if ! [ OBSDA_RUN(iof) .and. obs(iof)%nobs > 0 ]
    end do ! [ do iof = 1, OBS_IN_NUM ]

    deallocate (bsn, bsnext)


  rrtimer = MPI_WTIME()
  write (6,'(A,4x,F15.7)') '###### obsope_cal:bucket_sort:                  ', rrtimer-rrtimer00


  end if ! [ myrank_a == 0 ]

  deallocate ( obrank, obri, obrj )


  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


  call MPI_BCAST(bsna, (SLOT_END-SLOT_START+2)*nprocs_d, MPI_INTEGER, 0, MPI_COMM_a, ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,4x,F15.7)') '###### obsope_cal:sort_info_bcast:              ', rrtimer-rrtimer00
  rrtimer00 = rrtimer


  obsda%nobs = bsna(SLOT_END, myrank_d) - bsna(SLOT_START-1, myrank_d)
  call obs_da_value_allocate(obsda, 0)

  if (myrank_e == 0) then
    allocate (cnts(nprocs_d))
    allocate (dsps(nprocs_d))
    do ip = 0, nprocs_d-1
      dsps(ip+1) = bsna(SLOT_START-1, ip)
      cnts(ip+1) = bsna(SLOT_END, ip) - dsps(ip+1)
    end do

    call MPI_SCATTERV(obset_bufs, cnts, dsps, MPI_INTEGER, obsda%set, cnts(myrank_d+1), MPI_INTEGER, 0, MPI_COMM_d, ierr)
    call MPI_SCATTERV(obidx_bufs, cnts, dsps, MPI_INTEGER, obsda%idx, cnts(myrank_d+1), MPI_INTEGER, 0, MPI_COMM_d, ierr)
    call MPI_SCATTERV(ri_bufs2,   cnts, dsps, MPI_r_size,  obsda%ri,  cnts(myrank_d+1), MPI_r_size,  0, MPI_COMM_d, ierr)
    call MPI_SCATTERV(rj_bufs2,   cnts, dsps, MPI_r_size,  obsda%rj,  cnts(myrank_d+1), MPI_r_size,  0, MPI_COMM_d, ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,4x,F15.7)') '###### obsope_cal:mpi_scatterv:                 ', rrtimer-rrtimer00


    deallocate (cnts, dsps)
    deallocate (obset_bufs, obidx_bufs, ri_bufs2, rj_bufs2)
  end if ! [ myrank_e == 0 ]


  call MPI_BARRIER(MPI_COMM_a, ierr)
  rrtimer00 = MPI_WTIME()


  call MPI_BCAST(obsda%set, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_e, ierr)
  call MPI_BCAST(obsda%idx, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_e, ierr)
  call MPI_BCAST(obsda%ri,  obsda%nobs, MPI_r_size,  0, MPI_COMM_e, ierr)
  call MPI_BCAST(obsda%rj,  obsda%nobs, MPI_r_size,  0, MPI_COMM_e, ierr)


  rrtimer = MPI_WTIME()
  write (6,'(A,4x,F15.7)') '###### obsope_cal:mpi_broadcast:                ', rrtimer-rrtimer00
  rrtimer00 = rrtimer

!if (myrank_a == 37) then
!print *, obsda%nobs
!do n = 1, 100 !obsda%nobs
!print *, obsda%set(n), obsda%idx(n), obsda%ri(n), obsda%rj(n)
!end do
!end if


  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    if ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin) then
      write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is processing member ', &
            im, ', subdomain id #', proc2mem(2,it,myrank+1)

!write(6,*) '%%%%%%', MPI_WTIME(), 0

!!!      nobs = 0

      !!!!!!
      if (nobs_alldomain > 0) then
      !!!!!!

      obsda%qc = iqc_time

      do islot = SLOT_START, SLOT_END
!!!        slot_lb = (real(islot-SLOT_BASE,r_size) - 0.5d0) * SLOT_TINTERVAL
!!!        slot_ub = (real(islot-SLOT_BASE,r_size) + 0.5d0) * SLOT_TINTERVAL
        write (6,'(A,I3,A,F9.1,A,F9.1,A)') 'Slot #', islot-SLOT_START+1, ': time interval (', slot_lb, ',', slot_ub, '] sec'

        call read_ens_history_iter(it,islot,v3dg,v2dg)


  rrtimer = MPI_WTIME()
  write (6,'(A,I3,A,I3,A,4x,F15.7)') '###### obsope_cal:read_ens_history_iter:', it, ':', islot, ':', rrtimer-rrtimer00
  rrtimer00 = rrtimer


!!!        nn_0 = 0

!!!        do iof = 1, OBS_IN_NUM

!!!          if (.not. OBSDA_RUN(iof) .or. obs(iof)%nobs == 0) cycle

!!!          obs_idx_TCX = -1
!!!          obs_idx_TCY = -1
!!!          obs_idx_TCP = -1

!!!          nslot = 0
!!!          nobs_slot = 0

!!!!write(6,*) '%%%===', MPI_WTIME(), 'im:', im, 'islot:', islot, 'iof:', iof

!!!          IF(OBS_IN_FORMAT(iof) /= 3)THEN ! except H08 obs ! H08

!!!            ! do this small computation first, without OpenMP
!!!            nobs_0 = nobs
!!!            nn = nn_0
!!!            do n = 1, obs(iof)%nobs

!!!              select case (obs(iof)%elm(n))
!!!              case (id_tclon_obs)
!!!                obs_idx_TCX = n
!!!                cycle
!!!              case (id_tclat_obs)
!!!                obs_idx_TCY = n
!!!                cycle
!!!              case (id_tcmip_obs)
!!!                obs_idx_TCP = n
!!!                cycle
!!!              end select

!!!              if (obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub) then
!!!                nslot = nslot + 1
!!!                call phys2ij(obs(iof)%lon(n),obs(iof)%lat(n),rig,rjg)
!!!                call rij_g2l_auto(proc,rig,rjg,ritmp,rjtmp)

!!!                if (myrank_d == proc) then
!!!                  nobs = nobs + 1
!!!                  nobs_slot = nobs_slot + 1
!!!                  obsda%set(nobs) = iof
!!!                  obsda%idx(nobs) = n
!!!                  obsda%ri(nobs) = rig  ! obsda%ri: global grid coordinate
!!!                  obsda%rj(nobs) = rjg  !
!!!                  ri(nobs) = ritmp      ! ri: local grid coordinate
!!!                  rj(nobs) = rjtmp      !
!!!                end if ! [ myrank_d == proc ]
!!!              end if ! [ obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub ]

!!!!write (6, *) obsda%set(nobs), obsda%idx(nobs), obsda%ri(nobs), obsda%rj(nobs), ri(nobs), rj(nobs)

!!!                end if
!!!              end if


!!!            end do ! [ n = 1, obs%nobs ]

#ifdef H08
          ELSEIF( OBS_IN_FORMAT(iof) == 3) THEN ! for H08 obs (OBS_IN_FORMAT(iof) = 3) ! H08

            nprof_H08 = 0
            nobs_0 = nobs
            nallprof = obs(iof)%nobs/nch

            ALLOCATE(tmp_ri_H08(nallprof))
            ALLOCATE(tmp_rj_H08(nallprof))
            ALLOCATE(tmp_lon_H08(nallprof))
            ALLOCATE(tmp_lat_H08(nallprof))

            do n = 1, nallprof
              ns = (n - 1) * nch + 1
              if (obs(iof)%dif(ns) > slot_lb .and. obs(iof)%dif(ns) <= slot_ub) then
                nslot = nslot + 1
                call phys2ij(obs(iof)%lon(ns),obs(iof)%lat(ns),rig,rjg)
                call rij_g2l_auto(proc,rig,rjg,ritmp,rjtmp)

                if (myrank_d == proc) then
                  nprof_H08 = nprof_H08 + 1 ! num of prof in myrank node
                  tmp_ri_H08(nprof_H08) = ritmp
                  tmp_rj_H08(nprof_H08) = rjtmp
                  tmp_lon_H08(nprof_H08) = obs(iof)%lon(ns)
                  tmp_lat_H08(nprof_H08) = obs(iof)%lat(ns)

                  nobs = nobs + nch
                  nobs_slot = nobs_slot + 1
                  obsda%set(nobs-nch+1:nobs) = iof
                  obsda%ri(nobs-nch+1:nobs) = rig
                  obsda%rj(nobs-nch+1:nobs) = rjg
                  ri(nobs-nch+1:nobs) = ritmp
                  rj(nobs-nch+1:nobs) = rjtmp
                  do ch = 1, nch
                    obsda%idx(nobs-nch+ch) = ns + ch - 1
                  enddo

                end if ! [ myrank_d == proc ]
              end if ! [ obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub ]
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

            ENDIF

            DEALLOCATE(tmp_ri_H08,tmp_rj_H08)
            DEALLOCATE(tmp_lon_H08,tmp_lat_H08)

#endif
!!!          ENDIF ! end of nobs count [if (OBS_IN_FORMAT(iof) = 3)]


!!!  rrtimer = MPI_WTIME()
!!!  write (6,'(A,I3,A,I3,A,I3,A,F15.7)') '###### obsope_cal:obsope_step_1:        ', it, ':', islot, ':', iof, ':', rrtimer-rrtimer00
!!!  rrtimer00 = rrtimer


          ! then do this heavy computation with OpenMP

!!!          IF(OBS_IN_FORMAT(iof) /= 3)THEN ! H08


!write(6,*) '%%%===', MPI_WTIME(), nobs_0 + 1, nobs

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(nn,n,rk)

            n1 = bsna(islot-1, myrank_d) - bsna(SLOT_START-1, myrank_d) + 1
            n2 = bsna(islot, myrank_d) - bsna(SLOT_START-1, myrank_d)

            do nn = n1, n2
!!!            do nn = nobs_0 + 1, nobs
              iof = obsda%set(nn)
              n = obsda%idx(nn)

!              IF(OBS_IN_FORMAT(iof) /= 3)THEN ! H08 ????????????


              call rij_g2l(myrank_d, obsda%ri(nn), obsda%rj(nn), ril, rjl)


!if (myrank_a == 0) then
!print *, iof, n, obsda%ri(nn), obsda%rj(nn), ril, rjl
!end if


!if (mod(nn,50) == 0) then
!  write(6,*) '%%%%%%', MPI_WTIME(), nn
!end if

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
!!!                  call phys2ijkz(v3dg(:,:,:,iv3dd_hgt),ri(nn),rj(nn),obs(iof)%lev(n),rk,obsda%qc(nn))
                  call phys2ijkz(v3dg(:,:,:,iv3dd_hgt),ril,rjl,obs(iof)%lev(n),rk,obsda%qc(nn))
                end if
              else
!!!                call phys2ijk(v3dg(:,:,:,iv3dd_p),obs(iof)%elm(n),ri(nn),rj(nn),obs(iof)%lev(n),rk,obsda%qc(nn))
                call phys2ijk(v3dg(:,:,:,iv3dd_p),obs(iof)%elm(n),ril,rjl,obs(iof)%lev(n),rk,obsda%qc(nn))
              end if

              if (obsda%qc(nn) == iqc_good) then
                select case (OBS_IN_FORMAT(iof))
                case (1)
!!!                  call Trans_XtoY(obs(iof)%elm(n),ri(nn),rj(nn),rk, &
                  call Trans_XtoY(obs(iof)%elm(n),ril,rjl,rk, &
                                  obs(iof)%lon(n),obs(iof)%lat(n),v3dg,v2dg,obsda%val(nn),obsda%qc(nn))
                case (2)
!!!                  call Trans_XtoY_radar(obs(iof)%elm(n),obs(iof)%meta(1),obs(iof)%meta(2),obs(iof)%meta(3),ri(nn),rj(nn),rk, &
                  call Trans_XtoY_radar(obs(iof)%elm(n),obs(iof)%meta(1),obs(iof)%meta(2),obs(iof)%meta(3),ril,rjl,rk, &
                                        obs(iof)%lon(n),obs(iof)%lat(n),obs(iof)%lev(n),v3dg,v2dg,obsda%val(nn),obsda%qc(nn))
                  if (obsda%qc(nn) == iqc_ref_low) obsda%qc(nn) = iqc_good ! when process the observation operator, we don't care if reflectivity is too small

                !!!!!! may not need to do this at this stage...
                !if (obs(iof)%elm(n) == id_radar_ref_obs) then
                !  obsda%val(nn) = 10.0d0 * log10(obsda%val(nn))
                !end if
                !!!!!!

                end select
              end if

!              ENDIF ! H08 ????????????

!            end do ! [ nn = nobs_0 + 1, nobs ]
            end do ! [ nn = n1, n2 ]
!$OMP END PARALLEL DO

#ifdef H08
          ELSEIF((OBS_IN_FORMAT(iof) == 3).and.(nprof_H08 >=1 ))THEN ! H08
! -- Note: Trans_XtoY_H08 is called without OpenMP but it can use a parallel (with OpenMP) RTTOV routine
!
            !------
            if (.not. USE_OBS(23)) then
              obsda%qc(nobs_0+1:nobs) = iqc_otype
            else
            !------

            ALLOCATE(yobs_H08(nprof_H08*nch))
            ALLOCATE(yobs_H08_clr(nprof_H08*nch))
            ALLOCATE(plev_obs_H08(nprof_H08*nch))
            ALLOCATE(qc_H08(nprof_H08*nch))

            CALL Trans_XtoY_H08(nprof_H08,ri_H08,rj_H08,&
                                lon_H08,lat_H08,v3dg,v2dg,&
                                yobs_H08,plev_obs_H08,&
                                qc_H08,yobs_H08_clr=yobs_H08_clr)

! Clear sky yobs(>0)
! Cloudy sky yobs(<0)

            obsda%qc(nobs_0+1:nobs) = iqc_obs_bad

            ns = 0
            DO nn = nobs_0 + 1, nobs
              ns = ns + 1

              obsda%val(nn) = yobs_H08(ns)
              obsda%qc(nn) = qc_H08(ns)

              if(obsda%qc(nn) == iqc_good)then
                rig = obsda%ri(nn)
                rjg = obsda%rj(nn)

! -- tentative treatment around the TC center --
!                dist_MSLP_TC = sqrt(((rig - MSLP_TC_rig) * DX)**2&
!                                   +((rjg - MSLP_TC_rjg) * DY)**2)

!                if(dist_MSLP_TC <= dist_MSLP_TC_MIN)then
!                  obsda%qc(nn) = iqc_obs_bad
!                endif

! -- Rejecting Himawari-8 obs over the buffer regions. --
                if((rig <= bris) .or. (rig >= brie) .or.&
                   (rjg <= brjs) .or. (rjg >= brje))then
                  obsda%qc(nn) = iqc_obs_bad
                endif
              endif

!
!  NOTE: T.Honda (10/16/2015)
!  The original H08 obs does not inlcude the level information.
!  However, we have the level information derived by RTTOV (plev_obs_H08) here, 
!  so that we substitute the level information into obsda%lev.  
!  The substituted level information is used in letkf_tools.f90
!
              obsda%lev(nn) = plev_obs_H08(ns)
              obsda%val2(nn) = yobs_H08_clr(ns)

!              write(6,'(a,f12.1,i9)')'H08 debug_plev',obsda%lev(nn),nn

            END DO ! [ nn = nobs_0 + 1, nobs ]

            DEALLOCATE(ri_H08, rj_H08)
            DEALLOCATE(lon_H08, lat_H08)
            DEALLOCATE(yobs_H08, plev_obs_H08)
            DEALLOCATE(yobs_H08_clr)
            DEALLOCATE(qc_H08)

            !------
            end if ! [.not. USE_OBS(23)]
            !------

#endif
!!!          ENDIF ! H08

!!!          write (6,'(3A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot = ', nslot
!!!          write (6,'(3A,I6,A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot and processed by rank ' &
!!!                                    , myrank, ' = ', nobs_slot
          write (6,'(3A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot = ', sum(bsna(islot, :)) - sum(bsna(islot-1, :))
          write (6,'(3A,I6,A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot and processed by rank ' &
                                    , myrank, ' = ', bsna(islot, myrank_d) - bsna(islot-1, myrank_d)


  rrtimer = MPI_WTIME()
  write (6,'(A,I3,A,I3,A,I3,A,F15.7)') '###### obsope_cal:obsope_step_2:        ', it, ':', islot, ':', iof, ':', rrtimer-rrtimer00
  rrtimer00 = rrtimer



! ###  -- TC vital assimilation -- ###
          if (obs_idx_TCX > 0 .and. obs_idx_TCY > 0 .and. obs_idx_TCP > 0) then
          if (obs(iof)%dif(obs_idx_TCX) == obs(iof)%dif(obs_idx_TCY) .and. &
              obs(iof)%dif(obs_idx_TCY) == obs(iof)%dif(obs_idx_TCP)) then
           
            if (obs(iof)%dif(obs_idx_TCX) > slot_lb .and. &
              obs(iof)%dif(obs_idx_TCX) <= slot_ub) then
              nslot = nslot + 3 ! TC vital obs should have 3 data (i.e., lon, lat, and MSLP)

              !!! bTC(1,:) : lon, bTC(2,:): lat, bTC(3,:): mslp
              ! bTC(1,:) : tcx (m), bTC(2,:): tcy (m), bTC(3,:): mslp
              allocate(bTC(3,0:MEM_NP-1))

              bTC = 9.99d33

              ! Note: obs(iof)%dat(obs_idx_TCX) is not longitude (deg) but X (m).
              !       Units of the original TC vital position are converted in
              !       subroutine read_obs in common_obs_scale.f90.
              !
              call phys2ij(obs(iof)%lon(obs_idx_TCX),obs(iof)%lat(obs_idx_TCX),rig,rjg) 
              call rij_g2l_auto(proc,rig,rjg,ril,rjl)  
              call search_tc_subdom(rig,rjg,v2dg,bTC(1,myrank_d),bTC(2,myrank_d),bTC(3,myrank_d))
  
!              CALL MPI_BARRIER(MPI_COMM_d,ierr)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,bTC,3*MEM_NP,MPI_r_size,MPI_MIN,MPI_COMM_d,ierr)

              ! Assume MSLP of background TC is lower than 1100 (hPa). 
              bTC_mslp = 1100.0d2
              do n = 0, MEM_NP - 1
                write(6,'(3e20.5)')bTC(1,n),bTC(2,n),bTC(3,n) ! debug
                if (bTC(3,n) < bTC_mslp ) then
                  bTC_mslp = bTC(3,n)
                  bTC_proc = n
                endif
              enddo ! [ n = 0, MEM_NP - 1]

              if (myrank_d == proc) then
                do n = 1, 3
                  nobs = nobs + 1
                  nobs_slot = nobs_slot + 1
                  obsda%set(nobs) = iof
                  if(n==1) obsda%idx(nobs) = obs_idx_TCX
                  if(n==2) obsda%idx(nobs) = obs_idx_TCY
                  if(n==3) obsda%idx(nobs) = obs_idx_TCP
                  obsda%ri(nobs) = rig
                  obsda%rj(nobs) = rjg
                  ri(nobs) = ril
                  rj(nobs) = rjl

                  obsda%val(nobs) = bTC(n,bTC_proc)
                  obsda%qc(nobs) = iqc_good
                enddo ! [ n = 1, 3 ]

              endif
              deallocate(bTC)

            endif ! [ obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub ]
          endif ! [ obs_idx_TCX > 0 ...]
          endif !

!!!          nn_0 = nn_0 + obs(iof)%nobs

!!!        end do ! [ do iof = 1, OBS_IN_NUM ]

 
!      IF(NINT(elem(n)) == id_ps_obs .AND. odat(n) < -100.0d0) THEN
!        CYCLE
!      END IF
!      IF(NINT(elem(n)) == id_ps_obs) THEN
!        CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,dz)
!        rk = rlev(n) - dz
!        IF(ABS(rk) > threshold_dz) THEN ! pressure adjustment threshold
!!          WRITE(6,'(A)') '* PS obs vertical adjustment beyond threshold'
!!          WRITE(6,'(A,F10.2,A,F6.2,A,F6.2,A)') '*   dz=',rk,&
!!           & ', (lon,lat)=(',elon(n),',',elat(n),')'
!          CYCLE
!        END IF
!      END IF

      end do ! [ islot = SLOT_START, SLOT_END ]

!write(6,*) '%%%%%%', MPI_WTIME(), nobs

      !!!!!!
      end if ! [ nobs_alldomain > 0 ]
      !!!!!!

!!!      if (it == 1) then
!!!        obsda%nobs = nobs
!!!      else if (nobs /= obsda%nobs) then
!!!        write (6, '(A)') '[Error] numbers of observations found are different among members.'
!!!        stop
!!!      end if



      nobs = obsda%nobs




      write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' finishes processing member ', &
            im, ', subdomain id #', proc2mem(2,it,myrank+1)
      write (6,'(A,I8,A)') ' -- ', nobs, ' observations found'


  rrtimer00 = MPI_WTIME()


      if (OBSDA_OUT) then
        write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is writing observations for member ', &
              im, ', subdomain id #', proc2mem(2,it,myrank+1)
        if (im <= MEMBER) then
          call file_member_replace(im, OBSDA_OUT_BASENAME, obsdafile)
        else if (im == mmean) then
          obsdafile = OBSDA_MEAN_OUT_BASENAME
        else if (im == mmdet) then
          obsdafile = OBSDA_MDET_OUT_BASENAME
        end if
        write (obsda_suffix(2:7),'(I6.6)') proc2mem(2,it,myrank+1)
        call write_obs_da(trim(obsdafile)//obsda_suffix,obsda,0)


  rrtimer = MPI_WTIME()
  write (6,'(A,I3,A,8x,F15.7)') '###### obsope_cal:write_obs_da:         ', it, ':', rrtimer-rrtimer00
  rrtimer00 = rrtimer


      end if

      if (present(obsda_return)) then
        ! variables without an ensemble dimension
        if (it == 1) then
          if (present(nobs_extern)) then
            obsda_return%nobs = nobs + nobs_extern ! additional space for externally processed observations
          else
            obsda_return%nobs = nobs
          end if
          call obs_da_value_allocate(obsda_return, nensobs)
        end if
        if (nobs > 0) then
          if (it == 1) then
            obsda_return%set(1:nobs) = obsda%set(1:nobs)
            obsda_return%idx(1:nobs) = obsda%idx(1:nobs)
            obsda_return%ri(1:nobs) = obsda%ri(1:nobs)
            obsda_return%rj(1:nobs) = obsda%rj(1:nobs)
            obsda_return%qc(1:nobs) = obsda%qc(1:nobs)
#ifdef H08
            obsda_return%lev(1:nobs) = obsda%lev(1:nobs)
            obsda_return%val2(1:nobs) = obsda%val2(1:nobs)
#endif
          else
            obsda_return%qc(1:nobs) = max(obsda_return%qc(1:nobs), obsda%qc(1:nobs))
#ifdef H08
            if (im <= MEMBER) then ! only consider lev, val2 from members, not from the means
              obsda_return%lev(1:nobs) = obsda_return%lev(1:nobs) + obsda%lev(1:nobs)
              obsda_return%val2(1:nobs) = obsda_return%val2(1:nobs) + obsda%val2(1:nobs)
            end if
#endif
          end if

          ! variables with an ensemble dimension
          if (im == mmdetin) then
            obsda_return%ensval(mmdetobs,1:nobs) = obsda%val(1:nobs)
          else
            obsda_return%ensval(im,1:nobs) = obsda%val(1:nobs)
          end if
        end if ! [ nobs > 0 ]


  rrtimer = MPI_WTIME()
  write (6,'(A,I3,A,8x,F15.7)') '###### obsope_cal:obsda_return:         ', it, ':', rrtimer-rrtimer00
  rrtimer00 = rrtimer


      end if ! [ present(obsda_return) ]

    end if ! [ (im >= 1 .and. im <= MEMBER) .or. im == mmdetin ]

  end do ! [ it = 1, nitmax ]

  call obs_da_value_deallocate(obsda)

!  deallocate ( ri, rj, v3dg, v2dg )
  deallocate ( v3dg, v2dg )
  deallocate ( bsna )

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

  write (6,'(A,I6.6,A,I6.6)') 'MYRANK ',myrank,' is processing subdomain id #', proc2mem(2,1,myrank+1)

  allocate ( v3dg (nlevh,nlonh,nlath,nv3dd) )
  allocate ( v2dg (nlonh,nlath,nv2dd) )

  do iof = 1, OBS_IN_NUM
    obs(iof)%dat = 0.0d0
  end do

  nobs = 0
  do islot = SLOT_START, SLOT_END
    slot_lb = (real(islot-SLOT_BASE,r_size) - 0.5d0) * SLOT_TINTERVAL
    slot_ub = (real(islot-SLOT_BASE,r_size) + 0.5d0) * SLOT_TINTERVAL
    write (6,'(A,I3,A,F9.1,A,F9.1,A)') 'Slot #', islot-SLOT_START+1, ': time interval (', slot_lb, ',', slot_ub, '] sec'

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

!-------------------------------------------------------------------------------

  write (6,'(A,I6.6,A,I6.6)') 'MYRANK ', myrank, ' is processing subdomain id #', proc2mem(2,1,myrank+1)

  do j = 1, nlat
    rj = real(j + JHALO, r_size)

    do i = 1, nlon
      ri = real(i + IHALO, r_size)
      call MPRJ_xy2lonlat((ri-1.0d0) * DX + GRID_CX(1), (rj-1.0d0) * DY + GRID_CY(1), lon, lat)
      lon = lon * rad2deg
      lat = lat * rad2deg

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
!            case (id_H08IR_obs)               !!!!!! H08 as 2D observations ???
!              call Trans_XtoY_radar_H08(...)
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
