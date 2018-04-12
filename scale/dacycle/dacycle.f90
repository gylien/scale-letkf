program dacycle
!=======================================================================
!
! [PURPOSE:] Single program for SCALE/LETKF data assimilation cycles
!
!=======================================================================
!$USE OMP_LIB
  use common
  use common_mpi
  use common_nml
  use common_scale
  use common_scalerm
  use common_mpi_scale
  use common_obs_scale
  use letkf_obs
  use letkf_tools
  use obsope_tools, only: &
    obsope_cal

  use scale_stdio, only: &
    IO_L, &
    IO_FID_LOG
  use scale_time, only: &
    TIME_NOWDATE, &
    TIME_NOWMS, &
    TIME_NOWSTEP, &
    TIME_DTSEC, &
    TIME_NSTEP
  use scale_file_history, only: &
    FILE_HISTORY_write, &
    FILE_HISTORY_set_nowdate
  use scale_monitor, only: &
    MONIT_write
  use mod_admin_restart, only: &
#ifdef SCALEUV
    ADMIN_restart_write, &
    ADMIN_restart_write_additional
#else
    ADMIN_restart_write
#endif
  use mod_admin_time, only: &
    ADMIN_TIME_checkstate, &
    ADMIN_TIME_advance, &
    TIME_DOATMOS_step, &
    TIME_DOLAND_step, &
    TIME_DOURBAN_step, &
    TIME_DOOCEAN_step, &
    TIME_DOresume, &
    TIME_DOend, &
    TIME_DOATMOS_restart, &
    TIME_DTSEC_ATMOS_RESTART
  use mod_atmos_admin, only: &
    ATMOS_do
  use mod_atmos_driver, only: &
    ATMOS_driver, &
    ATMOS_driver_finalize
  use mod_ocean_admin, only: &
    OCEAN_do
  use mod_ocean_driver, only: &
    OCEAN_driver
  use mod_land_admin, only: &
    LAND_do
  use mod_land_driver, only: &
    LAND_driver
  use mod_urban_admin, only: &
    URBAN_do
  use mod_urban_driver, only: &
    URBAN_driver
  use mod_user, only: &
    USER_step
  implicit none

  real(r_size), allocatable :: gues3d(:,:,:,:)
  real(r_size), allocatable :: gues2d(:,:,:)
  real(r_size), allocatable:: anal3d(:,:,:,:)
  real(r_size), allocatable :: anal2d(:,:,:)

  character(len=7) :: stdoutf='-000000'
  character(len=6400) :: icmd

  logical :: anal_mem_out_now
  logical :: anal_mean_out_now
  logical :: anal_mdet_out_now
  logical :: gues_mean_out_now
  logical :: gues_sprd_out_now
  logical :: anal_sprd_out_now

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  call initialize_mpi_scale
  call mpi_timer('', 1)

  if (command_argument_count() >= 2) then
    call get_command_argument(2, icmd)
    if (trim(icmd) /= '') then
      write (stdoutf(2:7), '(I6.6)') myrank
!      WRITE (6,'(3A,I6.6)') 'STDOUT goes to ', trim(icmd)//stdoutf, ' for MYRANK ', myrank
      open (6, file=trim(icmd)//stdoutf)
      write (6,'(A,I6.6,2A)') 'MYRANK=', myrank, ', STDOUTF=', trim(icmd)//stdoutf
    end if
  end if

  write (6, '(A)') '============================================='
  write (6, '(A)') '  LOCAL ENSEMBLE TRANSFORM KALMAN FILTERING  '
  write (6, '(A)') '                                             '
  write (6, '(A)') '   LL      EEEEEE  TTTTTT  KK  KK  FFFFFF    '
  write (6, '(A)') '   LL      EE        TT    KK KK   FF        '
  write (6, '(A)') '   LL      EEEEE     TT    KKK     FFFFF     '
  write (6, '(A)') '   LL      EE        TT    KK KK   FF        '
  write (6, '(A)') '   LLLLLL  EEEEEE    TT    KK  KK  FF        '
  write (6, '(A)') '                                             '
  write (6, '(A)') '             WITHOUT LOCAL PATCH             '
  write (6, '(A)') '                                             '
  write (6, '(A)') '          Coded by Takemasa Miyoshi          '
  write (6, '(A)') '  Based on Ott et al (2004) and Hunt (2005)  '
  write (6, '(A)') '  Tested by Miyoshi and Yamane (2006)        '
  write (6, '(A)') '============================================='

!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------

  call set_common_conf(nprocs)

  call set_mem_node_proc(MEMBER_RUN)

  call scalerm_setup('DACYCLE')

  call mpi_timer('INITIALIZE', 1, barrier=MPI_COMM_WORLD)

  if (myrank_use) then

    icycle = 0
    lastcycle = int((TIME_DTSEC * TIME_NSTEP + 1.0d-6) / TIME_DTSEC_ATMOS_RESTART)

    write (6, '(A,I7)') 'Total cycle numbers:', lastcycle

!-----------------------------------------------------------------------
! Main loop
!-----------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START TIMESTEP ++++++'
    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_Loop', 0)

    do

      anal_mem_out_now = .false.
      anal_mean_out_now = .false.
      anal_mdet_out_now = .false.
      gues_mean_out_now = .false.
      gues_sprd_out_now = .false.
      anal_sprd_out_now = .false.

      ! report current time
      call ADMIN_TIME_checkstate

      if ( TIME_DOresume ) then
        ! resume state from restart files
        if (DIRECT_TRANSFER .and. icycle >= 1) then
          if (LOG_LEVEL >= 1) then
            write (6, '(A,I6,A)') '[Info] Cycle #', icycle, ': Use direct transfer; skip reading restart files'
          end if
          call resume_state(do_restart_read=.false.)
        else
          call resume_state(do_restart_read=.true.)
        end if

        ! history&monitor file output
        call MONIT_write('MAIN')
        call FILE_HISTORY_write ! if needed
      end if

      ! time advance
      call ADMIN_TIME_advance
      call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP )

      ! user-defined procedure
      call USER_step

      ! change to next state
      if( OCEAN_do .AND. TIME_DOOCEAN_step ) call OCEAN_driver
      if( LAND_do  .AND. TIME_DOLAND_step  ) call LAND_driver
      if( URBAN_do .AND. TIME_DOURBAN_step ) call URBAN_driver
      if( ATMOS_do .AND. TIME_DOATMOS_step ) call ATMOS_driver

      ! history&monitor file output
      call MONIT_write('MAIN')
      call FILE_HISTORY_write

      ! restart output before LETKF
      if (DIRECT_TRANSFER) then
        if (LOG_LEVEL >= 1 .and. TIME_DOATMOS_restart) then
          write (6, '(A,I6,A)') '[Info] Cycle #', icycle, ': Use direct transfer; skip writing restart files before LETKF'
        end if
      else
        call ADMIN_restart_write
      end if
      call ADMIN_restart_write_additional !!!!!! To do: control additional restart outputs for gues_mean, gues_sprd, and anal_sprd

      !-------------------------------------------------------------------------
      ! LETKF section start
      !-------------------------------------------------------------------------
      if (TIME_DOATMOS_restart) then

        call mpi_timer('SCALE', 1, barrier=MPI_COMM_a)

        icycle = icycle + 1

        if (ANAL_OUT_FREQ >= 1 .and. (mod(icycle, ANAL_OUT_FREQ) == 0 .or. icycle == lastcycle)) then
          anal_mem_out_now = .true.
        end if
        if (ANAL_MEAN_OUT_FREQ >= 1 .and. (mod(icycle, ANAL_MEAN_OUT_FREQ) == 0 .or. icycle == lastcycle)) then
          anal_mean_out_now = .true.
        end if
        if (ANAL_MDET_OUT_FREQ >= 1 .and. (mod(icycle, ANAL_MDET_OUT_FREQ) == 0 .or. icycle == lastcycle)) then
          anal_mdet_out_now = .true.
        end if

        if (GUES_MEAN_OUT_FREQ >= 1 .and. (mod(icycle, GUES_MEAN_OUT_FREQ) == 0 .or. icycle == lastcycle)) then
          gues_mean_out_now = .true.
        end if
        if (GUES_SPRD_OUT_FREQ >= 1 .and. (mod(icycle, GUES_SPRD_OUT_FREQ) == 0 .or. icycle == lastcycle)) then
          gues_sprd_out_now = .true.
        end if
        if (ANAL_SPRD_OUT_FREQ >= 1 .and. (mod(icycle, ANAL_SPRD_OUT_FREQ) == 0 .or. icycle == lastcycle)) then
          anal_sprd_out_now = .true.
        end if

        !-----------------------------------------------------------------------
        ! LETKF setups
        !-----------------------------------------------------------------------

        call timelabel_update(TIME_DTSEC_ATMOS_RESTART)

        if (icycle == 1) then
          call set_common_scale
          call set_common_mpi_scale
          call set_common_obs_scale

          allocate (obs(OBS_IN_NUM))

          call mpi_timer('INIT_LETKF', 1, barrier=MPI_COMM_a)
        end if

        !-----------------------------------------------------------------------
        ! Read observations
        !-----------------------------------------------------------------------

        call read_obs_all_mpi(obs)

        call mpi_timer('READ_OBS', 1, barrier=MPI_COMM_a)

        !-----------------------------------------------------------------------
        ! Observation operator
        !-----------------------------------------------------------------------

        if (OBSDA_IN) then
          call get_nobs_da_mpi(nobs_extern)
        else
          nobs_extern = 0
        end if

        !
        ! Compute observation operator, return the results in obsda
        ! with additional space for externally processed observations
        !
        call obsope_cal(obsda_return=obsda, nobs_extern=nobs_extern)

        call mpi_timer('OBS_OPERATOR', 1, barrier=MPI_COMM_a)

        !-----------------------------------------------------------------------
        ! Process observation data
        !-----------------------------------------------------------------------

        call set_letkf_obs

        call mpi_timer('PROCESS_OBS', 1, barrier=MPI_COMM_a)

        !-----------------------------------------------------------------------
        ! First guess ensemble
        !-----------------------------------------------------------------------

        !
        ! LETKF GRID setup
        !
        if (icycle == 1) then
          call set_common_mpi_grid

          allocate (gues3d(nij1,nlev,nens,nv3d))
          allocate (gues2d(nij1,nens,nv2d))
          allocate (anal3d(nij1,nlev,nens,nv3d))
          allocate (anal2d(nij1,nens,nv2d))

          call mpi_timer('SET_GRID', 1, barrier=MPI_COMM_a)
        end if

        !
        ! READ GUES
        !
        call read_ens_mpi(gues3d, gues2d)

        if (ENS_WITH_MDET .and. mmdetin /= mmdet) then
          gues3d(:,:,mmdet,:) = gues3d(:,:,mmdetin,:)
          gues2d(:,mmdet,:) = gues2d(:,mmdetin,:)
        end if

        call mpi_timer('READ_GUES', 1, barrier=MPI_COMM_a)

        !
        ! WRITE ENS MEAN and SPRD
        !
        if (DEPARTURE_STAT .and. LOG_LEVEL >= 1) then
          call write_ensmean(trim(GUES_MEAN_OUT_BASENAME) // trim(timelabel_anal), gues3d, gues2d, &
                             calced=.false., mean_out=gues_mean_out_now, monit_step=1)
        else
          call write_ensmean(trim(GUES_MEAN_OUT_BASENAME) // trim(timelabel_anal), gues3d, gues2d, &
                             calced=.false., mean_out=gues_mean_out_now)
        end if

        if (gues_sprd_out_now) then
          call write_enssprd(trim(GUES_SPRD_OUT_BASENAME) // trim(timelabel_anal), gues3d, gues2d)
        end if

        call mpi_timer('GUES_MEAN', 1, barrier=MPI_COMM_a)

        !-----------------------------------------------------------------------
        ! Data Assimilation
        !-----------------------------------------------------------------------

        !
        ! LETKF
        !
        call das_letkf(gues3d, gues2d, anal3d, anal2d)

        call mpi_timer('DAS_LETKF', 1, barrier=MPI_COMM_a)

        !-----------------------------------------------------------------------
        ! Analysis ensemble
        !-----------------------------------------------------------------------

        !
        ! COMPUTE ENS MEAN and SPRD
        !
        call ensmean_grd(MEMBER, nens, nij1, anal3d, anal2d)
        ! write analysis mean later in write_ens_mpi

        if (anal_sprd_out_now) then
          call write_enssprd(trim(ANAL_SPRD_OUT_BASENAME) // trim(timelabel_anal), anal3d, anal2d)
        end if

        call mpi_timer('ANAL_MEAN', 1, barrier=MPI_COMM_a)

        !
        ! WRITE ANAL and ENS MEAN
        !
        if (DEPARTURE_STAT .and. LOG_LEVEL >= 1) then
          call write_ens_mpi(anal3d, anal2d, monit_step=2)
        else
          call write_ens_mpi(anal3d, anal2d)
        end if

        call mpi_timer('WRITE_ANAL', 1, barrier=MPI_COMM_a)

        !-----------------------------------------------------------------------

      end if ! [ TIME_DOATMOS_restart ]
      !-------------------------------------------------------------------------
      ! LETKF section end
      !-------------------------------------------------------------------------

      ! restart output after LETKF
      if (DIRECT_TRANSFER .and. (anal_mem_out_now .or. anal_mean_out_now .or. anal_mdet_out_now)) then
        !!!!!! To do: control restart outputs separately for members, mean, and mdet
        if (LOG_LEVEL >= 1 .and. TIME_DOATMOS_restart) then
          write (6, '(A,I6,A)') '[Info] Cycle #', icycle, ': Use direct transfer; writing restart (analysis) files after LETKF'
        end if
        call ADMIN_restart_write
      end if

      if( TIME_DOend ) exit

      if( IO_L ) call flush(IO_FID_LOG)

    end do

    call PROF_rapend('Main_Loop', 0)

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

    ! LETKF finalize
    !-------------------------------------------------------------------

    deallocate (obs)
    deallocate (gues3d, gues2d, anal3d, anal2d)

    call unset_common_mpi_scale

    !-------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END TIMESTEP ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    call PROF_setprefx('FIN')

    call PROF_rapstart('All', 1)

    if( ATMOS_do ) call ATMOS_driver_finalize

  end if ! [ myrank_use ]

  call scalerm_finalize('DACYCLE')

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_WORLD)

  call finalize_mpi_scale

  stop
end program dacycle
