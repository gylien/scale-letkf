program dacycle
!=======================================================================
!
! [PURPOSE:] Single program for SCALE/LETKF data assimilation cycles
!
!=======================================================================
!$USE OMP_LIB
  use common
  use common_mpi
  use common_scale
  use common_mpi_scale
  use common_obs_scale
  use common_nml
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
    TIME_NOWSTEP
  use scale_file_history, only: &
    FILE_HISTORY_write, &
    FILE_HISTORY_set_nowdate
  use scale_monitor, only: &
    MONIT_write
  use mod_admin_restart, only: &
    ADMIN_restart_write
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

  integer :: icycle

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

  if (DET_RUN) then
    call set_mem_node_proc(MEMBER+2)
  else
    call set_mem_node_proc(MEMBER+1)
  end if
  call set_scalelib('DACYCLE')

  call mpi_timer('INITIALIZE', 1, barrier=MPI_COMM_WORLD)

  if (myrank_use) then

    icycle = 0

!-----------------------------------------------------------------------
! Main loop
!-----------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START TIMESTEP ++++++'
    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_Loop', 0)

    do

      ! report current time
      call ADMIN_TIME_checkstate

      if ( TIME_DOresume ) then
         ! resume state from restart files
         call resume_state

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

      ! restart output
      call ADMIN_restart_write

      !-------------------------------------------------------------------------
      ! LETKF section start
      !-------------------------------------------------------------------------
      if (TIME_DOATMOS_restart) then

        call mpi_timer('SCALE', 1, barrier=MPI_COMM_a)

        icycle = icycle + 1

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

        if (DET_RUN .and. mmdetin /= mmdet) then
          gues3d(:,:,mmdet,:) = gues3d(:,:,mmdetin,:)
          gues2d(:,mmdet,:) = gues2d(:,mmdetin,:)
        end if

        call mpi_timer('READ_GUES', 1, barrier=MPI_COMM_a)

        !
        ! WRITE ENS MEAN and SPRD
        !
        if (DEPARTURE_STAT .and. LOG_LEVEL >= 1) then
          call write_ensmean(trim(GUES_MEAN_INOUT_BASENAME) // trim(timelabel_anal), gues3d, gues2d, calced=.false., monit_step=1)
        else
          call write_ensmean(trim(GUES_MEAN_INOUT_BASENAME) // trim(timelabel_anal), gues3d, gues2d, calced=.false.)
        end if

        if (GUES_SPRD_OUT) then
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

        if (ANAL_SPRD_OUT) then
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

      if( TIME_DOend ) exit

      if( IO_L ) call flush(IO_FID_LOG)

    end do

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

    ! LETKF finalize
    !-------------------------------------------------------------------

    deallocate (obs)
    deallocate (gues3d, gues2d, anal3d, anal2d)

!    call unset_common_mpi_scale !!!!!! cause unknown MPI error in 'PROF_rapreport' in 'unset_scalelib' if enabled

    !-------------------------------------------------------------------

    call PROF_rapend('Main_Loop', 0)

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END TIMESTEP ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    call PROF_setprefx('FIN')

    call PROF_rapstart('All', 1)

    if( ATMOS_do ) call ATMOS_driver_finalize

  end if ! [ myrank_use ]

  call unset_scalelib('DACYCLE')

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_WORLD)

  call finalize_mpi_scale

  stop
end program dacycle



  subroutine resume_state
    use mod_atmos_driver, only: &
       ATMOS_driver_resume1, &
       ATMOS_driver_resume2, &
       ATMOS_SURFACE_SET
    use mod_ocean_driver, only: &
       OCEAN_driver_resume, &
       OCEAN_SURFACE_SET
    use mod_land_driver, only: &
       LAND_driver_resume, &
       LAND_SURFACE_SET
    use mod_urban_driver, only: &
       URBAN_driver_resume, &
       URBAN_SURFACE_SET
    use mod_atmos_vars, only: &
       ATMOS_vars_calc_diagnostics, &
       ATMOS_vars_history_setpres, &
       ATMOS_vars_restart_read
    use mod_ocean_vars, only: &
       OCEAN_vars_restart_read
    use mod_land_vars, only: &
       LAND_vars_restart_read
    use mod_urban_vars, only: &
       URBAN_vars_restart_read
    use mod_user, only: &
       USER_resume0, &
       USER_resume
    use mod_atmos_admin, only: &
       ATMOS_do
    use mod_ocean_admin, only: &
       OCEAN_do
    use mod_land_admin, only: &
       LAND_do
    use mod_urban_admin, only: &
       URBAN_do
    use mod_admin_restart, only: &
       ADMIN_restart_read
    implicit none
    !---------------------------------------------------------------------------

    ! read restart data
    call ADMIN_restart_read

    ! setup user-defined procedure before setup of other components
    call USER_resume0

    if ( ATMOS_do ) then
       ! calc diagnostics
       call ATMOS_vars_calc_diagnostics
       call ATMOS_vars_history_setpres
    endif

    ! setup surface condition
    if( ATMOS_do ) call ATMOS_SURFACE_SET( countup=.false. )
    if( OCEAN_do ) call OCEAN_SURFACE_SET( countup=.false. )
    if( LAND_do  ) call LAND_SURFACE_SET ( countup=.false. )
    if( URBAN_do ) call URBAN_SURFACE_SET( countup=.false. )

    ! setup submodel driver
    if( ATMOS_do ) call ATMOS_driver_resume1
    if( OCEAN_do ) call OCEAN_driver_resume
    if( LAND_do  ) call LAND_driver_resume
    if( URBAN_do ) call URBAN_driver_resume
    if( ATMOS_do ) call ATMOS_driver_resume2

    ! setup user-defined procedure
    call USER_resume

    return
  end subroutine resume_state
