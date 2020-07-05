PROGRAM letkf
!=======================================================================
!
! [PURPOSE:] Main program of LETKF
!
! [HISTORY:]
!   01/16/2009   Takemasa Miyoshi  created
!   October 2014 Guo-Yuan Lien     modified for SCALE model
!   ............ See git history for the following revisions
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
  USE common_obs_scale
  USE common_nml
  USE letkf_obs
  USE letkf_tools
  use obsope_tools, only: &
    obsope_cal
  IMPLICIT NONE

  REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: anal3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: anal2d(:,:,:)

  character(len=7) :: stdoutf='-000000'
  character(len=6400) :: cmd1, cmd2, icmd
  character(len=10) :: myranks
  integer :: iarg, iof

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  call initialize_mpi_scale
  call mpi_timer('', 1)

  if (command_argument_count() >= 3) then
    write (myranks, '(I10)') myrank
    call get_command_argument(3, icmd)
    cmd1 = 'bash ' // trim(icmd) // ' letkf_1' // ' ' // trim(myranks)
    cmd2 = 'bash ' // trim(icmd) // ' letkf_2' // ' ' // trim(myranks)
    do iarg = 4, command_argument_count()
      call get_command_argument(iarg, icmd)
      cmd1 = trim(cmd1) // ' ' // trim(icmd)
      cmd2 = trim(cmd2) // ' ' // trim(icmd)
    end do
  end if

  if (command_argument_count() >= 2) then
    call get_command_argument(2, icmd)
    if (trim(icmd) /= '') then
      WRITE(stdoutf(2:7), '(I6.6)') myrank
!      WRITE(6,'(3A,I6.6)') 'STDOUT goes to ',trim(icmd)//stdoutf,' for MYRANK ', myrank
      OPEN(6,FILE=trim(icmd)//stdoutf)
      WRITE(6,'(A,I6.6,2A)') 'MYRANK=',myrank,', STDOUTF=',trim(icmd)//stdoutf
    end if
  end if

  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '  LOCAL ENSEMBLE TRANSFORM KALMAN FILTERING  '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '   LL      EEEEEE  TTTTTT  KK  KK  FFFFFF    '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(6,'(A)') '   LL      EEEEE     TT    KKK     FFFFF     '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(6,'(A)') '   LLLLLL  EEEEEE    TT    KK  KK  FF        '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '             WITHOUT LOCAL PATCH             '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '          Coded by Takemasa Miyoshi          '
  WRITE(6,'(A)') '  Based on Ott et al (2004) and Hunt (2005)  '
  WRITE(6,'(A)') '  Tested by Miyoshi and Yamane (2006)        '
  WRITE(6,'(A)') '============================================='

!-----------------------------------------------------------------------
! Pre-processing scripts
!-----------------------------------------------------------------------

  if (command_argument_count() >= 3) then
    write (6,'(A)') 'Run pre-processing scripts'
    write (6,'(A,I6.6,3A)') 'MYRANK ',myrank,' is running a script: [', trim(cmd1), ']'
    call system(trim(cmd1))
  end if

  call mpi_timer('PRE_SCRIPT', 1, barrier=MPI_COMM_WORLD)

!-----------------------------------------------------------------------

  call set_common_conf

  if (DET_RUN) then
    call set_mem_node_proc(MEMBER+2)
  else
    call set_mem_node_proc(MEMBER+1)
  end if
  call set_scalelib('LETKF')

  if (myrank_use) then

    call set_common_scale
    call set_common_mpi_scale
    call set_common_obs_scale

    call mpi_timer('INITIALIZE', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Read observations
!-----------------------------------------------------------------------

    allocate (obs(OBS_IN_NUM))
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
    call set_common_mpi_grid

    allocate (gues3d(nij1,nlev,nens,nv3d))
    allocate (gues2d(nij1,nens,nv2d))
    allocate (anal3d(nij1,nlev,nens,nv3d))
    allocate (anal2d(nij1,nens,nv2d))

    call mpi_timer('SET_GRID', 1, barrier=MPI_COMM_a)

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
      call write_ensmean(GUES_MEAN_INOUT_BASENAME, gues3d, gues2d, calced=.false., monit_step=1)
    else
      call write_ensmean(GUES_MEAN_INOUT_BASENAME, gues3d, gues2d, calced=.false.)
    end if

    if (GUES_SPRD_OUT) then
      call write_enssprd(GUES_SPRD_OUT_BASENAME, gues3d, gues2d)
    end if

    call mpi_timer('GUES_MEAN', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------

    !
    ! LETKF
    !
    call das_letkf(gues3d,gues2d,anal3d,anal2d)

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
      call write_enssprd(ANAL_SPRD_OUT_BASENAME, anal3d, anal2d)
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

!!-----------------------------------------------------------------------
!! Monitor
!!-----------------------------------------------------------------------

    do iof = 1, OBS_IN_NUM
      call obs_info_deallocate(obs(iof))
    end do
    deallocate (obs)
    deallocate (gues3d, gues2d, anal3d, anal2d)

    call unset_common_mpi_scale

  end if ! [ myrank_use ]

  call unset_scalelib

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_WORLD)

!-----------------------------------------------------------------------
! Post-processing scripts
!-----------------------------------------------------------------------

  if (command_argument_count() >= 3) then
    write (6,'(A)') 'Run post-processing scripts'
    write (6,'(A,I6.6,3A)') 'MYRANK ',myrank,' is running a script: [', trim(cmd2), ']'
    call system(trim(cmd2))
  end if

  call mpi_timer('POST_SCRIPT', 1, barrier=MPI_COMM_WORLD)

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  call finalize_mpi_scale

  STOP
END PROGRAM letkf
