program obssim
!=======================================================================
!
! [PURPOSE:] Main program of model-to-observation simulator
!
!=======================================================================
!$use OMP_LIB
  use scale_precision, only: RP

  use common, only: r_size
  use common_mpi
  use common_scale
  use common_mpi_scale
  use common_obs_scale
  use common_nml
  use obsope_tools

  implicit none

  real(RP), allocatable :: v3dg(:,:,:,:)
  real(RP), allocatable :: v2dg(:,:,:)
  real(RP), allocatable :: topog(:,:)
  real(r_size), allocatable :: v3dgh(:,:,:,:)
  real(r_size), allocatable :: v2dgh(:,:,:)

  real(r_size), allocatable :: v3dgsim(:,:,:,:)
  real(r_size), allocatable :: v2dgsim(:,:,:)

  integer :: it

  real(r_size) :: rtimer00, rtimer
  integer :: ierr
  character(len=7) :: stdoutf = '-000000'
  character(len=11) :: timer_fmt = '(A30,F10.2)'
  character(len=6400) :: icmd

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  call initialize_mpi_scale
  rtimer00 = MPI_WTIME()

  if (command_argument_count() >= 2) then
    call get_command_argument(2, icmd)
    if (trim(icmd) /= '') then
      write (stdoutf(2:7), '(I6.6)') myrank
!      write (6,'(3A,I6.6)') 'STDOUT goes to ', trim(icmd)//stdoutf, ' for MYRANK ', myrank
      open (6, file=trim(icmd)//stdoutf)
      write (6,'(A,I6.6,2A)') 'MYRANK=', myrank, ', STDOUTF=', trim(icmd)//stdoutf
    end if
  end if

!-----------------------------------------------------------------------

  call set_common_conf(nprocs)

  call read_nml_obssim
  call read_nml_letkf_radar
  call read_nml_letkf_h08

  call set_mem_node_proc(1, NNODES, PPN, MEM_NODES, MEM_NP)
  call set_scalelib

  if (myrank_use) then

    call set_common_scale
    call set_common_mpi_scale
    call set_common_obs_scale

    call MPI_BARRIER(MPI_COMM_a, ierr)
    rtimer = MPI_WTIME()
    write (6, timer_fmt) '### TIMER(INITIALIZE):', rtimer-rtimer00
    rtimer00=rtimer

!-----------------------------------------------------------------------
! Read restart/history input data and run the model-to-obs simulator
!-----------------------------------------------------------------------

    allocate (v3dgh(nlevh, nlonh, nlath, nv3dd))
    allocate (v2dgh(       nlonh, nlath, nv2dd))
    allocate (v3dgsim(nlev, nlon, nlat, OBSSIM_NUM_3D_VARS))
    allocate (v2dgsim(      nlon, nlat, OBSSIM_NUM_2D_VARS))

    if (trim(OBSSIM_IN_TYPE) == 'restart') then

      allocate (v3dg(nlev,nlon,nlat,nv3d))
      allocate (v2dg(nlon,nlat,nv2d))
      allocate (topog(nlon,nlat))

      call read_restart(OBSSIM_RESTART_IN_BASENAME, v3dg, v2dg)
      call state_trans(v3dg)
!      call read_topo(OBSSIM_TOPO_IN_BASENAME, topog)
      topog = 0.0_r_size
      call state_to_history(v3dg, v2dg, topog, v3dgh, v2dgh)

      deallocate (v3dg, v2dg, topog)

      call obssim_cal(v3dgh, v2dgh, v3dgsim, v2dgsim, stggrd=1)

      call write_grd_mpi(OBSSIM_GRADS_OUT_NAME, OBSSIM_NUM_3D_VARS, OBSSIM_NUM_2D_VARS, &
                         1, v3dgsim, v2dgsim)

    else if (trim(OBSSIM_IN_TYPE) == 'history') then

      do it = OBSSIM_TIME_START, OBSSIM_TIME_END
        call read_history(OBSSIM_HISTORY_IN_BASENAME, it, v3dgh, v2dgh)

        call obssim_cal(v3dgh, v2dgh, v3dgsim, v2dgsim)

        call write_grd_mpi(OBSSIM_GRADS_OUT_NAME, OBSSIM_NUM_3D_VARS, OBSSIM_NUM_2D_VARS, &
                           it, v3dgsim, v2dgsim, obsout=OBSSIM_OBSOUT)
      end do

    else
      write (6, '(2A)') "[Error] Unsupported 'OBSSIM_IN_TYPE': ", trim(OBSSIM_IN_TYPE)
    end if

    deallocate (v3dgh, v2dgh)
    deallocate (v3dgsim, v2dgsim)

    call MPI_BARRIER(MPI_COMM_a, ierr)
    rtimer = MPI_WTIME()
    write (6, timer_fmt) '### TIMER(OBSSIM):', rtimer-rtimer00
    rtimer00=rtimer

    call unset_common_mpi_scale

    call unset_scalelib

  else ! [ myrank_use ]

    write (6, '(A,I6.6,A)') 'MYRANK=', myrank, ': This process is not used!'

  end if ! [ myrank_use ]

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  rtimer = MPI_WTIME()
  write (6, timer_fmt) '### TIMER(FINALIZE):', rtimer-rtimer00
  rtimer00=rtimer

  call finalize_mpi_scale

  stop
end program obssim
