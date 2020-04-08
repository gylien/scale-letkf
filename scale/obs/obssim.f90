program obssim
!=======================================================================
!
! [PURPOSE:] Main program of model-to-observation simulator
!
!=======================================================================
!$use OMP_LIB
  use common, only: r_size
  use common_mpi
  use common_scale
  use common_mpi_scale
  use common_obs_scale
  use common_nml
  use obsope_tools
  use scale_precision, only: RP
  implicit none

  real(RP), allocatable :: v3dg(:,:,:,:)
  real(RP), allocatable :: v2dg(:,:,:)
  real(RP), allocatable :: topog(:,:)
  real(r_size), allocatable :: v3dgh(:,:,:,:)
  real(r_size), allocatable :: v2dgh(:,:,:)

  real(r_size), allocatable :: v3dgsim(:,:,:,:)
  real(r_size), allocatable :: v2dgsim(:,:,:)

  integer :: it

  character(len=7) :: stdoutf = '-000000'
  character(len=6400) :: icmd

#ifdef GPR
  real(RP), allocatable :: rho(:,:,:)
  real(RP), allocatable :: v1dg(:,:)
  real(r_size), allocatable :: v1dgh(:,:)
#endif /* GPR */

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  call initialize_mpi_scale
  call mpi_timer('', 1)

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
  call read_nml_letkf_gpr

  call set_mem_node_proc(1)
  call set_scalelib

  if (myrank_use) then

    call set_common_scale
    call set_common_mpi_scale
    call set_common_obs_scale

    call mpi_timer('INITIALIZE', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Read restart/history input data and run the model-to-obs simulator
!-----------------------------------------------------------------------

    allocate (v3dgh(nlevh, nlonh, nlath, nv3dd))
    allocate (v2dgh(       nlonh, nlath, nv2dd))
    allocate (v3dgsim(nlev, nlon, nlat, OBSSIM_NUM_3D_VARS))
    allocate (v2dgsim(      nlon, nlat, OBSSIM_NUM_2D_VARS))
#ifdef GPR
    allocate (v1dgh(nlevh,               nv1dd))
#endif

    if (trim(OBSSIM_IN_TYPE) == 'restart') then

      allocate (v3dg(nlev,nlon,nlat,nv3d))
      allocate (v2dg(nlon,nlat,nv2d))
      allocate (topog(nlon,nlat))

      call read_restart(OBSSIM_RESTART_IN_BASENAME, v3dg, v2dg)
#ifdef GPR
      allocate (rho(nlev,nlon,nlat))
      rho(:,:,:)=v3dg(:,:,:,iv3d_rho)
#endif

      call state_trans(v3dg)

#ifdef GPR
      v3dg(:,:,:,iv3d_rho)=rho(:,:,:)
      deallocate(rho)
#endif

      call read_topo(OBSSIM_TOPO_IN_BASENAME, topog)
      call state_to_history(v3dg, v2dg, topog, v3dgh, v2dgh)

      deallocate (v3dg, v2dg, topog)

#ifdef GPR
      allocate (v1dg(nlev,nv1d))
      call read_restart_1d(OBSSIM_RESTART_IN_BASENAME, v1dg)
      call state_to_history_1d(v1dg, v1dgh)
      deallocate (v1dg)
#endif

#ifdef GPR
      call obssim_cal(v3dgh, v2dgh, v1dgh, v3dgsim, v2dgsim, stggrd=1)
#else
      call obssim_cal(v3dgh, v2dgh, v3dgsim, v2dgsim, stggrd=1)
#endif

      call write_grd_mpi(OBSSIM_GRADS_OUT_NAME, OBSSIM_NUM_3D_VARS, OBSSIM_NUM_2D_VARS, &
                         1, v3dgsim, v2dgsim)

    else if (trim(OBSSIM_IN_TYPE) == 'history') then

      do it = OBSSIM_TIME_START, OBSSIM_TIME_END
#ifdef GPR
        call read_history(OBSSIM_HISTORY_IN_BASENAME, it, v3dgh, v2dgh)
        call read_history_1d(OBSSIM_HISTORY_IN_BASENAME, 1, v1dgh)
        call obssim_cal(v3dgh, v2dgh, v1dgh, v3dgsim, v2dgsim)
#else
        call read_history(OBSSIM_HISTORY_IN_BASENAME, it, v3dgh, v2dgh)
        call obssim_cal(v3dgh, v2dgh, v3dgsim, v2dgsim)
#endif /* GPR */


        call write_grd_mpi(OBSSIM_GRADS_OUT_NAME, OBSSIM_NUM_3D_VARS, OBSSIM_NUM_2D_VARS, &
                           it, v3dgsim, v2dgsim)
      end do

    else
      write (6, '(2A)') "[Error] Unsupported 'OBSSIM_IN_TYPE': ", trim(OBSSIM_IN_TYPE)
    end if

    call mpi_timer('OBSSIM', 1, barrier=MPI_COMM_a)

    deallocate (v3dgh, v2dgh)
    deallocate (v3dgsim, v2dgsim)
#ifdef GPR
    deallocate (v1dgh)
#endif /* GPR */

    call unset_common_mpi_scale

  end if ! [ myrank_use ]

  call unset_scalelib

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_WORLD)

  call finalize_mpi_scale

  stop
end program obssim
