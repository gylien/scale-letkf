program obssim
!=======================================================================
!
! [PURPOSE:] Main program for creating the initial ensemble of an idealized OSE
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

  call read_nml_sc

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

    call  create_ens_mpi(SC_OUT_BASENAME,SC_NATURE_IN_BASENAME)

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
