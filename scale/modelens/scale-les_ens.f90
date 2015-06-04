!-------------------------------------------------------------------------------
!> Program SCALE-LES (a launcher of main routine)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
program scaleles_launcher
  !-----------------------------------------------------------------------------


  USE common
  USE common_mpi
  USE common_scale
  USE common_mpi_scale

  use common_nml


  !
  !++ used modules
  !
  use dc_log, only: &
     LogInit
  use gtool_file, only: &
     FileCloseAll
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_process, only: &
     PRC_setup,         &
     PRC_MPIstart,      &
     PRC_MPIsplit,      &
     PRC_MPIsplit_letkf,&
     PRC_MPIfinish,     &
     PRC_MPIstop,       &
     PRC_BULKsetup,     &
     MASTER_LOG,        &
     MASTER_COMM_WORLD, &
     MASTER_nmax,       &
     max_depth, &
     PRC_NUM_X, &
     PRC_NUM_Y, &
     LOCAL_COMM_WORLD
  use scale_fileio, only: &
     FILEIO_setup
  use scale_comm, only: &
     COMM_setup
  use mod_les_driver
  !
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

!  integer :: NUM_BULKJOB                  ! number of bulk jobs
!  integer :: NUM_DOMAIN                   ! number of domains
!  integer :: PRC_BULKJOB(max_depth)       ! # of total process in each bulk job
!  integer :: PRC_DOMAINS(max_depth)       ! # of total process in each domain
!  integer :: BULK_COMM_WORLD              ! split communicator for each bulk job
!  integer :: MY_COMM_WORLD                ! assigned local communicator
!  integer :: inter_parent                 ! inter communicator with parent
!  integer :: inter_child                  ! inter communicator with child

!  logical :: ABORT_ALL_JOBS = .false.     ! flag of abort all jobs or not
!  logical :: LOG_SPLIT = .false.          ! flag of log-output for mpi splitting

!  character(len=H_LONG) :: CONF_FILES(max_depth)  ! names of configulation files
!  character(len=H_LONG) :: CONF_BULKS(max_depth)  ! names of configulation files (dummy)
!  character(len=H_LONG) :: fname_launch           ! config file for launcher
!  character(len=H_LONG) :: fname_local            ! config file for local domain
!  character(len=H_LONG) :: fname_bulks            ! config file for dummy use
!  character(len=H_LONG) :: fname_scaleles         ! config file for scaleles

!  integer :: check, nprocs
!  integer :: ierr
!  integer :: LNC_FID_CONF


!  character(len=H_MID), parameter :: MODELNAME = "SCALE-LES"



  REAL(r_dble) :: rtimer00,rtimer
  INTEGER :: ierr, it, im
  CHARACTER(11) :: stdoutf='NOUT-000000'
  CHARACTER(11) :: timer_fmt='(A30,F10.2)'


  CHARACTER(len=H_LONG) :: confname='0000/run.conf'

  integer :: LOCAL_myrank, LOCAL_nmax

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  CALL initialize_mpi
  rtimer00 = MPI_WTIME()

  WRITE(stdoutf(6:11), '(I6.6)') myrank
  WRITE(6,'(3A,I6.6)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I6.6,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf

!-----------------------------------------------------------------------

  ! setup standard I/O
  call IO_setup( MODELNAME, .false.)

  call read_nml_letkf
  call read_nml_letkf_prc
  call read_nml_letkf_obs !!!!!!!!!!!!!!!!!!!!!! move outside of subroutine????
  call read_nml_letkf_obserr !!!!!!!!!!!!!!!!!!! move outside of subroutine????
  call read_nml_letkf_obs_radar !!!!!!!!!!!!!!!! move outside of subroutine????

  if (nprocs /= NNODES * PPN) then
    write(6,'(A,I10)') 'Number of MPI processes = ', nprocs
    write(6,'(A,I10)') 'NNODES = ', NNODES
    write(6,'(A,I10)') 'PPN    = ', PPN
    write(6,'(A)') 'Number of MPI processes should be equal to NNODES * PPN.'
    stop
  end if

  !
  ! Set up node and process distribution
  !
  call set_mem_node_proc(MEMBER+1,NNODES,PPN,MEM_NODES,MEM_NP)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  rtimer = MPI_WTIME()
  WRITE(6,timer_fmt) '### TIMER(INITIALIZE):        ',rtimer-rtimer00
  rtimer00=rtimer

!-----------------------------------------------------------------------
! Run SCALE-LES
!-----------------------------------------------------------------------

  if (scale_IO_use) then

!    if (MEM_NP /= PRC_NUM_X * PRC_NUM_Y) then
!      write(6,'(A,I10)') 'MEM_NP    = ', MEM_NP
!      write(6,'(A,I10)') 'PRC_NUM_X = ', PRC_NUM_X
!      write(6,'(A,I10)') 'PRC_NUM_Y = ', PRC_NUM_Y
!      write(6,'(A)') 'MEM_NP should be equal to PRC_NUM_X * PRC_NUM_Y.'
!      stop
!    end if


  !-----------------------------------------------------------------------------

    ! start SCALE MPI
    call PRC_MPIstart

    ! split MPI communicator for LETKF
    call PRC_MPIsplit_letkf(MEM_NP, nitmax, nprocs, proc2mem, myrank, &
                            LOCAL_myrank, LOCAL_nmax)

!  ! setup MPI
!  call PRC_MPIsetup( LOCAL_COMM_WORLD )


    do it = 1, nitmax
      im = proc2mem(1,it,myrank+1)
      if (im >= 1 .and. im <= MEMBER+1) then



        if (im == MEMBER+1) then
          WRITE(confname(1:4),'(A4)') 'mean'
        else
          WRITE(confname(1:4),'(I4.4)') proc2mem(1,it,myrank+1)
        end if

        WRITE(6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is running a model with configuration file: ', confname

        call scaleles ( LOCAL_COMM_WORLD,   &
                        MPI_COMM_NULL,    &
                        MPI_COMM_NULL,     &
                        confname )

      end if


!      CALL MPI_BARRIER(MPI_COMM_a,ierr)

    end do ! [ it = 1, nitmax ]






  else

    write (6, '(A,I6.6,A)') 'MYRANK=',myrank,': This process is not used!'

  end if

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  rtimer = MPI_WTIME()
  WRITE(6,timer_fmt) '### TIMER(SCALE_LES):          ',rtimer-rtimer00
  rtimer00=rtimer


!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  CALL finalize_mpi




  stop
  !=============================================================================
  !
end program scaleles_launcher
