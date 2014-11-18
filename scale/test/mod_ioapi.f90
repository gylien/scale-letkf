!-------------------------------------------------------------------------------
!> module fileio API
!!
!! @par Description
!!          fileio sample
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_ioapi
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  USE common
  USE common_mpi

  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: IOAPI
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> IO API
  subroutine IOAPI()
    use dc_log, only: &
       LogInit
    use gtool_file, only: &
       FileRead, &
       FileCloseAll

    use gtool_history, only: &
       HistoryInit

    use scale_precision
    use scale_stdio
    use scale_prof
    use scale_grid_index

    use scale_process, only: &
       PRC_setup,    &
       PRC_MPIstart, &
       PRC_MPIfinish, &
       PRC_master, &
       PRC_myrank, &
       PRC_myrank_world, &
       PRC_nu, &
       PRC_2Drank
    use scale_const, only: &
       CONST_setup
    use scale_calendar, only: &
       CALENDAR_setup
    use scale_random, only: &
       RANDOM_setup
    use scale_time, only: &
       TIME_setup
    use scale_grid, only: &
       GRID_setup
    use scale_grid_nest, only: &
       NEST_setup
!    use scale_land_grid_index, only: &
!       LAND_GRID_INDEX_setup
!    use scale_land_grid, only: &
!       LAND_GRID_setup
!    use scale_urban_grid_index, only: &
!       URBAN_GRID_INDEX_setup
!    use scale_urban_grid, only: &
!       URBAN_GRID_setup
    use scale_tracer, only: &
       TRACER_setup
    use scale_fileio, only: &
       FILEIO_setup, &
       FILEIO_write, &
       FILEIO_read
    use scale_history, only: &
       HIST_setup, &
       HIST_get
    use scale_comm, only: &
       COMM_setup, &
       COMM_vars8, &
       COMM_wait


    implicit none

    character(len=H_MID), parameter :: MODELNAME = "SCALE-LES"

    character(len=H_MID) :: DATATYPE = 'DEFAULT' !< REAL4 or REAL8

    real(RP), allocatable :: U(:,:,:), MOMX(:,:,:)

    integer :: k, i, j

    integer               :: dim1_max, dim1_S, dim1_E
    integer               :: dim2_max, dim2_S, dim2_E
    integer               :: dim3_max, dim3_S, dim3_E
!    integer               :: dim4_max, dim4_S, dim4_E
    real(RP), allocatable :: var3D(:,:,:)

    integer :: iolen

    character(len=100) :: basename
    character(len=100) :: varname
    integer :: step

    integer :: rankidx(2)

    basename = 'history'
    varname = 'U'
    step = 1
    !-----------------------------------------------------------------------------

    ! setup standard I/O
    call IO_setup( MODELNAME )

    ! start MPI
    call PRC_MPIstart

    ! setup process
    call PRC_setup

    ! setup Log
    call LogInit(IO_FID_CONF, IO_FID_LOG, IO_L)

    ! setup constants
    call CONST_setup

    ! setup time
    call TIME_setup( setup_TimeIntegration = .false. )

    call PROF_rapstart('Initialize')

    ! setup horizontal/vertical grid coordinates (cartesian,idealized)
    call GRID_INDEX_setup
    call GRID_setup

!    call LAND_GRID_INDEX_setup
!    call LAND_GRID_setup

!    call URBAN_GRID_INDEX_setup
!    call URBAN_GRID_setup

    ! setup file I/O
    call FILEIO_setup

    ! setup mpi communication
    call COMM_setup



    rankidx(1) = PRC_2Drank(PRC_myrank, 1)
    rankidx(2) = PRC_2Drank(PRC_myrank, 2)
    call HistoryInit('','','',IMAX*JMAX*KMAX,PRC_master,PRC_myrank,rankidx)


    call PROF_rapend('Initialize')

    !########## main ##########

    call PROF_rapstart('Main')
    if( IO_L ) write(IO_FID_LOG,*) '*** IOAPI test ***'





!print *, '######', PRC_nu, PRC_myrank_world, PRC_myrank, KA, IA, JA

    if (PRC_nu == 0) then
      basename = trim(basename) // '.u000000'
    else if (PRC_nu == 1) THEN
      basename = trim(basename) // '.u000001'
    end if

    allocate( U(KA,IA,JA) )
    allocate( MOMX(KA,IA,JA) )
    U = 0.0d0




    ! Read file

!    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(varname)

       dim1_max = IMAX !KMAX
       dim2_max = JMAX !IMAX
       dim3_max = KMAX !JMAX
       dim1_S   = IS !KS
       dim1_E   = IE !KE
       dim2_S   = JS !IS
       dim2_E   = JE !IE
       dim3_S   = KS !JS
       dim3_E   = KE !JE

    allocate( var3D(dim1_max,dim2_max,dim3_max) )


    call HIST_get(var3D, trim(basename), trim(varname), step=step)


    forall (i=1:IMAX, j=1:JMAX, k=1:KMAX) U(k+KHALO,i+IHALO,j+JHALO) = var3D(i,j,k)

    deallocate( var3D )


    if (PRC_nu == 0) then
      call FILEIO_read( MOMX(:,:,:),                          & ! [OUT]
                        'init.u000000', 'MOMX', 'ZXY', step=1 ) ! [IN]
    else if (PRC_nu == 1) then
      call FILEIO_read( MOMX(:,:,:),                          & ! [OUT]
                        'init.u000001', 'MOMX', 'ZXY', step=1 ) ! [IN]
    end if


    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j  = JS, JE
    do i  = IS, IE
       U(   1:KS-1,i,j) = U(KS,i,j)
       U(KE+1:KA,  i,j) = U(KE,i,j)
       MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
       MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
    enddo
    enddo


    call COMM_vars8( U   (:,:,:), 1 )
    call COMM_vars8( MOMX(:,:,:), 2 )
    call COMM_wait ( U   (:,:,:), 1 )
    call COMM_wait ( MOMX(:,:,:), 2 )





INQUIRE(IOLENGTH=iolen) iolen
OPEN(myrank+30,FORM='unformatted',ACCESS='direct',RECL=IA*JA*KA*iolen)
WRITE(myrank+30,REC=1) (((real(U(k,i,j),4),i=1,IA),j=1,JA),k=1,KA)
WRITE(myrank+30,REC=2) (((real(MOMX(k,i,j),4),i=1,IA),j=1,JA),k=1,KA)
CLOSE(myrank+30)


!    if( IO_L ) write(IO_FID_LOG,*) '*** IOAPI test end ***'
    call PROF_rapend('Main')

    !########## Finalize ##########

    call PROF_rapreport

    call FileCloseAll

    ! stop MPI
    call PRC_MPIfinish

    return
  end subroutine IOAPI

end module mod_ioapi
