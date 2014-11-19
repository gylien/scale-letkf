MODULE obsope_tools
!=======================================================================
!
! [PURPOSE:] Observation operator tools
!
! [HISTORY:]
!   November 2014  Guo-Yuan Lien  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
  USE common_obs_scale

  use common_letkf, only: nbv

  use letkf_nml

  use scale_process, only: &
       prc_myrank, &
       prc_myrank_world

  use scale_grid_index, only: &
    KHALO, IHALO, JHALO

  IMPLICIT NONE
  PUBLIC

!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: nslots=11 ! number of time slots for 4D-LETKF
  INTEGER,PARAMETER :: nbslot=6 ! basetime slot
  REAL(r_size),PARAMETER :: slotint=60.0d0 ! time interval between slots in second

  CHARACTER(7) :: obsfile='obs.dat' !IN
  CHARACTER(15) :: obsvalfile='obsval.0000.dat' !OUT


CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE obsope_cal
  IMPLICIT NONE



  REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)

  integer :: it, islot

  integer :: ierr
!integer :: i, j, k


!-----------------------------------------------------------------------


  if (scale_IO_group_n >= 0) then
    call set_scale_IO(MEM_NP)

!  print *, myrank, PRC_myrank


    ALLOCATE ( v3dg (nlev+2*KHALO,nlon+2*IHALO,nlat+2*JHALO,nv3dd) )
    ALLOCATE ( v2dg (nlon+2*IHALO,nlat+2*JHALO,nv2dd) )

    do it = 1, nitmax

      write (6,'(A,I)') '# Loop ', it

      do islot = 1, SLOT_NUM

        write (6,'(A,I)') '## Slot loop ', islot

!    WRITE(obsinfile(6:7),'(I2.2)') islot
!    CALL read_obs(obsinfile,nobslots(islot),&
!      & elem(1:nobslots(islot)),rlon(1:nobslots(islot)),&
!      & rlat(1:nobslots(islot)),rlev(1:nobslots(islot)),&
!      & odat(1:nobslots(islot)),oerr(1:nobslots(islot)),&
!      & otyp(1:nobslots(islot)))
!    tdif(1:nobslots(islot)) = REAL(islot-nbslot,r_size)*hourslot


        call read_ens_history_mpi('hist',it,islot,v3dg,v2dg)


!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!if (myrank ==0) then
!print *, v3dg(5,:,6,iv3dd_t)
!end if

!if (myrank<4) then
! OPEN(myrank+30,FORM='unformatted',ACCESS='direct',RECL=IA*JA*KA)
! WRITE(myrank+30,REC=1) (((real(v3dg(k,i,j,iv3dd_p),4),i=1,IA),j=1,JA),k=1,KA)
! WRITE(myrank+30,REC=2) (((real(v3dg(k,i,j,iv3dd_hgt),4),i=1,IA),j=1,JA),k=1,KA)
! CLOSE(myrank+30)
!end if




!    ohx=0.0d0
!    oqc=0

!    DO n=1,nobslots(islot)
!      CALL phys2ijk(v3d(:,:,:,iv3d_p),elem(n),rlon(n),rlat(n),rlev(n),ri,rj,rk)
!      IF(CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) THEN
!!        WRITE(6,'(A)') '* X-coordinate out of range'
!!        WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', rlon=',rlon(n)
!        CYCLE
!      END IF
!      IF(CEILING(rj) < 2 .OR. nlat < CEILING(rj)) THEN
!!        WRITE(6,'(A)') '* Y-coordinate out of range'
!!        WRITE(6,'(A,F6.2,A,F6.2)') '*   rj=',rj,', rlat=',rlat(n)
!        CYCLE
!      END IF
!      IF(CEILING(rk) > nlev) THEN
!!        CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,dz)
!!        WRITE(6,'(A)') '* Z-coordinate out of range'
!!        WRITE(6,'(A,F6.2,A,F10.2,A,F6.2,A,F6.2,A,F10.2)') &
!!         & '*   rk=',rk,', rlev=',rlev(n),&
!!         & ', (lon,lat)=(',rlon(n),',',rlat(n),'), phi0=',dz
!        CYCLE
!      END IF
!      IF(CEILING(rk) < 2 .AND. NINT(elem(n)) /= id_ps_obs) THEN
!        IF(NINT(elem(n)) > 9999) THEN
!          rk = 0.0d0
!        ELSE IF(NINT(elem(n)) == id_u_obs .OR. NINT(elem(n)) == id_v_obs) THEN
!          rk = 1.00001d0
!        ELSE
!!          CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,dz)
!!          WRITE(6,'(A)') '* Z-coordinate out of range'
!!          WRITE(6,'(A,F6.2,A,F10.2,A,F6.2,A,F6.2,A,F10.2)') &
!!           & '*   rk=',rk,', rlev=',rlev(n),&
!!           & ', (lon,lat)=(',rlon(n),',',rlat(n),'), phi0=',dz
!          CYCLE
!        END IF
!      END IF
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
!      !
!      ! observational operator
!      !
!      CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,ohx(n))
!      oqc(n) = 1
!    END DO
!    IF(firstwrite) THEN
!      CALL write_obs2(obsoutfile,nobslots(islot),&
!        & elem(1:nobslots(islot)),rlon(1:nobslots(islot)),&
!        & rlat(1:nobslots(islot)),rlev(1:nobslots(islot)),&
!        & odat(1:nobslots(islot)),oerr(1:nobslots(islot)),&
!        & otyp(1:nobslots(islot)),tdif(1:nobslots(islot)),&
!        & ohx(1:nobslots(islot)),oqc(1:nobslots(islot)),0)
!      firstwrite = .FALSE.
!    ELSE
!      CALL write_obs2(obsoutfile,nobslots(islot),&
!        & elem(1:nobslots(islot)),rlon(1:nobslots(islot)),&
!        & rlat(1:nobslots(islot)),rlev(1:nobslots(islot)),&
!        & odat(1:nobslots(islot)),oerr(1:nobslots(islot)),&
!        & otyp(1:nobslots(islot)),tdif(1:nobslots(islot)),&
!        & ohx(1:nobslots(islot)),oqc(1:nobslots(islot)),1)
!    END IF
!  END DO

!  IF(nobslots(nslots+1) > 0) THEN
!    v2d(:,:,iv2d_tprcp) = pp
!    CALL read_obs(obsinfile_mean,nobslots(nslots+1),&
!      & elem(1:nobslots(nslots+1)),rlon(1:nobslots(nslots+1)),&
!      & rlat(1:nobslots(nslots+1)),rlev(1:nobslots(nslots+1)),&
!      & odat(1:nobslots(nslots+1)),oerr(1:nobslots(nslots+1)),&
!      & otyp(1:nobslots(nslots+1)))
!    tdif(1:nobslots(nslots+1)) = 0.0d0
!    ohx=0.0d0
!    oqc=0
!    DO n=1,nobslots(nslots+1)
!      IF(elem(n) /= id_rain_obs) CYCLE
!      CALL phys2ij(rlon(n),rlat(n),ri,rj)
!      IF(CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) THEN
!!        WRITE(6,'(A)') '* X-coordinate out of range'
!!        WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', rlon=',rlon(n)
!        CYCLE
!      END IF
!      IF(CEILING(rj) < 2 .OR. nlat < CEILING(rj)) THEN
!!        WRITE(6,'(A)') '* Y-coordinate out of range'
!!       WRITE(6,'(A,F6.2,A,F6.2)') '*   rj=',rj,', rlat=',rlat(n)
!        CYCLE
!      END IF
!      !
!      ! observational operator
!      !
!      rk = 0.0d0
!      CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,ohx(n))
!      oqc(n) = 1
!    END DO
!    IF(firstwrite) THEN
!      CALL write_obs2(obsoutfile,nobslots(nslots+1),&
!        & elem(1:nobslots(nslots+1)),rlon(1:nobslots(nslots+1)),&
!        & rlat(1:nobslots(nslots+1)),rlev(1:nobslots(nslots+1)),&
!        & odat(1:nobslots(nslots+1)),oerr(1:nobslots(nslots+1)),&
!        & otyp(1:nobslots(nslots+1)),tdif(1:nobslots(nslots+1)),&
!        & ohx(1:nobslots(nslots+1)),oqc(1:nobslots(nslots+1)),0)
!      firstwrite = .FALSE.
!    ELSE
!      CALL write_obs2(obsoutfile,nobslots(nslots+1),&
!        & elem(1:nobslots(nslots+1)),rlon(1:nobslots(nslots+1)),&
!        & rlat(1:nobslots(nslots+1)),rlev(1:nobslots(nslots+1)),&
!        & odat(1:nobslots(nslots+1)),oerr(1:nobslots(nslots+1)),&
!        & otyp(1:nobslots(nslots+1)),tdif(1:nobslots(nslots+1)),&
!        & ohx(1:nobslots(nslots+1)),oqc(1:nobslots(nslots+1)),1)
!    END IF
!  END IF

!  DEALLOCATE( elem,rlon,rlat,rlev,odat,oerr,otyp,tdif,ohx,oqc )

      end do
    end do



    call unset_scale_IO
  end if ! [scale_IO_group_n >= 0]



!  REAL(r_size) :: v3d(nlon,nlat,nlev,nv3dx)
!  REAL(r_size) :: v2d(nlon,nlat,nv2dx)
!  REAL(r_size) :: pp(nlon,nlat)
!  REAL(r_size),PARAMETER :: threshold_dz=500.0d0
!  REAL(r_size) :: dz,tg,qg
!  INTEGER :: nobslots(nslots+1)
!  INTEGER :: nobs,maxnobs
!  REAL(r_size) :: ri,rj,rk
!  INTEGER :: n,islot
!  LOGICAL :: firstwrite

!  CALL set_common_gfs

!  DO islot=1,nslots
!    WRITE(obsinfile(6:7),'(I2.2)') islot
!    CALL get_nobs(obsinfile,7,nobslots(islot))
!    WRITE(6,'(2A,I9,A)') obsinfile, ':', nobslots(islot), ' OBSERVATIONS'
!  END DO
!  CALL get_nobs(obsinfile_mean,7,nobslots(nslots+1))
!  WRITE(6,'(2A,I9,A)') obsinfile_mean, ':', nobslots(nslots+1), ' OBSERVATIONS'
!  nobs = SUM(nobslots)
!  maxnobs = MAXVAL(nobslots)
!  WRITE(6,'(A,I9,A)') 'TOTAL:      ', nobs, ' OBSERVATIONS'
!  ALLOCATE( elem(maxnobs) )
!  ALLOCATE( rlon(maxnobs) )
!  ALLOCATE( rlat(maxnobs) )
!  ALLOCATE( rlev(maxnobs) )
!  ALLOCATE( odat(maxnobs) )
!  ALLOCATE( oerr(maxnobs) )
!  ALLOCATE( otyp(maxnobs) )
!  ALLOCATE( tdif(maxnobs) )
!  ALLOCATE( ohx(maxnobs) )
!  ALLOCATE( oqc(maxnobs) )
!  firstwrite = .TRUE.
!  pp = 0.0d0


end subroutine obsope_cal

END MODULE obsope_tools
