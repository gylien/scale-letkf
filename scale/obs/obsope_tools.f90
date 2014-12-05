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

!  use common_scalelib

  use common_nml

  use scale_process, only: &
       PRC_myrank, &
       PRC_myrank_world, &
       MPI_COMM_u

  use scale_grid_index, only: &
    KHALO, IHALO, JHALO

  IMPLICIT NONE
  PUBLIC

  !--- PARAM_LETKF_OBSMAKE
  real(r_size) :: OBSERR_U = 1.0d0
  real(r_size) :: OBSERR_V = 1.0d0
  real(r_size) :: OBSERR_T = 1.0d0
  real(r_size) :: OBSERR_Q = 0.001d0
  real(r_size) :: OBSERR_RH = 10.0d0
  real(r_size) :: OBSERR_PS = 100.0d0

!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: nslots=11 ! number of time slots for 4D-LETKF
  INTEGER,PARAMETER :: nbslot=6 ! basetime slot
  REAL(r_size),PARAMETER :: slotint=60.0d0 ! time interval between slots in second

  CHARACTER(7) :: obsfile='obs.dat' !IN
  CHARACTER(21) :: obsdafile='obsda.0000.000000.dat' !OUT


CONTAINS
!-----------------------------------------------------------------------
! Read namelist for obsope
!-----------------------------------------------------------------------
subroutine read_nml_obsope
  implicit none

  call read_nml_letkf_prc
  call read_nml_letkf_obs

  return
end subroutine read_nml_obsope
!-----------------------------------------------------------------------
! Read namelist for obsmake
!-----------------------------------------------------------------------
subroutine read_nml_obsmake
  implicit none

  call read_nml_letkf_prc
  call read_nml_letkf_obs
  call read_nml_letkf_obsmake

  return
end subroutine read_nml_obsmake
!-----------------------------------------------------------------------
! PARAM_LETKF_OBSMAKE
!-----------------------------------------------------------------------
subroutine read_nml_letkf_obsmake
  implicit none
  integer :: ierr

  namelist /PARAM_LETKF_OBSMAKE/ &
    OBSERR_U, &
    OBSERR_V, &
    OBSERR_T, &
    OBSERR_Q, &
    OBSERR_RH, &
    OBSERR_PS

  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_LETKF_OBSMAKE,iostat=ierr)
  if (ierr < 0) then !--- missing
    write(6,*) 'xxx Not found namelist. Check!'
    stop
  elseif (ierr > 0) then !--- fatal error
    write(6,*) 'xxx Not appropriate names in namelist LETKF_PARAM_OBSMAKE. Check!'
    stop
  endif

  return
end subroutine read_nml_letkf_obsmake

!-----------------------------------------------------------------------
! Observation operator calculation
!-----------------------------------------------------------------------
SUBROUTINE obsope_cal(obs)
  IMPLICIT NONE

  TYPE(obs_info),INTENT(IN) :: obs
  type(obs_da_value) :: obsda
  REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)

  integer :: it,islot,proc,im
  integer :: n,nslot,nproc,nprocslot,ierr
  real(r_size) :: rig,rjg,ri,rj,rk
  real(r_size) :: slot_lb,slot_ub

!-----------------------------------------------------------------------

  call set_scalelib(MEM_NP, nitmax, nprocs, proc2mem)

!  call set_mpi_along_domains

  obsda%nobs = obs%nobs
  call obs_da_value_allocate(obsda,0)

  allocate ( v3dg (nlevhalo,nlonhalo,nlathalo,nv3dd) )
  allocate ( v2dg (nlonhalo,nlathalo,nv2dd) )

  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    if (im >= 1 .and. im <= nbv) then
      write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is processing member ', &
            im, ', subdomain id #', proc2mem(2,it,myrank+1)

      nproc = 0
      do islot = SLOT_START, SLOT_END
        slot_lb = (real(islot-SLOT_BASE,r_size) - 0.5d0) * SLOT_TINTERVAL
        slot_ub = (real(islot-SLOT_BASE,r_size) + 0.5d0) * SLOT_TINTERVAL
        write (6,'(A,I3,A,F7.1,A,F7.1,A)') 'Slot #', islot-SLOT_START+1, ': time interval (', slot_lb, ',', slot_ub, '] sec'

        call read_ens_history_mpi('hist',it,islot,v3dg,v2dg)
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        nslot = 0
        nprocslot = 0
        do n = 1, obs%nobs

          if (obs%dif(n) > slot_lb .and. obs%dif(n) <= slot_ub) then
            nslot = nslot + 1

            call phys2ij(obs%lon(n),obs%lat(n),rig,rjg)
            call ijproc(rig,rjg,ri,rj,proc)

  !          if (PRC_myrank == 0) then
  !            print *, proc, rig, rjg, ri, rj
  !          end if

            if (PRC_myrank == proc) then
              nproc = nproc + 1
              nprocslot = nprocslot + 1
              obsda%idx(nproc) = n
              obsda%ri(nproc) = rig
              obsda%rj(nproc) = rjg

              call phys2ijk(v3dg(:,:,:,iv3dd_p),obs%elm(n),ri,rj,obs%lev(n),rk)

              if (rk == -1.0d0) then
                obsda%qc(nproc) = iqc_out_vhi
              else if (rk == -2.0d0) then
                obsda%qc(nproc) = iqc_out_vlo
              else if (rk == -3.0d0) then
                obsda%qc(nproc) = iqc_out_h
              else
                call Trans_XtoY(obs%elm(n),ri,rj,rk,v3dg,v2dg,obsda%val(nproc),obsda%qc(nproc))
              end if

            end if ! [ PRC_myrank == proc ]

          end if ! [ obs%dif(n) > slot_lb .and. obs%dif(n) <= slot_ub ]

        end do ! [ n = 1, obs%nobs ]

        write (6,'(A,I10)') ' -- nobs in the slot = ', nslot
        write (6,'(A,I6,A,I10)') ' -- nobs in the slot and processed by rank ', myrank, ' = ', nprocslot

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

      obsda%nobs = nproc

      write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' finishes processing member ', &
            im, ', subdomain id #', proc2mem(2,it,myrank+1)
      write (6,'(A,I8,A)') ' -- ', nproc, ' observations found'

      write (obsdafile(7:10),'(I4.4)') im
      write (obsdafile(12:17),'(I6.6)') proc2mem(2,it,myrank+1)
      call write_obs_da(obsdafile,obsda,0)

    end if

  end do ! [ it = 1, nitmax ]

  deallocate ( v3dg, v2dg )

  call unset_scalelib

end subroutine obsope_cal
!-----------------------------------------------------------------------
! Observation generator calculation
!-----------------------------------------------------------------------

!!! need to refine qc...

SUBROUTINE obsmake_cal(obs)
  IMPLICIT NONE

  TYPE(obs_info),INTENT(INOUT) :: obs
  REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)
  REAL(r_size),ALLOCATABLE :: error(:)

  integer :: islot,proc
  integer :: n,nslot,nproc,nprocslot,ierr,iqc
  real(r_size) :: rig,rjg,ri,rj,rk
  real(r_size) :: slot_lb,slot_ub
  real(r_size),allocatable :: bufr(:)

  CHARACTER(10) :: obsoutfile = 'obsout.dat'

!-----------------------------------------------------------------------

  call set_scalelib(MEM_NP, nitmax, nprocs, proc2mem)

  allocate ( v3dg (nlevhalo,nlonhalo,nlathalo,nv3dd) )
  allocate ( v2dg (nlonhalo,nlathalo,nv2dd) )

  write (6,'(A,I6.6,A,I6.6)') 'MYRANK ',myrank,' is processing subdomain id #', proc2mem(2,1,myrank+1)

  nproc = 0
  do islot = SLOT_START, SLOT_END
    slot_lb = (real(islot-SLOT_BASE,r_size) - 0.5d0) * SLOT_TINTERVAL
    slot_ub = (real(islot-SLOT_BASE,r_size) + 0.5d0) * SLOT_TINTERVAL
    write (6,'(A,I3,A,F7.1,A,F7.1,A)') 'Slot #', islot-SLOT_START+1, ': time interval (', slot_lb, ',', slot_ub, '] sec'

    call read_ens_history_mpi('hist',1,islot,v3dg,v2dg)

    nslot = 0
    nprocslot = 0
    do n = 1, obs%nobs

      if (obs%dif(n) > slot_lb .and. obs%dif(n) <= slot_ub) then
        nslot = nslot + 1

        call phys2ij(obs%lon(n),obs%lat(n),rig,rjg)
        call ijproc(rig,rjg,ri,rj,proc)

!          if (PRC_myrank == 0) then
!            print *, proc, rig, rjg, ri, rj
!          end if

        if (PRC_myrank == proc) then
          nproc = nproc + 1
          nprocslot = nprocslot + 1

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

          call phys2ijk(v3dg(:,:,:,iv3dd_p),obs%elm(n),ri,rj,obs%lev(n),rk)

          call Trans_XtoY(obs%elm(n),ri,rj,rk,v3dg,v2dg,obs%dat(n),iqc)
        end if ! [ PRC_myrank == proc ]

      end if ! [ obs%dif(n) > slot_lb .and. obs%dif(n) <= slot_ub ]

    end do ! [ n = 1, obs%nobs ]

    write (6,'(A,I10)') ' -- nobs in the slot = ', nslot
    write (6,'(A,I6,A,I10)') ' -- nobs in the slot and processed by rank ', myrank, ' = ', nprocslot

  end do ! [ islot = SLOT_START, SLOT_END ]

  deallocate ( v3dg, v2dg )

  if (PRC_myrank == 0) then
    allocate ( bufr (obs%nobs) )
  end if

  CALL MPI_REDUCE(obs%dat,bufr,obs%nobs,MPI_r_size,MPI_MAX,0,MPI_COMM_u,ierr)

  if (PRC_myrank == 0) then

    obs%dat = bufr
    deallocate ( bufr )

    ALLOCATE(error(obs%nobs))
    CALL com_randn(obs%nobs,error)
    do n = 1, obs%nobs
      select case(obs%elm(n))
      case(id_u_obs)
        obs%err(n) = OBSERR_U
      case(id_v_obs)
        obs%err(n) = OBSERR_V
      case(id_t_obs,id_tv_obs)
        obs%err(n) = OBSERR_T
      case(id_q_obs)
        obs%err(n) = OBSERR_Q
      case(id_rh_obs)
        obs%err(n) = OBSERR_RH
      case(id_ps_obs)
        obs%err(n) = OBSERR_PS
      case default
        write(6,'(A)') 'warning: skip assigning observation error (unsupported observation type)' 
      end select
      obs%dat(n) = obs%dat(n) + obs%err(n) * error(n)

print *, '######', obs%elm(n), obs%dat(n)

    end do ! [ n = 1, obs%nobs ]

    call write_obs(obsoutfile,obs)

  end if ! [ PRC_myrank == 0 ]

  call unset_scalelib

end subroutine obsmake_cal
!=======================================================================

END MODULE obsope_tools
