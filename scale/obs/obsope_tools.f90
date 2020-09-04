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

!  use common_scalelib

  use common_nml

!  use scale_process, only: &
!       PRC_myrank
!       MPI_COMM_d => LOCAL_COMM_WORLD

  use scale_grid_index, only: &
    KHALO, IHALO, JHALO
#ifdef H08
  use scale_grid, only: &
      DX, DY,           &
      BUFFER_DX,        &
      BUFFER_DY
#endif


  IMPLICIT NONE
  PUBLIC

!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------



CONTAINS
!!-----------------------------------------------------------------------
!! Read namelist for obsope
!!-----------------------------------------------------------------------
!subroutine read_nml_obsope
!  implicit none

!  call read_nml_letkf
!  call read_nml_letkf_prc

!  return
!end subroutine read_nml_obsope
!!-----------------------------------------------------------------------
!! Read namelist for obsmake
!!-----------------------------------------------------------------------
!subroutine read_nml_obsmake
!  implicit none

!  call read_nml_letkf
!  call read_nml_letkf_prc
!  call read_nml_letkf_obsmake

!  return
!end subroutine read_nml_obsmake
!-----------------------------------------------------------------------
! PARAM_LETKF_OBSMAKE
!-----------------------------------------------------------------------
!subroutine read_nml_letkf_obsmake
!  implicit none
!  integer :: ierr

!  namelist /PARAM_LETKF_OBSMAKE/ &
!    OBSERR_U, &
!    OBSERR_V, &
!    OBSERR_T, &
!    OBSERR_Q, &
!    OBSERR_RH, &
!    OBSERR_PS, &
!    OBSERR_RADAR_REF, &
!    OBSERR_RADAR_VR

!  rewind(IO_FID_CONF)
!  read(IO_FID_CONF,nml=PARAM_LETKF_OBSMAKE,iostat=ierr)
!  if (ierr < 0) then !--- missing
!    write(6,*) 'xxx Not found namelist. Check!'
!    stop
!  elseif (ierr > 0) then !--- fatal error
!    write(6,*) 'xxx Not appropriate names in namelist LETKF_PARAM_OBSMAKE. Check!'
!    stop
!  endif

!  return
!end subroutine read_nml_letkf_obsmake

!-----------------------------------------------------------------------
! Observation operator calculation
!-----------------------------------------------------------------------
SUBROUTINE obsope_cal(obs, obsda_return)
  IMPLICIT NONE

  TYPE(obs_info),INTENT(IN) :: obs(OBS_IN_NUM)
  type(obs_da_value),OPTIONAL,INTENT(INOUT) :: obsda_return
  type(obs_da_value) :: obsda
  REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)

  integer :: it,islot,proc,im,iof
  integer :: n,nn,nslot,nobs,nobs_0,nobs_slot,nobs_alldomain
!  real(r_size) :: rig,rjg,ri,rj,rk
  real(r_size) :: rig,rjg,rk
  real(r_size),allocatable :: ri(:),rj(:)
  real(r_size) :: ritmp,rjtmp
  real(r_size) :: slot_lb,slot_ub

#ifdef H08
! -- for Himawari-8 obs --
  real(r_size) :: yobs_H08(nlon,nlat,NIRB_HIM8)
  real(r_size) :: yobs_H08_clr(nlon,nlat,NIRB_HIM8)
  integer :: qc_H08(nlon,nlat,NIRB_HIM8)
  integer :: i8, j8
  
  real(r_size) :: plev_obs_H08(nlon,nlat,NIRB_HIM8)
!
#endif

! -- for TC vital assimilation --
  INTEGER :: obs_idx_TCX, obs_idx_TCY, obs_idx_TCP ! obs index
  INTEGER :: bTC_proc ! the process where the background TC is located.
! bTC: background TC in each subdomain
! bTC(1,:) : tcx (m), bTC(2,:): tcy (m), bTC(3,:): mslp (Pa)
  REAL(r_size),ALLOCATABLE :: bTC(:,:)
  REAL(r_size) :: bTC_mslp

  character(filelenmax) :: obsdafile
  character(11) :: obsda_suffix = '.000000.dat'

!-----------------------------------------------------------------------

  integer :: ierr
  REAL(r_dble) :: rrtimer00,rrtimer

  logical :: USE_HIM8 = .false. ! initialize

  real(r_size) :: fp_max, fp_max_lon, fp_max_lat

!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer00 = MPI_WTIME()


  nobs_alldomain = 0
  do iof = 1, OBS_IN_NUM
    if (OBSDA_RUN(iof)) then
      nobs_alldomain = nobs_alldomain + obs(iof)%nobs
    end if
  end do
  obsda%nobs = nobs_alldomain
  call obs_da_value_allocate(obsda,0)
  allocate ( ri(nobs_alldomain) )
  allocate ( rj(nobs_alldomain) )

  allocate ( v3dg (nlevh,nlonh,nlath,nv3dd) )
  allocate ( v2dg (nlonh,nlath,nv2dd) )

  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    if (im >= 1 .and. im <= MEMBER) then
      write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is processing member ', &
            im, ', subdomain id #', proc2mem(2,it,myrank+1)


!write(6,*) '%%%%%%', MPI_WTIME(), 0

      nobs = 0

      !!!!!!
      if (nobs_alldomain > 0) then
      !!!!!!

      obsda%qc = iqc_time

      do islot = SLOT_START, SLOT_END
        slot_lb = (real(islot-SLOT_BASE,r_size) - 0.5d0) * SLOT_TINTERVAL
        slot_ub = (real(islot-SLOT_BASE,r_size) + 0.5d0) * SLOT_TINTERVAL
        write (6,'(A,I3,A,F9.1,A,F9.1,A)') 'Slot #', islot-SLOT_START+1, ': time interval (', slot_lb, ',', slot_ub, '] sec'

        call read_ens_history_iter(HISTORY_IN_BASENAME,it,islot,v3dg,v2dg)
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,I3,A,I3,A,4x,F15.7)') '###### obsope_cal:read_ens_history_iter:',it,':',islot,':',rrtimer-rrtimer00
  rrtimer00=rrtimer



        if ( LT_2DLOC ) then
          call get_fpmax_loc2d( v3dg, fp_max, fp_max_lon, fp_max_lat ) 
        endif

        do iof = 1, OBS_IN_NUM

          if (.not. OBSDA_RUN(iof)) cycle

          obs_idx_TCX = -1
          obs_idx_TCY = -1
          obs_idx_TCP = -1

          nslot = 0
          nobs_slot = 0

!write(6,*) '%%%===', MPI_WTIME(), 'im:', im, 'islot:', islot, 'iof:', iof


            ! do this small computation first, without OpenMP
            nobs_0 = nobs
            do n = 1, obs(iof)%nobs

              select case (obs(iof)%elm(n))
              case (id_tclon_obs)
                obs_idx_TCX = n
                cycle
              case (id_tclat_obs)
                obs_idx_TCY = n
                cycle
              case (id_tcmip_obs)
                obs_idx_TCP = n
                cycle
              end select

              if (obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub) then
                nslot = nslot + 1
                call phys2ij(obs(iof)%lon(n),obs(iof)%lat(n),rig,rjg)
                call rij_g2l_auto(proc,rig,rjg,ritmp,rjtmp)

                if (myrank_d == proc) then
                  nobs = nobs + 1
                  nobs_slot = nobs_slot + 1
                  obsda%set(nobs) = iof
                  obsda%idx(nobs) = n
                  obsda%ri(nobs) = rig  ! obsda%ri: global grid coordinate
                  obsda%rj(nobs) = rjg  !
                  ri(nobs) = ritmp      ! ri: local grid coordinate
                  rj(nobs) = rjtmp      !
                end if ! [ myrank_d == proc ]
              end if ! [ obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub ]
            end do ! [ n = 1, obs%nobs ]




!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,I3,A,I3,A,I3,A,F15.7)') '###### obsope_cal:obsope_step_1:        ',it,':',islot,':',iof,':',rrtimer-rrtimer00
  rrtimer00=rrtimer


          ! then do this heavy computation with OpenMP

          IF(OBS_IN_FORMAT(iof) /= 3)THEN ! H08


!write(6,*) '%%%===', MPI_WTIME(), nobs_0 + 1, nobs

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(nn,n,rk)
            do nn = nobs_0 + 1, nobs
              n = obsda%idx(nn)

!if (mod(nn,50) == 0) then
!  write(6,*) '%%%%%%', MPI_WTIME(), nn
!end if

              if (.not. USE_OBS(obs(iof)%typ(n))) then
                obsda%qc(nn) = iqc_otype
                cycle
              end if

              if (obs(iof)%elm(n) == id_radar_ref_obs .or. obs(iof)%elm(n) == id_radar_ref_zero_obs .or. obs(iof)%elm(n) == id_radar_vr_obs) then
                if (obs(iof)%lev(n) > RADAR_ZMAX) then
                  obsda%qc(nn) = iqc_radar_vhi
#ifdef DEBUG
                  write(6,'(A,F8.1,A,I5)') 'warning: radar observation is too high: lev=', obs(iof)%lev(n), ', elem=', obs(iof)%elm(n)
#endif
                else
                  call phys2ijkz(v3dg(:,:,:,iv3dd_hgt),ri(nn),rj(nn),obs(iof)%lev(n),rk,obsda%qc(nn))
                end if
              else
                if (obs(iof)%lev(n) > LT_ZMAX .and. id_fp3d_obs) then
                  obsda%qc(nn) = iqc_radar_vhi
                else
                  call phys2ijkz(v3dg(:,:,:,iv3dd_hgt),ri(nn),rj(nn),obs(iof)%lev(n),rk,obsda%qc(nn)) ! z level 10/24/2018 by TH
                endif
              end if


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,I3,A,I3,A,I3,A,I8,A,F15.7)') '###### obsope_cal:obsope_step_2_phys2ijkz:',it,':',islot,':',iof,':',nn,':',rrtimer-rrtimer00
!  rrtimer00=rrtimer


              if (obsda%qc(nn) == iqc_good) then
                select case (OBS_IN_FORMAT(iof))
                case (1)
                  call Trans_XtoY(obs(iof)%elm(n),ri(nn),rj(nn),rk, &
                                  obs(iof)%lon(n),obs(iof)%lat(n),v3dg,v2dg,obsda%val(nn),obsda%qc(nn))
                case (2, 5)

                  call Trans_XtoY_radar(obs(iof)%elm(n),obs(iof)%meta(1),obs(iof)%meta(2),obs(iof)%meta(3),ri(nn),rj(nn),rk, &
                                        obs(iof)%lon(n),obs(iof)%lat(n),obs(iof)%lev(n),v3dg,v2dg,obsda%val(nn),obsda%qc(nn))
                  if (obsda%qc(nn) == iqc_ref_low) obsda%qc(nn) = iqc_good ! when process the observation operator, we don't care if reflectivity is too small

                !!!!!! may not need to do this at this stage...
                !if (obs(iof)%elm(n) == id_radar_ref_obs) then
                !  obsda%val(nn) = 10.0d0 * log10(obsda%val(nn))
                !end if
                !!!!!!

                case (4, 6) ! Lighting obs
                  if ( LT_2DLOC ) then
                    obsda%lev(nn) = 0.0d0
                    select case( obs(iof)%elm(n) )
                    case( id_fp2d_obs_max ) 
                      obsda%val(nn) = fp_max
                      obsda%qc(nn) = iqc_good
                    case( id_fp2d_obs_lon ) 
                      obsda%val(nn) = fp_max_lon
                      obsda%qc(nn) = iqc_good
                    case( id_fp2d_obs_lat ) 
                      obsda%val(nn) = fp_max_lat
                      obsda%qc(nn) = iqc_good
                    case default
                      obsda%qc(nn) = iqc_obs_bad
                    end select
                  else
                    call Trans_XtoY_LT( obs(iof)%elm(n), ri(nn), rj(nn), rk, &
                                        v3dg, v2dg, obsda%val(nn), obsda%qc(nn), obsda%lev(nn), myrank_d )
                  endif

                end select
              end if


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,I3,A,I3,A,I3,A,I8,A,F15.7)') '###### obsope_cal:obsope_step_2_Trans_XtoY_radar:',it,':',islot,':',iof,':',nn,':',rrtimer-rrtimer00
!  rrtimer00=rrtimer


            end do ! [ nn = nobs_0 + 1, nobs ]
!$OMP END PARALLEL DO

#ifdef H08
          ELSEIF((OBS_IN_FORMAT(iof) == 3) )THEN ! H08

            call Trans_XtoY_H08_allg(v3dg, v2dg, yobs_H08, yobs_H08_clr,&
                                     qc_H08)

            do nn = nobs_0 + 1, nobs
              n = obsda%idx(nn)
              i8 = nint( ri(nn) )
              j8 = nint( rj(nn) )
              obsda%val(nn) = yobs_H08(i8,j8,3)
              obsda%val2(nn) = yobs_H08_clr(i8,j8,3)
              obsda%qc(nn) = iqc_good ! test

              if ( mod( nint(obsda%ri(nn)), H08_THIN_NG ) /= 0 .or. &
                   mod( nint(obsda%rj(nn)), H08_THIN_NG ) /= 0 ) then
                obsda%qc(nn) = iqc_obs_thin
              endif
            enddo

! allg here

#endif
          ENDIF ! H08

          write (6,'(3A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot = ', nslot
          write (6,'(3A,I6,A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot and processed by rank ' &
                                    , myrank, ' = ', nobs_slot


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,I3,A,I3,A,I3,A,F15.7)') '###### obsope_cal:obsope_step_2:        ',it,':',islot,':',iof,':',rrtimer-rrtimer00
  rrtimer00=rrtimer



! ###  -- TC vital assimilation -- ###
          if (obs_idx_TCX > 0 .and. obs_idx_TCY > 0 .and. obs_idx_TCP > 0) then
          if (obs(iof)%dif(obs_idx_TCX) == obs(iof)%dif(obs_idx_TCY) .and. &
              obs(iof)%dif(obs_idx_TCY) == obs(iof)%dif(obs_idx_TCP)) then
           
            if (obs(iof)%dif(obs_idx_TCX) > slot_lb .and. &
              obs(iof)%dif(obs_idx_TCX) <= slot_ub) then
              nslot = nslot + 3 ! TC vital obs should have 3 data (i.e., lon, lat, and MSLP)

              !!! bTC(1,:) : lon, bTC(2,:): lat, bTC(3,:): mslp
              ! bTC(1,:) : tcx (m), bTC(2,:): tcy (m), bTC(3,:): mslp
              allocate(bTC(3,0:MEM_NP-1))

              bTC = 9.99d33

              ! Note: obs(iof)%dat(obs_idx_TCX) is not longitude (deg) but X (m).
              !       Units of the original TC vital position are converted in
              !       subroutine read_obs in common_obs_scale.f90.
              !
              call phys2ij(obs(iof)%lon(obs_idx_TCX),obs(iof)%lat(obs_idx_TCX),rig,rjg) 
              call rij_g2l_auto(proc,rig,rjg,ritmp,rjtmp)  
              call search_tc_subdom(rig,rjg,v2dg,bTC(1,myrank_d),bTC(2,myrank_d),bTC(3,myrank_d))
  
              CALL MPI_BARRIER(MPI_COMM_d,ierr)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,bTC,3*MEM_NP,MPI_r_size,MPI_MIN,MPI_COMM_d,ierr)

              ! Assume MSLP of background TC is lower than 1100 (hPa). 
              bTC_mslp = 1100.0d2
              do n = 0, MEM_NP - 1
                write(6,'(3e20.5)')bTC(1,n),bTC(2,n),bTC(3,n) ! debug
                if (bTC(3,n) < bTC_mslp ) then
                  bTC_mslp = bTC(3,n)
                  bTC_proc = n
                endif
              enddo ! [ n = 0, MEM_NP - 1]

              if (myrank_d == proc) then
                do n = 1, 3
                  nobs = nobs + 1
                  nobs_slot = nobs_slot + 1
                  obsda%set(nobs) = iof
                  if(n==1) obsda%idx(nobs) = obs_idx_TCX
                  if(n==2) obsda%idx(nobs) = obs_idx_TCY
                  if(n==3) obsda%idx(nobs) = obs_idx_TCP
                  obsda%ri(nobs) = rig
                  obsda%rj(nobs) = rjg
                  ri(nobs) = ritmp
                  rj(nobs) = rjtmp

                  obsda%val(nobs) = bTC(n,bTC_proc)
                  obsda%qc(nobs) = iqc_good
                enddo ! [ n = 1, 3 ]

              endif
              deallocate(bTC)

            endif ! [ obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub ]
          endif ! [ obs_idx_TCX > 0 ...]
          endif !

        end do ! [ do iof = 1, OBS_IN_NUM ]

 
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

!write(6,*) '%%%%%%', MPI_WTIME(), nobs

      !!!!!!
      end if ! [ nobs_alldomain > 0 ]
      !!!!!!

      if (it == 1) then
        obsda%nobs = nobs
      else if (nobs /= obsda%nobs) then
        write (6, '(A)') '[Error] numbers of observations found are different among members.'
        stop
      end if

      write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' finishes processing member ', &
            im, ', subdomain id #', proc2mem(2,it,myrank+1)
      write (6,'(A,I8,A)') ' -- ', nobs, ' observations found'


      if (OBSDA_OUT) then
        write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is writing observations for member ', &
              im, ', subdomain id #', proc2mem(2,it,myrank+1)
        call file_member_replace(im, OBSDA_OUT_BASENAME, obsdafile)
        write (obsda_suffix(2:7),'(I6.6)') proc2mem(2,it,myrank+1)
        call write_obs_da(trim(obsdafile)//obsda_suffix,obsda,0)
      end if

      if (present(obsda_return)) then
        ! variables without an ensemble dimension
        if (it == 1) then
          obsda_return%nobs = nobs + obsda_return%nobs ! obsda_return%nobs: additional space for externally processed observations
          call obs_da_value_allocate(obsda_return,MEMBER)
        end if
        if (nobs > 0) then
          if (it == 1) then
            obsda_return%set(1:nobs) = obsda%set(1:nobs)
            obsda_return%idx(1:nobs) = obsda%idx(1:nobs)
            obsda_return%ri(1:nobs) = obsda%ri(1:nobs)
            obsda_return%rj(1:nobs) = obsda%rj(1:nobs)
            obsda_return%qc(1:nobs) = obsda%qc(1:nobs)
#ifdef H08
            obsda_return%lev(1:nobs) = obsda%lev(1:nobs)
            obsda_return%val2(1:nobs) = obsda%val2(1:nobs)
#endif
          else
            obsda_return%qc(1:nobs) = max(obsda_return%qc(1:nobs), obsda%qc(1:nobs))
#ifdef H08
            obsda_return%lev(1:nobs) = obsda_return%lev(1:nobs) + obsda%lev(1:nobs)
            obsda_return%val2(1:nobs) = obsda_return%val2(1:nobs) + obsda%val2(1:nobs)
#endif
          end if

          ! variables with an ensemble dimension
          obsda_return%ensval(im,1:nobs) = obsda%val(1:nobs)
        end if ! [ nobs > 0 ]
      end if ! [ present(obsda_return) ]


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,I3,A,8x,F15.7)') '###### obsope_cal:write_obs_da:         ',it,':',rrtimer-rrtimer00
  rrtimer00=rrtimer


    end if ! [ im >= 1 .and. im <= MEMBER ]

  end do ! [ it = 1, nitmax ]

  call obs_da_value_deallocate(obsda)

  deallocate ( ri, rj, v3dg, v2dg )

!write(6,*) ri(1),rj(1)
!write(6,*) '$$$$ 0'
!!  deallocate ( ri )
!write(6,*) '$$$$ 1'
!!  deallocate ( rj )
!write(6,*) '$$$$ 2'
!!  deallocate ( v3dg )
!write(6,*) '$$$$ 3'
!!  deallocate ( v2dg )
!write(6,*) '$$$$ 4'

end subroutine obsope_cal

!-----------------------------------------------------------------------
! Observation generator calculation
!-----------------------------------------------------------------------
SUBROUTINE obsmake_cal(obs)
  IMPLICIT NONE

  TYPE(obs_info),INTENT(INOUT) :: obs(OBS_IN_NUM)
  REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)

  integer :: islot,proc
  integer :: n,nslot,nobs,nobs_slot,ierr,iqc,iof
  integer :: nobsmax,nobsall
  real(r_size) :: rig,rjg,ri,rj,rk
  real(r_size) :: slot_lb,slot_ub
  real(r_size),allocatable :: bufr(:)
  real(r_size),allocatable :: error(:)

  CHARACTER(10) :: obsoutfile = 'obsout.dat'
  INTEGER :: ns 

!-----------------------------------------------------------------------

  write (6,'(A,I6.6,A,I6.6)') 'MYRANK ',myrank,' is processing subdomain id #', proc2mem(2,1,myrank+1)

  allocate ( v3dg (nlevh,nlonh,nlath,nv3dd) )
  allocate ( v2dg (nlonh,nlath,nv2dd) )

  do iof = 1, OBS_IN_NUM
    obs(iof)%dat = 0.0d0
  end do

  nobs = 0
  do islot = SLOT_START, SLOT_END
    slot_lb = (real(islot-SLOT_BASE,r_size) - 0.5d0) * SLOT_TINTERVAL
    slot_ub = (real(islot-SLOT_BASE,r_size) + 0.5d0) * SLOT_TINTERVAL
    write (6,'(A,I3,A,F9.1,A,F9.1,A)') 'Slot #', islot-SLOT_START+1, ': time interval (', slot_lb, ',', slot_ub, '] sec'

    call read_ens_history_iter(HISTORY_IN_BASENAME,1,islot,v3dg,v2dg)

    do iof = 1, OBS_IN_NUM
      IF(OBS_IN_FORMAT(iof) /= 3)THEN ! except H08 obs
        nslot = 0
        nobs_slot = 0
        do n = 1, obs(iof)%nobs

          if (obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub) then
            nslot = nslot + 1

            call phys2ij(obs(iof)%lon(n),obs(iof)%lat(n),rig,rjg)
            call rij_g2l_auto(proc,rig,rjg,ri,rj)

  !          if (myrank_d == 0) then
  !            print *, proc, rig, rjg, ri, rj
  !          end if

            if (proc < 0 .and. myrank_d == 0) then ! if outside of the global domain, processed by myrank_d == 0
              obs(iof)%dat(n) = undef
            end if

            if (myrank_d == proc) then
              nobs = nobs + 1
              nobs_slot = nobs_slot + 1

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

!              if (obs(iof)%elm(n) == id_radar_ref_obs .or. obs(iof)%elm(n) == id_radar_ref_zero_obs .or. obs(iof)%elm(n) == id_radar_vr_obs) then
              call phys2ijkz(v3dg(:,:,:,iv3dd_hgt),ri,rj,obs(iof)%lev(n),rk,iqc) ! z level only 10/24/2018 by TH
!              else
!                call phys2ijk(v3dg(:,:,:,iv3dd_p),obs(iof)%elm(n),ri,rj,obs(iof)%lev(n),rk,iqc)
!              end if

              if (iqc /= iqc_good) then
                obs(iof)%dat(n) = undef
              else
                select case (OBS_IN_FORMAT(iof))
                case (1)
                  call Trans_XtoY(obs(iof)%elm(n),ri,rj,rk, &
                                  obs(iof)%lon(n),obs(iof)%lat(n),v3dg,v2dg,obs(iof)%dat(n),iqc)
                case (2, 5)
                  call Trans_XtoY_radar(obs(iof)%elm(n),obs(iof)%meta(1),obs(iof)%meta(2),obs(iof)%meta(3),ri,rj,rk, &
                                        obs(iof)%lon(n),obs(iof)%lat(n),obs(iof)%lev(n),v3dg,v2dg,obs(iof)%dat(n),iqc)

                case (4, 6)
                  call Trans_XtoY_LT(obs(iof)%elm(n),ri,rj,rk, &
                                     v3dg,v2dg,obs(iof)%dat(n),iqc,obs(iof)%lev(n),myrank_d)
                end select

 !!! For radar observation, when reflectivity value is too low, do not generate ref/vr observations
 !!! No consideration of the terrain blocking effects.....

                if (iqc /= iqc_good) then
                  obs(iof)%dat(n) = undef
                end if
              end if

            end if ! [ myrank_d == proc ]

          end if ! [ obs%dif(n) > slot_lb .and. obs%dif(n) <= slot_ub ]

        end do ! [ n = 1, obs%nobs ]

      ENDIF

    end do ! [ iof = 1, OBS_IN_NUM ]

    write (6,'(3A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot = ', nslot
    write (6,'(3A,I6,A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot and processed by rank ', myrank, ' = ', nobs_slot

  end do ! [ islot = SLOT_START, SLOT_END ]

  deallocate ( v3dg, v2dg )

  if (myrank_d == 0) then
    nobsmax = 0
    nobsall = 0
    do iof = 1, OBS_IN_NUM
      if (obs(iof)%nobs > nobsmax) nobsmax = obs(iof)%nobs
      nobsall = nobsall + obs(iof)%nobs
    end do

    allocate ( bufr(nobsmax) )
    allocate ( error(nobsall) )

    call com_randn(nobsall, error) ! generate all random numbers at the same time
    ns = 0
  end if

  do iof = 1, OBS_IN_NUM

    call MPI_REDUCE(obs(iof)%dat,bufr(1:obs(iof)%nobs),obs(iof)%nobs,MPI_r_size,MPI_SUM,0,MPI_COMM_d,ierr)

    if (myrank_d == 0) then
      obs(iof)%dat = bufr(1:obs(iof)%nobs)

      do n = 1, obs(iof)%nobs
        select case(obs(iof)%elm(n))
        case(id_u_obs)
          obs(iof)%err(n) = OBSERR_U
        case(id_v_obs)
          obs(iof)%err(n) = OBSERR_V
        case(id_t_obs,id_tv_obs)
          obs(iof)%err(n) = OBSERR_T
        case(id_q_obs)
          obs(iof)%err(n) = OBSERR_Q
        case(id_rh_obs)
          obs(iof)%err(n) = OBSERR_RH
        case(id_ps_obs)
          obs(iof)%err(n) = OBSERR_PS
        case(id_radar_ref_obs,id_radar_ref_zero_obs)
          obs(iof)%err(n) = OBSERR_RADAR_REF
        case(id_radar_vr_obs)
          obs(iof)%err(n) = OBSERR_RADAR_VR
!
! -- Not available (02/09/2015)
!        case(id_H08IR_obs) ! H08
!          obs(iof)%err(n) = OBSERR_H08(ch) !H08
!        case default
          write(6,'(A)') 'warning: skip assigning observation error (unsupported observation type)' 
        end select

        if (obs(iof)%dat(n) /= undef .and. obs(iof)%err(n) /= undef) then
          obs(iof)%dat(n) = obs(iof)%dat(n) + obs(iof)%err(n) * error(ns+n)
        end if

!print *, '######', obs%elm(n), obs%dat(n)
      end do ! [ n = 1, obs(iof)%nobs ]

      ns = ns + obs(iof)%nobs
    end if ! [ myrank_d == 0 ]

  end do ! [ iof = 1, OBS_IN_NUM ]

  if (myrank_d == 0) then
    deallocate ( bufr )
    deallocate ( error )

    call write_obs_all(obs, missing=.false., file_suffix='.out') ! only at the head node
  end if

end subroutine obsmake_cal

!-------------------------------------------------------------------------------
! Model-to-observation simulator calculation
!-------------------------------------------------------------------------------
subroutine obssim_cal(v3dgh, v2dgh, v3dgsim, v2dgsim, stggrd)
  use scale_grid, only: &
      GRID_CX, GRID_CY, &
      GRID_CZ, &
      DX, DY
  use scale_grid_index, only: &
      IHALO, JHALO, KHALO
  use scale_mapproj, only: &
      MPRJ_xy2lonlat

  implicit none

  real(r_size), intent(in) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  real(r_size), intent(in) :: v2dgh(nlonh,nlath,nv2dd)
  real(r_size), intent(out) :: v3dgsim(nlev,nlon,nlat,OBSSIM_NUM_3D_VARS)
  real(r_size), intent(out) :: v2dgsim(nlon,nlat,OBSSIM_NUM_2D_VARS)
  integer, intent(in), optional :: stggrd

  integer :: i, j, k, iv3dsim, iv2dsim
  real(r_size) :: ri, rj, rk
  real(r_size) :: lon, lat, lev
  real(r_size) :: tmpobs
  integer :: tmpqc
  real(r_size) :: tmplev

  integer :: ig, jg

#ifdef H08
! -- for Himawari-8 obs --
  real(r_size) :: yobs_H08(nlon,nlat,NIRB_HIM8),yobs_H08_clr(nlon,nlat,NIRB_HIM8)
  real(r_size) :: plev_obs_H08(nlon,nlat,NIRB_HIM8)
  integer :: qc_H08(nlon,nlat,NIRB_HIM8)
  integer :: ch

  logical :: USE_HIM8 = .false. ! initialize

#endif

!-------------------------------------------------------------------------------

  write (6,'(A,I6.6,A,I6.6)') 'MYRANK ', myrank, ' is processing subdomain id #', proc2mem(2,1,myrank+1)

  do j = 1, nlat
    rj = real(j + JHALO, r_size)

    do i = 1, nlon
      ri = real(i  + IHALO, r_size)

      call ij_l2g(myrank_d, i, j, ig, jg)

      if (RADAR_IDEAL) then
        ! Idealized experiment without map projections
        lon = (ri-1.0d0) * DX + GRID_CX(1)
        lat = (rj-1.0d0) * DY + GRID_CY(1)
      else
        call MPRJ_xy2lonlat((ri-1.0d0) * DX + GRID_CX(1), (rj-1.0d0) * DY + GRID_CY(1), lon, lat)
        lon = ri ! lon * rad2deg
        lat = rj ! lat * rad2deg
      endif


      do k = 1, nlev
        rk = real(k + KHALO, r_size)

        do iv3dsim = 1, OBSSIM_NUM_3D_VARS
          select case (OBSSIM_3D_VARS_LIST(iv3dsim))
          case (id_radar_ref_obs, id_radar_ref_zero_obs, id_radar_vr_obs, id_radar_prh_obs)
!            lev = v3dgh(k+KHALO, i+IHALO, j+JHALO, iv3dd_hgt)
! This is an idealized experiment without map projections
            lev = GRID_CZ(k+KHALO)
            call Trans_XtoY_radar(OBSSIM_3D_VARS_LIST(iv3dsim), OBSSIM_RADAR_LON, OBSSIM_RADAR_LAT, OBSSIM_RADAR_Z, ri, rj, rk, &
                                  lon, lat, lev, v3dgh, v2dgh, tmpobs, tmpqc, stggrd)
            if (tmpqc == iqc_ref_low) tmpqc = iqc_good ! when process the observation operator, we don't care if reflectivity is too small

          case (id_lt3d_obs, id_fp3d_obs)
            call Trans_XtoY_LT(OBSSIM_3D_VARS_LIST(iv3dsim), ri, rj, rk, &
                               v3dgh, v2dgh, tmpobs, tmpqc, tmplev, myrank_d)

          case (-999)
!            if (iv3dsim == 7) then
!              tmpobs = v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_qc) &
!                     + v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_qr) &
!                     + v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_qi) &
!                     + v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_qs) &
!                     + v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_qg)
!            elseif (iv3dsim == 8) then
!              tmpobs = v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_cc) &
!                     + v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_cr) &
!                     + v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_ci) &
!                     + v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_cs) &
!                     + v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_cg)
!              if (OBSSIM_IN_TYPE == 'restart') then
!                tmpobs = tmpobs * 1.e-6 ! fC(10^-15) => nC(10^-9)
!              endif
!            else
!            tmpobs = v3dgh(k+KHALO,i+IHALO,j+JHALO,abs(OBSSIM_3D_VARS_LIST(iv3dsim)))
!            tmpobs = v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dsim)
!            endif
            tmpqc = iqc_good 
          case default
            if (OBSSIM_3D_VARS_LIST(iv3dsim) < 0) then
              tmpobs = v3dgh(k+KHALO,i+IHALO,j+JHALO,abs(OBSSIM_3D_VARS_LIST(iv3dsim)))
              tmpqc = iqc_good 

            else
              call Trans_XtoY(OBSSIM_3D_VARS_LIST(iv3dsim), ri, rj, rk, &
                              lon, lat, v3dgh, v2dgh, tmpobs, tmpqc, stggrd)
            endif
          end select

          if (tmpqc == 0) then
            v3dgsim(k,i,j,iv3dsim) = real(tmpobs, r_sngl)
          else
            v3dgsim(k,i,j,iv3dsim) = real(undef, r_sngl)
          end if
        end do ! [ iv3dsim = 1, OBSSIM_NUM_3D_VARS ]

        ! 2D observations calculated when k = 1
        if (k == 1) then
          do iv2dsim = 1, OBSSIM_NUM_2D_VARS
            select case (OBSSIM_2D_VARS_LIST(iv2dsim))
!            case (id_H08IR_obs)               !!!!!! H08 as 2D observations ???
!              call Trans_XtoY_radar_H08(...)
!            case (id_tclon_obs, id_tclat_obs, id_tcmip_obs)
!              call ...

            case (id_lt2d_obs)
              call Trans_XtoY_LT(OBSSIM_2D_VARS_LIST(iv2dsim), ri, rj, rk, &
                                 v3dgh, v2dgh, tmpobs, tmpqc, tmplev, myrank_d)

            case (id_fp2d_obs)
              call Trans_XtoY_LT(OBSSIM_2D_VARS_LIST(iv2dsim), ri, rj, rk, &
                                 v3dgh, v2dgh, tmpobs, tmpqc, tmplev, myrank_d)
            case (id_H08IR_obs)
              USE_HIM8 = .true.
              cycle
            case default
              call Trans_XtoY(OBSSIM_2D_VARS_LIST(iv2dsim), ri, rj, rk, &
                              lon, lat, v3dgh, v2dgh, tmpobs, tmpqc, stggrd)
            end select

            if (tmpqc == 0) then
              v2dgsim(i,j,iv2dsim) = real(tmpobs, r_sngl)
            else
              v2dgsim(i,j,iv2dsim) = real(undef, r_sngl)
            end if
          end do ! [ iv2dsim = 1, OBSSIM_NUM_2D_VARS ]
        end if ! [ k == 1 ]

      end do ! [ k = 1, nlev ]

    end do ! [ i = 1, nlon ]

  end do ! [ j = 1, nlat ]

#ifdef H08
  if (USE_HIM8) then
    call Trans_XtoY_H08_allg(v3dgh,v2dgh,yobs_H08,yobs_H08_clr,&
                             qc_H08)
    ch = 0
    do iv2dsim = 1, OBSSIM_NUM_2D_VARS
      if(OBSSIM_2D_VARS_LIST(iv2dsim) /= id_H08IR_obs .or. ch > NIRB_HIM8)cycle

      ch = ch + 1
      v2dgsim(:,:,iv2dsim) = real(yobs_H08(:,:,ch), r_sngl)
    enddo ! [iv2dsim = 1, OBSSIM_NUM_2D_VARS]
  endif
#endif

!-------------------------------------------------------------------------------
  return
end subroutine obssim_cal

subroutine get_fpmax_loc2d( v3dgh, fp_max, fp_max_lon, fp_max_lat ) 
  use scale_grid, only: &
      GRID_CXG, GRID_CYG
  use scale_grid_index, only: &
      IHALO, JHALO, &
      KS, KE
  implicit none

  integer :: proc_i, proc_j
  integer :: ishift, jshift
  real(r_size), intent(in) :: v3dgh(nlevh,nlonh,nlath,nv3dd)

  real(r_size) :: fp2d(nlong,nlatg)

  integer :: i, j 
  integer :: ierr

  real(r_size), intent(out) :: fp_max 
  real(r_size), intent(out) :: fp_max_lon, fp_max_lat
  integer :: fp_maxi, fp_maxj
  real(r_size) :: fp_tmp


  call rank_1d_2d(myrank_d, proc_i, proc_j)
  ishift = proc_i * nlon
  jshift = proc_j * nlat

  fp2d(:,:) = 0.0d0
  do j = 1, nlat
  do i = 1, nlon
    fp2d(i+ishift, j+jshift) = real( sum(v3dgh(KS:KE,IHALO+i,JHALO+j,iv3dd_fp) ), r_size)
  enddo
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE, fp2d, nlong*nlatg, MPI_r_size, MPI_SUM, MPI_COMM_d, ierr)

  fp_max = -1.0d0
  fp_maxi = int( nlong / 2 )
  fp_maxj = int( nlatg / 2 )
  do j = 3, nlatg - 1
  do i = 3, nlong - 1
    if ( mod( i, XY_THINNING_LT ) /= 0 .or. &
         mod( j, XY_THINNING_LT ) /= 0 ) cycle

    ! should be consistent with Trans_XtoY_LT ( reduce into 8-km x 8-km grid )
    fp_tmp = sum( fp2d(i-2:i+1,j-2:j+1) )
    if ( fp_tmp > fp_max ) then
      fp_max = fp_tmp
      fp_maxi = i
      fp_maxj = j
    endif
  enddo
  enddo

  fp_max_lon = GRID_CXG( fp_maxi + IHALO )
  fp_max_lat = GRID_CYG( fp_maxj + JHALO )

  fp_max = max( fp_max, 0.0d0 )

  return
end subroutine get_fpmax_loc2d
!=======================================================================

END MODULE obsope_tools
