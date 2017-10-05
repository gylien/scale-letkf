module scale_H08_fwd12
implicit none

contains

subroutine SCALE_RTTOV12_fwd(nchannels,&
                             nlevels,&
                             nprof,&
                             tmp_p,&
                             tmp_t,&
                             tmp_qv,&
                             tmp_qc,&
                             tmp_qice,&
                             tmp_t2m,&
                             tmp_q2m,&
                             tmp_p2m,&
                             tmp_u2m,&
                             tmp_v2m,&
!                     & tmp_soze, tmp_soaz, tmp_saze, tmp_saaz,       &
                             tmp_elev,&
                             tmp_lon,&
                             tmp_lat,&
                             tmp_land,&
                             tmp_zenith,&
                             kidx_rlx, &
                             kadd_org, &
                             RD_presh, &
                             RD_temph, &
                             btall_out,&
                             btclr_out,&
                             mwgt_plev,&
                             ctop_out)
  !
  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2016, EUMETSAT, All Rights Reserved.
  !
  !     *************************************************************

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & platform_name,       &
       & inst_name,           &
       & qmin,                &
       & tmin

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit
  USE common, ONLY : r_size
  USE scale_const, ONLY: &
        Rdry    => CONST_Rdry, &
        Rvap    => CONST_Rvap, &
        Deg2Rad => CONST_D2R,  &
        CONST_GRAV
  USE common_nml, ONLY: &
        H08_RTTOV_CFRAC_CNST, &
        H08_RTTOV_MINQ_CTOP,  &
        H08_RTTOV_COEF_PATH,  &
        H08_RTTOV_PROF_SHIFT, &
        H08_RTTOV_KADD
  use scale_grid, only: &
      GRID_FZ
  use scale_grid_index, only: &
      KHALO
  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  INTEGER, INTENT(IN) :: nprof
  INTEGER, INTENT(IN) :: nlevels

  Real(r_size),INTENT(IN) :: tmp_p(nlevels, nprof)
  Real(r_size),INTENT(IN) :: tmp_t(nlevels, nprof)
  Real(r_size),INTENT(IN) :: tmp_qv(nlevels, nprof)
  Real(r_size),INTENT(IN) :: tmp_qc(nlevels, nprof)
  Real(r_size),INTENT(IN) :: tmp_qice(nlevels, nprof)
  Real(r_size),INTENT(IN) :: tmp_t2m(nprof)
  Real(r_size),INTENT(IN) :: tmp_q2m(nprof)
  Real(r_size),INTENT(IN) :: tmp_p2m(nprof)
  Real(r_size),INTENT(IN) :: tmp_u2m(nprof)
  Real(r_size),INTENT(IN) :: tmp_v2m(nprof)
  Real(r_size),INTENT(IN) :: tmp_elev(nprof)
  Real(r_size),INTENT(IN) :: tmp_lon(nprof)
  Real(r_size),INTENT(IN) :: tmp_lat(nprof)
  Real(r_size),INTENT(IN) :: tmp_land(nprof)
  Real(r_size),INTENT(IN) :: tmp_zenith(nprof)

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output
  INTEGER(KIND=jpim), PARAMETER :: mxchn = 9000 ! max number of channels

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls


  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=11)  :: NameOfRoutine = 'example_fwd'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename='./rtcoef_himawari_8_ahi.bin'
  CHARACTER(LEN=256) :: cld_coef_filename='./sccldcoef_himawari_8_ahi.bin'
  INTEGER(KIND=jpim) :: nthreads 
  INTEGER(KIND=jpim) :: dosolar = 0
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  ! loop variables
  INTEGER(KIND=jpim) :: j
  INTEGER(KIND=jpim) :: ich, jch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff

! by T.Honda
  REAL(Kind=r_size),INTENT(OUT) :: btall_out(nchannels,nprof)
  REAL(Kind=r_size),INTENT(OUT) :: btclr_out(nchannels,nprof)
  REAL(Kind=r_size),INTENT(OUT) :: ctop_out(nprof)
  REAL(Kind=r_size) :: ptmp

  REAL(Kind=r_size) :: rdp, max_wgt, tmp_wgt
  REAL(Kind=r_size),INTENT(OUT) :: mwgt_plev(nchannels,nprof) ! Max weight level (Pa)

  logical :: debug = .false.
  !logical :: debug = .true.

  real(kind=jprb) :: Rd 
  real(kind=jprb) :: Rv 
  real(kind=jprb) :: grav

  real(kind=jprb) :: epsb 
  real(kind=jprb) :: repsb 

  integer :: rdk, orgk

  integer(kind=jpim),intent(in) :: kidx_rlx, kadd_org
  real(kind=r_size),intent(in)  :: RD_presh(nlevs+kadd_org+1)
  real(kind=r_size),intent(in)  :: RD_temph(nlevs+kadd_org+1)
  integer(kind=jpim) :: kadd
  real(kind=jprb) :: rat, rdz, dz, tmp_dif

  if(debug) write(6,'(1x,a)')"hello from RTTOV"

  kadd = min(kadd_org,H08_RTTOV_KADD)

! -- set thermodynamic constants
  Rd = real(Rdry,kind=jprb)
  Rv = real(Rvap,kind=jprb)
  epsb = Rd / Rv 
  repsb = 1.0_jprb / epsb

  grav = real(CONST_GRAV,kind=jprb)

  !
  ! The usual steps to take when running RTTOV are as follows:
  !   1. Specify required RTTOV options
  !   2. Read coefficients
  !   3. Set up the chanprof array with the channels/profiles to simulate
  !   4. Allocate RTTOV profile, radiance and transmittance structures
  !   5. Read input profile(s)
  !   6. Set up surface emissivity and/or reflectance
  !   7. Call rttov_direct and store results
  !   8. Deallocate all structures and arrays

  errorstatus     = 0_jpim
  alloc_status(:) = 0_jpim

  dosolar=0

  do ich = 1, nchannels
    input_chan(ich)=ich
  end do
  input_ems(:)=0.0
  input_brdf(:)=0.0

  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  IF (dosolar == 1) THEN
    opts % rt_ir % addsolar = .TRUE.           ! Include solar radiation
  ELSE
    opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
  ENDIF
  opts % interpolation % addinterp  = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method
  opts % rt_all % addrefrac         = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addaerosl          = .FALSE. ! Don't include aerosol effects

  opts % rt_ir % addclouds          = .TRUE. ! Include cloud effects

  opts % rt_ir % ir_scatt_model      = 2       ! Scattering model for emission source term:
                                               !   1 => DOM; 2 => Chou-scaling
  opts % rt_ir % vis_scatt_model     = 1       ! Scattering model for solar source term:
                                               !   1 => DOM; 2 => single-scattering
  opts % rt_ir % dom_nstreams        = 8       ! Number of streams for Discrete Ordinates (DOM)

  opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
  opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
  opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
  opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
  opts % rt_ir % co_data             = .FALSE. !
  opts % rt_ir % so2_data            = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !

  opts % config % verbose            = .TRUE.  ! Enable printing of warnings

! added by T.Honda(2015/06/10)
  opts % interpolation % reg_limit_extrap  = .TRUE.  
  opts % config % apply_reg_limits  = .TRUE.


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------

  ! Read optical depth and cloud coefficient files together
  CALL rttov_read_coefs(errorstatus, coefs, opts, &
                        file_coef=trim(H08_RTTOV_COEF_PATH)//trim(coef_filename), &
                        file_sccld=trim(H08_RTTOV_COEF_PATH)//trim(sccoef_filename), &
                        form_coef='unformatted', form_sccld='unformatted')

  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF

  if(debug) write(6,'(1x,a)')"hello from RTTOV3"

  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nch * nprof

  write(6,'(a,3i9)')"DEBUG",nchanprof, nch, nprof

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        nprof,                   &
        nchanprof,               &
        nlevels+kadd,            &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance, &
        init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  ich = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nch
      ich = ich + 1_jpim
      chanprof(ich)%prof = j
      chanprof(ich)%chan = jch
!      chanprof(ich)%chan = channel_list(jch)
    ENDDO
  ENDDO


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== Read profiles == start =============

  if(debug) then
    WRITE(6,*) 'debug information'
    WRITE(6,*) 'max p, max t',maxval(tmp_p),maxval(tmp_t)
    WRITE(6,*) 'min p, min t',minval(tmp_p),minval(tmp_t)
    WRITE(6,*) 'max p2, max t2',maxval(tmp_p2m),maxval(tmp_t2m)
    WRITE(6,*) 'min p2, min t2',minval(tmp_p2m),minval(tmp_t2m)
    WRITE(6,*) 'max qv, min qv',maxval(tmp_qv),minval(tmp_qv)
    WRITE(6,*) '-- t prof --'
    WRITE(6,*) 'size t:',size(tmp_t(:,1)),size(tmp_p(:,1))
    do j = 1, nlevels
      WRITE(6,*) 'qv:',tmp_qv(j,2)
    enddo
    WRITE(6,*) '-- end t prof --'
    do j = 1, nlevels
      WRITE(6,*) 't:',tmp_t(j,2)
    enddo
  endif

  if(debug) WRITE(6,*) 'START SUBSTITUTE PROFILE'

  ! Note: Profiles are from top to surface.
  rdz = 1.0d3 / (GRID_FZ(nlevs+KHALO) - GRID_FZ(nlevs+KHALO-kidx_rlx))

  ! Loop over all profiles
  do iprof = 1, nprof

    if(H08_RTTOV_PROF_SHIFT)then
      tmp_dif = RD_temph(kadd+1) - tmp_t(1,iprof)
    else
      tmp_dif = 0.0d0
    endif ! H08_RTTOV_PROF_SHIFT

    do ilev = 1, kadd
      !rdk = (kadd_org + 1 + ilev) - kadd
      rdk = (kadd_org + ilev) - kadd
      if(debug .and. iprof == 1) print *,"ilev,rdk",ilev,rdk
      profiles(iprof)%p(ilev) = RD_presh(rdk) ! (hPa)
      profiles(iprof)%t(ilev) = max(RD_temph(rdk) - tmp_dif,tmin + tmin * 0.01_jprb) ! (K)

      profiles(iprof)%q(ilev) = qmin + qmin * 0.01_jprb

    enddo

    ! model grid values (> about [H08_RTTOV_RLX_HGT] m) are relaxed to climatology 
    do ilev = kadd + 1, kadd + kidx_rlx
      orgk = ilev - kadd ! k index for original profile
      rdk = (kadd_org + orgk) ! k index for rd (climatological) profile

      ! distance from the model top [GRID_FZ(nlevs+KHALO)]
      dz = abs(GRID_FZ(nlevs+KHALO) - GRID_FZ(nlevs+KHALO-orgk+1)) * 1.d-3 ! (km)

      if(kadd == 0 .or. H08_RTTOV_PROF_SHIFT)then
        rat = 1.0d0
      else
        rat = dz * rdz
      endif

      !if(debug .and. iprof == 1) print *,"CLIM:",tmp_p(orgk,iprof)*1.d-2,&
      !                                           tmp_t(orgk,iprof),&
      !                                           RD_presh(rdk),&
      !                                           RD_temph(rdk)

      profiles(iprof)%p(ilev) = real(tmp_p(orgk,iprof),kind=jprb) * 0.01_jprb ! (hPa)
      profiles(iprof)%t(ilev) = max(rat * real(tmp_t(orgk,iprof),kind=jprb) &
                                    + (1.0d0 - rat) * RD_temph(rdk), &
                                    tmin + tmin * 0.01_jprb) ! (K)

      profiles(iprof)%q(ilev) = qmin + qmin * 0.01_jprb

    enddo

    profiles(iprof)%s2m%t = real(tmp_t2m(iprof),kind=jprb)
    profiles(iprof)%s2m%q = real(tmp_q2m(iprof),kind=jprb) * q_mixratio_to_ppmv ! (ppmv)

    if(profiles(iprof)%s2m%t < tmin) profiles(iprof)%s2m%t = tmin + tmin * 0.01_jprb
    if(profiles(iprof)%s2m%q < qmin) profiles(iprof)%s2m%q = qmin + qmin * 0.01_jprb


    profiles(iprof)%s2m%p = real(tmp_p2m(iprof),kind=jprb) * 0.01_jprb ! (hPa)
    profiles(iprof)%s2m%u = real(tmp_u2m(iprof),kind=jprb)
    profiles(iprof)%s2m%v = real(tmp_v2m(iprof),kind=jprb)
    profiles(iprof)%s2m%wfetc = 100000.0_jprb

    profiles(iprof) % skin % t = max(real(tmp_t2m(iprof),kind=jprb), tmin + tmin * 0.01_jprb)

!
!   fastem is used for MW
!    profiles(iprof) % skin % fastem(1) = 3.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(2) = 5.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(3) =15.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(4) = 0.1 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(5) = 0.3 ! comment out (11/18/2015)
!

    profiles(iprof) % skin % surftype = int(tmp_land(iprof))
    profiles(iprof) % skin % watertype = 1 ! tentative (11/18/2015)

    profiles(iprof) % elevation = real(tmp_elev(iprof),kind=jprb) * 0.001_jprb ! (km)
    profiles(iprof) % latitude  = real(tmp_lat(iprof),kind=jprb)
    profiles(iprof) % longitude = real(tmp_lon(iprof),kind=jprb)


    profiles(iprof)% zenangle = real(tmp_zenith(iprof),kind=jprb) ! (11/18/2015)
    if(mod(iprof,1000) == 0 .and. debug)write(6,'(a,f15.10)'),'zenangle ',profiles(iprof)% zenangle
    if(mod(iprof,1000) == 0 .and. debug)write(6,'(a,2f10.5)'),' ',tmp_lon(iprof),tmp_lat(iprof)

! These are parameters for simple cloud.
! Not used.
    profiles(iprof) % ctp       = 500.0_jprb
    profiles(iprof) % cfraction = 0.0_jprb

!  --- 6.1 microphysical clouds 
!    if(icldprf==1) then
    ! --convert kg/kg into g/m3, make liq.cloud and ice.cloud amount, and upside-down

      ctop_out(iprof) = -1.0d0

      do ilev = kadd + 1, kadd + nlevs - 1
        orgk = ilev - kadd ! k index for original profile

        ! ilev
        tv = real(tmp_t(orgk,iprof) * (1.0d0+tmp_qv(orgk,iprof)) * repsb &
                   / (1.0d0 + tmp_qv(orgk,iprof)), kind=jprb)

        kgkg2gm3 = real(tmp_p(orgk,iprof),kind=jprb) / (Rd * tv) * 1000.0_jprb

        liqc1 = real(max(tmp_qc(orgk,iprof),0.0_r_size),kind=jprb) * kgkg2gm3
        icec1 = real(max(tmp_qice(orgk,iprof),0.0_r_size),kind=jprb) * kgkg2gm3

        ! ilev + 1
        tv = real(tmp_t(orgk+1,iprof) * (1.0d0+tmp_qv(orgk+1,iprof)) * repsb &
                   / (1.0d0 + tmp_qv(orgk+1,iprof)), kind=jprb)

        kgkg2gm3 = real(tmp_p(orgk+1,iprof),kind=jprb) / (Rd * tv) * 1000.0_jprb

        liqc2 = real(max(tmp_qc(orgk+1,iprof),0.0_r_size),kind=jprb) * kgkg2gm3
        icec2 = real(max(tmp_qice(orgk+1,iprof),0.0_r_size),kind=jprb) * kgkg2gm3

        profiles(iprof) % cloud(2,ilev) = & !stratus maritime (default)
                   (liqc1 + liqc2) * 0.5_jprb
        profiles(iprof) % cloud(6,ilev) = &
                   (icec1 + icec2) * 0.5_jprb
!
! cloud fraction & cloud top diagnosis
        ptmp = (tmp_p(orgk+1,iprof) + tmp_p(orgk,iprof))*0.5_jprb

        profiles(iprof) % cfrac(ilev) = min((profiles(iprof) % cloud(2,ilev) + &
                                             profiles(iprof) % cloud(6,ilev)) / H08_RTTOV_CFRAC_CNST, &
                                             1.0_jprb)
        ! Need to modify? if openmp
        if(profiles(iprof) % cloud(2,ilev) + &
           profiles(iprof) % cloud(6,ilev) >= H08_RTTOV_MINQ_CTOP)then
          if(ctop_out(iprof) < 0.0d0)then
            ctop_out(iprof) = ptmp
          endif
        endif

      end do ! ilev
    endif ! addclouds

    if(debug .and. mod(iprof,200)==0)then
      do ilev = 1, nlevs + kadd - 1
        write(6,'(a,i5,4f11.4)')"DEBUG PROF",ilev,profiles(iprof) % t(ilev),&
                                                  profiles(iprof) % q(ilev),&
                                                  profiles(iprof) % p(ilev),&
                                                  profiles(iprof) % cloud(6,ilev)
      end do ! ilev
    endif
  enddo ! prof









    ! Cloud part
     ! Note: profiles(:)%mmr_cldaer = .true. (default), so that units for cloud variables are kg/kg 

     ! Select the ice cloud parameterisation to use, in this case the Baran 2014 scheme:

    profiles(iprof) % ice_scheme = 2 

    profiles(iprof) % cloud(:,:) = 0._jprb
    profiles(iprof) % cfrac(:)   = 0._jprb


    ctop_out(iprof) = -1.0d0 
    do ilev=1,nlevels-1
      ! Stratus maritime 
      profiles(iprof) % cloud(2,ilev) = max(real((tmp_qc(ilev,iprof) + tmp_qc(ilev+1,iprof)) * 0.5_jprb,&
                                            kind=jprb),&
                                            0.0_jprb)

      ! Ice cloud
      profiles(iprof) % cloud(6,ilev) = max(real((tmp_qice(ilev,iprof) + tmp_qice(ilev+1,iprof)) * 0.5_jprb,&
                                            kind=jprb),&
                                            0.0_jprb)

      ! Cloud fraction diagnosis
      ptmp = (tmp_p(ilev+1,iprof) + tmp_p(ilev,iprof))*0.5_jprb

      profiles(iprof) % cfrac(ilev) = min((profiles(iprof) % cloud(2,ilev) + &
                                           profiles(iprof) % cloud(6,ilev)) / jcfrac_cnst, 1.0_jprb)

      ! Cloud top level (hPa) diagnosis
      if(profiles(iprof) % cloud(2,ilev) + &
         profiles(iprof) % cloud(6,ilev) >= H08_RTTOV_MINQ_CTOP)then
        if(ctop_out(iprof) < 0.0d0)then
          ctop_out(iprof) = ptmp
        endif
      endif

    enddo ! ilev
  enddo ! iprof
  if(debug) WRITE(6,*) '[RTTOV12]: END SUBSTITUTE PROFILE'

  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------
  
  ! In this example we have no values for input emissivities
  emissivity(:) % emis_in = 0._jprb
  
  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)
  
  ! In this example we have no values for input reflectances
  reflectance(:) % refl_in = 0._jprb
  
  ! Calculate BRDF within RTTOV where the input BRDF value is zero or less
  ! (all channels in this case)
  calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)
  
  ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
  reflectance(:) % refl_cloud_top = 0._jprb

  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------

  if(debug) write(6,*)"[RTTOV]: Enter direct"
  if(debug) write(6,*)"[RTTOV]: We assume that this job is excuted with 8 threads as in the K computer."
  nthreads = 8

  IF (nthreads <= 1) THEN
    CALL rttov_direct(                &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance) ! inout input/output BRDFs per channel
  ELSE
    CALL rttov_parallel_direct(       &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance,&! inout input/output BRDFs per channel
            nthreads    = nthreads)    ! in    number of threads to use
  ENDIF

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ELSE
    if(debug) write(6,*)"Pass direct"
  ENDIF

  ! --- Output the results --------------------------------------------------

  do iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    !
    !     OUTPUT RESULTS
    !
    btall_out(1:nchannels,iprof) = real(radiance%bt(1+joff:nchannels+joff), kind = r_size)
    btclr_out(1:nchannels,iprof) = real(radiance%bt_clear(1+joff:nchannels+joff), kind=r_size)

    do ich = 1, nchannels
      rdp = 1.0d0 / (abs(profiles(iprof)%p(1) - profiles(iprof)%p(2)) * 1.0d2) ! Pa
      max_wgt = abs(transmission % tau_levels(1,joff+ich) &
                  - transmission % tau_levels(2,joff+ich)) * rdp
      mwgt_plev(ich,iprof) = (profiles(iprof)%p(1) + profiles(iprof)%p(2)) * 0.5d2 ! Pa

      if(debug .and. mod(iprof,100)==0)then
        write(6,'(a,i6,a,i3,a,f11.5,a,i4,a,f8.2,a,f7.2,a,f10.4)')"WGT,",&
              iprof,",",ich+6,",",max_wgt*1.e6,",",1,",",&
              (profiles(iprof)%p(1) + profiles(iprof)%p(2)) * 0.5d0,",",&
              (profiles(iprof)%t(1) + profiles(iprof)%t(2)) * 0.5d0,",",&
              (profiles(iprof)%q(1) + profiles(iprof)%q(2)) * 0.5d0 / q_mixratio_to_ppmv

      endif

      ! TOA to the ground
      do ilev = 2, nlevs - 1 + kadd
        rdp = 1.0d0 / (abs(profiles(iprof)%p(ilev) - profiles(iprof)%p(ilev+1)) * 1.0d2) ! Pa
        tmp_wgt = abs(transmission % tau_levels(ilev,joff+ich) &
                    - transmission % tau_levels(ilev+1,joff+ich)) * rdp

        if(tmp_wgt > max_wgt)then
          max_wgt = tmp_wgt
          mwgt_plev(ich,iprof) = (profiles(iprof)%p(ilev) + profiles(iprof)%p(ilev+1)) * 0.5d2 ! Pa
        endif

        if(debug .and. mod(iprof,100)==0)then
          write(6,'(a,i6,a,i3,a,f11.5,a,i4,a,f8.2,a,f7.2,a,f10.4)')"WGT,",&
                iprof,",",ich+6,",",tmp_wgt*1.e6,",",ilev,",",&
                (profiles(iprof)%p(ilev) + profiles(iprof)%p(ilev+1)) * 0.5d0,",",&
                (profiles(iprof)%t(ilev) + profiles(iprof)%t(ilev+1)) * 0.5d0,",",&
                (profiles(iprof)%q(ilev) + profiles(iprof)%q(ilev+1)) * 0.5d0 / q_mixratio_to_ppmv
        endif

      enddo ! ilev
    enddo ! ich

  enddo ! iprof

  if(debug) write(6,'(a)')"End WGT calculation"

  ! --- End of output section -----------------------------------------------

  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  CALL rttov_alloc_direct( &
        errorstatus,                         &
        0_jpim,                              &  ! 0 => deallocate
        nprof,                               &
        nchanprof,                           &
        nlevels,                             &
        chanprof,                            &
        opts,                                &
        profiles,                            &
        coefs,                               &
        transmission,                        &
        radiance,                            &
        calcemis=calcemis,                   &
        emissivity=emissivity,               &
        calcrefl=calcrefl,                   &
        reflectance=reflectance,             &
        aer_maxnmom=opts%rt_ir%dom_nstreams, &
        aer_nphangle=nphangle,               &
        aer_opt_param=aer_opt_param)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF

! Format definitions for output
111  FORMAT(1X,10I8)
222  FORMAT(1X,10F8.2)
444  FORMAT(1X,10F8.3)
777  FORMAT(/,A,A9,I3)

  if(debug) write(*,*)'Successfully finished!!'

return
end subroutine SCALE_RTTOV12_fwd

end module scale_H08_fwd12

