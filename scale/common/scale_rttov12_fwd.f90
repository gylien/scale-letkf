module scale_rttov12_fwd
implicit none

contains

subroutine cld_ir_fwd(nchannels,&
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
                           tmp_sat_zangle, & ! satellite zenith angle
                           btall_out,& 
                           btclr_out)!,& 
!                           trans_out)

  ! rttov_const contains useful RTTOV constants
  use rttov_const, only :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         qmin,                &
         qmax

  ! rttov_types contains definitions of all RTTOV data types
  use rttov_types, only :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

  ! jpim, jprb and jplm are the RTTOV integer, real and logical kinds
  use parkind1, only : jpim, jprb, jplm

  use rttov_unix_env, only : rttov_exit

  use common, only : r_size
  use scale_const, only: &
        Rdry    => CONST_Rdry, &
        Rvap    => CONST_Rvap, &
        Deg2Rad => CONST_D2R
  use common_nml, only: &
        H08_RTTOV_COEF_PATH, &
        H08_RTTOV_MinQ, &
        H08_RTTOV_CFRAC_CNST, &
        H08_RTTOV_CLD,  &
        H08_RTTOV_NTHREAD, &
        NVIS_HIM8
  implicit none

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  ! RTTOV variables/structures
  !====================
  type(rttov_options)              :: opts                     ! Options structure
  type(rttov_coefs)                :: coefs                    ! Coefficients structure
  type(rttov_chanprof),    pointer :: chanprof(:)    => null() ! Input channel/profile list
  logical(kind=jplm),      pointer :: calcemis(:)    => null() ! Flag to indicate calculation of emissivity within RTTOV
  type(rttov_emissivity),  pointer :: emissivity(:)  => null() ! Input/output surface emissivity
  logical(kind=jplm),      pointer :: calcrefl(:)    => null() ! Flag to indicate calculation of BRDF within RTTOV
  type(rttov_reflectance), pointer :: reflectance(:) => null() ! Input/output surface BRDF
  type(rttov_profile),     pointer :: profiles(:)    => null() ! Input profiles
  type(rttov_transmission)         :: transmission             ! Output transmittances
  type(rttov_radiance)             :: radiance                 ! Output radiances

  integer(kind=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  integer(kind=jpim) :: alloc_status

  ! variables for input
  !====================
  character(len=256) :: coef_filename = './rtcoef_himawari_8_ahi.dat'
  character(len=256) :: sccldcoef_filename = './sccldcoef_himawari_8_ahi.dat'
!  integer(kind=jpim) :: nthreads = 1 ! OFP 
  integer(kind=jpim) :: dosolar
  integer(kind=jpim) :: nlevels
  integer(kind=jpim) :: nprof
  integer(kind=jpim) :: nchannels
  integer(kind=jpim) :: nchanprof
!  integer(kind=jpim), allocatable :: channel_list(:)           
  ! loop variables
  integer(kind=jpim) :: j, jch
  integer(kind=jpim) :: nch
  integer(kind=jpim) :: iprof, joff




!
! minQcfrac: Threshold for diagnosing cloud fraction
  real(kind=jprb) :: jcfrac_cnst
  real(kind=jprb) :: minQcfrac ! threshold water/ice contents (g m-3) for cloud fraction diagnosis

! ###

  real(r_size),intent(in) :: tmp_p(nlevels, nprof)
  real(r_size),intent(in) :: tmp_t(nlevels, nprof)
  real(r_size),intent(in) :: tmp_qv(nlevels, nprof)
  real(r_size),intent(in) :: tmp_qc(nlevels, nprof)
  real(r_size),intent(in) :: tmp_qice(nlevels, nprof)
  real(r_size),intent(in) :: tmp_t2m(nprof)
  real(r_size),intent(in) :: tmp_q2m(nprof)
  real(r_size),intent(in) :: tmp_p2m(nprof)
  real(r_size),intent(in) :: tmp_u2m(nprof)
  real(r_size),intent(in) :: tmp_v2m(nprof)
!  real(r_size),intent(in) :: tmp_soze(nprof)
!  real(r_size),intent(in) :: tmp_soaz(nprof)
!  real(r_size),intent(in) :: tmp_saze(nprof)
!  real(r_size),intent(in) :: tmp_saaz(nprof)
  real(r_size),intent(in) :: tmp_elev(nprof)
  real(r_size),intent(in) :: tmp_lon(nprof)
  real(r_size),intent(in) :: tmp_lat(nprof)
  real(r_size),intent(in) :: tmp_land(nprof)
  real(r_size),intent(in) :: tmp_sat_zangle(nprof)


  real(kind=jprb),allocatable :: kgkg2gm3(:) ! convert parameter [kg/kg] => [gm^-3]
  real(kind=jprb),allocatable :: icec(:) ! ice cloud content (ice + snow + graupel)
  real(kind=jprb),allocatable :: liqc(:) ! liquid cloud content (cloud water)
  real (kind=jprb):: tv ! virtual temp. (K)

  !--------------------------
  !
  integer(kind=jpim), parameter :: iup   = 20   ! unit for input profile file
  integer(kind=jpim), parameter :: ioout = 21   ! unit for output
  integer(kind=jpim), parameter :: mxchn = 9000 ! max number of channels



  ! variables for input
  !====================
  integer(kind=jpim) :: input_chan(mxchn)
  real(kind=jprb)    :: input_ems(mxchn), input_brdf(mxchn)
  integer(kind=jpim) :: ivch, ich
  integer(kind=jpim) :: ilev

! by T.Honda
  real(kind=r_size),intent(out) :: btall_out(nchannels,nprof)
  real(kind=r_size),intent(out) :: btclr_out(nchannels,nprof)
!  real(kind=r_size),intent(out) :: trans_out(nlevels,nchannels,nprof)

  logical :: debug = .false.
!  logical :: debug = .true.

  real(kind=jprb) :: Rd 
  real(kind=jprb) :: Rv 

  real(kind=jprb) :: epsb 
  real(kind=jprb) :: repsb 

  Real(Kind=jprb), Parameter :: q_mixratio_to_ppmv  = 1.60771704e+6_JPRB

  Rd = real(Rdry,kind=jprb)
  Rv = real(Rvap,kind=jprb)
  epsb = Rd / Rv 
  repsb = 1.0_jprb / epsb

  minQcfrac = real(H08_RTTOV_MinQ,kind=jprb)
  jcfrac_cnst = real(H08_RTTOV_CFRAC_CNST,kind=jprb)

  allocate(kgkg2gm3(nlevels))
  allocate(icec(nlevels))
  allocate(liqc(nlevels))

  if(debug) write(6,'(1x,a)')"hello from RTTOV12"
!  if(debug) nprof = 100 ! debug

  errorstatus     = 0_jpim

  dosolar = 0

  ! Read channel list including input surface emissivities and reflectances
  ! If the input emissivities/reflectances are zero or less, calcemis/calcrefl is set to true.

  do ich = 1, nchannels
    input_chan(ich)=ich
  end do
  input_ems(:)=0.0
  input_brdf(:)=0.0

  if(debug) write(6,'(1x,a)')"hello from RTTOV2"

  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  if (dosolar == 1) then
    opts % rt_ir % addsolar = .TRUE.           ! Include solar radiation
  else
    opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
  endif

  opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method
  opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects

  opts % rt_ir % addclouds           = .true.  ! Include cloud effects

  opts % rt_ir % ir_scatt_model      = 2       ! Scattering model for emission source term:
                                               !   1 => doM; 2 => Chou-scaling
  opts % rt_ir % vis_scatt_model     = 1       ! Scattering model for solar source term:
                                               !   1 => doM; 2 =>
                                               !   single-scattering; 3 =>
                                               !   MFASIS
  opts % rt_ir % dom_nstreams        = 8       ! Number of streams for Discrete Ordinates (doM)

  opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
  opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
  opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
  opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
  opts % rt_ir % co_data             = .FALSE. !
  opts % rt_ir % so2_data            = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !

  if (debug) then
    opts % config % verbose            = .true.  ! Enable printing of warnings
  else
    opts % config % verbose            = .false.  ! Enable printing of warnings
  endif


!  opts%config%apply_reg_limits       = .true.
  opts%config%do_checkinput          = .false.

!! added by T.Honda(2015/06/10)
!  opts % interpolation % reg_limit_extrap  = .TRUE.  ! see UG 7.3 (32pp)
!  opts % config % apply_reg_limits  = .TRUE.  ! see UG 7.3 (32pp)
!!  opts%config%do_checkinput = .FALSE. ! see UG 85pp

  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  if(debug) write(6,'(1x,a)')"hello from RTTOV3 IR"
  call rttov_read_coefs(errorstatus, coefs, opts, form_coef='formatted', &
                       file_coef=trim(H08_RTTOV_COEF_PATH)//"/"//coef_filename, &
                       file_sccld=trim(H08_RTTOV_COEF_PATH)//"/"//sccldcoef_filename)
!  CALL rttov_read_coefs(errorstatus, coefs, opts, form_coef='formatted', file_coef=coef_filename)
  if (errorstatus /= errorstatus_success) then
    write(*,*) 'fatal error reading coefficients'
    call rttov_exit(errorstatus)
  endif

!  ! Ensure input number of channels is not higher than number stored in coefficient file
!  if (nchannels > coefs % coef % fmv_chn) THEN
!    nchannels = coefs % coef % fmv_chn
!  endif

  ! Ensure the options and coefficients are consistent
  call rttov_user_options_checkinput(errorstatus, opts, coefs)
  if (errorstatus /= errorstatus_success) then
    write(*,*) 'error in rttov options'
    call rttov_exit(errorstatus)
  endif

  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for rttov_direct
  call rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
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
  if (errorstatus /= errorstatus_success) then
    write(*,*) 'allocation error for rttov_direct structures'
    call rttov_exit(errorstatus)
  endif

  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  do j = 1, nprof
    do jch = 1, nchannels
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = jch ! + NVIS_HIM8 !channel_list(jch)
    enddo
  enddo

  if(debug) write(6,'(1x,a)')"hello from RTTOV5"

  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------
  if(debug) then
    write(6,*) 'debug information'
    write(6,*) 'max p, max t',maxval(tmp_p),maxval(tmp_t)
    write(6,*) 'min p, min t',minval(tmp_p),minval(tmp_t)
    write(6,*) 'max p2, max t2',maxval(tmp_p2m),maxval(tmp_t2m)
    write(6,*) 'min p2, min t2',minval(tmp_p2m),minval(tmp_t2m)

    write(6,*) 'max qv, min qv',maxval(tmp_qv),minval(tmp_qv)
    write(6,*) '-- t prof --'
    write(6,*) 'size t:',size(tmp_t(:,1)),size(tmp_p(:,1))
    write(6,*) '-- zenith --'
    write(6,*) 'min za, min za',minval(tmp_sat_zangle(:)),minval(tmp_sat_zangle(:))

!    do j = 1, nlevels
!      write(6,*) 'qv:',tmp_qv(j,2)
!    enddo
!    write(6,*) '-- end t prof --'
!    do j = 1, nlevels
!      write(6,*) 't:',tmp_t(j,2)
!    enddo
  endif

  !===============================================
  !========== READ profiles == start =============
  if(debug) write(6,*) 'START SUBSTITUTE PROFILE'
  do iprof = 1, nprof

    profiles(iprof)% zenangle = real(tmp_sat_zangle(iprof),kind=jprb)

    !profiles(iprof) % gas_units = 1 ! kg/kg
    profiles(iprof) % gas_units = 2 ! ppmv 

    profiles(iprof)%p(:)=real(tmp_p(:,iprof),kind=jprb) * 0.01_jprb  ! (hpa)
    profiles(iprof)%t(:)=real(tmp_t(:,iprof),kind=jprb)
    !profiles(iprof)%q(:)=max(real(tmp_qv(:,iprof),kind=jprb), qmin) ! (kg/kg) 
    profiles(iprof)%q(:)=min(max(real(tmp_qv(:,iprof)*q_mixratio_to_ppmv,kind=jprb), qmin*1.01_jprb), qmax*0.99_jprb) ! (kg/kg) 
    profiles(iprof)%s2m%t=real(tmp_t2m(iprof),kind=jprb)
    !profiles(iprof)%s2m%q=max(real(tmp_q2m(iprof),kind=jprb), qmin) ! (kg/kg)
    profiles(iprof)%s2m%q=min(max(real(tmp_q2m(iprof),kind=jprb)*q_mixratio_to_ppmv, qmin*1.01_jprb), qmax*0.99_jprb) ! (kg/kg)

!    if(profiles(iprof)%s2m%q < qmin) profiles(iprof)%s2m%q = qmin + qmin * 0.01_jprb
!    do ilev=1,nlevels
!      if(profiles(iprof)%q(ilev) < qmin) profiles(iprof)%q(ilev) = qmin + qmin * 0.01_jprb
!    enddo

    profiles(iprof)%s2m%p=real(tmp_p2m(iprof),kind=jprb) * 0.01_jprb ! (hPa)
    profiles(iprof)%s2m%u=real(tmp_u2m(iprof),kind=jprb)
    profiles(iprof)%s2m%v=real(tmp_v2m(iprof),kind=jprb)
    profiles(iprof)%s2m%wfetc= 100000.0_jprb

    profiles(iprof) % skin % t = real(tmp_t2m(iprof),kind=jprb)
!    profiles(iprof) % skin % fastem(1) = 3.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(2) = 5.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(3) =15.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(4) = 0.1 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(5) = 0.3 ! comment out (11/18/2015)

    profiles(iprof) % skin % surftype = int(tmp_land(iprof))
    profiles(iprof) % skin % watertype = 1 ! tentative (11/18/2015)

    profiles(iprof) % elevation = real(tmp_elev(iprof),kind=jprb) * 0.001_jprb ! (km)
    profiles(iprof) % latitude  = real(tmp_lat(iprof),kind=jprb)
    profiles(iprof) % longitude = real(tmp_lon(iprof),kind=jprb)


!    profiles(iprof)% zenangle = 30.0_jprb ! tentative
!    profiles(iprof)%azangle=0.0_jprb     ! Not required for [opts % rt_ir % addsolar = .FALSE.] 
!    profiles(iprof)%sunzenangle=0.0_jprb ! Not required for [opts % rt_ir % addsolar = .FALSE.] 
!    profiles(iprof)%sunazangle=0.0_jprb  ! Not required for [opts % rt_ir % addsolar = .FALSE.]

!   These are parameters for simple cloud.
!   Not used.
    profiles(iprof) % ctp       = 500.0_jprb
    profiles(iprof) % cfraction = 0.0_jprb
 

  !-- 6 general cloud 
    if( opts % rt_ir % addclouds ) then

      ! Select the CLW and ice cloud properties:
      profiles(:) % clw_scheme = 1
      profiles(:) % ice_scheme = 2

      ! Set the ice Deff parameterisation to a suitable value
      profiles(:) % idg = 4

      profiles(iprof) % cloud(:,:) = 0._jprb
      profiles(iprof) % cfrac(:)   = 0._jprb

!  --- 6.1 microphysical clouds 
!    if(icldprf==1) then
    ! --convert kg/kg into g/m3, make liq.cloud and ice.cloud amount, and upside-down
      do ilev=1,nlevels
        tv = real(tmp_t(ilev,iprof),kind=jprb) * (1.0_jprb+real(tmp_qv(ilev,iprof),kind=jprb) * repsb) &
           / (1.0_jprb + real(tmp_qv(ilev,iprof),kind=jprb))
        kgkg2gm3(ilev) = real(tmp_p(ilev,iprof),kind=jprb) / (Rd * tv) * 1000.0_jprb 
        liqc(ilev) = real(max(tmp_qc(ilev,iprof),0.0_r_size),kind=jprb) * kgkg2gm3(ilev)
        icec(ilev) = real(max(tmp_qice(ilev,iprof),0.0_r_size),kind=jprb) * kgkg2gm3(ilev)
      end do !ilev

      do ilev=1,nlevels-1
        profiles(iprof) % cloud(2,ilev) = & !stratus maritime (default)
                   (liqc(ilev+1) + liqc(ilev)) * 0.5_jprb
        profiles(iprof) % cloud(6,ilev) = & 
                   (icec(ilev+1) + icec(ilev)) * 0.5_jprb
!
! cloud fraction diagnosis
        if(jcfrac_cnst <= 0.0_jprb)then
          if(((profiles(iprof) % cloud(2,ilev) + &
               profiles(iprof) % cloud(6,ilev)) > minQcfrac)) then
            profiles(iprof) % cfrac(ilev)   = 1.0_jprb  !cloud fraction 
          else
            profiles(iprof) % cfrac(ilev)   = 0.0_jprb  !cloud fraction 
          endif
        else
          profiles(iprof) % cfrac(ilev) = min((profiles(iprof) % cloud(2,ilev) + &
                                               profiles(iprof) % cloud(6,ilev)) / jcfrac_cnst, 1.0_jprb)
        endif
      end do ! ilev
    endif ! addclouds
    if (debug) write(6,*) 'qmax, qmin',maxval(profiles(iprof)%q(:)), minval(profiles(iprof)%q(:))
  enddo ! prof

  deallocate(kgkg2gm3, icec, liqc)


  if(debug) write(6,*) 'END SUBSTITUTE PROFILE'

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

  if(debug) write(6,*)"Enter direct"

  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  if (H08_RTTOV_NTHREAD <= 1) then
    call rttov_direct(                &
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
  else
    call rttov_parallel_direct(       &
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
            nthreads    = H08_RTTOV_NTHREAD)    ! in    number of threads to use
  endif

  if (errorstatus /= errorstatus_success) then
    write (*,*) 'rttov_direct error'
    call rttov_exit(errorstatus)
  endif

  ! --- Output the results --------------------------------------------------

  do iprof = 1, nprof 

    joff = (iprof-1_jpim) * nchannels

    !
    !     outPUT RESULTS
    !

    btall_out(1:nchannels,iprof) = real(radiance%bt(1+joff:nchannels+joff), kind=r_size)
    btclr_out(1:nchannels,iprof) = real(radiance%bt_clear(1+joff:nchannels+joff), kind=r_size)

    if (debug) print *,maxval(btall_out(:,iprof)), minval(btall_out(:,iprof))

!    do ilev = 1, nlevels
!      trans_out(ilev,1:nchannels,iprof) = real(transmission % tau_levels(ilev,1+joff:nchannels+joff), kind=r_size)
      !if(debug .and. mod(iprof,5) == 0) write(6,'(a,f10.7,i4)')"RTTOV debug trans",trans_out(ilev,1,iprof)
!    enddo

  enddo

  !============== Output results == end ==============
  !=====================================================

  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
!  deallocate (channel_list, stat=alloc_status)
!  if (alloc_status /= 0) THEN
!    write(*,*) 'mem dellocation error'
!  endif

  ! Deallocate structures for rttov_direct
  call rttov_alloc_direct( &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance)
  if (errorstatus /= errorstatus_success) then
    write(*,*) 'deallocation error for rttov_direct structures'
    call rttov_exit(errorstatus)
  endif

  if(debug) write(*,*)'Successfully finished IR!!'

return
end subroutine cld_ir_fwd

subroutine cld_mfasis_fwd(nchannels,&
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
                           tmp_elev,&
                           tmp_lon,&
                           tmp_lat,&
                           tmp_land,&
                           tmp_sat_zangle, & ! satellite zenith angle
                           tmp_sat_azm, & ! satellite azimuth angle
                           tmp_szangle,  & ! solar zenith angle
                           tmp_saangle,  &
                           reflall_out,& 
                           reflclr_out)!,& 
!                           trans_out)

  ! rttov_const contains useful RTTOV constants
  use rttov_const, only :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         qmin,                &
         qmax

  ! rttov_types contains definitions of all RTTOV data types
  use rttov_types, only :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

  ! jpim, jprb and jplm are the RTTOV integer, real and logical kinds
  use parkind1, only : jpim, jprb, jplm

  use rttov_unix_env, only : rttov_exit

  use common, only : r_size
  use scale_const, only: &
        Rdry    => CONST_Rdry, &
        Rvap    => CONST_Rvap, &
        Deg2Rad => CONST_D2R
  use common_nml, only: &
        H08_RTTOV_COEF_PATH, &
        H08_RTTOV_MinQ, &
        H08_RTTOV_CFRAC_CNST, &
        H08_RTTOV_CLD, &
        H08_RTTOV_NTHREAD
  implicit none

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_print_radiance_quality.interface"
#include "rttov_skipcommentline.interface"

  ! RTTOV variables/structures
  !====================
  type(rttov_options)              :: opts                     ! Options structure
  type(rttov_coefs)                :: coefs                    ! Coefficients structure
  type(rttov_chanprof),    pointer :: chanprof(:)    => null() ! Input channel/profile list
  logical(kind=jplm),      pointer :: calcemis(:)    => null() ! Flag to indicate calculation of emissivity within RTTOV
  type(rttov_emissivity),  pointer :: emissivity(:)  => null() ! Input/output surface emissivity
  logical(kind=jplm),      pointer :: calcrefl(:)    => null() ! Flag to indicate calculation of BRDF within RTTOV
  type(rttov_reflectance), pointer :: reflectance(:) => null() ! Input/output surface BRDF
  type(rttov_profile),     pointer :: profiles(:)    => null() ! Input profiles
  type(rttov_transmission)         :: transmission             ! Output transmittances
  type(rttov_radiance)             :: radiance                 ! Output radiances

  integer(kind=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  integer(kind=jpim) :: alloc_status

  ! variables for input
  !====================
!  character(len=256) :: sccldcoef_filename = './sccldcoef_himawari_8_ahi.dat'
  character(len=256) :: coef_filename = './rtcoef_himawari_8_ahi.dat'
  character(len=256) :: sccldcoef_filename = './sccldcoef_himawari_8_ahi.dat'
  character(len=256) :: mfasis_lut_filename = './rttov_mfasis_cld_himawari_8_ahi_opac.H5'
!  integer(kind=jpim) :: nthreads = 1 ! K
  integer(kind=jpim) :: nlevels
  integer(kind=jpim) :: nprof
  integer(kind=jpim) :: nchannels
  integer(kind=jpim) :: nchanprof
!  integer(kind=jpim), allocatable :: channel_list(:)           
  ! loop variables
  integer(kind=jpim) :: j, jch
  integer(kind=jpim) :: nch
  integer(kind=jpim) :: iprof, joff




!
! minQcfrac: Threshold for diagnosing cloud fraction
  real(kind=jprb) :: jcfrac_cnst
  real(kind=jprb) :: minQcfrac ! threshold water/ice contents (g m-3) for cloud fraction diagnosis

! ###

  real(r_size),intent(in) :: tmp_p(nlevels, nprof)
  real(r_size),intent(in) :: tmp_t(nlevels, nprof)
  real(r_size),intent(in) :: tmp_qv(nlevels, nprof)
  real(r_size),intent(in) :: tmp_qc(nlevels, nprof)
  real(r_size),intent(in) :: tmp_qice(nlevels, nprof)
  real(r_size),intent(in) :: tmp_t2m(nprof)
  real(r_size),intent(in) :: tmp_q2m(nprof)
  real(r_size),intent(in) :: tmp_p2m(nprof)
  real(r_size),intent(in) :: tmp_u2m(nprof)
  real(r_size),intent(in) :: tmp_v2m(nprof)
!  real(r_size),intent(in) :: tmp_soze(nprof)
!  real(r_size),intent(in) :: tmp_soaz(nprof)
!  real(r_size),intent(in) :: tmp_saze(nprof)
!  real(r_size),intent(in) :: tmp_saaz(nprof)
  real(r_size),intent(in) :: tmp_elev(nprof)
  real(r_size),intent(in) :: tmp_lon(nprof)
  real(r_size),intent(in) :: tmp_lat(nprof)
  real(r_size),intent(in) :: tmp_land(nprof)
  real(r_size),intent(in) :: tmp_sat_zangle(nprof)
  real(r_size),intent(in) :: tmp_sat_azm(nprof)
  real(r_size),intent(in) :: tmp_szangle(nprof)
  real(r_size),intent(in) :: tmp_saangle(nprof)


  real(kind=jprb),allocatable :: kgkg2gm3(:) ! convert parameter [kg/kg] => [gm^-3]
  real(kind=jprb),allocatable :: icec(:) ! ice cloud content (ice + snow + graupel)
  real(kind=jprb),allocatable :: liqc(:) ! liquid cloud content (cloud water)
  real (kind=jprb):: tv ! virtual temp. (K)

  !--------------------------
  !
  integer(kind=jpim), parameter :: iup   = 20   ! unit for input profile file
  integer(kind=jpim), parameter :: ioout = 21   ! unit for output
  integer(kind=jpim), parameter :: mxchn = 9000 ! max number of channels



  ! variables for input
  !====================
  integer(kind=jpim) :: input_chan(mxchn)
  real(kind=jprb)    :: input_ems(mxchn), input_brdf(mxchn)
  integer(kind=jpim) :: ivch, ich
  integer(kind=jpim) :: ilev

! by T.Honda
  real(kind=r_size),intent(out) :: reflall_out(nchannels,nprof)
  real(kind=r_size),intent(out) :: reflclr_out(nchannels,nprof)
!  real(kind=r_size),intent(out) :: trans_out(nlevels,nchannels,nprof)

  logical :: debug = .false.
!  logical :: debug = .true.

  real(kind=jprb) :: Rd 
  real(kind=jprb) :: Rv 

  real(kind=jprb) :: epsb 
  real(kind=jprb) :: repsb 

  Real(Kind=jprb), Parameter :: q_mixratio_to_ppmv  = 1.60771704e+6_JPRB

  logical :: unit_kgkg = .true.

  Rd = real(Rdry,kind=jprb)
  Rv = real(Rvap,kind=jprb)
  epsb = Rd / Rv 
  repsb = 1.0_jprb / epsb

  minQcfrac = real(H08_RTTOV_MinQ,kind=jprb)
  jcfrac_cnst = real(H08_RTTOV_CFRAC_CNST,kind=jprb)

  if ( unit_kgkg ) then
    ! g/kg => kg/kg
    minQcfrac = minQcfrac * 1.E-3_jprb
    jcfrac_cnst = jcfrac_cnst * 1.E-3_jprb
  endif

  allocate(kgkg2gm3(nlevels))
  allocate(icec(nlevels))
  allocate(liqc(nlevels))

  if(debug) write(6,'(1x,a)')"hello from RTTOV12"
!  if(debug) nprof = 100 ! debug

  errorstatus     = 0_jpim

  ! Read channel list including input surface emissivities and reflectances
  ! If the input emissivities/reflectances are zero or less, calcemis/calcrefl is set to true.

  do ich = 1, nchannels
    input_chan(ich)=ich
  end do
  input_ems(:)=0.0
  input_brdf(:)=0.0

  if(debug) write(6,'(1x,a)')"hello from RTTOV2"

  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  opts % rt_ir % addsolar = .TRUE.           ! Include solar radiation
  opts % rt_ir % addclouds           = .true.  ! Include cloud effects
  opts % rt_ir % vis_scatt_model     = 3       ! Scattering model for solar source term must be set
                                               !   to 3 (MFASIS) before reading
                                               !   coefficient files

  opts % rt_ir % ir_scatt_model      = 2       ! Scattering model for emission source term:
                                               !   1 => doM; 2 => Chou-scaling

  opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method

  if (debug) then
    opts % config % verbose            = .true.  ! Enable printing of warnings
  else
    opts % config % verbose            = .false.  ! Enable printing of warnings
  endif

  opts % rt_ir % dom_nstreams        = 8       ! Number of streams for Discrete Ordinates (doM)

  opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
  opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
  opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
  opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
  opts % rt_ir % co_data             = .FALSE. !
  opts % rt_ir % so2_data            = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !


!  opts%config%apply_reg_limits       = .true.
!  opts%config%do_checkinput          = .false.

!! added by T.Honda(2015/06/10)
!  opts % interpolation % reg_limit_extrap  = .TRUE.  ! see UG 7.3 (32pp)
!  opts % config % apply_reg_limits  = .TRUE.  ! see UG 7.3 (32pp)
!!  opts%config%do_checkinput = .FALSE. ! see UG 85pp

  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  if(debug) write(6,'(1x,a)')"hello from RTTOV3 MFASIS"
  call rttov_read_coefs(errorstatus, coefs, opts, & !form_coef='formatted', &
                       file_coef=trim(H08_RTTOV_COEF_PATH)//"/vis/"//coef_filename, &
                       file_sccld=trim(H08_RTTOV_COEF_PATH)//"/vis/"//sccldcoef_filename, &
                       file_mfasis_cld=trim(H08_RTTOV_COEF_PATH)//"/vis/"//mfasis_lut_filename)

!  CALL rttov_read_coefs(errorstatus, coefs, opts, form_coef='formatted', file_coef=coef_filename)
  if (errorstatus /= errorstatus_success) then
    write(*,*) 'fatal error reading coefficients'
    call rttov_exit(errorstatus)
  endif

!  ! Ensure input number of channels is not higher than number stored in coefficient file
!  if (nchannels > coefs % coef % fmv_chn) THEN
!    nchannels = coefs % coef % fmv_chn
!  endif

  ! Ensure the options and coefficients are consistent
  call rttov_user_options_checkinput(errorstatus, opts, coefs)
  if (errorstatus /= errorstatus_success) then
    write(*,*) 'error in rttov options'
    call rttov_exit(errorstatus)
  endif

  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for rttov_direct
  call rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
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
  if (errorstatus /= errorstatus_success) then
    write(*,*) 'allocation error for rttov_direct structures'
    call rttov_exit(errorstatus)
  endif

  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  do j = 1, nprof
    do jch = 1, nchannels ! Only for VIS
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = jch !channel_list(jch)
    enddo
  enddo

  if(debug) write(6,'(1x,a)')"hello from RTTOV5"

  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------
  if(debug) then
    write(6,*) 'debug information'
    write(6,*) 'max p, max t',maxval(tmp_p),maxval(tmp_t)
    write(6,*) 'min p, min t',minval(tmp_p),minval(tmp_t)
    write(6,*) 'max p2, max t2',maxval(tmp_p2m),maxval(tmp_t2m)
    write(6,*) 'min p2, min t2',minval(tmp_p2m),minval(tmp_t2m)

    write(6,*) 'max qv, min qv',maxval(tmp_qv),minval(tmp_qv)
    write(6,*) '-- t prof --'
    write(6,*) 'size t:',size(tmp_t(:,1)),size(tmp_p(:,1))

!    do j = 1, nlevels
!      write(6,*) 'qv:',tmp_qv(j,2)
!    enddo
!    write(6,*) '-- end t prof --'
!    do j = 1, nlevels
!      write(6,*) 't:',tmp_t(j,2)
!    enddo
  endif

  !===============================================
  !========== READ profiles == start =============
  if(debug) write(6,*) 'START SUBSTITUTE PROFILE'
  do iprof = 1, nprof

    profiles(iprof) % sunzenangle = real(tmp_szangle(iprof),kind=jprb)
    profiles(iprof) % sunazangle = real(tmp_saangle(iprof),kind=jprb)
    profiles(iprof)% zenangle = real(tmp_sat_zangle(iprof),kind=jprb)
    profiles(iprof)% azangle = real(tmp_sat_azm(iprof),kind=jprb)

    profiles(iprof)%p(:)=real(tmp_p(:,iprof),kind=jprb) * 0.01_jprb  ! (hpa)
    profiles(iprof)%t(:)=real(tmp_t(:,iprof),kind=jprb)
    profiles(iprof)%q(:)=max(real(tmp_qv(:,iprof),kind=jprb), qmin) ! (kg/kg) 
    !profiles(iprof)%q(:)=min(max(real(tmp_qv(:,iprof)*q_mixratio_to_ppmv,kind=jprb), qmin*1.01_jprb), qmax*0.99_jprb) ! (kg/kg) 
    profiles(iprof)%s2m%t=real(tmp_t2m(iprof),kind=jprb)
    profiles(iprof)%s2m%q=max(real(tmp_q2m(iprof),kind=jprb), qmin) ! (kg/kg)
    !profiles(iprof)%s2m%q=min(max(real(tmp_q2m(iprof),kind=jprb)*q_mixratio_to_ppmv, qmin*1.01_jprb), qmax*0.99_jprb) ! (kg/kg)

!    if(profiles(iprof)%s2m%q < qmin) profiles(iprof)%s2m%q = qmin + qmin * 0.01_jprb
!    do ilev=1,nlevels
!      if(profiles(iprof)%q(ilev) < qmin) profiles(iprof)%q(ilev) = qmin + qmin * 0.01_jprb
!    enddo

    profiles(iprof)%s2m%p=real(tmp_p2m(iprof),kind=jprb) * 0.01_jprb ! (hPa)
    profiles(iprof)%s2m%u=real(tmp_u2m(iprof),kind=jprb)
    profiles(iprof)%s2m%v=real(tmp_v2m(iprof),kind=jprb)
    profiles(iprof)%s2m%wfetc= 100000.0_jprb

    if ( unit_kgkg ) then
      profiles(iprof) % gas_units = 1 ! kg/kg
      profiles(iprof) % mmr_cldaer = .true. ! kg/kg
    else
      profiles(iprof) % gas_units = 2 ! ppmv 
      write(*, *) "unit_kgkg should be true"
      stop
    endif


    profiles(iprof) % skin % t = real(tmp_t2m(iprof),kind=jprb)
!    profiles(iprof) % skin % fastem(1) = 3.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(2) = 5.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(3) =15.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(4) = 0.1 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(5) = 0.3 ! comment out (11/18/2015)

    profiles(iprof) % skin % surftype = int(tmp_land(iprof))
    profiles(iprof) % skin % watertype = 1 ! tentative (11/18/2015)

    profiles(iprof) % elevation = real(tmp_elev(iprof),kind=jprb) * 0.001_jprb ! (km)
    profiles(iprof) % latitude  = real(tmp_lat(iprof),kind=jprb)
    profiles(iprof) % longitude = real(tmp_lon(iprof),kind=jprb)


!    profiles(iprof)%azangle=0.0_jprb     ! Not required for [opts % rt_ir % addsolar = .FALSE.] 
!    profiles(iprof)%sunzenangle=0.0_jprb ! Not required for [opts % rt_ir % addsolar = .FALSE.] 
!    profiles(iprof)%sunazangle=0.0_jprb  ! Not required for [opts % rt_ir % addsolar = .FALSE.]

!   These are parameters for simple cloud.
!   Not used.
    profiles(iprof) % ctp       = 500.0_jprb
    profiles(iprof) % cfraction = 0.0_jprb
 

  !-- 6 general cloud 
    if( opts % rt_ir % addclouds ) then

      ! Select the CLW and ice cloud properties:
      profiles(:) % clw_scheme = 1
      profiles(:) % ice_scheme = 1 ! SSEC ice properties

      ! Set the ice Deff parameterisation to a suitable value
      profiles(:) % idg = 4

      profiles(iprof) % cloud(:,:) = 0._jprb
      profiles(iprof) % cfrac(:)   = 0._jprb

!  --- 6.1 microphysical clouds 
!    if(icldprf==1) then
    ! --convert kg/kg into g/m3, make liq.cloud and ice.cloud amount, and upside-down
      do ilev=1,nlevels
        tv = real(tmp_t(ilev,iprof),kind=jprb) * (1.0_jprb+real(tmp_qv(ilev,iprof),kind=jprb) * repsb) &
           / (1.0_jprb + real(tmp_qv(ilev,iprof),kind=jprb))
        kgkg2gm3(ilev) = real(tmp_p(ilev,iprof),kind=jprb) / (Rd * tv) * 1000.0_jprb 
        liqc(ilev) = real(max(tmp_qc(ilev,iprof),0.0_r_size),kind=jprb) !* kgkg2gm3(ilev)
        icec(ilev) = real(max(tmp_qice(ilev,iprof),0.0_r_size),kind=jprb) !* kgkg2gm3(ilev)
      end do !ilev

      do ilev=1,nlevels-1
        profiles(iprof) % cloud(2,ilev) = & !stratus maritime (default)
                   (liqc(ilev+1) + liqc(ilev)) * 0.5_jprb 
        profiles(iprof) % cloud(6,ilev) = & 
                   (icec(ilev+1) + icec(ilev)) * 0.5_jprb 
!
! cloud fraction diagnosis
        if(jcfrac_cnst <= 0.0_jprb)then
          if(((profiles(iprof) % cloud(2,ilev) + &
               profiles(iprof) % cloud(6,ilev)) > minQcfrac)) then
            profiles(iprof) % cfrac(ilev)   = 1.0_jprb  !cloud fraction 
          else
            profiles(iprof) % cfrac(ilev)   = 0.0_jprb  !cloud fraction 
          endif
        else
          profiles(iprof) % cfrac(ilev) = min((profiles(iprof) % cloud(2,ilev) + &
                                               profiles(iprof) % cloud(6,ilev)) / jcfrac_cnst, 1.0_jprb)
        endif
      end do ! ilev
    endif ! addclouds
    if (debug .and. (mod(iprof,50)==0)) write(6,*) 'qmax, qmin',maxval(profiles(iprof)%q(:)), minval(profiles(iprof)%q(:))
  enddo ! prof

  deallocate(kgkg2gm3, icec, liqc)


  if(debug) write(6,*) 'END SUBSTITUTE PROFILE'

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

  if(debug) write(6,*)"Enter direct"

  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  if (H08_RTTOV_NTHREAD <= 1) then
    call rttov_direct(                &
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
  else
    call rttov_parallel_direct(       &
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
            nthreads    = H08_RTTOV_NTHREAD)    ! in    number of threads to use
  endif

  if (errorstatus /= errorstatus_success) then
    write (*,*) 'rttov_direct error'
    call rttov_exit(errorstatus)
  endif

  ! --- Output the results --------------------------------------------------

  do iprof = 1, nprof 

    joff = (iprof-1_jpim) * nchannels

    !
    !     outPUT RESULTS
    !

    ! VIS/NIR Him8
    reflall_out(1:nchannels,iprof) = real(radiance%refl(1+joff:nchannels+joff), kind=r_size)
    reflclr_out(1:nchannels,iprof) = real(radiance%refl_clear(1+joff:nchannels+joff), kind=r_size)

!    if (mod(iprof,200) == 0) then
!print *,"DEBUG refl ",radiance%refl(1+joff),radiance%refl(2+joff),radiance%refl(3+joff),radiance%refl(4+joff),radiance%refl(5+joff),radiance%refl(6+joff)!,radiance%refl(7+joff)
!print *,"DEBUG bt ",radiance%bt(1+joff),radiance%bt(2+joff),radiance%bt(3+joff),radiance%bt(4+joff),radiance%bt(5+joff),radiance%bt(6+joff)!,radiance%bt(7+joff)
!    endif

  enddo


  !============== Output results == end ==============
  !=====================================================

  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
!  deallocate (channel_list, stat=alloc_status)
!  if (alloc_status /= 0) THEN
!    write(*,*) 'mem dellocation error'
!  endif

  ! Deallocate structures for rttov_direct
  call rttov_alloc_direct( &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance)
  if (errorstatus /= errorstatus_success) then
    write(*,*) 'deallocation error for rttov_direct structures'
    call rttov_exit(errorstatus)
  endif

  if(debug) write(*,*)'Successfully finished!!'

return
end subroutine cld_mfasis_fwd


end module scale_rttov12_fwd

