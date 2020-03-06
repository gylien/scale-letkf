module scale_H08_fwd12
!$USE OMP_LIB

  USE common, ONLY : r_size
  USE scale_const, ONLY: &
        Rdry    => CONST_Rdry,  &
        Rvap    => CONST_Rvap,  &
        Deg2Rad => CONST_D2R,   &
        temp00  => CONST_TEM00, & 
        pres00  => CONST_PRE00, &
        Q_EPS   => CONST_EPS,   &
        CONST_GRAV
  implicit none

contains

subroutine SCALE_RTTOV_fwd12(nchannels,&
                           nlevs,&
                           nprof,&
                           prs,&
                           tk,&
                           qv,&
                           qc,&
                           qice,&
                           tk2m,&
                           q2m,&
                           prs2m,&
                           u2m,&
                           v2m,&
                           elev,&
                           lon,&
                           lat,&
                           land,&
                           zenith,&
!                           RD_presh, &
!                           RD_temph, &
                           btall_out,& 
                           btclr_out,& 
                           mwgt_plev,&
                           ctop_out)

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & platform_name,       &
       & inst_name,           &
!       & q_mixratio_to_ppmv,  &
       & qmin,                &
       & qmax,                &
       & tmin

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
       & rttov_options,       &
       & rttov_coefs,         &
       & rttov_profile,       &
       & rttov_transmission,  &
       & rttov_radiance,      &
       & rttov_chanprof,      &
       & rttov_emissivity,    &
       & rttov_reflectance

  USE rttov_unix_env, ONLY : rttov_exit
  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm
  USE common_nml, ONLY: &
        H08_RTTOV_CFRAC_CNST, &
        H08_RTTOV_MINQ_CTOP,  &
        H08_RTTOV_COEF_PATH,  &
!        H08_RTTOV_PROF_SHIFT, &
        H08_RTTOV_CFRAC,      &
        H08_RTTOV_KADD
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

  integer, intent(in) :: nprof
  integer, intent(in) :: nlevs

  real(r_size), intent(in) :: prs(nlevs, nprof)
  real(r_size), intent(in) :: tk(nlevs, nprof)
  real(r_size), intent(in) :: qv(nlevs, nprof)
  real(r_size), intent(in) :: qc(nlevs, nprof)
  real(r_size), intent(in) :: qice(nlevs, nprof)
  real(r_size), intent(in) :: tk2m(nprof)
  real(r_size), intent(in) :: q2m(nprof)
  real(r_size), intent(in) :: prs2m(nprof)
  real(r_size), intent(in) :: u2m(nprof)
  real(r_size), intent(in) :: v2m(nprof)
  real(r_size), intent(in) :: elev(nprof)
  real(r_size), intent(in) :: lon(nprof)
  real(r_size), intent(in) :: lat(nprof)
  real(r_size), intent(in) :: land(nprof)
  real(r_size), intent(in) :: zenith(nprof)


  real(kind=jprb) :: icec1, icec2 ! ice cloud content (ice + snow + graupel)
  real(kind=jprb) :: liqc1, liqc2 ! liquid cloud content (cloud water)

  ! variables for input
  !====================
  character(len=256) :: coef_filename='/rtcoef_himawari_8_ahi.dat'
  character(len=256) :: sccoef_filename='/sccldcoef_himawari_8_ahi.dat'
  integer(kind=jpim) :: nthreads = 1_jpim
  integer(kind=jpim) :: dosolar = 0_jpim
  integer(kind=jpim), intent(in) :: nchannels
  integer(kind=jpim) :: nchanprof
  integer(kind=jpim) :: ich
  ! loop variables
  integer(kind=jpim) :: j, jch
  integer(kind=jpim) :: nch
  integer(kind=jpim) :: ilev
  integer(kind=jpim) :: iprof, joff

! by T.Honda
  real(kind=r_size), intent(out) :: btall_out(nchannels,nprof)
  real(kind=r_size), intent(out) :: btclr_out(nchannels,nprof)
  real(kind=r_size), intent(out) :: ctop_out(nprof)
  real(kind=r_size) :: ptmp, tktmp, qvtmp

  real(kind=r_size) :: rdp, max_wgt, tmp_wgt
  real(kind=r_size), intent(out) :: mwgt_plev(nchannels,nprof) ! Max weight level (Pa)

  logical :: debug = .false.
!  logical :: debug = .true.

  real(kind=jprb) :: repsb 

  integer :: orgk

!  real(kind=r_size), intent(in)  :: RD_presh(nlevs+H08_RTTOV_KADD+1)
!  real(kind=r_size), intent(in)  :: RD_temph(nlevs+H08_RTTOV_KADD+1)
!  real(kind=jprb) :: tmp_dif
 
  real(Kind=jprb), parameter :: q_mixratio_to_ppmv  = 1.60771704e+6_JPRB
 
  if(debug) write(6,'(1x,a)')"hello from RTTOV"

! -- set thermodynamic constants

  repsb = real(Rvap / Rdry, kind=jprb)

  errorstatus     = 0_jpim

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

  opts%rt_ir%ir_sea_emis_model       = 2       ! IREMIS for IR emissivity

  if (debug) then
    opts % config % verbose            = .true.  ! Enable printing of warnings
  else
    opts % config % verbose            = .false.  ! Enable printing of warnings
  endif

  opts%config%apply_reg_limits       = .true.
  opts%config%do_checkinput          = .false.

!! added by T.Honda(2015/06/10)
!  opts % interpolation % reg_limit_extrap  = .TRUE.  ! see UG 7.3 (32pp)
!  opts % config % apply_reg_limits  = .TRUE.  ! see UG 7.3 (32pp)
!!  opts%config%do_checkinput = .FALSE. ! see UG 85pp

  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  if(debug) write(6,'(1x,a)')"hello from RTTOV3"
  call rttov_read_coefs(errorstatus, coefs, opts, form_coef='formatted', &
                       &file_coef=trim(H08_RTTOV_COEF_PATH)//trim(coef_filename), &
                       &file_sccld=trim(H08_RTTOV_COEF_PATH)//trim(sccoef_filename))
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    call rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
!  IF (nchannels > coefs % coef % fmv_chn) THEN
!    nchannels = coefs % coef % fmv_chn
!  ENDIF


  ! Ensure the options and coefficients are consistent
  call rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    call rttov_exit(errorstatus)
  ENDIF

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
        nlevs+H08_RTTOV_KADD,    &
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
      chanprof(nch)%chan = jch !+ NVIS_HIM8 !channel_list(jch)
    enddo
  enddo

  if(debug) write(6,'(1x,a)')"hello from RTTOV5"

  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------
  if(debug) then
    write(6,*) 'debug information'
    write(6,*) 'max p, max t',maxval(prs),maxval(tk)
    write(6,*) 'min p, min t',minval(prs),minval(tk)
    write(6,*) 'max p2, max t2',maxval(prs2m),maxval(tk2m)
    write(6,*) 'min p2, min t2',minval(prs2m),minval(tk2m)

    write(6,*) 'max qv, min qv',maxval(qv),minval(qv)
    write(6,*) '-- t prof --'
    write(6,*) 'size t:',size(tk(:,1)),size(prs(:,1))

!    do j = 1, nlevels
!      write(6,*) 'qv:',qv(j,2)
!    enddo
!    write(6,*) '-- end t prof --'
!    do j = 1, nlevels
!      write(6,*) 't:',tk(j,2)
!    enddo
  endif

  !===============================================
  !========== READ profiles == start =============
  if(debug) write(6,*) 'START SUBSTITUTE PROFILE'
  do iprof = 1, nprof

!    if(H08_RTTOV_PROF_SHIFT)then
!      tmp_dif = RD_temph(H08_RTTOV_KADD+1) - tk(1,iprof)
!    else
!      tmp_dif = 0.0_jprb
!    endif ! H08_RTTOV_PROF_SHIFT
!
!    ! Above the model top (p < RD_presh(H08_RTTOV_KADD+1))
!    do ilev = 1, H08_RTTOV_KADD
!      profiles(iprof)%p(ilev) = real(RD_presh(ilev),kind=jprb) ! (hPa)
!      profiles(iprof)%t(ilev) = max(RD_temph(ilev) - tmp_dif,tmin * 1.01_jprb) ! (K)
!
!      profiles(iprof)%q(ilev) = qmin * 1.01_jprb
!
!    enddo

    do ilev = H08_RTTOV_KADD + 1, H08_RTTOV_KADD + nlevs

      profiles(iprof)%p(ilev) = real(prs(ilev-H08_RTTOV_KADD,iprof),kind=jprb) * 0.01_jprb ! (hPa)
      profiles(iprof)%t(ilev) = real(tk(ilev-H08_RTTOV_KADD,iprof),kind=jprb) ! (K) 

      profiles(iprof)%q(ilev) = min(max(qv(ilev-H08_RTTOV_KADD,iprof) * q_mixratio_to_ppmv, qmin * 1.01_jprb), qmax*0.99) ! (ppmv)

    enddo

    profiles(iprof)%s2m%t = real(tk2m(iprof),kind=jprb)
    profiles(iprof)%s2m%q = real(q2m(iprof),kind=jprb) * q_mixratio_to_ppmv ! (ppmv)

    if(profiles(iprof)%s2m%t < tmin) profiles(iprof)%s2m%t = tmin + tmin * 0.01_jprb
    if(profiles(iprof)%s2m%q < qmin) profiles(iprof)%s2m%q = qmin + qmin * 0.01_jprb


    profiles(iprof)%s2m%p = real(prs2m(iprof),kind=jprb) * 0.01_jprb ! (hPa)
    profiles(iprof)%s2m%u = real(u2m(iprof),kind=jprb)
    profiles(iprof)%s2m%v = real(v2m(iprof),kind=jprb)
    profiles(iprof)%s2m%wfetc = 100000.0_jprb

    profiles(iprof) % skin % t = max(real(tk2m(iprof),kind=jprb), tmin + tmin * 0.01_jprb)



    profiles(iprof)% zenangle = real(zenith(iprof),kind=jprb)

    !profiles(iprof) % gas_units = 1 ! kg/kg
    profiles(iprof) % gas_units = 2 ! ppmv 

    profiles(iprof) % mmr_cldaer = .true. ! kg/kg

    profiles(iprof) % skin % t = real(tk2m(iprof),kind=jprb)

    profiles(iprof) % skin % surftype = int(land(iprof))
    profiles(iprof) % skin % watertype = 1 ! tentative (11/18/2015)

    profiles(iprof) % elevation = real(elev(iprof),kind=jprb) * 0.001_jprb ! (km)
    profiles(iprof) % latitude  = real(lat(iprof),kind=jprb)
    profiles(iprof) % longitude = real(lon(iprof),kind=jprb)


    profiles(iprof)% zenangle = real(zenith(iprof),kind=jprb) ! (11/18/2015)
    if(mod(iprof,1000) == 0 .and. debug)write(6,'(a,f15.10)'),'zenangle ',profiles(iprof)% zenangle
    if(mod(iprof,1000) == 0 .and. debug)write(6,'(a,2f10.5)'),' ',lon(iprof),lat(iprof)



!    profiles(iprof)% zenangle = 30.0_jprb ! tentative
!    profiles(iprof)%azangle=0.0_jprb     ! Not required for [opts % rt_ir %
!    addsolar = .FALSE.] 
!    profiles(iprof)%sunzenangle=0.0_jprb ! Not required for [opts % rt_ir %
!    addsolar = .FALSE.] 
!    profiles(iprof)%sunazangle=0.0_jprb  ! Not required for [opts % rt_ir %
!    addsolar = .FALSE.]

!   These are parameters for simple cloud.
!   Not used.
    profiles(iprof) % ctp       = 500.0_jprb
    profiles(iprof) % cfraction = 0.0_jprb

  !-- 6 general cloud 
    if( opts % rt_ir % addclouds ) then

      ! Select the CLW and ice cloud properties:
      profiles(:) % clw_scheme = 1
      profiles(:) % ice_scheme = 2 ! Baran


      profiles(iprof) % cloud(:,:) = 0._jprb
      profiles(iprof) % cfrac(:)   = 0._jprb


      ! Set the ice Deff parameterisation to a suitable value
      !profiles(:) % idg = 4
      !profiles(iprof) % ish = icecld_ish  !ice water shape

      !profiles(iprof) % idg = icecld_idg  !ice water effective diameter 
      !profiles(iprof) % icede(:)= 0._jprb !ice effective diameter, set non-zero if you give by yourself

      profiles(iprof) % cloud(:,:) = 0._jprb
      profiles(iprof) % cfrac(:)   = 0._jprb


      ctop_out(iprof) = -1.0d0

      do ilev = H08_RTTOV_KADD + 1, H08_RTTOV_KADD + nlevs - 1
        orgk = ilev - H08_RTTOV_KADD ! k index for original profile

        ! ilev
        liqc1 = real(max(qc(orgk,iprof),0.0_r_size),kind=jprb)
        icec1 = real(max(qice(orgk,iprof),0.0_r_size),kind=jprb) 

        ! ilev + 1

        liqc2 = real(max(qc(orgk+1,iprof),0.0_r_size),kind=jprb) 
        icec2 = real(max(qice(orgk+1,iprof),0.0_r_size),kind=jprb)

        profiles(iprof) % cloud(2,ilev) = & !stratus maritime (default)
                   (liqc1 + liqc2) * 0.5_jprb
        profiles(iprof) % cloud(6,ilev) = &
                   (icec1 + icec2) * 0.5_jprb

        ptmp = (prs(orgk+1,iprof) + prs(orgk,iprof))*0.5_jprb    ! (Pa)
        tktmp = (tk(orgk+1,iprof) + tk(orgk,iprof))*0.5_jprb ! (K)
        qvtmp = max((qv(orgk+1,iprof) + qv(orgk,iprof)) * 0.5_jprb, 0.0_r_size) ! (kgkg-1)

        !
        ! cloud fraction & cloud top diagnosis
        !
        select case (H08_RTTOV_CFRAC)
        case (0) ! use H08_RTTOV_MINQ_CTOP as in Honda et al. (2017a,b)
                 !
          profiles(iprof) % cfrac(ilev) = min((profiles(iprof) % cloud(2,ilev) + &
                                               profiles(iprof) % cloud(6,ilev)) / H08_RTTOV_CFRAC_CNST, &
                                               1.0_jprb)

        case (1) ! SCALE microphysics method with a minor modification
                 ! e.g.,
                 ! scalelib/src/atmos-physics/microphysics/scale_atmos_phy_mp_tomita08.F90 
                 !                                 /radiation/scale_atmos_phy_rd_mstrnx.F90
                 ! "subroutine ATMOS_PHY_MP_tomita08_CloudFraction"
                 ! 
          profiles(iprof) % cfrac(ilev) = 0.5_jprb + sign(0.5_jprb, profiles(iprof) % cloud(2,ilev) + &
                                                                    profiles(iprof) % cloud(6,ilev) - Q_EPS)

        case (2) ! Tompkins and Janiskova (2004QJRMS) method (as in Okamoto 2017QJRMS)
                 !
          profiles(iprof) % cfrac(ilev) = cldfrac_TJ04(ptmp,tktmp,qvtmp) ! Pa, K, kgkg-1

        end select

        ! Need to modify? if openmp
        if(profiles(iprof) % cloud(2,ilev) + &
           profiles(iprof) % cloud(6,ilev) >= H08_RTTOV_MINQ_CTOP)then
          if(ctop_out(iprof) < 0.0d0)then
            ctop_out(iprof) = ptmp
          endif
        endif

      end do ! ilev
    endif ! addclouds

    if(debug .and. mod(iprof,20)==0)then
      do ilev = 1, nlevs + H08_RTTOV_KADD - 1
        write(6,'(a,i5,4f11.4)')"DEBUG PROF",ilev,profiles(iprof) % t(ilev),&
                                                  profiles(iprof) % q(ilev),&
                                                  profiles(iprof) % p(ilev),&
                                                  profiles(iprof) % cloud(6,ilev)
      end do ! ilev
      print *,""
    endif
  enddo ! prof
!### OMP END PARALLEL DO

!   do ilev = 1, H08_RTTOV_KADD + nlevs - 1
!     print *,"DEBUG RD RTTOV",ilev,profiles(iprof) % p(ilev)
!   enddo

  if(debug) write(6,*)"ch2",nprof,nlevs

  if(debug) WRITE(6,*) 'END SUBSTITUTE PROFILE'


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

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

  if (nthreads <= 1) then
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
            nthreads    = nthreads)    ! in    number of threads to use
  endif

  if (errorstatus /= errorstatus_success) then
    write (*,*) 'rttov_direct error'
    call rttov_exit(errorstatus)
  endif

  ! --- Output the results --------------------------------------------------

  do iprof = 1, nprof 

    joff = (iprof-1_jpim) * nchannels

    !
    !     OUTPUT RESULTS
    !
    btall_out(1:nchannels,iprof) = real(radiance%bt(1+joff:nchannels+joff), kind = r_size)
    btclr_out(1:nchannels,iprof) = real(radiance%bt_clear(1+joff:nchannels+joff), kind=r_size)

    if(debug .and. iprof<=2) print *,"DEBUG HIM8 SCALE_RTTOV:",btall_out(3,iprof),btclr_out(3,iprof),iprof
    if(debug .and. mod(iprof,100)==0) print *,"DEBUG HIM8 SCALE_RTTOV:",btall_out(3,iprof),btclr_out(3,iprof)

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
      do ilev = 2, nlevs - 1 + H08_RTTOV_KADD - 1
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

  ! Deallocate structures for rttov_direct
  call rttov_alloc_direct( &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        nprof,                   &
        nchanprof,               &
        nlevs+H08_RTTOV_KADD,    &
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
end subroutine SCALE_RTTOV_fwd12

function cldfrac_TJ04(pres,temp,qv)
  implicit none

  real(r_size), intent(in) :: pres ! (Pa)
  real(r_size), intent(in) :: temp ! (K)
  real(r_size), intent(in) :: qv   ! (kgkg-1)

  real(r_size) :: rh, rh_c
  real(r_size) :: sigma
  real(r_size) :: kappa

  real(r_size) :: cldfrac_TJ04

  !
  ! cloud fraction diagnosis based on 
  !  Tompkins and Janiskova (2004QJRMS):
  !  A cloud scheme for data assimilation: Description and initial tests
  !  

  rh = get_RH(pres,temp,qv)

  sigma = pres / pres00 ! non-dimensinoal level

  kappa = max(0.0d0, 0.9d0 * (sigma - 0.2) ** 0.2) ! eq. (6) in TJ04
  rh_c = 0.86d0 - 0.7d0 * sigma * (1.0d0 - sigma) * (1.85d0 + 0.95d0 * (sigma - 0.5d0)) ! eq. (7) in TJ04

  cldfrac_TJ04 = 1.0d0 - dsqrt((1.0d0 - rh)/(1.0d0 - rh_c - kappa * (rh - rh_c)))

  return
end function cldfrac_TJ04

function get_RH(pres,temp,qv)
  implicit none

  real(r_size), intent(in) :: pres ! (Pa)
  real(r_size), intent(in) :: temp ! (K)
  real(r_size), intent(in) :: qv   ! (kgkg-1)

  real(r_size) :: get_RH ! (out)

  real(r_size) :: es_tmp
  real(r_size) :: e_tmp

  ! saturation vapor pressure
  !  Tetens' formula
  es_tmp = 6.112d0 * dexp(17.67d0 * (temp - temp00)/(temp - temp00 + 243.5d0)) * 100.0d0 ! (Pa)

  ! vapor pressure
  e_tmp = qv * pres / (qv + Rdry / Rvap)

  get_RH = e_tmp / es_tmp

  return
end function get_RH

end module scale_H08_fwd12

