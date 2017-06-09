module scale_H08_fwd12
implicit none

contains

subroutine SCALE_RTTOV12_fwd(nlevels,&
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
                           !tmp_ztop,&
                             btall_out,& 
                             btclr_out,& 
                             trans_out,&
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
        H08_RTTOV_EXTRA_US76, &
        H08_RTTOV_MINQ_CTOP,  &
        nch

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


!
! -  Added by T.Honda (11/18/2015)
! -- Note: Computation of the zenith angle in each obs point (P) is based on the formula in
!          LRIT/HRIT Global Specification.
!          http://www.cgms-info.org/index_.php/cgms/page?cat=publications&page=technical+publications
! 
  REAL(r_size),PARAMETER :: Rpol = 6356.7523d3 ! a polar radius of Earth (m) 
  REAL(r_size) :: Rl ! a local radius of Earth
  REAL(r_size),PARAMETER :: sub_lon_H08 = 140.7d0 ! longitude of Himawari-8 satellite
  REAL(r_size) :: rlon,rlat ! (lon,lat) (Radian)
!
!
! Vector components for a satellite coordinate frame
!
!
  REAL(r_size) :: rnps, rnep, c_lat ! auxiliary variables
  REAL(r_size) :: r1, r2, r3       ! components of location vector for point P 
  REAL(r_size) :: r1ps, r2ps, r3ps ! components of the vector from P to the satellite 
  REAL(r_size) :: r1ep, r2ep, r3ep  ! components of the vector from the center of Earth to P
  REAL(r_size) :: z_angle_H08 ! zenith angle of Himawari-8

!
! minQcfrac: Threshold for diagnosing cloud fraction
  REAL(kind=jprb) :: jcfrac_cnst
  REAL(kind=jprb) :: minQcfrac ! threshold water/ice contents (g m-3) for cloud fraction diagnosis

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
  !Real(r_size),INTENT(IN) :: tmp_ztop(nprof)



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
  Real(Kind=jprb),ALLOCATABLE :: tmp_btall_out(:,:)
  Real(Kind=jprb),ALLOCATABLE :: tmp_btclr_out(:,:)
  Real(Kind=jprb),ALLOCATABLE :: tmp_trans_out(:,:,:)
  REAL(Kind=r_size),INTENT(OUT) :: btall_out(nch,nprof)
  REAL(Kind=r_size),INTENT(OUT) :: btclr_out(nch,nprof)
  REAL(Kind=r_size),INTENT(OUT) :: trans_out(nlevels,nch,nprof)
  REAL(Kind=r_size),INTENT(OUT) :: ctop_out(nprof)
  REAL(Kind=r_size) :: ptmp

  logical :: debug = .false.
  !logical :: debug = .true.
  logical :: in_warning = .false.

  real(kind=jprb) :: Rd 
  real(kind=jprb) :: Rv 
  real(kind=jprb) :: grav

  real(kind=jprb) :: epsb 
  real(kind=jprb) :: repsb 

  integer :: slev, elev

  if(debug) write(6,'(1x,a)')"hello from RTTOV"


! -- set thermodynamic constants
  Rd = real(Rdry,kind=jprb)
  Rv = real(Rvap,kind=jprb)
  epsb = Rd / Rv 
  repsb = 1.0_jprb / epsb

  grav = real(CONST_GRAV,kind=jprb)

  jcfrac_cnst = real(H08_RTTOV_CFRAC_CNST,kind=jprb)

  ALLOCATE(tmp_btall_out(nch,nprof))
  ALLOCATE(tmp_btclr_out(nch,nprof))
  ALLOCATE(tmp_trans_out(nlevels,nch,nprof))



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

  if(debug) write(6,'(1x,a)')"hello from RTTOV2"


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
                        file_coef=coef_filename, file_sccld=cld_coef_filename, &
                        form_coef='unformatted', form_sccld='unformatted')
! default               file_coef=coef_filename, file_sccld=cld_coef_filename)

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
  slev = nlevels
  elev = 1

  ! Loop over all profiles
  do iprof = 1, nprof

    ! Pressure (hPa), temp (K), WV (kg/kg)
    profiles(iprof)%p(elev:slev)=real(tmp_p(1:nlevels,iprof) * 0.01,kind=jprb) ! (hPa)
    profiles(iprof)%t(elev:slev)=real(tmp_t(1:nlevels,iprof),kind=jprb)
    profiles(iprof)%q(elev:slev)=real(tmp_qv(1:nlevels,iprof),kind=jprb) ! (kg kg-1)

    ! Check T & Q inputs
    do ilev=elev,slev
      if(profiles(iprof)%t(ilev) < tmin)then
        if(.not.in_warning)then
          write(6,*)'!! WARNING !! T input for RTTOV has unphysical values!!'
          write(6,'(a,3f10.3)')'!! WARNING!!',profiles(iprof)%q(ilev),&
                                              profiles(iprof)%p(ilev),&
                                              profiles(iprof)%t(ilev)
          write(6,'(a,3i8)')'!! WARNING!!',ilev,elev,slev
          in_warning = .true.
        endif

        profiles(iprof)%t(ilev) = tmin + tmin * 0.01_jprb
      endif

      if(profiles(iprof)%q(ilev) < qmin)then
        if(.not.in_warning)then
          write(6,*)'!! WARNING !! Q input for RTTOV is too small!'
          write(6,'(a,3f10.3)')'!! WARNING!!',profiles(iprof)%q(ilev),&
                                              profiles(iprof)%p(ilev),&
                                              profiles(iprof)%t(ilev)
          in_warning = .true.
        endif

        profiles(iprof)%q(ilev) = qmin + qmin * 0.01_jprb
      endif
    enddo ! End of check T & Q inputs

    ! 2 meter air variables
    profiles(iprof)%s2m%p=real(tmp_p2m(iprof) * 0.01,kind=jprb) ! (hPa)
    profiles(iprof)%s2m%t=real(tmp_t2m(iprof),kind=jprb) ! (K)
    profiles(iprof)%s2m%q=real(tmp_q2m(iprof),kind=jprb) ! (kg/kg)
    if(profiles(iprof)%s2m%t < tmin) profiles(iprof)%s2m%t = tmin + tmin * 0.01_jprb
    if(profiles(iprof)%s2m%q < qmin) profiles(iprof)%s2m%q = qmin + qmin * 0.01_jprb

    profiles(iprof)%s2m%u=real(tmp_u2m(iprof),kind=jprb)
    profiles(iprof)%s2m%v=real(tmp_v2m(iprof),kind=jprb)

   ! Skin variables
    profiles(iprof) % skin % t = real(tmp_t2m(iprof),kind=jprb)
    profiles(iprof) % skin % surftype = int(tmp_land(iprof))
    profiles(iprof) % skin % watertype = 1 ! tentative (11/18/2015)
    profiles(iprof)%s2m%wfetc= 100000.0_jprb

    ! Elevation, latitude and longitude
    profiles(iprof) % elevation = real(tmp_elev(iprof)*0.001,kind=jprb) ! (km)
    profiles(iprof) % latitude  = real(tmp_lat(iprof),kind=jprb)
    profiles(iprof) % longitude = real(tmp_lon(iprof),kind=jprb)

    ! Satellite and solar angles

    rlat = tmp_lat(iprof)*Deg2Rad
    rlon = tmp_lon(iprof)*Deg2Rad

    c_lat = datan(0.993305616d0 * dtan(rlat))
    Rl = Rpol / dsqrt(1.0d0 - 0.00669438444d0 * dcos(c_lat)*dcos(c_lat))
    r1 = 42164.0d3 - Rl * dcos(c_lat) * dcos(rlon - sub_lon_H08*Deg2Rad)
    r2 = -Rl * dcos(c_lat) * dsin(rlon - sub_lon_H08*Deg2Rad)
    r3 = Rl * dsin(c_lat)
    rnps = dsqrt(r1*r1+r2*r2+r3*r3)

    r1ps = r1 * (-1.0d0)
    r2ps = r2 * (-1.0d0)
    r3ps = r3 * (-1.0d0)

    r1ep = r1 - 42164.0d3
    r2ep = r2 
    r3ep = r3
 
    rnep = dsqrt(r1ep*r1ep+r2ep*r2ep+r3ep*r3ep)
    z_angle_H08 = r1ps * r1ep + r2ps * r2ep + r3ps * r3ep ! internal product 
    z_angle_H08 = dacos(z_angle_H08/(rnps*rnep))/Deg2Rad  

    profiles(iprof)% zenangle = real(z_angle_H08,kind=jprb) ! (11/18/2015)
    if(mod(iprof,1000)==0.and.debug)write(6,'(a,f15.10)'),'zenangle ',profiles(iprof)% zenangle
    if(mod(iprof,1000)==0.and.debug)write(6,'(a,4f10.5)'),' ',rlon,rlat,tmp_lon(iprof),tmp_lat(iprof)

!    profiles(iprof)%azangle=0.0_jprb     ! Not required for [opts % rt_ir % addsolar = .FALSE.] 
!    profiles(iprof)%sunzenangle=0.0_jprb ! Not required for [opts % rt_ir % addsolar = .FALSE.] 
!    profiles(iprof)%sunazangle=0.0_jprb  ! Not required for [opts % rt_ir % addsolar = .FALSE.]


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


  DO iprof = 1, nprof 

    joff = (iprof-1_jpim) * nch

    !
    !     OUTPUT RESULTS
    !
    nprint = 1 + INT((nch-1)/10)

    tmp_btall_out(1:nch,iprof)=radiance%bt(1+joff:nch+joff) 
    tmp_btclr_out(1:nch,iprof)=radiance%bt_clear(1+joff:nch+joff)

    DO ilev = 1, nlevels
      tmp_trans_out(ilev,1:nch,iprof) = transmission % tau_levels(ilev,1+joff:nch+joff)
     if(debug .and. mod(iprof,50)==0) write(6,'(a,f10.7,i4)')"RTTOV debug trans",tmp_trans_out(ilev,1,iprof)
    ENDDO


  ENDDO

  trans_out = REAL(tmp_trans_out,kind=r_size)
  btall_out = REAL(tmp_btall_out,kind=r_size)
  btclr_out = REAL(tmp_btclr_out,kind=r_size)

  DEALLOCATE(tmp_btall_out,tmp_btclr_out,tmp_trans_out)

  ! --- End of output section -----------------------------------------------


  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  ! Deallocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        nprof,                   &
        nchanprof,               &
        nlevels,                   &
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
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF

  if(debug) write(*,*)'Successfully finished!!'

return
end subroutine SCALE_RTTOV12_fwd

end module scale_H08_fwd12

