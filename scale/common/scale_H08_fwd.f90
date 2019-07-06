module scale_H08_fwd
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

subroutine SCALE_RTTOV_fwd(nchannels,&
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
!                     & tmp_soze, tmp_soaz, tmp_saze, tmp_saaz,       &
                           elev,&
                           lon,&
                           lat,&
                           land,&
                           zenith,&
!                           kidx_rlx, &
                           RD_presh, &
                           RD_temph, &
                           btall_out,& 
                           btclr_out,& 
                           mwgt_plev,&
                           ctop_out)

  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2010, EUMETSAT, All Rights Reserved.
  !
  !     *************************************************************
  !
  !     TEST PROGRAM FOR RTTOV FORWARD MODEL
  !          RTTOV VERSION 11
  !
  ! To run this program you must have the following files
  ! either resident in the same directory or set up as a
  ! symbolic link:
  !   the file containing input profiles (e.g. prof.dat)
  !   the RTTOV coefficient file
  !
  ! The script run_example_fwd.sh may be used to run this program.
  ! The output is generated in a file called example_fwd_output.dat.
  !
  !
  ! To adapt the code to their own applications users should
  ! edit the code highlighted like this:
  !     !================================
  !     !======Read =====start===========
  !          code to be modified
  !     !======Read ===== end ===========
  !     !================================
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date        Comment
  ! -------   ----        -------
  !  1.0    27/04/2004   orginal (based on tstrad) P. Brunel
  !  1.1    09/08/2004   modified to allow for variable no. channels/per profile
  !                       R. Saunders
  !  1.2    13/04/2007   Modified for RTTOV-90
  !  1.3    31/07/2007   Modified for RTTOV-91 R Saunders
  !  1.4    11/10/2007   Parallel version P.Marguinaud
  !  2.0    25/06/2010   Modified for RTTOV-10 J Hocking
  !  2.1    23/11/2011   Updates for v10.2 J Hocking
  !  3.0    23/11/2012   Updated for v11 J Hocking
  !
  ! Code Description:
  !   Language:          Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & platform_name,       &
       & inst_name,           &
       & q_mixratio_to_ppmv,  &
       & qmin,                &
       & tmin

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
       & rttov_options,       &
       & rttov_coefs,         &
       & profile_type,        &
       & transmission_type,   &
       & radiance_type,       &
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
        H08_RTTOV_PROF_SHIFT, &
        H08_RTTOV_CFRAC,      &
        H08_RTTOV_KADD
  IMPLICIT NONE

#include "rttov_parallel_direct.interface"
#include "rttov_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  INTEGER, INTENT(IN) :: nprof
  INTEGER, INTENT(IN) :: nlevs
  integer(kind=jpim):: icecld_ish=4, icecld_idg=0  !improved Baran, new in rttov11.2 recommended in  Rttov11
! ###

  Real(r_size),INTENT(IN) :: prs(nlevs, nprof)
  Real(r_size),INTENT(IN) :: tk(nlevs, nprof)
  Real(r_size),INTENT(IN) :: qv(nlevs, nprof)
  Real(r_size),INTENT(IN) :: qc(nlevs, nprof)
  Real(r_size),INTENT(IN) :: qice(nlevs, nprof)
  Real(r_size),INTENT(IN) :: tk2m(nprof)
  Real(r_size),INTENT(IN) :: q2m(nprof)
  Real(r_size),INTENT(IN) :: prs2m(nprof)
  Real(r_size),INTENT(IN) :: u2m(nprof)
  Real(r_size),INTENT(IN) :: v2m(nprof)
!  Real(r_size),INTENT(IN) :: tmp_soze(nprof)
!  Real(r_size),INTENT(IN) :: tmp_soaz(nprof)
!  Real(r_size),INTENT(IN) :: tmp_saze(nprof)
!  Real(r_size),INTENT(IN) :: tmp_saaz(nprof)
  Real(r_size),INTENT(IN) :: elev(nprof)
  Real(r_size),INTENT(IN) :: lon(nprof)
  Real(r_size),INTENT(IN) :: lat(nprof)
  Real(r_size),INTENT(IN) :: land(nprof)
  Real(r_size),INTENT(IN) :: zenith(nprof)


  Real(Kind=jprb) :: kgkg2gm3 ! convert parameter [kg/kg] => [gm^-3]
  Real(Kind=jprb) :: icec1, icec2 ! ice cloud content (ice + snow + graupel)
  Real(Kind=jprb) :: liqc1, liqc2 ! liquid cloud content (cloud water)
  Real (Kind=jprb):: tv ! virtual temp. (K)

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output
  INTEGER(KIND=jpim), PARAMETER :: mxchn = 9000 ! max number of channels

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)                  :: opts           ! Options structure
  TYPE(rttov_coefs)                    :: coefs          ! Coefficients structure
  TYPE(rttov_chanprof),    ALLOCATABLE :: chanprof(:)    ! Input channel/profile list
  TYPE(profile_type),      ALLOCATABLE :: profiles(:)    ! Input profiles
  LOGICAL(KIND=jplm),      ALLOCATABLE :: calcemis(:)    ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  ALLOCATABLE :: emissivity(:)  ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      ALLOCATABLE :: calcrefl(:)    ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), ALLOCATABLE :: reflectance(:) ! Input/output surface BRDF
  TYPE(transmission_type)              :: transmission   ! Output transmittances
  TYPE(radiance_type)                  :: radiance       ! Output radiances

  INTEGER(KIND=jpim)                   :: errorstatus    ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status(10)
  CHARACTER(LEN=11)  :: NameOfRoutine = 'example_fwd'

  ! variables for input
  !====================
  INTEGER(KIND=jpim) :: input_chan(mxchn)
  REAL(KIND=jprb)    :: input_ems(mxchn), input_brdf(mxchn)
  CHARACTER(LEN=256) :: coef_filename='/rtcoef_himawari_8_ahi.dat'
  CHARACTER(LEN=256) :: sccoef_filename='/sccldcoef_himawari_8_ahi.dat'
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim),intent(in) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim) :: ich
  !REAL(KIND=jprb)    :: ems_val, brdf_val
  INTEGER(KIND=jpim) :: asw
  REAL(KIND=jprb), ALLOCATABLE :: emis(:), brdf(:)
  INTEGER(KIND=jpim), ALLOCATABLE :: nchan(:)
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch
  INTEGER(KIND=jpim) :: nch
  INTEGER(KIND=jpim) :: ilev
  INTEGER(KIND=jpim) :: iprof, joff

! by T.Honda
  REAL(Kind=r_size),INTENT(OUT) :: btall_out(nchannels,nprof)
  REAL(Kind=r_size),INTENT(OUT) :: btclr_out(nchannels,nprof)
  REAL(Kind=r_size),INTENT(OUT) :: ctop_out(nprof)
  REAL(Kind=r_size) :: ptmp, tktmp, qvtmp

  REAL(Kind=r_size) :: rdp, max_wgt, tmp_wgt
  REAL(Kind=r_size),INTENT(OUT) :: mwgt_plev(nchannels,nprof) ! Max weight level (Pa)

  logical :: debug = .false.
!  logical :: debug = .true.

  real(kind=jprb) :: grav

  real(kind=jprb) :: repsb 

  integer :: orgk

!  integer(kind=jpim),intent(in) :: kidx_rlx
  real(kind=r_size),intent(in)  :: RD_presh(nlevs+H08_RTTOV_KADD+1)
  real(kind=r_size),intent(in)  :: RD_temph(nlevs+H08_RTTOV_KADD+1)
  real(kind=jprb) :: tmp_dif
  
  if(debug) write(6,'(1x,a)')"hello from RTTOV"


! -- set thermodynamic constants

  grav = real(CONST_GRAV,kind=jprb)
  repsb = real(Rvap / Rdry, kind=jprb)

  !
  ! The usual steps to take when running RTTOV are as follows:
  !   1. Specify required RTTOV options
  !   2. Read coefficients
  !   3. Set up the chanprof array with the channels/profiles to simulate
  !   4. Allocate RTTOV profile, radiance and transmittance structures
  !   5. Read input profile(s)
  !   6. Set up surface emissivity and/or reflectance
  !   7. ll rttov_direct and store results
  !   8. Deallocate all structures and arrays

  errorstatus     = 0_jpim
  alloc_status(:) = 0_jpim

  dosolar=0


  do ich = 1, nchannels
    input_chan(ich)=ich
  end do
  input_ems(:) = 0.0_jprb
  input_brdf(:) = 0.0_jprb

  if(debug) write(6,'(1x,a)')"hello from RTTOV2"

  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  opts % rt_ir % addsolar = .FALSE.         ! Do not include solar radiation
  opts % interpolation % addinterp  = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method

  opts % rt_all % addrefrac         = .FALSE.  ! Include refraction in path calc
  opts % rt_ir % addclouds          = .TRUE. ! Include cloud effects
  opts % rt_ir % user_cld_opt_param   = .FALSE. 
  opts % rt_ir % addaerosl          = .FALSE. ! Don't include aerosol effects

  opts % rt_ir % ozone_data         = .FALSE.  ! We have an ozone profile
  opts % rt_ir % co2_data           = .FALSE. ! We do not have profiles
  opts % rt_ir % n2o_data           = .FALSE. !   for any other constituents
  opts % rt_ir % ch4_data           = .FALSE. !
  opts % rt_ir % co_data            = .FALSE. !
  opts % rt_mw % clw_data           = .FALSE. !

! added by T.Honda(2015/06/10)
  opts % interpolation % reg_limit_extrap  = .TRUE.  ! see UG 7.3 (32pp)
  opts % config % apply_reg_limits  = .TRUE.  ! see UG 7.3 (32pp)
  if(debug)then
    opts%config%do_checkinput = .FALSE. ! see UG 85pp 
  else
    opts%config%do_checkinput = .TRUE. ! see UG 85pp 
  endif
  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  if(debug) write(6,'(1x,a)')"hello from RTTOV3"
  CALL rttov_read_coefs(errorstatus, coefs, opts, form_coef='formatted', &
                       &file_coef=trim(H08_RTTOV_COEF_PATH)//trim(coef_filename), &
                       &file_sccld=trim(H08_RTTOV_COEF_PATH)//trim(sccoef_filename))
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
!  IF (nchannels > coefs % coef % fmv_chn) THEN
!    nchannels = coefs % coef % fmv_chn
!  ENDIF


  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! --------------------------------------------------------------------------
  ! 3. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------
  if(debug) write(6,'(1x,a)')"hello from RTTOV4"

  ! In general the number of channels simulated per profile can vary, but for
  ! this example we simulate all specified channels for each profile.

  ALLOCATE(nchan(nprof))     ! Number of channels per profile
  nchan(:) = nchannels
  nchanprof = SUM(nchan(:))  ! Size of chanprof array is total number of channels over all profiles

  ! Pack channels and input emissivity arrays
  ALLOCATE(chanprof(nchanprof))
  ALLOCATE(emis(nchanprof))
  ALLOCATE(brdf(nchanprof))

  nch = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nchan(j)
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = input_chan(jch)
      emis(nch)          = input_ems(jch)
      brdf(nch)          = input_brdf(jch)
    ENDDO
  ENDDO

  if(debug) write(6,'(1x,a)')"hello from RTTOV5"

  ! --------------------------------------------------------------------------
  ! 4. Allocate profiles, radiance and transmittance structures
  ! --------------------------------------------------------------------------

  if(debug) write(6,*)"np",nprof

  asw = 1 ! Switch for allocation passed into RTTOV subroutines

  ! Allocate input profile arrays
  ALLOCATE(profiles(nprof), stat=alloc_status(1))
  IF (ANY(alloc_status /= 0)) THEN
    WRITE(6,*) 'mem allocation error for profile array'
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_alloc_prof(         &
      & errorstatus,             &
      & nprof,                   &
      & profiles,                &
      & nlevs+H08_RTTOV_KADD,    &
      & opts,                    &
      & asw,                     &
      & coefs=coefs,             &
      & init=.TRUE._jplm         )
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(6,*) 'allocation error for profile arrays'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Allocate output radiance arrays
  CALL rttov_alloc_rad( &
      & errorstatus,    &
      & nchanprof,      &
      & radiance,       &
      & nlevs+H08_RTTOV_KADD-1_jpim, &
      & asw)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(6,*) 'allocation error for radiance arrays'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Allocate transmittance structure
  CALL rttov_alloc_transmission( &
      & errorstatus,             &
      & transmission,            &
      & nlevs+H08_RTTOV_KADD-1_jpim,  &
      & nchanprof,               &
      & asw,                     &
      & init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(6,*) 'allocation error for transmission arrays'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Allocate arrays for surface emissivity
  ALLOCATE(calcemis(nchanprof), stat=alloc_status(1))
  ALLOCATE(emissivity(nchanprof), stat=alloc_status(2))
  IF (ANY(alloc_status /= 0)) THEN
    WRITE(6,*) 'mem allocation error for emissivity arrays'
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  ! Allocate arrays for surface reflectance
  ALLOCATE(calcrefl(nchanprof), stat=alloc_status(1))
  ALLOCATE(reflectance(nchanprof), stat=alloc_status(2))
  IF (ANY(alloc_status /= 0)) THEN
    WRITE(6,*) 'memallocation error for reflectance arrays'
    CALL rttov_exit(errorstatus_fatal)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------
  if(debug) then
    WRITE(6,*) 'debug information'
    WRITE(6,*) 'max p, max t',maxval(prs),maxval(tk)
    WRITE(6,*) 'min p, min t',minval(prs),minval(tk)
    WRITE(6,*) 'max p2, max t2',maxval(prs2m),maxval(tk2m)
    WRITE(6,*) 'min p2, min t2',minval(prs2m),minval(tk2m)

    WRITE(6,*) 'max qv, min qv',maxval(qv),minval(qv)
    WRITE(6,*) '-- t prof --'
    WRITE(6,*) 'size t:',size(tk(:,1)),size(prs(:,1))

  endif

  !===============================================
  !========== READ profiles == start =============
  if(debug) WRITE(6,*) 'START SUBSTITUTE PROFILE'

  ! Note: Profiles are from top to surface.

!### $OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iprof,ilev,rat,rdk,orgk,ptmp,tktmp,qvtmp,tv,kgkg2gm3,icec1,icec2,liqc1,liqc2,dz,tmp_dif)
  DO iprof = 1, nprof ! iprof

    if(H08_RTTOV_PROF_SHIFT)then
      tmp_dif = RD_temph(H08_RTTOV_KADD+1) - tk(1,iprof)
    else
      tmp_dif = 0.0_jprb
    endif ! H08_RTTOV_PROF_SHIFT
   
    ! Above the model top (p < RD_presh(H08_RTTOV_KADD+1))
    do ilev = 1, H08_RTTOV_KADD
      profiles(iprof)%p(ilev) = real(RD_presh(ilev),kind=jprb) ! (hPa)
      profiles(iprof)%t(ilev) = max(RD_temph(ilev) - tmp_dif,tmin * 1.01_jprb) ! (K)

      profiles(iprof)%q(ilev) = qmin * 1.01_jprb

    enddo

    do ilev = H08_RTTOV_KADD + 1, H08_RTTOV_KADD + nlevs

      profiles(iprof)%p(ilev) = real(prs(ilev-H08_RTTOV_KADD,iprof),kind=jprb) * 0.01_jprb ! (hPa)
      profiles(iprof)%t(ilev) = real(tk(ilev-H08_RTTOV_KADD,iprof),kind=jprb) ! (K) 

      profiles(iprof)%q(ilev) = max(qv(ilev-H08_RTTOV_KADD,iprof) * q_mixratio_to_ppmv, qmin * 1.01_jprb) ! (ppmv)

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

!
!   fastem is used for MW
!    profiles(iprof) % skin % fastem(1) = 3.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(2) = 5.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(3) =15.0 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(4) = 0.1 ! comment out (11/18/2015)
!    profiles(iprof) % skin % fastem(5) = 0.3 ! comment out (11/18/2015)
!

    profiles(iprof) % skin % surftype = int(land(iprof))
    profiles(iprof) % skin % watertype = 1 ! tentative (11/18/2015)

    profiles(iprof) % elevation = real(elev(iprof),kind=jprb) * 0.001_jprb ! (km)
    profiles(iprof) % latitude  = real(lat(iprof),kind=jprb)
    profiles(iprof) % longitude = real(lon(iprof),kind=jprb)


    profiles(iprof)% zenangle = real(zenith(iprof),kind=jprb) ! (11/18/2015)
    if(mod(iprof,1000) == 0 .and. debug)write(6,'(a,f15.10)'),'zenangle ',profiles(iprof)% zenangle
    if(mod(iprof,1000) == 0 .and. debug)write(6,'(a,2f10.5)'),' ',lon(iprof),lat(iprof)

!    profiles(iprof)%azangle=0.0_jprb     ! Not required for [opts % rt_ir % addsolar = .FALSE.] 
!    profiles(iprof)%sunzenangle=0.0_jprb ! Not required for [opts % rt_ir % addsolar = .FALSE.] 
!    profiles(iprof)%sunazangle=0.0_jprb  ! Not required for [opts % rt_ir % addsolar = .FALSE.]

! These are parameters for simple cloud.
! Not used.
    profiles(iprof) % ctp       = 500.0_jprb
    profiles(iprof) % cfraction = 0.0_jprb
 

  !-- 6 general cloud 
    if( opts % rt_ir % addclouds ) then
      ! Cloud variables for cloud scattering scheme 

      profiles(iprof) % ish = icecld_ish  !ice water shape
      profiles(iprof) % idg = icecld_idg  !ice water effective diameter 
      profiles(iprof) % icede(:)= 0._jprb !ice effective diameter, set non-zero if you give by yourself

      profiles(iprof) % cloud(:,:) = 0._jprb
      profiles(iprof) % cfrac(:)   = 0._jprb

!  --- 6.1 microphysical clouds 
!    if(icldprf==1) then
    ! --convert kg/kg into g/m3, make liq.cloud and ice.cloud amount, and upside-down
 
      ctop_out(iprof) = -1.0d0 

      do ilev = H08_RTTOV_KADD + 1, H08_RTTOV_KADD + nlevs - 1
        orgk = ilev - H08_RTTOV_KADD ! k index for original profile

        ! ilev
        tv = real(tk(orgk,iprof) * (1.0d0+qv(orgk,iprof)) * repsb &
                   / (1.0d0 + qv(orgk,iprof)), kind=jprb)

        kgkg2gm3 = real(prs(orgk,iprof),kind=jprb) / (Rdry * tv) * 1000.0_jprb 

        liqc1 = real(max(qc(orgk,iprof),0.0_r_size),kind=jprb) * kgkg2gm3
        icec1 = real(max(qice(orgk,iprof),0.0_r_size),kind=jprb) * kgkg2gm3

        ! ilev + 1
        tv = real(tk(orgk+1,iprof) * (1.0d0+qv(orgk+1,iprof)) * repsb &
                   / (1.0d0 + qv(orgk+1,iprof)), kind=jprb)

        kgkg2gm3 = real(prs(orgk+1,iprof),kind=jprb) / (Rdry * tv) * 1000.0_jprb 

        liqc2 = real(max(qc(orgk+1,iprof),0.0_r_size),kind=jprb) * kgkg2gm3
        icec2 = real(max(qice(orgk+1,iprof),0.0_r_size),kind=jprb) * kgkg2gm3

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
                 ! e.g., scalelib/src/atmos-physics/microphysics/scale_atmos_phy_mp_tomita08.F90 
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

  ! Copy emissivities into RTTOV input emissivity array
  emissivity(:) % emis_in = emis(:)

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less
  calcemis(:) = emis(:) <= 0._jprb

  ! Copy BRDFs into RTTOV input reflectance array
  reflectance(:) % refl_in = brdf(:)

  ! Calculate BRDF within RTTOV where the input BRDF value is zero or less
  calcrefl(:) = brdf(:) <= 0._jprb

  ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
  reflectance(:) % refl_cloud_top = 0._jprb


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------

  if(debug) write(6,*)"Enter direct"

  CALL rttov_parallel_direct(                &
!  CALL rttov_direct(                &
        & errorstatus,              &! out   error flag
        & chanprof,                 &! in    channel and profile index structure
        & opts,                     &! in    options structure
        & profiles,                 &! in    profile array
        & coefs,                    &! in    coefficients strucutre
        & transmission,             &! inout computed transmittances
        & radiance,                 &! inout computed radiances
        & calcemis    = calcemis,   &! in    flag for internal emissivity calcs
        & emissivity  = emissivity, &! inout input/output emissivities per channel
        & calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
        & reflectance = reflectance ,& ! inout input/output BRDFs per channel
        & nthreads    = 8) ! Assume K computer
  IF (errorstatus /= errorstatus_success) THEN
    WRITE (6,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF

  if(debug) write(6,*)"Pass direct"

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
  DEALLOCATE (emis,        stat=alloc_status(1))
  DEALLOCATE (brdf,        stat=alloc_status(2))
  DEALLOCATE (nchan,       stat=alloc_status(3))
  DEALLOCATE (chanprof,    stat=alloc_status(4))
  DEALLOCATE (emissivity,  stat=alloc_status(5))
  DEALLOCATE (calcemis,    stat=alloc_status(6))
  DEALLOCATE (reflectance, stat=alloc_status(7))
  DEALLOCATE (calcrefl,    stat=alloc_status(8))


  IF (ANY(alloc_status /= 0)) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  if(debug) write(6,'(a)') 'CHK100'

  asw = 0 ! Switch for deallocation passed into RTTOV subroutines

  ! Deallocate radiance arrays
  CALL rttov_alloc_rad(errorstatus, nchannels, radiance, nlevs-1_jpim+H08_RTTOV_KADD, asw)
  IF(errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'radiance deallocation error'
  ENDIF

  if(debug) write(6,'(a)') 'CHK101'
  ! Deallocate transmission arrays
  CALL rttov_alloc_transmission(errorstatus, transmission, nlevs-1_jpim+H08_RTTOV_KADD, nchannels, asw)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'transmission deallocation error'
  ENDIF

  if(debug) write(6,'(a)') 'CHK102'
  ! Deallocate profile arrays
  CALL rttov_alloc_prof(errorstatus, nprof, profiles, nlevs+H08_RTTOV_KADD, opts, asw)
  if(debug) write(6,'(a)') 'CHK102.1'
  IF (errorstatus /= errorstatus_success .or. alloc_status(1) /= 0) THEN
    WRITE(*,*) 'profile deallocation error'
  ENDIF
  DEALLOCATE(profiles, stat=alloc_status(1))
  if(debug) write(6,'(a)') 'CHK102.2'
  IF (alloc_status(1) /= 0) THEN
    WRITE(*,*) 'mem deallocation error for profile array'
  ENDIF

  if(debug) write(6,'(a)') 'CHK103'
  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


  if(debug) write(*,*)'Successfully finished!!'

return
end subroutine SCALE_RTTOV_fwd

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

end module scale_H08_fwd

