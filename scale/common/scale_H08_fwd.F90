module scale_H08_fwd
implicit none

contains

SUBROUTINE SCALE_RTTOV_fwd(nchannels,&
                           nlevs,&
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
                           tmp_ztop,&
                           btall_out,& 
                           btclr_out,& 
                           trans_out)
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
       & qmin ! H08

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
        H08_RTTOV_MINQ, &
        H08_RTTOV_CFRAC_CNST, &
        H08_RTTOV_EXTRA_US76 
! for Obs sim.
!  use mod_net2g_vars, ONLY: & !
!        H08_RTTOV_MINQ, &
!        H08_RTTOV_CFRAC_CNST, &
!        H08_RTTOV_CFRAC_CNST_TC99
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

!## Only for Obs simulation
!  INTEGER,PARAMETER :: r_size=kind(0.0d0)
!  INTEGER,PARAMETER :: RP = r_size
!  real(RP) :: Rdry    =  287.04_RP          !< specific gas constant (dry air) [J/kg/K]
!  real(RP), parameter :: Rvap    = 461.46_RP           !< specific gas constant (water vapor)    [J/kg/K]
!  real(RP), parameter :: CONST_PI      = 3.14159265358979_RP !< pi
!  real(RP), parameter :: Deg2Rad     = CONST_PI / 180.0_RP
!  real(RP), parameter :: CONST_GRAV    = 9.80665_RP  !< standard acceleration of gravity [m/s2]

!  LOGICAL,INTENT(IN) :: CLD
!  LOGICAL :: H08_RTTOV_EXTRA_US76 = .true.

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
  INTEGER, INTENT(IN) :: nlevs
  integer(kind=jpim):: icecld_ish=4, icecld_idg=0  !improved Baran, new in rttov11.2 recommended in  Rttov11
! ###

  Real(r_size),INTENT(IN) :: tmp_p(nlevs, nprof)
  Real(r_size),INTENT(IN) :: tmp_t(nlevs, nprof)
  Real(r_size),INTENT(IN) :: tmp_qv(nlevs, nprof)
  Real(r_size),INTENT(IN) :: tmp_qc(nlevs, nprof)
  Real(r_size),INTENT(IN) :: tmp_qice(nlevs, nprof)
  Real(r_size),INTENT(IN) :: tmp_t2m(nprof)
  Real(r_size),INTENT(IN) :: tmp_q2m(nprof)
  Real(r_size),INTENT(IN) :: tmp_p2m(nprof)
  Real(r_size),INTENT(IN) :: tmp_u2m(nprof)
  Real(r_size),INTENT(IN) :: tmp_v2m(nprof)
!  Real(r_size),INTENT(IN) :: tmp_soze(nprof)
!  Real(r_size),INTENT(IN) :: tmp_soaz(nprof)
!  Real(r_size),INTENT(IN) :: tmp_saze(nprof)
!  Real(r_size),INTENT(IN) :: tmp_saaz(nprof)
  Real(r_size),INTENT(IN) :: tmp_elev(nprof)
  Real(r_size),INTENT(IN) :: tmp_lon(nprof)
  Real(r_size),INTENT(IN) :: tmp_lat(nprof)
  Real(r_size),INTENT(IN) :: tmp_land(nprof)
  Real(r_size),INTENT(IN) :: tmp_ztop(nprof)


  Real(Kind=jprb),allocatable :: kgkg2gm3(:) ! convert parameter [kg/kg] => [gm^-3]
  Real(Kind=jprb),allocatable :: icec(:) ! ice cloud content (ice + snow + graupel)
  Real(Kind=jprb),allocatable :: liqc(:) ! liquid cloud content (cloud water)
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
  CHARACTER(LEN=256) :: coef_filename='./rtcoef_himawari_8_ahi.dat'
  CHARACTER(LEN=256) :: sccoef_filename='./sccldcoef_himawari_8_ahi.dat'
  CHARACTER(LEN=256) :: prof_filename
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim),intent(in) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim) :: ivch, ich
  REAL(KIND=jprb)    :: ems_val, brdf_val
  INTEGER(KIND=jpim) :: asw
  REAL(KIND=jprb), ALLOCATABLE :: emis(:), brdf(:)
  INTEGER(KIND=jpim), ALLOCATABLE :: nchan(:)
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch
  INTEGER(KIND=jpim) :: nch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

! by T.Honda
  Real(Kind=jprb),ALLOCATABLE :: tmp_btall_out(:,:)
  Real(Kind=jprb),ALLOCATABLE :: tmp_btclr_out(:,:)
  Real(Kind=jprb),ALLOCATABLE :: tmp_trans_out(:,:,:)
  REAL(Kind=r_size),INTENT(OUT) :: btall_out(nchannels,nprof)
  REAL(Kind=r_size),INTENT(OUT) :: btclr_out(nchannels,nprof)
  REAL(Kind=r_size),INTENT(OUT) :: trans_out(nlevs,nchannels,nprof)

  logical :: debug = .false.
!  logical :: debug = .true.

  real(kind=jprb) :: Rd 
  real(kind=jprb) :: Rv 
  real(kind=jprb) :: grav

  real(kind=jprb) :: epsb 
  real(kind=jprb) :: repsb 

! -- U.S. standard atmosphere (only below the stratopause at 47 km)--
  integer,parameter :: US76_midx = 4 ! Maximum idx for US76
                                     ! Extrapolation will be done only below Zm_US76(US76_midx)
  real(kind=jprb) :: Zm_US76(5) = (/0.0d0, 11000.0d0, 20000.0d0, 32000.0d0, 47000.0d0/) ! (m)
  real(kind=jprb) :: dTdZm_US76(4) = (/-6.5d-3, 0.0d-3, 1.0d-3, 2.8d-3/) ! (K/m)
  integer :: nlevs_add
  real(kind=jprb),allocatable :: Zm_add(:) ! (m)
  real(kind=jprb),allocatable :: TK_add(:) ! (K)
  real(kind=jprb),allocatable :: Phpa_add(:) ! (hPa)
  real(kind=jprb),parameter :: dZm_add = 1000.0_jprb ! (m)
  real(kind=jprb) :: dz, ztop
  integer :: US76_sidx, US76_idx
  integer :: slev, elev, k
 
! - set arrays for extrapolation using the U.S. 1976 atmos. --
  if(H08_RTTOV_EXTRA_US76)then
    ! maximum height 
    ztop = maxval(tmp_ztop)
    !--  search model top location
    do j = 1, 4
      if((ztop >= Zm_US76(j)) .and. (ztop < Zm_US76(j+1)))then
        US76_sidx = j
        exit
      endif
    enddo

    nlevs_add = 5 - US76_sidx
    !nlevs_add = int((Zm_US76(US76_midx) - tmp_ztop) / dZm_add + 0.99d0)
    write(6,'(a,I10)'),"US76, nlevs_add ",nlevs_add

    allocate(Zm_add(1:nlevs_add))
    allocate(TK_add(1:nlevs_add))
    allocate(Phpa_add(1:nlevs_add))

    !--  set Z (m) levels
    k = 1
    do j = US76_sidx, 5
      Zm_add(k) = Zm_US76(j)
      k = k + 1
    enddo

!    Zm_add(nlevs_add) = real(int(hgt_top / 1000.0D0 + 0.99d0)*1000,kind=jprb)
!    do ilev = nlevs_add - 1, 1, -1
!      Zm_add(ilev) = Zm_add(ilev+1) + dZm_add
!    enddo 
!
!    US76_sidx = 0
!    do ilev = 1, US76_midx - 1
!      if((hgt_top >= Zm_US76(ilev)) .and. (hgt_top < Zm_US76(ilev+1)))then
!        US76_sidx = ilev
!      endif
!    enddo
!
!    if(US76_sidx==0)then
!      write(6,'(a)')'Error from RTTOV!!' 
!      stop
!    endif

  else
    nlevs_add = 0
  endif ! H08_RTTOV_EXTRA_US76

! -- set thermodynamic constants
  Rd = real(Rdry,kind=jprb)
  Rv = real(Rvap,kind=jprb)
  epsb = Rd / Rv 
  repsb = 1.0_jprb / epsb

  grav = real(CONST_GRAV,kind=jprb)

  minQcfrac = real(H08_RTTOV_MINQ,kind=jprb)
  jcfrac_cnst = real(H08_RTTOV_CFRAC_CNST,kind=jprb)

  ALLOCATE(tmp_btall_out(nchannels,nprof))
  ALLOCATE(tmp_btclr_out(nchannels,nprof))
  ALLOCATE(tmp_trans_out(nlevs,nchannels,nprof))


  allocate(kgkg2gm3(nlevs))
  allocate(icec(nlevs))
  allocate(liqc(nlevs))

  if(debug) write(6,'(1x,a)')"hello from RTTOV"

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
  input_ems(:)=0.0
  input_brdf(:)=0.0

  if(debug) write(6,'(1x,a)')"hello from RTTOV2"

  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  opts % rt_ir % addsolar = .FALSE.         ! Do not include solar radiation
  opts % interpolation % addinterp  = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method

  opts % rt_all % addrefrac         = .FALSE.  ! Include refraction in path calc
  opts % rt_ir % addclouds          = .TRUE. ! Include cloud effects
  opts % rt_ir % user_cld_opt_param   = .FALSE. ! include cloud effects
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
!  opts%config%do_checkinput = .FALSE. ! see UG 85pp

  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  if(debug) write(6,'(1x,a)')"hello from RTTOV3"
  CALL rttov_read_coefs(errorstatus, coefs, opts, form_coef='formatted', file_coef=coef_filename, &
                       &file_sccld=sccoef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
!    nchannels = coefs % coef % fmv_chn
  ENDIF


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
      & nlevs+nlevs_add,         &
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
      & nlevs-1_jpim+nlevs_add, &
      & asw)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(6,*) 'allocation error for radiance arrays'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Allocate transmittance structure
  CALL rttov_alloc_transmission( &
      & errorstatus,             &
      & transmission,            &
      & nlevs-1_jpim+nlevs_add,  &
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
    WRITE(6,*) 'max p, max t',maxval(tmp_p),maxval(tmp_t)
    WRITE(6,*) 'min p, min t',minval(tmp_p),minval(tmp_t)
    WRITE(6,*) 'max p2, max t2',maxval(tmp_p2m),maxval(tmp_t2m)
    WRITE(6,*) 'min p2, min t2',minval(tmp_p2m),minval(tmp_t2m)

    WRITE(6,*) 'max qv, min qv',maxval(tmp_qv),minval(tmp_qv)
    WRITE(6,*) '-- t prof --'
    WRITE(6,*) 'size t:',size(tmp_t(:,1)),size(tmp_p(:,1))

    do j = 1, nlevs
    WRITE(6,*) 'qv:',tmp_qv(j,2)

    enddo
    WRITE(6,*) '-- end t prof --'
    do j = 1, nlevs
    WRITE(6,*) 't:',tmp_t(j,2)

    enddo
  endif

  !===============================================
  !========== READ profiles == start =============
  if(debug) WRITE(6,*) 'START SUBSTITUTE PROFILE'
  ! Note: Profiles are from top to surface.
  slev = nlevs+nlevs_add
  elev = nlevs_add+1
  DO iprof = 1, nprof
    profiles(iprof)%p(elev:slev)=real(tmp_p(:,iprof),kind=jprb) * 0.01_jprb  ! (hpa)
    profiles(iprof)%t(elev:slev)=real(tmp_t(:,iprof),kind=jprb)
    profiles(iprof)%q(elev:slev)=real(tmp_qv(:,iprof),kind=jprb) * q_mixratio_to_ppmv ! (ppmv)
    profiles(iprof)%s2m%t=real(tmp_t2m(iprof),kind=jprb)
    profiles(iprof)%s2m%q=real(tmp_q2m(iprof),kind=jprb) * q_mixratio_to_ppmv ! (ppmv)

    if(profiles(iprof)%s2m%q < qmin) profiles(iprof)%s2m%q = qmin + qmin * 0.01_jprb
    do ilev=elev,slev
      if(profiles(iprof)%q(ilev) < qmin) profiles(iprof)%q(ilev) = qmin + qmin * 0.01_jprb
    enddo

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

! sattelite angle 

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

!    profiles(iprof)% zenangle = 30.0_jprb ! tentative
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
      do ilev=1,nlevs
        tv = real(tmp_t(ilev,iprof),kind=jprb) * (1.0_jprb+real(tmp_qv(ilev,iprof),kind=jprb) * repsb) &
           / (1.0_jprb + real(tmp_qv(ilev,iprof),kind=jprb))
        kgkg2gm3(ilev) = real(tmp_p(ilev,iprof),kind=jprb) / (Rd * tv) * 1000.0_jprb 
        liqc(ilev) = real(max(tmp_qc(ilev,iprof),0.0_r_size),kind=jprb) * kgkg2gm3(ilev)
        icec(ilev) = real(max(tmp_qice(ilev,iprof),0.0_r_size),kind=jprb) * kgkg2gm3(ilev)
      end do !ilev

      do ilev=1,nlevs-1
        profiles(iprof) % cloud(2,ilev+nlevs_add) = & !stratus maritime (default)
                   (liqc(ilev+1) + liqc(ilev)) * 0.5_jprb
        profiles(iprof) % cloud(6,ilev+nlevs_add) = & 
                   (icec(ilev+1) + icec(ilev)) * 0.5_jprb
!
! cloud fraction diagnosis
        if(jcfrac_cnst <= 0.0_jprb)then
          if(((profiles(iprof) % cloud(2,ilev+nlevs_add) + &
               profiles(iprof) % cloud(6,ilev+nlevs_add)) > minQcfrac)) then
            profiles(iprof) % cfrac(ilev+nlevs_add) = 1.0_jprb  !cloud fraction 
          else
            profiles(iprof) % cfrac(ilev+nlevs_add) = 0.0_jprb  !cloud fraction 
          endif
        else ! default!!
          profiles(iprof) % cfrac(ilev+nlevs_add) = min((profiles(iprof) % cloud(2,ilev+nlevs_add) + &
                                               profiles(iprof) % cloud(6,ilev+nlevs_add)) / jcfrac_cnst, 1.0_jprb)
        endif
      end do ! ilev
    endif ! addclouds

! -- extrapolation --
    if(H08_RTTOV_EXTRA_US76)then
      do k = 1, nlevs_add
        if(k == 1)then
          dz = Zm_add(k) - tmp_ztop(iprof)
          TK_add(1) = real(tmp_t(elev,iprof),kind=jprb) + dTdZm_US76(US76_sidx) * dz
          Phpa_add(1) = real(tmp_p(elev,iprof),kind=jprb) * dexp(-grav/Rd/TK_add(1) * dz)
        else
          dz = Zm_add(k) - Zm_add(k-1)
          TK_add(k) = TK_add(k-1) + dTdZm_US76(US76_sidx+k-1) * dz
          Phpa_add(k) = Phpa_add(k-1) * dexp(-grav/Rd/TK_add(k) * dz)
        endif
       
      enddo 
      ! Caution !
      ! Additional arrays (e.g., Zm_add) are stored from the surface to TOA,
      ! but input profiles are stored from TOA to surface.
      ! 
      profiles(iprof)%q(1:nlevs_add) = qmin + qmin * 0.01_jprb 
      profiles(iprof)%p(1:nlevs_add) = Phpa_add(nlevs_add:1)
      profiles(iprof)%t(1:nlevs_add) = TK_add(nlevs_add:1)

!      dz = Zm_add(nlevs_add) - hgt_top
!      TK_add(nlevs_add) = real(tmp_t(1,iprof),kind=jprb) + dTdZm_US76(US76_sidx) * dz
!      Phpa_add = real(tmp_p(1,iprof),kind=jprb) * dexp(-grav/Rd/TK_add(nlevs_add)*dz)
!
!      US76_idx = US76_sidx
!      DO ilev = nlevs_add - 1, 1, -1
!        dz = Zm_add(ilev) - Zm_add(ilev+1)
!        if(Zm_add(ilev+1) > Zm_US76(US76_idx+1))then
!          US76_idx = max(US76_idx + 1, US76_midx - 1)
!        endif
!        TK_add(ilev) = TK_add(ilev+1) + dTdZm_US76(US76_idx) * dz
!        Phpa_add(ilev) = Phpa_add(ilev+1) * dexp(-grav/Rd/TK_add(ilev)*dz)
!
        !write(6,'(a,3F10.2)'),"US76, ",Zm_add(ilev),TK_add(ilev),Phpa_add(ilev)
!      ENDDO
  
!      profiles(iprof)%p(1:nlevs_add) = Phpa_add(1:nlevs_add)
!      profiles(iprof)%t(1:nlevs_add) = TK_add(1:nlevs_add)
!      profiles(iprof)%q(1:nlevs_add) = qmin + qmin * 0.01_jprb 

    endif ! H08_RTTOV_EXTRA_US76
  ENDDO ! prof

  deallocate(kgkg2gm3)

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
        & reflectance = reflectance)!,& ! inout input/output BRDFs per channel
!        & nthreads    = 8) ! tentative (10/14/2015) T.Honda
  IF (errorstatus /= errorstatus_success) THEN
    WRITE (6,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF

  if(debug) write(6,*)"Pass direct"

  ! --- Output the results --------------------------------------------------



  DO iprof = 1, nprof 

    joff = (iprof-1_jpim) * nchannels

    !
    !     OUTPUT RESULTS
    !
    nprint = 1 + INT((nchannels-1)/10)

    tmp_btall_out(1:nchannels,iprof)=radiance%bt(1+joff:nchannels+joff) 
    tmp_btclr_out(1:nchannels,iprof)=radiance%bt_clear(1+joff:nchannels+joff)

    DO ilev = 1, nlevs
      tmp_trans_out(ilev,1:nchannels,iprof) = transmission % tau_levels(ilev+nlevs_add,1+joff:nchannels+joff)
     if(debug .and. mod(iprof,5)==0) write(6,'(a,f10.7,i4)')"RTTOV debug trans",tmp_trans_out(ilev,1,iprof)
     
    ENDDO

!    DO ich = 1, nchannels
      ! Select transmittance based on channel type (VIS/NIR or IR)
!      IF (coefs % coef % ss_val_chn(chanprof(j+joff) % chan) == 2) THEN
!        DO ilev = 1, nlevs
!         tmp_trans_out(ilev,ich,iprof)=transmission % tau_levels_path1(ilev,joff+ich) 
!        ENDDO
!      ELSE
!        DO ilev = 1, nlevs
!         tmp_trans_out(ilev,ich,iprof)=transmission % tau_levels(ilev,joff+ich)
!        ENDDO
!      ENDIF
!    ENDDO

  ENDDO


  trans_out = REAL(tmp_trans_out,kind=r_size)
  btall_out = REAL(tmp_btall_out,kind=r_size)
  btclr_out = REAL(tmp_btclr_out,kind=r_size)

  DEALLOCATE(tmp_btall_out,tmp_btclr_out,tmp_trans_out)


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

  asw = 0 ! Switch for deallocation passed into RTTOV subroutines

  ! Deallocate radiance arrays
  CALL rttov_alloc_rad(errorstatus, nchannels, radiance, nlevs-1_jpim+nlevs_add, asw)
  IF(errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'radiance deallocation error'
  ENDIF

  ! Deallocate transmission arrays
  CALL rttov_alloc_transmission(errorstatus, transmission, nlevs-1_jpim+nlevs_add, nchannels, asw)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'transmission deallocation error'
  ENDIF

  ! Deallocate profile arrays
  CALL rttov_alloc_prof(errorstatus, nprof, profiles, nlevs+nlevs_add, opts, asw)
  IF (errorstatus /= errorstatus_success .or. alloc_status(1) /= 0) THEN
    WRITE(*,*) 'profile deallocation error'
  ENDIF
  DEALLOCATE(profiles, stat=alloc_status(1))
  IF (alloc_status(1) /= 0) THEN
    WRITE(*,*) 'mem deallocation error for profile array'
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


  if(debug) write(*,*)'Successfully finished!!'

return
END SUBROUTINE SCALE_RTTOV_fwd

end module scale_H08_fwd

