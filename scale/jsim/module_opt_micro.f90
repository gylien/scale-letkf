  module module_opt_micro
    use common, only: r_size
    use common_geosatpr
    implicit none

    public :: opt_micro_radar, & ! compute microwave single-scattering properties
               watoptic          ! compute water dielectric constant for a given microwave frequency

    real(r_size),allocatable,dimension(:,:,:,:),public :: &
        kexttot, & ! total extinction coef
        salbtot, & ! total single scattering albedo
        asymtot, & ! total assymetry parameter
        sbacktot   ! total back scattering

    real(r_size),parameter :: q_min_micro = q_min_condensate ! min. detectable condensate amound [g/cm3]
    real(r_size),parameter :: min_bext = 1.d-10   ! minimum extinction for stability

    contains

  !>---------------------------------------------------------------------
  !> Masunaga, H., and C.D. Kummerow, 2005: Combined Radar and Radiometer Analysis of 
  !>     Precipitation Profiles for a Parametric Retrieval Algorithm. J. Atmos. Oceanic 
  !>     Technol., 22, 909-929.
  subroutine opt_micro_radar
    implicit none

    integer :: mxfreq   ! maximum num. of frequency
    real(r_size) :: atm_ext   ! atmos. gas extinction
    real(r_size) :: freqcy    ! frequency [GHz]
    real(r_size) :: tavg      ! temperature [K]
    real(r_size) :: tavg_c    ! temperature [K]
    real(r_size) :: pres      ! pressure
    real(r_size) :: melt_frac_bulk    ! melting fraction

    type( particle_gce ) :: &
        kext_gce, & ! extinction ceoff.
        salb_gce, & ! single scattering albedo
        asym_gce, & ! asymetry param.
        pbck_gce    ! backscattering ceoff.

    integer :: i, j ,k, nf


!    mxfreq = 1
!    nf = mxfreq
!    freqcy = freq_radar
!
!    ! allocate
!    if( allocated(kexttot) ) deallocate( kexttot,salbtot,asymtot,sbacktot ) !for total
!    if( .not. allocated(kexttot) )then
!      allocate( kexttot (is:ie,js:je,ks:ke,mxfreq), &
!                salbtot (is:ie,js:je,ks:ke,mxfreq), &
!                asymtot (is:ie,js:je,ks:ke,mxfreq), &
!                sbacktot(is:ie,js:je,ks:ke,mxfreq) )
!    endif
!
!    do k = 1, ke ; do j = js, je ; do i = is, ie
!      tavg = atmos(i,j,k)%t_air
!      tavg_c = tavg - 273.16 ! [degC]
!      pres = atmos(i,j,k)%press
!
!      ! Identify melting layer and fraction
!      if( atmos_stag(i,j,k-1)%t_air.ge.const_tmelt.and.atmos_stag(i,j,k)%t_air.le.const_tmelt )then
!        melt_frac_bulk = ( atmos_stag(i,j,k-1)%t_air - const_tmelt ) / &
!                         ( atmos_stag(i,j,k-1)%t_air - atmos_stag(i,j,k)%t_air )
!      elseif( atmos_stag(i,j,k-1)%t_air.lt.const_tmelt )then
!        melt_frac_bulk = 0.d0
!      elseif( atmos_stag(i,j,k)%t_air.gt.const_tmelt )then
!        melt_frac_bulk = 0.9d0
!      endif
!
!      ! Gas absorption
!      call gas_absorb( freqcy, tavg, pres, atmos(i,j,k)%rh, 0.0_r_size, atm_ext )
!
!      ! Condensate particle optical properties
!      call opt_gce_micro( nf, freqcy, tavg, melt_frac_bulk, q_gce(i,j,k), re_gce(i,j,k), &
!                          kext_gce, salb_gce, asym_gce, pbck_gce )
!
!      ! Total optical properties
!      call total_opt_gce( atm_ext, kext_gce, salb_gce, asym_gce, pbck_gce, &
!                          kexttot(i,j,k,nf), salbtot(i,j,k,nf), asymtot(i,j,k,nf), sbacktot(i,j,k,nf) )
!
!    enddo ; enddo ; enddo

    return

  end subroutine opt_micro_radar

  !>---------------------------------------------------------------------
  !> Compute single-scattering properties 
  subroutine opt_gce_micro( &
        nf, &
        freq, &
        tavg, &
        fmelt,&
        q, & 
        re, &
        kext_gce, &
        salb_gce, &
        asym_gce, &
        pbck_gce )
    implicit none

    integer,intent(in)  :: nf
    real(r_size),intent(in) :: tavg   ! layer temperature
    real(r_size),intent(in) :: freq   ! frequency [GHz]
    real(r_size),intent(in) :: fmelt  ! melt fraction

    type( particle_gce ),intent(in) :: &
        q, & ! mixing ratio [g/m3]
        re   ! effective radius [micron]
    type( particle_gce ),intent(out) :: &
        kext_gce, & ! extinction coefficient [kg-1]
        salb_gce, & ! single scattering albedo 
        asym_gce, & ! asymetry parameter
        pbck_gce    ! backscatter

    ! RAMS
    type particle_rams !rams paticle type
        real(r_size):: cloud1,cloud2,rain,ice1,ice2,snow,graupel,hail
    end type particle_rams

    type ( particle_rams ) :: & !scalar  parameters
        gnu_rams  ! shape parameter for RAMS generalized gamma distirbution

    integer :: ispc  !RAMS spicies index
    real(r_size) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation

    real(r_size),dimension(7,16) :: dstprms
    real(r_size),dimension(7,16) :: rams_dstprms !RAMS habit table
    data dstprms/ &
!  --------------------------------------------------------------------------------------
!   shape     cfmas   pwmas     cfvt    pwvt     dmb0      dmb1     lhcat - habit name
!  --------------------------------------------------------------------------------------
      .5,      524.,     3.,    3173.,     2.,   2.d-6,   40.d-6,  & !1 -cloud
      .5,      524.,     3.,     144.,   .497,   .1d-3,    5.d-3,  & !2 -rain
    .179,     110.8,   2.91,    1538.,   1.00,  15.d-6,  125.d-6,  & !3 -pris col
    .179,  2.739d-3,   1.74,     27.7,   .484,   .1d-3,   10.d-3,  & !4 -snow col
      .5,      .496,    2.4,     16.1,   .416,   .1d-3,   10.d-3,  & !5 -aggreg
      .5,      157.,     3.,     332.,   .786,   .1d-3,    5.d-3,  & !6 -graup
      .5,      471.,     3.,    152.1,   .497,   .8d-3,   10.d-3,  & !7 -hail 
    .429,     .8854,    2.5,   20801.,  1.377,      00,       00,  & !8 -pris hex
   .3183,   .377d-2,     2.,     56.4,   .695,      00,       00,  & !9 -pris den
   .1803,   1.23d-3,    1.8,   1617.9,   .983,      00,       00,  & !10-pris ndl
      .5,     .1001,  2.256,    6239.,   1.24,      00,       00,  & !11-pris ros
    .429,     .8854,    2.5,    30.08,   .563,      00,       00,  & !12-snow hex
   .3183,   .377d-2,     2.,     3.39,   .302,      00,       00,  & !13-snow den
   .1803,   1.23d-3,    1.8,     44.6,   .522,      00,       00,  & !14-snow ndl
      .5,     .1001,  2.256,    125.7,   .716,      00,       00,  & !15-snow ros
      .5,      524.,     3.,   1.26d7,   1.91,  65.d-6,  100.d-6/    !16-cloud2
!  --------------------------------------------------------------------------------------

    ! $(JointSimulator)/DATAFILES/rams.config.1moment
    gnu_rams%cloud1 = 2.0 ! small cloud
    gnu_rams%ice1   = 2.0 ! pristine ice

    rams_dstprms(1:7,1:16) = dstprms(1:7,1:16)

    ! Initialize
    kext_gce%cloud  = 0.d0; salb_gce%cloud  = 0.d0; asym_gce%cloud  = 0.d0; pbck_gce%cloud  = 0.d0  
    kext_gce%rain   = 0.d0; salb_gce%rain   = 0.d0; asym_gce%rain   = 0.d0; pbck_gce%rain   = 0.d0  
    kext_gce%ice    = 0.d0; salb_gce%ice    = 0.d0; asym_gce%ice    = 0.d0; pbck_gce%ice    = 0.d0
    kext_gce%snow   = 0.d0; salb_gce%snow   = 0.d0; asym_gce%snow   = 0.d0; pbck_gce%snow   = 0.d0
    kext_gce%graupel= 0.d0; salb_gce%graupel= 0.d0; asym_gce%graupel= 0.d0; pbck_gce%graupel= 0.d0  
    kext_gce%hail   = 0.d0; salb_gce%hail   = 0.d0; asym_gce%hail   = 0.d0; pbck_gce%hail   = 0.d0

    ! Mie
    ispc = 1 ; cfmas = rams_dstprms(2,ispc) ; pwmas = rams_dstprms(3,ispc)
    call mie_rams_micro('cloud1',freq,tavg,q%cloud,re%cloud,gnu_rams%cloud1, &  ! cloud
                    cfmas, pwmas , 1.0_r_size,  kext_gce%cloud,salb_gce%cloud,asym_gce%cloud,pbck_gce%cloud )
    call mie_gce_micro('qr',freq,tavg,q%rain,re%rain ,rho_gce%rain ,1.0_r_size,&   ! rain
                      kext_gce%rain,salb_gce%rain,asym_gce%rain,pbck_gce%rain)

    ispc = 3 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice -> RAMS ice1 col
    call mie_rams_micro('ice1',freq, tavg, q%ice ,re%ice , gnu_rams%ice1 , &
                   cfmas, pwmas , 0.0_r_size, kext_gce%ice,salb_gce%ice,asym_gce%ice,pbck_gce%ice )

    call mie_gce_micro('qs',freq,tavg,q%snow,re%snow,rho_gce%snow,fmelt, & !snow
                       kext_gce%snow,salb_gce%snow,asym_gce%snow,pbck_gce%snow )

    call mie_gce_micro('qg',freq,tavg,q%graupel,re%graupel, rho_gce%graupel, 0.0_r_size, & !graupel
                       kext_gce%graupel,salb_gce%graupel,asym_gce%graupel,pbck_gce%graupel  )

    call mie_gce_micro('qh', freq, tavg, q%hail, re%hail , rho_gce%hail , 0.0_r_size,& !hail
                        kext_gce%hail,salb_gce%hail,asym_gce%hail,pbck_gce%hail )

    return

  end subroutine opt_gce_micro
    
  !>---------------------------------------------------------------------
  !> Compute single-scattering properties 
  subroutine mie_gce_micro( &
        spc, &
        freq, &
        temp, &
        wc, &
        re, &
        dens, &
        frac_liq, &
        ksca, &
        asca, &
        gsca, &
        pbck )
    implicit none

    character(len=2),intent(in) :: spc !species character
    real(r_size),intent(in) :: freq  ! frequency of radiation [GHz]
    real(r_size),intent(in) :: temp  ! temperature of particles [K]
    real(r_size),intent(in) :: wc    ! water content of condensate [g/m**3]
    real(r_size),intent(in) :: re    ! drop effective radius [micron] 
    real(r_size),intent(in) :: dens  ! density of hydrometero [kg/m**3]
    real(r_size),intent(in) :: frac_liq ! fracion of liquid (1=rain or 0=ice condensates, 0~1=melting ice)
    real(r_size),intent(out) :: ksca  ! extinction coefficient [1/km]
    real(r_size),intent(out) :: asca  ! single-scatter albedo [-]
    real(r_size),intent(out) :: gsca  ! asymmetry factor [-]
    real(r_size),intent(out) :: pbck  ! backscatter phase function/(4*pi) [-]

!--------Local Parameters
    integer :: i    !looping 
    integer :: imax !max for looping
    
    real(r_size) :: n0   ! intercept for the exponetial DSD [1/m**4]
    real(r_size) ::  wave
    real(r_size) :: rad  ! radius of particle [mm]
    real(r_size) :: lam  ! the slope of the distribution [1/m]
    real(r_size) :: num  ! a particle number density per radius increment [1/m**4]
    real(r_size) :: faa  ! fraction of air in ice 
    real(r_size) :: eicere
    real(r_size) :: eiceim
    real(r_size) :: qext
    real(r_size) :: qsca
    real(r_size) :: asym
    real(r_size) :: qbsca
    real(r_size) :: bext
    real(r_size) :: bsca
    real(r_size) :: bsym
    real(r_size) :: bq11
    
    real(r_size) :: eimag, ereal
    
    real(r_size) :: &
         xd     ,&       ! = dble(2.*pi*rad/wave)
         freqd  ,&
         freqhzd,& ! frequency in double precision
         tempd  ,&       ! temperaure in doubpl precision
         sald   ,&
         eicered,& 
         ewatred,&
         eiceimd,&
         ewatimd,&
         qscad  ,&
         qextd  ,&
         asymd  ,&
         qbscad 
    
       
    complex(r_size) eice,ewat,eair !ice, water, and air permittivity
    complex(r_size) ei     !effective ice-air permittivity
    complex(r_size) emelt,emelt_wi,emelt_iw, emelt_em !effective water-ice-air permittivity
    
    complex(r_size) ewatd
    complex(r_size) eid, emeltd
    complex(r_size) crefd, crefd_liq, crefd_ice
    !complex*8 eice,ewat,eair !ice, water, and air permittivity
    !complex*8 ei     !effective ice-air permittivity
    !complex*8 emelt,emelt_wi,emelt_iw, emelt_em !effective water-ice-air permittivity
    !
    !complex*16 ewatd
    !complex*16 eid, emeltd
    !complex*16 crefd, crefd_liq, crefd_ice
    
    real(r_size),parameter :: dr = 0.05    ! increment bin of radius [mm]
    real(r_size),parameter :: densice = 0.917e+3  !solid ice density
    real(r_size),parameter :: densliq = 1.0e+3    !liquid density
    real(r_size),parameter :: wc_unit = 1.0  !unit water content 1 [g/m3]

!
! Assign some useful constants
!
    wave=300./freq

    if( wc.lt. q_min_micro )then
      ksca = 0.d0 ; asca = 0.d0 ; gsca = 0.d0 ; pbck = 0.d0
      return
    endif
!
! If hydrometeors are present, initialize the scattering parameters
!
    bext=0. ; bsca=0. ; bsym=0. ; bq11=0.
!
!     Loop over particle sizes:
!     increments of particle radius are 0.05 mm; the particle
!     size distribution is expressed as a particle number density,
!     num, per radius increment; the original psd is
!     n(D) = n0 * exp(-lam * D) where n is the number density
!     per diameter increment, n0 is the distribution intercept,
!     lam is the slope of the distribution of ln(n(D)), and D is
!     the particle diameter. It is assumed here that n0 is
!     and wc [g/m3] are prescribed, and lam [1/m] is therefore constrained
!     to yield the prescribed water content:
!     wc = integral {n(D) * pi * (D**3) * dens(D) * dD/6}
!     therefore:
!     lam=(n0*pi*dens/wc)**(0.25)

    imax = nint(10./dr)

    ! derive lambda and intercept of exponential DSD for unit water content
    lam = 3.0e+0/2.0e+0/(re*1d-6)  ![1/m]
    n0 = (lam**4)  * wc_unit * 1d-3 / (const_pi*dens) ![1/m4]
!debug
!write(6,*)'lam',lam
!write(6,*)'n0',n0
!debug

    ! size loop
    do i = 0, imax
      rad = dr*0.5_r_size+dr*float(i) ![mm]
      xd  = real(2.0_r_size*const_pi*rad/wave,kind=r_size)
      ! num is in terms of radius - factors of 2 account for conversion d(D) = 2d(r)
      num = 2.0_r_size*n0*exp(-lam*2.0_r_size*rad*(1.d-3))  ![1/m4]

      ! complex refractive index of ice 
      freqd=real(freq,kind=r_size)
      tempd=real(temp,kind=r_size)
      call iceoptic(freqd,tempd,eicered,eiceimd)
      !eicere=REAL(eicered)
      !eiceim=REAL(eiceimd)
      eicere=real(eicered,kind=r_size)
      eiceim=real(eiceimd,kind=r_size)
      eice=cmplx(eicere,eiceim,kind=r_size)
      eair=cmplx(1.0006_r_size,0.0_r_size,kind=r_size)
!debug
!write(6,*)'rad',rad
!write(6,*)'xd',xd
!write(6,*)'freqd',freqd
!write(6,*)'tempd',tempd
!write(6,*)'eicere',eicere
!write(6,*)'eiceim',eiceim
!write(6,*)'eice',eice
!write(6,*)'eair',eair
!debug

      ! Calculate effective dielectric constant of frozen particle using different methods/assumptions. 
      faa =1.d0 - (dens/densice) ! volume fraction of inclusions (air)
      ice_select: select case(ice_refraction_func)
      case(1)  ! air include ice 
        call mg_ellips(1.d0-faa,eair,eice,ei)
        eid=cmplx(ei,kind=r_size)
        crefd_ice=cdsqrt(eid)
        !okzk eid=dcmplx(ei)
        !okzk crefd_ice=cdsqrt(eid)

      case(2)  ! ice include air
        call mg_ellips(faa,eice,eair,ei)
        eid=dcmplx(ei)
        crefd_ice=cdsqrt(eid)

      case(3)  ! homogeneous mixing of ice and air
        call em_ellips( 1.d0-faa, eair, eice, ei)
        eid=dcmplx(ei)
        crefd_ice=cdsqrt(eid)

      end select ice_select

      ! complex refractive index of liquid water
      freqhzd=dble(freq*1.e+9)
      sald=0.d0
      call watoptic(freqhzd,tempd,sald,ewatred,ewatimd) 
      ewatd=dcmplx(ewatred,ewatimd) !complex*16
      crefd_liq=cdsqrt(ewatd)

      ewat=dcmplx(ewatred,ewatimd) !complex*8
!debug
!write(6,*)'ice_refraction_func',ice_refraction_func
!write(6,*)'faa',faa
!write(6,*)'eair',eair
!write(6,*)'eice',eice
!write(6,*)'ei',ei
!write(6,*)'eid',eid
!write(6,*)'cdsqrt(eid)',cdsqrt(eid)
!write(6,*)'ewatd',ewatd
!write(6,*)'crefd_liq',crefd_liq
!write(6,*)'crefd_ice',crefd_ice
!write(6,*)'ewat',ewat
!debug

      ! choose either liquid, ice, or melting
      if(frac_liq == 0.0_r_size) then    ! ice condensates
        crefd = crefd_ice
      elseif(frac_liq== 1.0_r_size) then ! liquid condensates
        crefd = crefd_liq
      else                       ! melting condensates 
        ! melting option
        melt_select: select case(melt_opt)
        case(0)   ! no melting scheme (just use ice dielectric const.)
          crefd = crefd_ice

        case(1)   ! water include ice
          call mg_ellips( 1.-frac_liq, ewat, ei, emelt_wi )
          emeltd=dcmplx(emelt_wi)
          crefd = cdsqrt(emeltd)

        case(2)   ! ice include water     
          call mg_ellips( frac_liq, ei, ewat, emelt_iw)
          emeltd=dcmplx(emelt_iw)
          crefd = cdsqrt(emeltd)

        case(3)   !averaging above two solutions
          call mg_ellips( 1.-frac_liq, ewat, ei, emelt_wi )
          call mg_ellips( frac_liq, ei, ewat, emelt_iw)
          eimag = 10.**( log10(aimag(emelt_wi))*(frac_liq) + log10(aimag(emelt_iw))*(1.-frac_liq)    )
          ereal = 10.**( log10(dble(emelt_wi))*frac_liq + log10(dble(emelt_iw))*(1.-frac_liq)    )
          emelt =dcmplx(ereal,eimag)
          emeltd=dcmplx(emelt)
          crefd =cdsqrt(emeltd)

        case(4)   ! homogeneous mixing of ice and water
          call em_ellips( 1-frac_liq, ewat, ei, emelt_em)
          emeltd=dcmplx(emelt_em)  !ewatd
          crefd =cdsqrt(emeltd)

        end select melt_select
      endif

      ! Mie
      call mie_sphere(xd,crefd,qscad,qextd,asymd,qbscad)
!debug
!write(6,*)'frac_liq',frac_liq
!write(6,*)'melt_opt',melt_opt
!write(6,*)'crefd',crefd
!write(6,*)'qscad',qscad
!write(6,*)'qextd',qextd
!write(6,*)'asymd',asymd
!write(6,*)'qbscad',qbscad
!debug

      qext=REAL(qextd,kind=r_size)
      qsca=REAL(qscad,kind=r_size)
      asym=REAL(asymd,kind=r_size)
      qbsca=REAL(qbscad,kind=r_size)

      ! integrate over particle size distribution; (Cross sectin extionction)*(cross sectional aerea)
      bext=bext+num*qext*const_pi*rad*rad*dr*1.d-6
      bsca=bsca+num*qsca*const_pi*rad*rad*dr*1.d-6
      bsym=bsym+num*qsca*asym*const_pi*rad*rad*dr*1.d-6
      bq11=bq11+num*qbsca*const_pi*rad*rad*dr*1.d-6
!debug
!write(6,*)'bext',bext
!write(6,*)'bsca',bsca
!write(6,*)'bsym',bsym
!write(6,*)'bq11',bq11
!debug

    enddo ! i

    !
    ! check for distribution with very small extinction;
    ! set parameters to zero to avoid numerical problems
    !
    if( bext .gt. min_bext) then
      ksca=bext      * wc  ! *wc is re-scaling to actual extinction
      asca=bsca/bext
      gsca=bsym/bsca
      pbck=bq11/bsca
    else
      ksca=0. ; asca=0. ; gsca=0. ; pbck=0.
    end if
!debug
!write(6,*)'ksca',bext
!write(6,*)'asca',bsca
!write(6,*)'gsca',bsym
!write(6,*)'pbck',bq11
!debug

    return

  end subroutine mie_gce_micro

  !>---------------------------------------------------------------------
 subroutine mg_ellips (fincl, ematrix, eincl, emg)
 implicit none
!
! Maxwell-Garnett dielctric function of 2-component media
! (elliptical inclusions) P. Bauer 1996 (core-shell mixture)
! See also Olson et al. 2001 JAM page 1153 equation (32)&(33)
! This version is closer to EM sollution, and slighly different 
! from Bohren and Battan 1980 (spherical inclusion). 
! 

 real(r_size),intent(in) ::  fincl  !volume fraction of inclusions
 complex(r_size),intent(in) :: ematrix ! permittivity of matrix
 complex(r_size),intent(in) :: eincl   ! permittivity of inclusions
 complex(r_size),intent(out) ::  emg   !  effective permittivity
 complex(r_size) :: gamma, q
 !complex*8,intent(in) :: ematrix ! permittivity of matrix
 !complex*8,intent(in) :: eincl   ! permittivity of inclusions
 !complex*8,intent(out) ::  emg   !  effective permittivity
 !complex*8 :: gamma, q


      q     = (eincl / (eincl - ematrix)) &
            * cdlog (eincl / ematrix) - 1.0
      gamma = 2.0 * ematrix * q / (eincl - ematrix)

      emg = ((1.0 - fincl) * ematrix + fincl * gamma * eincl) &
          / (1.0 - fincl + fincl * gamma)

    return
  end subroutine mg_ellips

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine em_ellips (fincl, ematrix, eincl, emg)
 implicit none

!
! Effective-Medium dielectric function  (Homogeneous mixture assumption)
! See eq. (2) in Bohren and Battan 1980
!
! Toshi Matsui@NASA GSFC; Initial
!

 real(r_size),intent(in) ::  fincl  !volume fraction of inclusions
 complex(r_size),intent(in) :: ematrix ! permittivity of matrix
 complex(r_size),intent(in) :: eincl   ! permittivity of inclusions
 complex(r_size),intent(out) :: emg   !  effective permittivity
 complex(r_size) :: a,b,c
 !complex*8,intent(in) :: ematrix ! permittivity of matrix
 !complex*8,intent(in) :: eincl   ! permittivity of inclusions
 !complex*8,intent(out) :: emg   !  effective permittivity
 !complex*8 :: a,b,c

!
! parameters for a quadratic equation
!
 a = cmplx(-2.,0.,kind=r_size)
 b = (3.*fincl-1.)*eincl - (3.*fincl-2.)*ematrix 
 c = ematrix*eincl

!
! quadratic solution for effective permittivity 
!
 emg = (-b - sqrt(b*b-4.*a*c) )  / (2.*a)

    return
  end subroutine em_ellips

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine mie_sphere (x, min, qscat, qexti, asym, qbscat)

!*
!*    Mie Routine P. Bauer
!*

      implicit none
!*
      integer limitx
      parameter (limitx = 1500)
      real(r_size) conv, pi, ghztohz
!*
      real(r_size) :: x
      real(r_size) :: mr, mi, n1, n2
      real(r_size) :: qscat, qexti, qabso, asym, qbscat
!*
      real(r_size) :: rfac1, rfac2
      real(r_size) :: rhelp1 (2), rhelp2 (2)
!*
      complex(r_size) m, mx, min
      complex(r_size) chelp1, chelp2, cfac1, cfac2, cbscat
!*
      complex(r_size) dn (0:limitx), wn (-1:limitx)
      complex(r_size) an   (limitx), bn    (limitx)
!*
      integer nend
      integer i100, i101
!*
      equivalence (chelp1, rhelp1 (1))
      equivalence (chelp2, rhelp2 (1))
!*
      real(r_size) ::  eps0, epshigh
!*
!*
      data pi      /  3.14159265358979323846d+00 /
      data conv    /  0.017453292d+00 /
      data ghztohz /  1.0d+09 /
!*
      data eps0     / 8.854d-12 /
      data epshigh  / 4.9d+00 /

!*
      m      = dconjg (min)
      chelp1 = m
      mr     =        rhelp1 (1)
      mi     = -1.0d0 * rhelp1 (2)
!*
      mx   = m  * x
      n1   = mr * x
      n2   = mi * x
!*
      if (x .le. 20000.0) nend = x + 4.00 * x ** (1.0 / 3.0) + 2.0
      if (x .le.  4200.0) nend = x + 4.05 * x ** (1.0 / 3.0) + 2.0
      if (x .le.     8.0) nend = x + 4.00 * x ** (1.0 / 3.0) + 1.0
      if (nend .le.    5) nend = 5
      if (nend .gt. limitx) nend = limitx
!*
      rfac1      = sin  (n1) * sin  (n1) + sinh (n2) * sinh (n2)
      rhelp1 (1) = sin  (n1) * cos  (n1) / rfac1
      rhelp1 (2) = sinh (n2) * cosh (n2) / rfac1
!*
      dn (0) = chelp1
!*
      rhelp1 (1) =             cos (x)
      rhelp1 (2) = -1.0d+00 * sin (x)
      rhelp2 (1) =             sin (x)
      rhelp2 (2) =             cos (x)
!*
      wn (-1) = chelp1
      wn ( 0) = chelp2
!*
      qexti  = 0.0
      qscat  = 0.0
      qbscat = 0.0
      qabso  = 0.0
      asym   = 0.0
      cbscat = cmplx (0.0,0.0,kind=r_size)
!*
      do 100 i100 = 1, nend
         dn (i100) = -1.0d0 * i100 / mx &
                   +  1.0d0 / (i100 / mx - dn (i100 - 1))
         wn (i100) = wn (i100 - 1) * (2.0d0 * i100 - 1.0d0) / x &
                   - wn (i100 - 2)
!*
         cfac1 = dn (i100) / m + i100 / x
         cfac2 = m * dn (i100) + i100 / x
!*
         chelp1 = wn (i100)
         chelp2 = wn (i100 - 1)
!*
         an (i100) = (cfac1 * rhelp1 (1) - rhelp2 (1)) &
                   / (cfac1 * chelp1     - chelp2    )
         bn (i100) = (cfac2 * rhelp1 (1) - rhelp2 (1)) &
                   / (cfac2 * chelp1     - chelp2    )
!*
         chelp1 = an (i100)
         chelp2 = bn (i100)
!*
         rfac1 = rhelp1 (1) + rhelp2 (1)
         rfac2 = cdabs (an (i100)) * cdabs (an (i100)) &
               + cdabs (bn (i100)) * cdabs (bn (i100))
!*
         qexti  = qexti  + (2.0d0 * i100 + 1.0d0) * rfac1
         qscat  = qscat  + (2.0d0 * i100 + 1.0d0) * rfac2
         cbscat = cbscat + (2.0d0 * i100 + 1.0d0) * (-1.0d0) ** i100 &
                * (an (i100) - bn (i100))
!*
         if (i100 .eq. 1) go to 100
!*
         chelp1 = an (i100 - 1) * dconjg (an (i100)) &
                + bn (i100 - 1) * dconjg (bn (i100))
         chelp2 = an (i100 - 1) * dconjg (bn (i100 - 1))
!*
         i101 = i100 - 1
         rfac1  = i101 * (i101 + 2) / (i101 + 1.0d0)
         rfac2  = (2.0d0 * i101 + 1.0d0) / (i101 * (i101 + 1.0d0))
!*
         asym = asym + rfac1 * rhelp1 (1) + rfac2 * rhelp2 (1)
100   continue
!*
      qexti  = qexti * 2.0d0 / (x * x)
      qscat  = qscat * 2.0d0 / (x * x)
      asym   = asym  * 4.0d0 / (x * x * qscat)
      qbscat = cdabs (cbscat) * cdabs (cbscat) / (x * x)
      if (qscat .gt. qexti) qscat = qexti
!*
   return
 end subroutine mie_sphere

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine gas_absorb(freq,temp,pres,relhum,water,kabs)
      implicit none

!----IO parameter----
      real(r_size), intent(in) :: freq   !frequency [GHz]
      real(r_size), intent(in) :: temp   !temperature [K]
      real(r_size), intent(in) :: relhum !relative humidity [%]
      real(r_size), intent(in) :: pres   !pressure [mb] 
      real(r_size), intent(in) :: water  !water content [g/m3]
      real(r_size), intent(out) :: kabs  !absorption coefficient
   
!----local parameter-------
      real(r_size),parameter:: a0=6.107799961d0
      real(r_size),parameter:: a1=4.436518521d-1
      real(r_size),parameter:: a2=1.428945805d-2
      real(r_size),parameter:: a3=2.650648471d-4
      real(r_size),parameter:: a4=3.031240396d-6
      real(r_size),parameter:: a5=2.034080948d-8
      real(r_size),parameter:: a6=6.136820929d-11

      real(r_size) :: tc !temperature [DegK]
      real(r_size) :: es !saturation papor pressure
      real(r_size) :: vapor_pressure !vapor pressure 
      real(r_size) :: rho            !water vapor density [g/m3] 
      real(r_size) :: kabs_h2o , kabs_o2 , kabs_n2 , kabs_clw  !absorption coefficient in each component. 

      tc = temp - 273.15
      es = a0+tc*(a1+tc*(a2+tc*(a3+tc*(a4+tc*(a5+a6*tc)))))
      if (es .lt. 0.) es = 0.
      vapor_pressure =  relhum*es/100.
      rho = vapor_pressure*100*18/(8.314*temp)

      call abs_h2o(temp, pres, rho, freq, kabs_h2o)

      call abs_o2(temp, pres, rho, freq, kabs_o2)

      call abs_n2(temp, pres, freq, kabs_n2)

      call absorb_clw(freq,temp,water,kabs_clw)

      kabs = kabs_h2o + kabs_o2 + kabs_n2 + kabs_clw

      return
      end subroutine gas_absorb

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine abs_h2o(t,p,rho,f,abh2o)
      implicit none
!
!     C. Kummerow, 8/2003.  Changed function to subroutine
!     Copyright (c) 2002 Massachusetts Institute of Technology
!
!  NAME- ABH2O    LANGUAGE- FORTRAN 77
!
! PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
!
!  CALLING SEQUENCE PARAMETERS-
!    SPECIFICATIONS
      real(r_size) :: t,p,rho,f,abh2o
!      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
!      T       KELVIN    I   TEMPERATURE
!      P       MILLIBAR  I   PRESSURE              .1 TO 1000
!      RHO     G/M**3    I   WATER VAPOR DENSITY
!      F       GHZ       I   FREQUENCY             0 TO 800
!      ABH2O   NEPERS/KM O   ABSORPTION COEFFICIENT
!
!   REFERENCES-
!   P.W. ROSENKRANZ, RADIO SCIENCE V.33, PP.919-928 (1998); V.34, P.1025 (1999).
!
!   LINE INTENSITIES SELECTION THRESHOLD=
!     HALF OF CONTINUUM ABSORPTION AT 1000 MB.
!   WIDTHS MEASURED AT 22,183,380 GHZ, OTHERS CALCULATED.
!     A.BAUER ET AL.ASA WORKSHOP (SEPT. 1989) (380GHz).
!     M. TRETYAKOV et al., J. MOLEC. SPEC. (2003)
!
!   REVISION HISTORY-
!    DATE- OCT.6, 1988  P.W.ROSENKRANZ - EQS AS PUBL. IN 1993.
!          OCT.4, 1995  PWR- USE CLOUGH'S DEFINITION OF LOCAL LINE
!                   CONTRIBUTION,  HITRAN INTENSITIES, ADD 7 LINES.
!          OCT. 24, 95  PWR -ADD 1 LINE.
!          JULY 7, 97   PWR -SEPARATE COEFF. FOR SELF-BROADENING,
!                       REVISED CONTINUUM.
!        Aug. 28, 2002  PWR - CORRECTED LINE INTENSITIES
!        Mar. 2, 2003   PWR - LINE SHIFT
!
!   LOCAL VARIABLES:
      integer nlines,i,j
      parameter (nlines=15)
      real(r_size) :: df(2),s1(nlines),b2(nlines),w3(nlines),fl(nlines),x(nlines), &
       ws(nlines),xs(nlines),sr(nlines)
      real(r_size) :: pvap,pda,den,ti,ti2,sum,width,wsq,s,base,res,con,shift
!     LINE FREQUENCIES:
      data fl/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, &
       443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360, &
       620.7008, 752.0332, 916.1712/
!     LINE INTENSITIES AT 300K:
      data s1/ .1314d-13, .2279d-11, .8058d-13, .2701d-11, .2444d-10, &
       .2185d-11, .4637d-12, .2568d-10, .8392d-12, .3272d-11, .6676d-12, &
       .1535d-08, .1711d-10, .1014d-08, .4238d-10/
!     T COEFF. OF INTENSITIES:
      data b2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, &
       3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
!     AIR-BROADENED WIDTH PARAMETERS AT 300K:
      data w3/.00281, .00287, .0023, .00278, .00287, .0021, .00186, &
       .00263, .00215, .00236, .0026, .00321, .00244, .00306, .00267/
!     T-EXPONENT OF AIR-BROADENING:
      data x/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69, &
       .71, .68, .70/
!     SELF-BROADENED WIDTH PARAMETERS AT 300K:
      data ws/.01349, .01491, .0108, .0135, .01541, .0090, .00788, &
       .01275, .00983, .01095, .01313, .01320, .01140, .01253, .01275/
!     T-EXPONENT OF SELF-BROADENING:
      data xs/ .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72, &
       1.0, .68, .84, .78/
!     RATIO OF SHIFT TO WIDTH
      data sr/ 0., -.017, 13*0./
!
      if(rho.le.0.) then
        abh2o = 0.
        return
      endif
      pvap = rho*t/217.
      pda = p -pvap
      den = 3.335d16*rho ! const includes isotopic abundance
      ti = 300./t
      ti2 = ti**2.5
!
!      CONTINUUM TERMS
      con = (5.43d-10*pda*ti**3 + 1.8d-8*pvap*ti**7.5)*pvap*f*f
!
!      ADD RESONANCES
      sum = 0.
      do 30 i=1,nlines
      width = w3(i)*pda*ti**x(i) + ws(i)*pvap*ti**xs(i)
      shift = sr(i)*width  ! unknown temperature dependence
      wsq = width*width
      s = s1(i)*ti2*exp(b2(i)*(1.-ti))
      df(1) = f - fl(i) - shift
      df(2) = f + fl(i) + shift
!  USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
      base = width/(562500. + wsq)
!  DO FOR POSITIVE AND NEGATIVE RESONANCES
      res = 0.
      do 20 j=1,2
      if(abs(df(j)).lt.750.) res = res + width/(df(j)**2+wsq) - base
20    continue
30    sum = sum + s*res*(f/fl(i))**2
      abh2o = .3183d-4*den*sum + con
      return
      end subroutine abs_h2o

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine abs_o2(temp,pres,vapden,freq,o2abs)
      implicit none
!
!     C. Kummerow, 8/2003.  Changed function to subroutine
!  Copyright (c) 2003 Massachusetts Institute of Technology
!
!     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
!              IN NEPERS/KM
!
!      5/1/95  P. Rosenkranz
!      11/5/97  P. Rosenkranz - 1- line modification.
!      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
!      8/21/02  pwr - revised width at 425
!      3/20/03  pwr - 1- line mixing and width revised
!
!     IMPLICIT NONE
!
!     ARGUMENTS:
    real(r_size) :: temp,pres,vapden,freq
    real(r_size),intent(out) :: o2abs
!
!     NAME    UNITS    DESCRIPTION        VALID RANGE
!
!     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
!                                          valid for atmosphere
!     PRES   MILLIBARS PRESSURE           3 TO 1000
!     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
!                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
!     FREQ    GHZ      FREQUENCY          0 TO 900
!
!     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
!     P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
!      BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
!     H.J. Liebe et al, JQSRT V.48, pp.629-643 (1992).
!     M.J. Schwartz, Ph.D. thesis, M.I.T. (1998).
!     A.F. Krupnov et al, J. Mol. Spect. v.215, pp.309-311 (2002).
!     M.Yu. Tretyakov et al, J. Mol. Spect. (2003 preprint).
!     SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
!
!     This version differs from Liebe's MPM92 in these significant respects:
!     1. The 1- line has the width and mixing coefficient measured by
!      Tretyakov et al.
!     2. It modifies the 1- line width temperature dependence to (1/T)**0.9
!     3. It uses the same temperature dependence (X) for submillimeter
!      line widths as in the 60 GHz band: (1/T)**0.8
!     4. The 425 GHz line width is from Krupnov et al.
!
!     Local variables:
      real(r_size) :: th,th1,b,preswv,presda,den,dens,dfnr,sum,str,y,sf1,sf2,fcen,df
      integer k
      real(r_size) :: x,wb300,w300(40),f(40),y300(40),s300(40),v(40),be(40)
      common /o2com/ x,wb300,w300,f,y300,s300,v,be
!      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
      data f/118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910, &
        59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002, &
        56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685, &
        55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241, &
        53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368, &
        52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7632, &
        487.2494, 715.3931, 773.8397, 834.1458/
        data s300/.2936d-14,.8079d-15, .2480d-14,.2228d-14, &
        .3351d-14,.3292d-14, .3721d-14,.3891d-14, &
        .3640d-14,.4005d-14, .3227d-14,.3715d-14, &
        .2627d-14,.3156d-14, .1982d-14,.2477d-14, &
        .1391d-14,.1808d-14, .9124d-15,.1230d-14, &
        .5603d-15,.7842d-15, .3228d-15,.4689d-15, &
        .1748d-15,.2632d-15, .8898d-16,.1389d-15, &
        .4264d-16,.6899d-16, .1924d-16,.3229d-16, &
        .8191d-17,.1423d-16, .6494d-15, .7083d-14, .3025d-14, &
        .1835d-14, .1158d-13, .3993d-14/
      data be/.009,.015, .083,.084, 2*.212, 2*.391, 2*.626, &
       2*.915, 2*1.260, 1.660,1.665, 2.119,2.115, 2.624,2.625, &
       2*3.194, 2*3.814, 2*4.484, 2*5.224, 2*6.004, 2*6.844, &
       2*7.744, .048, .044, .049, .145, .141, .145/
!      WIDTHS IN MHZ/MB
      data wb300/.56/, x/.8/
      data w300/1.67, 1.646, 1.468, 1.449, 1.382, 1.360, &
       1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171, &
       1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 2*1.05, &
       2*1.02,2*1.00,2*.97,2*.94,2*.92,2*.89, 3*1.64, 3*1.81/
      data y300/  -0.036,  0.2408, -0.3486,  0.5227, &
       -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311, &
        0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599, &
        0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246, &
        0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546, &
        0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529, 6*0./
      data v/  0.0079, -0.0978,  0.0844, -0.1273, &
        0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584, &
        0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675, &
        0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590, &
        0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091, &
        0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, 6*0./
!
      th = 300./temp
      th1 = th-1.
      b = th**x
      preswv = vapden*temp/217.
      presda = pres -preswv
      den = .001*(presda*b + 1.1*preswv*th)
      dens = .001*(presda*th**.9 + 1.1*preswv*th)
      dfnr = wb300*den
      sum = 1.6d-17*freq*freq*dfnr/(th*(freq*freq + dfnr*dfnr))
      do 32 k=1,40
      if(k.eq.1) then !exception for 1- line
        df = w300(1)*dens
      else
        df = w300(k)*den
      endif
      fcen = f(k)
      y = .001*pres*b*(y300(k)+v(k)*th1)
      str = s300(k)*exp(-be(k)*th1)
      sf1 = (df + (freq-fcen)*y)/((freq-fcen)**2 + df*df)
      sf2 = (df - (freq+fcen)*y)/((freq+fcen)**2 + df*df)
32    sum = sum + str*(sf1+sf2)*(freq/f(k))**2
      o2abs = .5034d12*sum*presda*th**3/3.14159
      o2abs = amax1(o2abs,0.)
      return
      end subroutine abs_o2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

   subroutine abs_n2(t,p,f,absn2)
   implicit none

   real(r_size),intent(in) :: t ! temperature [K]
   real(r_size),intent(in) :: p ! pressure [mb]
   real(r_size),intent(in) :: f ! frequency [GHz] (valid 0 - 1000GHz)
   real(r_size),intent(out) :: absn2 
   real(r_size) :: fdepen, bf, th

!
!     C. Kummerow, 8/2003.  Changed function to subroutine
!  Copyright (c) 2002 Massachusetts Institute of Technology
!     ABSN2 = COLLISION-INDUCED ABSORPTION COEFFICIENT (NEPER/KM)
!     5/22/02 P.Rosenkranz
!
!     Equations based on:
!      Borysow, A, and L. Frommhold,
!      Astrophysical Journal, v.311, pp.1043-1057 (1986)
!     with modification of 1.29 to account for O2-O2 and O2-N2
!     collisions, as suggested by
!      J.R. Pardo, E.Serabyn, J.Cernicharo, J. Quant. Spectros.
!      Radiat. Trans. v.68, pp.419-433 (2001).
!

      th = 300./t

      fdepen = .5 + .5/(1.+(f/450.)**2)

      bf = 6.5d-14*fdepen*p*p*f*f*th**3.6

      absn2 = 1.29*bf

      return

      end subroutine abs_n2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine absorb_clw(freq,temp,water,kabs_clw)
      implicit none

      real(r_size),intent(in) :: freq !frequency [GHz]
      real(r_size),intent(in) :: temp !temperature [K]
      real(r_size),intent(in) :: water !water content [g/m3]

!     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
!     ARGUMENTS (INPUT):
!     WATER IN G/M**3
!     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
!     TEMP IN KELVIN

      real(r_size),intent(out) :: kabs_clw
!
!     REFERENCES:
!     LIEBE, HUFFORD AND MANABE, INT. J. IR & MM WAVES V.12, pp.659-675
!      (1991);  Liebe et al, AGARD Conf. Proc. 542, May 1993.
!
!     REVISION HISTORY:
!        PWR 8/3/92   original version
!        PWR 12/14/98 temp. dependence of EPS2 eliminated to agree
!                     with MPM93
!        pwr 2/27/02  use exponential dep. on T, eq. 2b instead of eq. 4a
!
    real(r_size) :: eps0, eps1, eps2
    real(r_size) :: theta1, fp, fs
    complex eps,re
      if(water.le.0.) then
       kabs_clw = 0.
       return
      endif
      theta1 = 1.-300./temp
      eps0 = 77.66 - 103.3*theta1
      eps1 = .0671*eps0
      eps2 = 3.52                 ! from MPM93
      fp = 20.1*exp(7.88*theta1)  ! from eq. 2b
      fs = 39.8*fp
      eps = (eps0-eps1)/cmplx(1.,freq/fp,kind=r_size) + &
       (eps1-eps2)/cmplx(1.,freq/fs,kind=r_size) +eps2
      re = (eps-1.)/(eps+2.)
      kabs_clw = -.06286*aimag(re)*freq*water
      return
      end subroutine absorb_clw

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine iceoptic (fghz, t, epsri, epsii)
      implicit none
!
!     Compute  complex refractive index of cloud ice
!
!  Hufford (1991), see Brussard and Watson (1995), p.297
!  The radio frequency may range from 0 to 1000 GHz, and the temperature from 0° to –40°C.
!

      real(r_size) :: fghz, epsri, epsii, a, b, theta, t, tk
!
      epsri = 3.15d+00

      tk = t
      if (tk .gt. 273.16) tk = 273.16  !always below freezing 

      theta = 300.0 / tk
      a     = 1.0d-04 * (50.4 + 62.0 * (theta - 1.0)) &
            * exp (-22.1 * (theta - 1.0))
      b     = 1.0d-04 * (0.633 / theta - 0.131) &
            + (7.36d-04 * theta / (theta - 0.9927)) &
            * (7.36d-04 * theta / (theta - 0.9927))
      epsii = a / fghz + b * fghz
!
      return
      end subroutine iceoptic 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine watoptic (freqhz, ktemp, salinity, epsreal, epsimag)
      implicit none
!
! complex refractive index (permittivity) of liquid water
!

      real(r_size) :: fac1
      real(r_size) :: conv, pi, ghztohz
      real(r_size) :: ktemp, ctemp, salinity
      real(r_size) :: eps0, epsstat, epshigh, epsreal, epsimag
      real(r_size) :: trelax, freqhz, omega
      data pi      /  3.141592654d+00 /
      data conv    /  0.017453292d+00 /
      data ghztohz /  1.0d+09 /
      data eps0     / 8.854d-12 /  !permittivity in vacume
      data epshigh  / 4.9d+00 /  !high-frequency permittivity


!*
      ctemp = ktemp - 273.15d+00
      omega = (pi + pi) * freqhz
!
! parameterized static permittivity (v->0)
!
      epsstat = (87.134d+00 &
              - 1.949d-01 * ctemp &
              - 1.276d-02 * ctemp * ctemp &
              + 2.491d-04 * ctemp * ctemp * ctemp) &
              * (1.0d+00 &
              + 1.613d-05 * salinity * ctemp &
              - 3.656d-03 * salinity &
              + 3.210d-05 * salinity * salinity &
              - 4.232d-07 * salinity * salinity * salinity)
!
! parameterized effective relaxation time
!
      trelax = (1.768d-11 &
             - 6.086d-13 * ctemp &
             + 1.104d-14 * ctemp * ctemp &
             - 8.111d-17 * ctemp * ctemp * ctemp) &
             * (1.0d+00 &
             + 2.282d-05  * salinity * ctemp &
             - 7.638d-04  * salinity &
             - 7.760d-06  * salinity * salinity &
             + 1.105d-08  * salinity * salinity * salinity)

      fac1    = 1.0d+00 + omega * omega * trelax * trelax

!
! real and imaginery components of permittivity
!
      epsreal = epshigh + (epsstat - epshigh) / fac1
      epsimag = ((epsstat - epshigh) * omega * trelax) / fac1

      return
      end subroutine watoptic

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine total_opt_gce(atm_ext, kext_gce, salb_gce, asym_gce, pbck_gce, &
                            kexttot, salbtot, asymtot, sbacktot )
 implicit none

 real(r_size),intent(in) :: atm_ext  !atmosphere extinction 
 type ( particle_gce ),intent(in) :: & ! particle_gce is defined in module_simulater
     kext_gce,  & ! extinction coefficient
     salb_gce,  & ! single scattering albedo [-]
     asym_gce,  & ! asymetry parameter [-]
     pbck_gce     ! backscattering coefficient

 real(r_size),intent(out) ::  kexttot, salbtot, asymtot, sbacktot !total optical parameters

!
! total extinction 
!
   kexttot    = atm_ext         + &
                kext_gce%cloud  + & 
                kext_gce%rain   + &
                kext_gce%ice    + &
                kext_gce%snow   + &
                kext_gce%graupel+ &
                kext_gce%hail 

!
! total single scattering albedo
!
   if ( kexttot <= 0.d0 ) then
        salbtot = 0.d0
   else
       salbtot = ( &
                  salb_gce%cloud  * kext_gce%cloud   + &
                  salb_gce%rain   * kext_gce%rain    + &
                  salb_gce%ice    * kext_gce%ice     + &
                  salb_gce%snow   * kext_gce%snow    + &
                  salb_gce%graupel* kext_gce%graupel + &
                  salb_gce%hail   * kext_gce%hail      &
                 ) / kexttot

   endif

!
! total asymetry parameter
!
   if ( salbtot <= 0.d0 ) then
       asymtot = 0.d0
   else
       asymtot = ( &
                  asym_gce%cloud  * salb_gce%cloud  * kext_gce%cloud   + &
                  asym_gce%rain   * salb_gce%rain   * kext_gce%rain    + &
                  asym_gce%ice    * salb_gce%ice    * kext_gce%ice     + &
                  asym_gce%snow   * salb_gce%snow   * kext_gce%snow    + &
                  asym_gce%graupel* salb_gce%graupel* kext_gce%graupel + &
                  asym_gce%hail   * salb_gce%hail   * kext_gce%hail      &
                 ) / ( salbtot * kexttot )

   endif

!
! total backscattering coef.  (particle only, does not include atmosphere)
!
   sbacktot = &
                  pbck_gce%cloud  * salb_gce%cloud  * kext_gce%cloud   + &
                  pbck_gce%rain   * salb_gce%rain   * kext_gce%rain    + &
                  pbck_gce%ice    * salb_gce%ice    * kext_gce%ice     + &
                  pbck_gce%snow   * salb_gce%snow   * kext_gce%snow    + &
                  pbck_gce%graupel* salb_gce%graupel* kext_gce%graupel + &
                  pbck_gce%hail   * salb_gce%hail   * kext_gce%hail      


 return
 end subroutine total_opt_gce


  !>---------------------------------------------------------------------
  subroutine mie_rams_micro(spc,freq, temp, wc, re, gnu, cfmas, pwmas, frac_liq, ksca, asca, gsca, pbck)
  implicit none

!--------------------------------------------------------------------------------
! Purpose : Compute the extinction, absorption, asymmetry parameter and
!           backscatter for a given water content and effective radius of condensates , and
!           a particle size distribution parameter for RAMS 1moment & 2moment microphyiscs
!
!  03/2010  Toshi Matsui@NASA GSFC ; Additional options for effective refractive index for solid particles.
! Toshi Matsui @ NASA GSFC: Initial
!
!--------------------------------------------------------------------------------
 character*(*),intent(in) :: spc !species character
 real(r_size),intent(in) :: freq  ! frequency of radiation [GHz]
 real(r_size),intent(in) :: temp  ! temperature of particles [K]
 real(r_size),intent(in) :: wc    ! water content of condensate [g/m3]
 real(r_size),intent(in) :: re    ! drop effective radius [micron] 
 real(r_size),intent(in) :: gnu   ! PSD shape parameter for generalized gamma distribution

 real(r_size),intent(in) :: cfmas  ! alpha_m [kg /m**beta_m]
 real(r_size),intent(in) :: pwmas  ! beta_m  [-]
 real(r_size),intent(in) :: frac_liq ! fracion of liquid (1=rain or 0=ice condensates, 0~1=melting snow)
 real(r_size),intent(out) :: ksca  ! extinction coefficient [1/km]
 real(r_size),intent(out) :: asca  ! single-scatter albedo [-]
 real(r_size),intent(out) :: gsca  ! asymmetry factor [-]
 real(r_size),intent(out) :: pbck  ! backscatter phase function/(4*pi) [-]

!--------Local Parameters
 integer :: i    !looping 
 integer :: imax !max for looping


 real(r_size) ::& 
   wave  ,& ! wavelength [mm]
   rad   ,& ! radius of particle [mm]
   d     ,& !particle diameter 
   dens  ,& !particle density [kg/m3]
   num   ,& ! a particle number density per radius increment [1/m4]
   faa   ,& ! fraction of air in ice 
   eicere,& ! real component of ice complex refractive index
   eiceim,& !imaginary component of ice complex refractive index
   qext  ,& !extinction
   qsca  ,&!scattering
   asym  ,& !asymetry 
   qbsca ,& !backscattering
   bext  ,& !extinction
   bsca  ,& !scattering
   bsym  ,& !asymetry
   bq11  ,& !backscattering
   eimag ,& !complex refractive index
   ereal ,& !
   gamfac,& ! gamma PSD factors [-]
   gfac1 ,& !
   gfac2 ,&   !
   dn         ,& ! charactristic diameter [m]
   d_increment,& ! diameter increment over size loops [mm]
   d_max      ,& ! maximum diameter of particle
   d_m        ,& ! particle diameter [mm]
   d_mm       ,& ! particle diameter [mm]
   rad_new    ,& ! updated particle radius [mm]
   fmelt      ,& ! melt fraction
   fmelt_new  ,& ! updated melt fraction [-]
   mean_mass  ,& ! mean mass [kg]
   fgamma     ,& ! gamma function
   ntot       ,& ! total particle number concentrations [#/m3]
   factor        ! factor = num * pi * 0.25 * d_mm * d_mm * d_increment * 1.e-6

 real(r_size) :: &
     xd      ,&      ! size parameter [-] in double precision 
     freqd   ,&
     freqhzd ,& ! frequency in double precision [GHz]
     tempd   ,&      ! temperaure in doubpl precision
     sald    ,&      ! salinity of water
     eicered ,&
     ewatred ,& !real component of ice/water complex refractive index
     eiceimd ,&
     ewatimd ,& !imaginary component of ice/water complex refractive index
     qscad   ,&
     qextd   ,&
     asymd   ,&
     qbscad
    
 !complex *8 eice,ewat,eair !ice, water, and air permittivity
 !complex *8 ei     !effective ice-air permittivity
 !complex *8 emelt,emelt_wi,emelt_iw, emelt_em !effective water-ice-air permittivity
 complex(r_size) eice,ewat,eair !ice, water, and air permittivity
 complex(r_size) ei     !effective ice-air permittivity
 complex(r_size) emelt,emelt_wi,emelt_iw, emelt_em !effective water-ice-air permittivity

 !complex*16 ewatd
 !complex*16 eid, emeltd
 !complex*16 crefd, crefd_liq, crefd_ice
 complex(r_size) ewatd
 complex(r_size) eid, emeltd
 complex(r_size) crefd, crefd_liq, crefd_ice

 real(r_size),parameter :: dr = 0.05          ! increment bin of radius [mm]
 real(r_size),parameter :: densice = 0.917e+3 ! solid ice density [kg/m3]
 real(r_size),parameter :: densliq = 1.0e+3   ! liquid water density [kg/m3]
 real(r_size),parameter :: wc_unit = 1.0      ! unit water content [g/m3]


!
! Begin by checking if hydrometeors of this species have enough water content.
! If not, set scattering parameters to zero and return.
!
  if(wc .lt. q_min_micro) then
     ksca=0. ; asca=0. ; gsca=0. ; pbck=0.
     return
  endif

!
! Assign some useful constants
!

  wave=300./freq !wavelgnth [mm] <- frequency [GHz]


!
! Get characteristic diameter from effective radius (re) and PSD shape parameter (gnu)
!
  call gamma_reff(gnu+3.0_r_size,gfac1)
  call gamma_reff(gnu+2.0_r_size,gfac2)    
  gamfac = gfac1/gfac2
  dn = 2.0 * re / gamfac * 1d-6  !charactristic diameter [m]


!
! get total number concentration for a given characteristic diameter, mass mixing ratio, PSD shape parameter
! , and gamma factors.
!
  call gamma_reff(gnu+pwmas,gfac1)
  call gamma_reff(gnu      ,gfac2)
  gamfac = gfac1/gfac2
  mean_mass = cfmas * ( dn ** pwmas ) * gamfac  !mean mass [kg]
  ntot = wc * 1d-3 / mean_mass                   !total particle number concentration [#/m3]


!
! initialize the optical parameters
!
  bext=0. ; bsca=0. ; bsym=0. ; bq11=0.

  d_increment = 0.10  ! [mm]
  d_max = 20.
! define imax for a givne radius increment 20. is maximum diameter
  imax = nint(d_max/d_increment)

!
! loop over particle size
!
  SIZE_LOOP: do i=0,imax

!
! Compute diameter, size parameter particle number density, and particle density
!
     d_mm = d_increment * 0.5 + d_increment * FLOAT(i)  !diameter [mm]
     d_m  = d_mm * 1d-3                        ! diameter but using different unit [m]
     xd  = dble(const_pi*d_mm/wave)                  ! size parameter [-]
     dens = (6.0 * cfmas / const_pi) * (d_m ** (pwmas-3.d0)) ! particle density [kg/m3]

     fmelt = frac_liq  !melt fraction

!
! Snow compatction impact (Li-Matsui scheme)
!

     if(trim(spc) == 'snow' .and. frac_liq > 0.0 .and. melt_opt /= 0) then
        call snow_compaction( 0.5*d_mm , dens*1d-3 , frac_liq , rad_new, fmelt_new )
        d_mm = 2.*rad_new         ! compacted diameter [mm]
        d_m  = d_mm * 1d-3        ! diameter but using different unit [m]
        xd  = dble(const_pi*d_mm/wave)  ! size parameter [-]
        fmelt = fmelt_new         ! update liquid fraction
     endif 

!
! Gamma function
!
     call gamma_function(gnu,d_m,dn, fgamma)   ! get gamma function
     num = ntot * fgamma                       ! particle number density[1/m4]


     tempd=dble(temp)  !temperature [K] in double precision

     ICE: if(fmelt < 1.0) then !ice or melting particles
!
! complex refractive index of pure ice 
!
       freqd=dble(freq)  !frequency [GHz] in double precision
       call iceoptic(freqd,tempd,eicered,eiceimd)  !get complex refractive index
       !eicere=REAL(eicered)      !real      component of complex refractive index
       !eiceim=REAL(eiceimd)      !imaginary component of complex refractive index
       eicere=dble(eicered)      !real      component of complex refractive index
       eiceim=dble(eiceimd)      !imaginary component of complex refractive index
       eice=dcmplx(eicere,eiceim) !complex refractive index of ice
       eair=dcmplx(1.0006,0.0)    !complex refractive index of air

!
! Calculate effective dielectric constant of frozen particle using different methods/assumptions. 
!
        faa =1.d0 -  (dens/densice)  ! volume fraction of inclusions (air)

        ice_select: select case(ice_refraction_func)
        case(1)  ! air include ice 

          call mg_ellips(1.d0-faa,eair,eice,ei)
          eid=dcmplx(ei)
          crefd_ice=cdsqrt(eid)

        case(2)  ! ice include air

          call mg_ellips(faa,eice,eair,ei)
          eid=dcmplx(ei)
          crefd_ice=cdsqrt(eid)

        case(3)  ! homogeneous mixing of ice and air

          call em_ellips( 1.d0-faa, eair, eice, ei)
          eid=dcmplx(ei)
          crefd_ice=cdsqrt(eid)

        end select ice_select


     endif ICE

     LIQ: if(fmelt > 0.0) then !liquid or melting particles
!
! complex refractive index of liquid water
!
       freqhzd=dble(freq*1.e+9)      !frequency
       sald=0.d0                     !salinity
       call watoptic(freqhzd,tempd,sald,ewatred,ewatimd)  !get complex refractive index 
       ewatd=dcmplx(ewatred,ewatimd) !double precision complex*16
       crefd_liq=cdsqrt(ewatd)       !square root of complex refractive index
       ewat=cmplx(ewatred,ewatimd)   !complex*8

     endif LIQ


!
! choose either liquid, ice, or melting
!
   if(fmelt == 0.) then    ! ice condensates
      crefd = crefd_ice
   elseif(fmelt== 1.) then ! liquid condensates
      crefd = crefd_liq
   else                       ! mixed-phase condensates 

     melt_select: select case(melt_opt)
     case(0)   ! no melting scheme (just use ice dielectric const.)
         crefd = crefd_ice
     case(1)   ! water include ice
        call mg_ellips( 1.-fmelt, ewat, ei, emelt_wi )
        emeltd=dcmplx(emelt_wi)
        crefd = cdsqrt(emeltd)
     case(2)   ! ice include water     
         call mg_ellips( fmelt, ei, ewat, emelt_iw)
         emeltd=dcmplx(emelt_iw)
         crefd = cdsqrt(emeltd)
     case(3)   !averaging above two solutions
         call mg_ellips( 1.-fmelt, ewat, ei, emelt_wi )
         call mg_ellips( fmelt, ei, ewat, emelt_iw)
         eimag = 10.**( log10(aimag(emelt_wi))*(fmelt) + log10(aimag(emelt_iw))*(1.-fmelt)    )
         ereal = 10.**( log10(dble(emelt_wi))*fmelt + log10(dble(emelt_iw))*(1.-fmelt)    )
         emelt =dcmplx(ereal,eimag)
         emeltd=dcmplx(emelt)
         crefd =cdsqrt(emeltd)
     case(4)   ! homogeneous mixing of ice and water
         call em_ellips( 1-fmelt, ewat, ei, emelt_em)
         emeltd=dcmplx(emelt_em)  !ewatd
         crefd =cdsqrt(emeltd)
     end select melt_select

   endif
 
!
! call Mie program to get optical properties
!
     call mie_sphere(xd,crefd,qscad,qextd,asymd,qbscad)

!
!double to single precision
!
     qext = REAL(qextd,kind=r_size)   !extinction efficiency 
     qsca = REAL(qscad,kind=r_size)   !scattering 
     asym = REAl(asymd,kind=r_size)   !asymetry
     qbsca= REAL(qbscad,kind=r_size)  !backscattering

!
! integrate over particle size distribution; (Cross sectin extionction)*(cross sectional aerea)
! toshii (make sure the unit...)
!

     factor = num * const_pi * 0.25 * d_mm * d_mm * d_increment * 1.d-6  

     bext = bext + qext      * factor   ![km-1]
     bsca = bsca + qsca      * factor   ![km-1]
     bsym = bsym + qsca*asym * factor   ![km-1]
     bq11 = bq11 + qbsca     * factor   ![km-1]

  enddo SIZE_LOOP

!
! check for distribution with very small extinction;
! set parameters to zero to avoid numerical problems
!

  if( bext .gt. min_bext) then
      ksca=bext        
      asca=bsca/bext
      gsca=bsym/bsca
      pbck=bq11/bsca
   else
      ksca=0. ; asca=0. ; gsca=0. ; pbck=0.
   end if


 return
 end subroutine mie_rams_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine snow_compaction( rad , den , fmelt , rad_new , fmelt_new )
 implicit none
!--------------------------------------------------------------------------------
! Purpose : Compact snow aggregate radius as a function of melting fraction. 
! 
! Toshi Matsui @ NASA GSFC: initial
!--------------------------------------------------------------------------------
 real(r_size),intent(in) :: rad        ! initial particle radius [mm]
 real(r_size),intent(in) :: den        ! particle density [g/cm3]
 real(r_size),intent(in) :: fmelt      ! 0~1 [-]
 real(r_size),intent(out) :: rad_new   ! compacted particle radius  [mm]
 real(r_size),intent(out) :: fmelt_new ! 0~1 [-]
 real(r_size) :: rad_liq ![mm]
 real(r_size),parameter :: den_liq = 1.0 !liquid density [g/cm3]


  rad_liq = ( (rad**3.)*den*fmelt / den_liq  )**(1./3.)  !radius for liquid water [mm]

  rad_new = ( ( (rad**3.)*den + (rad_liq**3.)*(den-den_liq) )/den  )**(1./3.)  ! [mm]

  fmelt_new = (rad_liq/rad_new)**3  !this is actually volumetric fraction of water within melting ice-air inclusion

  return
 end subroutine snow_compaction

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gamma_reff(x,ga_out)

!---------------------------------------------------------------------------------------------------
! Comments:                             
!   compute the gamma function a(x) for single precision floating point. 
!       input :  x  --- argument of a(x)
!                       ( x is not equal to 0,-1,-2,... )
!       output:  ga --- a(x)
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Adapted to SDSU
!           
! References: 
!----------------------------------------------------------------------------------------------------
 implicit none
 integer :: k, m, m1
 real(r_size),dimension(26) :: g
 real(r_size) :: pi, ga, gr, r, z
 data g/1.0d0,0.5772156649015329d0, &
       -0.6558780715202538d0, -0.420026350340952d-1, &
        0.1665386113822915d0,-.421977345555443d-1, &
        -.96219715278770d-2, .72189432466630d-2, &
        -.11651675918591d-2, -.2152416741149d-3, &
        .1280502823882d-3, -.201348547807d-4, &
        -.12504934821d-5, .11330272320d-5, &
        -.2056338417d-6, .61160950d-8, &
         .50020075d-8, -.11812746d-8, &
        .1043427d-9, .77823d-11, &
        -.36968d-11, .51d-12, &
        -.206d-13, -.54d-14, .14d-14, .1d-15/
 real(r_size) :: x, ga_out
 pi=3.141592653589793d0
 if (x.eq.int(x)) then
     if (x.gt.0.0d0) then
         ga=1.0d0
         m1=x-1
        do k=2,m1
           ga=ga*k
        enddo
! <<< T. Hashino added 2011/02/07
     elseif(x.eq.0.0d0) then
        ga=1.0d0
! >>> T. Hashino added 2011/02/07
     else
        ga=1.0d+300
     endif
  else
     if (dabs(dble(x)).gt.1.0d0) then
         z=dabs(dble(x))
         m=int(z)
         r=1.0d0
        do k=1,m
           r=r*(z-k)
        enddo
        z=z-m
     else
        z=dble(x)
     endif
     gr=g(26)
     do k=25,1,-1
        gr=gr*z+g(k)
     enddo
     ga=1.0d0/(gr*z)
     if (dabs(dble(x)).gt.1.0d0) then
         ga=ga*r
         if (x.lt.0.0d0) ga=-pi/(x*ga*dsin(pi*x))
     endif
  endif

  ga_out = real(ga,kind=r_size)

  return
 end subroutine gamma_reff

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gamma_function(gnu,d,dn, fgamma)
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   compute the gamma function for single-precision floating pint
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Adapted to SDSU
!           
! References: 
!----------------------------------------------------------------------------------------------------

 implicit none
 real(r_size) ,intent(in)  :: gnu    ! PSD shape parameter for generalized gamma distribution
 real(r_size) ,intent(in)  :: d      ! diameter [m] 
 real(r_size) ,intent(in)  :: dn     ! characteristic diameter [m]
 real(r_size) ,intent(out) :: fgamma ! gamma function [1/m]

 real(r_size) :: gfac

 call gamma_reff(gnu,gfac)

 fgamma = (1.0_r_size/gfac) * ( (d/dn)**(gnu-1.0_r_size) ) * (1.0_r_size/dn) * exp(-d/dn)  ![1/m]
  
 return
 end subroutine gamma_function

 end module module_opt_micro
