  module module_simulator
    use common, only: r_size
    use common_geosatpr
    implicit none

  contains

  !>---------------------------------------------------------------------
  !moved to Trans_XtoY_GPR subroutine get_others
  !moved to Trans_XtoY_GPR   implicit none
  !moved to Trans_XtoY_GPR   
  !moved to Trans_XtoY_GPR   integer :: i,j,k,ic

  !moved to Trans_XtoY_GPR   do k = 1, ke ; do j = js, je ; do i = is, ie

  !moved to Trans_XtoY_GPR     ! specific humidity 
  !moved to Trans_XtoY_GPR     ! exner function
  !moved to Trans_XtoY_GPR     atmos(i,j,k)%exner = (atmos(i,j,k)%press * 0.001d0)**(const_Rd/1004.d0)

  !moved to Trans_XtoY_GPR     ! layer thickness [km]
  !moved to Trans_XtoY_GPR     atmos(i,j,k)%dhgt = atmos_stag(i,j,k)%hgt-atmos_stag(i,j,k-1)%hgt

  !moved to Trans_XtoY_GPR     ! layer height [km]
  !moved to Trans_XtoY_GPR     atmos(i,j,k)%hgt = 0.5d0 * (atmos_stag(i,j,k)%hgt+atmos_stag(i,j,k-1)%hgt)

  !moved to Trans_XtoY_GPR     ! vertical velocity [m/s]
  !moved to Trans_XtoY_GPR     atmos(i,j,k)%omega = 0.5d0 * (atmos_stag(i,j,k)%omega+atmos_stag(i,j,k-1)%omega)

  !moved to Trans_XtoY_GPR   enddo ; enddo ; enddo

  !moved to Trans_XtoY_GPR   do j = js, je ; do i = is, ie
  !moved to Trans_XtoY_GPR   
  !moved to Trans_XtoY_GPR     ! column-integrated condensate [kg/m2]
  !moved to Trans_XtoY_GPR     qcol_gce(i,j)%cloud   = sum(q_gce(i,j,1:ke)%cloud   * atmos(i,j,1:ke)%dhgt)
  !moved to Trans_XtoY_GPR     qcol_gce(i,j)%rain    = sum(q_gce(i,j,1:ke)%rain    * atmos(i,j,1:ke)%dhgt)
  !moved to Trans_XtoY_GPR     qcol_gce(i,j)%ice     = sum(q_gce(i,j,1:ke)%ice     * atmos(i,j,1:ke)%dhgt)
  !moved to Trans_XtoY_GPR     qcol_gce(i,j)%snow    = sum(q_gce(i,j,1:ke)%snow    * atmos(i,j,1:ke)%dhgt)
  !moved to Trans_XtoY_GPR     qcol_gce(i,j)%graupel = sum(q_gce(i,j,1:ke)%graupel * atmos(i,j,1:ke)%dhgt)
  !moved to Trans_XtoY_GPR     qcol_gce(i,j)%hail    = sum(q_gce(i,j,1:ke)%hail    * atmos(i,j,1:ke)%dhgt)

  !moved to Trans_XtoY_GPR   enddo ; enddo

  !moved to Trans_XtoY_GPR   return

  !moved to Trans_XtoY_GPR end subroutine get_others

  !>---------------------------------------------------------------------
  subroutine bulk_dsd
    implicit none
    
    real(r_size) :: re
    real(r_size) :: temp,rc

    n0_gce%cloud = 50.0d+6   ! number concentration assumed for cloud [1/m3] (=50[1/cm3]) for maritime
    n0_gce%rain = 0.08d+8   ! intercept for rain [1/m4] (=0.08[1/cm4])
    n0_gce%snow = 0.03d+8   ! intercept for snow [1/m4] (=0.03[1/cm4])
    n0_gce%graupel= 0.04d+8   ! intercept for graupel [1/m4] (=0.04[1/cm4])
    n0_gce%hail = 0.002d+8  ! does not exist

    rho_gce%cloud = 1.0d+3  ! density of cloud [kg/m3] (=1.0[g/cm3])
    rho_gce%rain = 1.0d+3  ! density of rain [kg/m3] (=1.0[g/cm3])
    rho_gce%ice = 0.91668d+3  ! density of ice [kg/m3] (=0.91668[g/cm3])
    rho_gce%snow = 0.1d+3  ! density of snow [kg/m3] (=0.1[g/cm3]) 
    rho_gce%graupel= 0.4d+3  ! density of graupel [kg/m3] (=0.4[g/cm3])
    rho_gce%hail = 0.9d+3  ! does not exist

    call re_LUT_Heymsfield_Platt_1984('init',270.0_r_size,1.0_r_size, re)

    return

  end subroutine bulk_dsd


  !>---------------------------------------------------------------------
! MOVE TO Trans_XtoY_GPR in common/common_obs_scale.f90
!
!  ! compute effective radius of condensates
!  subroutine re_all
!    implicit none
!    integer :: i,j,k,ic
!    real(r_size) :: n0 ! intercept [1/m**4]
!
!    do k=ks,ke ; do j=js,je ; do i=is,ie
!      call re_bulk_mono(0, rho_gce%cloud  , n0_gce%cloud ,q_gce(i,j,k)%cloud , re_gce(i,j,k)%cloud )
!!tmp      call re_bulk_mono(1, rho_gce%ice    , n0_gce%ice, q_gce(i,j,k)%ice , re_gce(i,j,k)%ice )
!      call re_LUT_Heymsfield_Platt_1984('proc', atmos(i,j,k)%t_air, q_gce(i,j,k)%ice, re_gce(i,j,k)%ice )
!      call re_bulk_exp(rho_gce%rain  , n0_gce%rain  , q_gce(i,j,k)%rain , re_gce(i,j,k)%rain)
!      call re_bulk_exp(rho_gce%snow  , n0_gce%snow  , q_gce(i,j,k)%snow , re_gce(i,j,k)%snow)
!      call re_bulk_exp(rho_gce%graupel , n0_gce%graupel , q_gce(i,j,k)%graupel , re_gce(i,j,k)%graupel)
!      call re_bulk_exp(rho_gce%hail  , n0_gce%hail  , q_gce(i,j,k)%hail , re_gce(i,j,k)%hail)
!    enddo ; enddo ; enddo
!    
!    return
!
!  end subroutine re_all

  !>---------------------------------------------------------------------
  subroutine re_LUT_Heymsfield_Platt_1984(proc,temp,wc,re)
    implicit none

    character(len=4),intent(in) :: proc !'init' or 'proc'
    real(r_size),intent(in)  :: temp  !temperature [K]
    real(r_size)  :: wc    ! water content [g/m3]
    real(r_size),intent(out) :: re    !effective radius [micron]

    integer :: t ! loop indice
    integer,parameter :: tmin = 210 , tmax = 260
    real(r_size),save :: re_lut(tmin:tmax)  !LUT of re as a function of temperature [micron]
    real(r_size) :: wgt1,wgt2 !weighting function

    proc_select: select case(proc)
    case('init')
      do t = tmin, tmax
        call re_Heymsfield_Platt_1984(REAL(t,kind=r_size),1.0_r_size, re)
        re_lut(t) = re
      enddo !t

    case('proc')
      if( wc < q_min_condensate ) then
        re = 0.d0 ; wc = 0.d0;
        return
      endif

    end select proc_select

    return

  end subroutine

  !>---------------------------------------------------------------------
  ! Compute drop effective radius from observed DSD in Heymsfield Platt 1984]  
  ! Heymsfield, A.J., and C. Platt, 1984: A Parameterization of the Particle Size Spectrum of Ice 
  ! Clouds in Terms of the Ambient Temperature and the Ice Water Content. J. Atmos. Sci., 41, 846-855.
  subroutine re_Heymsfield_Platt_1984(temp,lwc,re)
    implicit none
    real(r_size),intent(in)  :: temp  !temperature [K]
    real(r_size),intent(in)  :: lwc   !liquid water content [g/m3]
    real(r_size),intent(out) :: re    !effective radius [micron]
    
    integer :: imax
    integer :: i
    integer :: nopt
    real(r_size) :: rad
    real(r_size) :: densice
    real(r_size) :: densliq
    real(r_size) :: density
    real(r_size) :: densi
    real(r_size) :: iwctest
    real(r_size) :: norm
    real(r_size) :: num
    real(r_size) :: faa
    real(r_size) :: dr
    
    data densliq/1.0d+3/
    data densice/0.917d+3/
    real(r_size) :: third_moment  ![m3]  
    real(r_size) :: sec_moment     ![m2]
    
    third_moment = 0.  ![m3]  
    sec_moment    =0.  ![m2]
!
!     Begin by checking if hydrometeors of this species are present.
!     If not, set scattering parameters to zero and return.
!
    if(lwc .lt. q_min_condensate) then
       return
    endif

!
!     Loop over particle sizes:

!     increments of particle radius are 0.005 mm; the particle
!     size distribution is expressed as a particle number density,
!     num, per radius increment.  This distribution is given
!     by the fit to observed cloud ice distributions by
!     Heymsfield and Platt (1984)
!     two options are available:
!     nopt=0:   ice particle mass distributed in spherical
!               volume with diameter equal to maximum crystal
!               dimension (l).
!     nopt=1:   ice particle described as pure ice sphere
!               with same mass as elongated crystal.
    nopt=1
    dr = 0.005d0
    imax = nint(2.5d0/dr)
!
!     first compute normalization factor for particle size distribution
    norm=1.
    iwctest=0.

    do i=0,imax
      rad=dr*0.5d0+dr*float(i)
      call heymplatt(nopt,temp,lwc,rad,densi,num,norm)  !num [1/m**4]
      density=densi
      iwctest=iwctest+ &
      num*(1.d-3)*density*4.d0*const_pi*((rad*1.d0-1.d0)**3.)/3.d0*.005d0
    enddo
    
    norm=iwctest/lwc
    
    do i=0,imax
      rad=dr*0.5d0+dr*float(i)
    
      call heymplatt(nopt,temp,lwc,rad,densi,num,norm) !num [1/m**4]
    
!
! compute 2nd 3rd moment
!
      third_moment = third_moment + ( rad**3) * num *dr ! [mm3]  
      sec_moment   = sec_moment   + ( rad**2) * num *dr ! [mm2]
    end do

!
! compute effective radius
!
    re = third_moment / sec_moment *1d+3

    return
  end subroutine re_Heymsfield_Platt_1984

  !>---------------------------------------------------------------------
  ! Return particle density and number density
  ! of ice crystals, given temperature, ice water content of
  ! ice crystal distribution, and maximum particle
  ! dimension.  Follows the empirical relation of
  ! Heymsfield and Platt (1984).
  !
  ! Input
  ! Two options are available:
  ! nopt=0:   ice particle mass distributed in spherical
  !             volume with diameter equal to maximum crystal
  !            dimension (l).
  ! nopt=1:      ice particle described as pure ice sphere
  !            with same mass as elongated crystal.
  ! t            temperature [K]
  ! iwc      equivalent water content of cloud ice [g/m**3]
  ! r            particle radius [mm]
  ! norm     normalization factor required to contrain
  !            integrated distribution to equal the equivalent
  !            ice water content.  heymplatt must be called
  !            first with norm set equal to 1 to give
  !           proper normalization
  !
  ! Output
  ! den      particle density [g/cm**3]
  ! numden   particle number density [1/m**4]
  !
  ! References: 
  !   Heymsfield, A.J., and C. Platt, 1984: A Parameterization of the Particle Size Spectrum of Ice 
  !     Clouds in Terms of the Ambient Temperature and the Ice Water Content. J. Atmos. Sci., 41, 846-855.
  subroutine heymplatt(nopt,t,iwc,r,den,numden,norm)
    implicit none
    
    integer :: n
    integer :: nopt
    real(r_size):: t
    real(r_size):: iwc
    real(r_size):: r
    real(r_size):: den
    real(r_size):: numden
    real(r_size):: a(0:4)
    real(r_size):: b(0:4)
    real(r_size):: c(0:4)
    real(r_size):: d(0:4)
    real(r_size):: e(0:4)
    real(r_size):: suma
    real(r_size):: sumb
    real(r_size):: sumc
    real(r_size):: sumd
    real(r_size):: sume
    real(r_size):: tc
    real(r_size):: mass
    real(r_size):: l
    real(r_size):: l0
    real(r_size):: b1
    real(r_size):: b2
    real(r_size):: a1
    real(r_size):: a2
    real(r_size):: n100diwc
    real(r_size):: n1000diwc
    real(r_size):: lc
    real(r_size):: dc
    real(r_size):: densice
    real(r_size):: norm
    real(r_size):: dldr
    
    data a/-1.1430d+1,-7.3892d-1,-1.8647d-2,-1.4045d-4,0.0d+0/
    data b/1.8940d+1,2.7658d0,1.2833d-1,2.7750d-3,2.2994d-5/
    data c/-1.0159d+1,-1.4538d0,-1.3511d-2,1.1318d-3,2.2360d-5/
    data d/1.6764d+1,-1.5072d-1,-1.9713d-2,-3.5051d-4,-1.6727d-6/
    data e/1.5508d+2,1.8377d+1,8.5312d-1,1.6879d-2,1.1873d-4/
    data densice/0.917/
    
    
! Temperature [C]
  tc=t-const_Kel2Cel

! Note: empirical formulae only apply to range
! -60 < tc < -20
  if(tc .gt. -20.) tc=-20.
  if(tc .lt. -60.) tc=-60.

! Calculate maximum particle dimension [microns]
  if(nopt .eq. 0) then
     l=2.*r*1.d+3
  else
     mass=densice*4.*const_pi*((r*1.d-1)**3.)/3.
    if(2.*r .le. 0.3) then
       l=((mass/(densice*(3.*sqrt(3.)/2.)* &
         ((0.5/2.)**2.)*1.d-3))**.33333)* &
        1.d+3
       dldr=(4./3.)*(1.d-3)*densice*const_pi*(mass**(-.66666))*r*r/ &
            ((densice*(3.*sqrt(3.)/2.)*((0.5/2.)**2.)*1.d-3)**.33333)
    else
       l=((mass/(densice*(3.*sqrt(3.)/2.)*((.2/2.)**2.)*1.d-3)) &
         **.55249)*1.d+3
       dldr=(4./1.82)*(1.d-3)*densice*const_pi*(mass**(1./1.82-1.))*r*r/ &
            ((densice*(3.*sqrt(3.)/2.)*((0.2/2.)**2.)*1.d-3)**(1./1.82))
    endif
 endif

 suma=0.
 sumb=0.
 sumc=0.
 sumd=0.
 sume=0.
 do n=0,4
    suma=suma+a(n)*(tc**n)
    sumb=sumb+b(n)*(tc**n)
    sumc=sumc+c(n)*(tc**n)
    sumd=sumd+d(n)*(tc**n)
    sume=sume+e(n)*(tc**n)
 enddo
 b1=suma
 b2=sumb

!  Liou's fit of b2 fails at low temperature;
!  since Heymsfield and Platt (1984) data indicate
!  a nearly constant value of -4., we use it here.
 b2=-4.

 if(tc .ge. -37.5) then
    n100diwc=exp(sumc)
 else
    n100diwc=exp(sumd)
 endif
 n1000diwc=sume

 a1=n100diwc/(100.**b1)
 a2=n1000diwc/(1000.**b2)
 l0=(a2/a1)**(1./(b1-b2))

 if(nopt .eq. 0) then
    if(l .le. l0) then
       numden=(1.d+6)*2.*a1*(l**b1)*iwc/norm
    else
       numden=(1.d+6)*2.*a2*(l**b2)*iwc/norm
    endif
 else
    if(l .le. l0) then
       numden=(1.d+6)*dldr*a1*(l**b1)*iwc/norm
    else
       numden=(1.d+6)*dldr*a2*(l**b2)*iwc/norm
    endif
 endif

! Compute particle density
! Assume randomly oriented hexagonal column with
! width/length relationship determined
! from Heymsfield's empirical relation

 if(nopt .eq. 0) then
    if(l*1.d-3 .le. 0.3) then
       dc=0.5*(l*1.d-3)*1.d-1
       lc=l*1.d-4
    else
       dc=0.2*((l*1.d-3)**0.41)*1.d-1
       lc=l*1.d-4
    endif
       mass=densice*(3.*sqrt(3.)/2.)*((dc/2.)**2.)*lc
       den=mass/(4.*const_pi*((r*1.d-1)**3.)/3.)
 else
    den=0.917
 end if

 return
 end subroutine heymplatt

  !>---------------------------------------------------------------------
  ! Compute drop effective radius for one-moment bulk scheme using analytic solution (fast). 
  subroutine re_bulk_mono(switch, den, n0,  wc, re )
    implicit none

    integer,intent(in) :: switch  ! switch for diagnosis
    real(r_size),intent(in) :: den  !density [kg/m3]
    real(r_size),intent(inout) :: n0   !intercept of exponential DSD [1/m4]
    real(r_size),intent(inout) :: wc   !liquid/ice water content [g/m3] 
    real(r_size),intent(out) :: re  !effective radius [micron] 
    
    real(r_size), parameter :: amice=4.19d-13 ! mass of one cloud ice [kg]

    ! for no particles. 
    if( wc < q_min_condensate ) then
      re = 0.d0 ; wc = 0.d0
    else
      ! compute drop effective radius for mono disperse distribution
      select case(switch)
      case(0)
        ! here n0 is actually the total number concentration of the specie [1/m3]
        re = 1.0d+6*(1.0d-3*wc/n0/den/(4.0d0/3.0d0*const_pi))**(1.0/3.0) ![micron] <- [m]
      case(1)
        re = 1.0d+6*(amice/den/(4.0d0/3.0d0*const_pi))**(1.0/3.0) ![micron] <- [m]
        n0 = wc*1.0d-3/amice ! [1/m3]
      end select
    endif

    return

  end subroutine re_bulk_mono

  !>---------------------------------------------------------------------
  ! Compute drop effective radius for one-moment bulk scheme using analytic solution (fast). 
  subroutine re_bulk_exp(den, n0 , wc, re)
    implicit none

    real(r_size),intent(in) :: den  !density [kg/m3]
    real(r_size),intent(in) :: n0   !intercept of exponential DSD [1/m4]
    real(r_size),intent(inout) :: wc   !liquid/ice water content [g/m3] 
    real(r_size),intent(out) :: re  !effective radius [micron] 
    
    real(r_size) :: lam   !intercept [1/m] ! 2010/09 T. Hashino added 

    ! for no particles. 
    if( wc < q_min_condensate ) then
      re = 0.d0 ; wc = 0.d0; lam = 0.0d0
    else
      ! compute drop effective radius for exponential distribution N(D) = N0*exp(-lam*D)
      lam = (n0*const_pi*den/(wc*(1.d-3)))**(0.25d0)  ![1/m]
      re = 3.0d0/2.0d0/lam *1.0d+6  ![micron] <- [m] ! 2011/02/21 T. Hashino bug found and fixed.
    endif

    return
  end subroutine re_bulk_exp


  end module module_simulator
