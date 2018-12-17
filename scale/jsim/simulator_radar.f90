 subroutine simulator_radar
 use common, only: r_size
 use common_geosatpr
 use module_opt_micro
 implicit none

!--------------------------------------------------------------------------------------------------
!
!               Joint-Simulator (Joint Simulator for Satellite Sensors)
!
! The ‘Joint-Simulator’ is developed by Japan Aerospace Exploration Agency (JAXA) EarthCARE mission 
! and Joint-Simulator team, however JAXA and Joint-Simulator team are not responsible for the safety of 
! the Joint-Simulator nor the reliability of information provided by the Joint-Simulator.  User(s) shall 
! not make any claim against JAXA and Joint-Simulator team with respect to any damages or loss of 
! any kind that may be caused by the use of the Joint-Simulator.
!
!
!
!
!                              TERMS AND CONDITIONS
!
! - Do not distribute the Joint-Simulator to any third party.
! - Please report any bugs you have found.
! - Please give us a feedback on the usability of the Joint-Simulator so as to improve it.
! - Please let us know the achievements you made with the Joint-Simulator, so that we can make a 
!   list on the Joint-Simulator website (https://sites.google.com/site/jointsimulator/home).
! - When you present or publish the results obtained with the Joint-Simulator, please cite the logo 
!   of the Joint-Simulator (http://www.eorc.jaxa.jp/EARTHCARE/Joint-Sim_logo.png) on your presentations, 
!   acknowledge in your papers as follows: "The Joint-Simulator is developed by Japan Aerospace 
!   Exploration Agency (JAXA) EarthCARE mission and Joint-Simulator team", and cite the publications 
!   describing the sensor simulators you used. The relevant publications for a sensor can be 
!   found in the codes or on the Joint-Simulator website.
! - If a Joint-Simulator team member is heavily involved with your project, please acknowledge or 
!   consider adding him or her as a co-author.
!
!
!                             Enjoy using Joint-Simulator!
!
!--------------------------------------------------------------------------------------------
!              = Goddard Satellite Data Simulator Unit =
!
!
! NASA GSFC makes no representations about the suitability of software for any purpose. 
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US 
! government) nor Principal Developers (their organizations) shall be liable for any 
! damages suffered by the user of this software. In addition, please do not distribute 
! the software to third party.
!
!
! Comments: 
!  The methodology to simulate radar echoes and path-integrated attenuation
!  is described in the Appendix of Masunaga and Kummerow (2005). The geophysical
!  parameters that affect radar simulation are summarized in the following.
!
!
! History:
!  09/2010  Tempei Hashino@Univ. of Tokyo AORI ; 
!               - beam convolution is implemented for each process, so no need to gather
!                 all data into the root process.
!  08/2008  Toshi Matsui@NASA GSFC ; Remove optical properties part into module_opt_micowave
!  07/2007  Toshi Matsui@NASA GSFC ; Add dynamic allocation and better-memory coding.  
!  03/2007  Toshi Matsui@NASA GSFC ; Add options for attenuation and ground-based radar options. 
!  03/2007  Toshi Matsui@NASA GSFC ; Adapted for Goddard SDSU.
!
! References:
!  Masunaga, H., and C.D. Kummerow, 2005: Combined Radar and Radiometer Analysis of 
!    Precipitation Profiles for a Parametric Retrieval Algorithm. J. Atmos. Oceanic 
!    Technol., 22, 909-929.
!
!-----------------------------------------------------------------------------------------------------
 integer :: ii, jj, i,j,k,nf ! looping indice
 integer :: ierr             ! alloc statistics

!for FOV convolution
 integer :: n_ct, n_dt  !number of pixels from the center points
 integer :: i_dt, j_ct  !loop indices
 real(r_size) :: xa,ya       !physical distance from the center points [km]
 real(r_size),parameter :: lne = 0.43429448d0 ! = log10(2.71828183)

!for radar constant
 real(8) :: freqhzd, epsreal, epsimag  !dielectric constant parameters
 real(r_size) :: kradar ! radar constant |K^2| 

!
! 1D optical properties
!
 real(r_size),dimension(nln) :: &
       kext1D,sback1D
!
! 1D atmospheric variables
!
 real(r_size),dimension(nln) :: &
       dh1D


!
! sensor angle (umu is global parameter)
!
 umu = cos(view_angle_radar*const_degrad) !sensor view angle

 nf = mxfreq_radar

!
! Derive radar constant for different frequency, if not known.  
!
     freqhzd=dble(freq_radar*1.d+9)
     call watoptic(freqhzd, dble(273.16d0), dble(0.d0), epsreal, epsimag)
     kradar = ((epsreal-1.d0)**2 + epsimag**2)/ ((epsreal+2.d0)**2 + epsimag**2)

     if( k2 /= undefined ) kradar = k2

!move to Trans_XtoY_GPR !
!move to Trans_XtoY_GPR ! Estimate radar echo (Only nadir profiles are currently simulated. -> now 1D driver)
!move to Trans_XtoY_GPR !
!move to Trans_XtoY_GPR      do j = js, je ; do i = is, ie
!move to Trans_XtoY_GPR 
!move to Trans_XtoY_GPR         kext1D(1:ke)=kexttot(i,j,1:ke,nf)
!move to Trans_XtoY_GPR         sback1D(1:ke)=sbacktot(i,j,1:ke,nf)
!move to Trans_XtoY_GPR         dh1D(1:ke)=atmos(i,j,1:ke)%dhgt
!move to Trans_XtoY_GPR 
!move to Trans_XtoY_GPR         call radar_echo( nln, dh1D(1:nln), kext1D(1:nln), sback1D(1:nln), &
!move to Trans_XtoY_GPR                          2.997925e2/freq_radar, kradar, z_out(i,j,1:nln,nf) )
!move to Trans_XtoY_GPR !
!move to Trans_XtoY_GPR ! shift the data back over the topography
!move to Trans_XtoY_GPR !
!move to Trans_XtoY_GPR         do k=ke-(surface(i,j)%k1layer-1),1,-1
!move to Trans_XtoY_GPR            z_out(i,j,k+surface(i,j)%k1layer-1,nf)=z_out(i,j,k,nf)
!move to Trans_XtoY_GPR         enddo
!move to Trans_XtoY_GPR         do k=1,surface(i,j)%k1layer-1
!move to Trans_XtoY_GPR            z_out(i,j,k,nf)=undefined
!move to Trans_XtoY_GPR         enddo
!move to Trans_XtoY_GPR 
!move to Trans_XtoY_GPR      enddo ; enddo

 return
 end subroutine simulator_radar 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine radar_echo( nln, dhgt, kexttot, sbacktot, lambda, &
                         kradar, Z)
  use common_geosatpr, only : r_size, undefined, attenuation, const_pi
  use module_simulator, only : umu
  implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Simulating attenuating (or non-attenuating) radar echo for either ground-based or space-based. 
!  This does not account multiple scattering. 
! 
! History:
! Toshi Matsui@NASA GSFC :: Add kradar as input Dec 2007
! Toshi Matsui@NASA GSFC :: Add ground-based path Jun 2007
! Toshi Matsui@NASA GSFC :: Adapted for Goddard SDSU.
!           
! References: 
!-----------------------------------------------------------------------------------------------------
  integer,intent(in) :: nln ! number of layers
  real(r_size),intent(in) :: kradar          !radar constant |k^2|
  real(r_size),intent(in) :: dhgt(1:nln)   !layer thickness [km] 
  real(r_size),intent(in) :: kexttot(nln)  !total extinction
  real(r_size),intent(in) :: sbacktot(nln) !total backscatter [1/km]
  real(r_size),intent(in) :: lambda !Radar echo wavelength [mm]
  real(r_size),intent(out) :: Z(nln)      !radar reflectivity
  real(r_size) :: pia    !path-integrated attenuation
  real(r_size) :: atten  !attenuation (z=surface to TOA) : attenuation measured from TOA

  integer :: k 
  real(r_size) :: tauext(0:nln+1)  !(z=surface to TOA): Path-integrated Extinction optical thickness
  real(r_size) :: dtau  ! optical depth for each layer
  real(r_size),parameter :: lne = 0.43429448d0 ! = log10(2.71828183)
  real(r_size) :: dbz ![dBZ]
!
  tauext = 0.

!<<<okzk
    do k = nln, 1, -1 !space-born satellite
       dtau = kexttot(k)* dhgt(k)/umu
       tauext(k-1) = tauext(k)  + dtau
       if (sbacktot(k).le.0.d0) then
           dbz = undefined
             Z(k) = 0.d0
       else
!
! Units : dbz in [mm6/m3], lambda in [mm], sbacktot in [km^-1]
!
           if(attenuation) then  ! compute attenuating radar reflectivity
                     dbz = 10.d0*(4.d0*log10(lambda) + &
                                     log10(sbacktot(k)) - &
                                     log10(kradar) - 5.d0*log10(const_pi) - &
                                     2.d0*tauext(k)*lne + 3.d0 + &
                                     log10(1.d0-exp(-2.d0*dtau)) - &
                                     log10(2.d0*dtau))
            else !compute non-attenuating radar reflectivity
                     dbz = 10.d0*(4.d0*log10(lambda) + &
                                     log10(sbacktot(k)) - &
                                     log10(kradar) - 5.d0*log10(const_pi) + 3.d0)
            endif !attenuation
            Z(k) = 10.d0**(dbz/10.d0)  !back to Z
        endif !sbacktot > 0
     enddo !k
!
! Compute PIA and total attentuation
!
     pia = 20.d0*tauext(0)*lne
     atten = tauext(0)


!!!       if (sbacktot(k).le.0.e0) then
!!!           dbz = undefined
!!!             Z(k) = 0.e0
!!!       else
!!!!
!!!! Units : dbz in [mm6/m3], lambda in [mm], sbacktot in [km^-1]
!!!!
!!!!           if(attenuation) then  ! compute attenuating radar reflectivity
!!!!                     dbz = 10.e0*(4.e0*log10(lambda) + &
!!!!                                     log10(sbacktot(k)) - &
!!!!                                     log10(kradar) - 5.e0*log10(const_pi) - &
!!!!                                     2.e0*tauext(k)*lne + 3.e0 + &
!!!!                                     log10(1.e0-exp(-2.e0*dtau)) - &
!!!!                                     log10(2.e0*dtau))
!!!!            else !compute non-attenuating radar reflectivity
!!!                     dbz = 10.e0*(4.e0*log10(lambda) + &
!!!                                     log10(sbacktot(k)) - &
!!!                                     log10(kradar) - 5.e0*log10(const_pi) + 3.e0)
!!!!            endif !attenuation
!!!            Z(k) = 10.e0**(dbz/10.e0)  !back to Z
!!!        endif !sbacktot > 0
!!!     enddo !k

!>>>okzk



 return
 end subroutine radar_echo

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

