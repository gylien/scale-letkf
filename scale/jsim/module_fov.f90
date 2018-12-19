  module module_fov
    use common
    use common_nml
    use common_scale
    use common_geosatpr

    implicit none

    contains

  !>---------------------------------------------------------------------
  !> beam convolution for PR on geostationary satellite
  !subroutine radar_fov
  !  implicit none


  !  ! Interpolate
  !  select case(beamcnv_switch)
  !  case(1)
  !  case(2)
  !  case default
  !    write(*,*) 'no interpolation'
  !    isv=is ; iev=ie ; jsv=js ; jev=je ; ksv=ks ; kev=ke
  !    dbz_out(is:ie,js:je,ks:ke,1) = z_out(is:ie,js:je,ks:ke,1)
  !  end select
  !    

  !  ! convert output from Z to dBZ
  !  do k = ksv, kev ; do j = js, je ; do i = is, ie
  !    if( dbz_out(i,j,k,1).eq.undefined ) then
  !      dbz_out(i,j,k,1) = undefined
  !    elseif( dbz_out(i,j,k,1).eq.0.d0 ) then
  !      dbz_out(i,j,k,1) = undefined
  !    else
  !      dbz_out(i,j,k,1) = 10.d0 * log10( dbz_out(i,j,k,1) ) !back to dBZ
  !    endif
  !    if( dbz_out(i,j,k,1).le.min_echo ) dbz_out(i,j,k,1) = undefined
  !  enddo ; enddo; enddo

  !  ! Undef the ground

  !end subroutine radar_fov

  !>---------------------------------------------------------------------
  !> convert model coordinate --> satellite coordinate
  !> reference: Coordination Group for Meteorological Satellites,
  !>            LRIT/HRIT Global Specification
  subroutine phys_scale2gpr(lon_e,lat_e,lev_e,lon_s,lat_s,lev_s)
    implicit none

    real(r_size),intent(out) :: lon_s ! satellite coordinate (rad)
    real(r_size),intent(out) :: lat_s ! satellite coordinate (rad)
    real(r_size),intent(out) :: lev_s ! satellite coordinate (m)
    real(r_size),intent(in) :: lon_e ! earth coordinate (rad)
    real(r_size),intent(in) :: lat_e ! earth coordinate (rad)
    real(r_size),intent(in) :: lev_e ! earth coordinate (m)

    ! internal work
    real(r_size) :: r1,r2,r3

    r1 = (sat_lev_e+Rearth)-(Rearth+lev_e)*cos(lat_e)*cos(lon_e-sat_lon_e)
    r2 =                    (Rearth+lev_e)*cos(lat_e)*sin(lon_e-sat_lon_e)
    r3 =                    (Rearth+lev_e)*sin(lat_e)

    lev_s = sqrt(r1**2+r2**2+r3**2)
    lon_s = atan(r2/r1)
    lat_s = asin(r3/lev_s)

    return
  end subroutine phys_scale2gpr
  !----------------------------------------------------------------------
  !> convert satellite coordinate --> model coordinate
  subroutine phys_gpr2scale(lon_s,lat_s,lev_s,lon_e,lat_e,lev_e)
    implicit none

    real(r_size),intent(out) :: lon_e ! earth coordinate (rad)
    real(r_size),intent(out) :: lat_e ! earth coordinate (rad)
    real(r_size),intent(out) :: lev_e ! earth coordinate (m)
    real(r_size),intent(in) :: lon_s ! satellite coordinate (rad)
    real(r_size),intent(in) :: lat_s ! satellite coordinate (rad)
    real(r_size),intent(in) :: lev_s ! satellite coordinate (m)

    ! internal work
    real(r_size) :: r1,r2,r3

    r1 = (sat_lev_e+Rearth)-lev_s*cos(lat_s)*cos(lon_s)
    r2 =                    lev_s*cos(lat_s)*sin(lon_s)
    r3 =                    lev_s*sin(lat_s)

    lev_e = sqrt(r1**2+r2**2+r3**2)
    lon_e = atan(r2/r1)+sat_lon_e
    lat_e = asin(r3/lev_e)
    lev_e = lev_e-Rearth

    return
  end subroutine phys_gpr2scale


  end module module_fov
