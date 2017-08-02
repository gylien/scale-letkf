!
!  module common_precip
!
!  module for precipitation assimilation
!
!  created  Jan. 2012, Guo-Yuan Lien, UMD
!  adopted to GFS-LETKF and modified, May 2013, Guo-Yuan Lien, UMD
!  modified, Spetember 2013, Guo-Yuan Lien, UMD
!  adopted to SCALE-LETKF, Aug 1, 2017, Cheng Da, UMD
!
!  function dinvnorm(p) modified from Ren-Raw Chen, 
!    Rutgers University in New Brunswick, New Jersey
!    http://home.online.no/~pjacklam/notes/invnorm/
!
!-------------------------------------------------------------------------------
!
!  subroutine read_ppcdf     (cdffile_m, cdffile_o, ppcdf_m, ppcdf_o, ppzero_m, ppzero_o)
!  function   pptrans_normal (pp, ppcdf, ppzero)
!  function   pptrans_log    (pp)
!  subroutine pptrans_normal_mdzero_def (pp_ens, ppcdf, ppzero,           zero_mem, ym, sigma)
!  function   pptrans_normal_mdzero     (pp,     ppcdf, ppzero, ppzero_m, zero_mem, ym, sigma)
!  function   compact_tail   (pos_cdf)
!  function   dinvnorm       (p)
!
!-------------------------------------------------------------------------------
module common_precip

  use common,       only : r_size, pi
  use common_nml,   only : MEMBER, &
                           use_precip, ncdf, ppzero_thres, gausstail_thres, &
                           opt_pptrans, opt_ppobserr, log_trans_tiny, &
                           min_obserr_rain, obserr_rain_percent, obserr_rain_lt, &
                           obserr_rain_gt


  !use common_scale, only : nlon, nlat, nlonh, nlath, nlong, nlatg
  implicit none

!-------------------------------------------------------------------------------
  public 

  real(r_size),allocatable,save :: ppcdf_m(:,:,:) ! nlon*nlat*(0:ncdf)
  real(r_size),allocatable,save :: ppcdf_o(:,:,:) ! nlon*nlat*(0:ncdf)
  real(r_size),allocatable,save :: ppzero_m(:,:,:) ! nlon*nlat
  real(r_size),allocatable,save :: ppzero_o(:,:,:) ! nlon*nlat

!-------------------------------------------------------------------------------

  !integer, parameter :: pp_bg_nlev = 2
  !integer, parameter :: pp_bg_levs(pp_bg_nlev-1) = &
  !                      (/24/)
  !integer, parameter :: pp_ob_nlev = 2
  !real(r_size), parameter :: pp_ob_levs(pp_ob_nlev-1) = &
  !                      (/ppzero_thres/)
  !logical, parameter :: pp_criterion(pp_bg_nlev,pp_ob_nlev) = reshape((/ &
! !          bg1   , bg2
  !         .false.,.true., &  ! ob1
  !         .false.,.true.  &  ! ob2
  !         /), (/pp_bg_nlev,pp_ob_nlev/))

contains

!-------------------------------------------------------------------------------

!subroutine read_ppcdf (cdffile_m, cdffile_o, ppcdf_m, ppcdf_o, ppzero_m, ppzero_o)
!
!  implicit none
!
!  character(len=*), intent(in) :: cdffile_m
!  character(len=*), intent(in) :: cdffile_o
!  real(r_size), intent(out) :: ppcdf_m(nlon,nlat,0:ncdf)
!  real(r_size), intent(out) :: ppcdf_o(nlon,nlat,0:ncdf)
!  real(r_size), intent(out) :: ppzero_m(nlon,nlat)
!  real(r_size), intent(out) :: ppzero_o(nlon,nlat)
!
!  real(r_sngl) :: ppcdf_ms(nlon,nlat,0:ncdf)
!  real(r_sngl) :: ppzero_ms(nlon,nlat)
!  real(r_sngl) :: ppcdf_os(nlon,nlat,0:ncdf)
!  real(r_sngl) :: ppzero_os(nlon,nlat)
!  integer :: i, j, b, iolen
!  logical :: ex
!
!end subroutine read_ppcdf

!-------------------------------------------------------------------------------
function pptrans_normal (pp, ppcdf, ppzero)

  implicit none

  real(r_size) :: pptrans_normal

  real(r_size), intent(in) :: pp
  real(r_size), intent(in) :: ppcdf(0:ncdf)
  real(r_size), intent(in) :: ppzero

  real(r_size) :: pos_cdf, rr
  integer :: b

  if (ppcdf(0) < -1.0d0 .or. ppzero < -1.0d0) then
    write (*, *) '[Error] Wrong input CDF.'
    stop
  end if

  if (pp < ppzero_thres) then

    if (opt_pptrans == 2) then
      pos_cdf = ppzero * 0.5d0
    else
      write (*, *) '[Error] Unsupported transformation method.'
      stop
    end if

  else ! [pp >= ppzero_thres]

    if (pp < ppcdf(0)) then
      pos_cdf = 0.0d0
    else
      do b = 1, ncdf+1
        if (b > ncdf) then
          pos_cdf = 1.0d0
          exit
        else if (pp < ppcdf(b)) then
          rr = (pp - ppcdf(b-1)) / (ppcdf(b) - ppcdf(b-1))
          pos_cdf = ((1.0d0-rr) * real(b-1,r_size) + rr * real(b,r_size)) / real(ncdf,r_size)
          exit
        end if
      end do
    end if

  end if

  pptrans_normal = dinvnorm(compact_tail(pos_cdf))

end function pptrans_normal

!-------------------------------------------------------------------------------
function pptrans_log (pp)

  implicit none

  real(r_size) :: pptrans_log
  real(r_size), intent(in) :: pp

  if (pp < ppzero_thres) then
    pptrans_log = log(log_trans_tiny)
  else
    pptrans_log = log(pp + log_trans_tiny)
  end if

end function pptrans_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pptrans_normal_mdzero_def (pp_ens, ppcdf, ppzero, zero_mem, ym, sigma)

  implicit none

  real(r_size), intent(inout) :: pp_ens(MEMBER)
  real(r_size), intent(in)    :: ppcdf(0:ncdf)
  real(r_size), intent(in)    :: ppzero
  integer, intent(out)        :: zero_mem
  real(r_size), intent(out)   :: ym
  real(r_size), intent(out)   :: sigma

  real(r_size) :: pos_cdf
  real(r_size) :: ppzero_b, pprain_b
  real(r_size) :: y_trace, y_trace_b
  real(r_size) :: alpha, beta
  logical :: zero(MEMBER)
  integer :: n

!------
!  real(r_size) :: pp_ens_ori(MEMBER)
!  pp_ens_ori = pp_ens
!------

  if (ppcdf(0) < -1.0d0 .or. ppzero < -1.0d0) then
    write (*, *) '[Error] Wrong input CDF.'
    stop
  end if
  if (opt_pptrans /= 3) then
    write (*, *) '[Error] Unsupported transformation method.'
    stop
  end if

  beta = 0.0d0
  zero_mem = 0
  zero = .false.
  do n = 1, MEMBER
    if (pp_ens(n) < ppzero_thres) then
      zero_mem = zero_mem + 1
      zero(n) = .true.
    else
      pp_ens(n) = pptrans_normal(pp_ens(n), ppcdf, ppzero)
      beta = beta + pp_ens(n)
    end if
  end do
  beta = beta / real(MEMBER, r_size)
  ppzero_b = real(zero_mem, r_size) / real(MEMBER, r_size)
  pprain_b = 1.0d0 - ppzero_b

  y_trace = dinvnorm(compact_tail(ppzero))
  y_trace_b = dinvnorm(compact_tail(ppzero_b))

  alpha = 0.0d0 - exp(0.0d0 - 0.5d0*y_trace_b*y_trace_b) / sqrt(2.0d0*pi)
  ym = (alpha * y_trace + beta * y_trace_b) / (alpha + pprain_b * y_trace_b)
  sigma = (pprain_b * y_trace - beta) / (alpha + pprain_b * y_trace_b)

  do n = 1, MEMBER
    if (zero(n)) then
      pos_cdf = ppzero_b * 0.5d0
      pp_ens(n) = ym + sigma * dinvnorm(compact_tail(pos_cdf))
    end if
  end do

!------
!  print *, '----'
!  print *, ppzero, ppzero_b, dinvnorm(pos_cdf)
!  print *, y_trace, ym, sigma
!  do n = 1, MEMBER
!    print *, pp_ens_ori(n), pp_ens(n), zero(n)
!  end do
!------

end subroutine pptrans_normal_mdzero_def

!-------------------------------------------------------------------------------
function pptrans_normal_mdzero (pp, ppcdf, ppzero, ppzero_m, zero_mem, ym, sigma)

  implicit none

  real(r_size) :: pptrans_normal_mdzero

  real(r_size), intent(in) :: pp
  real(r_size), intent(in) :: ppcdf(0:ncdf)
  real(r_size), intent(in) :: ppzero
  real(r_size), intent(in) :: ppzero_m
  integer, intent(in)      :: zero_mem
  real(r_size), intent(in) :: ym
  real(r_size), intent(in) :: sigma

  real(r_size) :: pos_cdf, rr
  integer :: b

  if (ppcdf(0) < -1.0d0 .or. ppzero < -1.0d0) then
    write (*, *) '[Error] Wrong input CDF.'
    stop
  end if
  if (opt_pptrans /= 3) then
    write (*, *) '[Error] Unsupported transformation method.'
    stop
  end if

  if (pp < ppzero_thres) then

    pos_cdf = ppzero * 0.5d0

  else ! [pp >= ppzero_thres]

    if (pp < ppcdf(0)) then
      pos_cdf = 0.0d0
    else
      do b = 1, ncdf+1
        if (b > ncdf) then
          pos_cdf = 1.0d0
          exit
        else if (pp < ppcdf(b)) then
          rr = (pp - ppcdf(b-1)) / (ppcdf(b) - ppcdf(b-1))
          pos_cdf = ((1.0d0-rr) * real(b-1,r_size) + rr * real(b,r_size)) / real(ncdf,r_size)
          exit
        end if
      end do
    end if

  end if

!------
!  print *, '---###'
!  print *, ppzero, ppzero_m, zero_mem, ym, sigma
!  print *, pos_cdf
!------

  if (pos_cdf < ppzero_m) then
    pos_cdf = (pos_cdf / ppzero_m) * (real(zero_mem, r_size) / real(MEMBER, r_size))
    pptrans_normal_mdzero = ym + sigma * dinvnorm(compact_tail(pos_cdf))
  else
    pptrans_normal_mdzero = dinvnorm(compact_tail(pos_cdf))
  end if

!------
!  print *, pos_cdf
!  print *, pptrans_normal_mdzero
!------

end function pptrans_normal_mdzero

!-------------------------------------------------------------------------------
function compact_tail (pos_cdf)

  implicit none

  real(r_size) :: compact_tail
  real(r_size), intent(in) :: pos_cdf

  compact_tail = pos_cdf
  if (compact_tail < gausstail_thres        ) compact_tail = gausstail_thres
  if (compact_tail > 1.0d0 - gausstail_thres) compact_tail = 1.0d0 - gausstail_thres

end function compact_tail

!-------------------------------------------------------------------------------
! ren-raw chen, rutgers business school
! normal inverse
! translate from http://home.online.no/~pjacklam/notes/invnorm
! a routine written by john herrero
!-------------------------------------------------------------------------------
! ren-raw chen, rutgers business school
real*8 function dinvnorm(p)
      real*8 p,p_low,p_high
      real*8 a1,a2,a3,a4,a5,a6
      real*8 b1,b2,b3,b4,b5
      real*8 c1,c2,c3,c4,c5,c6
      real*8 d1,d2,d3,d4
      real*8 z,q,r
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=dsqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=dsqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   dinvnorm=z
      return
end function dinvnorm


end module common_precip

