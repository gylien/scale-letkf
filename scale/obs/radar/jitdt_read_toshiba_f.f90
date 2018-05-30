module jitdt_read_toshiba_f
  use iso_c_binding
  use read_toshiba_f

  interface
     integer(kind=c_int) function jitdt_read_toshiba_c(n_type, jitdt_place, hd, az, el, rtdat) bind(C, name="jitdt_read_toshiba")
       use iso_c_binding
       import c_pawr_header
       import RDIM, AZDIM, ELDIM

       integer(kind=c_int), value :: n_type
       character(kind=c_char) :: jitdt_place(*)
       type(c_pawr_header) :: hd(n_type)
       real(kind=c_float) :: az(AZDIM, ELDIM, n_type)
       real(kind=c_float) :: el(AZDIM, ELDIM, n_type)
       real(kind=c_float) :: rtdat(RDIM, AZDIM, ELDIM, n_type)
     end function jitdt_read_toshiba_c
  end interface

  public

contains

  function jitdt_read_toshiba(n_type, jitdt_place, hd, az, el, rtdat)
    integer :: jitdt_read_toshiba
    integer, intent(in) :: n_type
    character(*), intent(in) :: jitdt_place
    type(c_pawr_header), intent(out) :: hd(n_type)
    real(kind=c_float), intent(out) :: az(AZDIM, ELDIM, n_type)
    real(kind=c_float), intent(out) :: el(AZDIM, ELDIM, n_type)
    real(kind=c_float), intent(out) :: rtdat(RDIM, AZDIM, ELDIM, n_type)
    character(kind=c_char) :: c_jitdt_place*1025

    !write(*, *) "jitdt_read_toshiba_f#jitdt_read_toshiba"
    c_jitdt_place = trim(jitdt_place) // c_null_char

    jitdt_read_toshiba = jitdt_read_toshiba_c(int(n_type, c_int), c_jitdt_place, hd, az, el, rtdat)
  end function jitdt_read_toshiba
end module jitdt_read_toshiba_f

