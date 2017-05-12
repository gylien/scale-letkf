module common_rand
!===============================================================================

  implicit none
  private

  public :: Knuth_Shuffle

  logical :: random_seed_init = .false.

contains
!===============================================================================

subroutine Knuth_Shuffle(num, a)
  implicit none
  integer, intent(in) :: num
  integer, intent(out) :: a(num)
  integer :: i, randpos, tmp
  real(8) :: r
 
  do i = 1, num
    a(i) = i
  end do

  if (.not. random_seed_init) then
    call random_seed()
    random_seed_init = .true.
  end if

  do i = num, 2, -1
    call random_number(r)
    randpos = int(r * i) + 1
    tmp = a(randpos)
    a(randpos) = a(i)
    a(i) = tmp
  end do
end subroutine Knuth_Shuffle

!===============================================================================
end module common_rand
