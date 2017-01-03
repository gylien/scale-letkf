module common_sort
  use common, only: r_size

  implicit none
  public

contains

subroutine swap(x, y)
  implicit none
  real(r_size), intent(inout) :: x, y
  real(r_size) :: tmp
  tmp = x
  x = y
  y = tmp
end subroutine swap

subroutine swap_i(x, y)
  implicit none
  integer, intent(inout) :: x, y
  integer :: tmp
  tmp = x
  x = y
  y = tmp
end subroutine swap_i

subroutine partition(A, left, right, pivot, B, C, I)
  implicit none
  real(r_size), intent(inout) :: A(:)
  integer, intent(in) :: left, right
  integer, intent(inout) :: pivot
  real(r_size), intent(inout), optional :: B(:), C(:)
  integer, intent(inout), optional :: I(:)
  real(r_size) :: A_pivot
  integer :: idx, store_idx

  A_pivot = A(pivot)
  call swap(A(pivot), A(right))
  if (present(B)) call swap(B(pivot), B(right))
  if (present(C)) call swap(C(pivot), C(right))
  if (present(I)) call swap_i(I(pivot), I(right))

  store_idx = left
  do idx = left, right-1
    if (A(idx) < A_pivot) then
      call swap(A(store_idx), A(idx))
      if (present(B)) call swap(B(store_idx), B(idx))
      if (present(C)) call swap(C(store_idx), C(idx))
      if (present(I)) call swap_i(I(store_idx), I(idx))
      store_idx = store_idx + 1
    end if
  end do

  call swap(A(right), A(store_idx))
  if (present(B)) call swap(B(right), B(store_idx))
  if (present(C)) call swap(C(right), C(store_idx))
  if (present(I)) call swap_i(I(right), I(store_idx))
  pivot = store_idx
end subroutine partition

subroutine middle_of_three(A, i1, i2, i3, i)
  implicit none
  real(r_size), intent(in) :: A(:)
  integer, intent(in) :: i1, i2, i3
  integer, intent(out) :: i

  if (A(i1) < A(i2)) then
    if (A(i2) < A(i3)) then
      i = i2
    else if (A(i1) < A(i3)) then
      i = i3
    else
      i = i1
    end if
  else
    if (A(i1) < A(i3)) then
      i = i1
    else if (A(i2) < A(i3)) then
      i = i3
    else
      i = i2
    end if
  end if
end subroutine middle_of_three

recursive subroutine QUICKSELECT(A, left, right, K, B, C, I)
  implicit none
  real(r_size), intent(inout) :: A(:)
  integer, intent(in) :: left, right
  integer, intent(in) :: K
  real(r_size), intent(inout), optional :: B(:), C(:)
  integer, intent(inout), optional :: I(:)
  integer :: middle, pivot

  if (left < right) then
    middle = (left + right) / 2
    call middle_of_three(A, left, middle, right, pivot)
    call partition(A, left, right, pivot, B, C, I)
    if (K < pivot) then
      call QUICKSELECT(A, left, pivot-1, K, B, C, I)
    else if (K > pivot) then
      call QUICKSELECT(A, pivot+1, right, K, B, C, I)
    end if
  end if
end subroutine QUICKSELECT

end module common_sort
