module common_sort
!===============================================================================

  use common, only: r_size

  implicit none
  private

  public :: QUICKSELECT
  public :: QUICKSELECT_arg
  public :: QUICKSELECT_desc
  public :: QUICKSELECT_desc_arg

contains
!===============================================================================

subroutine swap(x, y)
  implicit none
  real(r_size), intent(inout) :: x, y
  real(r_size) :: tmp
  tmp = x
  x = y
  y = tmp
end subroutine swap

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

subroutine swap_i(x, y)
  implicit none
  integer, intent(inout) :: x, y
  integer :: tmp
  tmp = x
  x = y
  y = tmp
end subroutine swap_i

!-------------------------------------------------------------------------------

subroutine partition(A, left, right, pivot)
  implicit none
  real(r_size), intent(inout) :: A(:)
  integer, intent(in) :: left, right
  integer, intent(inout) :: pivot
  real(r_size) :: A_pivot
  integer :: idx, store_idx

  A_pivot = A(pivot)
  call swap(A(pivot), A(right))

  store_idx = left
  do idx = left, right-1
    if (A(idx) < A_pivot) then
      call swap(A(store_idx), A(idx))
      store_idx = store_idx + 1
    end if
  end do

  call swap(A(right), A(store_idx))
  pivot = store_idx
end subroutine partition

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

subroutine partition_arg(A, X, left, right, pivot)
  implicit none
  real(r_size), intent(in) :: A(:)
  integer, intent(inout) :: X(:)
  integer, intent(in) :: left, right
  integer, intent(inout) :: pivot
  real(r_size) :: A_pivot
  integer :: idx, store_idx

  A_pivot = A(X(pivot))
  call swap_i(X(pivot), X(right))

  store_idx = left
  do idx = left, right-1
    if (A(X(idx)) < A_pivot) then
      call swap_i(X(store_idx), X(idx))
      store_idx = store_idx + 1
    end if
  end do

  call swap_i(X(right), X(store_idx))
  pivot = store_idx
end subroutine partition_arg

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

subroutine partition_desc(A, left, right, pivot)
  implicit none
  real(r_size), intent(inout) :: A(:)
  integer, intent(in) :: left, right
  integer, intent(inout) :: pivot
  real(r_size) :: A_pivot
  integer :: idx, store_idx

  A_pivot = A(pivot)
  call swap(A(pivot), A(right))

  store_idx = left
  do idx = left, right-1
    if (A(idx) > A_pivot) then
      call swap(A(store_idx), A(idx))
      store_idx = store_idx + 1
    end if
  end do

  call swap(A(right), A(store_idx))
  pivot = store_idx
end subroutine partition_desc

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

subroutine partition_desc_arg(A, X, left, right, pivot)
  implicit none
  real(r_size), intent(in) :: A(:)
  integer, intent(inout) :: X(:)
  integer, intent(in) :: left, right
  integer, intent(inout) :: pivot
  real(r_size) :: A_pivot
  integer :: idx, store_idx

  A_pivot = A(X(pivot))
  call swap_i(X(pivot), X(right))

  store_idx = left
  do idx = left, right-1
    if (A(X(idx)) > A_pivot) then
      call swap_i(X(store_idx), X(idx))
      store_idx = store_idx + 1
    end if
  end do

  call swap_i(X(right), X(store_idx))
  pivot = store_idx
end subroutine partition_desc_arg

!-------------------------------------------------------------------------------

subroutine median_of_three(A, i1, i2, i3, i)
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
end subroutine median_of_three

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

subroutine median_of_three_arg(A, X, i1, i2, i3, i)
  implicit none
  real(r_size), intent(in) :: A(:)
  integer, intent(in) :: X(:)
  integer, intent(in) :: i1, i2, i3
  integer, intent(out) :: i

  if (A(X(i1)) < A(X(i2))) then
    if (A(X(i2)) < A(X(i3))) then
      i = i2
    else if (A(X(i1)) < A(X(i3))) then
      i = i3
    else
      i = i1
    end if
  else
    if (A(X(i1)) < A(X(i3))) then
      i = i1
    else if (A(X(i2)) < A(X(i3))) then
      i = i3
    else
      i = i2
    end if
  end if
end subroutine median_of_three_arg

!-------------------------------------------------------------------------------

subroutine sample_second_min(A, left, right, K, i_2)
  implicit none
  real(r_size), intent(in) :: A(:)
  integer, intent(in) :: left, right
  integer, intent(in) :: K
  integer, intent(out) :: i_2
  real(r_size) :: A_min, A_min_2
  integer :: i, j

  i = left
  i_2 = left
  A_min = huge(A)
  A_min_2 = huge(A)
  do j = left, right, K
    if (A(j) < A_min) then
      A_min_2 = A_min
      A_min = A(j)
      i_2 = i
      i = j
    else if (A(j) < A_min_2) then
      A_min_2 = A(j)
      i_2 = j
    end if
  end do
end subroutine sample_second_min

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

subroutine sample_second_min_arg(A, X, left, right, K, i_2)
  implicit none
  real(r_size), intent(in) :: A(:)
  integer, intent(in) :: X(:)
  integer, intent(in) :: left, right
  integer, intent(in) :: K
  integer, intent(out) :: i_2
  real(r_size) :: A_min, A_min_2
  integer :: i, j

  i = left
  i_2 = left
  A_min = huge(A)
  A_min_2 = huge(A)
  do j = left, right, K
    if (A(X(j)) < A_min) then
      A_min_2 = A_min
      A_min = A(X(j))
      i_2 = i
      i = j
    else if (A(X(j)) < A_min_2) then
      A_min_2 = A(X(j))
      i_2 = j
    end if
  end do
end subroutine sample_second_min_arg

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

subroutine sample_second_max(A, left, right, K, i_2)
  implicit none
  real(r_size), intent(in) :: A(:)
  integer, intent(in) :: left, right
  integer, intent(in) :: K
  integer, intent(out) :: i_2
  real(r_size) :: A_max, A_max_2
  integer :: i, j

  i = left
  i_2 = left
  A_max = 0.0d0 - huge(A)
  A_max_2 = 0.0d0 - huge(A)
  do j = left, right, K
    if (A(j) > A_max) then
      A_max_2 = A_max
      A_max = A(j)
      i_2 = i
      i = j
    else if (A(j) > A_max_2) then
      A_max_2 = A(j)
      i_2 = j
    end if
  end do
end subroutine sample_second_max

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

subroutine sample_second_max_arg(A, X, left, right, K, i_2)
  implicit none
  real(r_size), intent(in) :: A(:)
  integer, intent(in) :: X(:)
  integer, intent(in) :: left, right
  integer, intent(in) :: K
  integer, intent(out) :: i_2
  real(r_size) :: A_max, A_max_2
  integer :: i, j

  i = left
  i_2 = left
  A_max = 0.0d0 - huge(A)
  A_max_2 = 0.0d0 - huge(A)
  do j = left, right, K
    if (A(X(j)) > A_max) then
      A_max_2 = A_max
      A_max = A(X(j))
      i_2 = i
      i = j
    else if (A(X(j)) > A_max_2) then
      A_max_2 = A(X(j))
      i_2 = j
    end if
  end do
end subroutine sample_second_max_arg

!-------------------------------------------------------------------------------

recursive subroutine QUICKSELECT(A, left, right, K)
  implicit none
  real(r_size), intent(inout) :: A(:)
  integer, intent(in) :: left, right
  integer, intent(in) :: K
  integer :: middle, pivot

  if (left < right) then
    if ((right-left)/K >= 2) then
      call sample_second_min(A, left, right, K, pivot)
#ifdef DEBUG
      write (6, '(4I8,F10.4,A)') K, left, right, pivot, A(pivot), 'sample_second_min'
#endif
    else
      middle = (left + right) / 2
      call median_of_three(A, left, middle, right, pivot)
#ifdef DEBUG
      write (6, '(4I8,F10.4,A)') K, left, right, pivot, A(pivot), 'median_of_three'
#endif
    end if
    call partition(A, left, right, pivot)
    if (K < pivot) then
      call QUICKSELECT(A, left, pivot-1, K)
    else if (K > pivot) then
      call QUICKSELECT(A, pivot+1, right, K)
    end if
  end if
end subroutine QUICKSELECT

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

recursive subroutine QUICKSELECT_arg(A, X, left, right, K)
  implicit none
  real(r_size), intent(in) :: A(:)
  integer, intent(inout) :: X(:)
  integer, intent(in) :: left, right
  integer, intent(in) :: K
  integer :: middle, pivot

  if (left < right) then
    if ((right-left)/K >= 2) then
      call sample_second_min_arg(A, X, left, right, K, pivot)
#ifdef DEBUG
      write (6, '(4I8,F10.4,A)') K, left, right, pivot, A(X(pivot)), 'sample_second_min'
#endif
    else
      middle = (left + right) / 2
      call median_of_three_arg(A, X, left, middle, right, pivot)
#ifdef DEBUG
      write (6, '(4I8,F10.4,A)') K, left, right, pivot, A(X(pivot)), 'median_of_three'
#endif
    end if
    call partition_arg(A, X, left, right, pivot)
    if (K < pivot) then
      call QUICKSELECT_arg(A, X, left, pivot-1, K)
    else if (K > pivot) then
      call QUICKSELECT_arg(A, X, pivot+1, right, K)
    end if
  end if
end subroutine QUICKSELECT_arg

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

recursive subroutine QUICKSELECT_desc(A, left, right, K)
  implicit none
  real(r_size), intent(inout) :: A(:)
  integer, intent(in) :: left, right
  integer, intent(in) :: K
  integer :: middle, pivot

  if (left < right) then
    if ((right-left)/K >= 2) then
      call sample_second_max(A, left, right, K, pivot)
#ifdef DEBUG
      write (6, '(4I8,F10.4,A)') K, left, right, pivot, A(pivot), 'sample_second_max'
#endif
    else
      middle = (left + right) / 2
      call median_of_three(A, left, middle, right, pivot)
#ifdef DEBUG
      write (6, '(4I8,F10.4,A)') K, left, right, pivot, A(pivot), 'median_of_three'
#endif
    end if
    call partition_desc(A, left, right, pivot)
    if (K < pivot) then
      call QUICKSELECT_desc(A, left, pivot-1, K)
    else if (K > pivot) then
      call QUICKSELECT_desc(A, pivot+1, right, K)
    end if
  end if
end subroutine QUICKSELECT_desc

!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

recursive subroutine QUICKSELECT_desc_arg(A, X, left, right, K)
  implicit none
  real(r_size), intent(in) :: A(:)
  integer, intent(inout) :: X(:)
  integer, intent(in) :: left, right
  integer, intent(in) :: K
  integer :: middle, pivot

  if (left < right) then
    if ((right-left)/K >= 2) then
      call sample_second_max_arg(A, X, left, right, K, pivot)
#ifdef DEBUG
      write (6, '(4I8,F10.4,A)') K, left, right, pivot, A(X(pivot)), 'sample_second_max'
#endif
    else
      middle = (left + right) / 2
      call median_of_three_arg(A, X, left, middle, right, pivot)
#ifdef DEBUG
      write (6, '(4I8,F10.4,A)') K, left, right, pivot, A(X(pivot)), 'median_of_three'
#endif
    end if
    call partition_desc_arg(A, X, left, right, pivot)
    if (K < pivot) then
      call QUICKSELECT_desc_arg(A, X, left, pivot-1, K)
    else if (K > pivot) then
      call QUICKSELECT_desc_arg(A, X, pivot+1, right, K)
    end if
  end if
end subroutine QUICKSELECT_desc_arg

!===============================================================================
end module common_sort
