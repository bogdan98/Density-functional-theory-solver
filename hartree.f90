module hartree
use int_interp, only : norm00
implicit none

contains
subroutine VHartree(n, r, rab, nr, V)
  !calculates Hartree potential at all points of r, given input
  !density (n) and radial grid. Adds the Hartree potential to V, which
  !is given as an input
  integer :: nr, i;
  real(kind  = kind(1.d0)) :: n(nr), r(nr), V(nr), vec(nr), rab(nr);
  do i = 2, nr-1
    vec(1:i) = (r(1:i)**2)*n(1:i)/r(i);
    vec(i+1:nr) = r(i+1:nr)*n(i+1:nr);
    V(i) = V(i) + norm00(vec, r, rab, 1, nr);
  end do
  V(1) = V(1) + norm00(n(:)*r(:), r, rab, 1, nr);
  V(nr) = V(nr) + norm00((r(:)**2)*n(:)/r(nr), r, rab, 1, nr);
  !4*3.14159265359*
end subroutine
subroutine GSimpson(x, y, n, sum)
  !juas a subroutine for integration by 1/3 Simpson's rule
  integer n, i;
  REAL (kind = kind(1.d0)) :: x(n), y(n), sum
  do i = 3, n, 2
    sum = sum + y(i-2)*(x(i) -x(i-2)) + (x(i) - x(i-2))**2/2*sqb2(i-2, i-1, x, y, n) + &
    & 0.5d0*sqb3(i-2, i-1, i, x, y, n)*((x(i) - x(i-2))*(x(i) - x(i-1))**2 - (x(i) - &
    &x(i-1))**3/3 + (x(i-2) - x(i-1))**3/3);
  end do
end subroutine
REAL function sqb2(i1, i2, x, y, n)
  REAL (kind = kind(1.d0)) :: x(n), y(n)
  INTEGER ::i1, i2, n
  sqb2 = (y(i2) - y(i1))/(x(i2) - x(i1))
  return
end function
REAL function sqb3(i1, i2, i3, x, y, n)
  REAL (kind = kind(1.d0)) :: x(n), y(n)
  INTEGER ::i1, i2, i3, n
  sqb3 = (sqb2(i2, i3, x, y, n) - sqb2(i1, i2, x, y, n))/(x(i3) - x(i1));
  return
end function

end module hartree
