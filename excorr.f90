module excorr
implicit none
contains
subroutine vxc(density, nr, v_xc)
  !takes the array of density and number of points in the grid
  !adds exchange-correlation potential to the potential v_xc,
  !which is given as an input
  integer nr, i;
  real(kind = kind(1.d0)) density(nr), a0, rs(nr), pi, vex(nr), &
  &vc(nr), v_xc(nr), te, be;
  pi = 3.14159265359;
  a0 = (4./9./pi)**(1./3.);
  vex = - 1./pi/a0*(density/3.)**(1./3.);
  do i = 1, nr
    if (density(i).ge.(1.0E-44)) then
      rs(i) = (3./density(i))**(1./3.);
    if (rs(i).ge.1.d0) then
      !LDA parametrization
      te = 1.0+(7.0/6.0)*1.0529*sqrt(rs(i))+(4.0/3.0)*0.3334*rs(i);
      be  = 1.0+1.0529*sqrt(rs(i))+0.3334*rs(i);
      vc(i) = -0.2846*te/2.0/be**2; !in hartrees
    else
      vc(i) = ((0.0311 + 2.0/3.0*0.002*rs(i))*log(rs(i)) - &
      (0.0480 + 0.0311/3.0) - &
      (2.0/3.0*0.0116 + 0.002/3.0)*rs(i)); !in hartrees
    end if
    v_xc(i) = vex(i) + v_xc(i) + vc(i);
  end if
  end do
end subroutine

end module excorr
