module correction
use hartree
use excorr
use int_interp
implicit none;
contains

subroutine energycorrf(density, r, rab, nr, E)
  !given the input total charge density and the radial grid, adds the exchange
  !and correlation contributions to the total energy density E
  integer nr, i;
  real(kind = kind(1.d0)) r(nr), rab(nr), E(nr), density(nr),&
  &a0, rs(nr), te, be, ex(nr), ec(nr), pi
  pi = 3.14159265359;
  a0 = (4./9./pi)**(1./3.);
  do i = 1, nr
    if (density(i).le.(1.0E-44)) then
      E(i) = 0.d0; !to avoid NaN
    else
      rs(i) = (3./density(i))**(1./3.);
      ex(i) = -3./4./(pi*rs(i)*a0); !exchange energy in hartrees
      if (rs(i).ge.1.d0) then
        be  = 1.0+1.0529*sqrt(rs(i))+0.3334*rs(i);
        ec(i) = -0.2846/2.0/be; !in hartrees
      else
        ec(i) =  ((0.0311+0.002*rs(i))*log(rs(i))-&
        0.0480-0.01160*rs(i)); !in hartrees
      end if
      E(i) = ec(i) + ex(i); !add hartree and excorr contributions to total
      !energy density
    end if
  end do
end subroutine


subroutine totalcorr(density, r, rab, nr, corrE, deltaE, Z, ec, enuc, exc)
  !given the total charge density, nuclear charge and the radial grid, calculates
  !Coulomb energe, nuclear energy, exchange-correlation energy and the correction to total
  !energy
  integer nr, i;
  real(kind = kind(1.d0)) r(nr), rab(nr), v0(nr), v1(nr), v2(nr), deltaE(nr), &
  & density(nr), corrE, ec, enuc, exc, Z;
  ec = 0.d0;
  enuc = 0.d0;
  exc= 0.d0;
  v0 = 0.d0; v1 = 0.d0; v2 = 0.d0; !these are arrays containing energy/potential
  !due to hartree and exchange-correlation, which are used to calculate
  !correction to total energy
  deltaE = 0.d0;
  call energycorrf(density, r, rab, nr, v0);
  call VHartree(density, r, rab, nr, v1);
  call vxc(density, nr, v2);
  ec = 0.5d0*norm00(v1*density*r**2, r, rab, 1, nr);
  enuc = - Z*norm00(density*r, r, rab, 1, nr);
  exc = norm00(v0*density*r**2, r, rab, 1, nr);
  deltaE = v0 - v1/2.0 - v2;
  corrE = exc - ec - norm00(v2*density*r**2, r, rab, 1, nr);
end subroutine
end module correction
