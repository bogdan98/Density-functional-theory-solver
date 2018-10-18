module LSDAexc
use int_interp
use hartree
implicit none
contains
subroutine spinexc(nup, ndown, nr, vexu, vexd)
  !exchange correlation potential for up and down spins
  !nup - up spin density, ndown - down spin density
  integer nr, i;
  real(kind = kind(1.d0)) nup(nr), ndown(nr), a0, rs, pi, vexu(nr), &
  &vexd(nr);
  pi = 3.14159265359;
  a0 = (4./9./pi)**(1./3.);
  vexu = vexu - 1./pi/a0*(2.*nup/3.)**(1./3.);
  vexd = vexd - 1./pi/a0*(2.*ndown/3.)**(1./3.);
end subroutine
subroutine spincorr(nup, ndown, nr, vcup, vcdown)
  !correlation energy and potentials for up and down spin densities
  !same principle as spinexc subroutine
  integer nr, i;
  real(kind = kind(1.d0)) nup(nr), ndown(nr), ntot(nr), rs, &
  &z(nr), f(nr), f1(nr), ecp(nr), ecu(nr), be, be1, vcup(nr), vcdown(nr), &
  vcu(nr), vcp(nr), te, te1;
  do i = 1, nr
    if ((nup(i) + ndown(i)).le.(1.0E-44)) then
      z(i) = 0.d0;
    else
      z(i) = (nup(i) - ndown(i))/(nup(i) + ndown(i));
    end if
  end do
  ntot = nup + ndown;
  f = ((1 + z)**(4./3.) + (1 - z)**(4./3.) - 2)/(2**(4./3.) - 2);
  f1 = 4./3.*((1 + z)**(1./3.) - (1 - z)**(1./3.))/(2**(4./3.) - 2);
  do i = 1, nr
    if (ntot(i).gt.(1.0E-44)) then
      rs = (3./ntot(i))**(1./3.); !to avoid NaN
      if (rs.ge.1.d0) then
        !parametrization of excorr in LSDA
        be  = 1.0+1.0529*sqrt(rs)+0.3334*rs;
        be1 = 1.0+1.3981*sqrt(rs)+0.2611*rs;
        te = 1.0+(7.0/6.0)*1.0529*sqrt(rs)+(4.0/3.0)*0.3334*rs;
        vcu(i) = -0.2846*te/2.0/be**2; !in hartrees
        te1 = 1.0+(7.0/6.0)*1.3981*sqrt(rs)+(4.0/3.0)*0.2611*rs;
        vcp(i) = -0.1686*te1/2.0/be1**2; !in hartrees
        ecu(i) = -0.2846/2.0/be; !in hartrees
        ecp(i) = -0.1686/2.0/be1;  !in hartrees
      else
        vcu(i) = ((0.0311 + 2.0/3.0*0.002*rs)*log(rs) - &
        (0.0480 + 0.0311/3.0) - &
        (2.0/3.0*0.0116 + 0.002/3.0)*rs); !in hartrees
        vcp(i) = ((0.01555 + 2.0/3.0*0.0007*rs)*log(rs) - &
        (0.0269 + 0.01555/3.0) - &
        (2.0/3.0*0.0048 + 0.0007/3.0)*rs); !in hartrees
        ecu(i) =  ((0.0311+0.002*rs)*log(rs)-&
        0.0480-0.01160*rs); !in hartrees
        ecp(i) =  ((0.01555+0.0007*rs)*log(rs)-&
        0.0269-0.0048*rs); !in hartrees
      end if
    end if
  end do
  vcup = vcup + (1 - f)*vcu + f*vcp + (ecp - ecu)*(+1 - z)*f1;
  vcdown = vcdown + (1 - f)*vcu + f*vcp + (ecp - ecu)*(-1 - z)*f1;
end subroutine
subroutine LSDAtotalcorr(nup, ndown, r, rab, nr, corrE, Ec, En, Exc, AtNum)
  !given up and down spin densities, atomic number and the radial grid,
  !calculates correction to total energy, Coulomb, nuclear and exchange-correlation energies
  integer nr, i;
  real(kind = kind(1.d0)) r(nr), rab(nr), vexu(nr), vexd(nr), eex(nr), ecu(nr), ecp(nr),&
  & nup(nr), ndown(nr), ntot(nr), corrE, vha(nr), ecorr(nr), vcup(nr), &
  & vcdown(nr), z(nr), f(nr), rs, be, be1, Ek, Ec, En, Exc, AtNum;
  vha = 0.d0;
  vcup = 0.d0;
  vcdown = 0.d0;
  vexu = 0.d0;
  vexd = 0.d0;
  do i = 1, nr
    if (nup(i).le.(1.0E-44)) then
      nup(i) = 0.d0; !to avoid NaN
    end if
    if (ndown(i).le.(1.0E-44)) then
      ndown(i) = 0.d0; !to avoid NaN
    end if
  end do
  do i = 1, nr
    if ((nup(i) + ndown(i)).eq.(0.d0)) then
      z(i) = 0.d0; !to prevent NaN coming put
    else
      z(i) = (nup(i) - ndown(i))/(nup(i) + ndown(i));
    end if
  end do
  ntot = nup + ndown;
  f = ((1 + z)**(4./3.) + (1 - z)**(4./3.) - 2)/(2**(4./3.) - 2);
  corrE = 0.d0;
  do i = 1, nr
    if ((ntot(i)).ge.(1.0E-44)) then
      rs = (3./ntot(i))**(1./3.);
      if (rs.ge.1.d0) then
        !LSDA parametrization
        be  = 1.0+1.0529*sqrt(rs)+0.3334*rs;
        ecu(i) = -0.2846/2.0/be; !in hartrees
        be1 = 1.0+1.3981*sqrt(rs)+0.2611*rs;
        ecp(i) = -0.1686/2.0/be1;  !in hartrees
      else
        ecu(i) =  ((0.0311+0.002*rs)*log(rs)-&
        0.0480-0.01160*rs); !in hartrees
        ecp(i) =  ((0.01555+0.0007*rs)*log(rs)-&
        0.0269-0.0048*rs); !in hartrees
      end if
    end if
  end do
  ecorr = (1 - f)*ecu + f*ecp; !correlation energy correction
  call spinexc(nup, ndown, nr, vexu, vexd);!exchange energy part
  call VHartree(ntot, r, rab, nr, vha); !hartree energy part
  call spincorr(nup, ndown, nr, vcup, vcdown); !correlation energy part
  Ec = 0.5d0*norm00(vha*ntot*r**2, r, rab, 1, nr); !Coulomb energy
  En = -AtNum*norm00(ntot*r, r, rab, 1, nr); !nuclear energy
  Exc = norm00((ecorr )*ntot*r**2, r, rab, 1, nr) + norm00((3./4.*vexu)*nup*r**2, &
  & r, rab, 1, nr) + norm00((3./4.*vexd)*ndown*r**2, r, rab, 1, nr); !exchange-correlation energy
  corrE = Exc - Ec - &
  &norm00((vcup + vexu)*nup*r**2, r, rab, 1, nr)&
  & - norm00((vcdown + vexd)*ndown*r**2, r, rab, 1, nr);
end subroutine
end module LSDAexc
