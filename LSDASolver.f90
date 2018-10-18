module LSDASolver
use hartree
use LSDAexc
use ode
use int_interp
use parameters
implicit none
contains

  subroutine spiniterates(ns, ls, spins, number, nup, ndown, vext, r, rab, nr, u, ud, a, b,  &
    & rmax, energies, Et, Z, Ecoul, Enuc, Exc, iterations)
    !spins is an array of 1s and 0s, 1 for up spin and zero for down spin
    !ns and ls are arrays containing energy level numbers and angular momenta
    !respectively. vext is external potential, -Z/r in case of an atom. Radial grid
    !is provided as well. u and ud are just variables used in calculations
    !Total(taking correction into account), nuclear, Coulomb and exchange-correlation
    !energies are calculated, as well as eigenenergies
    integer nr, i, lp, l, ls(number), ns(number), spins(number), number, nlps, &
    & iterations;
    real(kind = kind(1.d0)) nup(nr), ndown(nr), r(nr), rab(nr), u(nr), &
    &ud(nr), a, b, rmax, el1(number), vext(nr), v2(nr), corrE, ediff, ein,&
    &deltaE(nr), ang(nr), energies(:), Z, d2u(nr), d2d(nr), v3(nr), &
    & Et, Ecoul, Enuc, Exc;
    iterations = 0;
    Et = 0.d0;
    Ecoul = 0.d0;
    Enuc = 0.d0;
    Exc = 0.d0;
    corrE = 0.d0;
    ein = 0.d0;
    ediff = 100.d0;
    d2u = nup;
    d2d = ndown;
    do while(ediff > 1.0E-8)
      iterations = iterations + 1;
      nup = d2u;
      ndown = d2d;
      v2 = vext;
      v3 = vext;
      call VHartree(nup + ndown, r, rab, nr, v2); !add local Hartree contributions
      !to both spin up and down potentials
      call VHartree(nup + ndown, r, rab, nr, v3);
      call spinexc(nup, ndown, nr, v2, v3) !adding exchange potential
      call spincorr(nup, ndown, nr, v2, v3) !adding correlation potential
      !v2 for vcup, v3 for vcdown
      nup = 0.d0;
      ndown = 0.d0;
      Et = 0.d0;
      do lp = 1, number
        u=0.d0 ; ud=0.d0 ; l=ls(lp) ; energies(lp)=-0.1d0 !starting guess for eigenenergy, can be altered;
        ang(:) = 0.5d0*dble(l*(l+1))/r(:)**2
        if (spins(lp).eq.(1)) then !solving for eigenfunctions for spin up or
          !down electrons
          call difnrl(lp,v2(:)+ang(:),              &
           &         0.d0,u,ud,nr,a,b,r,rab,number,ns,ls,Z,                 &
           &         rmax,energies);
           nup(2:nr) = nup(2:nr) + (u(2:nr)/r(2:nr))**2;
           nup(1) = nup(2);
           Et = Et + energies(lp);
         else
           call difnrl(lp,v3(:)+ang(:),              &
            &         0.d0,u,ud,nr,a,b,r,rab,number,ns,ls,Z,                 &
            &         rmax,energies);
            ndown(2:nr) = ndown(2:nr) + (u(2:nr)/r(2:nr))**2;
            ndown(1) = ndown(2);
            Et = Et + energies(lp);
          end if
       end do
       d2u = 0.5d0*nup + 0.5d0*d2u; !starting density for the next iteration
       d2d = 0.5d0*ndown + 0.50d0*d2d;
       ediff = abs(Et - ein); !calculating sum of eigenenergiues difference
       !compared to previous iteration to monitor self-consistency
      ein = Et;
    end do
    call LSDAtotalcorr(nup, ndown, r, rab, nr, corrE, Ecoul, Enuc, Exc, Z);
    !calculate correction tot total energy, as well as other energies of interest
    Et = Et + corrE;
  end subroutine

end module LSDASolver
