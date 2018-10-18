MODULE ode
USE int_interp, only : dinterp
IMPLICIT NONE

CONTAINS


!--- Solve SE with supplied potential and grid.
 SUBROUTINE difnrl(iorb,v,v0,ar,br,nr,a,b,r,rab,norb,no,lo, &
&                  znuc,rwell,ev)
 IMPLICIT NONE
!- integrate the Schroedinger equation:
!- find eigenvalue ev, wavefunction ar
!- and its derivative br = d(ar)/dr
!- $Id: difnrl.x,v 1.1 89/10/26 19:53:18 sverre Exp $

!- $Log:	difnrl.x,v $
!- Revision 1.1  89/10/26  19:53:18  sverre
!- Initial revision
 
 integer iorb,nr,norb,no(:),lo(:),itmax
 real(kind=kind(1.d0)) v(:),ar(:),br(:),r(:),rab(:),ev(:)
 real(kind=kind(1.d0)) v0,a,b,znuc,rwell

 integer i,j,k,jj,jm1,ll,lp,nwell,ninf,nctp,idone,nodes,istart,istop
 real(kind=kind(1.d0)) r2(nr)
 real(kind=kind(1.d0)) dlnto2,zeff,aa,bb,vzero,var0,emax,emin,evold,devmax, &
  &                    fa0,fa1,fa2,fa3,fa4,fa5,fb0,fb1,fb2,fb3,fb4,fb5,     &
  &                    arp,brp,arc,brc,arctp,brctp,alf,dnwell,rj,rabj,        &
  &                    factor,dev
 real(kind=kind(1.d0)) abc1,abc2,abc3,abc4,abc5,amc0,amc1,amc2,amc3,amc4,tol
 parameter(				&
& abc1 = 1901.D0/720.D0,		&
& abc2 = -1387.D0/360.D0,		&
& abc3 = 109.D0/30.D0,			&
& abc4 = -637.D0/360.D0,		&
& abc5 = 251.D0/720.D0,			&
& amc0 = 251.D0/720.D0,			&
& amc1 = 323.D0/360.D0,			&
& amc2 = -11.D0/30.D0,			&
& amc3 = 53.D0/360.D0,			&
& amc4 = -19.D0/720.D0,			&
& tol = 1.D-13,				&
& itmax = 500)				!     eigenvalue tolerance and max number of iterations

!-    determine effective charge and vzero for
!-    startup of outward integration
!-    
!-    ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
!-    
!-    aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)

 ev=2.0d0*ev ; v=2.0d0*v ; v0=2.d0*v0			! hartree to rydbergs

 do i=1,nr
  r2(i)=r(i)*r(i)
 enddo
!---  Want a larger radius than current ninf definition provides, so modifiy this v.
!dlnto2=log(tol)**2
 dlnto2=log(tol/1000d0)**2

 lp = lo(iorb)+1
 zeff = znuc
 aa = -zeff/lp
 vzero = -2.d0*zeff*aa + v0				! used to make bb later
 var0 = 0.D0						! u''(r=0) for l>1
 if (lo(iorb) .eq. 0) var0 = -2.d0*zeff			!              l=0
 if (lo(iorb) .eq. 1) var0 = 2.D0			!              l=1

!-    these are used to bracket eigenvalue

 emax = +1.D+20
 emin = -1.D+20

!-    max step size for eigenvalue changes

 devmax = -ev(iorb) * 0.2d0
 if (devmax .lt. 0.3D0) devmax = 0.3D0

!-    begin iteration loop

 do 190 i=1,itmax

!-    find closest point inside rwell - nwell,
!-    practical infinity ninf, and
!-    classical turning point nctp

   nwell = nr
   ninf = nr
   nctp = nr
   do jj=2,nr
     j = nr-jj+2
     ar(j) = 0.D0
     br(j) = 0.D0
     idone = 1
     if (r(j) .gt. rwell) then
       nwell = j - 1
       idone = 0
     end if
     if (r2(j)*(v(j)-ev(iorb)) .gt. dlnto2) then
       ninf = j
       idone = 0
     end if
     if (v(j) .gt. ev(iorb)) then
       nctp = j
       idone = 0
     end if
     if (idone .eq. 1) goto 110
   enddo

!-    three possibilities (nwell is normally equal to nr)

!     nctp < ninf < nwell  -- normal case, exponetial inward startup
!     nctp < nwell < ninf  -- bounded case, power series inward startup
!     nwell < nctp         -- bounded case, no inward integration

!     reset ninf and nctp to allow at least two inward startup points

 110 if (ninf .gt. nwell) ninf = nwell
   if (nctp .gt. nwell - 1) nctp = nwell - 1


!     outward integration from 1 to nctp -- startup

   bb = (vzero-ev(iorb))/(4.d0*lp+2.d0)
   ar(1) = 0.D0
   br(1) = 0.D0
   if (lo(iorb) .eq. 0) br(1) = b*a
   do j=2,5
     ar(j) = r(j)**lp * (1.d0+(aa+bb*r(j))*r(j))
     br(j) = rab(j)*r(j)**lo(iorb)                       &
    &        * (lp+(aa*(lp+1)+bb*(lp+2.d0)*r(j))*r(j))
   enddo
   fa5 = br(1)
   fb5 = b*br(1) + rab(1)*rab(1)*var0
   fa4 = br(2)
   fb4 = b*br(2) + rab(2)*rab(2)*(v(2)-ev(iorb))*ar(2)
   fa3 = br(3)
   fb3 = b*br(3) + rab(3)*rab(3)*(v(3)-ev(iorb))*ar(3)
   fa2 = br(4)
   fb2 = b*br(4) + rab(4)*rab(4)*(v(4)-ev(iorb))*ar(4)
   fa1 = br(5)
   fb1 = b*br(5) + rab(5)*rab(5)*(v(5)-ev(iorb))*ar(5)

!     outward integration loop

   nodes = 0
   do j=6,nctp

!     predictor (Adams-Bashforth)

     arp = ar(j-1) + abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5
     brp = br(j-1) + abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5
     fa0 = brp
     fb0 = b*brp + rab(j)*rab(j)*(v(j)-ev(iorb))*arp

!     corrector (Adams-Moulton)

     jm1=j-1
     arc = ar(jm1) + amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4
     brc = br(jm1) + amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4
     fa5 = fa4
     fb5 = fb4
     fa4 = fa3
     fb4 = fb3
     fa3 = fa2
     fb3 = fb2
     fa2 = fa1
     fb2 = fb1
     fa1 = brc
     fb1 = b*brc + rab(j)*rab(j)*(v(j)-ev(iorb))*arc
     ar(j) = arc + amc0*(fa1-fa0)
     br(j) = brc + amc0*(fb1-fb0)
     fa1 = br(j)
     fb1 = b*br(j) + rab(j)*rab(j)*(v(j)-ev(iorb))*ar(j)

!     count nodes

     if (ar(j)*ar(jm1) .le. 0) nodes = nodes + 1
   enddo

!     end outward integration

!     if incorrect number of nodes modify energy stepwise

   if (nodes .gt. no(iorb)-lo(iorb)-1) then 		! too many nodes -- decrease ev
     if (ev(iorb) .lt. emax) emax = ev(iorb)
     if (devmax .gt. 0.D0) devmax = -devmax * 0.5d0
     ev(iorb) = ev(iorb) + devmax
     goto 190
   elseif (nodes .lt. no(iorb)-lo(iorb)-1) then		! too few nodes -- increase ev
     if (ev(iorb) .gt. emin) emin = ev(iorb)
     if (devmax .lt. 0.D0) devmax = -devmax * 0.5d0
     ev(iorb) = ev(iorb) + devmax
     goto 190
   endif

!     correct number of nodes

   arctp = ar(nctp)
   brctp = br(nctp)

!     inward integration from ninf to nctp -- startup

   if (ninf .lt. nwell) then

!     normal startup

     istart = ninf - nctp + 1
     if (istart .gt. 5) istart = 5
     do jj=1,istart
       j = ninf-jj+1
       alf = v(j) - ev(iorb)
       if (alf .lt. 0.D0) alf = 0.D0
       alf = sqrt(alf)
       ar(j) = exp(-alf*r(j))
       br(j) = -rab(j)*alf*ar(j)
     enddo
     fa5 = br(ninf)
     fb5 = b*br(ninf) + rab(ninf)*rab(ninf)*(v(ninf)-ev(iorb))*ar(ninf)	
     fa4 = br(ninf-1)
     fb4 = b*br(ninf-1) + rab(ninf-1)*rab(ninf-1)*(v(ninf-1)-ev(iorb))*ar(ninf-1)
     fa3 = br(ninf-2)
     fb3 = b*br(ninf-2) + rab(ninf-2)*rab(ninf-2)*(v(ninf-2)-ev(iorb))*ar(ninf-2)
     fa2 = br(ninf-3)
     fb2 = b*br(ninf-3) + rab(ninf-3)*rab(ninf-3)*(v(ninf-3)-ev(iorb))*ar(ninf-3)
     fa1 = br(ninf-4)
     fb1 = b*br(ninf-4) + rab(ninf-4)*rab(ninf-4)*(v(ninf-4)-ev(iorb))*ar(ninf-4)
     istart = 5
   else

!     power series startup

!     find (v(rwell) - ev(iorb) / 6

     if (nwell .lt. nr) then
       alf = (v(nwell)*(r(nwell+1)-rwell)           &
        &     - v(nwell+1)*(r(nwell)-rwell))        &
        &     / (r(nwell+1) - r(nwell))
     else
       alf = v(nr)
     endif
     alf = (alf - ev(iorb)) * .16666666666666666666666666666666d0

     dnwell=dble(nwell)
     rj = a * (exp(b*(dnwell-1.d0)) - 1.d0)
     rabj = (rj + a) * b
     rj = rj - rwell
     fa1 = -(1.d0 + 3.d0 * alf * rj**2) * rabj
     fb1 = fa1 * b - 6.d0 * alf * rj * rabj * rabj
     ar(nwell) = -rj * (1.d0 + alf * rj**2)
     br(nwell) = fa1

     rj = a * (exp(b*(dnwell)) - 1.d0)
     rabj = (rj + a) * b
     rj = rj - rwell
     fa2 = -(1.d0 + 3.d0 * alf * rj**2) * rabj
     fb2 = fa2 * b - 6.d0 * alf * rj * rabj * rabj

     rj = a * (exp(b*(dnwell+1.d0)) - 1.d0)
     rabj = (rj + a) * b
     rj = rj - rwell
     fa3 = -(1.d0 + 3.d0 * alf * rj**2) * rabj
     fb3 = fa3 * b - 6.d0 * alf * rj * rabj * rabj

     rj = a * (exp(b*(dnwell+2.d0)) - 1.d0)
     rabj = (rj + a) * b
     rj = rj - rwell
     fa4 = -(1.d0 + 3.d0 * alf * rj**2) * rabj
     fb4 = fa4 * b - 6.d0 * alf * rj * rabj * rabj

     rj = a * (exp(b*(dnwell+3.d0)) - 1.d0)
     rabj = (rj + a) * b
     rj = rj - rwell
     fa5 = -(1.d0 + 3.d0 * alf * rj**2) * rabj
     fb5 = fa5 * b - 6.d0 * alf * rj * rabj * rabj

     istart = 1
   endif

!     integration loop

   istop = ninf - nctp
   do jj=istart,istop
     j = ninf - jj

!     predictor (Adams-Bashforth)

     arp = ar(j+1) - (abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5)
     brp = br(j+1) - (abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5)
     fa0 = brp
     fb0 = b*brp + rab(j)*rab(j)*(v(j)-ev(iorb))*arp

!     corrector (Adams-Moulton)

     arc = ar(j+1) - (amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4)
     brc = br(j+1) - (amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4)
     fa5 = fa4
     fb5 = fb4
     fa4 = fa3
     fb4 = fb3
     fa3 = fa2
     fb3 = fb2
     fa2 = fa1
     fb2 = fb1
     fa1 = brc
     fb1 = b*brc + rab(j)*rab(j)*(v(j)-ev(iorb))*arc
     ar(j) = arc - amc0*(fa1-fa0)
     br(j) = brc - amc0*(fb1-fb0)
     fa1 = br(j)
     fb1 = b*br(j) + rab(j)*rab(j)*(v(j)-ev(iorb))*ar(j)
   enddo

!     end inward integration

!     rescale ar and br outside nctp to match
!     ar(nctp) from outward integration

   factor = arctp/ar(nctp)
   do j=nctp,ninf
     ar(j) = factor * ar(j)
     br(j) = factor * br(j)
   enddo

!     find normalization

   factor = 0.D0
   ll = 4
   do j=2,ninf
     factor = factor + ll*ar(j)*ar(j)*rab(j)
     ll = 6 - ll
   enddo
   factor = factor *0.3333333333333333333333333333333333333333d0

!     modify eigenvalue ev

   dev = arctp * (brctp-br(nctp)) / (factor * rab(nctp))

!     resort to bisection if dev too large

   if (abs(dev) .gt. abs(devmax)) then
     if (devmax*dev .lt. 0.D0) devmax = -devmax * 0.5d0
     dev = devmax
   endif
   evold = ev(iorb)
   ev(iorb) = ev(iorb) + dev
   if (ev(iorb) .gt. emax) ev(iorb) = (evold + emax) * 0.5d0
   if (ev(iorb) .lt. emin) ev(iorb) = (evold + emin) * 0.5d0
   if (abs(dev) .lt. tol*(1+abs(ev(iorb)))) goto 220

 190  continue

!     eigenpar not converged in itmax iterations

!     if missing -- find normalization

   if (nodes .ne. no(iorb)-lo(iorb)-1) then
     factor = 0.D0
     ll = 4
     do j=2,ninf
       factor = factor + ll*ar(j)*ar(j)*rab(j)
       ll = 6 - ll
     enddo
   endif
   factor = factor * 0.33333333333333333333333333333333333333d0

!     error message

 210 format(' orb #',i3,' did not converge',/,            &
    &     ' ev =',e18.10,' nodes =',i2,' dev =',e18.10)

!     normalize wavefunction and change br from d(ar)/dj to d(ar)/dr

 220 factor = 1 / sqrt(factor)
   do j=1,ninf
     ar(j) = factor*ar(j)
     br(j) = factor*br(j) / rab(j)
   enddo

   ev=0.5d0*ev ; v=0.5d0*v ; v0=0.5d0*v0		! rydbergs to hartrees


 return

 END SUBROUTINE difnrl


!--- Solve GENERALISED SE with supplied potential and grid.
!--- Generalised SE solved for is (IN HARTREES):
![0) -nabla.A(r).nabla.psi + B(r)/r^2.L^2.psi - 1/2.C(r).nabla.psi - 1/2.nabla.C(r).psi + v(r).psi = e.psi]
!--- Generalised SE solved for is (IN RYDBERGS):
! 1) -nabla.A(r).nabla.psi + B(r)/r^2.L^2.psi -     C(r).nabla.psi -     nabla.C(r).psi + v(r).psi = e.psi
!--- Spherical coords (and removing Ylm from psi) reduces this to:
! 2) -1/r^2.d/dr[r^2.A(r).dpsi/dr] + B(r).l(l+1)/r^2.psi - C(r).d/dr[psi] - d/dr[C(r).psi] + v(r).psi = e.psi
!--- Differentiating to move A,B,C fields to LHS reduces this to:
! 3) -A.1/r^2.d/dr[r^2.dpsi/dr] - A'.dpsi/dr + B(r).l(l+1)/r^2.psi - 2.C.d/dr[psi] - C'.psi + v(r).psi = e.psi
!--- In terms of u, with psi=u/r:
! 4) - A.u" - [ A' + 2.C ].u' + [ A' + 2.C ].1/r.u + B.l(l+1)/r^2.u - C'.u + (v-e).u = 0
!--- Finally to form psi = r^{n-1}.y or u=r^n.y with n fractional and y polynomial in limit r->0. [ For A0=B0 n=l+1 ].
! 5) - A.y" - [ 2.A.n.1/r + A'+2C ].y' - [ A.n.(n-1) - B.l(l+1) ].1/r^2.y - [ A'+2C ].(n-1).1/r.y - C'.y + (v-e).y = 0
!--- 	  [ 2.A.n.1/r + A'+2C ]						term1
!--- 	- [ A.n.(n-1) - B.l(l+1) ].1/r^2 - [ A'+2C ].(n-1).1/r		term0
!---	- C'.y + (v-e).y						in place
 SUBROUTINE difnrl_gen(iorb,am,amd,bm,cm,fm,v,ar,br,cr,nr,a,b,r,rab,norb,no,lo, &
&                  znuc,rwell,ev,y,yd,ydd)
 IMPLICIT NONE
!- integrate the generalised Schroedinger equation:
!- find eigenvalue ev, wavefunction ar
!- and its derivative br = d(ar)/dr
!- Based on difnrl above.

 integer iorb,nr,norb,no(:),lo(:),itmax
 real(kind=kind(1.d0)) v(:),ar(:),br(:),cr(:),r(:),rab(:),ev(:)
 real(kind=kind(1.d0)) y(:),yd(:),ydd(:)
 real(kind=kind(1.d0)) v0,a,b,znuc,rwell
 real(kind=kind(1.d0)) am(:),amd(:),bm(:),cm(:),fm(:)		! Provided on full grid.
 real(kind=kind(1.d0)) amdd, bmd, bmdd, cmd			! Evaluated here at r=0 only.

 integer i,j,k,jj,jm1,ll,lp,nwell,ninf,nctp,idone,nodes,istart,istop
 real(kind=kind(1.d0)) r2(nr)
 real(kind=kind(1.d0)) dlnto2,zeff,aa,bb,vzero,var0,emax,emin,evold,devmax, &
  &                    fa0,fa1,fa2,fa3,fa4,fa5,fb0,fb1,fb2,fb3,fb4,fb5,     &
  &                    arp,brp,arc,brc,arctp,brctp,alf,dnwell,rj,rabj,        &
  &                    factor,dev
 real(kind=kind(1.d0)) abc1,abc2,abc3,abc4,abc5,amc0,amc1,amc2,amc3,amc4,tol
 real(kind=kind(1.d0)) nn,t1,t2,t3							! JRT
 parameter(				&
& abc1 = 1901.D0/720.D0,		&
& abc2 = -1387.D0/360.D0,		&
& abc3 = 109.D0/30.D0,			&
& abc4 = -637.D0/360.D0,		&
& abc5 = 251.D0/720.D0,			&
& amc0 = 251.D0/720.D0,			&
& amc1 = 323.D0/360.D0,			&
& amc2 = -11.D0/30.D0,			&
& amc3 = 53.D0/360.D0,			&
& amc4 = -19.D0/720.D0,			&
& tol = 1.D-13,				&
& itmax = 500)				!     eigenvalue tolerance and max number of iterations

!-    determine effective charge and vzero for
!-    startup of outward integration
!-    
!-    ar = (1 + aa r + bb r**2 + ... )
!-    

 ev=2.0d0*ev ; v=2.0d0*v 				! hartree to rydbergs
 v0=v(1)

 do i=1,nr
  r2(i)=r(i)*r(i)
 enddo
!---  Want a larger radius than current ninf definition provides, so modifiy this v.
!dlnto2=log(tol)**2
 dlnto2=log(tol/1000d0)**2

 lp = lo(iorb)+1
 zeff = znuc
 nn = 0.5d0 + dsqrt( 0.25d0 + (lp-1)*lp*bm(1)/am(1) )					! JRT
 call dinterp(0.d0,1,r,am,nr,t1,t2,amdd)						! JRT
 call dinterp(0.d0,1,r,bm,nr,t1,bmd,bmdd)						! JRT
 call dinterp(0.d0,1,r,am,nr,t1,cmd,t3)							! JRT
 aa = amd(1)*(nn+1)*(nn-1) - bmd*(lp-1)*lp + cm(1)*(nn-1) + 2.d0*zeff			! JRT
 aa = -aa/( 2.0d0*am(1)*nn )								! JRT
 vzero = ( amd(1)*nn*(nn+2) - bmd*(lp-1)*lp + cm(1)*nn + 2.d0*zeff )			! JRT
 vzero = vzero*aa + ( amdd*(nn-1)*(nn+2) - bmdd*(lp-1)*lp + cmd*(nn-1) - v0 )		! JRT
 vzero = -vzero										! JRT used to make bb later

!-    these are used to bracket eigenvalue

 emax = +1.D+20
 emin = -1.D+20

!-    max step size for eigenvalue changes

 devmax = -ev(iorb) * 0.2d0
 if (devmax .lt. 0.3D0) devmax = 0.3D0

!-    begin iteration loop

 do 190 i=1,itmax

!-    find closest point inside rwell - nwell,
!-    practical infinity ninf, and
!-    classical turning point nctp

   nwell = nr
   ninf = nr
   nctp = nr
   do jj=2,nr
     j = nr-jj+2
     ar(j) = 0.D0
     br(j) = 0.D0
     idone = 1
     if (r(j) .gt. rwell) then
       nwell = j - 1
       idone = 0
     end if
     if (r2(j)*(bm(j)*lp*(lp-1)/r2(j)+v(j)-ev(iorb)) .gt. dlnto2) then			! JRT ***
       ninf = j
       idone = 0
     end if
     if (bm(j)*lp*(lp-1)/r2(j)+v(j) .gt. ev(iorb)) then					! JRT ***
       nctp = j
       idone = 0
     end if
     if (idone .eq. 1) goto 110
   enddo

!-    three possibilities (nwell is normally equal to nr)

!     nctp < ninf < nwell  -- normal case, exponetial inward startup
!     nctp < nwell < ninf  -- bounded case, power series inward startup
!     nwell < nctp         -- bounded case, no inward integration

!     reset ninf and nctp to allow at least two inward startup points

 110 if (ninf .gt. nwell) ninf = nwell
   if (nctp .gt. nwell - 1) nctp = nwell - 1


!     outward integration from 1 to nctp -- startup

   bb = (vzero-ev(iorb))/( am(1)*(4.d0*nn + 2.d0) )					! JRT
   var0=2.0d0*bb									! JRT
   ar(1) = 1.D0										! JRT
   br(1) = rab(1)*aa									! JRT
   cr(1) = rab(1)*rab(1)*2.d0*bb							! JRT
   do j=2,5
     ar(j) = (1.d0+(aa+bb*r(j))*r(j))
     br(j) = rab(j)*(aa+2.d0*bb*r(j))
     cr(j) = rab(j)*rab(j)*2.d0*bb
   enddo
   fa5 = br(1)
   fa4 = br(2)
   fa3 = br(3)
   fa2 = br(4)
   fa1 = br(5)
  fb5 = b*br(1) + rab(1)*rab(1)*var0							! JRT
  fb4 = b*br(2) + term1(2,am,amd,bm,cm,r,nn,lp,rab)*br(2) + &				! JRT
   &              term0(2,am,amd,bm,cm,r,nn,lp,rab)*ar(2) + &				! JRT
   &              rab(2)*rab(2)*(v(2)-ev(iorb))/am(2)*ar(2)				! JRT
  fb3 = b*br(3) + term1(3,am,amd,bm,cm,r,nn,lp,rab)*br(3) + &				! JRT
   &              term0(3,am,amd,bm,cm,r,nn,lp,rab)*ar(3) + &				! JRT
   &              rab(3)*rab(3)*(v(3)-ev(iorb))/am(3)*ar(3)				! JRT
  fb2 = b*br(4) + term1(4,am,amd,bm,cm,r,nn,lp,rab)*br(4) + &				! JRT
   &              term0(4,am,amd,bm,cm,r,nn,lp,rab)*ar(4) + &				! JRT
   &              rab(4)*rab(4)*(v(4)-ev(iorb))/am(4)*ar(4)				! JRT
  fb1 = b*br(5) + term1(5,am,amd,bm,cm,r,nn,lp,rab)*br(5) + &				! JRT
   &              term0(5,am,amd,bm,cm,r,nn,lp,rab)*ar(5) + &				! JRT
   &              rab(5)*rab(5)*(v(5)-ev(iorb))/am(5)*ar(5)				! JRT
! A.y" + [ 2.A.n.1/r + A'+2C ].y' + [ A.n(n-1)-B.l(l+1) ].1/r^2.y + [A'+2C].(n-1).1/r.y - [C'+v-e].y = 0

!     outward integration loop

   nodes = 0
   do j=6,nctp

!     predictor (Adams-Bashforth)

     arp = ar(j-1) + abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5
     brp = br(j-1) + abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5
     fa0 = brp
     fb0 = b*brp + term1(j,am,amd,bm,cm,r,nn,lp,rab)*brp + &				! JRT
   &               term0(j,am,amd,bm,cm,r,nn,lp,rab)*arp + &				! JRT
   &               rab(j)*rab(j)*(v(j)-ev(iorb))/am(j)*arp     				! JRT

!     corrector (Adams-Moulton)

     jm1=j-1
     arc = ar(jm1) + amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4
     brc = br(jm1) + amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4
     fa5 = fa4
     fb5 = fb4
     fa4 = fa3
     fb4 = fb3
     fa3 = fa2
     fb3 = fb2
     fa2 = fa1
     fb2 = fb1
     fa1 = brc
     fb1 = b*brc + term1(j,am,amd,bm,cm,r,nn,lp,rab)*brc + &				! JRT
   &               term0(j,am,amd,bm,cm,r,nn,lp,rab)*arc + &				! JRT
   &               rab(j)*rab(j)*(v(j)-ev(iorb))/am(j)*arc				! JRT
     ar(j) = arc + amc0*(fa1-fa0)
     br(j) = brc + amc0*(fb1-fb0)
     fa1 = br(j)
     fb1 = b*br(j) + term1(j,am,amd,bm,cm,r,nn,lp,rab)*br(j) + &			! JRT
   &                 term0(j,am,amd,bm,cm,r,nn,lp,rab)*ar(j) + &			! JRT
   &                 rab(j)*rab(j)*(v(j)-ev(iorb))/am(j)*ar(j)				! JRT

     cr(j) =         term1(j,am,amd,bm,cm,r,nn,lp,rab)*br(j) + &			! JRT
   &                 term0(j,am,amd,bm,cm,r,nn,lp,rab)*ar(j) + &			! JRT
   &                 rab(j)*rab(j)*(v(j)-ev(iorb))/am(j)*ar(j)				! JRT

!     count nodes

     if (ar(j)*ar(jm1) .le. 0) nodes = nodes + 1
   enddo

!     end outward integration

!     if incorrect number of nodes modify energy stepwise

   if (nodes .gt. no(iorb)-lo(iorb)-1) then 		! too many nodes -- decrease ev
     if (ev(iorb) .lt. emax) emax = ev(iorb)
     if (devmax .gt. 0.D0) devmax = -devmax * 0.5d0
     ev(iorb) = ev(iorb) + devmax
     goto 190
   elseif (nodes .lt. no(iorb)-lo(iorb)-1) then		! too few nodes -- increase ev
     if (ev(iorb) .gt. emin) emin = ev(iorb)
     if (devmax .lt. 0.D0) devmax = -devmax * 0.5d0
     ev(iorb) = ev(iorb) + devmax
     goto 190
   endif

!     correct number of nodes

   arctp = ar(nctp)
   brctp = br(nctp)

!     inward integration from ninf to nctp -- startup

   if (ninf .lt. nwell) then

!     normal startup

     istart = ninf - nctp + 1
     if (istart .gt. 5) istart = 5
     do jj=1,istart
       j = ninf-jj+1
       alf = bm(i)*lp*(lp-1)/r2(j) + v(j) - ev(iorb)					! JRT ***
       if (alf .lt. 0.D0) alf = 0.D0
       alf = sqrt(alf)
       ar(j) = exp(-alf*r(j))
       br(j) = -rab(j)*alf*ar(j)
       cr(j) = +rab(j)*rab(j)*alf**2*ar(j)						! JRT
     enddo
     fa5 = br(ninf)
     fa4 = br(ninf-1)
     fa3 = br(ninf-2)
     fa2 = br(ninf-3)
     fa1 = br(ninf-4)
     fb5 = b*br(ninf) + term1(ninf,am,amd,bm,cm,r,nn,lp,rab)*br(ninf) +         &		! JRT
   &                    term0(ninf,am,amd,bm,cm,r,nn,lp,rab)*ar(ninf) +         &		! JRT
   &                    rab(ninf)*rab(ninf)*(v(ninf)-ev(iorb))/am(ninf)*ar(ninf)			! JRT
     fb4 = b*br(ninf-1) + term1(ninf-1,am,amd,bm,cm,r,nn,lp,rab)*br(ninf-1) +   &		! JRT
   &                      term0(ninf-1,am,amd,bm,cm,r,nn,lp,rab)*ar(ninf-1) +   &		! JRT
   &                      rab(ninf-1)*rab(ninf-1)*(v(ninf-1)-ev(iorb))/am(ninf-1)*ar(ninf-1)	! JRT
     fb3 = b*br(ninf-2) + term1(ninf-2,am,amd,bm,cm,r,nn,lp,rab)*br(ninf-2) +   &		! JRT
   &                      term0(ninf-2,am,amd,bm,cm,r,nn,lp,rab)*ar(ninf-2) +   &		! JRT
   &                      rab(ninf-2)*rab(ninf-2)*(v(ninf-2)-ev(iorb))/am(ninf-2)*ar(ninf-2)	! JRT
     fb2 = b*br(ninf-3) + term1(ninf-3,am,amd,bm,cm,r,nn,lp,rab)*br(ninf-3) +   &		! JRT
   &                      term0(ninf-3,am,amd,bm,cm,r,nn,lp,rab)*ar(ninf-3) +   &		! JRT
   &                      rab(ninf-3)*rab(ninf-3)*(v(ninf-3)-ev(iorb))/am(ninf-3)*ar(ninf-3)	! JRT
     fb1 = b*br(ninf-4) + term1(ninf-4,am,amd,bm,cm,r,nn,lp,rab)*br(ninf-4) +   &		! JRT
   &                      term0(ninf-4,am,amd,bm,cm,r,nn,lp,rab)*ar(ninf-4) +   &		! JRT
   &                      rab(ninf-4)*rab(ninf-4)*(v(ninf-4)-ev(iorb))/am(ninf-4)*ar(ninf-4)	! JRT
! A.y" + [ 2.A.n.1/r + A'+2C ].y' + [ A.n(n-1)-B.l(l+1) ].1/r^2.y + [A'+2C].(n-1).1/r.y - [C'+v-e].y = 0

     istart = 5
   else

!     power series startup

!     find (v(rwell) - ev(iorb) / 6

     if (nwell .lt. nr) then
       alf = (v(nwell)*(r(nwell+1)-rwell)           &
        &     - v(nwell+1)*(r(nwell)-rwell))        &
        &     / (r(nwell+1) - r(nwell))
     else
       alf = v(nr)
     endif
     alf = (alf - ev(iorb)) * .16666666666666666666666666666666d0

     dnwell=dble(nwell)
     rj = a * (exp(b*(dnwell-1.d0)) - 1.d0)
     rabj = (rj + a) * b
     rj = rj - rwell
     fa1 = -(1.d0 + 3.d0 * alf * rj**2) * rabj
     fb1 = fa1 * b - 6.d0 * alf * rj * rabj * rabj
     ar(nwell) = -rj * (1.d0 + alf * rj**2)
     br(nwell) = fa1

     rj = a * (exp(b*(dnwell)) - 1.d0)
     rabj = (rj + a) * b
     rj = rj - rwell
     fa2 = -(1.d0 + 3.d0 * alf * rj**2) * rabj
     fb2 = fa2 * b - 6.d0 * alf * rj * rabj * rabj

     rj = a * (exp(b*(dnwell+1.d0)) - 1.d0)
     rabj = (rj + a) * b
     rj = rj - rwell
     fa3 = -(1.d0 + 3.d0 * alf * rj**2) * rabj
     fb3 = fa3 * b - 6.d0 * alf * rj * rabj * rabj

     rj = a * (exp(b*(dnwell+2.d0)) - 1.d0)
     rabj = (rj + a) * b
     rj = rj - rwell
     fa4 = -(1.d0 + 3.d0 * alf * rj**2) * rabj
     fb4 = fa4 * b - 6.d0 * alf * rj * rabj * rabj

     rj = a * (exp(b*(dnwell+3.d0)) - 1.d0)
     rabj = (rj + a) * b
     rj = rj - rwell
     fa5 = -(1.d0 + 3.d0 * alf * rj**2) * rabj
     fb5 = fa5 * b - 6.d0 * alf * rj * rabj * rabj

     istart = 1
   endif

!     integration loop

   istop = ninf - nctp
   do jj=istart,istop
     j = ninf - jj

!     predictor (Adams-Bashforth)

     arp = ar(j+1) - (abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5)
     brp = br(j+1) - (abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5)
     fa0 = brp
     fb0 = b*brp + term1(j,am,amd,bm,cm,r,nn,lp,rab)*brp + &				! JRT
    &              term0(j,am,amd,bm,cm,r,nn,lp,rab)*arp + &				! JRT
    &              rab(j)*rab(j)*(v(j)-ev(iorb))/am(j)*arp				! JRT

!     corrector (Adams-Moulton)

     arc = ar(j+1) - (amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4)
     brc = br(j+1) - (amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4)
     fa5 = fa4
     fb5 = fb4
     fa4 = fa3
     fb4 = fb3
     fa3 = fa2
     fb3 = fb2
     fa2 = fa1
     fb2 = fb1
     fa1 = brc
     fb1 = b*brc + term1(j,am,amd,bm,cm,r,nn,lp,rab)*brc + &				! JRT
    &              term0(j,am,amd,bm,cm,r,nn,lp,rab)*arc + &				! JRT
    &              rab(j)*rab(j)*(v(j)-ev(iorb))/am(j)*arc				! JRT

     ar(j) = arc - amc0*(fa1-fa0)
     br(j) = brc - amc0*(fb1-fb0)
     fa1 = br(j)
     fb1 = b*br(j) + term1(j,am,amd,bm,cm,r,nn,lp,rab)*br(j) + &			! JRT
   &                 term0(j,am,amd,bm,cm,r,nn,lp,rab)*ar(j) + &			! JRT
   &                 rab(j)*rab(j)*(v(j)-ev(iorb))/am(j)*ar(j)				! JRT

     cr(j) =         term1(j,am,amd,bm,cm,r,nn,lp,rab)*br(j) + &			! JRT
   &                 term0(j,am,amd,bm,cm,r,nn,lp,rab)*ar(j) + &			! JRT
   &                 rab(j)*rab(j)*(v(j)-ev(iorb))/am(j)*ar(j)				! JRT

   enddo

!     end inward integration

!     rescale ar and br outside nctp to match
!     ar(nctp) from outward integration

   factor = arctp/ar(nctp)
   do j=nctp,ninf
     ar(j) = factor * ar(j)
     br(j) = factor * br(j)
     cr(j) = factor * cr(j)								! JRT
   enddo

!     find normalization

   factor = 0.D0
   ll = 4
   do j=2,ninf
     factor = factor + ll*ar(j)*ar(j)*rab(j)*r(j)**( 2*nn )				! JRT
     ll = 6 - ll
   enddo
   factor = factor *0.3333333333333333333333333333333333333333d0

!     modify eigenvalue ev

   dev = r(nctp)**(2*nn-1)*ar(nctp) * (r(nctp)*(brctp-br(nctp)) + nn*(arctp-ar(nctp)))	! JRT
   dev = dev /( factor*rab(nctp) )							! JRT

!     resort to bisection if dev too large

   if (abs(dev) .gt. abs(devmax)) then
     if (devmax*dev .lt. 0.D0) devmax = -devmax * 0.5d0
     dev = devmax
   endif
   evold = ev(iorb)
   ev(iorb) = ev(iorb) + dev
   if (ev(iorb) .gt. emax) ev(iorb) = (evold + emax) * 0.5d0
   if (ev(iorb) .lt. emin) ev(iorb) = (evold + emin) * 0.5d0
   if (abs(dev) .lt. tol*(1+abs(ev(iorb)))) goto 220

 190  continue

!     eigenpar not converged in itmax iterations

!     if missing -- find normalization

   if (nodes .ne. no(iorb)-lo(iorb)-1) then
     factor = 0.D0
     ll = 4
     do j=2,ninf
       factor = factor + ll*ar(j)*ar(j)*rab(j)*r(j)**( 2*nn )				! JRT
       ll = 6 - ll
     enddo
   endif
   factor = factor * 0.33333333333333333333333333333333333333d0

!     error message

 210 format(' orb #',i3,' did not converge',/,            &
    &     ' ev =',e18.10,' nodes =',i2,' dev =',e18.10)

!     normalize wavefunction and change br from d(ar)/dj to d(ar)/dr

 220 factor = 1 / sqrt(factor)
   do j=1,ninf
     ar(j) = factor*ar(j)
     br(j) = factor*br(j) / rab(j)
     cr(j) = factor*cr(j) / (rab(j)*rab(j))						! JRT
   enddo

!     convert y to u=r^n.(1 + aa.r + bb.r^2 + ...) [ = r^{l+1}.(1 + aa.r + bb.r^2 + ...) for std. SE ]

   y=ar ; yd=br ; ydd=cr
   do j=1,ninf										! JRT
     cr(j) = r(j)**(nn-2.d0)*(r(j)**2*cr(j)+2.d0*nn*r(j)*br(j)+nn*(nn-1.d0)*ar(j))	! JRT
     br(j) = r(j)**(nn-1.d0)*(r(j)*br(j)+nn*ar(j))					! JRT
     ar(j) = ar(j)*r(j)**nn								! JRT
   enddo										! JRT

   ev=0.5d0*ev ; v=0.5d0*v				! rydbergs to hartrees

!write(*,*) 'l=                        :',lp-1						! VTEMP
!write(*,*) 'lp=                       :',lp						! VTEMP
!write(*,*) 'aa=                       :',aa						! VTEMP
!write(*,*) 'bb=                       :',bb						! VTEMP
!write(*,*) 'n=                        :',nn						! VTEMP
!write(*,*) 'nctp , ninf , nwell , nr  :',nctp,ninf,nwell,nr				! VTEMP
!write(*,*) 'dev                       :',abs(dev)					! VTEMP
!write(*,*) 'iterations                :',i						! VTEMP

 return

 END SUBROUTINE difnrl_gen




 REAL(kind=kind(1.d0)) FUNCTION term1(i,am,amd,bm,cm,r,nn,lp,rab)
 IMPLICIT NONE
 INTEGER i,lp
 REAL(kind=kind(1.d0)) am(:),amd(:),bm(:),cm(:),r(:),nn,rab(:)

!- A.y" + [ 2.A.n.1/r + A'+2C ].y' + [ A.n(n-1)-B.l(l+1) ].1/r^2.y + [A'+2C].(n-1).1/r.y - [C'+v-e].y = 0

 term1 = -rab(i)*( 2.d0*am(i)*nn/r(i) + amd(i) + cm(i) )/am(i)

 END FUNCTION


 REAL(kind=kind(1.d0)) FUNCTION term0(i,am,amd,bm,cm,r,nn,lp,rab)
 IMPLICIT NONE
 INTEGER i,lp
 REAL(kind=kind(1.d0)) am(:),amd(:),bm(:),cm(:),r(:),nn,rab(:)
 REAL(kind=kind(1.d0)) t1,t2,t3,t4

!- A.y" + [ 2.A.n.1/r + A'+2C ].y' + [ A.n(n-1)-B.l(l+1) ].1/r^2.y + [A'+2C].(n-1).1/r.y - [C'+v-e].y = 0

 t1= ( am(i)*nn*(nn-1.d0)-bm(i)*(lp-1)*lp )/r(i)**2
 t2= ( amd(i) + cm(i) )*(nn-1.d0)/r(i)

 term0= -rab(i)*rab(i)*(t1+t2)/am(i)

 END FUNCTION


END MODULE ode

