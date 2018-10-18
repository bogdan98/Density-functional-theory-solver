MODULE int_interp
IMPLICIT NONE

CONTAINS


!--- Find f=0 of tabulated f.
!--- Supplied index i0 points at region (i0,i0+1) that brackets f(r)=0.
!--- Supplied index i1 points at index about which interpolation is performed.
!--- i1=i0 or i1=i0+1 are equally justified...
 SUBROUTINE fsolve(rmid,i0,i1,f,r,nr)
 IMPLICIT NONE
 INTEGER :: i0,i1,nr
 INTEGER :: j
 REAL(kind=kind(1.d0)), PARAMETER :: xacc=1.d-14
 REAL(kind=kind(1.d0)) :: r(:),f(:),rmid
 REAL(kind=kind(1.d0)) :: r1,r2,f1,f2,rtbis,dr
 REAL(kind=kind(1.d0)) :: t1,t2,t3

 r1=r(i0) ; r2=r(i0+1)
 call dinterp(r1,i1,r,f,nr,f1,t2,t3) ; call dinterp(r2,i1,r,f,nr,f2,t2,t3)
 if(f1.lt.0.)then
   rtbis=r1
   dr=r2-r1
 else
   rtbis=r2
   dr=r1-r2
 endif
 do j=1,50
   dr=dr*.5
   rmid=rtbis+dr
   call dinterp(rmid,i1,r,f,nr,f2,t2,t3)
   if(f2.le.0.)rtbis=rmid
   if(abs(dr).lt.xacc .or. f2.eq.0.) exit
 enddo

 END SUBROUTINE fsolve


!--- Integrate from ir1->ir2.
 REAL(kind=kind(1.d0)) FUNCTION norm00(ar,r,rab,ir1,ir2)
 IMPLICIT NONE
 !----------------------------------------------------------------------!
 !  Integrate ar r(n1) to r(n2).                                        !
 !----------------------------------------------------------------------!
 real(kind=kind(1.d0)) ar(:),r(:),rab(:)
 integer ir1,ir2
 real(kind=kind(1.d0)) valp
 integer i,j,k,irr,ll,none,ntwo,nthree

 ll = 2
 valp = -ar(ir2)*rab(ir2)
 irr = ir2 - ir1 + 1
 if (irr .ne. 2*(irr/2)) then
   do k = ir2, ir1, -1
     valp  = valp + ll*ar(k)*rab(k)
     ll = 6 - ll
   enddo
   valp  = valp - ar(ir1)*rab(ir1)
 else
   nthree = ir1 + 3
   ntwo = ir1 + 2
   none = ir1 + 1
   do k = ir2, nthree, -1
     valp = valp + ll*ar(k)*rab(k)
     ll = 6 - ll
   enddo
   valp = valp - ar(nthree)*rab(nthree)
   valp = valp + 9*( ar(ir1)*rab(ir1)+  &
  &       3*ar(none)*rab(none)+        &
  &       3*ar(ntwo)*rab(ntwo)+        &
  &       ar(nthree)*rab(nthree))/8
 endif
 norm00=valp/3.d0

 END FUNCTION norm00


!--- Wrapper for polintd. Interpolates y using pn points centered on i.
!--- Provides value at supplied r0, so can extrapolate or interpolate.
!--- Adjust for i at start/end of array.
 SUBROUTINE dinterp(r0,i,r,x,nr,y,yd,ydd)
 IMPLICIT NONE
 INTEGER :: i,ii,nr
 INTEGER, PARAMETER :: pn=5					! MUST be odd.
 REAL(kind=kind(1.d0)) :: r0,r(:),x(:),y,yd,ydd

 if(i.le.pn/2)then
  ii=1
 elseif(i.ge.nr+1-pn/2)then
  ii=nr+1-pn
 else
  ii=i-pn/2
 endif
 call polintd( r(ii:ii+pn-1), x(ii:ii+pn-1), pn, r0, y,yd,ydd)

 END SUBROUTINE dinterp


 SUBROUTINE polintd(xa,ya,n,x,y,y_d,y_dd)
 IMPLICIT NONE
 INTEGER nmax,n,i,m,ns
 PARAMETER (nmax=10)
 REAL(kind=kind(1.d0)) dy,dy_d,dy_dd,x,y,y_d,y_dd,xa(n),ya(n), &
  &                    den,dif,dift,ho,hp,w,w_d,w_dd,          &
  &                    c(nmax),d(nmax),c_d(nmax),d_d(nmax),    &
  &                    c_dd(nmax),d_dd(nmax)

 ns = 1
 dif = abs(x-xa(1))

!---  ns is index of closest table entry.

 do i = 1,n

   dift = abs(x-xa(I))
   if (dift.lt.dif) then
     ns  = i
     dif = dift
   endif
!--- Initializing the tableau of entries of c's and d's.
   c   (i) = ya(i)
   d   (i) = ya(i)
   c_d (i) = 0.d0
   d_d (i) = 0.d0
   c_dd(i) = 0.d0
   d_dd(i) = 0.d0

 enddo


!---   This is the initial approximation to y.
 y   = ya(ns)
 y_d = 0.d0
 y_dd= 0.d0
 ns  = ns - 1

!---  For each column of the tableu, loop over the current c's and d's
!---  and update them.

 do m = 1,n - 1

   do i = 1,n - m
     ho   = xa(i) - x
     hp   = xa(i+m) - x
     w    = c   (i+1) - d   (i)
     w_d  = c_d (i+1) - d_d (i)
     w_dd = c_dd(i+1) - d_dd(i)
     den  = ho - hp
!---  Here the c's and d's are updated.
     d   (i) = hp*w/den
     c   (i) = ho*w/den
     d_d (i) = (-w + hp*w_d)/den
     c_d (i) = (-w + ho*w_d)/den
     d_dd(i) = (-2.d0*w_d + hp* w_dd)/den
     c_dd(i) = (-2.d0*w_d + ho* w_dd)/den
   enddo

!---  After each column in the tableau is completed, decide which
!---  correction, c or d needs to be added to the accumulating value
!---  of y i.e. which path to take through the tableau-forking  up or
!---  down. This is done in such a way as to take the most "stright-line"
!---  route through the tableau to its apex, updating ns accordingly, to
!---  keep track of where we are. This route keeps the partial approx
!---  centered on the target x. The last dy added is the error indication.

   if (2*ns.lt.n-m) then
     dy   = c   (ns+1)
     dy_d = c_d (ns+1)
     dy_dd= c_dd(ns+1)
   else
     dy    = d   (ns)
     dy_d  = d_d (ns)
     dy_dd = d_dd(ns)
     ns = ns - 1
   endif

   y    = y    + dy
   y_d  = y_d  + dy_d
   y_dd = y_dd + dy_dd

 enddo

 END SUBROUTINE polintd



END MODULE int_interp
