MODULE orbital
USE int_interp, only : dinterp , norm00
IMPLICIT NONE

CONTAINS


 SUBROUTINE wf_ytou(nn,r,ppy,dppy,ddppy,ppwf,dppwf,ddppwf,nr,lpmax)
 IMPLICIT NONE
 REAL(kind=kind(1.d0)) :: n,nn(3),r(:)
 REAL(kind=kind(1.d0)) :: ppy(:,:),dppy(:,:),ddppy(:,:),ppwf(:,:),dppwf(:,:),ddppwf(:,:)
 INTEGER :: i,j,nr,lpmax
!--- Construct u orbitals from y orbitals, using supplied details.
!--- Note it is better to get u" from the SE or GSE.
!---  u =                                              r^n.y
!---  u'= r^n.y' + n.r^{n-1}.y                       = r^{n-1}.[ r.y' + n.y ]
!---  u"= r^n.y" + 2.n.r^{n-1).y' + n(n-1).r^{n-2}.y = r^{n-2}.[ r^2.y" + 2.n.r.y' + n(n-1).y ]
 do j=1,lpmax
   n=nn(j)
   do i=1,nr
      ppwf(i,j) = r(i)**n*ppy(i,j)
     dppwf(i,j) = r(i)**(n-1.d0)*( r(i)*dppy(i,j) + n*ppy(i,j) )
    ddppwf(i,j) = r(i)**(n-2.d0)*( r(i)**2*ddppy(i,j) + 2.d0*n*r(i)*dppy(i,j) + n*(n-1.d0)*ppy(i,j) )
   enddo
 enddo
 END SUBROUTINE wf_ytou

 SUBROUTINE wf_utoy(nn,r,ppy,dppy,ddppy,ppwf,dppwf,ddppwf,nr,lpmax)
 IMPLICIT NONE
 REAL(kind=kind(1.d0)) :: n,nn(3)
 REAL(kind=kind(1.d0)) :: ppy(:,:),dppy(:,:),ddppy(:,:),ppwf(:,:),dppwf(:,:),ddppwf(:,:),r(:)
 INTEGER :: i,j,nr,lpmax
!--- Construct y orbitals from u orbitals, using supplied details.
!--- Note it is better to get y" from the SE or GSE.
!---  y  =                                                       u/r^n
!---  y' = u'/r^{n} - n.u.1/r^{n+1}	                       = [ r.u' - n.u ]/r^{n+1}
!---  y" = u"/r^{n} - 2.n.u'.1/r^{n+1} - n.u.(-(n+1))1/r^{n+2} = [r^2.u" - 2.n.r.u' + n(n+1).u]/r^{n+2}
 do j=1,lpmax
   n=nn(j)
   do i=1,nr
      ppy(i,j) = ppwf(i,j)/r(i)**n
     dppy(i,j) = ( r(i)*dppwf(i,j) - n*ppwf(i,j) )/r(i)**(n+1)
    ddppy(i,j) = ( r(i)**2*ddppwf(i,j) - 2.d0*n*r(i)*dppwf(i,j) + n*(n+1.d0)*ppwf(i,j) )/r(i)**(n+2)
   enddo
 enddo
 END SUBROUTINE wf_utoy


END MODULE orbital

