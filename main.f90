MODULE commonvars
USE parameters
IMPLICIT NONE

!--- For radial grid and tabulation:
REAL(kind=kind(1.d0)) :: a,b,aaa,bbb,rmax
REAL(kind=kind(1.d0)),allocatable :: u(:),ud(:),udd(:), density(:), V2(:), V3(:)
REAL(kind=kind(1.d0)),allocatable :: r(:),rab(:), energies(:)
REAL(kind=kind(1.d0)),allocatable :: v1(:), ediff, enguess
REAL(kind=kind(1.d0)),allocatable :: el(:),el1(:)
REAL(kind=kind(1.d0)),allocatable :: ppwf1(:,:),dppwf1(:,:),ddppwf1(:,:)
REAL(kind=kind(1.d0)),allocatable :: wk(:)
INTEGER :: nr

!--- For parameterised ppots:
REAL(kind=kind(1.d0)) :: a_coef(maxlp,maxpower,maxterm),b_coef(maxlp,maxpower,maxterm)
REAL(kind=kind(1.d0)) :: zeff, aalpha, rrl, alpha_coef, rl_coef
REAL(kind=kind(1.d0)) :: rl0(maxlp), rc(maxlp)
INTEGER :: nc(maxlp),no(maxlp),lo(maxlp),nog,lpmax
INTEGER :: v_inc(maxlp,maxpower,maxterm),v_var(maxlp,maxpower,maxterm),alpha_var,rl_var	! Unused - from flags in pp.start
CHARACTER(24) :: line1,char1,input,atm

INTEGER :: ic,i0

END MODULE commonvars



PROGRAM main
USE parameters
USE commonvars
USE io
USE ode
USE int_interp
USE orbital
USE hartree
USE excorr
USE correction
USE DFTSolver
IMPLICIT NONE

!--- For general use:
REAL(kind=kind(1.d0)) :: t0,t1,t2, ekin, ecoul, enuc, exc
INTEGER :: i,j,k,l,lp,kk,ii, iterations


!--- Define gradial grid, read in gaussian-parameterised norm-conserving pseudopotential, and tabulate.
!--- Given here as an example on how to define the radial grid, parameters of which are used as inputs in 
!--- DFT solver subroutines. The part below that reads from the file can be deleted, since it is irrelevant to 
!--- DFT calculations. However, radial grid parameters still need to be defined
 input='pp.start'
 write(*,*) 'Reading parameterised ppot...'
 call readpp(input,line1,a_coef,b_coef,alpha_coef,rl_coef,v_inc,v_var,alpha_var,  &
   &         rl_var,rc,rl0,nog,lpmax)
 read(line1,*) char1,atm
 call reportpp(line1,a_coef,b_coef,alpha_coef,rl_coef,v_inc,v_var,alpha_var,rl_var,nog,lpmax)
 do i=1,lpmax
  lo(i)=i-1 ; no(i)=lo(i)+1
 enddo
 zeff=a_coef(lpmax,1,1)
 write(*,*)

 rmax=100.d0
 aaa=8.d0 ; bbb=75.d0 								! Grid parameters.
 a=bbb*exp(-aaa*log(10.d0))/10.0 !divide by nuclear charge of the element of interest 
 b=1.d0/bbb
 write(*,*) 'Grid spacing at r=0:    ',a*b
 write(*,*) 'Grid spacing at r=rmax: ',(rmax+a)*b
 do i=1,10000
  t1=a*(exp(b*dble(i-1))-1.d0)
  if(t1>rmax)exit
 enddo
 nr=i

 allocate(r(nr),rab(nr), V2(nr), V3(nr), energies(20))
 allocate(u(nr),ud(nr),udd(nr),wk(nr),el(maxlp))
 allocate(el1(maxlp),v1(nr), density(nr))
 allocate(ppwf1(nr,lpmax),dppwf1(nr,lpmax),ddppwf1(nr,lpmax))

 do i=1,nr
   r(i)=a*(exp(b*dble(i-1))-1.d0) ; rab(i)=(r(i)+a)*b		! r & dr/di
 enddo

 do i=1,lpmax									! Adjust rc(s) to grid.
   do j=1,nr ; if(r(j)>rc(i))exit ; enddo ; nc(i)=j-1 ; rc(i)=r(nc(i))
 enddo

!--- Solve SE. [ Centrifugal 'potential' is included in v supplied. ]
do lp=1,lpmax
  u=0.d0 ; ud=0.d0 ; l=lo(lp) ; el1(lp)=0.5d0*v1(1) ; t1=v1(1)
  call difnrl(lp,v1(:)+0.5d0*dble(l*(l+1))/r(:)**2,              &
   &         0.d0,u,ud,nr,a,b,r,rab,lpmax,no,lo,1.d0,                 &
   &         rmax,el1)
 ppwf1(:,lp)=u(:) ; dppwf1(:,lp)=ud(:)
enddo

!Starting from here, DFT calculations are done


enguess = 0.d0;
ekin = 0.d0;
ecoul = 0.d0;
enuc = 0.d0;
exc= 0.d0;
ediff = 1.d0;
density(:) = (4.0/(sqrt(3.14159265359)))*exp(-r(:)**2);
v1(1) = 0.d0;
v1(2:nr) = -10.0/r(2:nr);
V2 = 0.d0;
V3 = 0.d0;
open (unit = 254, file = "densitytest")
do i = 1, nr
  write(254, *) r(i), density(i)
end do
!this is an example for Neon
call iterates((/1, 1, 2, 2, 2, 2, 2, 2, 2, 2/), (/0, 0, 0, 0, 1, 1, 1, 1, 1, 1/),&
& 10, density, &
&v1, r, rab, nr, u, ud, a, b, rmax, energies, enguess, 10.d0,  ecoul, enuc, exc, &
& iterations);
open (unit =  257, file = "AtomicData")
do i = 1, 10
  write(257,*), energies(i)
end do
write(257,*), enguess
write(257,*), enguess - ecoul - enuc - exc ;
write(257,*), ecoul
write(257,*), enuc
write(257,*), exc
print*, iterations

!DFT Calculations end here
 deallocate(r,rab,u,ud,udd,wk,el)
 deallocate(el1,v1, density, V2)
 deallocate(ppwf1,dppwf1,ddppwf1)

END PROGRAM main
