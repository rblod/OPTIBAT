module m_zeroin
contains
subroutine zeroin(func,zeropkt,ax,bx,tol,length,dx,fval,n1,n2)
! Finds zero of function f.
! A zero of the function  $func(x,length,dx,n1,n2)$
! is computed in the interval $[ax,bx]$.
! Zeroin| returns a zero $x$ in the given interval
! to within a tolerance  $4*macheps*abs(x) + tol$, where macheps
! is the relative machine precision.
! NOTA: do not trap underflow during compiling (error)

! This function subprogram is a slightly  modified  translation  of
! the algol 60 procedure  zero  given in  richard brent, algorithms for
! minimization without derivatives, prentice - hall, inc. (1973).

   real, external :: func

   integer n1,n2
   real zeropkt,length,dx
   real ax   ! left endpoint of initial interval
   real bx   ! right endpoint of initial interval
   real tol  !  desired length of the interval of uncertainty of the
   real  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
   real  abs,sign,fval

!  compute eps, the relative machine precision

#ifdef DEBUG
   write(*,*)'A=',ax,bx,tol,length,dx
#endif
   icorr=0

   eps = 1.0
 10 eps = eps/2.0
   tol1 = 1.0 + eps
   if (tol1 .gt. 1.0) go to 10

#ifdef DEBUG
   print*, 'Zeroin: done computing eps=epsilon() '
#endif

! initialization
 77 a = ax
   b = bx
   fa = func(a,length,dx,n1,n2)
   fb = func(b,length,dx,n1,n2)


   if (fa*fb.gt.0.0) then
#ifdef DEBUG
      write(*,*)'fa=',fa
      write(*,*)'fb=',fb
      write(*,*)'fa*fb =',fa*fb,'is greater than zero'
#endif
      ax=0.1*ax
      bx=10.0*bx
      icorr=icorr+1
      if (icorr < 20) then
         goto 77
      else
         write(*,'(2(a,g13.5))')'zeroin: No convergence, ax=',ax,' bx=',bx
         stop
      endif
   endif

! begin step

 20 c = a
   fc = fa
   d = b - a
   e = d
 30 if (abs(fc) .ge. abs(fb)) go to 40
   a = b
   b = c
   c = a
   fa = fb
   fb = fc
   fc = fa

! convergence test

 40 tol1 = 2.0*eps*abs(b) + 0.5*tol
   xm = .5*(c - b)
   if (abs(xm) .le. tol1) go to 90
   if (fb .eq. 0.0) go to 90

! is bisection necessary

   if (abs(e) .lt. tol1) go to 70
   if (abs(fa) .le. abs(fb)) go to 70

! is quadratic interpolation possible

   if (a .ne. c) go to 50

! linear interpolation

   s = fb/fa
   p = 2.0*xm*s
   q = 1.0 - s
   go to 60

! inverse quadratic interpolation

 50 q = fa/fc
   r = fb/fc
   s = fb/fa
   p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
   q = (q - 1.0)*(r - 1.0)*(s - 1.0)

! adjust signs

 60 if (p .gt. 0.0) q = -q
   p = abs(p)

! is interpolation acceptable

   if ((2.0*p) .ge. (3.0*xm*q - abs(tol1*q))) go to 70
   if (p .ge. abs(0.5*e*q)) go to 70
   e = d
   d = p/q
   go to 80

! bisection

 70 d = xm
   e = d

! complete step

 80 a = b
   fa = fb
   if (abs(d) .gt. tol1) b = b + d
   if (abs(d) .le. tol1) b = b + sign(tol1, xm)
   fb = func(b,length,dx,n1,n2)
   if ((fb*(fc/abs(fc))) .gt. 0.0) go to 20
   go to 30

! done

 90 zeropkt = b
   fval=func(b,length,dx,n1,n2)
end subroutine zeroin
end module m_zeroin
