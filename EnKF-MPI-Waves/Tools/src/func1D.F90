real function func1D(sigg,length,dx,n1,n2)
implicit none
real sum1,sum2,sigg,length
real sigg2,pi2,kappa,kappa2,lambda,lambda2,dx
integer l,p,n1,n2
real, parameter :: pi=3.141592653589

sigg2=sigg**2

pi2=2.0*pi
kappa=pi2/(float(2*n1)*dx)
kappa2=kappa**2

!      /* Calculate sum1 */
sum1=0.0
do l=-n1/2+1,n1/2
   sum1=sum1+&
   exp(-2.0*(kappa2*float(l*l))/sigg2)*cos(kappa*float(l)*length)
enddo

!      /* Calculate sum2 */
sum2=0.0
do l=-n1/2+1,n1/2
   sum2=sum2+exp(-2.0*(kappa2*float(l*l))/sigg2)
enddo

func1D = sum1/sum2 - exp(-1.0)

end function func1D
