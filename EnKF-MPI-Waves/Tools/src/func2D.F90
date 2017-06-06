real function func2D(sigma,length,dx,n1,n2)
! Function used to calculate $sigma$ and $c$.
   implicit none
   real sum1,sum2,sigma,length
   real sigma2,pi2,kappa,kappa2,lambda,lambda2,dx
   integer l,p,n1,n2
   real, parameter :: pi=3.141592653589
   real :: tmpexponent

   sigma2=sigma**2

   pi2=2.0*pi
   kappa=pi2/(float(n1)*dx)
   kappa2=kappa**2

   lambda=pi2/(float(n2)*dx)
   lambda2=lambda**2


   ! Calculate sum1
   sum1=0.0
   do p=-n2/2+1,n2/2
   do l=-n1/2+1,n1/2
      !sum1=sum1+exp(-2.0*(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)&
      !*cos(kappa*float(l)*length)

      tmpexponent=-2.0*(kappa2*float(l*l)+lambda2*float(p*p))/sigma2
      tmpexponent=max(-1e3,tmpexponent)
      sum1=sum1+exp(tmpexponent)*cos(kappa*float(l)*length)
   enddo
   enddo

   ! Calculate sum2
   sum2=0.0
   do p=-n2/2+1,n2/2
   do l=-n1/2+1,n1/2
      !sum2=sum2+exp(-2.0*(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)

      tmpexponent=-2.0*(kappa2*float(l*l)+lambda2*float(p*p))/sigma2
      tmpexponent=max(-1e3,tmpexponent)
      sum2=sum2+exp(tmpexponent)
   enddo
   enddo

   func2D = sum1/sum2 - exp(-1.0)

end function func2D
