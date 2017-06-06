module m_fixsample
! This routines removes the ensemble mean and scales its variance to 1

contains
subroutine fixsample(A,nx,ny,nrens)

   implicit none
   integer, intent(in) :: nx, ny, nrens
   real, intent(inout) :: A(nx,ny,nrens)

   real mean(nx,ny), variance        ! mean field, global variance
   integer n,k

!Sums
   mean = sum(A,dim=3)
   variance =  sum(A**2.) - sum(mean**2.)/real(nx*ny*nrens)

!Scales
   mean = mean / real(nrens)
   if (nrens.LE.1) stop 'm_fixsample: less than 1 member'
   variance  = variance / real(nx*ny*nrens-1)
   !print*, ' Sample variance', variance

!Adjusts
   do n = 1, nrens
      A(:,:,n) = A (:,:,n)- mean
   enddo
   A = A / sqrt(variance)

end subroutine fixsample
end module m_fixsample
