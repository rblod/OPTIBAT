module m_sort2

contains
subroutine sort2(n,arr,brr)

! Sorts an array arr(1:n) into ascending numerical order, by
! straight insertion, while making the corresponding rearrangement of the array brr(1:n). 


   implicit none

   integer, intent(in)  :: n 
   real, intent(inout),dimension(n) :: arr,brr

  integer i,j 
  real a,b 

  do j=2,n
    a=arr(j) 
    b=brr(j) 
    do i=j-1,1,-1 
      if(arr(i).le.a) goto 10
      arr(i+1)=arr(i) 
      brr(i+1)=brr(i) 
    enddo 
    i=0
 10 arr(i+1)=a 
    brr(i+1)=b 
  enddo 

end subroutine sort2

RECURSIVE SUBROUTINE quick_sort(list, order)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.

IMPLICIT NONE
REAL, DIMENSION (:), INTENT(IN OUT)  :: list
INTEGER, DIMENSION (:), INTENT(OUT)  :: order

! Local variable
INTEGER :: i

DO i = 1, SIZE(list)
  order(i) = i
END DO

CALL quick_sort_1(1, SIZE(list))

CONTAINS

RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL                :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END SUBROUTINE quick_sort_1


SUBROUTINE interchange_sort(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL                :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort


   subroutine search_value(z,x,n,ind)
   implicit none
   integer::n,ind
   real,dimension(n)::z
   real::x
   integer::deb,fin,mil
   
   deb=1
   fin=n

   do while(deb.le.fin)
     mil=(deb+fin)/2
     if(x==z(mil))then
       ind=mil
       return
     endif
     if(x.lt.z(mil))then
       fin=mil-1
     else
       deb=mil+1
     endif
     ind=fin
   enddo
  
   end subroutine


end module m_sort2
