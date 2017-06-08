module m_ensemble_stats
   implicit none
   contains
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! fields to be analyzed!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer function get_nrfields_stats(infile)
      implicit none 
      character(len=*):: infile
      integer :: ios,first,last
      logical :: ex
      character(len=3) :: char3

      inquire(exist=ex,file=infile)
      if (.not. ex) then
        print *,'Could not find '//infile
      end if

      open(103,status='old',form='formatted',file=infile)
      ios=0
      get_nrfields_stats=0
      do while (ios==0)
        read(103,100,iostat=ios) char3
        if (ios==0) get_nrfields_stats=get_nrfields_stats+1
      end do
      close(103)
  100 format (a3)
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
   
      subroutine get_fields_stats(infile,numfields,fieldsnames)
      implicit none    
      character(len=*):: infile 
      integer :: numfields     
      character(len=3), dimension(numfields):: fieldsnames
      
      integer :: first,last,k,nfld,ios
      logical :: ex
      character(len=3) :: char3


      inquire(exist=ex,file=infile)
      if (.not. ex) then
        print *,'Could not find '//infile
      end if

      open(103,status='old',form='formatted',file=infile)
      ios=0
      nfld=0
      do while (ios==0)
        read(103,100,iostat=ios) char3
        if (ios==0) then
          fieldsnames(nfld+1)=char3
          nfld=nfld+1
        end if
      end do
      close(103)
  100 format (a3)

      if (nfld/=numfields) then
        print *,'An error occured when reading '//infile
      end if

      ! List fields used in analysis
      do k=1,numfields
        print *,fieldsnames(k)
      end do

      end subroutine      

end module m_ensemble_stats 
