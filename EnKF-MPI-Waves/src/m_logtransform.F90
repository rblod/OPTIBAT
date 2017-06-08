MODULE m_logtransform
    IMPLICIT NONE
!


! !PUBLIC MEMBER FUNCTIONS:
!
! !PUBLIC DATA MEMBERS:
!
    type anamorphosis
!        character(len=3) smooth !nature of smoothness of the empirical anamorphosis!
        character(len=3) id !nature of the physical variable!
!        integer::n !size of the sample!
!        real, dimension(:), allocatable::y  !gaussian values invG(i/N)!
!        real, dimension(:), allocatable::z  !sorted physical values!
!        real, dimension(:), allocatable::phi !anamorphosis coefficients!
        real::Gmin,Gmax !bounds in gaussian space!
        real::Zmin,Zmax !bounds in real space!
!        real:: alpha !linear interpolation: translation coefficient!
!        integer:: m !0:linear tail, 1:nonlinear tail! 
!        integer::nsamp,samp !size of the sampling after clustering!
!        real::eps !threshold (lowest relevant value)!
!	integer::ifirst
!	real::a1,a2 !param of the left tail: nonlinear function!
!	logical::obs
    end type anamorphosis
    
    character(len=*), parameter, public :: infile_ana='analysisfields_ana.in'
    integer,save :: numfields_ana
    character(len=3), dimension(:), save, allocatable:: fieldnames_ana
    type(anamorphosis),dimension(:),save,allocatable::ana_enkf
    
CONTAINS    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer function ana_get_ind(ana,nana,char3) 
      implicit none
      integer :: nana
      type(anamorphosis),dimension(nana)::ana
      character(len=3) :: char3
      integer::k

      do k=1,nana
        if(trim(ana(k)%id)==trim(char3))then
          ana_get_ind=k
          return
       endif  
      enddo
      ana_get_ind=-1
   
      end function
   
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!!!!         reading of analysisfields_ana.in     !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer function get_nrfields_ana()
      implicit none
      integer :: ios,first,last
      logical :: ex
      character(len=3) :: char3

      inquire(exist=ex,file=infile_ana)
      if (.not. ex) then
        print *,'Could not find '//infile_ana
      end if

      open(103,status='old',form='formatted',file=infile_ana)
      ios=0
      get_nrfields_ana=0
      do while (ios==0)
        read(103,100,iostat=ios) char3
        if (ios==0) get_nrfields_ana=get_nrfields_ana+1
      end do
      close(103)
  100 format (a3,2i3)
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
   
      subroutine get_analysisfields_ana()
      implicit none
      integer :: first,last,k,nfld,ios
      logical :: ex
      character(len=3) :: char3

      numfields_ana=get_nrfields_ana()
      print *,'numfields_ana is ',numfields_ana
      if (numfields_ana<=0 .or.numfields_ana > 20) then !
         print *,'numfields_ana is higher than max allowed setting or = 0'
      end if
      allocate(fieldnames_ana(numfields_ana))

      inquire(exist=ex,file=infile_ana)
      if (.not. ex) then
        print *,'Could not find '//infile_ana
      end if

      open(103,status='old',form='formatted',file=infile_ana)
      ios=0
      nfld=0
      do while (ios==0)
        read(103,100,iostat=ios) char3
        if (ios==0) then
          fieldnames_ana (nfld+1)=char3
          nfld=nfld+1
        end if
      end do
      close(103)
  100 format (a3,2i3)

      if (nfld/=numfields_ana) then
        print *,'An error occured when reading '//infile_ana
      end if

      ! List fields used in analysis
      do k=1,numfields_ana
        print *,fieldnames_ana(k)
      end do

      end subroutine   
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
         
   subroutine ana_log_param(ana)
  ! reading the input files : one per variable to transform
   implicit none
   type(anamorphosis)::ana
   character(len=80) ::memfile
   logical::ex
   
   memfile='infile_ana.'//trim(ana%id)//'.in'
   
   inquire(exist=ex,file=memfile)
   if (.not. ex) then
      print *,'Could not find '//memfile
   end if
   
   open(105,status='old',form='formatted',file=memfile)
   read(105,*) ana%Zmin   
   read(105,*) ana%Zmax
   close(105)
   
   ana%Gmin=log(ana%Zmin)
   ana%Gmax=log(ana%Zmax)
   
   end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! Building the structure

    subroutine init_anamorphosis()
    implicit none
!    type(anamorphosis), dimension(:),allocatable::ana  
    integer::k,m

    
    !reading of the data which will be anamorphosed: analysisfields_ana.in!
    call get_analysisfields_ana()
   
    !allocation of the anamorphosis structure!
    allocate(ana_enkf(numfields_ana))
       
    do k=1,numfields_ana
      ana_enkf(k)%id=trim(fieldnames_ana(k))
      print *,'----------   ',ana_enkf(k)%id,'   ----------'
      
      !reading of the paramaters of the transformation!
      call ana_log_param(ana_enkf(k))
    enddo
    
    end   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subrotuines called to transform variables

  subroutine ana_log(y,z,ana)
  use m_sort2
  implicit none
  real :: y,z
  type(anamorphosis)::ana
  integer::ind,i
  real::b,tmp
  real,dimension(:),allocatable::yt
  real::d
  
    if(y.le.log(ana%Zmin))then
      z=ana%Zmin
      !print*,'ana_log',trim(ana%id),z 
      return
    elseif(y.ge.log(ana%Zmax))then 
      z=ana%Zmax
      !print*,'ana_log',trim(ana%id),z
      return
    endif  
    
    z=exp(y)

  end subroutine  
   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine ana_log_inv(z,y,ana)
  use m_sort2
  implicit none
  real :: y,z
  type(anamorphosis)::ana
  integer::ind,i
  real::b,tmp
  real,dimension(:),allocatable::zt
  real::d


    if(z.le.ana%Zmin)then
      y=log(ana%Zmin)    
      return
    elseif(z.ge.ana%Zmax)then 
      y=log(ana%Zmax)     
      return
    endif

    y=log(z)
  
  end subroutine
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine ana_log_inv_2D(fld,nx,ny,depth,ana)   
   implicit none
   integer :: nx,ny   
   type(anamorphosis)::ana
   real, dimension(nx,ny) :: fld,depth 
   integer :: i,j
   real:: d
  
    do j=1,ny
      do i=1,nx
        if(depth(i,j).gt.0.)then
          call ana_log_inv(fld(i,j),d,ana)    
          fld(i,j)=d
	endif
      enddo
    enddo
        
   end subroutine
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
   subroutine ana_log_2D(fld,nx,ny,depth,ana)   
   implicit none
   integer :: nx,ny   
   type(anamorphosis)::ana
   real, dimension(nx,ny) :: fld,depth
   integer :: i,j
   real ::d
   
    do j=1,ny
      do i=1,nx
        if(depth(i,j).gt.0.)then
          call ana_log(fld(i,j),d,ana)       
          fld(i,j)=d
	endif  
      enddo
    enddo
       
   end subroutine 
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    subroutine ana_chi_exp_S(S,nrobs,nrens,iens,ana)   
   implicit none
   integer :: nrobs,nrens,iens        ! Number of measurements
   real, dimension(1:nrobs,1:nrens) :: S
   type(anamorphosis)::ana
   integer::k
   real:: d
     
    do k=1,nrobs
      call ana_log_inv(S(k,iens),d,ana) 
      S(k,iens)=d
    enddo 
       
   end subroutine  
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine ana_var_obs(obs,ana,nrens)
      !compute the observation error variance in log-space and transform osbervations
      use mod_measurement
      use m_random
      implicit none
      integer::nrens
      type(measurement):: obs
      type(anamorphosis)::ana
      real::varo!,meandx
      integer::k,iens
      real::md,d
      real,dimension(nrens)::work,work2
                  
      call randn(nrens,work2)
      
      do iens=1,nrens
        if(obs%d==0.)then
          print*,'obs null',iens
          obs%d=ana%Zmin
        endif
        work(iens)=obs%d*exp(sqrt(obs%var)*&
	                 (work2(iens)-sqrt(obs%var)/2.))
	
        call ana_log_inv(work(iens),d,ana)
        work(iens)=d	
      enddo	
      
      md=sum(work)/real(nrens)
      varo=0.
      do iens=1,nrens
        varo=varo+(work(iens)-md)**2
      enddo	
      obs%var=varo/real(nrens-1)
    
      
      !print*,'avant',obs%d
      call ana_log_inv(obs%d,d,ana)
      obs%d=d
      !print*,'apres',obs%d

      end 
      
       
end module m_logtransform
