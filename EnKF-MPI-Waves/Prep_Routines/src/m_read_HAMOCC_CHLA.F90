module m_read_HAMOCC_CHLA
! Ehouarn !
! Reads CHLA1 data and grid from hamocc output files !
contains


  subroutine get_CHLA_mask(fname,nx,ny,nwnd,depths,mask,cpt)
  use nfw_mod
  implicit none
  character(80) :: fname
  integer::nx,ny,nwnd,cpt
  integer, dimension(nx,ny,nwnd)::mask
  real::att(4)
  character(80) :: memfile
  integer::i,j,k,vFIELD_ID,ncid
  character(3)::cmonth
  real,dimension(nx,ny)::fld,fld2,depths
  
  mask(:,:,:)=0
  cpt=0
  do k=1,nwnd
      write(cmonth,'(i3.3)')k
      memfile=trim(fname)//cmonth
      
      call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'swa',vFIELD_ID)
      call nfw_get_var_double(trim(memfile)//'.nc', ncid, vFIELD_ID, fld)
      call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'add_offset', att(2))
      call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'scale_factor', att(1))
    
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'fice',vFIELD_ID)
      call nfw_get_var_double(trim(memfile)//'.nc', ncid, vFIELD_ID, fld2)
      call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'add_offset', att(4))
      call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'scale_factor', att(3))
     
      call nfw_close(trim(memfile)//'.nc', ncid)
       
      do j=1,ny
        do i=1,nx	  
	  if((depths(i,j).gt.0.))then
	  if(((fld(i,j)*att(1)+att(2)).gt. 40.).and.((fld2(i,j)*att(3)+att(4))==0.))then
	    mask(i,j,k)=1
	    cpt=cpt+1	  
	  endif 
	  endif	  
	enddo
      enddo	 
  enddo
  
  end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   subroutine read_hamocc_ACHLA(fname,obs,mask,nx,ny,nwnd,nobs,nrobs,modlon,modlat)
   use mod_measurement
   use m_set_random_seed2
   use m_random
   use nfw_mod
   implicit none
   integer::nx,ny,nwnd,nobs,nrobs
   type(measurement):: obs(nobs)
   character(80) :: fname,memfile
   integer::d0,vFIELD_ID,ncid
   character(3)::cmonth
   real,dimension(nx,ny)::fld,modlon,modlat
   integer,dimension(nx,ny,nwnd)::mask
   real::att(2)
   integer::i,j,k,l,cpt
   real,dimension(nx*ny)::noise
   integer,dimension(4)::nc,ns
   real::varb,c2chl
   
   
   call set_random_seed3
  
   ns(1:4)=1
   nc(1)=nx
   nc(2)=ny
   nc(3:4)=1
   varb=0.09
   c2chl=60.
   
   cpt=1
   do l=1,nwnd
      k=1
      call randn(nx*ny,noise)
      
      write(cmonth,'(i3.3)')l
      memfile=trim(fname)//cmonth
      
      call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'phyc',vFIELD_ID)
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vFIELD_ID, ns, nc, fld)
      call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'add_offset', att(2))
      call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'scale_factor', att(1))
      call nfw_close(trim(memfile)//'.nc', ncid)
      
      do j=1,ny
        do i=1,nx
          if(mask(i,j,l)==1)then
	    obs(cpt)%id='ACHLA'
	    obs(cpt)%d=(fld(i,j)*att(1)+att(2))*exp(sqrt(varb)*noise(k)-0.5*varb)*1.e6/c2chl
	    obs(cpt)%var=0.09
	    obs(cpt)%lat=modlat(i,j)
	    obs(cpt)%lon=modlon(i,j)
	    obs(cpt)%ipiv=i
	    obs(cpt)%jpiv=j
	    obs(cpt)%status=.true.
	    !if(l.gt.1)obs(cpt)%status=.false.	    
	    obs(cpt)%depth=0.
	    obs(cpt)%h=1
	    obs(cpt)%ns=1
	    obs(cpt)%a1=1.
	    obs(cpt)%i_orig_grid = -1
            obs(cpt)%j_orig_grid = -1
	    obs(cpt)%date=l
	    
	    obs(cpt)%a1 = 1
            obs(cpt)%a2 = 0
            obs(cpt)%a3 = 0
            obs(cpt)%a4 = 0
     
	    cpt=cpt+1
	  endif  
	  k=k+1  
	enddo
      enddo	
   enddo
   nrobs=nobs
    
   end subroutine 
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    
  subroutine get_SST_mask(memfile,nx,ny,depths,mask,cpt)
  use nfw_mod
  implicit none
  character(80) :: fname
  integer::nx,ny,cpt
  integer, dimension(nx,ny)::mask
  real::att(2)
  character(80) :: memfile
  integer::i,j,k,vFIELD_ID,ncid
  character(3)::cmonth
  real,dimension(nx,ny)::fld,depths
  
  mask(:,:)=0
  cpt=0

  call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
  call nfw_inq_varid(trim(memfile)//'.nc', ncid,'fice',vFIELD_ID)
  call nfw_get_var_double(trim(memfile)//'.nc', ncid, vFIELD_ID, fld)
  call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'add_offset', att(2))
  call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'scale_factor', att(1))
  call nfw_close(trim(memfile)//'.nc', ncid)
       
  do j=1,ny
    do i=1,nx	  
      if((depths(i,j).gt.0.).and.((fld(i,j)*att(1)+att(2))==0.))then
	 mask(i,j)=1
	 cpt=cpt+1	  
      endif 
    enddo
  enddo	 
  
  end subroutine
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine read_hamocc_SST(memfile,obs,mask,nx,ny,nobs,nrobs,modlon,modlat)
   use mod_measurement
   use m_set_random_seed2
   use m_random
   use nfw_mod
   implicit none
   integer::nx,ny,nobs,nrobs
   type(measurement):: obs(nobs)
   character(80) :: fname,memfile
   integer::d0,vFIELD_ID,ncid
   character(3)::cmonth
   real,dimension(nx,ny)::fld,modlon,modlat
   integer,dimension(nx,ny)::mask
   real::att(2)
   integer::i,j,k,l,cpt
   real,dimension(nx*ny)::noise
   real::varb
   
   
   call set_random_seed3
  
   varb=0.01
   
   call randn(nx*ny,noise)    
      
   call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
   call nfw_inq_varid(trim(memfile)//'.nc', ncid,'sst',vFIELD_ID)    
   call nfw_get_var_double(trim(memfile)//'.nc', ncid, vFIELD_ID, fld)
   call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'add_offset', att(2))
   call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'scale_factor', att(1))     
   call nfw_close(trim(memfile)//'.nc', ncid)
  
   cpt=1 
   do j=1,ny
      do i=1,nx
        if(mask(i,j)==1)then
	   obs(cpt)%id='SST'
	   obs(cpt)%d=max((fld(i,j)*att(1)+att(2))+sqrt(varb)*noise(cpt),-1.81618)
	   obs(cpt)%var=0.01
	   obs(cpt)%lat=modlat(i,j)
	   obs(cpt)%lon=modlon(i,j)
	   obs(cpt)%ipiv=i
	   obs(cpt)%jpiv=j
	   obs(cpt)%status=.true.
	   obs(cpt)%depth=0.
	   obs(cpt)%h=1
	   obs(cpt)%ns=1
	   obs(cpt)%a1=1.
	   obs(cpt)%i_orig_grid = -1
           obs(cpt)%j_orig_grid = -1
	   obs(cpt)%date=0
	   cpt=cpt+1
	endif  
	  
       enddo
    enddo	

   nrobs=nobs
    
   end subroutine 
   
  

end module m_read_HAMOCC_CHLA
