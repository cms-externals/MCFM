
! Modified FeB 11 rescale_pjet is now a driver routine for various nprocs
! Still takes pjet as an argument

c---- Subroutine for rescaling momentum by z_frag
      subroutine rescale_pjet(pjet)
      include 'constants.f'
      include 'process.f'
      double precision pjet(mxpart,4)
        
     

      if((case.eq.'Wgamma').or.(case.eq.'Zgamma')) then 
         call rescale_1phot(pjet,5)
      elseif(case.eq.'gamgam') then
         call rescale_1phot(pjet,4)
      elseif(case.eq.'dirgam') then 
         call rescale_1phot(pjet,3) 
      elseif(case.eq.'Z_2gam') then 
         call rescale_1phot(pjet,6) 
      elseif(case.eq.'Zgajet') then 
         call rescale_1phot(pjet,5) 
      else

         write(6,*) 'Error: tried to rescale unknown quantity' 
         stop
      endif     

!--- Old routine rescales all photons by z_frag -----     
!      do i=1,mxpart  
!         if(plabel(i) .eq. 'ga')  then
!            do j=1,4
!               pjet(i,j) = z_frag*pjet(i,j)
!            enddo
!         endif
!      enddo
      
      end 
      
    
      

c---- Subroutine for returning original momenta 

      subroutine return_pjet(pjet)
      include 'constants.f'
      include 'process.f'
      double precision pjet(mxpart,4)
      

      
      if((case.eq.'Wgamma').or.(case.eq.'Zgamma')) then 
         call return_1phot(pjet,5)
      elseif(case.eq.'gamgam') then
         call return_1phot(pjet,4)
      elseif(case.eq.'dirgam') then 
         call return_1phot(pjet,3)
      elseif(case.eq.'Z_2gam') then 
         call return_1phot(pjet,6) 
      elseif(case.eq.'Zgajet') then 
         call return_1phot(pjet,5) 
      else
         write(6,*) 'Error: tried to rescale unknown quantity' 
         stop
      endif     


!     Old routine
!      do i=1,mxpart  
!         if(plabel(i) .eq. 'ga')  then
!            do j=1,4
!               pjet(i,j) = (1d0/z_frag)*pjet(i,j)
!            enddo
!         endif
!      enddo
      
      end 

      
!---- New routines 

      subroutine rescale_1phot(p,i) 
      implicit none
      include 'constants.f'
      include 'frag.f'
      integer i,j 
      double precision p(mxpart,4)
      
      do j=1,4
         p(i,j)=z_frag*p(i,j)
      enddo
      
      return 
      end 

      subroutine return_1phot(p,i) 
      implicit none
      include 'constants.f'
      include 'frag.f'
      integer i,j 
      double precision p(mxpart,4)
      
      do j=1,4
         p(i,j)=(1d0/z_frag)*p(i,j)
      enddo
      
      return 
      end 

      
     
