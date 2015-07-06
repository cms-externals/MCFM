!--------------------------------------------------------------- 
!   This subroutine checks the number of external dipoles----
!---absorbing the correct number into the fragmenation functions 
!   it then returns the finite (msq_qcd*dip) ---------------------
!--------------------------------------------------------------- 
!--- Author C. Williams April 2012 
!-----------------------------------------------------------------



      subroutine qqb_zaa_fragdips(p,qcd_tree,msq_out) 
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'frag.f'
      double precision p(mxpart,4)
      double precision msq_qcd(-nf:nf,-nf:nf),msq_out(-nf:nf,-nf:nf)
      integer j,k
      double precision virt_dips(2),xl(2),dot,fsq 
      double precision aewo2pi,fi_gaq
      external qcd_tree
            

      aewo2pi=esq/(fourpi*twopi)      
      
      fsq=frag_scale**2
     
      
!---- Integrated dipoles are functions of p_gamma = z * pjet so need to rescale pjet

     
      call rescale_pjet(p) 
      
     

      do j=1,2
         xl(j)=dlog(-two*dot(p,j,6)/fsq)
      enddo
      
      do j=1,2
         virt_dips(j)=+aewo2pi*(fi_gaq(z_frag,p,xl(j),6,j,2))
      enddo
      

!---- Matrix elements conserve momenta thro pjet = sum of rest so return orignal pjet

      call return_pjet(p)
     
     



      do j=-nf,nf
         do k=-nf,nf
            msq_qcd(j,k)=0d0
            msq_out(j,k)=0d0
         enddo
      enddo

  
      
      call qcd_tree(p,msq_qcd) 


      do j=-nf,nf
         do k=-nf,nf
            
!   factor of two cancelled by statistical factor because two photons
            if((j.eq.0).and.(k.gt.0)) then
                  msq_out(j,k)=msq_qcd(j,k)*Q(k)**2*virt_dips(2)
            elseif((j.eq.0).and.(k.lt.0)) then 
                 msq_out(j,k)=msq_qcd(j,k)*Q(abs(k))**2*virt_dips(2)
            elseif((j.gt.0).and.(k.eq.0)) then
                  msq_out(j,k)=msq_qcd(j,k)*Q(j)**2*virt_dips(1)              
            elseif((j.lt.0).and.(k.eq.0)) then          
                 msq_out(j,k)=msq_qcd(j,k)*Q(abs(j))**2*virt_dips(1)
            
            endif
            
         enddo
      enddo
     

      return 
      end subroutine


        
