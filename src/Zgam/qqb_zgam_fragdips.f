!--------------------------------------------------------------- 
!   This subroutine checks the number of external dipoles----
!---absorbing the correct number into the fragmenation functions 
!   it then returns the finite (msq_qcd*dip) ---------------------
!--------------------------------------------------------------- 
!--- Author C. Williams Dec 2010 
!-----------------------------------------------------------------


      subroutine qqb_zgam_fragdips(p,qcd_tree,msq_out) 
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

      call rescale_pjet(p)

      do j=1,2
         xl(j)=dlog(-two*dot(p,j,5)/fsq)
      enddo
      
      do j=1,2
         virt_dips(j)=+aewo2pi*(fi_gaq(z_frag,p,xl(j),5,j,2))
      enddo
      
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
            
           if((j.eq.0).and.(k.ne.0)) then 
              msq_out(j,k)=Q(k)**2*msq_qcd(j,k)*virt_dips(2)
           elseif((j.ne.0).and.(k.eq.0)) then 
              msq_out(j,k)=Q(j)**2*msq_qcd(j,k)*virt_dips(1)
           endif
            
         enddo
      enddo
     
     
      return 
      end subroutine


        
