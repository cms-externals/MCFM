!--------------------------------------------------------------- 
! Subroutine for pp->gamma + jet integrated frag dipoles-------------------
!--------------------------------------------------------------- 
!--- Author C. Williams Feb 2011 
!-----------------------------------------------------------------



      subroutine qqb_dirgam_fragdips(p,qcd_tree,qcd_tree_s,msq_out) 
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'frag.f'
      double precision p(mxpart,4)
      double precision msq_qcd(-nf:nf,-nf:nf),msq_out(-nf:nf,-nf:nf)
      double precision msq_qcd_s(-nf:nf,-nf:nf)
      integer j,k
      double precision virt_dips(3),xl(3),dot,fsq 
      double precision aewo2pi,fi_gaq,ff_gaq
      external qcd_tree
      external qcd_tree_s

      aewo2pi=esq/(fourpi*twopi)      
      
      fsq=frag_scale**2
      

     
!---- Integrated dipoles are functions of p_gamma = z * pjet so need to rescale pjet

     
      call rescale_pjet(p) 
      
     

      do j=1,2
         xl(j)=dlog(-two*dot(p,j,3)/fsq)
      enddo

      xl(3)=dlog(two*dot(p,3,4)/fsq)
     

      do j=1,2
         virt_dips(j)=+aewo2pi*(fi_gaq(z_frag,p,xl(j),3,j,2))
      enddo

      
      virt_dips(3)=+aewo2pi*(ff_gaq(z_frag,xl(3),2))


!---- Matrix elements conserve momenta thro pjet = sum of rest so return orignal pjet

      call return_pjet(p)
     
     



      do j=-nf,nf
         do k=-nf,nf
            msq_qcd(j,k)=0d0
            msq_out(j,k)=0d0
            msq_qcd_s(j,k)=0d0
         enddo
      enddo

      
      
      call qcd_tree(p,msq_qcd) 
      call qcd_tree_s(p,msq_qcd_s)
    
      do j=-nf,nf
         do k=-nf,nf
            
           

            if((j.eq.0).and.(k.gt.0)) then
               msq_out(j,k)=msq_qcd(j,k)*Q(k)**2*virt_dips(2)
            elseif((j.eq.0).and.(k.lt.0)) then 
               msq_out(j,k)=msq_qcd(j,k)*Q(abs(k))**2*virt_dips(2)
            elseif((j.gt.0).and.(k.eq.0)) then
               msq_out(j,k)=msq_qcd(j,k)*Q(j)**2*virt_dips(1)              
            elseif((j.lt.0).and.(k.eq.0)) then          
               msq_out(j,k)=msq_qcd(j,k)*Q(abs(j))**2*virt_dips(1)
            elseif((j.eq.0).and.(k.eq.0)) then 
!---- gg->qq*int_dips 
               msq_out(j,k)=two*msq_qcd(j,k)*(virt_dips(3))*
     &              (dfloat(nf-2)*Q(1)**2+2d0*Q(2)**2)   
         
            elseif((j.ne.-k).and.(j.ne.0).and.(k.ne.0)) then 
               msq_out(j,k)=msq_qcd(j,k)*(virt_dips(1)*Q(abs(j))**2+
     &          virt_dips(2)*Q(abs(k))**2)
            elseif((j.eq.-k).and.(j.ne.0)) then 
!---- t-channel ones
               msq_out(j,k)=msq_qcd(j,k)*(virt_dips(1)*Q(abs(j))**2
     &          +virt_dips(2)*Q(abs(k))**2)
               msq_out(j,k)=msq_out(j,k)+msq_qcd_s(j,k)
     &                 *two*(virt_dips(3))*
     &              (dfloat(nf-2)*Q(1)**2+2d0*Q(2)**2)   
            endif
            
         enddo
      enddo
     
     
      return 
      end subroutine


        
