      subroutine qqb_gamgam_gs(p,msq)
************************************************************************
*     Author: J.M. Campbell                                            *
*     October, 2002.                                                   *
*     Modified by CW to Include photon fragmentation dipoles Feb 11    *
*    Matrix element SUBTRACTION squared averag'd over init'l colors    *
*    and spins                                                         *
*     f(-p1) + f(-p2) -->  gamma(p3) + gamma(p4) + parton(p5)          *
************************************************************************

      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'frag.f'
      include 'ewcharge.f'
      integer j,k,nd
      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),sub15_2v,sub25_1v,
     & msq15_2v(-nf:nf,-nf:nf),msq25_1v(-nf:nf,-nf:nf),
     & msq35_1(-nf:nf,-nf:nf),msq35_2(-nf:nf,-nf:nf),
     & msq45_1(-nf:nf,-nf:nf),msq45_2(-nf:nf,-nf:nf), 
     & sub35_1,sub35_2,sub45_1,sub45_2
      
      logical phot_dip(mxpart)
      common/phot_dip/phot_dip 

      external qqb_gamgam,qqb_gamgam_gvec
      external qqb_dirgam,qqb_dirgam_swap 
      


      if(frag) then 
         ndmax=6
      else
         ndmax=2
      endif

!---- Intialise Photon dipoles 
      do j=1,mxpart
        phot_dip(j)=.false.
      enddo

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_gamgam,qqb_gamgam_gvec)
      call dips(2,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & qqb_gamgam,qqb_gamgam_gvec)

!---- Extra photon dipoles 
      if (frag) then 
!---- gamma(3) dipoles
         call dipsfrag(3,p,3,5,1,sub35_1,msq35_1,qqb_dirgam_swap)
         call dipsfrag(4,p,3,5,2,sub35_2,msq35_2,qqb_dirgam_swap)         
!----- gamma(4) dipoles 
         call dipsfrag(5,p,4,5,1,sub45_1,msq45_1,qqb_dirgam) 
         call dipsfrag(6,p,4,5,2,sub45_2,msq45_2,qqb_dirgam)
         do j=3,6 
           phot_dip(j)=.true. 
         enddo
      endif


      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0d0
      enddo
      
      if(frag) then 
!      g-q frag dipoles 
         if((j.eq.0) .and. (k.gt.0)) then            
            msq(4,j,k)=Q(k)**2*msq35_2(j,k)*sub35_2*half
            msq(6,j,k)=Q(k)**2*msq45_2(j,k)*sub45_2*half
         elseif((j.eq.0).and.(k.lt.0)) then 
!     g-qbar frag dipoles
            msq(4,j,k)=Q(abs(k))**2*msq35_2(j,k)*sub35_2*half 
            msq(6,j,k)=Q(abs(k))**2*msq45_2(j,k)*sub45_2*half
          elseif((j.gt.0).and.(k.eq.0)) then 
!     q-g frag dipoles 
             msq(3,j,k)=Q(j)**2*msq35_1(j,k)*sub35_1*half
             msq(5,j,k)=Q(j)**2*msq45_1(j,k)*sub45_1*half
          elseif((j.lt.0).and.(k.eq.0)) then 
!     qbar-g frag dipoles
             msq(3,j,k)=Q(abs(j))**2*msq35_1(j,k)*sub35_1*half
             msq(5,j,k)=Q(abs(j))**2*msq45_1(j,k)*sub45_1*half
          endif
       endif
          
       if ((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) goto 19    

c--- do only q-qb and qb-q cases      
      if (  ((j .gt. 0).and.(k .lt. 0))
     . .or. ((j .lt. 0).and.(k .gt. 0))) then
C-----half=statistical factor ????
         msq(1,j,k)=2d0*cf*sub15_2(qq)*msq15_2(j,k)
         msq(2,j,k)=2d0*cf*sub25_1(qq)*msq25_1(j,k)
      elseif ((j .ne. 0) .and. (k .eq. 0)) then
         msq(2,j,k)=2d0*tr*sub25_1(qg)*msq25_1(j,-j)
      elseif ((j .eq. 0) .and. (k .ne. 0)) then
         msq(1,j,k)=2d0*tr*sub15_2(qg)*msq15_2(-k,k)
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
c--- Comment out the following four lines to remove gg contribution
         msq(1,j,k)=2d0*xn*(sub15_2(gg)*msq15_2(j,k)
     &                     +sub15_2v*msq15_2v(j,k))
         msq(2,j,k)=2d0*xn*(sub25_1(gg)*msq25_1(j,k)
     &                     +sub25_1v*msq25_1v(j,k))
      endif

 19   continue
      enddo
      enddo

      return
      end
      
