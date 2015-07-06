************************************************************************ 
*     This subroutine calculates dipoles with an  identified           *
*     Final state photon                                               *
*     C. Williams Dec 2010                                             * 
*     Returns the dipoles in sub,subv and matrix elements in msq,msqv  *
*     nd labels the dipole configurations                              *
*     ip labels the emitter PHOTON                                     *
*     jp labels the emitted parton                                     *
*     kp labels the spectator parton                                   *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *
************************************************************************

c--- Extension to dipolesfrag that also returns msqx,
c--- 4 dimensional array indexed by initial and final parton labels

      subroutine dipsfragx(nd,p,ip,jp,kp,sub,msq,msqx,subr_born) 
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'ptilde.f'
      include 'dynamicscale.f'
      include 'initialscales.f' 
      include 'dipolescale.f'
      include 'facscale.f'
      include 'betacut.f'
      double precision p(mxpart,4),ptrans(mxpart,4),sub
      double precision z,omz,sij,sik,sjk,dot,u
      double precision msq(-nf:nf,-nf:nf)
      double precision mdum1(0:2,fn:nf,fn:nf) ! mqq
      double precision msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision mdum2(0:2,-nf:nf,-nf:nf) !msqx_cs
      integer nd,ip,jp,kp,j,k
      logical incldip(0:maxd),check_nv,phot_pass
      common/incldip/incldip
      external subr_born

      z=0d0
      omz=1d0
      u=0d0
      sub=0d0
      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=0d0
         enddo
      enddo
      
      incldip(nd)=.true.

      sij=two*dot(p,ip,jp)
      sik=two*dot(p,ip,kp) 
      sjk=two*dot(p,jp,kp)

     
******************************************************************************* 
************************ INITIAL - INITIAL ************************************
*******************************************************************************

**** I === I not implemented (need photon PDFS rather than frags) need rapidity cuts to remove collinear sing 
      if ((ip .le. 2) .and. (kp .le. 2)) return  

******************************************************************************* 
************************ INITIAL - FINAL **************************************
*******************************************************************************

**** I === F not implemented (need photon PDFS rather than frags) need rapidity cuts to remove collinear sing
      if ((ip .le. 2) .and. (kp .gt. 2)) return  

******************************************************************************* 
************************ FINAL - INITIAL  *************************************
*******************************************************************************

      
      if ((ip .gt. 2) .and. (kp .le. 2)) then
         
     
         if(check_nv(p,ip,jp,kp).eqv..false.) then 
            incldip(nd)=.false.
            return 
         endif
         call transformfrag(p,ptrans,z,ip,jp,kp)
         call storeptilde(nd,ptrans) 
         call store_zdip(nd,z)
         omz=one-z
c----  Check that photon will still pass cuts 
	 if (phot_pass(ptrans,ip,z) .eqv. .false.) return 
c--- if using a dynamic scale, set that scale with dipole kinematics	
	if (dynamicscale) then
           call rescale_z_dip(ptrans,nd,ip)
           call scaleset(initscale,initfacscale,ptrans)
           call return_z_dip(ptrans,nd,ip)
	   dipscale(nd)=facscale
	endif

         call subr_born(ptrans,msq,mdum1,msqx,mdum2)
         
         sub=two*(esq/sij)*((one+omz**2)/z)
        
         
******************************************************************************* 
************************ FINAL - FINAL ****************************************
*******************************************************************************
         
      elseif ((ip .gt. 2) .and. (kp .gt. 2)) then 

         z=(sij+sik)/(sik+sjk+sij) 
         omz=one-z
         
         u=sij/(sij+sik) 
       
         if(u.gt.bff) then
            incldip(nd)=.false.
            return 
         endif
            
c--- Calculate the ptrans-momenta 
         call transformfrag(p,ptrans,z,ip,jp,kp)
         call storeptilde(nd,ptrans) 
         call store_zdip(nd,z)

c----  Check that photon will still pass cuts 
	 if (phot_pass(ptrans,ip,z) .eqv. .false.) return 
c--- if using a dynamic scale, set that scale with dipole kinematics	
	if (dynamicscale) then
           call rescale_z_dip(ptrans,nd,ip)
           call scaleset(initscale,initfacscale,ptrans)
           call return_z_dip(ptrans,nd,ip)
	   dipscale(nd)=facscale
	endif


         call subr_born(ptrans,msq,mdum1,msqx,mdum2)
      
         sub=two*(esq/sij)*((one+omz**2)/z)
      
      endif

c      write(6,*) 'dipole ',nd,' momenta'
c      call writeout(ptrans)

      return 
      end

