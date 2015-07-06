      subroutine gen_photons_jets(r,nphots,njets,p,wt,*)
c----  WARNING: although mostly written in a generic manner, this
c----           routine has been tailored for direct photon (dirgam)
c----           and photon+two jets (gamjet) production

c---- generate phase space for 2-->nphots+njets process
c----   and 3,..,3+nphots the photons
c----   and 4+nphots,..,4+nphots+njets the jets
c----
c---- This routine is adapted from gen_phots_jets.f
c----
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(4+2n))
c----
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'leptcuts.f'
      include 'reset.f'
      include 'xmin.f'
      include 'part.f'
      double precision r(mxdim),p(mxpart,4),psumjet(4),wt
      double precision hmin,hmax,delh,h,sqrts,pt,etamax,xx(2)
      double precision y,sinhy,coshy,phi
      double precision ptjetmin,etajetmin,etajetmax,pbreak
      double precision ptmin_part,etamax_part,pbreak_part,swap
      integer j,nu,nphots,njets,nphotsjets,ijet,icount
      logical first,flatreal
      common/energy/sqrts
      common/x1x2/xx
      parameter(flatreal=.false.)
      integer perm 
      data perm /1/ 
      data first/.true./
      save first,ptjetmin,etajetmin,etajetmax,pbreak
      
      if (first .or. reset) then
        first=.false.
        reset=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
        if (part .eq. 'real') then
c--- if we're generating phase space for real emissions, then we need
c--- to produce partons spanning the whole phase space pt>0,eta<10;
c--- in this case, pbreak=ptjetmin simply means that we
c--- generate pt approx. 1/x for pt > pbreak and
c--- pt approx. uniformly for pt < pbreak
          pbreak=ptjetmin
          ptjetmin=0d0
          etajetmax=20d0
        else
c--- for lord and virt, the partons produced here can be generated
c--- right up to the jet cut boundaries and there is no need for pbreak
          pbreak=0d0
        endif
c--- in case this routine is used for very small values of ptjetmin
        if ((ptjetmin .lt. 5d0) .and. (part .ne. 'real')) pbreak=5d0
c--- for processes in which it is safe to jet ptmin to zero at NLO
	if ((part .eq. 'real') .and. (pbreak .lt. 1d-8)) pbreak=5d0
      endif        


c--- total number of photons and jets
      nphotsjets=nphots+njets

      do nu=1,4
        do j=1,4+nphotsjets
          p(j,nu)=0d0
        enddo
        psumjet(nu)=0d0
      enddo 

      wt=2d0*pi
                       
     
      icount=0
      do ijet=1,nphotsjets-1
c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
        wt=wt/16d0/pi**3
c        xmin=2d0/sqrts
c        xmax=1d0/ptjetmin

        if (ijet .le. nphots) then
          ptmin_part=gammpt
          etamax_part=gammrap
          pbreak_part=0d0
          if (part .eq. 'real') then
c--- cannot generate exactly to match, since dipoles transform photon
	  ptmin_part=0d0
          etamax_part=20d0
          pbreak_part=gammpt
	  endif
        else
	  if (part .eq. 'lord') then
c--- (gamjet) generate jets according to phase-space boundaries
            ptmin_part=ptjetmin
            etamax_part=etajetmax
            pbreak_part=pbreak
	  else
c--- (dirgam) attempt to generate jets to balance photon
            ptmin_part=0d0
            etamax_part=20d0
            pbreak_part=gammpt
	  endif
        endif

        if ((flatreal) .and. (part .eq. 'real')) then
c--- generate flat pt for the real contribution
          pt=r(icount+1)*((sqrts/2d0)-ptmin_part)+ptmin_part
          wt=wt*((sqrts/2d0)-ptmin_part)*pt
        else
c--- favour small pt region 
          hmin=1d0/dsqrt((sqrts/2d0)**2+pbreak_part**2)
          hmax=1d0/dsqrt(ptmin_part**2+pbreak_part**2)
          delh=hmax-hmin
          h=hmin+r(icount+1)*delh        
          pt=dsqrt(1d0/h**2-pbreak_part**2)
          wt=wt*delh/h**3
        endif

        etamax=sqrts/2d0/pt
        if (etamax**2 .le. 1d0) then
          write(6,*) 'etamax**2 .le. 1d0 in gen_phots_jets.f',etamax**2 
          wt=0d0
          return 1
        endif
        etamax=dlog(etamax+dsqrt(etamax**2-1d0))
        
        etamax=min(etamax,etamax_part)
        y=etamax*(2d0*r(icount+2)-1d0)
        wt=wt*2d0*etamax
        
        sinhy=dsinh(y)
        coshy=dsqrt(1d0+sinhy**2)
        
        p(2+ijet,4)=pt*coshy
        
        phi=2d0*pi*r(icount+3)
        wt=wt*2d0*pi
        
        p(2+ijet,1)=pt*dcos(phi)
        p(2+ijet,2)=pt*dsin(phi)
        p(2+ijet,3)=pt*sinhy
       
        do nu=1,4
          psumjet(nu)=psumjet(nu)+p(2+ijet,nu)
        enddo
	
      icount=icount+3
      enddo

!----- have generated nphotsjets-1 massless momenta, last pt is fixed, rapidity is unconstrained       
      etamax=sqrts/2d0/pt
      if (etamax**2 .le. 1d0) then
         write(6,*) 'etamax**2 .le. 1d0 in gen_phots_jets.f',etamax**2 
         wt=0d0
         return 1
      endif
      etamax=dlog(etamax+dsqrt(etamax**2-1d0))
        
     
      etamax=min(etamax,etamax_part)
      
      y=etamax*(2d0*r(3*nphotsjets-2)-1d0)
      wt=wt*2d0*etamax
        
      sinhy=dsinh(y)
      coshy=dsqrt(1d0+sinhy**2)

      pt=dsqrt(psumjet(1)**2+psumjet(2)**2)
      p(2+nphotsjets,4)=pt*coshy
               
      p(2+nphotsjets,1)=-psumjet(1)
      p(2+nphotsjets,2)=-psumjet(2)
      p(2+nphotsjets,3)=pt*sinhy
       
      do nu=1,4
      psumjet(nu)=psumjet(nu)+p(2+nphotsjets,nu)
      enddo

c--- now make the initial state momenta
      xx(1)=(psumjet(4)+psumjet(3))/sqrts
      xx(2)=(psumjet(4)-psumjet(3))/sqrts
      
      if   ((xx(1) .gt. 1d0) .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin).or. (xx(2) .lt. xmin)) then
         wt=0d0
         return 1
      endif 
      
      p(1,4)=-0.5d0*xx(1)*sqrts
      p(1,1)=0d0
      p(1,2)=0d0
      p(1,3)=-0.5d0*xx(1)*sqrts

      p(2,4)=-0.5d0*xx(2)*sqrts
      p(2,1)=0d0
      p(2,2)=0d0
      p(2,3)=+0.5d0*xx(2)*sqrts

      wt=wt/sqrts**2
      
c      call writeout(p) 
c      pause
     
c--- symmetrize phase-space by switching jet momenta every other point    
      if(perm.eq.1) then 
         do nu=1,4
            swap=p(4,nu)
            p(4,nu)=p(5,nu)
            p(5,nu)=swap
         enddo
         perm=2
      elseif(perm.eq.2) then 
         perm=1 
      endif

      return
      end
      
      
      
      
      
      
      
      
      
      
      
