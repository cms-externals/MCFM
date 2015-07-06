      subroutine GGdR_frag(x,i,D,str) 
c--- LO fragmentation functions from:
c---   A.~Gehrmann-De Ridder, E.~W.~N.~Glover,
c---   %``Final state photon production at LEP,''
c---   Eur.\ Phys.\ J.\  {\bf C7}, 29-48 (1999).  [hep-ph/9806316].
      implicit none 
      include 'frag.f'
      include 'qcdcouple.f'
      include 'ewcouple.f' 
      include 'constants.f'
      include 'ewcharge.f'
      double precision x,mu_zero,omx
      double precision pref_ew,pref_ewas
      double precision Pqga,Pqq,P1loop
      double precision aewon2pi,log_mu,ddilog
      double precision lx,dix,lmx,D,Dnp,PxP,PxD
      double precision GGdR_NonP,GGdR_PxP,GGdR_DxP
      integer pow,str,i

c--- Safety cuts::
c---    x < 0.999 to avoid x->1 singularity
c---    x > 0.1 to limit region of extrapolation
      if((x .gt. 0.999d0) .or. (x .lt. 0.1d0)) then 
         D=0d0 
         return 
      endif

c--- set up logarithms 
      omx=one-x
      lx=dlog(x)
      dix=ddilog(1-x) 
      lmx=dlog(omx) 
      pow=2

c--- no gluon fragmentation
      if (i.eq.0) then 
        D=0d0 
        return 
      endif

      aewon2pi=esq/(fourpi*twopi) 
  
c--- prefactors 
      pref_ew=aewon2pi*Q(i)**2
      if     (str.eq.0) then 
c----- LO only, set alpha_s terms= 0 
         pref_ewas = 0d0 
         mu_zero=0.14d0 
      elseif (str.eq.1) then 
         write(6,*) 'NLO not currently implemented.'
         stop 
!         pref_ewas=pref_ew*ason2pi
!         mu_zero=0.64d0 
      else
         write(6,*) 'Unrecognised order in Frag'  
      endif

      log_mu=dlog(frag_scale**2/mu_zero**2)

      
c-- splitting functions      
      Pqga=(one+omx**2)/x
c      P1loop=Cf*(-1d0/2d0+9d0/2d0*x+(-8d0+x/2d0)*lx+2d0*x*lmx
c     &     +(1d0-x/2d0)*lx**2+(lmx**2+4d0*lx*lmx+8d0*dix
c     &     -4d0*pi**2/3d0)*Pqga)

    
!----- Non perturbative pieces 
      Dnp=GGdR_NonP(x,str)

!      Pxp=GGdR_PxP(x) 
!      PxD=GGdR_DxP(x)      

      D=pref_ew*Dnp+log_mu*pref_ew*Pqga

      return 
      end
       
      
      
      double precision function GGdR_NonP(x,str) 
!-----returns the non-pertubrbative pieces str=0,LO 1,NLO 
      implicit none 
      include 'constants.f' 
      double precision x,Pqga
      double precision omx,lmx
      integer str

      omx=one-x
      lmx=dlog(omx**2)
      Pqga=(one+omx**2)/x

      if     (str.eq.0) then 
c-- LO Non Pert eq. 2.13 with prefactor removed 
         GGdR_NonP=-Pqga*lmx-13.26d0 
      elseif (str.eq.1) then
         GGdR_NonP=-Pqga*lmx+20.8d0*omx-11.07d0 
      else
         write(6,*) 'Unrecognised Non P input in GGdR' 
         stop
      endif

      return 
      end



!      double precision function GGdR_PxP(z) 
!      implicit none 
!      include 'constatns.f'
!      include 'grid_ggdrConv.f'
! !---- interpolating P conv P piece of 2.12 
!      double precision z 
!      integer i,jz,it_p
!      parameter(it_p=4)
!      double precision lz(num_g),lpxp(num_g) 
!      double precision Z1(it_p),LP1(it_p)
!      double precision z_min,z_max,a,out_2
!      logical first
!      data first/.true./


!      if(first) then 
!         first=.false.
!         do i=1,num_g
!     lz(i)=dlog10(z_grid(i)) 
!            lpxp(i)=PxPconv(i)
!     enddo
!      endif
      
     
!      z_min=z_grid(1)
!      z_max=z_grid(num_g)

!      GGdR_PxP=0d0
!      if((z.gt.z_min).and.(z.lt.z_max)) then
!------find z in grid 
!         a=dlog10(z) 
!         call locate (lz,num_g,a,jz)
 !        do i=1,it_p
!------Special case when jz= 1
c$$$            if(jz.eq.1) then 
c$$$               Z1(i)=lz(i) 
c$$$               LP1(i)=lpxp(i) 
c$$$            elseif(jz.eq.(num_g-1)) then   
c$$$               Z1(i)=lz(jz-3+i)
c$$$               LP1(i)=lpxp(jz-3+i)
c$$$            elseif(jz.eq.num_g) then 
c$$$               Z1(i)=lz(jz-4+i) 
c$$$               LP1(i)=lpxp(jz-4+i) 
c$$$            else
c$$$               Z1(i)=lz(jz-2+i) 
c$$$               LP1(i)=lpxp(jz-2+i)
c$$$            endif
c$$$         enddo
c$$$      endif
c$$$!------- interpolate 
c$$$           
c$$$      
c$$$      call dpolint(Z1,LP1,4,a,GGdR_PxP,out_2)
c$$$ 
c$$$!-------- return 
c$$$      return 
c$$$      end function 
c$$$
c$$$        double precision function GGdR_DxP(z) 
c$$$      implicit none 
c$$$      include 'constatns.f'
c$$$      include 'grid_ggdrConv.f'
c$$$ !---- interpolating P conv P piece of 2.12 
c$$$      double precision z 
c$$$      integer i,jz,it_p
c$$$      parameter(it_p=4)
c$$$      double precision lz(num_g),lpxp(num_g) 
c$$$      double precision Z1(it_p),LP1(it_p)
c$$$      double precision z_min,z_max,a,out_2
c$$$      logical first
c$$$      data first/.true./
c$$$
c$$$
c$$$      if(first) then 
c$$$         first=.false.
c$$$         do i=1,num_g
c$$$            lz(i)=dlog10(z_grid(i)) 
c$$$            lpxp(i)=DxPconv(i)
c$$$         enddo
c$$$      endif
c$$$      
c$$$     
c$$$      z_min=z_grid(1)
c$$$      z_max=z_grid(num_g)
c$$$
c$$$      GGdR_DxP=0d0
c$$$      if((z.gt.z_min).and.(z.lt.z_max)) then
c$$$!------find z in grid 
c$$$         a=dlog10(z) 
c$$$         call locate (lz,num_g,a,jz)
c$$$         do i=1,it_p
c$$$!------Special case when jz= 1
c$$$            if(jz.eq.1) then 
c$$$               Z1(i)=lz(i) 
c$$$               LP1(i)=lpxp(i) 
c$$$            elseif(jz.eq.(num_g-1)) then   
c$$$               Z1(i)=lz(jz-3+i)
c$$$               LP1(i)=lpxp(jz-3+i)
c$$$            elseif(jz.eq.num_g) then 
c$$$               Z1(i)=lz(jz-4+i) 
c$$$               LP1(i)=lpxp(jz-4+i) 
c$$$            else
c$$$               Z1(i)=lz(jz-2+i) 
c$$$               LP1(i)=lpxp(jz-2+i)
c$$$            endif
c$$$         enddo
c$$$      endif
c$$$!------- interpolate 
c$$$           
c$$$      
c$$$      call dpolint(Z1,LP1,4,a,GGdR_DxP,out_2)
c$$$ 
c$$$!-------- return 
c$$$      return 
c$$$      end function 
