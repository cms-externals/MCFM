      subroutine gg_ZZ(p,msqgg)
      implicit none
c--- Author: J. M. Campbell, March 2011
c--- For now, work in the approximation of BDK, i.e. that we
c--- retain only 1/mt^2 for top quark loops
c--- Box contributions are then complete (terms of order 1/mt^4 discarded)
c--- No triangle contributions

      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer h1,h2,h34,h56,up,dn
      double precision p(mxpart,4),msqgg,fac,cvec(2),cax(2),
     & pttwo,ptWsafetycut_massless,cl1(2),cl2(2)
      double complex Avec(2,2,2,2),prop34,prop56
      double complex a64v
      parameter(up=2,dn=1)

c--- omit massless loops for pt(W) < "ptWsafetycut_massless" (for num. stability)
c      ptWsafetycut_massless=0.05d0
c--- omit massless loops for pt(W) < "ptWsafetycut_massless" (for num. stability)
      ptWsafetycut_massless=2d0

c--- set up spinor products      
      call spinoru(6,p,za,zb)

      
c--- fill amplitudes
c--- labels are: helicity of gluons, lepton 3 and lepton 5

      if (pttwo(3,4,p) .lt. ptWsafetycut_massless) then
c--- ensure numerical stability: set massless loops to zero
c--- for pt(W) < "ptWsafetycut_massless" GeV
        do h1=1,2
        do h2=1,2
        do h34=1,2
        do h56=1,2
          Avec(h1,h2,h34,h56)=czip
        enddo
        enddo
        enddo
        enddo
      else
      Avec(2,2,1,1)=a64v('q+qb-g-g-',3,4,1,2,6,5,zb,za)*(-im)
      Avec(2,1,1,1)=a64v('q+qb-g-g+',3,4,1,2,6,5,zb,za)*(-im)
      Avec(1,2,1,1)=a64v('q+qb-g+g-',3,4,1,2,6,5,zb,za)*(-im)
      Avec(1,1,1,1)=a64v('q+qb-g+g+',3,4,1,2,6,5,zb,za)*(-im)

      Avec(2,2,2,1)=a64v('q+qb-g+g+',3,4,1,2,5,6,za,zb)*(-im)
      Avec(2,1,2,1)=a64v('q+qb-g+g-',3,4,1,2,5,6,za,zb)*(-im)
      Avec(1,2,2,1)=a64v('q+qb-g-g+',3,4,1,2,5,6,za,zb)*(-im)
      Avec(1,1,2,1)=a64v('q+qb-g-g-',3,4,1,2,5,6,za,zb)*(-im)

      Avec(2,2,1,2)=a64v('q+qb-g-g-',3,4,1,2,5,6,zb,za)*(-im)
      Avec(2,1,1,2)=a64v('q+qb-g-g+',3,4,1,2,5,6,zb,za)*(-im)
      Avec(1,2,1,2)=a64v('q+qb-g+g-',3,4,1,2,5,6,zb,za)*(-im)
      Avec(1,1,1,2)=a64v('q+qb-g+g+',3,4,1,2,5,6,zb,za)*(-im)

      Avec(2,2,2,2)=a64v('q+qb-g+g+',3,4,1,2,6,5,za,zb)*(-im)
      Avec(2,1,2,2)=a64v('q+qb-g+g-',3,4,1,2,6,5,za,zb)*(-im)
      Avec(1,2,2,2)=a64v('q+qb-g-g+',3,4,1,2,6,5,za,zb)*(-im)
      Avec(1,1,2,2)=a64v('q+qb-g-g-',3,4,1,2,6,5,za,zb)*(-im)
      endif


c--- propagator factors
      prop34=s(3,4)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=s(5,6)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)

c--- left and right-handed lepton couplings as an array
c--- cl1 associated with Z(3+4), cl2 associated with Z(5+6)
      cl1(1)=l1 
      cl1(2)=r1
      cl2(1)=l2
      cl2(2)=r2

c--- vector and axial couplings as an array for up/down quarks
      cvec(up)=half*(L(up)+R(up))
      cvec(dn)=half*(L(dn)+R(dn))
      cax(up)=half*(L(up)-R(up))
      cax(dn)=half*(L(dn)-R(dn))
      
c--- assume 5 massless flavors in the loop, (2 up, 3 down)
      msqgg=0d0
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2
      msqgg=msqgg+cdabs(Avec(h1,h2,h34,h56)
c--- vector couplings in the loop
     & *(2d0*(Qu*q1+cvec(up)*cl1(h34)*prop34)
     &      *(Qu*q1+cvec(up)*cl2(h56)*prop56)
     &  +3d0*(Qd*q1+cvec(dn)*cl1(h34)*prop34)
     &      *(Qd*q1+cvec(dn)*cl2(h56)*prop56)
c--- axial couplings in the loop
     &  +2d0*(cax(up)*cl1(h34)*prop34)*(cax(up)*cl2(h56)*prop56)
     &  +3d0*(cax(dn)*cl1(h34)*prop34)*(cax(dn)*cl2(h56)*prop56)
     & ))**2
      enddo
      enddo
      enddo
      enddo
      
c--- overall factor from diagrams
      fac=avegg*V*(2d0*4d0*esq*gsq/(16d0*pisq)*esq)**2
      
      msqgg=msqgg*fac
      
      return
      end
      
      
