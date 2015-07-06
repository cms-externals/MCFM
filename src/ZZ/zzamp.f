      subroutine zzamp(p1,p2,p3,p4,p5,p6,za,zb,zaa,z34,z56)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'srdiags.f'
      double precision tinv,s12,s34,s56
      double complex zaa(2,2,2),z34(2,2,2),z56(2,2,2),zab2
      integer p1,p2,p3,p4,p5,p6
C--   order of indices is gluon helicity 
      tinv(p1,p2,p3)=1d0/(s(p1,p2)+s(p2,p3)+s(p1,p3))
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      zaa(2,2,2)= + s34**(-1)*s56**(-1) * (  - za(p2,p4)*zb(p1,p5)*
     &    zab2(p6,p1,p5,p3)*tinv(p1,p5,p6) - za(p2,p6)*zb(p1,p3)*zab2(
     &    p4,p1,p3,p5)*tinv(p1,p3,p4) )

      zaa(1,2,2)= + s34**(-1)*s56**(-1) * (  - za(p1,p4)*zb(p2,p5)*
     &    zab2(p6,p1,p4,p3)*tinv(p1,p3,p4) - za(p1,p6)*zb(p2,p3)*zab2(
     &    p4,p1,p6,p5)*tinv(p1,p5,p6) )

      zaa(2,1,2)= + s34**(-1)*s56**(-1) * (  - za(p2,p3)*zb(p1,p5)*
     &    zab2(p6,p1,p5,p4)*tinv(p1,p5,p6) - za(p2,p6)*zb(p1,p4)*zab2(
     &    p3,p1,p4,p5)*tinv(p1,p3,p4) )

      zaa(2,2,1)= + s34**(-1)*s56**(-1) * (  - za(p2,p4)*zb(p1,p6)*
     &    zab2(p5,p1,p6,p3)*tinv(p1,p5,p6) - za(p2,p5)*zb(p1,p3)*zab2(
     &    p4,p1,p3,p6)*tinv(p1,p3,p4) )

      zaa(1,1,1)= + s34**(-1)*s56**(-1) * (  - za(p1,p3)*zb(p2,p6)*
     &    zab2(p5,p1,p3,p4)*tinv(p1,p3,p4) - za(p1,p5)*zb(p2,p4)*zab2(
     &    p3,p1,p5,p6)*tinv(p1,p5,p6) )

      zaa(2,1,1)= + s34**(-1)*s56**(-1) * (  - za(p2,p3)*zb(p1,p6)*
     &    zab2(p5,p1,p6,p4)*tinv(p1,p5,p6) - za(p2,p5)*zb(p1,p4)*zab2(
     &    p3,p1,p4,p6)*tinv(p1,p3,p4) )

      zaa(1,2,1)= + s34**(-1)*s56**(-1) * (  - za(p1,p4)*zb(p2,p6)*
     &    zab2(p5,p1,p4,p3)*tinv(p1,p3,p4) - za(p1,p5)*zb(p2,p3)*zab2(
     &    p4,p1,p5,p6)*tinv(p1,p5,p6) )

      zaa(1,1,2)= + s34**(-1)*s56**(-1) * (  - za(p1,p3)*zb(p2,p5)*
     &    zab2(p6,p1,p3,p4)*tinv(p1,p3,p4) - za(p1,p6)*zb(p2,p4)*zab2(
     &    p3,p1,p6,p5)*tinv(p1,p5,p6) )


c--- do not calculate single resonant diagrams unnecessarily
      if (srdiags .eqv. .false.) return

      z34(2,2,2)= + s34**(-1)*s12**(-1) * (  - za(p2,p6)*zb(p3,p5)*
     &    zab2(p4,p3,p5,p1)*tinv(p3,p4,p5) - za(p4,p6)*zb(p1,p5)*zab2(
     &    p2,p1,p5,p3)*tinv(p1,p2,p5) )

      z34(1,2,2)= + s34**(-1)*s12**(-1) * (  - za(p1,p6)*zb(p3,p5)*
     &    zab2(p4,p3,p5,p2)*tinv(p3,p4,p5) - za(p4,p6)*zb(p2,p5)*zab2(
     &    p1,p2,p5,p3)*tinv(p1,p2,p5) )

      z34(2,1,2)= + s34**(-1)*s12**(-1) * (  - za(p2,p6)*zb(p4,p5)*
     &    zab2(p3,p4,p5,p1)*tinv(p3,p4,p5) - za(p3,p6)*zb(p1,p5)*zab2(
     &    p2,p1,p5,p4)*tinv(p1,p2,p5) )

      z34(2,2,1)= + s34**(-1)*s12**(-1) * (  - za(p2,p5)*zb(p3,p6)*
     &    zab2(p4,p2,p5,p1)*tinv(p1,p2,p5) - za(p4,p5)*zb(p1,p6)*zab2(
     &    p2,p4,p5,p3)*tinv(p3,p4,p5) )

      z34(1,1,1)= + s34**(-1)*s12**(-1) * (  - za(p1,p5)*zb(p4,p6)*
     &    zab2(p3,p1,p5,p2)*tinv(p1,p2,p5) - za(p3,p5)*zb(p2,p6)*zab2(
     &    p1,p3,p5,p4)*tinv(p3,p4,p5) )

      z34(2,1,1)= + s34**(-1)*s12**(-1) * (  - za(p2,p5)*zb(p4,p6)*
     &    zab2(p3,p2,p5,p1)*tinv(p1,p2,p5) - za(p3,p5)*zb(p1,p6)*zab2(
     &    p2,p3,p5,p4)*tinv(p3,p4,p5) )

      z34(1,2,1)= + s34**(-1)*s12**(-1) * (  - za(p1,p5)*zb(p3,p6)*
     &    zab2(p4,p1,p5,p2)*tinv(p1,p2,p5) - za(p4,p5)*zb(p2,p6)*zab2(
     &    p1,p4,p5,p3)*tinv(p3,p4,p5) )

      z34(1,1,2)= + s34**(-1)*s12**(-1) * (  - za(p1,p6)*zb(p4,p5)*
     &    zab2(p3,p4,p5,p2)*tinv(p3,p4,p5) - za(p3,p6)*zb(p2,p5)*zab2(
     &    p1,p2,p5,p4)*tinv(p1,p2,p5) )

      z56(2,2,2)= + s56**(-1)*s12**(-1) * ( za(p2,p4)*zb(p3,p5)*zab2(p6
     &    ,p3,p5,p1)*tinv(p3,p5,p6) + za(p4,p6)*zb(p1,p3)*zab2(p2,p1,p3
     &    ,p5)*tinv(p1,p2,p3) )

      z56(1,2,2)= + s56**(-1)*s12**(-1) * ( za(p1,p4)*zb(p3,p5)*zab2(p6
     &    ,p3,p5,p2)*tinv(p3,p5,p6) + za(p4,p6)*zb(p2,p3)*zab2(p1,p2,p3
     &    ,p5)*tinv(p1,p2,p3) )

      z56(2,1,2)= + s56**(-1)*s12**(-1) * ( za(p2,p3)*zb(p4,p5)*zab2(p6
     &    ,p2,p3,p1)*tinv(p1,p2,p3) + za(p3,p6)*zb(p1,p4)*zab2(p2,p3,p6
     &    ,p5)*tinv(p3,p5,p6) )

      z56(2,2,1)= + s56**(-1)*s12**(-1) * ( za(p2,p4)*zb(p3,p6)*zab2(p5
     &    ,p3,p6,p1)*tinv(p3,p5,p6) + za(p4,p5)*zb(p1,p3)*zab2(p2,p1,p3
     &    ,p6)*tinv(p1,p2,p3) )

      z56(1,1,1)= + s56**(-1)*s12**(-1) * ( za(p1,p3)*zb(p4,p6)*zab2(p5
     &    ,p1,p3,p2)*tinv(p1,p2,p3) + za(p3,p5)*zb(p2,p4)*zab2(p1,p3,p5
     &    ,p6)*tinv(p3,p5,p6) )

      z56(2,1,1)= + s56**(-1)*s12**(-1) * ( za(p2,p3)*zb(p4,p6)*zab2(p5
     &    ,p2,p3,p1)*tinv(p1,p2,p3) + za(p3,p5)*zb(p1,p4)*zab2(p2,p3,p5
     &    ,p6)*tinv(p3,p5,p6) )

      z56(1,2,1)= + s56**(-1)*s12**(-1) * ( za(p1,p4)*zb(p3,p6)*zab2(p5
     &    ,p3,p6,p2)*tinv(p3,p5,p6) + za(p4,p5)*zb(p2,p3)*zab2(p1,p2,p3
     &    ,p6)*tinv(p1,p2,p3) )

      z56(1,1,2)= + s56**(-1)*s12**(-1) * ( za(p1,p3)*zb(p4,p5)*zab2(p6
     &    ,p1,p3,p2)*tinv(p1,p2,p3) + za(p3,p6)*zb(p2,p4)*zab2(p1,p3,p6
     &    ,p5)*tinv(p3,p5,p6) )

      return
      end
