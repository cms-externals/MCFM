      subroutine zzgamp(p1,p2,p3,p4,p5,p6,p7,za,zb,zaa,z34,z56)
c--- This is the new code for the amplitudes for ZZ+gluon production
c--- (singly-resonant diagrams are included)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'srdiags.f'
      double precision tinv,s34,s56,tinv346,tinv356,tinv127,
     & tinv456,tinv345,tinv156,tinv256,tinv134,tinv234
      double complex zaa(2,2,2,2),z34(2,2,2,2),z56(2,2,2,2),zab2
      double complex zaa1727,zbb1727
      integer p1,p2,p3,p4,p5,p6,p7
C--   order of indices is gluon helicity 
      tinv(p1,p2,p3)=1d0/(s(p1,p2)+s(p2,p3)+s(p1,p3))
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s34=s(p3,p4)
      s56=s(p5,p6)
      zaa1727=za(p1,p7)*za(p2,p7)
      zbb1727=zb(p1,p7)*zb(p2,p7)
      tinv127=tinv(p1,p2,p7)
      tinv134=tinv(p1,p3,p4)
      tinv234=tinv(p2,p3,p4)
      tinv156=tinv(p1,p5,p6)
      tinv256=tinv(p2,p5,p6)
      tinv356=tinv(p3,p5,p6)
      tinv456=tinv(p4,p5,p6)
      tinv346=tinv(p3,p4,p6)
      tinv345=tinv(p3,p4,p5)
      zaa(2,2,2,2)= + zaa1727**(-1)*s34**(-1)*s56**(-1) * (  - za(p1,p7
     &    )*za(p2,p4)**2*zb(p1,p5)*zb(p3,p4)*zab2(p6,p1,p5,p7)*tinv156*
     &    tinv234 - za(p1,p7)*za(p2,p6)**2*zb(p1,p3)*zb(p5,p6)*zab2(p4,
     &    p1,p3,p7)*tinv256*tinv134 - za(p2,p4)*zab2(p2,p1,p7,p5)*zab2(
     &    p6,p2,p4,p3)*tinv234 - za(p2,p6)*zab2(p2,p1,p7,p3)*zab2(p4,p2
     &    ,p6,p5)*tinv256 )

      zaa(1,1,1,1)= + zbb1727**(-1)*s34**(-1)*s56**(-1) * ( za(p1,p3)*
     &    za(p5,p6)*zb(p1,p7)*zb(p2,p6)**2*zab2(p7,p1,p3,p4)*tinv256*
     &    tinv134 + za(p1,p5)*za(p3,p4)*zb(p1,p7)*zb(p2,p4)**2*zab2(p7,
     &    p1,p5,p6)*tinv156*tinv234 + zb(p2,p4)*zab2(p3,p2,p4,p6)*zab2(
     &    p5,p1,p7,p2)*tinv234 + zb(p2,p6)*zab2(p3,p1,p7,p2)*zab2(p5,p2
     &    ,p6,p4)*tinv256 )

      zaa(1,2,2,2)= + zaa1727**(-1)*s34**(-1)*s56**(-1) * (  - za(p1,p4
     &    )**2*za(p2,p7)*zb(p2,p5)*zb(p3,p4)*zab2(p6,p2,p5,p7)*tinv256*
     &    tinv134 - za(p1,p4)*zab2(p1,p2,p7,p5)*zab2(p6,p1,p4,p3)*
     &    tinv134 - za(p1,p6)**2*za(p2,p7)*zb(p2,p3)*zb(p5,p6)*zab2(p4,
     &    p2,p3,p7)*tinv156*tinv234 - za(p1,p6)*zab2(p1,p2,p7,p3)*zab2(
     &    p4,p1,p6,p5)*tinv156 )

      zaa(2,1,2,2)= + zaa1727**(-1)*s34**(-1)*s56**(-1) * ( za(p1,p7)*
     &    za(p2,p3)**2*zb(p1,p5)*zb(p3,p4)*zab2(p6,p1,p5,p7)*tinv156*
     &    tinv234 - za(p1,p7)*za(p2,p6)**2*zb(p1,p4)*zb(p5,p6)*zab2(p3,
     &    p1,p4,p7)*tinv256*tinv134 - za(p2,p3)*zab2(p2,p1,p7,p5)*zab2(
     &    p6,p2,p3,p4)*tinv234 - za(p2,p6)*zab2(p2,p1,p7,p4)*zab2(p3,p2
     &    ,p6,p5)*tinv256 )

      zaa(2,2,1,2)= + zaa1727**(-1)*s34**(-1)*s56**(-1) * (  - za(p1,p7
     &    )*za(p2,p4)**2*zb(p1,p6)*zb(p3,p4)*zab2(p5,p1,p6,p7)*tinv156*
     &    tinv234 + za(p1,p7)*za(p2,p5)**2*zb(p1,p3)*zb(p5,p6)*zab2(p4,
     &    p1,p3,p7)*tinv256*tinv134 - za(p2,p4)*zab2(p2,p1,p7,p6)*zab2(
     &    p5,p2,p4,p3)*tinv234 - za(p2,p5)*zab2(p2,p1,p7,p3)*zab2(p4,p2
     &    ,p5,p6)*tinv256 )

      zaa(2,2,2,1)= + zbb1727**(-1)*s34**(-1)*s56**(-1) * (  - za(p2,p4
     &    )*za(p5,p6)*zb(p1,p5)**2*zb(p2,p7)*zab2(p7,p2,p4,p3)*tinv156*
     &    tinv234 - za(p2,p6)*za(p3,p4)*zb(p1,p3)**2*zb(p2,p7)*zab2(p7,
     &    p2,p6,p5)*tinv256*tinv134 + zb(p1,p3)*zab2(p4,p1,p3,p5)*zab2(
     &    p6,p2,p7,p1)*tinv134 + zb(p1,p5)*zab2(p4,p2,p7,p1)*zab2(p6,p1
     &    ,p5,p3)*tinv156 )

      zaa(1,1,2,2)= + zaa1727**(-1)*s34**(-1)*s56**(-1) * ( za(p1,p3)**
     &    2*za(p2,p7)*zb(p2,p5)*zb(p3,p4)*zab2(p6,p2,p5,p7)*tinv256*
     &    tinv134 - za(p1,p3)*zab2(p1,p2,p7,p5)*zab2(p6,p1,p3,p4)*
     &    tinv134 - za(p1,p6)**2*za(p2,p7)*zb(p2,p4)*zb(p5,p6)*zab2(p3,
     &    p2,p4,p7)*tinv156*tinv234 - za(p1,p6)*zab2(p1,p2,p7,p4)*zab2(
     &    p3,p1,p6,p5)*tinv156 )

      zaa(1,2,1,2)= + zaa1727**(-1)*s34**(-1)*s56**(-1) * (  - za(p1,p4
     &    )**2*za(p2,p7)*zb(p2,p6)*zb(p3,p4)*zab2(p5,p2,p6,p7)*tinv256*
     &    tinv134 - za(p1,p4)*zab2(p1,p2,p7,p6)*zab2(p5,p1,p4,p3)*
     &    tinv134 + za(p1,p5)**2*za(p2,p7)*zb(p2,p3)*zb(p5,p6)*zab2(p4,
     &    p2,p3,p7)*tinv156*tinv234 - za(p1,p5)*zab2(p1,p2,p7,p3)*zab2(
     &    p4,p1,p5,p6)*tinv156 )

      zaa(1,2,2,1)= + zbb1727**(-1)*s34**(-1)*s56**(-1) * (  - za(p1,p4
     &    )*za(p5,p6)*zb(p1,p7)*zb(p2,p5)**2*zab2(p7,p1,p4,p3)*tinv256*
     &    tinv134 - za(p1,p6)*za(p3,p4)*zb(p1,p7)*zb(p2,p3)**2*zab2(p7,
     &    p1,p6,p5)*tinv156*tinv234 + zb(p2,p3)*zab2(p4,p2,p3,p5)*zab2(
     &    p6,p1,p7,p2)*tinv234 + zb(p2,p5)*zab2(p4,p1,p7,p2)*zab2(p6,p2
     &    ,p5,p3)*tinv256 )

      zaa(2,1,1,2)= + zaa1727**(-1)*s34**(-1)*s56**(-1) * ( za(p1,p7)*
     &    za(p2,p3)**2*zb(p1,p6)*zb(p3,p4)*zab2(p5,p1,p6,p7)*tinv156*
     &    tinv234 + za(p1,p7)*za(p2,p5)**2*zb(p1,p4)*zb(p5,p6)*zab2(p3,
     &    p1,p4,p7)*tinv256*tinv134 - za(p2,p3)*zab2(p2,p1,p7,p6)*zab2(
     &    p5,p2,p3,p4)*tinv234 - za(p2,p5)*zab2(p2,p1,p7,p4)*zab2(p3,p2
     &    ,p5,p6)*tinv256 )

      zaa(2,1,2,1)= + zbb1727**(-1)*s34**(-1)*s56**(-1) * (  - za(p2,p3
     &    )*za(p5,p6)*zb(p1,p5)**2*zb(p2,p7)*zab2(p7,p2,p3,p4)*tinv156*
     &    tinv234 + za(p2,p6)*za(p3,p4)*zb(p1,p4)**2*zb(p2,p7)*zab2(p7,
     &    p2,p6,p5)*tinv256*tinv134 + zb(p1,p4)*zab2(p3,p1,p4,p5)*zab2(
     &    p6,p2,p7,p1)*tinv134 + zb(p1,p5)*zab2(p3,p2,p7,p1)*zab2(p6,p1
     &    ,p5,p4)*tinv156 )

      zaa(2,2,1,1)= + zbb1727**(-1)*s34**(-1)*s56**(-1) * ( za(p2,p4)*
     &    za(p5,p6)*zb(p1,p6)**2*zb(p2,p7)*zab2(p7,p2,p4,p3)*tinv156*
     &    tinv234 - za(p2,p5)*za(p3,p4)*zb(p1,p3)**2*zb(p2,p7)*zab2(p7,
     &    p2,p5,p6)*tinv256*tinv134 + zb(p1,p3)*zab2(p4,p1,p3,p6)*zab2(
     &    p5,p2,p7,p1)*tinv134 + zb(p1,p6)*zab2(p4,p2,p7,p1)*zab2(p5,p1
     &    ,p6,p3)*tinv156 )

      zaa(2,1,1,1)= + zbb1727**(-1)*s34**(-1)*s56**(-1) * ( za(p2,p3)*
     &    za(p5,p6)*zb(p1,p6)**2*zb(p2,p7)*zab2(p7,p2,p3,p4)*tinv156*
     &    tinv234 + za(p2,p5)*za(p3,p4)*zb(p1,p4)**2*zb(p2,p7)*zab2(p7,
     &    p2,p5,p6)*tinv256*tinv134 + zb(p1,p4)*zab2(p3,p1,p4,p6)*zab2(
     &    p5,p2,p7,p1)*tinv134 + zb(p1,p6)*zab2(p3,p2,p7,p1)*zab2(p5,p1
     &    ,p6,p4)*tinv156 )

      zaa(1,2,1,1)= + zbb1727**(-1)*s34**(-1)*s56**(-1) * ( za(p1,p4)*
     &    za(p5,p6)*zb(p1,p7)*zb(p2,p6)**2*zab2(p7,p1,p4,p3)*tinv256*
     &    tinv134 - za(p1,p5)*za(p3,p4)*zb(p1,p7)*zb(p2,p3)**2*zab2(p7,
     &    p1,p5,p6)*tinv156*tinv234 + zb(p2,p3)*zab2(p4,p2,p3,p6)*zab2(
     &    p5,p1,p7,p2)*tinv234 + zb(p2,p6)*zab2(p4,p1,p7,p2)*zab2(p5,p2
     &    ,p6,p3)*tinv256 )

      zaa(1,1,2,1)= + zbb1727**(-1)*s34**(-1)*s56**(-1) * (  - za(p1,p3
     &    )*za(p5,p6)*zb(p1,p7)*zb(p2,p5)**2*zab2(p7,p1,p3,p4)*tinv256*
     &    tinv134 + za(p1,p6)*za(p3,p4)*zb(p1,p7)*zb(p2,p4)**2*zab2(p7,
     &    p1,p6,p5)*tinv156*tinv234 + zb(p2,p4)*zab2(p3,p2,p4,p5)*zab2(
     &    p6,p1,p7,p2)*tinv234 + zb(p2,p5)*zab2(p3,p1,p7,p2)*zab2(p6,p2
     &    ,p5,p4)*tinv256 )

      zaa(1,1,1,2)= + zaa1727**(-1)*s34**(-1)*s56**(-1) * ( za(p1,p3)**
     &    2*za(p2,p7)*zb(p2,p6)*zb(p3,p4)*zab2(p5,p2,p6,p7)*tinv256*
     &    tinv134 - za(p1,p3)*zab2(p1,p2,p7,p6)*zab2(p5,p1,p3,p4)*
     &    tinv134 + za(p1,p5)**2*za(p2,p7)*zb(p2,p4)*zb(p5,p6)*zab2(p3,
     &    p2,p4,p7)*tinv156*tinv234 - za(p1,p5)*zab2(p1,p2,p7,p4)*zab2(
     &    p3,p1,p5,p6)*tinv156 )


c--- do not calculate single resonant diagrams unnecessarily
      if (srdiags .eqv. .false.) return

      z34(2,2,2,2)= + zaa1727**(-1)*s34**(-1) * (  - za(p1,p2)*za(p2,p4
     &    )*za(p4,p6)*zb(p1,p5)*zb(p3,p4)*tinv127*tinv346 - za(p1,p2)*
     &    za(p2,p6)*za(p3,p4)*zb(p1,p3)*zb(p3,p5)*tinv127*tinv345 + za(
     &    p1,p2)*za(p2,p6)*za(p4,p5)*zb(p1,p5)*zb(p3,p5)*tinv127*
     &    tinv345 - za(p1,p2)*za(p2,p6)*za(p4,p6)*zb(p1,p5)*zb(p3,p6)*
     &    tinv127*tinv346 - za(p2,p4)*za(p2,p7)*za(p4,p6)*zb(p3,p4)*zb(
     &    p5,p7)*tinv127*tinv346 - za(p2,p6)*za(p2,p7)*za(p3,p4)*zb(p3,
     &    p5)*zb(p3,p7)*tinv127*tinv345 + za(p2,p6)*za(p2,p7)*za(p4,p5)
     &    *zb(p3,p5)*zb(p5,p7)*tinv127*tinv345 - za(p2,p6)*za(p2,p7)*
     &    za(p4,p6)*zb(p3,p6)*zb(p5,p7)*tinv127*tinv346 )

      z34(1,1,1,1)= + zbb1727**(-1)*s34**(-1) * ( za(p1,p3)*za(p3,p5)*
     &    zb(p1,p2)*zb(p2,p6)*zb(p3,p4)*tinv127*tinv345 + za(p1,p5)*za(
     &    p3,p4)*zb(p1,p2)*zb(p2,p4)*zb(p4,p6)*tinv127*tinv346 - za(p1,
     &    p5)*za(p3,p5)*zb(p1,p2)*zb(p2,p6)*zb(p4,p5)*tinv127*tinv345
     &     + za(p1,p5)*za(p3,p6)*zb(p1,p2)*zb(p2,p6)*zb(p4,p6)*tinv127*
     &    tinv346 + za(p3,p4)*za(p5,p7)*zb(p2,p4)*zb(p2,p7)*zb(p4,p6)*
     &    tinv127*tinv346 + za(p3,p5)*za(p3,p7)*zb(p2,p6)*zb(p2,p7)*zb(
     &    p3,p4)*tinv127*tinv345 - za(p3,p5)*za(p5,p7)*zb(p2,p6)*zb(p2,
     &    p7)*zb(p4,p5)*tinv127*tinv345 + za(p3,p6)*za(p5,p7)*zb(p2,p6)
     &    *zb(p2,p7)*zb(p4,p6)*tinv127*tinv346 )

      z34(1,2,2,2)= + zaa1727**(-1)*s34**(-1) * (  - za(p1,p2)*za(p1,p4
     &    )*za(p4,p6)*zb(p2,p5)*zb(p3,p4)*tinv127*tinv346 - za(p1,p2)*
     &    za(p1,p6)*za(p3,p4)*zb(p2,p3)*zb(p3,p5)*tinv127*tinv345 + za(
     &    p1,p2)*za(p1,p6)*za(p4,p5)*zb(p2,p5)*zb(p3,p5)*tinv127*
     &    tinv345 - za(p1,p2)*za(p1,p6)*za(p4,p6)*zb(p2,p5)*zb(p3,p6)*
     &    tinv127*tinv346 + za(p1,p4)*za(p1,p7)*za(p4,p6)*zb(p3,p4)*zb(
     &    p5,p7)*tinv127*tinv346 + za(p1,p6)*za(p1,p7)*za(p3,p4)*zb(p3,
     &    p5)*zb(p3,p7)*tinv127*tinv345 - za(p1,p6)*za(p1,p7)*za(p4,p5)
     &    *zb(p3,p5)*zb(p5,p7)*tinv127*tinv345 + za(p1,p6)*za(p1,p7)*
     &    za(p4,p6)*zb(p3,p6)*zb(p5,p7)*tinv127*tinv346 )

      z34(2,1,2,2)= + zaa1727**(-1)*s34**(-1) * ( za(p1,p2)*za(p2,p3)*
     &    za(p3,p6)*zb(p1,p5)*zb(p3,p4)*tinv127*tinv346 + za(p1,p2)*za(
     &    p2,p6)*za(p3,p4)*zb(p1,p4)*zb(p4,p5)*tinv127*tinv345 + za(p1,
     &    p2)*za(p2,p6)*za(p3,p5)*zb(p1,p5)*zb(p4,p5)*tinv127*tinv345
     &     - za(p1,p2)*za(p2,p6)*za(p3,p6)*zb(p1,p5)*zb(p4,p6)*tinv127*
     &    tinv346 + za(p2,p3)*za(p2,p7)*za(p3,p6)*zb(p3,p4)*zb(p5,p7)*
     &    tinv127*tinv346 + za(p2,p6)*za(p2,p7)*za(p3,p4)*zb(p4,p5)*zb(
     &    p4,p7)*tinv127*tinv345 + za(p2,p6)*za(p2,p7)*za(p3,p5)*zb(p4,
     &    p5)*zb(p5,p7)*tinv127*tinv345 - za(p2,p6)*za(p2,p7)*za(p3,p6)
     &    *zb(p4,p6)*zb(p5,p7)*tinv127*tinv346 )

      z34(2,2,1,2)= + zaa1727**(-1)*s34**(-1) * ( za(p1,p2)*za(p2,p4)*
     &    za(p4,p5)*zb(p1,p6)*zb(p3,p4)*tinv127*tinv345 + za(p1,p2)*za(
     &    p2,p5)*za(p3,p4)*zb(p1,p3)*zb(p3,p6)*tinv127*tinv346 + za(p1,
     &    p2)*za(p2,p5)*za(p4,p5)*zb(p1,p6)*zb(p3,p5)*tinv127*tinv345
     &     - za(p1,p2)*za(p2,p5)*za(p4,p6)*zb(p1,p6)*zb(p3,p6)*tinv127*
     &    tinv346 + za(p2,p4)*za(p2,p7)*za(p4,p5)*zb(p3,p4)*zb(p6,p7)*
     &    tinv127*tinv345 + za(p2,p5)*za(p2,p7)*za(p3,p4)*zb(p3,p6)*zb(
     &    p3,p7)*tinv127*tinv346 + za(p2,p5)*za(p2,p7)*za(p4,p5)*zb(p3,
     &    p5)*zb(p6,p7)*tinv127*tinv345 - za(p2,p5)*za(p2,p7)*za(p4,p6)
     &    *zb(p3,p6)*zb(p6,p7)*tinv127*tinv346 )

      z34(2,2,2,1)= + zbb1727**(-1)*s34**(-1) * ( za(p2,p4)*za(p4,p6)*
     &    zb(p1,p2)*zb(p1,p5)*zb(p3,p4)*tinv127*tinv346 + za(p2,p6)*za(
     &    p3,p4)*zb(p1,p2)*zb(p1,p3)*zb(p3,p5)*tinv127*tinv345 - za(p2,
     &    p6)*za(p4,p5)*zb(p1,p2)*zb(p1,p5)*zb(p3,p5)*tinv127*tinv345
     &     + za(p2,p6)*za(p4,p6)*zb(p1,p2)*zb(p1,p5)*zb(p3,p6)*tinv127*
     &    tinv346 - za(p3,p4)*za(p6,p7)*zb(p1,p3)*zb(p1,p7)*zb(p3,p5)*
     &    tinv127*tinv345 + za(p4,p5)*za(p6,p7)*zb(p1,p5)*zb(p1,p7)*zb(
     &    p3,p5)*tinv127*tinv345 - za(p4,p6)*za(p4,p7)*zb(p1,p5)*zb(p1,
     &    p7)*zb(p3,p4)*tinv127*tinv346 - za(p4,p6)*za(p6,p7)*zb(p1,p5)
     &    *zb(p1,p7)*zb(p3,p6)*tinv127*tinv346 )

      z34(1,1,2,2)= + zaa1727**(-1)*s34**(-1) * ( za(p1,p2)*za(p1,p3)*
     &    za(p3,p6)*zb(p2,p5)*zb(p3,p4)*tinv127*tinv346 + za(p1,p2)*za(
     &    p1,p6)*za(p3,p4)*zb(p2,p4)*zb(p4,p5)*tinv127*tinv345 + za(p1,
     &    p2)*za(p1,p6)*za(p3,p5)*zb(p2,p5)*zb(p4,p5)*tinv127*tinv345
     &     - za(p1,p2)*za(p1,p6)*za(p3,p6)*zb(p2,p5)*zb(p4,p6)*tinv127*
     &    tinv346 - za(p1,p3)*za(p1,p7)*za(p3,p6)*zb(p3,p4)*zb(p5,p7)*
     &    tinv127*tinv346 - za(p1,p6)*za(p1,p7)*za(p3,p4)*zb(p4,p5)*zb(
     &    p4,p7)*tinv127*tinv345 - za(p1,p6)*za(p1,p7)*za(p3,p5)*zb(p4,
     &    p5)*zb(p5,p7)*tinv127*tinv345 + za(p1,p6)*za(p1,p7)*za(p3,p6)
     &    *zb(p4,p6)*zb(p5,p7)*tinv127*tinv346 )

      z34(1,2,1,2)= + zaa1727**(-1)*s34**(-1) * ( za(p1,p2)*za(p1,p4)*
     &    za(p4,p5)*zb(p2,p6)*zb(p3,p4)*tinv127*tinv345 + za(p1,p2)*za(
     &    p1,p5)*za(p3,p4)*zb(p2,p3)*zb(p3,p6)*tinv127*tinv346 + za(p1,
     &    p2)*za(p1,p5)*za(p4,p5)*zb(p2,p6)*zb(p3,p5)*tinv127*tinv345
     &     - za(p1,p2)*za(p1,p5)*za(p4,p6)*zb(p2,p6)*zb(p3,p6)*tinv127*
     &    tinv346 - za(p1,p4)*za(p1,p7)*za(p4,p5)*zb(p3,p4)*zb(p6,p7)*
     &    tinv127*tinv345 - za(p1,p5)*za(p1,p7)*za(p3,p4)*zb(p3,p6)*zb(
     &    p3,p7)*tinv127*tinv346 - za(p1,p5)*za(p1,p7)*za(p4,p5)*zb(p3,
     &    p5)*zb(p6,p7)*tinv127*tinv345 + za(p1,p5)*za(p1,p7)*za(p4,p6)
     &    *zb(p3,p6)*zb(p6,p7)*tinv127*tinv346 )

      z34(1,2,2,1)= + zbb1727**(-1)*s34**(-1) * ( za(p1,p4)*za(p4,p6)*
     &    zb(p1,p2)*zb(p2,p5)*zb(p3,p4)*tinv127*tinv346 + za(p1,p6)*za(
     &    p3,p4)*zb(p1,p2)*zb(p2,p3)*zb(p3,p5)*tinv127*tinv345 - za(p1,
     &    p6)*za(p4,p5)*zb(p1,p2)*zb(p2,p5)*zb(p3,p5)*tinv127*tinv345
     &     + za(p1,p6)*za(p4,p6)*zb(p1,p2)*zb(p2,p5)*zb(p3,p6)*tinv127*
     &    tinv346 + za(p3,p4)*za(p6,p7)*zb(p2,p3)*zb(p2,p7)*zb(p3,p5)*
     &    tinv127*tinv345 - za(p4,p5)*za(p6,p7)*zb(p2,p5)*zb(p2,p7)*zb(
     &    p3,p5)*tinv127*tinv345 + za(p4,p6)*za(p4,p7)*zb(p2,p5)*zb(p2,
     &    p7)*zb(p3,p4)*tinv127*tinv346 + za(p4,p6)*za(p6,p7)*zb(p2,p5)
     &    *zb(p2,p7)*zb(p3,p6)*tinv127*tinv346 )

      z34(2,1,1,2)= + zaa1727**(-1)*s34**(-1) * (  - za(p1,p2)*za(p2,p3
     &    )*za(p3,p5)*zb(p1,p6)*zb(p3,p4)*tinv127*tinv345 - za(p1,p2)*
     &    za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)*tinv127*tinv346 + za(
     &    p1,p2)*za(p2,p5)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*tinv127*
     &    tinv345 - za(p1,p2)*za(p2,p5)*za(p3,p6)*zb(p1,p6)*zb(p4,p6)*
     &    tinv127*tinv346 - za(p2,p3)*za(p2,p7)*za(p3,p5)*zb(p3,p4)*zb(
     &    p6,p7)*tinv127*tinv345 - za(p2,p5)*za(p2,p7)*za(p3,p4)*zb(p4,
     &    p6)*zb(p4,p7)*tinv127*tinv346 + za(p2,p5)*za(p2,p7)*za(p3,p5)
     &    *zb(p4,p5)*zb(p6,p7)*tinv127*tinv345 - za(p2,p5)*za(p2,p7)*
     &    za(p3,p6)*zb(p4,p6)*zb(p6,p7)*tinv127*tinv346 )

      z34(2,1,2,1)= + zbb1727**(-1)*s34**(-1) * (  - za(p2,p3)*za(p3,p6
     &    )*zb(p1,p2)*zb(p1,p5)*zb(p3,p4)*tinv127*tinv346 - za(p2,p6)*
     &    za(p3,p4)*zb(p1,p2)*zb(p1,p4)*zb(p4,p5)*tinv127*tinv345 - za(
     &    p2,p6)*za(p3,p5)*zb(p1,p2)*zb(p1,p5)*zb(p4,p5)*tinv127*
     &    tinv345 + za(p2,p6)*za(p3,p6)*zb(p1,p2)*zb(p1,p5)*zb(p4,p6)*
     &    tinv127*tinv346 + za(p3,p4)*za(p6,p7)*zb(p1,p4)*zb(p1,p7)*zb(
     &    p4,p5)*tinv127*tinv345 + za(p3,p5)*za(p6,p7)*zb(p1,p5)*zb(p1,
     &    p7)*zb(p4,p5)*tinv127*tinv345 + za(p3,p6)*za(p3,p7)*zb(p1,p5)
     &    *zb(p1,p7)*zb(p3,p4)*tinv127*tinv346 - za(p3,p6)*za(p6,p7)*
     &    zb(p1,p5)*zb(p1,p7)*zb(p4,p6)*tinv127*tinv346 )

      z34(2,2,1,1)= + zbb1727**(-1)*s34**(-1) * (  - za(p2,p4)*za(p4,p5
     &    )*zb(p1,p2)*zb(p1,p6)*zb(p3,p4)*tinv127*tinv345 - za(p2,p5)*
     &    za(p3,p4)*zb(p1,p2)*zb(p1,p3)*zb(p3,p6)*tinv127*tinv346 - za(
     &    p2,p5)*za(p4,p5)*zb(p1,p2)*zb(p1,p6)*zb(p3,p5)*tinv127*
     &    tinv345 + za(p2,p5)*za(p4,p6)*zb(p1,p2)*zb(p1,p6)*zb(p3,p6)*
     &    tinv127*tinv346 + za(p3,p4)*za(p5,p7)*zb(p1,p3)*zb(p1,p7)*zb(
     &    p3,p6)*tinv127*tinv346 + za(p4,p5)*za(p4,p7)*zb(p1,p6)*zb(p1,
     &    p7)*zb(p3,p4)*tinv127*tinv345 + za(p4,p5)*za(p5,p7)*zb(p1,p6)
     &    *zb(p1,p7)*zb(p3,p5)*tinv127*tinv345 - za(p4,p6)*za(p5,p7)*
     &    zb(p1,p6)*zb(p1,p7)*zb(p3,p6)*tinv127*tinv346 )

      z34(2,1,1,1)= + zbb1727**(-1)*s34**(-1) * ( za(p2,p3)*za(p3,p5)*
     &    zb(p1,p2)*zb(p1,p6)*zb(p3,p4)*tinv127*tinv345 + za(p2,p5)*za(
     &    p3,p4)*zb(p1,p2)*zb(p1,p4)*zb(p4,p6)*tinv127*tinv346 - za(p2,
     &    p5)*za(p3,p5)*zb(p1,p2)*zb(p1,p6)*zb(p4,p5)*tinv127*tinv345
     &     + za(p2,p5)*za(p3,p6)*zb(p1,p2)*zb(p1,p6)*zb(p4,p6)*tinv127*
     &    tinv346 - za(p3,p4)*za(p5,p7)*zb(p1,p4)*zb(p1,p7)*zb(p4,p6)*
     &    tinv127*tinv346 - za(p3,p5)*za(p3,p7)*zb(p1,p6)*zb(p1,p7)*zb(
     &    p3,p4)*tinv127*tinv345 + za(p3,p5)*za(p5,p7)*zb(p1,p6)*zb(p1,
     &    p7)*zb(p4,p5)*tinv127*tinv345 - za(p3,p6)*za(p5,p7)*zb(p1,p6)
     &    *zb(p1,p7)*zb(p4,p6)*tinv127*tinv346 )

      z34(1,2,1,1)= + zbb1727**(-1)*s34**(-1) * (  - za(p1,p4)*za(p4,p5
     &    )*zb(p1,p2)*zb(p2,p6)*zb(p3,p4)*tinv127*tinv345 - za(p1,p5)*
     &    za(p3,p4)*zb(p1,p2)*zb(p2,p3)*zb(p3,p6)*tinv127*tinv346 - za(
     &    p1,p5)*za(p4,p5)*zb(p1,p2)*zb(p2,p6)*zb(p3,p5)*tinv127*
     &    tinv345 + za(p1,p5)*za(p4,p6)*zb(p1,p2)*zb(p2,p6)*zb(p3,p6)*
     &    tinv127*tinv346 - za(p3,p4)*za(p5,p7)*zb(p2,p3)*zb(p2,p7)*zb(
     &    p3,p6)*tinv127*tinv346 - za(p4,p5)*za(p4,p7)*zb(p2,p6)*zb(p2,
     &    p7)*zb(p3,p4)*tinv127*tinv345 - za(p4,p5)*za(p5,p7)*zb(p2,p6)
     &    *zb(p2,p7)*zb(p3,p5)*tinv127*tinv345 + za(p4,p6)*za(p5,p7)*
     &    zb(p2,p6)*zb(p2,p7)*zb(p3,p6)*tinv127*tinv346 )

      z34(1,1,2,1)= + zbb1727**(-1)*s34**(-1) * (  - za(p1,p3)*za(p3,p6
     &    )*zb(p1,p2)*zb(p2,p5)*zb(p3,p4)*tinv127*tinv346 - za(p1,p6)*
     &    za(p3,p4)*zb(p1,p2)*zb(p2,p4)*zb(p4,p5)*tinv127*tinv345 - za(
     &    p1,p6)*za(p3,p5)*zb(p1,p2)*zb(p2,p5)*zb(p4,p5)*tinv127*
     &    tinv345 + za(p1,p6)*za(p3,p6)*zb(p1,p2)*zb(p2,p5)*zb(p4,p6)*
     &    tinv127*tinv346 - za(p3,p4)*za(p6,p7)*zb(p2,p4)*zb(p2,p7)*zb(
     &    p4,p5)*tinv127*tinv345 - za(p3,p5)*za(p6,p7)*zb(p2,p5)*zb(p2,
     &    p7)*zb(p4,p5)*tinv127*tinv345 - za(p3,p6)*za(p3,p7)*zb(p2,p5)
     &    *zb(p2,p7)*zb(p3,p4)*tinv127*tinv346 + za(p3,p6)*za(p6,p7)*
     &    zb(p2,p5)*zb(p2,p7)*zb(p4,p6)*tinv127*tinv346 )

      z34(1,1,1,2)= + zaa1727**(-1)*s34**(-1) * (  - za(p1,p2)*za(p1,p3
     &    )*za(p3,p5)*zb(p2,p6)*zb(p3,p4)*tinv127*tinv345 - za(p1,p2)*
     &    za(p1,p5)*za(p3,p4)*zb(p2,p4)*zb(p4,p6)*tinv127*tinv346 + za(
     &    p1,p2)*za(p1,p5)*za(p3,p5)*zb(p2,p6)*zb(p4,p5)*tinv127*
     &    tinv345 - za(p1,p2)*za(p1,p5)*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*
     &    tinv127*tinv346 + za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p3,p4)*zb(
     &    p6,p7)*tinv127*tinv345 + za(p1,p5)*za(p1,p7)*za(p3,p4)*zb(p4,
     &    p6)*zb(p4,p7)*tinv127*tinv346 - za(p1,p5)*za(p1,p7)*za(p3,p5)
     &    *zb(p4,p5)*zb(p6,p7)*tinv127*tinv345 + za(p1,p5)*za(p1,p7)*
     &    za(p3,p6)*zb(p4,p6)*zb(p6,p7)*tinv127*tinv346 )

      z56(2,2,2,2)= + zaa1727**(-1)*s56**(-1) * ( za(p1,p2)*za(p2,p4)*
     &    za(p3,p6)*zb(p1,p3)*zb(p3,p5)*tinv127*tinv356 - za(p1,p2)*za(
     &    p2,p4)*za(p4,p6)*zb(p1,p3)*zb(p4,p5)*tinv127*tinv456 + za(p1,
     &    p2)*za(p2,p4)*za(p5,p6)*zb(p1,p5)*zb(p3,p5)*tinv127*tinv356
     &     + za(p1,p2)*za(p2,p6)*za(p4,p6)*zb(p1,p3)*zb(p5,p6)*tinv127*
     &    tinv456 + za(p2,p4)*za(p2,p7)*za(p3,p6)*zb(p3,p5)*zb(p3,p7)*
     &    tinv127*tinv356 - za(p2,p4)*za(p2,p7)*za(p4,p6)*zb(p3,p7)*zb(
     &    p4,p5)*tinv127*tinv456 + za(p2,p4)*za(p2,p7)*za(p5,p6)*zb(p3,
     &    p5)*zb(p5,p7)*tinv127*tinv356 + za(p2,p6)*za(p2,p7)*za(p4,p6)
     &    *zb(p3,p7)*zb(p5,p6)*tinv127*tinv456 )

      z56(1,1,1,1)= + zbb1727**(-1)*s56**(-1) * (  - za(p1,p3)*za(p3,p5
     &    )*zb(p1,p2)*zb(p2,p4)*zb(p3,p6)*tinv127*tinv356 + za(p1,p3)*
     &    za(p4,p5)*zb(p1,p2)*zb(p2,p4)*zb(p4,p6)*tinv127*tinv456 - za(
     &    p1,p3)*za(p5,p6)*zb(p1,p2)*zb(p2,p6)*zb(p4,p6)*tinv127*
     &    tinv456 - za(p1,p5)*za(p3,p5)*zb(p1,p2)*zb(p2,p4)*zb(p5,p6)*
     &    tinv127*tinv356 - za(p3,p5)*za(p3,p7)*zb(p2,p4)*zb(p2,p7)*zb(
     &    p3,p6)*tinv127*tinv356 - za(p3,p5)*za(p5,p7)*zb(p2,p4)*zb(p2,
     &    p7)*zb(p5,p6)*tinv127*tinv356 + za(p3,p7)*za(p4,p5)*zb(p2,p4)
     &    *zb(p2,p7)*zb(p4,p6)*tinv127*tinv456 - za(p3,p7)*za(p5,p6)*
     &    zb(p2,p6)*zb(p2,p7)*zb(p4,p6)*tinv127*tinv456 )

      z56(1,2,2,2)= + zaa1727**(-1)*s56**(-1) * ( za(p1,p2)*za(p1,p4)*
     &    za(p3,p6)*zb(p2,p3)*zb(p3,p5)*tinv127*tinv356 - za(p1,p2)*za(
     &    p1,p4)*za(p4,p6)*zb(p2,p3)*zb(p4,p5)*tinv127*tinv456 + za(p1,
     &    p2)*za(p1,p4)*za(p5,p6)*zb(p2,p5)*zb(p3,p5)*tinv127*tinv356
     &     + za(p1,p2)*za(p1,p6)*za(p4,p6)*zb(p2,p3)*zb(p5,p6)*tinv127*
     &    tinv456 - za(p1,p4)*za(p1,p7)*za(p3,p6)*zb(p3,p5)*zb(p3,p7)*
     &    tinv127*tinv356 + za(p1,p4)*za(p1,p7)*za(p4,p6)*zb(p3,p7)*zb(
     &    p4,p5)*tinv127*tinv456 - za(p1,p4)*za(p1,p7)*za(p5,p6)*zb(p3,
     &    p5)*zb(p5,p7)*tinv127*tinv356 - za(p1,p6)*za(p1,p7)*za(p4,p6)
     &    *zb(p3,p7)*zb(p5,p6)*tinv127*tinv456 )

      z56(2,1,2,2)= + zaa1727**(-1)*s56**(-1) * ( za(p1,p2)*za(p2,p3)*
     &    za(p3,p6)*zb(p1,p4)*zb(p3,p5)*tinv127*tinv356 - za(p1,p2)*za(
     &    p2,p3)*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*tinv127*tinv456 - za(p1,
     &    p2)*za(p2,p3)*za(p5,p6)*zb(p1,p5)*zb(p4,p5)*tinv127*tinv456
     &     - za(p1,p2)*za(p2,p6)*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*tinv127*
     &    tinv356 + za(p2,p3)*za(p2,p7)*za(p3,p6)*zb(p3,p5)*zb(p4,p7)*
     &    tinv127*tinv356 - za(p2,p3)*za(p2,p7)*za(p4,p6)*zb(p4,p5)*zb(
     &    p4,p7)*tinv127*tinv456 - za(p2,p3)*za(p2,p7)*za(p5,p6)*zb(p4,
     &    p5)*zb(p5,p7)*tinv127*tinv456 - za(p2,p6)*za(p2,p7)*za(p3,p6)
     &    *zb(p4,p7)*zb(p5,p6)*tinv127*tinv356 )

      z56(2,2,1,2)= + zaa1727**(-1)*s56**(-1) * ( za(p1,p2)*za(p2,p4)*
     &    za(p3,p5)*zb(p1,p3)*zb(p3,p6)*tinv127*tinv356 - za(p1,p2)*za(
     &    p2,p4)*za(p4,p5)*zb(p1,p3)*zb(p4,p6)*tinv127*tinv456 - za(p1,
     &    p2)*za(p2,p4)*za(p5,p6)*zb(p1,p6)*zb(p3,p6)*tinv127*tinv356
     &     - za(p1,p2)*za(p2,p5)*za(p4,p5)*zb(p1,p3)*zb(p5,p6)*tinv127*
     &    tinv456 + za(p2,p4)*za(p2,p7)*za(p3,p5)*zb(p3,p6)*zb(p3,p7)*
     &    tinv127*tinv356 - za(p2,p4)*za(p2,p7)*za(p4,p5)*zb(p3,p7)*zb(
     &    p4,p6)*tinv127*tinv456 - za(p2,p4)*za(p2,p7)*za(p5,p6)*zb(p3,
     &    p6)*zb(p6,p7)*tinv127*tinv356 - za(p2,p5)*za(p2,p7)*za(p4,p5)
     &    *zb(p3,p7)*zb(p5,p6)*tinv127*tinv456 )

      z56(2,2,2,1)= + zbb1727**(-1)*s56**(-1) * (  - za(p2,p4)*za(p3,p6
     &    )*zb(p1,p2)*zb(p1,p3)*zb(p3,p5)*tinv127*tinv356 + za(p2,p4)*
     &    za(p4,p6)*zb(p1,p2)*zb(p1,p3)*zb(p4,p5)*tinv127*tinv456 - za(
     &    p2,p4)*za(p5,p6)*zb(p1,p2)*zb(p1,p5)*zb(p3,p5)*tinv127*
     &    tinv356 - za(p2,p6)*za(p4,p6)*zb(p1,p2)*zb(p1,p3)*zb(p5,p6)*
     &    tinv127*tinv456 + za(p3,p6)*za(p4,p7)*zb(p1,p3)*zb(p1,p7)*zb(
     &    p3,p5)*tinv127*tinv356 - za(p4,p6)*za(p4,p7)*zb(p1,p3)*zb(p1,
     &    p7)*zb(p4,p5)*tinv127*tinv456 + za(p4,p6)*za(p6,p7)*zb(p1,p3)
     &    *zb(p1,p7)*zb(p5,p6)*tinv127*tinv456 + za(p4,p7)*za(p5,p6)*
     &    zb(p1,p5)*zb(p1,p7)*zb(p3,p5)*tinv127*tinv356 )

      z56(1,1,2,2)= + zaa1727**(-1)*s56**(-1) * ( za(p1,p2)*za(p1,p3)*
     &    za(p3,p6)*zb(p2,p4)*zb(p3,p5)*tinv127*tinv356 - za(p1,p2)*za(
     &    p1,p3)*za(p4,p6)*zb(p2,p4)*zb(p4,p5)*tinv127*tinv456 - za(p1,
     &    p2)*za(p1,p3)*za(p5,p6)*zb(p2,p5)*zb(p4,p5)*tinv127*tinv456
     &     - za(p1,p2)*za(p1,p6)*za(p3,p6)*zb(p2,p4)*zb(p5,p6)*tinv127*
     &    tinv356 - za(p1,p3)*za(p1,p7)*za(p3,p6)*zb(p3,p5)*zb(p4,p7)*
     &    tinv127*tinv356 + za(p1,p3)*za(p1,p7)*za(p4,p6)*zb(p4,p5)*zb(
     &    p4,p7)*tinv127*tinv456 + za(p1,p3)*za(p1,p7)*za(p5,p6)*zb(p4,
     &    p5)*zb(p5,p7)*tinv127*tinv456 + za(p1,p6)*za(p1,p7)*za(p3,p6)
     &    *zb(p4,p7)*zb(p5,p6)*tinv127*tinv356 )

      z56(1,2,1,2)= + zaa1727**(-1)*s56**(-1) * ( za(p1,p2)*za(p1,p4)*
     &    za(p3,p5)*zb(p2,p3)*zb(p3,p6)*tinv127*tinv356 - za(p1,p2)*za(
     &    p1,p4)*za(p4,p5)*zb(p2,p3)*zb(p4,p6)*tinv127*tinv456 - za(p1,
     &    p2)*za(p1,p4)*za(p5,p6)*zb(p2,p6)*zb(p3,p6)*tinv127*tinv356
     &     - za(p1,p2)*za(p1,p5)*za(p4,p5)*zb(p2,p3)*zb(p5,p6)*tinv127*
     &    tinv456 - za(p1,p4)*za(p1,p7)*za(p3,p5)*zb(p3,p6)*zb(p3,p7)*
     &    tinv127*tinv356 + za(p1,p4)*za(p1,p7)*za(p4,p5)*zb(p3,p7)*zb(
     &    p4,p6)*tinv127*tinv456 + za(p1,p4)*za(p1,p7)*za(p5,p6)*zb(p3,
     &    p6)*zb(p6,p7)*tinv127*tinv356 + za(p1,p5)*za(p1,p7)*za(p4,p5)
     &    *zb(p3,p7)*zb(p5,p6)*tinv127*tinv456 )

      z56(1,2,2,1)= + zbb1727**(-1)*s56**(-1) * (  - za(p1,p4)*za(p3,p6
     &    )*zb(p1,p2)*zb(p2,p3)*zb(p3,p5)*tinv127*tinv356 + za(p1,p4)*
     &    za(p4,p6)*zb(p1,p2)*zb(p2,p3)*zb(p4,p5)*tinv127*tinv456 - za(
     &    p1,p4)*za(p5,p6)*zb(p1,p2)*zb(p2,p5)*zb(p3,p5)*tinv127*
     &    tinv356 - za(p1,p6)*za(p4,p6)*zb(p1,p2)*zb(p2,p3)*zb(p5,p6)*
     &    tinv127*tinv456 - za(p3,p6)*za(p4,p7)*zb(p2,p3)*zb(p2,p7)*zb(
     &    p3,p5)*tinv127*tinv356 + za(p4,p6)*za(p4,p7)*zb(p2,p3)*zb(p2,
     &    p7)*zb(p4,p5)*tinv127*tinv456 - za(p4,p6)*za(p6,p7)*zb(p2,p3)
     &    *zb(p2,p7)*zb(p5,p6)*tinv127*tinv456 - za(p4,p7)*za(p5,p6)*
     &    zb(p2,p5)*zb(p2,p7)*zb(p3,p5)*tinv127*tinv356 )

      z56(2,1,1,2)= + zaa1727**(-1)*s56**(-1) * ( za(p1,p2)*za(p2,p3)*
     &    za(p3,p5)*zb(p1,p4)*zb(p3,p6)*tinv127*tinv356 - za(p1,p2)*za(
     &    p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)*tinv127*tinv456 + za(p1,
     &    p2)*za(p2,p3)*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*tinv127*tinv456
     &     + za(p1,p2)*za(p2,p5)*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*tinv127*
     &    tinv356 + za(p2,p3)*za(p2,p7)*za(p3,p5)*zb(p3,p6)*zb(p4,p7)*
     &    tinv127*tinv356 - za(p2,p3)*za(p2,p7)*za(p4,p5)*zb(p4,p6)*zb(
     &    p4,p7)*tinv127*tinv456 + za(p2,p3)*za(p2,p7)*za(p5,p6)*zb(p4,
     &    p6)*zb(p6,p7)*tinv127*tinv456 + za(p2,p5)*za(p2,p7)*za(p3,p5)
     &    *zb(p4,p7)*zb(p5,p6)*tinv127*tinv356 )

      z56(2,1,2,1)= + zbb1727**(-1)*s56**(-1) * (  - za(p2,p3)*za(p3,p6
     &    )*zb(p1,p2)*zb(p1,p4)*zb(p3,p5)*tinv127*tinv356 + za(p2,p3)*
     &    za(p4,p6)*zb(p1,p2)*zb(p1,p4)*zb(p4,p5)*tinv127*tinv456 + za(
     &    p2,p3)*za(p5,p6)*zb(p1,p2)*zb(p1,p5)*zb(p4,p5)*tinv127*
     &    tinv456 + za(p2,p6)*za(p3,p6)*zb(p1,p2)*zb(p1,p4)*zb(p5,p6)*
     &    tinv127*tinv356 + za(p3,p6)*za(p3,p7)*zb(p1,p4)*zb(p1,p7)*zb(
     &    p3,p5)*tinv127*tinv356 - za(p3,p6)*za(p6,p7)*zb(p1,p4)*zb(p1,
     &    p7)*zb(p5,p6)*tinv127*tinv356 - za(p3,p7)*za(p4,p6)*zb(p1,p4)
     &    *zb(p1,p7)*zb(p4,p5)*tinv127*tinv456 - za(p3,p7)*za(p5,p6)*
     &    zb(p1,p5)*zb(p1,p7)*zb(p4,p5)*tinv127*tinv456 )

      z56(2,2,1,1)= + zbb1727**(-1)*s56**(-1) * (  - za(p2,p4)*za(p3,p5
     &    )*zb(p1,p2)*zb(p1,p3)*zb(p3,p6)*tinv127*tinv356 + za(p2,p4)*
     &    za(p4,p5)*zb(p1,p2)*zb(p1,p3)*zb(p4,p6)*tinv127*tinv456 + za(
     &    p2,p4)*za(p5,p6)*zb(p1,p2)*zb(p1,p6)*zb(p3,p6)*tinv127*
     &    tinv356 + za(p2,p5)*za(p4,p5)*zb(p1,p2)*zb(p1,p3)*zb(p5,p6)*
     &    tinv127*tinv456 + za(p3,p5)*za(p4,p7)*zb(p1,p3)*zb(p1,p7)*zb(
     &    p3,p6)*tinv127*tinv356 - za(p4,p5)*za(p4,p7)*zb(p1,p3)*zb(p1,
     &    p7)*zb(p4,p6)*tinv127*tinv456 - za(p4,p5)*za(p5,p7)*zb(p1,p3)
     &    *zb(p1,p7)*zb(p5,p6)*tinv127*tinv456 - za(p4,p7)*za(p5,p6)*
     &    zb(p1,p6)*zb(p1,p7)*zb(p3,p6)*tinv127*tinv356 )

      z56(2,1,1,1)= + zbb1727**(-1)*s56**(-1) * (  - za(p2,p3)*za(p3,p5
     &    )*zb(p1,p2)*zb(p1,p4)*zb(p3,p6)*tinv127*tinv356 + za(p2,p3)*
     &    za(p4,p5)*zb(p1,p2)*zb(p1,p4)*zb(p4,p6)*tinv127*tinv456 - za(
     &    p2,p3)*za(p5,p6)*zb(p1,p2)*zb(p1,p6)*zb(p4,p6)*tinv127*
     &    tinv456 - za(p2,p5)*za(p3,p5)*zb(p1,p2)*zb(p1,p4)*zb(p5,p6)*
     &    tinv127*tinv356 + za(p3,p5)*za(p3,p7)*zb(p1,p4)*zb(p1,p7)*zb(
     &    p3,p6)*tinv127*tinv356 + za(p3,p5)*za(p5,p7)*zb(p1,p4)*zb(p1,
     &    p7)*zb(p5,p6)*tinv127*tinv356 - za(p3,p7)*za(p4,p5)*zb(p1,p4)
     &    *zb(p1,p7)*zb(p4,p6)*tinv127*tinv456 + za(p3,p7)*za(p5,p6)*
     &    zb(p1,p6)*zb(p1,p7)*zb(p4,p6)*tinv127*tinv456 )

      z56(1,2,1,1)= + zbb1727**(-1)*s56**(-1) * (  - za(p1,p4)*za(p3,p5
     &    )*zb(p1,p2)*zb(p2,p3)*zb(p3,p6)*tinv127*tinv356 + za(p1,p4)*
     &    za(p4,p5)*zb(p1,p2)*zb(p2,p3)*zb(p4,p6)*tinv127*tinv456 + za(
     &    p1,p4)*za(p5,p6)*zb(p1,p2)*zb(p2,p6)*zb(p3,p6)*tinv127*
     &    tinv356 + za(p1,p5)*za(p4,p5)*zb(p1,p2)*zb(p2,p3)*zb(p5,p6)*
     &    tinv127*tinv456 - za(p3,p5)*za(p4,p7)*zb(p2,p3)*zb(p2,p7)*zb(
     &    p3,p6)*tinv127*tinv356 + za(p4,p5)*za(p4,p7)*zb(p2,p3)*zb(p2,
     &    p7)*zb(p4,p6)*tinv127*tinv456 + za(p4,p5)*za(p5,p7)*zb(p2,p3)
     &    *zb(p2,p7)*zb(p5,p6)*tinv127*tinv456 + za(p4,p7)*za(p5,p6)*
     &    zb(p2,p6)*zb(p2,p7)*zb(p3,p6)*tinv127*tinv356 )

      z56(1,1,2,1)= + zbb1727**(-1)*s56**(-1) * (  - za(p1,p3)*za(p3,p6
     &    )*zb(p1,p2)*zb(p2,p4)*zb(p3,p5)*tinv127*tinv356 + za(p1,p3)*
     &    za(p4,p6)*zb(p1,p2)*zb(p2,p4)*zb(p4,p5)*tinv127*tinv456 + za(
     &    p1,p3)*za(p5,p6)*zb(p1,p2)*zb(p2,p5)*zb(p4,p5)*tinv127*
     &    tinv456 + za(p1,p6)*za(p3,p6)*zb(p1,p2)*zb(p2,p4)*zb(p5,p6)*
     &    tinv127*tinv356 - za(p3,p6)*za(p3,p7)*zb(p2,p4)*zb(p2,p7)*zb(
     &    p3,p5)*tinv127*tinv356 + za(p3,p6)*za(p6,p7)*zb(p2,p4)*zb(p2,
     &    p7)*zb(p5,p6)*tinv127*tinv356 + za(p3,p7)*za(p4,p6)*zb(p2,p4)
     &    *zb(p2,p7)*zb(p4,p5)*tinv127*tinv456 + za(p3,p7)*za(p5,p6)*
     &    zb(p2,p5)*zb(p2,p7)*zb(p4,p5)*tinv127*tinv456 )

      z56(1,1,1,2)= + zaa1727**(-1)*s56**(-1) * ( za(p1,p2)*za(p1,p3)*
     &    za(p3,p5)*zb(p2,p4)*zb(p3,p6)*tinv127*tinv356 - za(p1,p2)*za(
     &    p1,p3)*za(p4,p5)*zb(p2,p4)*zb(p4,p6)*tinv127*tinv456 + za(p1,
     &    p2)*za(p1,p3)*za(p5,p6)*zb(p2,p6)*zb(p4,p6)*tinv127*tinv456
     &     + za(p1,p2)*za(p1,p5)*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*tinv127*
     &    tinv356 - za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p3,p6)*zb(p4,p7)*
     &    tinv127*tinv356 + za(p1,p3)*za(p1,p7)*za(p4,p5)*zb(p4,p6)*zb(
     &    p4,p7)*tinv127*tinv456 - za(p1,p3)*za(p1,p7)*za(p5,p6)*zb(p4,
     &    p6)*zb(p6,p7)*tinv127*tinv456 - za(p1,p5)*za(p1,p7)*za(p3,p5)
     &    *zb(p4,p7)*zb(p5,p6)*tinv127*tinv356 )

      return
      end
