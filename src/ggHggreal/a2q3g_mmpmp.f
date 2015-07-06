*********************************************************
* Based on eq.(3.15) of hep-th/0412275 which gives q1- 2g- g3+ g4- qb5+
* Explicit result obtained with the Maple file NMHV_qq (m2=2,m3=4)
      DOUBLE COMPLEX FUNCTION  A2q3g_mmpmp(I1,I2,I3,I4,I5,za,zb)
* ---------------------------------------------------------------------
*                            q1- 2g- 3g+ 4g- qb5+
* ---------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER I1,I2,I3,I4,I5
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
C      include 'zprods_com.f'
      double complex zab,zab2,zab3
      double precision s3,s4
      integer j1,j2,j3,j4,j5


C---statement functions
       zab(j1,j2,j3)           = za(j1,j2)*zb(j2,j3)
       zab2(j1,j2,j3,j4)       = zab(j1,j2,j4)+zab(j1,j3,j4)
       zab3(j1,j2,j3,j4,j5)    = zab2(j1,j2,j3,j5)+zab(j1,j4,j5)

       S3(j1,j2,j3) = S(j1,j2)+S(j1,j3)+S(j2,j3)
       S4(j1,j2,j3,j4) = S(j1,j2)+S(j1,j3)+S(j1,j4)+S(j2,j3)+
     .      S(j2,j4)+S(j3,j4)

       
C---   Simplified expression with specific choice of q=i5 
       A2q3g_mmpmp = 
     . (za(i2,i1)**3*zab2(i4,i1,i2,i5)**3*za(i4,i5)/zab2(i5,i1,i2,i5)
     . /zab2(i3,i1,i2,i5)/zab(i1,i2,i5)/S(i1,i2)*za(i5,i1)*za(i2,i3)+
     . za(i2,i1)**3*zab2(i2,i1,i3,i5)*zab3(i4,i1,i2,i3,i5)**2*za(i4,i5)
     . /zab3(i5,i1,i2,i3,i5)/zab2(i1,i2,i3,i5)/zab2(i3,i1,i2,i5)/
     . S3(i1,i2,i3)*za(i5,i1)*za(i3,i4))/za(i1,i2)/
     .  za(i2,i3)/za(i3,i4)/za(i4,i5)/za(i5,i1)+(za(i2,i1)**3*
     .  zab(i4,i3,i5)**3*za(i4,i5)/zab2(i1,i3,i4,i5)/
     .  zab(i3,i4,i5)/zab2(i5,i3,i4,i5)/S3(i3,i4,i5)
     .  *za(i2,i3)*za(i5,i1))/za(i1,i2)/
     .  za(i2,i3)/za(i3,i4)/za(i4,i5)/za(i5,i1)
     .    +za(i2,i1)**3*
     .  za(i2,i5)*zab(i4,i3,i5)**3/zab2(i2,i3,i4,i5)/zab2(i5,i3,i4,i5)
     .  /zab(i3,i4,i5)/S(i3,i4)/za(i1,i2)/za(i3,i4)/za(i5,i1)+
     .  (za(i2,i1)**3*za(i2,i5)*zab2(i4,i1,i2,i5)**3/zab2(i3,i1,i2
     .  ,i5)/zab2(i5,i1,i2,i5)/zab(i2,i1,i5)/S3(i1,i2,i5)
     .  *za(i4,i5)*za(i2,i3)+za(i2,i1)**3*za(i2,i5)*
     .  zab3(i4,i1,i2,i3,i5)**2/zab3(i5,i1,i2,i3,i5)/
     .  zab2(i3,i1,i2,i5)/S4(i1,i2,i3,i5)
     .  *za(i4,i5)*za(i3,i4))/za(i1,i2)/za(i2,i3)/
     .  za(i3,i4)/za(i4,i5)/za(i5,i1)+zab3(i1,i2,i3,i4,i5)**2*
     .  za(i2,i4)**4/zab2(i2,i3,i4,i5)/zab2(i4,i2,i3,i5)/
     .  S3(i2,i3,i4)/za(i2,i3)/za(i3,i4)/za(i5,i1)
     .  +za(i4,i1)**3*
     .  zab2(i2,i2,i3,i5)**3/zab2(i1,i2,i3,i5)/zab2(i4,i2,i3,i5)/
     .  zab(i3,i2,i5)/S(i2,i3)/za(i2,i3)/za(i5,i1)+(za(i4,i1)**3*
     .  za(i4,i5)*zab3(i2,i1,i3,i4,i5)**2/zab2(i3,i1,i4,i5)/
     .  zab2(i1,i3,i4,i5)/S4(i1,i3,i4,i5)
     .  *za(i2,i3)*za(i1,i2)+za(i4,i1)**3*za(i4,i5)*
     .  zab2(i2,i1,i4,i5)**3/zab2(i3,i1,i4,i5)/zab(i4,i1,i5)/
     .  zab(i1,i4,i5)/S3(i1,i4,i5)*za(i3,i4)*
     .  za(i1,i2))/za(i1,i2)/za(i2,i3)/za(i3,i4)/za(i4,i5)/za(i5,i1)

       ! NEED TO UNDERSTAND THE SIGN! 
       A2q3g_mmpmp = - A2q3g_mmpmp - zb(i3,i1)*zb(i3,i5)**3/
     .  (zb(i1,i2)*zb(i2,i3)*zb(i3,i4)*zb(i4,i5)*zb(i5,i1))

      
        END


CC---   Not simplified 
C       A2q3g_mmpmp2 = 
C     . (za(i2,i1)**3*zab2(i4,i1,i2,q)**3*za(i4,i5)/zab2(i5,i1,i2,q)
C     . /zab2(i3,i1,i2,q)/zab2(i1,i1,i2,q)/S(i1,i2)*za(i5,i1)*za(i2,i3)+
C     . za(i2,i1)**3*zab3(i2,i1,i2,i3,q)*zab3(i4,i1,i2,i3,q)**2*za(i4,i5)
C     . /zab3(i5,i1,i2,i3,q)/zab3(i1,i1,i2,i3,q)/zab3(i3,i1,i2,i3,q)/
C     . (S(i1,i2)+S(i1,i3)+S(i2,i3))*za(i5,i1)*za(i3,i4))/za(i1,i2)/
C     .  za(i2,i3)/za(i3,i4)/za(i4,i5)/za(i5,i1)+(za(i2,i1)**3*
C     .  zab3(i4,i3,i4,i5,q)**3*za(i4,i5)/zab3(i1,i3,i4,i5,q)/
C     .  zab3(i3,i3,i4,i5,q)/zab3(i5,i3,i4,i5,q)/(S(i3,i4)+S(i3,i5)+
C     .  S(i4,i5))*za(i2,i3)*za(i5,i1)+za(i2,i1)**3*zab2(i2,i4,i5,q)*
C     .  zab2(i4,i4,i5,q)**2*za(i4,i5)/zab2(i3,i4,i5,q)/zab2(i1,i4,i5,q)
C     .  /zab2(i5,i4,i5,q)/S(i4,i5)*za(i3,i4)*za(i5,i1))/za(i1,i2)/
C     .  za(i2,i3)/za(i3,i4)/za(i4,i5)/za(i5,i1)+za(i2,i1)**3*
C     .  za(i2,i5)*zab2(i4,i3,i4,q)**3/zab2(i2,i3,i4,q)/zab2(i5,i3,i4,q)
C     .  /zab2(i3,i3,i4,q)/S(i3,i4)/za(i1,i2)/za(i3,i4)/za(i5,i1)+
C     .  (za(i2,i1)**3*za(i2,i5)*zab3(i4,i1,i2,i5,q)**3/zab3(i3,i1,i2
C     .  ,i5,q)/zab3(i5,i1,i2,i5,q)/zab3(i2,i1,i2,i5,q)/(S(i1,i2)+
C     .  S(i1,i5)+S(i2,i5))*za(i4,i5)*za(i2,i3)+za(i2,i1)**3*za(i2,i5)*
C     .  zab4(i4,i1,i2,i3,i5,q)**2/zab4(i5,i1,i2,i3,i5,q)/
C     .  zab4(i3,i1,i2,i3,i5,q)/(S(i1,i2)+S(i1,i3)+S(i1,i5)+S(i2,i3)+
C     .  S(i2,i5)+S(i3,i5))*za(i4,i5)*za(i3,i4))/za(i1,i2)/za(i2,i3)/
C     .  za(i3,i4)/za(i4,i5)/za(i5,i1)+zab3(i1,i2,i3,i4,q)**2*
C     .  za(i2,i4)**4/zab3(i2,i2,i3,i4,q)/zab3(i4,i2,i3,i4,q)/
C     .  (S(i2,i3)+S(i2,i4)+S(i3,i4))/za(i2,i3)/za(i3,i4)/za(i5,i1)+
C     .  zab2(i1,i1,i5,q)**2*za(i2,i4)**4/zab2(i4,i1,i5,q)/zab2(i2,i1,i5,
C     .  q)/S(i1,i5)/za(i2,i3)/za(i3,i4)/za(i5,i1)+za(i4,i1)**3*
C     .  zab2(i2,i2,i3,q)**3/zab2(i1,i2,i3,q)/zab2(i4,i2,i3,q)/
C     .  zab2(i3,i2,i3,q)/S(i2,i3)/za(i2,i3)/za(i5,i1)+(za(i4,i1)**3*
C     .  za(i4,i5)*zab4(i2,i1,i3,i4,i5,q)**2/zab4(i3,i1,i3,i4,i5,q)/
C     .  zab4(i1,i1,i3,i4,i5,q)/(S(i1,i3)+S(i1,i4)+S(i1,i5)+S(i3,i4)+
C     .  S(i3,i5)+S(i4,i5))*za(i2,i3)*za(i1,i2)+za(i4,i1)**3*za(i4,i5)*
C     .  zab3(i2,i1,i4,i5,q)**3/zab3(i3,i1,i4,i5,q)/zab3(i4,i1,i4,i5,q)/
C     .  zab3(i1,i1,i4,i5,q)/(S(i1,i4)+S(i1,i5)+S(i4,i5))*za(i3,i4)*
C     .  za(i1,i2))/za(i1,i2)/za(i2,i3)/za(i3,i4)/za(i4,i5)/za(i5,i1)
C
C       if (abs(A2q3g_mmpmp-A2q3g_mmpmp2) > 0.0001d0) then 
C          write(*,*) 'A2q3g_mmpmp',A2q3g_mmpmp,A2q3g_mmpmp2
C          stop
C       endif

