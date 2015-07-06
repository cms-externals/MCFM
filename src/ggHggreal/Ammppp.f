*********************************************************
*  Code by Giulia Zanderighi
*  Based on eq.(3.27) of hep-th/0411092 which gives 1+ 2+ 3- 4- 5-
*  Taking the complex conjugate gives then 1- 2- 3+ 4+ 5+
*  Explicit result obtained with the Maple file NMHV 
      DOUBLE COMPLEX FUNCTION  Ammppp(I1,I2,I3,I4,I5)                         
* ---------------------------------------------------------------------
*                            1- 2- 3+ 4+ 5+                            
* ---------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER I1,I2,I3,I4,I5
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex zab,zab2,zab3,za_4,zab2_3,zab3_3
      double precision s3,s4
      integer j1,j2,j3,j4,j5

c      double complex Ammppp2 

C---statement functions
C       zab(j1,j2,j3)        = za(j1,j2)*zb(j2,j3)
C       zab2(j1,j2,j3,j4)    = zab(j1,j2,j4)+zab(j1,j3,j4)
C       zab3(j1,j2,j3,j4,j5) = zab(j1,j2,j5)+zab(j1,j3,j5)+zab(j1,j4,j5)
C       zab4(j1,j2,j3,j4,j5,j6) = zab(j1,j2,j6)+zab(j1,j3,j6)
C     .      +zab(j1,j4,j6)+zab(j1,j5,j6)

       zab(j1,j2,j3)           = za(j1,j2)*zb(j2,j3)
       zab2(j1,j2,j3,j4)       = zab(j1,j2,j4)+zab(j1,j3,j4)
       zab3(j1,j2,j3,j4,j5)    = zab2(j1,j2,j3,j5)+zab(j1,j4,j5)

       za_4(j1,j2)               = za(j1,j2)**4
       zab2_3(j1,j2,j3,j4)       = zab2(j1,j2,j3,j4)**3
       zab3_3(j1,j2,j3,j4,j5)    = zab3(j1,j2,j3,j4,j5)**3

       S3(i1,i2,i3) = S(i1,i2)+S(i1,i3)+S(i2,i3)

       S4(i1,i2,i3,i4) = S(i1,i2)+S(i1,i3)+S(i1,i4)+S(i2,i3)+
     .      S(i2,i4)+S(i3,i4)


C---   Fully simplified with specific choice of q=i2 
       Ammppp = 
     . (
     . za_4(i4,i5)*zab2_3(i3,i4,i5,i2)*za(i3,i4)*za(i5,i1)/        
     .(zab2(i1,i4,i5,i2)*zab(i4,i5,i2)*zab(i5,i4,i2)*S(i4,i5))
     .+za_4(i4,i5)*zab3_3(i3,i1,i4,i5,i2)*za(i3,i4)*za(i1,i2)/
     .(zab3(i2,i1,i4,i5,i2)*zab2(i4,i1,i5,i2)*zab2(i1,i4,i5,i2)*
     .(S3(i1,i4,i5)))
     .+za_4(i4,i5)*zab3(i3,i1,i4,i5,i2)**2*za(i3,i4)*za(i2,i3)/
     .(zab2(i4,i1,i5,i2)*zab3(i2,i1,i4,i5,i2)*
     .(S4(i1,i2,i4,i5)))
     .+za_4(i5,i3)*zab3(i4,i1,i3,i5,i2)**2*za(i4,i5)*za(i3,i4)/
     .(zab2(i5,i1,i3,i2)*zab2(i3,i1,i5,i2)*
     .(S4(i1,i2,i3,i5)))
     .+za_4(i3,i4)*zab3(i5,i1,i3,i4,i2)**2*za(i5,i1)*za(i4,i5)/
     .(zab2(i1,i3,i4,i2)*zab2(i4,i1,i3,i2)
     .*(S4(i1,i2,i3,i4)))
     .+za_4(i3,i4)*zab2_3(i5,i3,i4,i2)*za(i1,i2)*za(i4,i5)/
     .(zab2(i1,i3,i4,i2)*zab2(i2,i3,i4,i2)*zab(i4,i3,i2)*
     .(S3(i2,i3,i4)))
     .+za_4(i3,i4)*zab2_3(i5,i3,i4,i2)*za(i2,i3)*za(i4,i5)/
     .(zab2(i2,i3,i4,i2)*zab(i3,i4,i2)*zab(i4,i3,i2)*S(i3,i4))
     .+za_4(i4,i5)*zab(i3,i1,i2)**3*za(i3,i4)*za(i5,i1)
     ./(zab(i1,i3,i2)*zab2(i4,i1,i3,i2)*zab2(i5,i1,i3,i2)
     .*(S3(i1,i2,i3)))
     .+za_4(i3,i4)*zab2_3(i5,i1,i5,i2)*za(i1,i2)*za(i4,i5)/
     .(zab(i1,i5,i2)*zab2(i2,i1,i5,i2)*zab2(i4,i1,i5,i2)*S(i1,i5))
     .+za_4(i3,i4)*zab(i5,i1,i2)**3*za(i2,i3)*za(i4,i5)
     ./(zab2(i2,i1,i5,i2)*zab2(i3,i1,i5,i2)*zab2(i4,i1,i5,i2)
     .*(S3(i1,i2,i5)))
     .)/(za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i4,i5)*za(i5,i1))


!!C---   Fully simplified 
!       Ammppp2 = 
!     . (
!     . za_4(i4,i5)*zab2_3(i3,i4,i5,q)*za(i3,i4)*za(i5,i1)/        
!     .(zab2(i1,i4,i5,q)*zab(i4,i5,q)*zab(i5,i4,q)*S(i4,i5))
!     .+za_4(i4,i5)*zab3_3(i3,i1,i4,i5,q)*za(i3,i4)*za(i1,i2)/
!     .(zab3(i2,i1,i4,i5,q)*zab2(i4,i1,i5,q)*zab2(i1,i4,i5,q)*
!     .(S3(i1,i4,i5)))
!     .+za_4(i4,i5)*zab4_2(i3,i1,i2,i4,i5,q)*za(i3,i4)*za(i2,i3)/
!     .(zab3(i4,i1,i2,i5,q)*zab3(i2,i1,i4,i5,q)*
!     .(S4(i1,i2,i4,i5)))
!     .+za_4(i5,i3)*zab4_2(i4,i1,i2,i3,i5,q)*za(i4,i5)*za(i3,i4)/
!     .(zab3(i5,i1,i2,i3,q)*zab3(i3,i1,i2,i5,q)*
!     .(S4(i1,i2,i3,i5)))
!     .+za_4(i3,i4)*zab4_2(i5,i1,i2,i3,i4,q)*za(i5,i1)*za(i4,i5)/
!     .(zab3(i1,i2,i3,i4,q)*zab3(i4,i1,i2,i3,q)
!     .*(S4(i1,i2,i3,i4)))
!     .+za_4(i3,i4)*zab3_3(i5,i2,i3,i4,q)*za(i1,i2)*za(i4,i5)/
!     .(zab3(i1,i2,i3,i4,q)*zab2(i2,i3,i4,q)*zab2(i4,i2,i3,q)*
!     .(S3(i2,i3,i4)))
!     .+za_4(i3,i4)*zab2_3(i5,i3,i4,q)*za(i2,i3)*za(i4,i5)/
!     .(zab2(i2,i3,i4,q)*zab(i3,i4,q)*zab(i4,i3,q)*S(i3,i4))
!     .+za_4(i4,i5)*zab2_3(i3,i1,i2,q)*za(i3,i4)*za(i5,i1)
!     ./(zab2(i1,i2,i3,q)*zab3(i4,i1,i2,i3,q)*zab3(i5,i1,i2,i3,q)
!     .*(S3(i1,i2,i3)))
!     .+za_4(i4,i5)*zab2_3(i3,i2,i3,q)*za(i3,i4)*za(i1,i2)/
!     .(zab(i2,i3,q)*zab2(i4,i2,i3,q)*zab2(i1,i2,i3,q)*S(i2,i3))
!     .+za_4(i3,i4)*zab2_3(i5,i1,i5,q)*za(i1,i2)*za(i4,i5)/
!     .(zab(i1,i5,q)*zab2(i2,i1,i5,q)*zab2(i4,i1,i5,q)*S(i1,i5))
!     .+za_4(i3,i4)*zab2_3(i5,i1,i2,q)*za(i2,i3)*za(i4,i5)
!     ./(zab2(i2,i1,i5,q)*zab3(i3,i1,i2,i5,q)*zab3(i4,i1,i2,i5,q)
!     .*(S3(i1,i2,i5)))
!     .)/(za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i4,i5)*za(i5,i1))
!
!
!       if (abs(Ammppp-Ammppp2) > 0.0001d0 ) then 
!          write(*,*) Ammppp,Ammppp2
!          stop 
!       endif
!
      Ammppp=Ammppp-zb(i1,i2)**4/(zb(i1,i2)*zb(i2,i3)*zb(i3,i4)*
     .zb(i4,i5)*zb(i5,i1))

      Ammppp= dconjg(Ammppp)   
      !write(*,*) 'Ammppp',Ammppp
      !pause 

      END


! ORIG 
!       Ammppp = 
!     . (
!     . za(i4,i5)**4*zab2(i3,i4,i5,q)**3*za(i3,i4)*za(i5,i1)/        
!     .(zab2(i1,i4,i5,q)*zab2(i4,i4,i5,q)*zab2(i5,i4,i5,q)*S(i4,i5))
!     .+za(i4,i5)**4*zab3(i3,i1,i4,i5,q)**3*za(i3,i4)*za(i1,i2)/
!     .(zab3(i2,i1,i4,i5,q)*zab3(i4,i1,i4,i5,q)*zab3(i1,i1,i4,i5,q)*
!     .(S(i1,i4)+S(i1,i5)+S(i4,i5)))
!     .+za(i4,i5)**4*zab4(i3,i1,i2,i4,i5,q)**2*za(i3,i4)*za(i2,i3)/
!     .(zab4(i4,i1,i2,i4,i5,q)*zab4(i2,i1,i2,i4,i5,q)*
!     .(S(i1,i2)+S(i1,i4)+S(i1,i5)+S(i2,i4)+S(i2,i5)+S(i4,i5)))
!     .+za(i5,i3)**4*zab4(i4,i1,i2,i3,i5,q)**2*za(i4,i5)*za(i3,i4)/
!     .(zab4(i5,i1,i2,i3,i5,q)*zab4(i3,i1,i2,i3,i5,q)*
!     .(S(i1,i2)+S(i1,i3)+S(i1,i5)+S(i2,i3)+S(i2,i5)+S(i3,i5)))
!     .+za(i3,i4)**4*zab4(i5,i1,i2,i3,i4,q)**2*za(i5,i1)*za(i4,i5)/
!     .(zab4(i1,i1,i2,i3,i4,q)*zab4(i4,i1,i2,i3,i4,q)
!     .*(S(i1,i2)+S(i1,i3)+S(i1,i4)+S(i2,i3)+S(i2,i4)+S(i3,i4)))
!     .+za(i3,i4)**4*zab3(i5,i2,i3,i4,q)**3*za(i1,i2)*za(i4,i5)/
!     .(zab3(i1,i2,i3,i4,q)*zab3(i2,i2,i3,i4,q)*zab3(i4,i2,i3,i4,q)*
!     .(S(i2,i3)+S(i2,i4)+S(i3,i4)))
!     .+za(i3,i4)**4*zab2(i5,i3,i4,q)**3*za(i2,i3)*za(i4,i5)/
!     .(zab2(i2,i3,i4,q)*zab2(i3,i3,i4,q)*zab2(i4,i3,i4,q)*S(i3,i4))
!     .+za(i4,i5)**4*zab3(i3,i1,i2,i3,q)**3*za(i3,i4)*za(i5,i1)
!     ./(zab3(i1,i1,i2,i3,q)*zab3(i4,i1,i2,i3,q)*zab3(i5,i1,i2,i3,q)
!     .*(S(i1,i2)+S(i1,i3)+S(i2,i3)))
!     .+za(i4,i5)**4*zab2(i3,i2,i3,q)**3*za(i3,i4)*za(i1,i2)/
!     .(zab2(i2,i2,i3,q)*zab2(i4,i2,i3,q)*zab2(i1,i2,i3,q)*S(i2,i3))
!     .+za(i3,i4)**4*zab2(i5,i1,i5,q)**3*za(i1,i2)*za(i4,i5)/
!     .(zab2(i1,i1,i5,q)*zab2(i2,i1,i5,q)*zab2(i4,i1,i5,q)*S(i1,i5))
!     .+za(i3,i4)**4*zab3(i5,i1,i2,i5,q)**3*za(i2,i3)*za(i4,i5)
!     ./(zab3(i2,i1,i2,i5,q)*zab3(i3,i1,i2,i5,q)*zab3(i4,i1,i2,i5,q)
!     .*(S(i1,i2)+S(i1,i5)+S(i2,i5)))
!     .)/(za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i4,i5)*za(i5,i1))

! trivial semplifications 
!       Ammppp = 
!     . (
!     . za(i4,i5)**4*zab2(i3,i4,i5,q)**3*za(i3,i4)*za(i5,i1)/        
!     .(zab2(i1,i4,i5,q)*zab(i4,i5,q)*zab(i5,i4,q)*S(i4,i5))
!     .+za(i4,i5)**4*zab3(i3,i1,i4,i5,q)**3*za(i3,i4)*za(i1,i2)/
!     .(zab3(i2,i1,i4,i5,q)*zab2(i4,i1,i5,q)*zab2(i1,i4,i5,q)*
!     .(S(i1,i4)+S(i1,i5)+S(i4,i5)))
!     .+za(i4,i5)**4*zab4(i3,i1,i2,i4,i5,q)**2*za(i3,i4)*za(i2,i3)/
!     .(zab3(i4,i1,i2,i5,q)*zab3(i2,i1,i4,i5,q)*
!     .(S(i1,i2)+S(i1,i4)+S(i1,i5)+S(i2,i4)+S(i2,i5)+S(i4,i5)))
!     .+za(i5,i3)**4*zab4(i4,i1,i2,i3,i5,q)**2*za(i4,i5)*za(i3,i4)/
!     .(zab3(i5,i1,i2,i3,q)*zab3(i3,i1,i2,i5,q)*
!     .(S(i1,i2)+S(i1,i3)+S(i1,i5)+S(i2,i3)+S(i2,i5)+S(i3,i5)))
!     .+za(i3,i4)**4*zab4(i5,i1,i2,i3,i4,q)**2*za(i5,i1)*za(i4,i5)/
!     .(zab3(i1,i2,i3,i4,q)*zab3(i4,i1,i2,i3,q)
!     .*(S(i1,i2)+S(i1,i3)+S(i1,i4)+S(i2,i3)+S(i2,i4)+S(i3,i4)))
!     .+za(i3,i4)**4*zab3(i5,i2,i3,i4,q)**3*za(i1,i2)*za(i4,i5)/
!     .(zab3(i1,i2,i3,i4,q)*zab2(i2,i3,i4,q)*zab2(i4,i2,i3,q)*
!     .(S(i2,i3)+S(i2,i4)+S(i3,i4)))
!     .+za(i3,i4)**4*zab2(i5,i3,i4,q)**3*za(i2,i3)*za(i4,i5)/
!     .(zab2(i2,i3,i4,q)*zab(i3,i4,q)*zab(i4,i3,q)*S(i3,i4))
!     .+za(i4,i5)**4*zab2(i3,i1,i2,q)**3*za(i3,i4)*za(i5,i1)
!     ./(zab2(i1,i2,i3,q)*zab3(i4,i1,i2,i3,q)*zab3(i5,i1,i2,i3,q)
!     .*(S(i1,i2)+S(i1,i3)+S(i2,i3)))
!     .+za(i4,i5)**4*zab(i3,i2,q)**3*za(i3,i4)*za(i1,i2)/
!     .(zab(i2,i3,q)*zab2(i4,i2,i3,q)*zab2(i1,i2,i3,q)*S(i2,i3))
!     .+za(i3,i4)**4*zab(i5,i1,q)**3*za(i1,i2)*za(i4,i5)/
!     .(zab(i1,i5,q)*zab2(i2,i1,i5,q)*zab2(i4,i1,i5,q)*S(i1,i5))
!     .+za(i3,i4)**4*zab2(i5,i1,i2,q)**3*za(i2,i3)*za(i4,i5)
!     ./(zab2(i2,i1,i5,q)*zab3(i3,i1,i2,i5,q)*zab3(i4,i1,i2,i5,q)
!     .*(S(i1,i2)+S(i1,i5)+S(i2,i5)))
!     .)/(za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i4,i5)*za(i5,i1))



!C---   further Simplified 
!       Ammppp2 = 
!     . (
!     . za_4(i4,i5)*zab2_3(i3,i4,i5,q)*za(i3,i4)*za(i5,i1)/        
!     .(zab2(i1,i4,i5,q)*zab(i4,i5,q)*zab(i5,i4,q)*S(i4,i5))
!     .+za_4(i4,i5)*zab3_3(i3,i1,i4,i5,q)*za(i3,i4)*za(i1,i2)/
!     .(zab3(i2,i1,i4,i5,q)*zab2(i4,i1,i5,q)*zab2(i1,i4,i5,q)*
!     .(S(i1,i4)+S(i1,i5)+S(i4,i5)))
!     .+za_4(i4,i5)*zab4_2(i3,i1,i2,i4,i5,q)*za(i3,i4)*za(i2,i3)/
!     .(zab3(i4,i1,i2,i5,q)*zab3(i2,i1,i4,i5,q)*
!     .(S(i1,i2)+S(i1,i4)+S(i1,i5)+S(i2,i4)+S(i2,i5)+S(i4,i5)))
!     .+za_4(i5,i3)*zab4_2(i4,i1,i2,i3,i5,q)*za(i4,i5)*za(i3,i4)/
!     .(zab3(i5,i1,i2,i3,q)*zab3(i3,i1,i2,i5,q)*
!     .(S(i1,i2)+S(i1,i3)+S(i1,i5)+S(i2,i3)+S(i2,i5)+S(i3,i5)))
!     .+za_4(i3,i4)*zab4_2(i5,i1,i2,i3,i4,q)*za(i5,i1)*za(i4,i5)/
!     .(zab3(i1,i2,i3,i4,q)*zab3(i4,i1,i2,i3,q)
!     .*(S(i1,i2)+S(i1,i3)+S(i1,i4)+S(i2,i3)+S(i2,i4)+S(i3,i4)))
!     .+za_4(i3,i4)*zab3_3(i5,i2,i3,i4,q)*za(i1,i2)*za(i4,i5)/
!     .(zab3(i1,i2,i3,i4,q)*zab2(i2,i3,i4,q)*zab2(i4,i2,i3,q)*
!     .(S(i2,i3)+S(i2,i4)+S(i3,i4)))
!     .+za_4(i3,i4)*zab2_3(i5,i3,i4,q)*za(i2,i3)*za(i4,i5)/
!     .(zab2(i2,i3,i4,q)*zab(i3,i4,q)*zab(i4,i3,q)*S(i3,i4))
!     .+za_4(i4,i5)*zab2_3(i3,i1,i2,q)*za(i3,i4)*za(i5,i1)
!     ./(zab2(i1,i2,i3,q)*zab3(i4,i1,i2,i3,q)*zab3(i5,i1,i2,i3,q)
!     .*(S(i1,i2)+S(i1,i3)+S(i2,i3)))
!     .+za_4(i4,i5)*zab2_3(i3,i2,i3,q)*za(i3,i4)*za(i1,i2)/
!     .(zab(i2,i3,q)*zab2(i4,i2,i3,q)*zab2(i1,i2,i3,q)*S(i2,i3))
!     .+za_4(i3,i4)*zab2_3(i5,i1,i5,q)*za(i1,i2)*za(i4,i5)/
!     .(zab(i1,i5,q)*zab2(i2,i1,i5,q)*zab2(i4,i1,i5,q)*S(i1,i5))
!     .+za_4(i3,i4)*zab2_3(i5,i1,i2,q)*za(i2,i3)*za(i4,i5)
!     ./(zab2(i2,i1,i5,q)*zab3(i3,i1,i2,i5,q)*zab3(i4,i1,i2,i5,q)
!     .*(S(i1,i2)+S(i1,i5)+S(i2,i5)))
!     .)/(za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i4,i5)*za(i5,i1))
