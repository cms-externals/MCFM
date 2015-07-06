      subroutine zzamps(j1,j2,j3,j4,j5,j6,j7,za,zb,f)
c--- This is the old code for the amplitudes for ZZ+gluon production
c--- (no singly-resonant diagrams are included)
c---labels on array f are polarizations of partons:
c--- 1=quark, 2=lepton1, 3=lepton2, 4=gluon

      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      double complex A7trees
      double complex f(2,2,2,2)

      f(1,1,1,2)=+A7trees(j1,j2,j3,j4,j5,j6,j7,za,zb)  
      f(1,1,2,2)=+A7trees(j1,j2,j4,j3,j5,j6,j7,za,zb)  
      f(1,2,1,2)=+A7trees(j1,j2,j3,j4,j6,j5,j7,za,zb)  
      f(1,2,2,2)=+A7trees(j1,j2,j4,j3,j6,j5,j7,za,zb)  

      f(1,1,1,1)=-A7trees(j2,j1,j5,j6,j3,j4,j7,zb,za)  
      f(1,1,2,1)=-A7trees(j2,j1,j5,j6,j4,j3,j7,zb,za)  
      f(1,2,1,1)=-A7trees(j2,j1,j6,j5,j3,j4,j7,zb,za)  
      f(1,2,2,1)=-A7trees(j2,j1,j6,j5,j4,j3,j7,zb,za)  

      f(2,2,2,1)=-A7trees(j1,j2,j3,j4,j5,j6,j7,zb,za)  
      f(2,2,1,1)=-A7trees(j1,j2,j4,j3,j5,j6,j7,zb,za)  
      f(2,1,2,1)=-A7trees(j1,j2,j3,j4,j6,j5,j7,zb,za)  
      f(2,1,1,1)=-A7trees(j1,j2,j4,j3,j6,j5,j7,zb,za)  

      f(2,2,2,2)=+A7trees(j2,j1,j5,j6,j3,j4,j7,za,zb)  
      f(2,2,1,2)=+A7trees(j2,j1,j5,j6,j4,j3,j7,za,zb)  
      f(2,1,2,2)=+A7trees(j2,j1,j6,j5,j3,j4,j7,za,zb)  
      f(2,1,1,2)=+A7trees(j2,j1,j6,j5,j4,j3,j7,za,zb)  


c      f(mplus,minus,minus,mplus)=+A7trees(j1,j2,j4,j3,j6,j5,j7,zb,za)  
c      f(mplus,minus,mplus,mplus)=+A7trees(j1,j2,j3,j4,j6,j5,j7,zb,za)  
c      f(mplus,mplus,minus,mplus)=+A7trees(j1,j2,j4,j3,j5,j6,j7,zb,za)  
c      f(mplus,mplus,mplus,mplus)=+A7trees(j1,j2,j3,j4,j5,j6,j7,zb,za)  

c      f(mplus,minus,minus,minus)=-A7trees(j2,j1,j6,j5,j4,j3,j7,za,zb)  
c      f(mplus,minus,mplus,minus)=-A7trees(j2,j1,j6,j5,j3,j4,j7,za,zb)  
c      f(mplus,mplus,minus,minus)=-A7trees(j2,j1,j5,j6,j4,j3,j7,za,zb)  
c      f(mplus,mplus,mplus,minus)=-A7trees(j2,j1,j5,j6,j3,j4,j7,za,zb)  

      return
      end
