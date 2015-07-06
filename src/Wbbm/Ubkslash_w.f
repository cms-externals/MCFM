c-----multiplication of a barred spinor with k-slash from the right
C-----and return resultant spinor f. Weyl representation
      subroutine Ubkslash(sp,k,f) 
      implicit none
      double complex sp(4),k(4),f(4),czip,im,kslash(4,4),E,kx,ky,kz
      integer i,j
      parameter(czip=(0d0,0d0),im=(0d0,1d0))

C----create kslash
      E=k(4)
      kx=+k(3)
      ky=-k(2)
      kz=+k(1)

      kslash(1,1)=czip
      kslash(1,2)=czip
      kslash(1,3)=E+kz
      kslash(1,4)=kx-im*ky

      kslash(2,1)=czip
      kslash(2,2)=czip
      kslash(2,3)=kx+im*ky
      kslash(2,4)=E-kz

      kslash(3,1)=E-kz
      kslash(3,2)=-kx+im*ky
      kslash(3,3)=czip
      kslash(3,4)=czip

      kslash(4,1)=-kx-im*ky
      kslash(4,2)=E+kz
      kslash(4,3)=czip
      kslash(4,4)=czip

      do i=1,4
      f(i)=czip
      do j=1,4
      f(i)=f(i)+sp(j)*kslash(j,i)
      enddo
      enddo  
      return
      end



