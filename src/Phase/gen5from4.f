      subroutine gen5from4(q,z,rtalpha,phit,p,jac,*)
c----jac is the total wt of the whole business (4) and (5 from 4)
      implicit none
      include 'constants.f'
      include 'debug.f'
      integer nmin,nmax,j,iseed,i1(8),i2(8),k
      double precision p(mxpart,4),q(mxpart,4),z,rtalpha,phit,
     . wt5_4,msq(-nf:nf,-nf:nf)
      double precision sum(0:8),wtc(8),apweight(8),jac,ran0,myran
      common/apwt/apweight
      common/nmin/nmin
      common/nmax/nmax
      data iseed/1768/
      data i1/1,2,1,1,2,2,5,6/
      data i2/2,1,5,6,5,6,6,5/
      
c      call writeout(p)     

      do j=1,7
      do k=1,4
      p(j,k)=q(j,k)
      enddo
      enddo
 
      sum(nmin-1)=0d0

      do j=nmin,nmax
      apweight(j)=1d0/dfloat(nmax-nmin+1)
      sum(j)=sum(j-1)+apweight(j)
      if (debug) then
      write(6,*) 'j',j 
      write(6,*) 'apweight(j)',apweight(j) 
      write(6,*) 'sum(j)',sum(j) 
      endif
      enddo
  
      myran=ran0(iseed)

      do j=nmin,nmax
      if ((myran .gt. sum(j-1)) .and. (myran .lt. sum(j))) then
c---genrad is a switchyard routine routing to genrii,genrif,genrff
      call genrad(p,i1(j),i2(j),7,z,rtalpha,phit,wt5_4,*999)
c---although genrad returns wt5_4 we shall not use it 
c---in this step we have generated the new p's (only one set)
c---only one option is pursued in this do-loop 
      endif      
      enddo

c---Sum over channels
c---Initialize jac
      jac=0d0
      do j=nmin,nmax
         if ((j .eq. 1) .or. (j .eq. 2)) 
     .      call genii(j,p,wtc(j),msq)
         if ((j .eq. 3) .or. (j .eq. 4)) 
     .      call genff(j-2,p,wtc(j),msq)
         if ((j .eq. 5) .or. (j .eq. 6) .or. (j .eq. 7) .or. (j .eq. 8))
     .      call genif(j-4,p,wtc(j),msq)
        jac=jac+apweight(j)/wtc(j)
      enddo
      jac=1d0/jac


      return 

 999  jac=0d0
      return 1
      end

