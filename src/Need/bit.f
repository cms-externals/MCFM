***************************** Quark-Quark *****************************
*     (case when the emitter is massive and spectator is massless)    *
***********************************************************************
      double precision function ff_mqq0(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,ddilog
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for a
c--- final-final quark-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C --MSbar
c g zipffqq0=colfac*del(omx)*(
c epinv^2/2+epinv*(1-L/2+log(mbarsq)/2-log(1-mbarsq))-L+3+L**2/4-pisq*5/12
c  -dilog(1-mbarsq)-log(mbarsq)*(log(mbarsq)/2-1+L)/2
c  -log(1-mbarsq)*(log(mbarsq,1-mbarsq)-L+2)
c  -mbarsq/(1-mbarsq)*log(mbarsq)
c  )-ffqq0;

c--- Note: this uses the symmetric form of the eikonal integrals, so
c---       that it can match with ff_mgg

      if (aff .lt. 1d0) then
        write(6,*) 'Integrated dipole routine ff_mqq0 does not yet'
	write(6,*) ' support values of alpha < 1.'
	stop
      endif

      ff_mqq0=0d0 
      if (vorz .eq. 1) then
        if (mbar .gt. 1d0) then
            write(6,*) 'Problem with mbar in ff_mqq0, mbar=',mbar
            stop
        endif

        mbarsq=mbar**2
        ff_mqq0=epinv*epinv2/2d0
     &         +epinv*(1d0-L/2d0+dlog(mbarsq)/2d0-dlog(1d0-mbarsq))
     &         -L+3d0+L**2/4d0-pisq*5d0/12d0
     &         -ddilog(1d0-mbarsq)
     &         -dlog(mbarsq)*(dlog(mbarsq)/2d0-1d0+L)/2d0
     &         -dlog(1d0-mbarsq)*(dlog(mbarsq/(1d0-mbarsq))-L+2d0)
     &         -mbarsq/(1d0-mbarsq)*dlog(mbarsq)
	
      endif
      
      return
      end

