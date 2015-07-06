c--- This subroutine is a basic jet clustering algorithim for 
c--- f(p1)+f(p2) --> gamma(p3)+f(p4)+f(p5) 


c---- Takes in pin, which has passed frixione cuts will cluster partons and determine whether 
c---- any partons lie in cone Rij < delta_0
      subroutine photgenclust(pin,Rmin,pfinal,isub,ipow) 
      implicit none
      include 'constants.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'npart.f'
      include 'frag.f'
      double precision pin(mxpart,4),Rmin,pjet(mxpart,4),
     .     pfinal(mxpart,4)
      double precision dijmin,dkmin,aetarap,pt,Rgen
      integer isub,i,nu,iter,ipow,nmin1,nmin2,maxjet,jetindex(mxpart) 
      integer nk,ajet
      logical jetmerge,insideacone
      integer photindex(npart),nphotons,j
      integer softjet
      common/jetmerge/jetmerge

c--- debug - need to generalise to arb photon 
      jets=0
      nphotons=0
      maxjet=0
      jetmerge =.false.

c      pjet=0d0

c---- Pick out jets 
      do i=3,npart+2-isub
         if ( (plabel(i) .eq. 'pp') .or. (plabel(i) .eq. 'pj') 
     &    .or.(plabel(i) .eq. 'bq') .or. (plabel(i) .eq. 'ba') 
     &    .or.(plabel(i) .eq. 'qj') )then
            maxjet=maxjet+1
            jetindex(maxjet)=i
            jetlabel(maxjet)=plabel(i) 
            do nu=1,4
               pjet(maxjet,nu)=pin(i,nu)
            enddo
         elseif (plabel(i) .eq. 'ga') then 
            nphotons=nphotons+1
            photindex(nphotons)=i
         endif
      enddo

      
      if (maxjet .eq. 0 ) then 
         do i =1,mxpart 
            do nu=1,4
               pfinal(i,nu)=pin(i,nu)
            enddo 
         enddo
         jets=0
         return 
      endif

      if (maxjet .eq. 1) goto 2

      iter = 0
           
 1    iter =iter +1 
      
      call findmind(pin,pjet,iter,maxjet,dijmin,nmin1,nmin2,ipow)
      
      call findminet(pin,pjet,iter,maxjet,dkmin,nk,ipow) 
      dkmin=dkmin*Rmin

      if (dijmin .lt. dkmin) then 
         jetmerge = .true. 
         call combine(pjet,nmin1,nmin2) 
         
         call swapjet(pjet,jetindex,nmin2,maxjet)
         maxjet=maxjet-1
         iter =iter-1

      else
         jets=jets+1
         call swapjet(pjet,jetindex,jets,nk)
      endif

      if (iter .lt. maxjet-1) goto 1


 2    continue 
      jets=jets+1

      do i=1,2
         do nu=1,4
            pfinal(i,nu)=pin(i,nu) 
         enddo 
      enddo

      do i=3,npart+2
         do nu=1,4
            pfinal(i,nu)=0d0
             if ( (plabel(i) .ne. 'pp') .and. (plabel(i) .ne. 'pj') 
     &      .and. (plabel(i) .ne. 'bq') .and. (plabel(i) .ne. 'ba') 
     &      .and. (plabel(i) .ne. 'qj') )then
                pfinal(i,nu)=pin(i,nu)
             endif
          enddo
       enddo
       
      

c       if(isub .eq. 0) then 
c      write(*,*) 'AFTER CLUSTERING: Obtained ',jets,' jets'
c      endif
       
      
       ajet=0
       softjet=0

c--- loop over all partons
       do i=1,jets 

c--- check to see whether parton is inside one of the photon cones
         insideacone=.false.
         do j=1,nphotons           

c           write(6,*) i,j,Rgen(pjet,i,pin,photindex(j)),cone_ang
           if(Rgen(pjet,i,pin,photindex(j)) .lt. cone_ang) then 
	     insideacone=.true.
	   endif
         enddo

         if (insideacone) then	       
c--- if passed frix and is in isolation cone then do not apply cuts
c--- will add to jet tally
           softjet=softjet+1
	 else
c--- if jet doesnt lie within photon cone apply cuts             
           if ((pt(i,pjet) .ge. ptjetmin) .and. 
     &         (aetarap(i,pjet) .ge. etajetmin) .and.
     &         (aetarap(i,pjet) .le. etajetmax)) then                  
             ajet=ajet+1
             do nu=1,4
               pfinal(jetindex(ajet),nu)=pjet(i,nu)
             enddo
           endif
         endif

       enddo

c       write(6,*) 'isub,jets,ajet,softjet',isub,jets,ajet,softjet
c       if (softjet .eq. 1) pause
         
c--- if no jets are removed by eta and pt cuts, then jets=ajet
       if (ajet .lt. jets) then
          do i=ajet+1,jets
             do nu=1,4
                pfinal(jetindex(i),nu)=0d0
             enddo
          enddo        
       endif

       jets=ajet
       
c      write(6,*) 'input momenta' 
c      call writeout(pin)
c      write(6,*) '*****************************************'
c      write(6,*) ' output momenta' 
c      call writeout(pfinal) 
c      pause
      
      return 
      end
