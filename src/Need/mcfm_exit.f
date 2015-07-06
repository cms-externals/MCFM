      subroutine mcfm_exit(itmx,xinteg,xinteg_err)
************************************************************************
*                                                                      *
*  This routine should perform the final processing and print-outs     *
*                                                                      *
************************************************************************
      implicit none
      include 'efficiency.f'
      include 'process.f'
      include 'PDFerrors.f'
      include 'part.f'
      integer j,k,itmx
      double precision xinteg,xinteg_err,minPDFxsec,maxPDFxsec
      double precision PDFerror,PDFperror,PDFnerror
      double precision lord_bypart(-1:1,-1:1),lordnorm,rescale
      double precision ggpart,gqpart,qgpart,qqpart,qqbpart,
     . gqbpart,qbgpart,qbqbpart,qbqpart
      logical creatent,dswhisto
      common/outputflags/creatent,dswhisto
      common/bypart/lord_bypart
      common/finalpart/ggpart,gqpart,qgpart,qqpart,qqbpart

      double precision PDFMCav, PDFMCer, sum1,sum2
      
c--- Print-out the value of the integral and its error
      write(6,*) 
      if (xinteg .lt. 5d7) then 
        write(6,53)'Value of final ',part,' integral is',
     .   xinteg,' +/-',xinteg_err, ' fb'
      else 
        write(6,53)'Value of final ',part,' integral is',
     .   xinteg/1d6,' +/-',xinteg_err/1d6, ' nb'
        write(6,*) '(WARNING: result in nanobarns)'
      endif

c--- for gg->H+X processes, also write out the cross section
c---  normalized by sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)
      if (((case(1:5) .eq. 'ggfus') .or. (case(1:3) .eq. 'HWW')
     & .or.(case(1:3) .eq. 'HZZ')) .and. (case .ne. 'HWWint')
     &  .and. (case .ne. 'HWW_tb') .and. (case .ne. 'HZZint')
     &  .and. (case .ne. 'HZZ_tb') ) then
        call finitemtcorr(rescale)
        write(6,*)
	write(6,*) 'Cross section normalized by the ratio'
	write(6,*) 'sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)'
	write(6,*) '(i.e. exact for gg->H process, but '//
     .               'approx. for gg->H+n jets, n=1,2,3)'
        write(6,*)
        write(6,53)' Rescaled ',part,' integral is',
     .   xinteg*rescale,' +/-',xinteg_err*rescale, ' fb'   
        write(6,'(a25,f7.3,a2)') '   (Rescaling factor is ',rescale,')'  
      endif
     
   53 format(a15,a4,a12,f13.3,a4,f10.3,a3)

c--- Print-out a summary of the effects of jets and cuts
      write(6,*) 
      write(6,*) 'Total number of shots       : ',ntotshot
      write(6,*) 'Total no. failing cuts      : ',ntotzero
      write(6,*) 'Number failing jet cuts     : ',njetzero
      write(6,*) 'Number failing process cuts : ',ncutzero
      write(6,*) 
      call flush(6)

c--- Calculate the actual number of shots that were passed
c--- through the jet and cut routines
      ntotshot=ntotshot-(ntotzero-njetzero-ncutzero)
      write(6,54) 'Jet efficiency : ',
     .  100d0-100d0*dfloat(njetzero)/dfloat(ntotshot)
      write(6,54) 'Cut efficiency : ',
     .  100d0-100d0*dfloat(ncutzero)/dfloat((ntotshot-njetzero))
      write(6,54) 'Total efficiency : ',
     .  100d0-100d0*dfloat((njetzero+ncutzero))/dfloat(ntotshot)
      write(6,*) 
      
      lordnorm=0d0
      do j=-1,1
      do k=-1,1
        lordnorm=lordnorm+lord_bypart(j,k)
      enddo
      enddo
      ggpart=lord_bypart( 0, 0)/lordnorm
      gqpart=lord_bypart( 0,+1)/lordnorm
      gqbpart=lord_bypart( 0,-1)/lordnorm
      qgpart=lord_bypart(+1, 0)/lordnorm
      qbgpart=lord_bypart(-1, 0)/lordnorm
      qqpart=lord_bypart(+1,+1)/lordnorm
      qbqbpart=lord_bypart(-1,-1)/lordnorm
      qqbpart=lord_bypart(+1,-1)/lordnorm
      qbqpart=lord_bypart(-1,+1)/lordnorm
      write(6,*) 'Contribution from parton sub-processes:'
      write(6,*) '---------------------------------------'      
      write(6,55) '   GG    ',ggpart*xinteg,ggpart*100d0
      write(6,55) '   GQ    ',gqpart*xinteg,gqpart*100d0
      write(6,55) '   GQB   ',gqbpart*xinteg,gqbpart*100d0
      write(6,55) '   QG    ',qgpart*xinteg,qgpart*100d0
      write(6,55) '   QBG   ',qbgpart*xinteg,qbgpart*100d0
      write(6,55) '   QQ    ',qqpart*xinteg,qqpart*100d0
      write(6,55) '   QBQB  ',qbqbpart*xinteg,qbqbpart*100d0
      write(6,55) '   QQB   ',qqbpart*xinteg,qqbpart*100d0
      write(6,55) '   QBQ   ',qbqpart*xinteg,qbqpart*100d0
      write(6,*) '---------------------------------------'
      call flush(6)

   54 format(a20,f6.2,'%')
   55 format(4x,a9,' |',f18.5,f8.2,'%')

c--- If we've calculated PDF errors, present results.   
c--- Note that asymmetric errors are calculated according to
c--- Eq. (43) of "Hard Interactions of Quarks and Gluons",
c---  J.Campbell, J.Huston, W.J. Stirling, Rep. Prog. Phys. 70 (2007) 89
c--- (called "HEPDATA" method below)
      if (PDFerrors) then
        open(unit=91,status='unknown',file='pdferrors.res')
        write(6,*)
        write(6,58) '************ PDF error analysis ************'
        write(6,58) '*                                          *'
!     Compute PDF errors for sets which provide eigenvectors
        minPDFxsec=PDFxsec(0)
        maxPDFxsec=PDFxsec(0)
        PDFerror=0d0
        do j=0,maxPDFsets
          write(6,56) j,PDFxsec(j)
          if     (PDFxsec(j) .lt. minPDFxsec) then
            minPDFxsec=PDFxsec(j)
          elseif (PDFxsec(j) .gt. maxPDFxsec) then
            maxPDFxsec=PDFxsec(j)
          endif
          if ( (j .gt. 0) .and. (j/2 .eq. (j-1)/2) ) then
            PDFerror=PDFerror+(PDFxsec(j)-PDFxsec(j+1))**2
          endif
        enddo
	PDFperror=0d0
	PDFnerror=0d0
	do j=1,maxPDFsets-1,2
	  PDFperror=PDFperror+max(
     .     PDFxsec(j)-PDFxsec(0),PDFxsec(j+1)-PDFxsec(0),0d0)**2
	  PDFnerror=PDFnerror+max(
     .     PDFxsec(0)-PDFxsec(j),PDFxsec(0)-PDFxsec(j+1),0d0)**2
        enddo
        PDFerror=0.5d0*dsqrt(PDFerror)
        PDFperror=dsqrt(PDFperror)
        PDFnerror=dsqrt(PDFnerror)

c--- Compute PDF errors with the MC prescription
c---  (see Appendix B of arXiv:0808.1231 [hep-ph])
        sum1=0d0
        sum2=0d0

        do j=1,maxPDFsets
           sum1=sum1+PDFxsec(j)
           sum2=sum2+PDFxsec(j)**2d0
        enddo

        PDFMCav = sum1/maxPDFsets
        PDFMCer =dsqrt( sum2/maxPDFsets -  PDFMCav**2d0 )

        write(6,58) '*                                          *'
        write(6,58) '* --------------- SUMMARY ---------------- *'
        write(6,58) '*                                          *'
        write(6,58) '*            HEPDATA prescription          *'
	write(6,58) '*     (see, for example Eqn. (43) of       *'
	write(6,58) '*      J.Campbell, J.Huston, W.J.Stirling, *'
	write(6,58) '*      Rep. Prog. Phys. 70 (2007) 89)      *'
        write(6,58) '*                                          *'
        write(6,57) 'Minimum value',minPDFxsec
        write(6,57) 'Central value',PDFxsec(0)
        write(6,57) 'Maximum value',maxPDFxsec
        write(6,58) '*                                          *'
        write(6,57) 'Err estimate +/-',PDFerror
        write(6,57) '   +ve direction',PDFperror
        write(6,57) '   -ve direction',PDFnerror
        write(6,59) 'Fractional error',PDFerror/PDFxsec(0)
        write(6,58) '*                                          *'
        write(6,58) '*              MC prescription             *'
        write(6,58) '*       (for details and references,       *'
        write(6,58) '*        see Eqn. (158) in Appendix B      *'
        write(6,58) '*        of arXiv:0808.1231 [hep-ph])      *'
        write(6,58) '*                                          *'
        write(6,57) 'Central value',PDFMCav
        write(6,57) 'Err estimate +/-',PDFMCer
        write(6,59) 'Fractional error',PDFMCer/PDFMCav
        write(6,58) '********************************************'

        write(91,58) '* --------------- SUMMARY ---------------- *'
        write(91,58) '*                                          *'
        write(91,58) '*            HEPDATA prescription          *'
	write(91,58) '*     (see, for example Eqn. (43) of       *'
	write(91,58) '*      J.Campbell, J.Huston, W.J.Stirling, *'
	write(91,58) '*      Rep. Prog. Phys. 70 (2007) 89)      *'
        write(91,58) '*                                          *'
        write(91,57) 'Minimum value',minPDFxsec
        write(91,57) 'Central value',PDFxsec(0)
        write(91,57) 'Maximum value',maxPDFxsec
        write(91,58) '*                                          *'
        write(91,57) 'Err estimate +/-',PDFerror
        write(91,57) '   +ve direction',PDFperror
        write(91,57) '   -ve direction',PDFnerror
        write(91,59) 'Fractional error',PDFerror/PDFxsec(0)
        write(91,58) '*                                          *'
        write(91,58) '*              MC prescription             *'
        write(91,58) '*       (for details and references,       *'
        write(91,58) '*        see Eqn. (158) in Appendix B      *'
        write(91,58) '*        of arXiv:0808.1231 [hep-ph])      *'
        write(91,58) '*                                          *'
        write(91,57) 'Central value',PDFMCav
        write(91,57) 'Err estimate +/-',PDFMCer
        write(91,59) 'Fractional error',PDFMCer/PDFMCav
        write(91,58) '********************************************'

   
      endif
      close(91)
      

   56 format('* PDF error set ',i3,' -->',f15.3,' fb  *')
   57 format('*   ',a16,f14.3,' fb      *')
   58 format(a44)
   59 format('*   ',a16,f14.3,'         *')
 
c--- Finalize the histograms, if we're not filling ntuples instead
      if (creatent .eqv. .false.) then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call histofin(xinteg,xinteg_err,0,itmx)
        else
c--- DSW histograms - store the information
          call dswhbook(200,'Sigma   ',1.0d0,0.0d0,10.0d0)
          call dswhfill(200,0.5d0,xinteg)
          call dswhfill(200,1.5d0,xinteg_err)
c--- DSW histograms - output and close file
          call NTfinalize
        endif
      else
c--- ADDED - to produce normal histograms as well
        call histofin(xinteg,xinteg_err,0,itmx)
        call NTfinalize
      endif

      return
      
      end
