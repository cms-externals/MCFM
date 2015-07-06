      logical function gencuts(pjet,njets)
************************************************************************
*   Author: J.M. Campbell, 5th December 2001                           *
*                                                                      *
*   This routine imposes a generic set of cuts that can be applied     *
*   to all processes in process.DAT, using the parton momenta in pjet  *
*   which have already passed through the jet clustering algorithm     *
*                                                                      *
*   Only a basic set of variables is tested:                           *
*     pt(lepton) > leptpt, eta(lepton) < leptrap, missing Et > misspt  *
*                                                                      *
*   If leptpt2 or leptrap2 is not zero, then any leptons beyond the    *
*   leading-pt one must instead satisfy:                               *
*     pt(lepton) > leptpt2, eta(lepton) < leptrap2                     *
*                                                                      *
*   For processes where one would like to apply jet-like cuts, but     *
*   no clustering has been performed, additional cuts apply:           *
*     pt(jet) > jetpt, eta(jet) < jetrap, R(jet1,jet2) > Rcut          *
*                                                                      *
*   Finally, if further process-specific cuts are necessary,           *
*   an appropriate second routine may be called                        *
*                                                                      *
*   Return TRUE if this point FAILS the cuts                           *
*                                                                      *
************************************************************************
      implicit none
      include 'bbproc.f'
      include 'constants.f'
      include 'jetcuts.f'
      include 'leptcuts.f'
      include 'process.f'
      include 'plabel.f'
      logical first,passedlept
      integer njets,j,k,countb,bindex(mxpart),jindex,kindex,ib1,ib2
      integer countlept,leptindex(mxpart),countgamm,
     & countjet,jetindex(mxpart),pntr,maxparts,notag
      double precision pjet(mxpart,4),etvec(4)
      double precision pt,etarap,etmiss,evtmisset,R,Rcut,etaj,etak,
     & etalept,mll,jetpt,jetrap
      double precision etabuffer,m3456
      logical passed,hwwjetcuts
      double precision ht,qeta,mlbnu,merecon,reconcorr
      double precision dphi_ll,m_ll,mtrans,scut1,scut2
      double precision pttwo,phill,phillcut,etajet2cut,mllcut
      double precision deltaeta,deltaphi,pt5sq,pt6sq,ptthree
      integer ij,i,ilep,igam,inu
      character*30 runstring
      common/runstring/runstring
      common/stopvars/ht,qeta,mlbnu,merecon,reconcorr
      common/hwwvars/dphi_ll,m_ll,mtrans,scut1,scut2
      common/rcut/Rcut
      common/notag/notag
************************************************************************
*     Set-up the jet-like cut parameters here                          *
      parameter (jetpt=15d0,jetrap=2d0)
************************************************************************
      parameter(phillcut=1.2d0,etajet2cut=2.5d0,mllcut=75d0)
      data first/.true./
      
      gencuts=.false.
      
      hwwjetcuts=.false.
      
      
c--- THIS FUNCTIONALITY REMOVED, TO AVOID CONFUSION
c      if (case .eq. 'Zgamma') then
cc--- DEBUG: cut on pt of Z
c      if (pttwo(3,4,pjet) .lt. 200d0) then
c        gencuts=.true.
c	return
c      endif
cc--- DEBUG: cut on maximum m345
c      m3456=
c     &    +2d0*(pjet(3,4)*pjet(4,4)-pjet(3,1)*pjet(4,1)
c     &         -pjet(3,2)*pjet(4,2)-pjet(3,3)*pjet(4,3)) 
c     &    +2d0*(pjet(3,4)*pjet(5,4)-pjet(3,1)*pjet(5,1)
c     &         -pjet(3,2)*pjet(5,2)-pjet(3,3)*pjet(5,3)) 
c     &    +2d0*(pjet(4,4)*pjet(5,4)-pjet(4,1)*pjet(5,1)
c     &         -pjet(4,2)*pjet(5,2)-pjet(4,3)*pjet(5,3)) 
c      m3456=dsqrt(m3456)
c      if (m3456 .gt. 4000d0) then
c        gencuts=.true.
c	return
c      endif
c      endif
            
c--- perform extra H(->WW)+jet search cuts if there is a minimum
c---  jet rapidity, as is usually done
c      if (abs(etajetmin) .gt. 1d-6) then
c        hwwjetcuts=.true.
c      else
c        hwwjetcuts=.false.
c      endif 
     

c--- CMS Higgs search cuts
c      if(runstring(1:8).eq.'CMSHiggs') then 
c         call cms_higgsWW(pjet,gencuts) 
c         return 
c      endif
      
c--- ATLAS HIGGS SEARCH CUTS 
c      if(runstring(1:10).eq.'ATLASHiggs') then 
c         call ATLAS_hww(pjet,gencuts) 
c         return 
c      endif


c--- Special runstring for CDF W+2j cuts
c      if(runstring(1:10).eq.'cdf_Wdijet') then 
c         call cdf_Wdijet_v3(pjet,gencuts)
c         return 
c      endif

cc--- additional cuts for CDF Mjj analysis
c      if (runstring(1:3) .eq. 'mjj') then
cc--- transverse mass of the W > 30
c        mtrans=
c     .   (pjet(3,1)*pjet(4,1)+pjet(3,2)*pjet(4,2))
c     .   /dsqrt((pjet(3,1)**2+pjet(3,2)**2)
c     .         *(pjet(4,1)**2+pjet(4,2)**2))
c        mtrans=2d0*dsqrt(pjet(3,1)**2+pjet(3,2)**2)
c     .   *dsqrt(pjet(4,1)**2+pjet(4,2)**2)*(1d0-mtrans)
c        mtrans=dsqrt(max(mtrans,0d0))
c        if (mtrans .lt. 30d0) then
c	  gencuts=.true.
c	  return
c	endif       
c--- DeltaEta(J1,J2)| < 2.5
c        deltaeta=etarap(5,pjet)-etarap(6,pjet)
c	if (abs(deltaeta) .gt. 2.5d0) then
c	  gencuts=.true.
c	  return
c	endif
cc--- Pt(J1,J2) > 40 (Pt of vectorial sum of the two jets > 40) 
cc--- (this cut may be avoided)
c	if (runstring(4:9) .ne. 'noptjj') then
c	  if (pttwo(5,6,pjet) .lt. 40d0) then
c	    gencuts=.true.
c	    return
c	  endif
c	endif
cc--- DeltaPhi(J1, MET) >0.4
c	pt5sq=pjet(5,1)**2+pjet(5,2)**2
c	pt6sq=pjet(6,1)**2+pjet(6,2)**2
c	if (pt5sq .gt. pt6sq) then
c          ij=5
c	else
c	  ij=6
c	endif
c        deltaphi=atan2(pjet(3,1),pjet(3,2))-atan2(pjet(ij,1),pjet(ij,2))
c        if (deltaphi .gt. pi) deltaphi=twopi-deltaphi
c        if (deltaphi .lt. -pi) deltaphi=twopi+deltaphi
c	if (abs(deltaphi) .lt. 0.4d0) then
c	  gencuts=.true.
c	  return
c	endif
c      endif

cc--- extra transverse mass cut in W+jets for CDF
c      if (runstring(1:7) .eq. 'cdfjoey') then
c        mtrans=
c     .   (pjet(3,1)*pjet(4,1)+pjet(3,2)*pjet(4,2))
c     .   /dsqrt((pjet(3,1)**2+pjet(3,2)**2)
c     .         *(pjet(4,1)**2+pjet(4,2)**2))
cc---    transverse mass calculation
c        mtrans=2d0*dsqrt(pjet(3,1)**2+pjet(3,2)**2)
c     .   *dsqrt(pjet(4,1)**2+pjet(4,2)**2)*(1d0-mtrans)
c        mtrans=dsqrt(max(mtrans,0d0))
c        if (mtrans .lt. 20d0) then
c	  gencuts=.true.
c	  return
c	endif       
c      endif
      
c      if (runstring(1:4) .eq. 'stop') then
cc--- do single-top search cuts instead
c        maxparts=4+njets
c        if ((case .eq. 'tt_bbl') .or. (case .eq. 'tt_bbh')
c     .  .or.(case .eq. 'tt_bbu')) maxparts=6+njets
c        call stopcuts(pjet,maxparts,ht,qeta,mlbnu,merecon,reconcorr)  
c        if (ht .lt. 0d0) gencuts=.true.
c        return
c      endif
       
c      if (runstring(1:6) .eq. 'cdfhww') then
cc--- do CDF H->WW search cuts instead
c        maxparts=6+njets
c        call cdfhwwcuts(pjet,maxparts,passed)  
c        if (passed .eqv. .false.) gencuts=.true.
c        return
c      endif
       
c      if (runstring(1:3) .eq. 'hww') then
cc--- do H->WW search cuts instead
c        maxparts=6+njets
c        call hwwcuts(pjet,maxparts,dphi_ll,m_ll,mtrans,scut1,scut2)  
c        if (mtrans .lt. 0d0) gencuts=.true.
c        return
c      endif
       
c      if (runstring(1:3) .eq. 'wbf') then
cc--- do WBF search cuts instead
c        maxparts=4+njets
c	if (runstring(4:8) .eq. 'jeppe') then
c          call wbfcuts_jeppe(pjet,maxparts,passed)  
c        else
c	  call wbfcuts(pjet,maxparts,passed)  
c	endif
c        if (passed .eqv. .false.) gencuts=.true.
c        return
c      endif
 
       
c--- Look for particles that should be treated as jets,
c--- so far only b decays from Z-bosons and
c--- hadronic decay of the W in diboson processes
c--- If these are present, we will do additional cuts
      countb=0
      do j=3,mxpart
         if ((plabel(j) .eq. 'qb') .or. (plabel(j) .eq. 'ab')
     .  .or. (plabel(j) .eq. 'qq') .or. (plabel(j) .eq. 'qa')) then
           countb=countb+1
           bindex(countb)=j
         endif
      enddo

c--- write-out the cuts we are using
      if (first) then
      first=.false.
      
      write(6,*)
      write(6,*)  '****************** Generic cuts ********************'
      write(6,*)  '*                                                  *'
      write(6,99) '*        pt(lepton)      >   ',leptpt,
     .                ' GeV            *'
      write(6,99) '*      |eta(lepton)|     <   ',leptrap,
     .                '                *'
      write(6,99) '*       pt(missing)      >   ',misspt,
     .                ' GeV            *'
      if ((leptpt2 .ne. 0d0) .or. (leptrap2 .ne. 0d0)) then
      write(6,99) '*     pt(2nd+ lepton)    >   ',leptpt2,
     .                ' GeV            *'
      write(6,99) '*   |eta(2nd+ lepton)|   <   ',leptrap2,
     .                '                *'
      endif
      if     (case.eq.'Wgamma') then
      if (mtrans34cut .lt. 0d0) then
      write(6,99) '*   (3,4,5) trans. mass  >   ',abs(mtrans34cut),
     .                ' GeV            *'
      else
      write(6,99) '* (e-gam,nu) trans. mass >   ',mtrans34cut,
     .                ' GeV            *'
      endif
      elseif (case.eq.'Zgamma') then
      write(6,99) '*    (3,4,5) inv. mass   >   ',mtrans34cut,
     .                ' GeV            *'
      else
      write(6,99) '*  (3,4) transverse mass >   ',mtrans34cut,
     .                ' GeV            *'
      endif
      write(6,99) '*      R(jet,lepton)     >   ',Rjlmin,
     .                '                *'
      write(6,99) '*     R(lepton,lepton)   >   ',Rllmin,
     .                '                *'
      write(6,99) '* |eta(jet1)-eta(jet2)|  >   ',delyjjmin,
     .                '                *'
      if (jetsopphem) then
      write(6,*) '*           eta(jet1) . eta(jet2)  <  0            *' 
      endif
      if     (lbjscheme .eq. 1) then
      write(6,*) '*        eta(jet1)  <  eta(lept)  <  eta(jet2)     *' 
      elseif (lbjscheme .ge. 2) then
      write(6,*) '*  eta(jet1)+Rcut  <  eta(lept)  <  eta(jet2)-Rcut *' 
      lbjscheme=2
      endif
c      write(6,99) '*   pt(photon)           >   ',gammpt,
c     .                '                *'
c      write(6,99) '*   eta(photon)          <   ',gammrap,
c     .                '                *'
c      write(6,99) '*     R(photon,lepton)   >   ',Rgalmin,
c     .                '                *'
c      write(6,99) '*   pt(hadronic)         <   ',gammcut,
c     .                '    pt(photon)  *'
c      write(6,99) '*   (in cone around photon of',gammcone,
c     .                ')               *'
      if (countb .gt. 0) then
      write(6,*)  '*                                                  *'
      write(6,99) '*      pt(jet)       >   ',jetpt,
     .                ' GeV                *'
      write(6,99) '*    |eta(jet)|      <   ',jetrap,
     .                '                    *'
      write(6,99) '*   R(jet1,jet2)     >   ',Rcut,
     .                '                    *'
      endif
      if (hwwjetcuts) then
      write(6,99) '*   phi(lepton,lepton)   <   ',phillcut,
     .                '                *'
      write(6,99) '*  |eta(2nd jet)|        >   ',etajet2cut,
     .                '                *'
      write(6,99) '*  m(lepton,lepton)      <   ',mllcut,
     .                '                *'
      endif
      write(6,*)  '****************************************************'
      endif
c-- debug compare against JETPHOX
c      if(case.eq.'dirgam') then 
c         if((pt(3,pjet).lt.14.5d0).or.(pt(3,pjet).gt.15.5d0)) then 
c            gencuts=.true. 
c            return 
c         endif
c      endif

      countlept=0
      countgamm=0
c--- lepton pt and rapidity cuts
      do j=3,mxpart
        if (     (plabel(j) .eq. 'el') .or. (plabel(j) .eq. 'ea')
     .      .or. (plabel(j) .eq. 'ml') .or. (plabel(j) .eq. 'ma')
     .      .or. (plabel(j) .eq. 'tl') .or. (plabel(j) .eq. 'ta')) then
          countlept=countlept+1
          leptindex(countlept)=j
        endif
        if (plabel(j) .eq. 'ga') then
          countgamm=countgamm+1
        endif
      enddo
      
C     Basic pt and rapidity cuts for lepton
      if     (countlept .eq. 1) then
          if (     (pt(leptindex(1),pjet) .lt. leptpt) .or.
     .      (abs(etarap(leptindex(1),pjet)) .gt. leptrap)) then
            gencuts=.true.
            return
          endif
      elseif (countlept .gt. 1) then
c--- loop over all the lepton possibilities for lepton 1 (j)
          j=0
  77      continue
          j=j+1
          passedlept=.true.
          if (     (pt(leptindex(j),pjet) .lt. leptpt) .or.
     .      (abs(etarap(leptindex(j),pjet)) .gt. leptrap)) then
            passedlept=.false.
            goto 78
          endif
          do k=1,countlept
              if (k .ne. j) then
                  if (     (pt(leptindex(k),pjet) .lt. leptpt2) .or.
     .              (abs(etarap(leptindex(k),pjet)) .gt. leptrap2)) then
                    passedlept=.false.
                   endif
              endif 
          enddo         
  78      continue
c--- return to beginning if we failed and there are more leptons to try
          if ((passedlept .eqv. .false.).and.(j .lt. countlept)) goto 77
          gencuts=.not.(passedlept)
      endif

c--- These cuts are now performed in photoncuts.f
C     Basic pt and rapidity cuts for photon
c      if (countgamm .gt. 0) then
c            do j=1,countgamm
c              if (     (pt(gammindex(j),pjet) .lt. gammpt) .or.
c     .          (abs(etarap(gammindex(j),pjet)) .gt. gammrap)) then
c                gencuts=.true.
c                return
c              endif
c             enddo
c      endif

c--- This cut is now performed in photoncuts.f
c--- lepton-photon separation 
c      if ((countlept .ge. 1) .and. (countgamm .ge. 1)) then
c        do j=1,countgamm
c        do k=1,countlept
c          if (R(pjet,gammindex(j),leptindex(k)) .lt. Rgalmin) then
c            gencuts=.true.
c            return
c          endif
c        enddo
c        enddo
c      endif
 
c--- missing energy cut
      evtmisset=etmiss(pjet,etvec)
      if ((evtmisset .lt. misspt) .and. (evtmisset .ne. 0d0)) then
        gencuts=.true.
        return
      endif
      
c--- mtrans34cut is used for three roles:
c---  1) Wgamma    --> mtrans34cut<0: transverse mass cut on (3,4,5) system
c---               --> mtrans34cut<0: transverse mass cut on (e-gam,nu) system
c---  2) Zgamma    --> invariant mass cut on (3,4,5) system
c---  3) otherwise --> transverse mass cut on (3,4) system
c---
      if (case.eq.'Wgamma') then 
c--- cut on transverse mass of (3,4,5) system for Wgamma
        if (mtrans34cut .lt. 0d0) then
        mtrans=0d0
        do j=3,5
           mtrans=mtrans+dsqrt(pjet(j,1)**2+pjet(j,2)**2)
        enddo
        mtrans=mtrans**2
        do j=1,2 
           mtrans=mtrans-(pjet(3,j)+pjet(4,j)+pjet(5,j))**2
        enddo
        mtrans=dsqrt(max(mtrans,0d0)) 
        if (mtrans .lt. abs(mtrans34cut)) then
           gencuts=.true.
           return
        endif
	else
c--- cut on (e-gam,nu) transverse mass for Wgamma
        if (plabel(3) .eq. 'nl') then
	  inu=3
	  ilep=4
	else
	  inu=4
	  ilep=3
	endif
	igam=5
	mtrans=(pjet(ilep,4)+pjet(igam,4))**2
     &        -(pjet(ilep,1)+pjet(igam,1))**2
     &        -(pjet(ilep,2)+pjet(igam,2))**2
     &        -(pjet(ilep,3)+pjet(igam,3))**2
        mtrans=mtrans+(pjet(ilep,1)+pjet(igam,1))**2
     &               +(pjet(ilep,2)+pjet(igam,2))**2
	mtrans=dsqrt(max(mtrans,0d0))
     &        +dsqrt(pjet(inu,1)**2+pjet(inu,2)**2)
	mtrans=mtrans**2
        do j=1,2 
           mtrans=mtrans-(pjet(3,j)+pjet(4,j)+pjet(5,j))**2
        enddo
        mtrans=dsqrt(max(mtrans,0d0)) 
        if (mtrans .lt. mtrans34cut) then
           gencuts=.true.
           return
        endif
	endif
c--- cut on invariant mass of (3,4,5) system for Zgamma
      elseif (case.eq.'Zgamma') then 
	mtrans=(pjet(3,4)+pjet(4,4)+pjet(5,4))**2
        do j=1,3 
          mtrans=mtrans-(pjet(3,j)+pjet(4,j)+pjet(5,j))**2 
        enddo
        mtrans=dsqrt(max(mtrans,0d0)) 
        if(mtrans.lt.mtrans34cut) then
           gencuts=.true.
           return
        endif
      else
c--- cut on transverse mass of (3,4) pair otherwise
        mtrans=
     &   (pjet(3,1)*pjet(4,1)+pjet(3,2)*pjet(4,2))
     &   /dsqrt((pjet(3,1)**2+pjet(3,2)**2)
     &         *(pjet(4,1)**2+pjet(4,2)**2))
        mtrans=2d0*dsqrt(pjet(3,1)**2+pjet(3,2)**2)
     &   *dsqrt(pjet(4,1)**2+pjet(4,2)**2)*(1d0-mtrans)
        mtrans=dsqrt(max(mtrans,0d0))
        if (mtrans .lt. mtrans34cut) then
          gencuts=.true.
        return
        endif
      endif
        
c--- lepton-lepton separation (if there are 2 or more leptons)
      if ((countlept .gt. 1)) then
        do j=1,countlept
        do k=j+1,countlept
          if (R(pjet,leptindex(j),leptindex(k)) .lt. Rllmin) then
            gencuts=.true.
            return
          endif
        enddo
        enddo
c--- extra cut on phi(lept,lept) for H(->WW)+jet search
	if (hwwjetcuts) then
          phill=
     .       (pjet(leptindex(1),1)*pjet(leptindex(2),1)
     .       +pjet(leptindex(1),2)*pjet(leptindex(2),2))
     .       /dsqrt((pjet(leptindex(1),1)**2+pjet(leptindex(1),2)**2)
     .             *(pjet(leptindex(2),1)**2+pjet(leptindex(2),2)**2))
          if (phill .lt. -0.999999999D0) phill=-1d0
          phill=dacos(phill) 
	  if (phill .gt. phillcut) then
	    gencuts=.true.
	    return
	  endif
	endif
c--- extra cut on m(lept,lept) for H(->WW)+jet search
	if (hwwjetcuts) then
          mll=dsqrt(2d0*(
     .       +pjet(leptindex(1),4)*pjet(leptindex(2),4)
     .       -pjet(leptindex(1),1)*pjet(leptindex(2),1)
     .       -pjet(leptindex(1),2)*pjet(leptindex(2),2)
     .       -pjet(leptindex(1),3)*pjet(leptindex(2),3)))
	  if (mll .gt. mllcut) then
	    gencuts=.true.
	    return
	  endif
	endif
      endif
 
c--- if there are no cuts on the jets - or no jets - we are done      
      if ((Rjlmin .le. 0d0) .and. (delyjjmin .le. 0d0)) return
      if ((njets .eq. 0) .and. (countb .eq. 0)) return
      
c--- identify the jets
      countjet=0      
      do j=3,mxpart
        if (     (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     .      .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
          countjet=countjet+1
          jetindex(countjet)=j
        endif
      enddo

c--- countjet will pick up the extra 'pp' needed for the real piece,
c--- therefore we should subtract 1 from this number     
      if (countjet .gt. njets) countjet=countjet-1

      if ((njets .ne. countjet) .and. (notag .eq. 0)) then
        write(6,*) 'Something is wrong in gencuts.f -'
        write(6,*) 'countjet = ',countjet,' BUT njets = ',njets
        stop
      endif
      
c--- extra cut on eta(2nd jet) for H(->WW)+jet search
      if ((hwwjetcuts) .and. (countjet .ge. 2)) then
        if (pt(jetindex(1),pjet) .gt. pt(jetindex(2),pjet)) then
          if (abs(etarap(jetindex(2),pjet)) .lt. etajet2cut) then
            gencuts=.true.
            return
          endif
        else
          if (abs(etarap(jetindex(1),pjet)) .lt. etajet2cut) then
            gencuts=.true.
            return
          endif
        endif
      endif

c--- jet-lepton separation (if there are 1 or more jets and leptons)
      if ((njets .gt. 0) .and. (countlept .gt. 0)) then
        do j=1,countlept
        do k=1,njets
          if (R(pjet,leptindex(j),jetindex(k)) .lt. Rjlmin) then
            gencuts=.true.
            return
          endif
        enddo
        enddo
      endif
      
c--- jet-photon separation (if there are 1 or more jets and photons)
c      if ((njets .gt. 0) .and. (countgamm .gt. 0)) then
c	do j=1,countgamm
c	do k=1,njets
c	  if (R(pjet,gammindex(j),jetindex(k)) .lt. Rjlmin) then
c	    gencuts=.true.
c	    return
c	  endif
c	enddo
c	enddo
c      endif

c--- DEBUG: removed all isolation      
cc--- photon/hadron isolation     
c      if ((njets .gt. 0) .and. (countgamm .gt. 0)) then
c        do j=1,countgamm
cc--- Frixione cut, hep-ph/9801442
c          if (njets .gt. 2) then
c            write(6,*) 'Photon-hadron isolation not coded for njets > 2'
c            stop
c          endif
c          do k=1,njets
c            delta(k)=R(pjet,gammindex(j),jetindex(k))
c            ptjet(k)=pt(jetindex(k),pjet)
c          enddo
c          pntr=1
c          if (delta(2) .lt. delta(1)) pntr=2
c          if (njets .eq. 1) delta(2)=gammcone   
c          discr=(1d0-dcos(delta(pntr)))/(1d0-dcos(gammcone))
c          if (ptjet(pntr) .gt. discr*pt(gammindex(j),pjet)) then
c            gencuts=.true.
c            return
c          endif 
c          if (njets .ge. 2) then
c            discr=(1d0-dcos(delta(3-pntr)))/(1d0-dcos(gammcone))
c            if (ptjet(1)+ptjet(2) .gt. discr*pt(gammindex(j),pjet))then
c              gencuts=.true.
c              return
c            endif
c          endif
	  
c--- this block was already removed
c--- optional jet-veto
c          if (   (pt(4+countgamm+1,pjet) .gt. 50d0)    
c     .     .and. (abs(etarap(4+countgamm+1,pjet)) .lt. 2.5d0)) then
c            gencuts=.true.
c          endif                     
c--- de-Florian,Signer cut
c          do nu=1,2
c            sumjetpt(nu)=0d0
c          enddo
c          do k=1,njets
c            if (R(pjet,gammindex(j),4+countgamm+k) .lt. gammcone) then
c              do nu=1,2
c              sumjetpt(nu)=sumjetpt(nu)+pjet(4+countgamm+k,nu)
c              enddo
c            endif
c          enddo
c          if ( dsqrt(sumjetpt(1)**2+sumjetpt(2)**2) 
c     .    .gt. gammcut*pt(gammindex(j),pjet) ) then
c            gencuts=.true.
c          endif
c--- this block was already removed

c        enddo
c      endif  
c--- DEBUG: removed all isolation      
      
c--- WBF-style cuts (if there are 2 or more jets)
      if ((njets .gt. 1)) then
c--- jet-jet rapidity separation 
c--- j and k point to the two highest pt ('tagging') jets
        j=1
        k=2
        if (njets .eq. 3) then
          if     ( pt(jetindex(1),pjet) .lt.
     .      min(pt(jetindex(2),pjet),pt(jetindex(3),pjet)) ) then
            j=2
            k=3
          elseif ( pt(jetindex(2),pjet) .lt.
     .      min(pt(jetindex(1),pjet),pt(jetindex(3),pjet)) ) then
            j=1
            k=3
          endif
        endif
        if (abs(etarap(jetindex(j),pjet)-etarap(jetindex(k),pjet))
     .         .lt. delyjjmin) then
          gencuts=.true.
          return
        endif

c--- Requirement that jets be in opposite hemispheres
        if (jetsopphem) then
          if(etarap(jetindex(j),pjet)*etarap(jetindex(k),pjet) .gt. 0d0)
     .       then
            gencuts=.true.
            return
          endif
        endif

        if (lbjscheme .gt. 0) then
c--- Cut to require lepton to be between jets
          etabuffer=dble(lbjscheme-1)*Rcut
          etaj=etarap(jetindex(j),pjet)
          etak=etarap(jetindex(k),pjet)
          do pntr=1,countlept
            etalept=etarap(leptindex(pntr),pjet)
            if ( (etalept .lt. min(etaj,etak)+etabuffer) .or.
     .           (etalept .gt. max(etaj,etak)-etabuffer) ) then 
              gencuts=.true.
              return
            endif

c--- these lines impose an alternative WBF selection
c             if ( ((etaj .gt. 1.6d0) .and. (etaj .lt. 4.4d0))
c     .        .or.((etak .gt. 1.6d0) .and. (etak .lt. 4.4d0)) ) then
c             else
c               gencuts=.true.
c             endif             
          enddo
          
          
        endif

      endif
      
c-- cuts on b-quarks
      if (bbproc) then
        call getbs(pjet,ib1,ib2)
        if ( (abs(etarap(ib1,pjet)) .gt. etabjetmax)
     .  .or. (pt(ib1,pjet) .lt. ptbjetmin) ) gencuts=.true. 
        if ( (abs(etarap(ib2,pjet)) .gt. etabjetmax)
     .  .or. (pt(ib2,pjet) .lt. ptbjetmin) ) gencuts=.true. 
      endif

c--- completed basic cuts
C--- if there are jet-like particles (see above), do more cuts
      if (countb .gt. 0) then
        do j=1,countb
          jindex=bindex(j)
          if (          (pt(jindex,pjet) .lt. jetpt) .or.
     .           (abs(etarap(jindex,pjet)) .gt. jetrap)) then
            gencuts=.true.
            return
          endif
          do k=j+1,countb
            kindex=bindex(k)
            if ((r(pjet,jindex,kindex) .lt. rcut)) then
              gencuts=.true.
              return
            endif
          enddo
        enddo
      endif
 
      return

   99 format(1x,a29,f6.2,a17)
      

      end
 
 
 
 
