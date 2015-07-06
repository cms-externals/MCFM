      logical function includedipole(nd,ptrans)
c--- This function returns TRUE if the specified point ptrans,
c--- corresponding to dipole nd (nd=0 => real radiation),
c--- should be included 
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'process.f'
      include 'frag.f'
      double precision ptrans(mxpart,4),pjet(mxpart,4),rcut
      integer i,j,nd,nqcdjets,nqcdstart,notag,isub
      logical gencuts,failedgencuts,photoncuts,makecuts,filterWbbmas,
     & photonfailed,filterW_bjet
      integer ij,count_photo,nphotons
      character*30 runstring
      double precision y32,y43,z3,z4,z5,z6
      double precision dphizj,pt5sq,pt6sq,pt7sq
      logical passed_frix,iso
      logical phot_dip(mxpart)

      common/phot_dip/phot_dip     
      common/runstring/runstring
      common/nqcdjets/nqcdjets,nqcdstart
      common/rcut/rcut
      common/makecuts/makecuts
      common/notag/notag

c--- default: include this contribution
      includedipole=.true.

c--- isub=1 for dipole subtractions, isub=0 for real radiation
      if (nd .gt. 0) then
        isub=1
      else
        isub=0
      endif

      nphotons=count_photo(ptrans) 
      if (nphotons .gt. 0) then 
c--- Photons: Frixione isolation cuts if no fragmentation included 
         if (frag .eqv. .false.) then
            do j=3,mxpart 
               if(plabel(j).eq.'ga') then 
                  call frix(ptrans,passed_frix,j,isub)
                  if(passed_frix.eqv..false.) then 
                     includedipole=.false.
                     return 
                  endif
               endif
            enddo
            call genclustphotons(ptrans,rcut,pjet,isub)
c--- The block of code commented out below is now deprecated in favour
c--- of the single routine genclustphotons() called above
c          if    ((case .eq. 'dirgam') .or. (case .eq. 'gamjet'))then 
c            call basic_frix(ptrans,passed_frix,isub)
c          elseif((case .eq. 'Wgamma') .or. (case .eq. 'Zgamma')
c     &       .or.(case .eq. 'Wgajet') .or. (case .eq. 'Zgajet')) then
c            call basic_Vfrix(ptrans,passed_frix,isub)
c          elseif((case .eq. 'gamgam') .or. (case .eq. 'Higaga')
c     &     .or.  (case .eq. 'gmgmjt')) then 
c            call basic_di_frix(ptrans,passed_frix,isub)
c         elseif((case.eq.'WHgaga').or.(case.eq.'ZHgaga')) then 
c            call basic_VHgaga_frix(ptrans,passed_frix,isub)
c         elseif(case.eq.'qq_Hgg') then 
c            call basic_2jet_frix(ptrans,passed_frix,isub)
c	  else
c	    write(6,*) 'Unforeseen photon process in includedipole.f'
c	    write(6,*) 'Needs a new Frixione routine.'
c	    stop
c         endif
c             if (passed_frix .eqv. .false.) then 
c                includedipole=.false.
c                return
c          else
c            call genclustphotons(ptrans,rcut,pjet,isub)
c         endif
        else 
c--- Photons: not Frixione, need fragmentation and isolation
          call genclust2(ptrans,rcut,pjet,isub)
c---  Isolate photon
          do j=3,mxpart
            if (plabel(j).eq.'ga') then 
               if(nd.eq.0) then 
                  if (iso(ptrans,isub,.false.,j,nd) .eqv. .false.) then
                     includedipole=.false.
                     return 
                  endif
               else
             if (iso(ptrans,isub,phot_dip(nd),j,nd) .eqv. .false.) then
                     includedipole=.false.
                     return 
                  endif
               endif
            endif
          enddo
        endif
      else
c--- No photons: the usual case
         call genclust2(ptrans,rcut,pjet,isub)
      endif

c--- perform mass cuts
      call masscuts(pjet,*999)
          
c--- fill ptilde array as persistent storage for the jet momenta
      do j=1,4
        do i=1,npart+2
          ptildejet(nd,i,j)=pjet(i,j)
        enddo
      enddo
     
c--- for the Wbb process, we divide up contributions in a specific way;
c--- therefore filter events using special code and skip normal jet testing   
      if (case .eq. 'Wbbmas') then
        includedipole=filterWbbmas()
	if (includedipole .eqv. .false.) return
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)
          if (failedgencuts) includedipole=.false.
        endif
	goto 99
      endif
      

c--- for the Wb+X process, we divide up contributions in a specific way;
c--- therefore filter events using special code and skip normal jet testing   
      if (case .eq. 'W_bjet') then
        includedipole=filterW_bjet()
	if (includedipole .eqv. .false.) return
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)
          if (failedgencuts) includedipole=.false.
        endif
	goto 99
      endif
     
c--- if the number of jets is not correct, then do not include dipole
      if ((clustering .and. (jets .ne. nqcdjets-notag)
     &       .and. (inclusive .eqv. .false.)) .or.
     &    (clustering .and. (jets .lt. nqcdjets-notag)
     &       .and. (inclusive .eqv. .true.))) then
          includedipole=.false.
          return
      else
c--- otherwise, if it is correct, check the photon cuts if appropriate
c--- rescale photons for fragmentation processes
c	  if (rescale) call rescale_pjet(pjet)
        if (nphotons .gt. 0) then 
	  if (rescale) call rescale_pjet(pjet)
          photonfailed=photoncuts(pjet,isub,phot_dip(nd),nd)
c--- return photons to original values
	  if (rescale) call return_pjet(pjet)
	  if (photonfailed) then
	    includedipole=.false.
	    return
	  endif
	endif
c--- ... and then the lepton cuts, if necessary
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)	  
          if (failedgencuts) then
	    includedipole=.false.
	    return
	  endif
c          call writeout(ptrans)
c          pause
        endif
      endif
      
c      call Wgamcuts(ptrans,failedgencuts)
c      if (failedgencuts .eqv. .false.) includedipole=.false.
      
     
   99 continue
      
      return

  999 continue
      includedipole=.false.
      return    
      
      end
            




c--- This block of code may be reinstated above if required
cc--- special procedure for comparison with C. Oleari's
cc---  e+e- --> QQbg calculation; run Durham algorithm
c      if (runstring(1:5) .eq. 'carlo') then
c        call durhamalg(ptrans,npart-isub,y32,y43,z3,z4,z5,z6)
cc---     set up momentum array just in case it's needed elsewhere
c        do j=1,4
c          do i=1,npart+2
c          ptildejet(nd,i,j)=ptrans(i,j)
c          enddo
c        enddo
c	jets=1
c        if (y43 .gt. 0d0) jets=2
c	if     ((case .eq. 'qq_tbg') .or. (case .eq. 'epem3j')) then
cc---     only keep events that have 3 jets when ycut=rcut
c	  if ((y32 .lt. rcut) .or. (y43 .gt. rcut)) then
c	    includedipole=.false.
c	    jets=-1
c	  endif
c	  return
c	elseif (case .eq. 'qqtbgg') then
cc---     only keep events that have 4 jets when ycut=rcut
c	  if (y43 .lt. rcut) then
c	    includedipole=.false.
c	    jets=-1
c	  endif
c	  return
c	else
c	  write(6,*) 'Unexpected case in includedipole.f!'
c	  stop
c	endif
c      endif
c

c--- This block of code may be reinstated above if required
cc--- added extra check here, to allow for analysis of G. Hesketh et al.
cc--- that requires Z+2 jets with only one jet within cuts, to obtain
cc--- prediction for Delta_phi(Z,jet) at NLO
c      if (runstring(1:6) .eq. 'dphizj') then
c        if    (nqcdjets .eq. 1) then
c	  ij=5
c	elseif (nqcdjets .eq. 2) then
c	  pt5sq=pjet(5,1)**2+pjet(5,2)**2
c	  pt6sq=pjet(6,1)**2+pjet(6,2)**2
c	  if (pt5sq .gt. pt6sq) then
c            ij=5
c	  else
c	    ij=6
c	  endif
c	elseif (nqcdjets .eq. 3) then
c	  pt5sq=pjet(5,1)**2+pjet(5,2)**2
c	  pt6sq=pjet(6,1)**2+pjet(6,2)**2
c	  pt7sq=pjet(7,1)**2+pjet(7,2)**2
c	  if     (pt5sq .gt. max(pt6sq,pt7sq)) then
c            ij=5
c	  elseif (pt6sq .gt. max(pt5sq,pt7sq)) then
c	    ij=6
c	  else
c	    ij=7
c	  endif
c	endif
c        dphizj=atan2(pjet(3,1)+pjet(4,1),pjet(3,2)+pjet(4,2))
c     .        -atan2(pjet(ij,1),pjet(ij,2))
c        if (dphizj .gt. pi) dphizj=twopi-dphizj
c        if (dphizj .lt. -pi) dphizj=twopi+dphizj
c	if (abs(dphizj) .gt. pi-1d-1) includedipole=.false.
c      endif
c      
