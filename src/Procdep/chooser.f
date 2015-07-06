      subroutine chooser
c---- Note added 4/21/03
c---- plabel set to 'ig' (for 'ignore') means that this
c---- particle should not be subject to any cuts, so that the
c---- total cross-section comes out correctly when the BR is removed
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'vegas_common.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'bbproc.f'
      include 'nwz.f'
      include 'process.f'
      include 'flags.f'
      include 'heavyflav.f'
      include 'nflav.f'
      include 'nodecay.f'
      include 'stopscales.f'
      include 'scale.f'
      include 'facscale.f'
      include 'nlooprun.f'
      include 'b0.f'
      include 'colstruc.f'
      include 'stopbmass.f'
      include 'fourthgen.f'
      include 'anomcoup.f'
      include 'srdiags.f'
      include 'clustering.f'
      include 'frag.f'
      include 'plabel.f'
      include 'interference.f'
      include 'couple.f'
      include 'part.f'
      include 'hdecaymode.f'
      include 'breit.f'
      double precision wwbr,zzbr,tautaubr,gamgambr,zgambr,Rcut,Rbbmin,
     . alphas,cmass,bmass
      double precision br,BrnRat,brwen,brzee,brznn,brtau,brtop,brcharm
      integer nproc,mproc,j,nqcdjets,nqcdstart,isub,notag,imhq
      character*83 pname
      character*1 order
      character*72 string
      character*30 runstring
      double precision f0q,f2q,f4q
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/bitflags/f0q,f2q,f4q
      common/Rbbmin/Rbbmin
      common/Rcut/Rcut
      common/runstring/runstring
      common/nproc/nproc
      common/BrnRat/BrnRat
      common/nqcdjets/nqcdjets,nqcdstart
      common/isub/isub
      common/notag/notag
      common/qmass/cmass,bmass
      data hdecaymode/'xxxx'/
      do j=1,mxpart
      plabel(j)=''
      enddo

      string='process.DAT' 
      open(unit=21,file=string,status='old',err=43)
      call checkversion(21,string)
      
      write(6,*) 'Chooser:process chosen by nproc=',nproc

      do j=1,600
      read(21,*,err=44) mproc,pname,order
      
      if (nproc .lt. 0) then 
      write(6,*) mproc,pname 
      endif

      if (mproc .eq. nproc) go to 42
      if (pname .eq. 'EOF') go to 44
      enddo
      goto 44

 42   write(6,*)
      write(6,*) '*************************** f(p1)+f(p2) --> *****'//
     . '********************'
      write(6,*) '* ',pname(19:83),' *'
      write(6,*) '*************************************************'//
     . '********************'
      write(6,*)

      close(unit=21)

      plabel(1)='pp'
      plabel(2)='pp'

c--- the default behaviour is to remove no branching ratio
      BrnRat=1d0

c--- set up most parameters
      call coupling

      notag=0
      nqcdjets=0
      isub=0
      bbproc=.false.
      nodecay=.false.
      nfonly=.false.
      caonly=.false.
      fourthgen=.false.
      rescale=.false.

c-- Rbbmin is an additional variable, added so that the separation
c-- between two b jets can be controlled separately from the Delta_R
c-- cut between any other types of jet
c-- Default behaviour: the same value as for the other jets
      Rbbmin=Rcut
      
c-----------------------------------------------------------------------

      if ((nproc .eq. 1) .or. (nproc .eq. 6)) then
        case='W_only'
        mass3=wmass
        width3=wwidth
        n3=1
        ndim=4
        nqcdjets=0
c---W^+
        if     (nproc .eq. 1) then
C-- 1  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))'
C--    '  f(p1)+f(p2) --> W^+ (for total Xsect)' (removebr=.true.)
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='pp'
          nwz=1
c---W^-
        elseif (nproc .eq. 6) then
c-- 6  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))'
c--    '  f(p1)+f(p2) --> W^- (for total Xsect)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='pp'
          nwz=-1
        else
          call nprocinvalid()
        endif

c--- total cross-section
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 11) .or. (nproc .eq. 16)) then
        case='W_1jet'
        nqcdjets=1
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        plabel(5)='pp'
        plabel(6)='pp'

        if     (nproc .eq. 11) then
C-- 11 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + f(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 16) then
c-- 16 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+f(p5)'
c--    '  f(p1)+f(p2) --> W^- (no BR) + f(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 12) .or. (nproc .eq. 17)) then
        case='Wbfrmc'
        nqcdjets=1
        ndim=7
        n2=0
        n3=1
        nflav=4
        plabel(5)='bq'
        plabel(6)='pp'

        if     (nproc .eq. 12) then
c-- 12 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+bbar(p5)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + bbar(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 17) then
c-- 17 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+b(p5)'
c--    '  f(p1)+f(p2) --> W^- (no BR) + b(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
        mass3=wmass
        width3=wwidth
        mass2=mb
        if (Vcb .eq. 0d0) Vcb=0.041d0
        write(6,*) 'Setting Vcb=0.041 for this process'
        if (Vub .eq. 0d0) Vub=0.00347d0
        write(6,*) 'Setting Vub=0.00347 for this process'

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 13) .or. (nproc .eq. 18)) then
        case='W_cjet'
        nqcdjets=1
        ndim=7
        mb=0
        n2=0
        n3=1
        nflav=3
        plabel(5)='bq'
        plabel(6)='pp'

        if     (nproc .eq. 13) then
c-- 13 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+cbar(p5)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + cbar(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 18) then
c-- 18 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5)'
c--    '  f(p1)+f(p2) --> W^- (no BR) + c(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c--- change W mass (after couplings and BRs already calculated)
c	if     (runstring(4:8) .eq. 'mw_80') then
c	  wmass=80.4d0
c	elseif (runstring(4:8) .eq. 'mw200') then
c	  wmass=200d0
c	elseif (runstring(4:8) .eq. 'mw400') then
c	  wmass=400d0
c	endif
	
c--- change charm mass
c	if     (runstring(9:13) .eq. 'mc1.3') then
c	  mc=1.3d0
c	  mcsq=mc**2
c	elseif (runstring(9:13) .eq. 'mc4.5') then
c	  mc=4.5d0
c	  mcsq=mc**2
c	elseif (runstring(9:13) .eq. 'mc20.') then
c	  mc=20d0
c	  mcsq=mc**2
c	endif
	
        mass3=wmass
        width3=wwidth
        mass2=mc

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 14) .or. (nproc .eq. 19)) then
        case='Wcjet0'
        nqcdjets=1
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        mass2=0d0
        nflav=3
        plabel(5)='bq'
        plabel(6)='pp'

        if     (nproc .eq. 14) then
c-- 13 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+cbar(p5) [massless]'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + cbar(p5) [massless]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 19) then
c-- 18 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5) [massless]'
c--    '  f(p1)+f(p2) --> W^- (no BR) + c(p5) [massless]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 20) .or. (nproc .eq. 25)) then
        case='Wbbmas'
        write(6,*) 'mb=',mb
        nqcdjets=2
        flav=5
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
 
        if     (nproc .eq. 20) then
c-- 20 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+b~(p6) [massive]'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6) [massive]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 25) then
c-- 25 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5)+b~(p6) [massive]'
c--    '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6) [massive]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 21) .or. (nproc .eq. 26)) then
        case='Wbbbar'
        write(6,*) 'mb=',mb
        nqcdjets=2
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        mb=0
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 21) then
c-- 21 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+b~(p6)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 26) then
c-- 26 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5)+b~(p6)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6)' (removebr=.true.)
          nwz=-1 
          plabel(3)='el'
          plabel(4)='na'
        endif
             
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 22) .or. (nproc .eq. 27)) then
        case='W_2jet'
        nqcdjets=2
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 22) then
c-- 22 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)+f(p6)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +f(p5)+f(p6)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 27) then
c-- 27 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + f(p5)+f(p6)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +f(p5)+f(p6)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 23) .or. (nproc .eq. 28)) then
        case='W_3jet'
        nqcdjets=3
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=13
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 23) then
c-- 23 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)+f(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +f(p5)+f(p6)+f(p7)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 28) then
c-- 28 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + f(p5)+f(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +f(p5)+f(p6)+f(p7)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 24) .or. (nproc .eq. 29)) then
        case='Wbbjet'
        write(6,*) 'mb=',mb
        nqcdjets=3
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        mb=0d0
        ndim=13
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 24) then
c-- 24 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+b~(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6)+f(p7)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 29) then
c-- 29 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5)+b~(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6)+f(p7)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .gt. 30) .and. (nproc .le. 35)) then
        case='Z_only'
        nqcdjets=0
        nwz=0
        mass3=zmass
        width3=zwidth
        n3=1
        ndim=4
        plabel(5)='pp'

        if     (nproc .eq. 31) then
c-- 31 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))'
c--    '  f(p1)+f(p2) --> Z^0 (for total Xsect)' (removebr=.true.)
          call checkminzmass(1)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1d0
          l1=le
          r1=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee
          endif
        elseif (nproc .eq. 32) then
c-- 32 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4)))'
          plabel(3)='nl'
          plabel(4)='na'
          q1=0d0
          l1=ln*dsqrt(3d0)
          r1=rn*dsqrt(3d0)
        elseif (nproc .eq. 33) then
c-- 33 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))'
          call checkminzmass(1)
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        elseif (nproc .eq. 34) then
c-- 34 '  f(p1)+f(p2) --> Z^0(-->3*(d(p3)+d~(p4)))'
          call checkminzmass(1)
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(1)*dsqrt(3d0*xn)
          l1=l(1)*dsqrt(3d0*xn)
          r1=r(1)*dsqrt(3d0*xn)
        elseif (nproc .eq. 35) then
c-- 35 '  f(p1)+f(p2) --> Z^0(-->2*(u(p3)+u~(p4)))'
          call checkminzmass(1)
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(2)*dsqrt(2d0*xn)
          l1=l(2)*dsqrt(2d0*xn)
          r1=r(2)*dsqrt(2d0*xn)
        else
          call nprocinvalid()
        endif 

c-----------------------------------------------------------------------

        elseif (nproc .eq. 36) then
        case='ttZbbl'
        nwz=1
        ndim=16
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth
        bbproc=.true.
        
c--  36 '  f(p1)+f(p2) -> Z -> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+e^-(p7)+nu~(p8)'
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brwen*brtop)**2
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            plabel(7)='ig'
            plabel(8)='ig'
            nqcdjets=0
            bbproc=.false.
          endif


c-----------------------------------------------------------------------
      
      elseif ((nproc .ge. 41) .and. (nproc .le. 43)) then
        case='Z_1jet'
        nqcdjets=1
        nwz=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        plabel(5)='pp'
        plabel(6)='pp'

        if     (nproc .eq. 41) then
c-- 41 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +f(p5)' (removebr=.true.)
          call checkminzmass(1)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1d0
          l1=le
          r1=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee
          endif
        elseif (nproc .eq. 42) then
c-- 42 '  f(p1)+f(p2) --> Z_0(-->3*(nu(p3)+nu~(p4)))-(sum over 3 nu)+f(p5)'
            plabel(3)='nl'
            plabel(4)='na'
            q1=0d0
            l1=ln*dsqrt(3d0)
            r1=rn*dsqrt(3d0)
        elseif (nproc .eq. 43) then
c-- 43 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))+f(p5)'
          call checkminzmass(1)
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        endif
      
c-----------------------------------------------------------------------
      
      elseif (nproc .eq. 44) then
c-- 44 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)+f(p6)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +f(p5)+f(p6)' (removebr=.true.)
        case='Z_2jet'
        call checkminzmass(1)
        ndim=10
        n2=0
        n3=1
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        q1=-1d0
        l1=le
        r1=re
        nwz=0   
        mass3=zmass
        width3=zwidth
       
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
        endif

c--- added extra check here, to allow for analysis of G. Hesketh et al.
c--- that requires Z+2 jets with only one jet within cuts, to obtain
c--- prediction for Delta_phi(Z,jet) at NLO
c        if (runstring(1:6) .eq. 'dphizj') notag=1
	
      elseif (nproc .eq. 45) then
c-- 45 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)+f(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +f(p5)+f(p6)+f(p7)' (removebr=.true.)
        case='Z_3jet'
        call checkminzmass(1)
        ndim=13
        n2=0
        n3=1
        nqcdjets=3
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        q1=-1d0
        l1=le
        r1=re
        nwz=0   
        mass3=zmass
        width3=zwidth
       
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
        endif

c--- added extra check here, to allow for analysis of G. Hesketh et al.
c--- that requires Z+2 jets with only one jet within cuts, to obtain
c--- prediction for Delta_phi(Z,jet) at NLO
c        if (runstring(1:6) .eq. 'dphizj') notag=1
	
c-----------------------------------------------------------------------
      
      elseif (nproc .eq. 46) then
c-- 46 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))+f(p5)+f(p6)'
        case='Z_2jet'
        ndim=10
        n2=0
        n3=1
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='na'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        q1=0d0
        l1=ln*dsqrt(3d0)
        r1=rn*dsqrt(3d0)
        nwz=0   
        mass3=zmass
        width3=zwidth
       
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brznn
        endif

      elseif (nproc .eq. 47) then
c-- 47 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))+f(p5)+f(p6)+f(p7)'
        case='Z_3jet'
        ndim=13
        n2=0
        n3=1
        nqcdjets=3
        plabel(3)='nl'
        plabel(4)='na'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        q1=0d0
        l1=ln*dsqrt(3d0)
        r1=rn*dsqrt(3d0)
        nwz=0   
        mass3=zmass
        width3=zwidth
       
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brznn
        endif

c-----------------------------------------------------------------------
          
      elseif (nproc .eq. 50) then
c-- 50 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b~(p5)+b(p6) (massive)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +b~(p5)+b(p6) (massive)' (removebr=.true.)
        case='Zbbmas'
        call checkminzmass(1)
        write(6,*) 'mb=',mb
        bbproc=.true.
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        q1=-1d0
        l1=le
        r1=re
        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
        endif

      elseif ((nproc .ge. 51) .and. (nproc .le. 53)
     .   .or. (nproc .ge. 56) .and. (nproc .le. 58)) then
        case='Zbbbar'
        call checkminzmass(1)
        bbproc=.true.
        mb=0d0
        nqcdjets=2
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        q1=-1d0
        l1=le
        r1=re
        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        if     (nproc .le. 53) then
          flav=5
          nflav=4
        else
          flav=4
          nflav=3
        endif

        if     ((nproc .eq. 51) .or. (nproc .eq. 56)) then
c-- 51 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+b~(p6)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +b(p5)+b~(p6)' (removebr=.true.)
c-- 56 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)+c~(p6)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +c(p5)+c~(p6)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1d0
          l1=le
          r1=re
c--- total cross-section             
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee
          endif
        elseif (nproc .eq. 52) then
c-- 52 '  f(p1)+f(p2) --> Z_0(-->3*(nu(p3)+nu~(p4)))+b(p5)+b~(p6)'
          plabel(3)='nl'
          plabel(4)='na'
          q1=0d0
          l1=ln*dsqrt(3d0)
          r1=rn*dsqrt(3d0)
        elseif (nproc .eq. 53) then
c-- 53 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))+b(p5)+b~(p6)'
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        endif
        
      elseif (nproc .eq. 54) then
c-- 54 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+b~(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +b(p5)+b~(p6)+f(p7)' (removebr=.true.)
        case='Zbbjet'
        ndim=13
        bbproc=.true.
        mb=0d0
        nqcdjets=3
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        q1=-1d0
        l1=le
        r1=re
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        flav=5
        nflav=4

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
        endif

c-----------------------------------------------------------------------
          
      elseif (nproc/10 .eq. 6) then
        case='WWqqbr'
        call readcoup
        nqcdjets=0  
        plabel(7)='pp'
        nwz=1
        ndim=10
        mb=0d0
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

c--- include singly resonant diagrams if zerowidth=.false. , but only
c---  as long as anomtgc=.false. too
        srdiags=((zerowidth .eqv. .false.)
     &     .and. ( anomtgc  .eqv. .false.))
      
c--- zero srdiags for CDFdijet calculation 
c        if(runstring(1:10).eq.'cdf_Wdijet') then 
c           srdiags=.false.
c        endif

        if    ((nproc .eq. 61) .or. (nproc .eq. 69)) then
c--  61 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6))'
c--     '  f(p1)+f(p2) --> W^+ + W^- (for total Xsect)' (removebr=.true.)
          if (nproc .eq. 69) then
c--  69 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6)) [no pol]'
            case='WWnpol'
          endif
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='el'
          plabel(6)='na'
          l1=1d0
c--- total cross-section             
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen**2
          endif
        elseif (nproc .eq. 62) then
c--  62 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->q(p5)+q~(p6))'
c--- note: scattering diagrams are NOT included, only couplings change
          case='WWqqbr'
          nqcdjets=2
          notag=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='qj'
          plabel(6)='qj'
          plabel(7)='pp'
          l1=dsqrt(xn*2d0)

        elseif (nproc .eq. 63) then 
c--  63 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->q(p5)+q~(p6)) [rad.in.dk]'
          case='WWqqdk'
          nqcdjets=2
          notag=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='qj'
          plabel(6)='qj'
          plabel(7)='pp'
          l1=dsqrt(xn*2d0)

        elseif (nproc .eq. 64) then
c--  64 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+W^+(--> q(p5)+ q~(p6))'
c--- note: scattering diagrams are NOT included, only couplings change
          case='WWqqbr'
          nqcdjets=2
          notag=2
          plabel(5)='qj'
          plabel(6)='qj'
          plabel(3)='el'
          plabel(4)='na'
          plabel(7)='pp'
          l1=dsqrt(xn*2d0)
        elseif (nproc .eq. 65) then
c--  65 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+W^+(--> q(p5)+ q~(p6)),[rad.in.dk]'
c--- note: scattering diagrams are NOT included, only couplings change
          case='WWqqdk'
          nqcdjets=2
          notag=2
          plabel(5)='qj'
          plabel(6)='qj'
          plabel(3)='el'
          plabel(4)='na'
          plabel(7)='pp'
          l1=dsqrt(xn*2d0)
        elseif (nproc .eq. 66) then
c--  66 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+W^-(-->e^-(p5)+nu~(p6))+f(p7)'
c--     '  f(p1)+f(p2) --> W^+ + W^- + f(p7) (for total Xsect)' (removebr=.true.)
          case='WW_jet'
          nflav=4
	  nqcdjets=1
	  ndim=13
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='el'
          plabel(6)='na'
          plabel(7)='pp'
          plabel(8)='pp'
          l1=1d0
c--- total cross-section             
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen**2
          endif
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .ge. 70) .and. (nproc .le. 80)) then
        case='WZbbar'
        call checkminzmass(2)
        call readcoup
        plabel(7)='pp'
        ndim=10
        mb=0
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=wmass
        width3=wwidth

c--- include singly resonant diagrams if zerowidth=.false. , but only
c---  as long as anomtgc=.false. too
        srdiags=((zerowidth .eqv. .false.)
     &     .and. ( anomtgc  .eqv. .false.))
      
c--- Zero srdiags for CDFdijet calculation 
c        if(runstring(1:10).eq.'cdf_Wdijet') then 
c           srdiags=.false.
c        endif

        if (nproc .le. 75) then
c-- W^+Z
          nwz=+1

          if     (nproc .eq. 71) then             
c--  71 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->e^-(p5)+e^+(p6))'
c--     '  f(p1)+f(p2) --> W^+ (for total Xsect) + Z^0 ' (removebr=.true.)
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='ml'
            plabel(6)='ma'
            q1=-1d0
            l1=le
            r1=re
c--- total cross-section             
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              plabel(5)='ig'
              plabel(6)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brwen*brzee
            endif
          elseif (nproc .eq. 72) then
c--  72 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->nu_e(p5)+nu~_e(p6))'
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='nl'
            plabel(6)='na'
            plabel(7)='pp'
            q1=0d0
            l1=ln*dsqrt(3d0)
            r1=rn*dsqrt(3d0)
          elseif (nproc .eq. 73) then
c--  73 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->b(p5)+b~(p6))'
            bbproc=.true.
            nqcdjets=2
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='bq'
            plabel(6)='ba'
            plabel(7)='pp'
            q1=Q(5)*dsqrt(xn)
            l1=l(5)*dsqrt(xn)
            r1=r(5)*dsqrt(xn)
          elseif (nproc .eq. 74) then
c--  74 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->3*(d(p5)+d~(p6)))'
            nqcdjets=2
            notag=2
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='qj'
            plabel(6)='qj'
            plabel(7)='pp'
            q1=Q(5)*dsqrt(3d0*xn)
            l1=l(5)*dsqrt(3d0*xn)
            r1=r(5)*dsqrt(3d0*xn)
          elseif (nproc .eq. 75) then
c--  75 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->2*(u(p5)+u~(p6)))'
            nqcdjets=2
            notag=2 
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='qj'
            plabel(6)='qj'
            plabel(7)='pp'
            q1=Q(4)*dsqrt(2d0*xn)
            l1=l(4)*dsqrt(2d0*xn)
            r1=r(4)*dsqrt(2d0*xn)
          else
            call nprocinvalid()
          endif 

        elseif (nproc .ge. 76) then
c-- W^-Z
          nwz=-1

          if     (nproc .eq. 76) then
c--  76 '  f(p1)+f(p2) --> W^-(-->mu^-(p3)+nu~(p4))+Z^0(-->e^-(p5)+e^+(p6))'
c--     '  f(p1)+f(p2) --> W^- + Z^0 (for total Xsect)' (removebr=.true.)
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='ml'
            plabel(6)='ma'
            plabel(7)='pp'
            q1=-1d0
            l1=le
            r1=re
c--- total cross-section             
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              plabel(5)='ig'
              plabel(6)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brwen*brzee
            endif
          elseif (nproc .eq. 77) then
c--  77 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->nu(p5)+nu~(p6))'
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='nl'
            plabel(6)='na'
            plabel(7)='pp'
            q1=0d0
            l1=ln*dsqrt(3d0)
            r1=rn*dsqrt(3d0)
          elseif (nproc .eq. 78) then
c--  78 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->b(p5)+b~(p6))'
            bbproc=.true.
            nqcdjets=2
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='bq'
            plabel(6)='ba'
            plabel(7)='pp'
            q1=Q(5)*dsqrt(xn)
            l1=l(5)*dsqrt(xn)
            r1=r(5)*dsqrt(xn)
          elseif (nproc .eq. 79) then
c--  79 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->3*(d(p5)+d~(p6))'
            nqcdjets=0
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='qj'
            plabel(6)='qj'
            plabel(7)='pp'
            q1=Q(5)*dsqrt(3d0*xn)
            l1=l(5)*dsqrt(3d0*xn)
            r1=r(5)*dsqrt(3d0*xn)
          elseif (nproc .eq. 80) then
c--  80 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->2*(u(p5)+u~(p6)))'
            nqcdjets=0
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='qj'
            plabel(6)='qj' 
            plabel(7)='pp'
            q1=Q(4)*dsqrt(2d0*xn)
            l1=l(4)*dsqrt(2d0*xn)
            r1=r(4)*dsqrt(2d0*xn)
          else
            call nprocinvalid()
          endif 

        endif
            
c-----------------------------------------------------------------------

      elseif ((nproc .gt. 80) .and. (nproc .le. 90)) then
        vsymfact=1d0
        case='ZZlept'
        call checkminzmass(1)
        if ((nproc .eq. 81) .or. (nproc .eq. 83)) call checkminzmass(2)
        call readcoup
        plabel(7)='pp'
        nqcdjets=0
        nwz=0
        ndim=10
        mb=0
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth
        q1=-1d0
        l1=le
        r1=re
        interference=.false.

c--- only include singly-resonant diagrams when not in zerowidth approx.        
	if (zerowidth) then
	  srdiags=.false.
	else
	  srdiags=.true.
	endif
	
        if (nproc .eq. 81 .or. nproc .eq. 86) then
c--  81 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
c--     '  f(p1)+f(p2) --> Z^0 + Z^0 (for total Xsect)' (removebr=.true.)
c--  86 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->mu^-(p3)+mu^+(p4)) (NO GAMMA*)'
c--     '  f(p1)+f(p2) --> Z^0 + Z^0 (for total Xsect) (NO GAMMA*)' (removebr=.true.)
          plabel(3)='el' 
          plabel(4)='ea'
          plabel(5)='ml'
          plabel(6)='ma'
          q2=-1d0
          l2=le
          r2=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2d0*brzee**2  ! factor of 2 for identical particles
          endif
        elseif (nproc .eq. 82 .or. nproc .eq. 87) then
c--  82 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->3*(nu(p5)+nu~(p6)))'
c--  87 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->3*(nu(p5)+nu~(p6))) [no gamma^*]'
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='nl'
          plabel(6)='na'
          q2=0d0
          l2=ln*dsqrt(3d0)
          r2=rn*dsqrt(3d0)
        elseif (nproc .eq. 83 .or. nproc .eq. 88) then
c--  83 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->b(p5)+b~(p6))'
c--  88 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+Z^0(-->b(p5)+b~(p6)) [no gamma^*]'
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          q2=Q(5)*dsqrt(xn)
          l2=l(5)*dsqrt(xn)
          r2=r(5)*dsqrt(xn)
        elseif ((nproc .eq. 84) .or. (nproc .eq. 89))  then
c--  84 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + Z^0(-->3*(nu(p5)+nu~(p6)))'
c--  89 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + Z^0(-->3*(nu(p5)+nu~(p6))) [no gamma^*]'
          bbproc=.true.
          nqcdjets=2
          plabel(3)='bq'
          plabel(4)='ba'
          plabel(5)='nl'
          plabel(6)='na'
          q2=0d0
          l2=ln*dsqrt(3d0)
          r2=rn*dsqrt(3d0)
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
	elseif (nproc .eq. 85) then
c---  85 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->3*(nu(p5)+nu~(p6)))+f(p7)'
	  case='ZZ_jet'
	  nqcdjets=1
	  ndim=13
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='nl'
          plabel(6)='na'
          plabel(7)='pp'
          q2=0d0
          l2=ln*dsqrt(3d0)
          r2=rn*dsqrt(3d0)
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2d0*brzee**2  ! factor of 2 for identical particles
          endif
        elseif (nproc .eq. 90) then
c--  90 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
          interference=.true.
          vsymfact=0.25d0
          plabel(3)='el' 
          plabel(4)='ea'
          plabel(5)='el'
          plabel(6)='ea'
          q2=-1d0
          l2=le
          r2=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2d0*brzee**2  ! factor of 2 for identical particles
          endif
	  
        else
          call nprocinvalid()
        endif 

c-- remove gamma^* if necessary
        if ((nproc .gt. 85) .and. (nproc .lt. 90)) then
          q1=0d0
          q2=0d0
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 91) .or. (nproc .eq. 96)) then
        case='WHbbar'
        hdecaymode='bqba'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=2
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        
        ndim=10
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 91) then
c--  91 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) H(-->b(p5)+b~(p6)) '
c--     '  f(p1)+f(p2) --> W+ + H (for total Xsect)' (removebr=.true.)
          plabel(3)='nl'
          plabel(4)='ea'
          nwz=1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*br
            bbproc=.false.
            nqcdjets=0
          endif
        elseif (nproc .eq. 96) then
c--  96 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+ H(-->b(p5)+b~(p6))' 
c--     '  f(p1)+f(p2) --> W- + H (for total Xsect)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='na'
          nwz=-1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*br
            bbproc=.false.
            nqcdjets=0
          endif
        else
            call nprocinvalid()
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 92) .or. (nproc .eq. 97)) then
        case='WH__WW'
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=0
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='pp'
        
        ndim=16
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 92) then
c--- 92 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + H(-->W^+(nu(p3),e^+(p4))W^-(e^-(p5),nub(p6)))'
          plabel(3)='nl'
          plabel(4)='ea'
          nwz=+1
        elseif (nproc .eq. 97) then
          plabel(3)='el'
          plabel(4)='na'
          nwz=-1
        endif

        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'	       
          plabel(7)='ig'
          plabel(8)='ig'	       
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen**3*wwbr
        endif

c--- print warning if we're below threshold
        if (hmass .lt. 2d0*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
c        stop
	endif
	if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
                 
c-----------------------------------------------------------------------

        elseif     ((nproc .eq. 93) .or. (nproc .eq. 98)) then
C---93  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + H(-->Z(e^-(p5),e^+(p6))+Z(mu^-(p7),mu(p8)))' 
C---98  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + H(-->Z(e^-(p5),e^+(p6))+Z(mu^-(p7),mu(p8)))'
        case='WH__ZZ'
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=0
        plabel(5)='el'
        plabel(6)='ea'
        plabel(7)='el'
        plabel(8)='ea'


        ndim=16
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth

        l1=le
        r1=re
        l2=le
        r2=re

        if     (nproc .eq. 93) then
C---93  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + H(-->Z(e^-(p5),e^+(p6))+Z(mu^-(p7),mu(p8)))' 
          plabel(3)='nl'
          plabel(4)='ea'
          nwz=+1
        elseif (nproc .eq. 98) then
C--98  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + H(-->Z(e^-(p5),e^+(p6))+Z(mu^-(p7),mu(p8)))'
          plabel(3)='el'
          plabel(4)='na'
          nwz=-1
        endif

        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'	       
          plabel(7)='ig'
          plabel(8)='ig'	       
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=2d0*brwen*brzee**2*zzbr  ! factor of 2 for identical particles
        endif

c--- print warning if we're below threshold
        if (hmass .lt. 2d0*zmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
c        stop
	endif
	if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 94) .or. (nproc .eq. 99)) then
        case='WHgaga'
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=0
        hdecaymode='gaga'
        plabel(5)='ga'
        plabel(6)='ga'
        
        ndim=10
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 94) then
c '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + H(-->gamma(p5)+gamma(p6)' 'N'
          plabel(3)='nl'
          plabel(4)='ea'
          nwz=+1
        elseif (nproc .eq. 99) then
c '  f(p1)+f(p2) --> W^-(-->nu(p3)+e^+(p4)) + H(-->gamma(p5)+gamma(p6)' 'N'
          plabel(3)='el'
          plabel(4)='na'
          nwz=-1
        endif

        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'	       
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*gamgambr
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .ge. 101) .and. (nproc .le. 105)) then
        case='ZHbbar'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=2
        bbproc=.true.
        plabel(7)='pp'

        ndim=10
        nwz=0
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=zmass
        width3=zwidth

        if (nproc .eq. 101) then
c--  101 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->b(p5)+b~(p6))'
c--      '  f(p1)+f(p2) --> H + Z0 (for total Xsect)' (removebr=.true.)
          hdecaymode='bqba'
          call checkminzmass(1)
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          q1=-1d0
          l1=le
          r1=re
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee*br
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
            bbproc=.false.
            nqcdjets=0           
          endif
        elseif (nproc .eq. 102) then
c--  102 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + H(-->b(p5)+b~(p6))'
          hdecaymode='bqba'
          plabel(3)='nl'
          plabel(4)='na'
          plabel(5)='bq'
          plabel(6)='ba'
          q1=0d0
          l1=ln*dsqrt(3d0)
          r1=rn*dsqrt(3d0)
        elseif (nproc .eq. 103) then
c--  103 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + H(-->b(p5)+b~(p6))'     
          hdecaymode='bqba'
          call checkminzmass(1)
          nqcdjets=4
          plabel(3)='bq'
          plabel(4)='ba'
          plabel(5)='bq'
          plabel(6)='ba'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        elseif (nproc .eq. 104) then
c--  104 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->gamma(p5)+gamma(p6))' 'N'
          call checkminzmass(1)
          case='ZHgaga'
          nqcdjets=0
          bbproc=.false.
          q1=-1d0
          l1=le
          r1=re
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='ga'
          plabel(6)='ga'
          hdecaymode='gaga'
           if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee*gamgambr          
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
          endif
        elseif (nproc .eq. 105) then
c--  105 '  f(p1)+f(p2) --> Z^0(-->-->3*(nu(p3)+nu~(p4))) + H(-->gamma(p5)+gamma(p6))' 'N'
          case='ZHgaga'
          nqcdjets=0
          bbproc=.false.
          plabel(3)='nl'
          plabel(4)='na'
          plabel(5)='ga'
          plabel(6)='ga'          
          hdecaymode='gaga'
          q1=0d0
          l1=ln*dsqrt(3d0)
          r1=rn*dsqrt(3d0)
           if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)           
            BrnRat=brznn*gamgambr
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
           endif
        else
          call nprocinvalid()
        endif 

c-----------------------------------------------------------------------

      elseif ((nproc .ge. 106) .and. (nproc .le. 108)) then
        case='ZH__WW'
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=0
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='pp'
        
        ndim=16
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=zmass
        width3=zwidth

        if (nproc .eq. 106) then
c--  106 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8)))'
c--      '  f(p1)+f(p2) --> H + Z0 (for total Xsect)' (removebr=.true.)
          call checkminzmass(1)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1d0
          l1=le
          r1=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
            plabel(7)='ig'
            plabel(8)='ig'	       
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee*brwen**2*wwbr
          endif
        elseif (nproc .eq. 107) then
c--  107 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8)))'
          plabel(3)='nl'
          plabel(4)='na'
          q1=0d0
          l1=ln*dsqrt(3d0)
          r1=rn*dsqrt(3d0)
        elseif (nproc .eq. 108) then
c--  108 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8)))'     
          call checkminzmass(1)
          nqcdjets=2
          plabel(3)='bq'
          plabel(4)='ba'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        else
          call nprocinvalid()
        endif 

c--- print warning if we're below threshold
        if (hmass .lt. 2d0*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
                 
c-----------------------------------------------------------------------

        elseif (nproc .eq. 109) then
c--  109 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->Z(e^-(p5),e^+(p6))+Z(mu^-(p7),mu(p8)))'
        case='ZH__ZZ'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=0
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='ea'
        plabel(7)='el'
        plabel(8)='ea'
        plabel(9)='pp'
        
        ndim=16
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=zmass
        width3=zwidth
          q1=-1d0
          l1=le
          r1=re
          l2=le
          r2=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
            plabel(7)='ig'
            plabel(8)='ig'	       
            plabel(9)='ig'	       
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2d0*brzee**3*zzbr
          endif
c--- print warning if we're below threshold
        if (hmass .lt. 2d0*zmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
                 
c-----------------------------------------------------------------------

       elseif ((nproc .eq. 111) .or. (nproc .eq. 112)) then
        case='ggfus0'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(5)='pp'
        ndim=4
      
        n2=0
        n3=1
        mass3=hmass
        width3=hwidth

        if     (nproc .eq. 111) then
c--  111 '  f(p1)+f(p2) --> H(-->b(p3)+bbar(p4))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)       
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=0
            BrnRat=br
          endif
          
        elseif (nproc .eq. 112) then
c--  112 '  f(p1)+f(p2) --> H(-->tau^-(p3)+tau^+(p4))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)       
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=0
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            BrnRat=tautaubr
          endif
        endif

      elseif ((nproc .eq. 113) .or. (nproc .eq. 126)
     &   .or. (nproc .eq. 127)) then
c--  113 '  f(p1)+f(p2) --> H (--> W^+(nu(p3)+e^+(p4)) + W^-(e^-(p5)+nu~(p6)))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
        if     (nproc .eq. 113) then
          case='HWW_4l'  
	elseif (nproc .eq. 126) then
          case='HWW_tb'  
	elseif (nproc .eq. 127) then
	  case='HWWint'
	endif           
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
c--- widths according to Kauer et al., for comparison with gg2WW
c        if (abs(hmass-140d0) .lt. 1d-4) hwidth=0.008235d0
c        if (abs(hmass-170d0) .lt. 1d-4) hwidth=0.3837d0
c        if (abs(hmass-200d0) .lt. 1d-4) hwidth=1.426d0
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='pp'
        nqcdjets=0
        ndim=10
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

c--- print warning if we're below threshold
        if (hmass .lt. 2d0*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen**2*wwbr
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'               
        endif
        
      elseif (nproc .eq. 114) then
      case='HWW2lq'  
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
c--- widths according to Kauer et al., for comparison with gg2WW
c        if (abs(hmass-140d0) .lt. 1d-4) hwidth=0.008235d0
c        if (abs(hmass-170d0) .lt. 1d-4) hwidth=0.3837d0
c        if (abs(hmass-200d0) .lt. 1d-4) hwidth=1.426d0
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        nqcdjets=2
        ndim=10
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

c--- print warning if we're below threshold
        if (hmass .lt. 2d0*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=2d0*xn*brwen**2*wwbr
c	  if (part .eq. 'todk') BrnRat=BrnRat*(1d0+as/pi)
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'        
	  nqcdjets=0       
        endif
        
      elseif (nproc .eq. 115) then
      case='HWWdkW'  
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        nqcdjets=2
        ndim=10
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=2d0*xn*brwen**2*wwbr
c	  if (part .eq. 'todk') BrnRat=BrnRat*(1d0+as/pi)
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'   
	  nqcdjets=0            
        endif
        
c--- print warning if we're below threshold
        if (hmass .lt. 2d0*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
        
      elseif ((nproc .ge. 116) .and. (nproc .le. 118)) then
        case='HZZ_4l'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(7)='pp'
        nqcdjets=0
        ndim=10
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth

c--- print warning if we're below threshold
        if (hmass .lt. 2d0*zmass) then
          write(6,*)
          write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
          write(6,*) 'may not yield sensible results - check the number'
          write(6,*) 'of integration points'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
c        stop
	endif
          if (zerowidth) then
          write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c          stop
          endif
        endif
        
        if     (nproc .eq. 116) then
c--  116 '  f(p1)+f(p2) --> H(-->Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='ml'
          plabel(6)='ma'
          l1=le
          r1=re
          l2=le
          r2=re
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2d0*brzee**2*zzbr  ! factor of 2 for identical particles
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'             
          endif
        elseif (nproc .eq. 117) then
c--  117 '  f(p1)+f(p2) --> H(-->Z^0(3*(nu(p3)+nu~(p4)))+ Z^0(e^-(p5)+e^+(p6))'
          plabel(3)='nl'
          plabel(4)='na'
          plabel(5)='ml'
          plabel(6)='ma'
          l1=le
          r1=re
          l2=ln*dsqrt(3d0)
          r2=rn*dsqrt(3d0)      
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2d0*brzee*brznn*zzbr  ! factor of 2 for identical particles
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'             
          endif
        elseif (nproc .eq. 118) then
c--  118 '  f(p1)+f(p2) --> H(-->Z^0(mu^-(p3)+mu^+(p4)) + Z^0(b(p5)+b~(p6))'
          nqcdjets=2
          plabel(3)='ml'
          plabel(4)='ma'
          plabel(5)='bq'
          plabel(6)='ba'
          l1=le
          r1=re
          l2=l(5)*dsqrt(xn)
          r2=r(5)*dsqrt(xn)
        else
          call nprocinvalid()
        endif

        elseif (nproc .eq. 119) then
        case='Higaga'
        plabel(5)='pp'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        ndim=4
      
        n2=0
        n3=1
        mass3=hmass
        width3=hwidth

c--  119 '  f(p1)+f(p2) --> H(-->gamma^-(p3)+gamma^+(p4))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)       
        plabel(3)='ga'
        plabel(4)='ga'
        hdecaymode='gaga'
        nqcdjets=0
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          BrnRat=gamgambr
        endif

        elseif (nproc .eq. 120) then
        case='Hi_Zga'
        call checkminzmass(1)
        nqcdjets=0
        plabel(6)='pp'
        ndim=7
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
      
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

c--  120 '  f(p1)+f(p2) --> H(-->Z^0(mu^-(p3)+mu^+(p4)) + gamma(p5)')'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ga'
        l1=le
        r1=re
        q1=-1d0
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee*zgambr 
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
        endif

        elseif (nproc .eq. 121) then
        case='Hi_Zga'
        nqcdjets=0
        plabel(6)='pp'
        ndim=7
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
      
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

c--  121 '  f(p1)+f(p2) --> H(-->Z^0(3*(nu(p3)+nu~(p4))) + gamma(p5)')'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
        plabel(3)='nl'
        plabel(4)='na'
        q1=0d0
        l1=ln*dsqrt(3d0)
        r1=rn*dsqrt(3d0)      
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brznn*zgambr
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .ge. 131) .and. (nproc .le. 133)) then
        case='H_1jet'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        call setmb_msbar
        ndim=7
        plabel(3)='qb'
        plabel(4)='ab'
        plabel(5)='bq'
        hdecaymode='bqba'
        nqcdjets=1

        n2=0
        n3=1
        mass3=hmass
        width3=hwidth
        
        if ( (nproc .eq. 132) .and.  
     .       ((part .eq. 'virt') .or. (part .eq. 'tota')) ) then
          write(6,*) 'This process number is not suitable for the'
          write(6,*) 'NLO calculation. Please run processes'
          write(6,*) '131 (virtual+real) and 132 (real) separately.'
c          stop
        endif
        if ( (nproc .eq. 133) .and. (part .ne. 'real') ) then
          write(6,*) 'This process number is not suitable for such a'
          write(6,*) 'calculation. Please run process 133 (real) only.'
c          stop
        endif
             
        if     (nproc .eq. 131) then
c--  131 '  f(p1)+f(p2) --> H (no BR) + b(p5) [+g(p6)]'
          isub=1
          plabel(6)='pp'
        elseif (nproc .eq. 132) then
c--  132 '(p1)+f(p2) --> H (no BR) + b~(p5) [+b(p6)]'
          isub=2
          plabel(5)='ba'
          plabel(6)='bq'
        elseif (nproc .eq. 133) then
c--  133 '  f(p1)+f(p2) --> H (no BR) + b(p5) + b~(p6) [both observed]'
          isub=2
          plabel(5)='ba'
          plabel(6)='bq'
          nqcdjets=2
        endif
        
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          BrnRat=br
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 141) 
     &   .or. (nproc .eq. 142) 
     &   .or. (nproc .eq. 145) 
     &   .or. (nproc .eq. 146)
     &   .or. (nproc .eq. 147)
     &   .or. (nproc .eq. 148)
     &   .or. (nproc .eq. 149)
     &   .or. (nproc .eq. 150)
     &   .or. (nproc .eq. 151)) then
        ndim=16
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth
        bbproc=.true.
        
        if (nproc .eq. 141) then
c--  141 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+e^-(p7)+nu~(p8)'
c--      '  f(p1)+f(p2) --> t t~ (with BR for total Xsect)' (removebr=.true.)
          case='tt_bbl'
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brwen*brtop)**2
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            plabel(7)='ig'
            plabel(8)='ig'
            nqcdjets=0
            bbproc=.false.
          endif
        elseif (nproc .eq. 142) then
c--  142 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+e^-(p7)+nu~(p8) [radiation in top decay]'
          case='tt_ldk'
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brwen*brtop)**2
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            plabel(7)='ig'
            plabel(8)='ig'
            nqcdjets=0
            bbproc=.false.
          endif
        elseif (nproc .eq. 145) then
c--  145 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->b~(p6)+e^-(p7)+nu~(p8))'
c--           (uncorrelated)'
          case='tt_bbu'
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brwen*brtop)**2
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            plabel(7)='ig'
            plabel(8)='ig'
            nqcdjets=0
            bbproc=.false.
          endif
        elseif (nproc .eq. 146) then
c--  146 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+q(p7)+q~(p8)'
          case='tt_bbh'
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='pp'
          nqcdjets=4
        elseif (nproc .eq. 147) then
c--  147 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+q(p7)+q~(p8) [radiation in top decay]'
          case='tt_hdk'
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='pp'
          nqcdjets=4
        elseif (nproc .eq. 148) then
c--  148 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+q(p7)+q~(p8) [radiation in hadronic W decay]'
          case='tthWdk'
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='pp'
          nqcdjets=4
        elseif (nproc .eq. 149) then
c---  149 '  f(p1)+f(p2) --> t(-->q(p3)+q~(p4)+b(p5))+t~(-->b~(p6)+e-(p7)+nu~(p8))'
          case='tt_bbh'
          plabel(3)='pp'
          plabel(4)='pp'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          nqcdjets=4
        elseif (nproc .eq. 150) then
c---  150 '  f(p1)+f(p2) --> t(-->q(p3)+q~(p4)+b(p5))+t~(-->b~(p6)+e-(p7)+nu~(p8)) [radiation in top decay]'
          case='tt_hdk'
          plabel(3)='pp'
          plabel(4)='pp'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          nqcdjets=4
        elseif (nproc .eq. 151) then
c---  151 '  f(p1)+f(p2) --> t(-->q(p3)+q~(p4)+b(p5))+t~(-->b~(p6)+e-(p7)+nu~(p8)) [radiation in hadronic W decay]'
          case='tthWdk'
          plabel(3)='pp'
          plabel(4)='pp'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          nqcdjets=4
        endif 
c-----------------------------------------------------------------------

      elseif (nproc .eq. 143) then
c--  143 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->nu~(p7)+e^-(p8)+b~(p6))+g(p9)'
c--      '  f(p1)+f(p2)-->t(p345)+t~(p678)+g(p9)' (removebr=.true.)
        case='qq_ttg'
        nwz=1
        ndim=19
        nqcdjets=3
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='pp'

c--- total cross-section             
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=(brwen*brtop)**2
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'               
          plabel(7)='ig'
          plabel(8)='ig'
          nqcdjets=1
        endif

c-----------------------------------------------------------------------

      elseif (nproc .eq. 157) then
c--  157 '  f(p1)+f(p2) --> t t~ (for total Xsect)'
        case='tt_tot'
        nqcdjets=0
        ndim=4
        mass2=mt
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
      elseif (nproc .eq. 158) then
c--  158 '  f(p1)+f(p2) --> b b~ (for total Xsect)'
        case='bb_tot'
	nflav=4
        nqcdjets=0
        ndim=4
        mass2=mb
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
      elseif (nproc .eq. 159) then
c--  159 '  f(p1)+f(p2) --> c c~ (for total Xsect)'
        case='cc_tot'
	nflav=3
        nqcdjets=0
        ndim=4
        mass2=mc
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
 
      elseif (nproc .eq. 160) then
      if  ((part .eq. 'tota')
     . .or.(part .eq. 'virt')
     . .or.(part .eq. 'real')) then
          write(6,*) 'This process number is available only at LO'
          write(6,*) 'Please set part = lord and rerun'
c             stop
      endif
c--  160 '  f(p1)+f(p2) --> t t~ +jet (for total Xsect)'
        case='tt_glu'
        nqcdjets=1
        ndim=7
        mass2=mt
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 161) .or. (nproc .eq. 163)) then
c--  161 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+q(p6) [t-channel]'
c--      '  f(p1)+f(p2) --> t(no BR) + q(p6)' (removebr=.true.)
        case='bq_tpq'
        isub=1
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=+1
	
	if (nproc .eq. 161) then ! usual approach, mb=0 
c--- extra b that can appear at NLO is massless
          masslessb=.true.	
	else                     ! proper ACOT, mb>0 (must run 231 LO)
	  masslessb=.false.
	endif

c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0d0
	write(6,*) 'Enforcing mb=0 for this process!'
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1   
c--- the default is now that the parameter notag=1, so that
c--- the calculation is inclusive of all additional jets; it may
c--- be commented out to explicitly require an additional jet at LO
	  notag=1 
        endif

      elseif (nproc .eq. 162) then
c--  162 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+q(p6) [decay]'
        case='ttdkay'
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=+1

        if (part .eq. 'lord') then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '161 (lord) or process 162 (virt+real).'
c          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0d0
	write(6,*) 'Enforcing mb=0 for this process!'
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
c--- the default is now that the parameter notag=1, so that
c--- the calculation is inclusive of all additional jets; it may
c--- be commented out to explicitly require an additional jet at LO
	  notag=1 
        endif

      elseif ((nproc .eq. 166) .or. (nproc .eq. 168)) then
c--  166 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+q(p6) [t-channel]''
c--      '  f(p1)+f(p2) --> t~(no BR) + q(p6)' (removebr=.true.)
        case='bq_tpq'
        isub=1
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=-1

	if (nproc .eq. 166) then ! usual approach, mb=0 
c--- extra b that can appear at NLO is massless
          masslessb=.true.	
	else                     ! proper ACOT, mb>0 (must run 236 LO)
	  masslessb=.false.
	endif

c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0d0
	write(6,*) 'Enforcing mb=0 for this process!'
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
c--- the default is now that the parameter notag=1, so that
c--- the calculation is inclusive of all additional jets; it may
c--- be commented out to explicitly require an additional jet at LO
	  notag=1 
        endif
        
      elseif (nproc .eq. 167) then
c--  167 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+q(p6) [decay]'
        case='ttdkay'
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=-1

        if (part .eq. 'lord') then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '166 (lord) or process 167 (virt+real).'
c          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0d0
	write(6,*) 'Enforcing mb=0 for this process!'
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
c--- the default is now that the parameter notag=1, so that
c--- the calculation is inclusive of all additional jets; it may
c--- be commented out to explicitly require an additional jet at LO
	  notag=1 
        endif

c-----------------------------------------------------------------------

      elseif (nproc .eq. 171) then
c--  171 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)) [s-channel]'
c--      '  f(p1)+f(p2) --> t(no BR) + b~(p6)' (removebr=.true.)
        case='t_bbar'
        isub=2
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nwz=1
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif
        
      elseif (nproc .eq. 172) then
c--  172 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)) [decay]'
        case='tdecay'
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nwz=1
        
        if (part .eq. 'lord') then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '171 (lord) or process 172 (virt+real).'
c          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif
        
      elseif (nproc .eq. 176) then
c--  176 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+b(p6)) [s-channel]'
c--      '  f(p1)+f(p2) --> t~(no BR) + b(p6)' (removebr=.true.)
        case='t_bbar'
        isub=2
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='bq'
        plabel(7)='pp'
        nwz=-1
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
        
         if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif
             
      elseif (nproc .eq. 177) then
c--  177 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+b(p6)) [decay]'
        case='tdecay'
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='bq'
        plabel(7)='pp'
        nwz=-1
        
        if (part .eq. 'lord') then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '176 (lord) or process 177 (virt+real).'
c          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif
        
c-----------------------------------------------------------------------

      elseif (nproc .eq. 180) then
c--  180 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(p5)'
        case='W_tndk'
        nqcdjets=0
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ig'
        plabel(6)='pp'
        mass2=mt
        nflav=5
        nwz=-1

        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
             
      elseif (nproc .eq. 181) then
c--  181 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(nu(p5)+e^+(p6)+b(p7))'
        case='W_twdk'
        nqcdjets=1
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='bq'
        plabel(8)='pp'
        nflav=5
        nwz=-1

        ndim=13
        mb=0d0
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=0
        endif
             
      elseif (nproc .eq. 182) then
c--  182 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(nu(p5)+e^+(p6)+b(p7)) [decay]'
        
        case='Wtdkay'
        nqcdjets=1
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='bq'
        plabel(8)='pp'
        nflav=5
        nwz=-1

        if (part .eq. 'lord') then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '181 (lord) or process 182 (virt+real).'
c          stop
        endif
        
        ndim=13
        mb=0d0
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=0
        endif
             
      elseif (nproc .eq. 183) then
c--  183 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(nu(p5)+e^+(p6)+b(p7))+b(p8)'
        case='Wtbwdk'
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='bq'
        plabel(8)='pp'
        nflav=5
        nwz=-1

        ndim=16
c--- (this process can also be used for non-zero mb)
        mb=0d0
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=1
        endif
             
      elseif (nproc .eq. 184) then
c--  184 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(p5)+b(p6) [massive b]'
        case='Wtbndk'
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='bq'
        plabel(6)='ba'
        mass2=mt
        nflav=5
        nwz=-1

        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
             
      elseif (nproc .eq. 185) then
c--  185 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+tbar(p5)'
        case='W_tndk'
        nqcdjets=0
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='ig'
        plabel(6)='pp'
        mass2=mt
        nflav=5
        nwz=+1
        
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
             
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif (nproc .eq. 186) then
c--  186 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+t~(e^-(p5)+nu~(p6)+bbar(p7))'
        case='W_twdk'
        nqcdjets=1
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='ba'
        plabel(8)='pp'
        nflav=5
        nwz=+1
        
        ndim=13
        mb=0
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=0           
        endif
             
      elseif (nproc .eq. 187) then
c--  182 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+t~(e^-(p5)+nu~(p6)+bbar(p7)) [decay]'
        
        case='Wtdkay'
        nqcdjets=1
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='bq'
        plabel(8)='pp'
        nflav=5
        nwz=+1

        if (part .eq. 'lord') then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '186 (lord) or process 187 (virt+real).'
c          stop
        endif
        
        ndim=13
        mb=0
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=0           
        endif
             
c-----------------------------------------------------------------------

      elseif (nproc .eq. 190) then
c--  190 '  f(p1)+f(p2)-->t(p3)+t~(p4)+H(p5)'
        case='tottth'
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        nwz=1
        n2=0
        n3=0
        ndim=7
        
      elseif (nproc .eq. 191) then
c--  191 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))
c--         +t~(-->nu~(p7)+e^-(p8)+b~(p6))+H(p9+p10)'
c--      '  f(p1)+f(p2)-->t(p3+p4+p5)+t~(p6+p7+p8)+H(p9+p10)' (removebr=.true.)
        case='qq_tth'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='bq'
        plabel(10)='ba'
        hdecaymode='bqba'
	nqcdjets=4

        nwz=1
        ndim=22
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=(brtop*brwen)**2*br
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          plabel(8)='ig'
          plabel(9)='ig'
          plabel(10)='ig'
	  nqcdjets=0
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .ge. 200) .and. (nproc .le. 210)) then
        case='httjet'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=1
        plabel(5)='pp'
 
        ndim=7
        n2=0
        n3=1
        mass3=hmass
        width3=hwidth
        
        if     (nproc .eq. 201) then
c--  201 '  f(p1)+f(p2)--> H(-->b(p3)+b~(p4)) + f(p5) [full mt dep.]'
c--      '  f(p1)+f(p2)--> H(p3+p4) + f(p5) (for total Xsect)' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'  
          nqcdjets=3
          if (removebr) then        
            BrnRat=br
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
          endif

        elseif (nproc .eq. 202) then
c--  202 '  f(p1)+f(p2)--> H (-> tau(p3) tau~(p4)) + f(p5) [full mt dep.]'
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          if (removebr) then        
            BrnRat=tautaubr
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
          endif

        elseif ((nproc .eq. 203) .or. (nproc .eq. 204)) then
          case='ggfus1'
          nqcdjets=1
          mb=0
          plabel(5)='pp'
          plabel(6)='pp'
          ndim=7
      
          n2=0
          n3=1

          if     (nproc .eq. 203) then
c--  203 '  f(p1)+f(p2) -->H(-->b(p3)+b~(p4)) + f(p5)'
c--      '  f(p1)+f(p2)--> H(p3+p4) + f(p5) (for total Xsect)' (removebr=.true.)
            hdecaymode='bqba'
            plabel(3)='bq'
            plabel(4)='ba'  
            nqcdjets=3
            if (removebr) then        
              BrnRat=br
              plabel(3)='ig'
              plabel(4)='ig'
              nqcdjets=1
            endif
          elseif (nproc .eq. 204) then
c--  204 '  f(p1)+f(p2) -->H(-->tau^-(p3)+tau^+(p4)) + f(p5)'
            hdecaymode='tlta'
            plabel(3)='tl'
            plabel(4)='ta'
            nqcdjets=1
            if (removebr) then        
              plabel(3)='ig'
              plabel(4)='ig'
              Brnrat=tautaubr
            endif
          endif
        
        elseif (nproc .eq. 206) then
c--  206 '  f(p1)+f(p2)--> A(-->b(p3)+b~(p4)) + f(p5) [full mt dep.]'
c--      '  f(p1)+f(p2)--> A(p3+p4) + f(p5) (for total Xsect)' (removebr=.true.)
          case='attjet'
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'  
          nqcdjets=3
          if (removebr) then        
            BrnRat=br
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
          endif

        elseif (nproc .eq. 207) then
c--  207 '  f(p1)+f(p2)--> A (--> tau(p3) tau~(p4)) + f(p5) [full mt dep.]'
         case='attjet'
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          if (removebr) then        
            BrnRat=tautaubr
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
          endif
        endif

        if     (nproc .eq. 208) then
c-- 208 '  f(p1)+f(p2) --> H(-->W^+(p3,p4)W^-(p5,p6)) + f(p7)'
          case='HWWjet'
          ndim=13
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='el'
          plabel(6)='na'
          plabel(7)='pp'
          plabel(8)='pp'
          nqcdjets=1
          n2=1
          n3=1
          mass2=wmass
          width2=wwidth
          mass3=wmass
          width3=wwidth
c--- print warning if we're below threshold
          if (hmass .lt. 2d0*wmass) then
          write(6,*)
          write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
          write(6,*) 'may not yield sensible results - check the number'
          write(6,*) 'of integration points and the value of zerowidth'
	  if (removebr) then
	  write(6,*)
	write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
c          stop
	  endif
          if (zerowidth) then
          write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c          stop
          endif
          endif
        
          if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=wwbr*brwen**2
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          endif
        endif

       if ((nproc .eq. 209)) then
          case='HZZjet'
          l1=le
          l2=le
          r1=re
          r2=re
          ndim=13
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='ml'
          plabel(6)='ma'
          plabel(7)='pp'
          plabel(8)='pp'
          nqcdjets=1
          n2=1
          n3=1
          mass2=zmass
          width2=zwidth
          mass3=zmass
          width3=zwidth
          call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

c--- print warning if we're below threshold
        if (hmass .lt. 2d0*zmass) then
          write(6,*)
          write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
          write(6,*) 'may not yield sensible results - check the number'
          write(6,*) 'of integration points'
          if (zerowidth) then
          write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c          stop
          endif
        endif
        
          l1=le
          r1=re
          l2=le
          r2=re
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2d0*brzee**2*zzbr  ! factor of 2 for identical particles
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'             
          endif

        elseif     (nproc .eq. 210) then
c--  210 '  f(p1)+f(p2) -->H(-->gamma(p3)+gamma(p4)) + f(p5)'
            case='Hgagaj'
            plabel(3)='ga'
            plabel(4)='ga'  
            plabel(5)='pp'  
            plabel(6)='pp'  
            hdecaymode='gaga'
          nqcdjets=1
            if (removebr) then        
              BrnRat=gamgambr
              plabel(3)='ig'
              plabel(4)='ig'
            endif
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 211) .or. (nproc .eq. 212)) then
        case='qq_Hqq'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nwz=2
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1
	
        mass3=hmass
        width3=hwidth

        if     (nproc .eq. 211) then
c--  211 '  f(p1)+f(p2)--> H(-->b(p3)+b~(p4))+f(p5)+f(p6) [WBF]'
c--      '  f(p1)+f(p2)--> H(p3+p4)+f(p5)+f(p6) [WBF]' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=4
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=2
            BrnRat=br
          endif
        elseif (nproc .eq. 212) then
c--  212 '  f(p1)+f(p2)--> H(-->tau-(p3)+tau+(p4))+f(p5)+f(p6) [WBF]'
c--      '  f(p1)+f(p2)--> H(p3+p4)+f(p5)+f(p6) [WBF]' (removebr=.true.)
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            Brnrat=tautaubr
          endif
        endif
          
      elseif (nproc .eq. 213) then
        case='qq_HWW'
        mb=0d0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nwz=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        ndim=16
        nqcdjets=2
c	notag=1 ! If only one jet is required
c	notag=0 ! FOR CHECKING VS 211
c        nqcdjets=3 ! DEBUG
	
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth
c--- print warning if we're below threshold
        if (hmass .lt. 2d0*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=wwbr*brwen**2
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif
          
      elseif (nproc .eq. 214) then
        case='qq_HZZ'
        l1=le
        r1=re
        l2=le
        r2=re
        mb=0d0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nwz=2
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='ea'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        ndim=16
        nqcdjets=2
	
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth
c--- print warning if we're below threshold
        if (hmass .lt. 2d0*zmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=2d0*zzbr*brzee**2  ! factor of 2 for identical particles
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif
          
      elseif (nproc .eq. 215) then
        case='qq_Hgg'
	hdecaymode='gaga'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='ga'
        plabel(4)='ga'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
        nqcdjets=2
	
        n2=0
        n3=1
        mass3=hmass
        width3=hwidth
       
        
	if (removebr) then
	BrnRat=gamgambr
	plabel(3)='ig'
	plabel(4)='ig'
	endif
          
      elseif ((nproc .eq. 216) .or. (nproc .eq. 217)) then
        case='qqHqqg'
        mb=0d0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nwz=2
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=13
        n2=0
        n3=1

        mass3=hmass
        width3=hwidth

        if     (nproc .eq. 216) then
c-- 216 '  f(p1)+f(p2)--> H(-->b(p3)+b~(p4))+f(p5)+f(p6)+f(p7) [WBF+jet]'
c--     '  f(p1)+f(p2)--> H(p3+p4)+f(p5)+f(p6)+f(p7) [WBF+jet]' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=5
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=3
            BrnRat=br
          endif
        elseif (nproc .eq. 217) then
c-- 217 '  f(p1)+f(p2)--> H(-->tau-(p3)+tau+(p4))+f(p5)+f(p6)+f(p7) [WBF+jet]'
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=3
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
            Brnrat=tautaubr
          endif
        endif

c-----------------------------------------------------------------------

      elseif (nproc .eq. 221) then
        case='tautau'
c--  221 '  f(p1)+f(p2)--> tau^-(-->e^-(p3)+nu~_e(p4)+nu_tau(p5))+tau^+(-->nu~_tau(p6)+nu_e(p7)+e^+(p8))'
c--      '  f(p1)+f(p2)--> tau tau~ [for total Xsect]' (removebr=.true.)
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='nl'
        plabel(6)='na'
        plabel(7)='nl'
        plabel(8)='ea'
        nqcdjets=0
        nwz=1
        ndim=16
        n2=1
        n3=1
        mass2=mtau
        width2=tauwidth
        mass3=mtau
        width3=tauwidth

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brtau**2
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          plabel(8)='ig'
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 231) .or. (nproc .eq. 236)) then
c--  231 '  f(p1)+f(p2) --> t(p3)+b~(p4)+q(p5) [t-channel]'
c--  236 '  f(p1)+f(p2) --> t~(p3)+b(p4)+q(p5) [t-channel]'
        case='qg_tbq'
        nqcdjets=1
c--- the default is now that the parameter notag=0, so that
c--- the calculation requires the presence of the light jet that
c--- is present at LO; to compute an inclusive cross section, one
c--- can set notag=1, or use jet cuts that have no effect
	notag=0
        ndim=7
        mass2=mt
        mass3=mb
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
        plabel(6)='pp'
	if (nproc .eq. 231) then
	  nwz=+1
	else
	  nwz=-1
	endif

c---  in the SM, the logical fourthgen should be false
c---  for BSM calculations, it should be true and it indicates that
c---   5 flavours should be used in the PDF and alpha-s
        fourthgen=.false.

	if (fourthgen) then
c--- BSM: full 5 light flavours
          nflav=5
	  bmass=4.7d0  !  set b-mass to its usual value
	else
c--- SM: only 4 light flavours
          nflav=4
          bmass=1001d0 !  enforce 4-flavour running in alfamz.f
        endif
	
c--- set up correct scales and as on heavy and light quark lines
        facscale_H=initfacscale_H
        facscale_L=initfacscale_L
        renscale_H=initrenscale_H
        renscale_L=initrenscale_L
c--- make sure it works even if not specifying separate scales
        if (initrenscale_L .eq. 0d0) then 
	  facscale_H=facscale
	  facscale_L=facscale
	  renscale_H=scale
	  renscale_L=scale
	endif
	
        b0=(xn*11d0-2d0*nflav)/6d0
	as_H=alphas(abs(renscale_H),amz,nlooprun)
	as_L=alphas(abs(renscale_L),amz,nlooprun)

	
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 232) .or. (nproc .eq. 237)) then
c--  232 '  f(p1)+f(p2) --> t(p3)+b~(p4)+q(p5)+q(p6) [t-channel]'
c--  237 '  f(p1)+f(p2) --> t~(p3)+b(p4)+q(p5)+q(p6) [t-channel]'
        case='qgtbqq'
        nqcdjets=2
c--- the default is now that the parameter notag=0, so that
c--- the calculation requires the presence of two light jets that
c--- are present; to compute an inclusive cross section, one
c--- can set notag=1, or use jet cuts that have no effect
	notag=0
        ndim=10
        mass2=mt
        mass3=mb
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
        plabel(6)='pp'
	if (nproc .eq. 232) then
	  nwz=+1
	else
	  nwz=-1
	endif

c---  in the SM, the logical fourthgen should be false
c---  for BSM calculations, it should be true and it indicates that
c---   5 flavours should be used in the PDF and alpha-s
        fourthgen=.false.

	if (fourthgen) then
c--- BSM: full 5 light flavours
          nflav=5
	  bmass=4.7d0  !  set b-mass to its usual value
	else
c--- SM: only 4 light flavours
          nflav=4
          bmass=1001d0 !  enforce 4-flavour running in alfamz.f
        endif
	
c--- set up correct scales and as on heavy and light quark lines
        facscale_H=initfacscale_H
        facscale_L=initfacscale_L
        renscale_H=initrenscale_H
        renscale_L=initrenscale_L
c--- make sure it works even if not specifying separate scales
        if (initrenscale_L .eq. 0d0) then 
	  facscale_H=facscale
	  facscale_L=facscale
	  renscale_H=scale
	  renscale_L=scale
	endif
	
        b0=(xn*11d0-2d0*nflav)/6d0
	as_H=alphas(abs(renscale_H),amz,nlooprun)
	as_L=alphas(abs(renscale_L),amz,nlooprun)

	
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 233) .or. (nproc .eq. 238)
     &   .or. (nproc .eq. 234) .or. (nproc .eq. 239)) then
c--  233 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)+q(p7) [t-channel]'
c--  234 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)+q(p7) [t-channel, rad. in decay]'
c--  238 '  f(p1)+f(p2) --> t~(-->e-(p3)+nu~(p4)+b~(p5))+b(p6)+q(p7) [t-channel]'
c--  239 '  f(p1)+f(p2) --> t~(-->e-(p3)+nu~(p4)+b~(p5))+b(p6)+q(p7) [t-channel, rad. in decay]'
c--      '  f(p1)+f(p2) --> t(no BR) + b~(p6) + q(p7)' (removebr=.true.)
        if ((nproc .eq. 233) .or. (nproc .eq. 238)) then
	  case='4ftwdk'
	else
	  case='dk_4ft'
	endif
        nqcdjets=3
c--- the default is now that the parameter notag=0, so that
c--- the calculation requires the presence of the light jet that
c--- is present at LO; to compute an inclusive cross section, one
c--- can set notag=1, or use jet cuts that have no effect
	notag=0
c--- NB: for the studies in arXiv:1204.1513, FERMILAB-PUB-12-078-T we have set notag=1
c---     notag=1	
	if ((nproc .eq. 233) .or. (nproc .eq. 234)) then
	  nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
	else
	  nwz=-1
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='ba'
          plabel(6)='bq'
	endif
        plabel(7)='pp'
        plabel(8)='pp'
		
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=12
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=2   
        endif

c---  in the SM, the logical fourthgen should be false
c---  for BSM calculations, it should be true and it indicates that
c---   5 flavours should be used in the PDF and alpha-s
        fourthgen=.false.

	if (fourthgen) then
c--- BSM: full 5 light flavours
          nflav=5
	  bmass=4.7d0  !  set b-mass to its usual value
	else
c--- SM: only 4 light flavours
          nflav=4
          bmass=1001d0 !  enforce 4-flavour running in alfamz.f
        endif
	
c--- set up correct scales and as on heavy and light quark lines
        facscale_H=initfacscale_H
        facscale_L=initfacscale_L
        renscale_H=initrenscale_H
        renscale_L=initrenscale_L
c--- make sure it works even if not specifying separate scales
        if (initrenscale_L .eq. 0d0) then 
	  facscale_H=facscale
	  facscale_L=facscale
	  renscale_H=scale
	  renscale_L=scale
	endif
	
        b0=(xn*11d0-2d0*nflav)/6d0
	as_H=alphas(abs(renscale_H),amz,nlooprun)
	as_L=alphas(abs(renscale_L),amz,nlooprun)
	
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 235) .or. (nproc .eq. 240)) then
c--  235 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)+q(p7)+f(p8) [t-channel]'
c--  240 '  f(p1)+f(p2) --> t~(-->e-(p3)+nu~(p4)+b~(p5))+b(p6)+q(p7)+f(p8) [t-channel]'
c--      '  f(p1)+f(p2) --> t(no BR) + b~(p6) + q(p7) + f(p8)' (removebr=.true.)
        case='4ftjet'
        nqcdjets=4
c--- the default is now that the parameter notag=0, so that
c--- the calculation requires the presence of four light jets that
c--- are present at LO; to compute an inclusive cross section, one
c--- can set notag=1, or use jet cuts that have no effect
	notag=0
	if (nproc .eq. 235) then
	  nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
	else
	  nwz=-1
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='ba'
          plabel(6)='bq'
	endif
        plabel(7)='pp'
        plabel(8)='pp'
		
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=15
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=3   
        endif

c---  in the SM, the logical fourthgen should be false
c---  for BSM calculations, it should be true and it indicates that
c---   5 flavours should be used in the PDF and alpha-s
        fourthgen=.false.

	if (fourthgen) then
c--- BSM: full 5 light flavours
          nflav=5
	  bmass=4.7d0  !  set b-mass to its usual value
	else
c--- SM: only 4 light flavours
          nflav=4
          bmass=1001d0 !  enforce 4-flavour running in alfamz.f
        endif
	
c--- set up correct scales and as on heavy and light quark lines
        facscale_H=initfacscale_H
        facscale_L=initfacscale_L
        renscale_H=initrenscale_H
        renscale_L=initrenscale_L
c--- make sure it works even if not specifying separate scales
        if (initrenscale_L .eq. 0d0) then 
	  facscale_H=facscale
	  facscale_L=facscale
	  renscale_H=scale
	  renscale_L=scale
	endif
	
        b0=(xn*11d0-2d0*nflav)/6d0
	as_H=alphas(abs(renscale_H),amz,nlooprun)
	as_L=alphas(abs(renscale_L),amz,nlooprun)
	
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 241) .or. (nproc .eq. 246)
     .   .or. (nproc .eq. 242) .or. (nproc .eq. 247)) then
c--  241 '  f(p1)+f(p2) --> t(p3)+b~(p4)+f(p5) [s-channel]'
c--  246 '  f(p1)+f(p2) --> t~(p3)+b(p4)+f(p5) [s-channel]'

cc--- for comparison with C. Oleari's e+e- --> QQbg calculation
c	if (runstring(1:5) .eq. 'carlo') then
cc---     heavy quark mass passed via chars 6 and 7 of runstring
c	  read(runstring(6:7),67) imhq
c	  mt=dfloat(imhq)
c	  mb=mt
c          wmass=0d0
c	  write(6,*)
c	  write(6,*) ' >>> HEAVY QUARK MASS = ',mt,' GeV <<<'
c	endif
c   67   format(i2)
	
	if     ((nproc .eq. 241) .or. (nproc .eq. 246)) then
          case='qq_tbg'
          nqcdjets=1
          ndim=7
	elseif ((nproc .eq. 242) .or. (nproc .eq. 247)) then
          case='qqtbgg'
          nqcdjets=2
          ndim=10
	else
	  write(6,*) 'Unexpected value of nproc in chooser.f!'
c	  stop
	endif
        mass2=mt
        mass3=mb
	nflav=4
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
        plabel(6)='pp'
	if ((nproc .eq. 241) .or. (nproc .eq. 242)) then
	  nwz=+1
	else
	  nwz=-1
	endif

c--- set up correct scales and as on heavy and light quark lines
        facscale_H=initfacscale_H
        facscale_L=initfacscale_L
        renscale_H=initrenscale_H
        renscale_L=initrenscale_L
c--- make sure it works even if not specifying separate scales
        if (initrenscale_L .eq. 0d0) then 
	  facscale_H=facscale
	  facscale_L=facscale
	  renscale_H=scale
	  renscale_L=scale
	endif
	
	bmass=1001d0 ! since nflav=4
        b0=(xn*11d0-2d0*nflav)/6d0
	as_H=alphas(abs(renscale_H),amz,nlooprun)
	as_L=alphas(abs(renscale_L),amz,nlooprun)

	
c-----------------------------------------------------------------------

      elseif (nproc .eq. 249) then
c--  e+ e^- -> 3 jets as check of massless limit of process 241

cc--- for comparison with C. Oleari's e+e- --> QQbg calculation
c	if (runstring(1:5) .eq. 'carlo') then
cc---     heavy quark mass passed via chars 6 and 7 of runstring
c	  read(runstring(6:7),67) imhq
c	  mt=0d0
c	  mb=0d0
c          wmass=0d0
c	  write(6,*)
c	  write(6,*) ' >>> HEAVY QUARK MASS = ',mt,' GeV <<<'
c	endif
	
        case='epem3j'
        nqcdjets=1
        ndim=7

        mass2=0d0
        mass3=0d0
	nflav=4
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
        plabel(6)='pp'
	nwz=+1
	
	bmass=1001d0 ! since nflav=4
        b0=(xn*11d0-2d0*nflav)/6d0
	
c-----------------------------------------------------------------------

      elseif (nproc .eq. 251) then
c-- 251 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + W^+(-->nu(p5)+e^+(p6))+f(p7)+f(p8)'

        case='WpWp2j'
	nqcdjets=2
	ndim=16
	mb=0d0
	plabel(3)='nl'
	plabel(4)='ea'
	plabel(5)='nl'
	plabel(6)='ea'
	plabel(7)='pp'
	plabel(8)='pp'
	plabel(9)='pp'
	l1=1d0

	n2=1
	n3=1
	mass2=wmass
	width2=wwidth
	mass3=wmass
	width3=wwidth

	write(*,*)'Setting zerowidth to true for process 131'
	zerowidth = .true.
	write(*,*)'Setting removebr to false for process 131'
	removebr = .false.


c-----------------------------------------------------------------------

      elseif (nproc .eq. 252) then
c-- 252 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + W^+(-->nu(p5)+e^+(p6))+f(p7)+f(p8)+f(p9)'

        case='WpWp3j'
	nqcdjets=3
	ndim=19
	mb=0d0
	plabel(3)='nl'
	plabel(4)='ea'
	plabel(5)='nl'
	plabel(6)='ea'
	plabel(7)='pp'
	plabel(8)='pp'
	plabel(9)='pp'
	plabel(10)='pp'
	l1=1d0

	n2=1
	n3=1
	mass2=wmass
	width2=wwidth
	mass3=wmass
	width3=wwidth
        
	write(*,*)'Setting zerowidth to true for process 132'
	zerowidth = .true.
	write(*,*)'Setting removebr to false for process 132'
	removebr = .false.


c-----------------------------------------------------------------------

      elseif ((nproc .eq. 261) .or. (nproc .eq. 266)) then
c--  261 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)'
c--  266 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)[+b~(p6)]'
        case='gQ__ZQ'
        nqcdjets=1
        flav=5
        nwz=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        isub=1+(nproc-261)/5
        if (nproc .eq. 261) then
          plabel(6)='pp'
        else
          plabel(6)='ba'
        endif
        q1=-1d0
        l1=le
        r1=re

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif ((nproc .eq. 262) .or. (nproc .eq. 267)) then
c--  262 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)'
c--  267 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)[+c~(p6)]'
        case='gQ__ZQ'
        nqcdjets=1
        flav=4
        nwz=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        isub=1+(nproc-262)/5
        if (nproc .eq. 262) then
          plabel(6)='pp'
        else
          plabel(6)='ba'
        endif
        q1=-1d0
        l1=le
        r1=re
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif (nproc .eq. 263) then
c--  263 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b~(p5)+b(p6) (1 b-tag)'
        case='Zbbmas'
        nqcdjets=2
        notag=1

        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        write(6,*) 'mb=',mb
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        q1=-1d0
        l1=le
        r1=re
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif (nproc .eq. 264) then
c--  264 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c~(p5)+c(p6) (1 c-tag)'
        case='Zccmas'
        nqcdjets=2
        notag=1
        
        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        
        mb=mc
        write(6,*) 'mc=',mb
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        q1=-1d0
        l1=le
        r1=re

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif      
           
c-----------------------------------------------------------------------
          
      elseif ((nproc .eq. 270) .or. (nproc .eq. 271) 
     &   .or. (nproc .eq. 272)) then
      
c--- turn off Higgs decay, for speed
c        nodecay=.true.      
c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one
      
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
      
        n2=0
        n3=1

        mass3=hmass
        width3=hwidth
        
        if     (nproc .eq. 270) then
c-- 270 '  f(p1)+f(p2) --> H(gamma(p3)+gamma(p4))+f(p5)+f(p6)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)[in heavy top limit]' (removebr=.true.)
          hdecaymode='gaga'
          plabel(3)='ga'
          plabel(4)='ga'
          case='gagajj'
          nqcdjets=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            BrnRat=gamgambr
          endif
          
        elseif     (nproc .eq. 271) then
c-- 271 '  f(p1)+f(p2) --> H(b(p3)+b~(p4))+f(p5)+f(p6)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)[in heavy top limit]' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          case='ggfus2'
          nqcdjets=4
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=2
            BrnRat=br
          endif
          
        elseif (nproc .eq. 272) then
c-- 272 '  f(p1)+f(p2) --> H(tau-(p3)+tau+(p4))+f(p5)+f(p6)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)[in heavy top limit]' (removebr=.true.)
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          case='ggfus2'
          nqcdjets=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            Brnrat=tautaubr
          endif
        endif
                
c-----------------------------------------------------------------------

      elseif     (nproc .eq. 273) then
c-- 273 '  f(p1)+f(p2) -->` H(-->W^+(p3,p4)W^-(p5,p6)) + f(p7) + f(p8)'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one
      
        case='HWW2jt'
        ndim=16
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        nqcdjets=2
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth
c--- print warning if we're below threshold
        if (hmass .lt. 2d0*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=wwbr*brwen**2
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif

c-----------------------------------------------------------------------

      elseif     (nproc .eq. 274) then
c-- 274 f(p1)+f(p2)->H(Z^+(e^-(p3),e^+(p4))Z(mu^-(p5),mu^+(p6)))+f(p7)+f(p8)
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one
      
        case='HZZ2jt'
        l1=le
        r1=re
        l2=le
        r2=re
        ndim=16
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ml'
        plabel(6)='ma'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        nqcdjets=2
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth
c--- print warning if we're below threshold
        if (hmass .lt. 2d0*zmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=2d0*brzee**2*zzbr  ! factor of 2 for identical particles
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 275) .or. (nproc .eq. 276)) then

c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one

        case='ggfus3'
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=13
      
        n2=0
        n3=1

        mass3=hmass
        width3=hwidth
        
        if     (nproc .eq. 275) then
c-- 275 '  f(p1)+f(p2) --> H(b(p3)+b~(p4))+f(p5)+f(p6)+f(p7)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)+f(p7)[in heavy top limit]' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=5
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=3
            BrnRat=br
          endif
          
        elseif (nproc .eq. 276) then
c-- 276 '  f(p1)+f(p2) --> H(tau-(p3)+tau+(p4))+f(p5)+f(p6)+f(p7)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)+f(p7)[in heavy top limit]' (removebr=.true.)
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=3
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            BrnRat=tautaubr
          endif
        endif

c-----------------------------------------------------------------------

      elseif     (nproc .eq. 278) then
c-- 278 '  f(p1)+f(p2) --> H(-->W^+(p3,p4)W^-(p5,p6)) + f(p7) + f(p8) + f(p9)'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one
      
        case='HWW3jt'
        ndim=19
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        nqcdjets=3
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth
c--- print warning if we're below threshold
        if (hmass .lt. 2d0*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=wwbr*brwen**2
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif

c-----------------------------------------------------------------------

      elseif     (nproc .eq. 279) then
c-- 279 f(p1)+f(p2)->H(Z^+(e^-(p3),e^+(p4))Z(mu^-(p5),mu^+(p6)))+f(p7)+f(p8)+f(p9)
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one
      
        case='HZZ3jt'
        l1=le
        r1=re
        l2=le
        r2=re
        ndim=19
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ml'
        plabel(6)='ma'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        nqcdjets=3
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth
c--- print warning if we're below threshold
        if (hmass .lt. 2d0*zmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
	if (removebr) then
	write(6,*)
	write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
c        stop
	endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
c        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=2d0*brzee**2*zzbr  ! factor of 2 for identical particles
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif

c-----------------------------------------------------------------------

      elseif (nproc .eq. 280) then
             ndim=4 
             case='dirgam'
             plabel(3)='ga'
             plabel(4)='pp'
             plabel(5)='pp'
             nqcdjets=1
             n3=0
c--- Jet not required, just inclusive photon production
             notag=1
             inclusive=.true.
	     write(6,*)
	     write(6,*) 'Setting inclusive = .true. '//
     &                  'for direct photon production.'
	     
      elseif (nproc .eq. 282) then
             ndim=7 
             case='gamjet'
             plabel(3)='ga'
             plabel(4)='pp'
             plabel(5)='pp'
             plabel(6)='pp'
             nqcdjets=2
             n3=0
                
                
      elseif (nproc .eq. 285) then
             ndim=4 
             case='gamgam'
             plabel(3)='ga'
             plabel(4)='ga'
             plabel(5)='pp'
             nqcdjets=0
             n3=0

      elseif (nproc .eq. 286) then
             ndim=7 
             case='gmgmjt'
             plabel(3)='ga'
             plabel(4)='ga'
             plabel(5)='pp'
             plabel(6)='pp'
             nqcdjets=1
             n3=0
                
c-----------------------------------------------------------------------

c--- These two processes need to be moved to other numbers
      elseif (nproc .eq. 9280) then
c--  280      '  f(p1)+f(p2)--> f(p3)+f(p4)'
             ndim=4 
             case='twojet'
             plabel(3)='pp'
             plabel(4)='pp'
             plabel(5)='pp'
             nqcdjets=2
             n3=0
                
      elseif (nproc .eq. 9281) then
c--  281      '  f(p1)+f(p2)--> f(p3)+f(p4)+f(p5)'
             ndim=7 
             case='thrjet'
             plabel(3)='pp'
             plabel(4)='pp'
             plabel(5)='pp'
             plabel(6)='pp'
             nqcdjets=3
             n3=0
                
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 290) .or. (nproc .eq. 295)) then
        case='Wgamma'
        nqcdjets=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        plabel(5)='ga'
        plabel(6)='pp'

        if     (nproc .eq. 290) then
c-- 290 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+gamma(p5)'
c--     '  f(p1)+f(p2) --> W^+ (no BR) + gamma(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 295) then
c-- 295 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~+(p4))+gamma(p5)'
c--     '  f(p1)+f(p2) --> W^- (no BR) + gamma(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
      
        if (zerowidth .eqv. .false.) then
	  write(6,*)
          write(6,*) 'Setting removebr to .false. in order to ensure'
	  write(6,*) 'lepton-photon singularity can be removed'
	  removebr=.false.
	endif
      
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 292) .or. (nproc .eq. 297)) then 
       
         case='Wgajet'
         nqcdjets=1
         ndim=10
         mb=0
         rescale=.false.
         n2=0
         n3=1
         mass3=wmass
         width3=wwidth
         plabel(5)='ga'
         plabel(6)='pp'
        
         
        if     (nproc .eq. 292) then
c 292 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+(f(p5) --> gamma(p5))'
c     '  f(p1)+f(p2) --> W^+(No BR)+(f(p5) --> gamma(p5)) (removebr =.true.)'
           nwz=1
           plabel(3)='nl'
           plabel(4)='ea'
	   
        elseif (nproc .eq. 297) then
c 297 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+(f(p5) -->gamma(p5))'
c     '  f(p1)+f(p2) --> W^-(no BR)+(f(p5) -->gamma(p5)) (removebr=.true.)'
           nwz=-1
           plabel(3)='el'
           plabel(4)='na'
       endif
        
       if (zerowidth .eqv. .false.) then
	 write(6,*)
         write(6,*) 'Setting removebr to .false. in order to ensure'
	 write(6,*) 'lepton-photon singularity can be removed'
	 removebr=.false.
       endif
      
c---  total cross-section             
       if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
       endif

      
c--------------------------------------------------------------------------------------------------

        elseif ((nproc .eq. 300) .or. (nproc .eq. 305)) then
          case='Zgamma'
          nqcdjets=0
          ndim=7
          n2=0
          n3=1
          mass3=zmass
          width3=zwidth
          nwz=0
          plabel(5)='ga'
          plabel(6)='pp'
          
          if     (nproc .eq. 300) then
c-- 300 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+gamma(p5)'
c--     '  f(p1)+f(p2) --> Z^0 (no BR) +gamma(p5)' (removebr=.true.)
            call checkminzmass(1)
            plabel(3)='el'
            plabel(4)='ea'
            q1=-1d0
            l1=le
            r1=re
            if (zerowidth .eqv. .false.) then
	      write(6,*)
              write(6,*)'Setting removebr to .false. in order to ensure'
	      write(6,*)'lepton-photon singularity can be removed'
	      removebr=.false.
	    endif
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brzee
            endif
          elseif (nproc .eq. 305) then
c-- 305 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4)))-(sum over 3 nu)+gamma(p5)'
            plabel(3)='nl'
            plabel(4)='na'
            q1=0d0
            l1=ln*dsqrt(3d0)
            r1=rn*dsqrt(3d0)
          endif

c-----------------------------------------------------------------------

        elseif ((nproc .eq. 302) .or. (nproc .eq. 307)) then
          case='Zgajet'
          nqcdjets=1
          ndim=10
          n2=0
          n3=0
          mass3=zmass
          width3=zwidth
          nwz=0
          plabel(5)='ga'
          plabel(6)='pp'
          plabel(7)='pp'
          
          if     (nproc .eq. 302) then
c-- 302 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+gamma(p5)+f(p6)'
c--     '  f(p1)+f(p2) --> Z^0 (no BR) +gamma(p5)+jet(p6)' (removebr=.true.)
            call checkminzmass(1)
            plabel(3)='el'
            plabel(4)='ea'
            q1=-1d0
            l1=le
            r1=re
            if (zerowidth .eqv. .false.) then
	      write(6,*)
              write(6,*)'Setting removebr to .false. in order to ensure'
	      write(6,*)'lepton-photon singularity can be removed'
	      removebr=.false.
	    endif
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brzee
            endif
          elseif  (nproc .eq. 307) then
c-- 307 '  f(p1)+f(p2) --> Z^0(-->nu(p3)+nu~(p4))+gamma(p5)+f(p6)'
            plabel(3)='nl'
            plabel(4)='na'
            q1=0d0
            l1=ln*dsqrt(3d0)
            r1=rn*dsqrt(3d0)
          endif
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 301) .or. (nproc .eq. 306)) then
c-- 301 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4))+gamma(p5)+gamma(p6)'
c-- 306 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))+gamma(p5)+gamma(p6)'
        case='Z_2gam'
        ndim=10
        n2=0
        n3=0
        if (nproc .eq. 301) then
           call checkminzmass(1)
           plabel(3)='el'
           plabel(4)='ea'
           q1=-1d0
           l1=le
           r1=re
           if (zerowidth .eqv. .false.) then
	     write(6,*)
             write(6,*)'Setting removebr to .false. in order to ensure'
	     write(6,*)'lepton-photon singularity can be removed'
	     removebr=.false.
	   endif
c--- total cross-section             
           if (removebr) then
             plabel(3)='ig'
             plabel(4)='ig'
             call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
             BrnRat=brzee
           endif
        elseif (nproc .eq. 306) then
           plabel(3)='nl'
           plabel(4)='na'
           q1=0d0
           l1=ln*dsqrt(3d0)
           r1=rn*dsqrt(3d0)
        endif
        plabel(5)='ga'
        plabel(6)='ga'
        plabel(7)='pp'
        nwz=0   
        mass3=zmass
        width3=zwidth
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 303) .or. (nproc .eq. 308)) then
c-- 303 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4))+gamma(p5)+gamma(p6)+f(p7)'
c-- 308 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))+gamma(p5)+gamma(p6)+f(p7)'
        case='Z2gajt'
        nqcdjets=1
        ndim=13
        n2=0
        n3=0
        if (nproc .eq. 303) then
           call checkminzmass(1)
           plabel(3)='el'
           plabel(4)='ea'
           q1=-1d0
           l1=le
           r1=re
           if (zerowidth .eqv. .false.) then
	     write(6,*)
             write(6,*)'Setting removebr to .false. in order to ensure'
	     write(6,*)'lepton-photon singularity can be removed'
	     removebr=.false.
	   endif
c--- total cross-section             
           if (removebr) then
             plabel(3)='ig'
             plabel(4)='ig'
             call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
             BrnRat=brzee
           endif
        elseif (nproc .eq. 308) then
           plabel(3)='nl'
           plabel(4)='na'
           q1=0d0
           l1=ln*dsqrt(3d0)
           r1=rn*dsqrt(3d0)
        endif
        plabel(5)='ga'
        plabel(6)='ga'
        plabel(7)='pp'
        nwz=0   
        mass3=zmass
        width3=zwidth
                 
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 304) .or. (nproc .eq. 309)) then
c-- 304 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4))+gamma(p5)+f(p6)+f(p7)'
c-- 309 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))+gamma(p5)+f(p6)+f(p7)'
        case='Zga2jt'
        nqcdjets=2
        ndim=13
        n2=0
        n3=0
        if (nproc .eq. 304) then
           call checkminzmass(1)
           plabel(3)='el'
           plabel(4)='ea'
           q1=-1d0
           l1=le
           r1=re
           if (zerowidth .eqv. .false.) then
	     write(6,*)
             write(6,*)'Setting removebr to .false. in order to ensure'
	     write(6,*)'lepton-photon singularity can be removed'
	     removebr=.false.
	   endif
c--- total cross-section             
           if (removebr) then
             plabel(3)='ig'
             plabel(4)='ig'
             call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
             BrnRat=brzee
           endif
        elseif (nproc .eq. 309) then
           plabel(3)='nl'
           plabel(4)='na'
           q1=0d0
           l1=ln*dsqrt(3d0)
           r1=rn*dsqrt(3d0)
        endif
        plabel(5)='ga'
        plabel(6)='pp'
        plabel(7)='pp'
        nwz=0   
        mass3=zmass
        width3=zwidth
                         
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 311) .or. (nproc .eq. 316)) then
        case='W_bjet'
        nqcdjets=2
        flav=5
        isub=1
        
        nflav=5
        mb=0d0
        plabel(5)='bq'
        plabel(6)='pp'
        plabel(7)='pp'
        
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 311) then
c--  311 '  f(p1)+b(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+f(p6)'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 316) then
c--  316 '  f(p1)+b(p2) --> W^-(-->e^-(p3)+nu~(p4))+b(p5)+f(p6)'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+b(p2) --> W(no BR)+b(p5)+f(p6)' (removebr=.true.)
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 321) .or. (nproc .eq. 326)) then
        case='W_bjet'
        nqcdjets=2
        flav=4
        isub=1
        
        nflav=4
        mb=0d0
        plabel(5)='bq'
        plabel(6)='pp'
        plabel(7)='pp'
        
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 321) then
c--  321 '  f(p1)+b(p2) --> W^+(-->nu(p3)+e^+(p4))+c(p5)+f(p6)'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 326) then
c--  326 '  f(p1)+b(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5)+f(p6)'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+b(p2) --> W(no BR)+c(p5)+f(p6)' (removebr=.true.)
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 331) .or. (nproc .eq. 336)) then
        case='Wcjetg'
        nqcdjets=2
        nflav=3
        
        plabel(5)='bq'
        plabel(6)='pp'

        ndim=10
        mb=0
        n2=0
        n3=1
        mass2=0d0
        mass3=wmass
        width3=wwidth
        
        if     (nproc .eq. 331) then
c--  331 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+c(p5)+f(p6) [c-s interaction]'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 336) then
c--  336 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+c(p5)+f(p6) [c-s interaction]'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+f(p2) --> W(no BR)+c(p5)+f(p6) [c-s interaction]' (removebr=.true.)
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 341) .or. (nproc .eq. 351) 
     .   .or. (nproc .eq. 342) .or. (nproc .eq. 352)) then
        case='Z_bjet'
        call checkminzmass(1)
        ndim=10
        n2=0
        n3=1
        nqcdjets=2
        
        mb=0d0
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        
        if     ((nproc .eq. 341) .or. (nproc .eq. 351)) then
          isub=1        
          plabel(6)='pp'
          plabel(7)='pp'
        elseif ((nproc .eq. 342) .or. (nproc .eq. 352)) then
          isub=2        
          plabel(6)='ba'
          plabel(7)='pp'
        endif
        
        q1=-1d0
        l1=le
        r1=re
        nwz=0   
        mass3=zmass
        width3=zwidth

        if     ((nproc .eq. 341) .or. (nproc .eq. 342)) then
c--  341 '  f(p1)+b(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+f(p6)'
          flav=5
          nflav=5
        elseif ((nproc .eq. 351) .or. (nproc .eq. 352)) then
c--  351 '  f(p1)+c(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)+f(p6)'
          flav=4
          nflav=4
        endif
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 346) .or. (nproc .eq. 356) 
     .   .or. (nproc .eq. 347) .or. (nproc .eq. 357)) then
        case='Zbjetg'
        call checkminzmass(1)
        ndim=13
        n2=0
        n3=1
        nqcdjets=3
        
        mb=0d0
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        
        if     ((nproc .eq. 346) .or. (nproc .eq. 356)) then
c--  346 '  f(p1)+b(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+f(p6)+f(p7)'
          isub=1        
          plabel(6)='pp'
          plabel(7)='pp'
        elseif ((nproc .eq. 347) .or. (nproc .eq. 357)) then
c--  347 '  f(p1)+b(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+f(p6)+b~(p7)'
          isub=2        
          plabel(6)='ba'
          plabel(7)='pp'
        endif
        
        q1=-1d0
        l1=le
        r1=re
        nwz=0   
        mass3=zmass
        width3=zwidth

        if     ((nproc .eq. 346) .or. (nproc .eq. 347)) then
c--  346 '  f(p1)+b(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+f(p6)+f(p7)'
          flav=5
          nflav=5
        elseif ((nproc .eq. 356) .or. (nproc .eq. 357)) then
c--  356 '  f(p1)+c(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)+f(p6)+f(p7)'
          flav=4
          nflav=4
        endif
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif (nproc/10 .eq. 36) then
        case='W_only'
        n3=1
        ndim=4
        nqcdjets=0
C-- 361 '  c(p1)+sbar(p2) --> W^+(-->nu(p3)+e^+(p4))'
c--- This can be used to calculate cs->W at NLO, in the ACOT-like
c--- scheme where the charm quark appears in the initial state but
c--- the real corrections use a massless charm quark in the final state
c--- (c.f. processes 362 and 363 below)
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='pp'
        nwz=1

c--- total cross-section
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
	
c--- change W mass (after couplings and BRs already calculated)
c	if     (runstring(4:8) .eq. 'mw_80') then
c	  wmass=80.4d0
c	elseif (runstring(4:8) .eq. 'mw200') then
c	  wmass=200d0
c	elseif (runstring(4:8) .eq. 'mw400') then
c	  wmass=400d0
c	endif
	
c--- change charm mass
c	if     (runstring(9:13) .eq. 'mc1.3') then
c	  mc=1.3d0
c	  mcsq=mc**2
c	elseif (runstring(9:13) .eq. 'mc4.5') then
c	  mc=4.5d0
c	  mcsq=mc**2
c	elseif (runstring(9:13) .eq. 'mc20.') then
c	  mc=20d0
c	  mcsq=mc**2
c	endif
	
        mass3=wmass
        width3=wwidth
c--- set CKM matrix to remove all elements except for Vcs
        Vud=0d0
        Vus=0d0
        Vub=0d0
        Vcd=0d0
        Vcs=1d0
        Vcb=0d0

c--- To obtain a complete prediction for cs->W at NLO, including the
c---  the effect of the charm quark mass (as one should, according to
c---  ACOT) processes 362 and 363 should be summed (using 'tota')

c--- real matrix elements W+s and W+g, with corresponding (massless)
c---  integrated counterterms and virtual matrix elements
	if (nproc .eq. 362) case='Wcsbar'

c--- real matrix elements W+c including the charm quark mass,
c---  with the virtual contribution representing the counterterm
c---  consisting of the logarithm convolution	
	if (nproc .eq. 363) case='Wcs_ms'

        if ( (part .eq. 'lord') .and.  
     .       ((nproc .eq. 362) .or. (nproc .eq. 363)) ) then
          write(6,*) 'This process number is not suitable for the'
          write(6,*) 'LO calculation. Please run process 361'
          write(6,*) 'for the corresponding Born contribution.'
c          stop
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .ge. 401) .and. (nproc .le. 408)) then
        case='Wbbmas'
        write(6,*) 'mb=',mb
        flav=5

        bbproc=.false.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nqcdjets=2
	notag=1
	
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth           

        if     (nproc .le. 403) then
c--- 401  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5) [massive]'
c---      '  f(p1)+f(p2) --> W^+ (no BR) +b(p5) [massive]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .ge. 406) then
c--- 406  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5) [massive]'
c---      '  f(p1)+f(p2) --> W^- (no BR) +b(p5) [massive]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 411) .or. (nproc .eq. 416)) then
        case='W_bjet'
        nqcdjets=2
        flav=5
        isub=1
        
c--- check for Wb+X flag and allow one jet to be untagged in that case
	notag=1
	write(6,*)
        write(6,*)'****************************************************'
	write(6,*)'* WARNING: cuts allow final state of Wb+X          *'
        write(6,*)'****************************************************'
	
        nflav=5
        mb=0d0
        plabel(5)='bq'
        plabel(6)='pp'
        plabel(7)='pp'
        
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 411) then
c--  411 '  f(p1)+b(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+f(p6)'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 416) then
c--  416 '  f(p1)+b(p2) --> W^-(-->e^-(p3)+nu~(p4))+b(p5)+f(p6)'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+b(p2) --> W(no BR)+b(p5)+f(p6)' (removebr=.true.)
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 421) .or. (nproc .eq. 426)) then
        case='Wbbmas'
        write(6,*) 'mb=',mb
        flav=5

        bbproc=.false.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nqcdjets=2
	notag=1
	
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth           

        if     (nproc .eq. 421) then
c--- 421  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5) [massive]'
c---      '  f(p1)+f(p2) --> W^+ (no BR) +b(p5) [massive]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .ge. 426) then
c--- 426  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5) [massive]'
c---      '  f(p1)+f(p2) --> W^- (no BR) +b(p5) [massive]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 431) .or. (nproc .eq. 436)) then
        case='Wbbjem'
        write(6,*) 'mb=',mb
	nqcdjets=3
        flav=5
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        ndim=13
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
 
        if     (nproc .eq. 431) then
c-- 431 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +b(p5)+b~(p6)+f(p7) [massive]'
c--     '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6)+f(p7) [massive]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 436) then
c-- 436 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) +b(p5)+b~(p6)+f(p7) [massive]'
c--     '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6)+f(p7) [massive]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section	     
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
	     
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 500) .or. (nproc .eq. 510)) then
        case='Wttmas'
        write(6,*) 'mt=',mt
        flav=6
        plabel(5)='ig'
        plabel(6)='ig'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
 
        if     (nproc .eq. 500) then
C-- 500 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +t(p5)+t~(p6) [massive]' 'N'
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 510) then
C-- 510 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(p5)+t~(p6) [massive]' 'N'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 501) .or. (nproc .eq. 502)
     &   .or. (nproc .eq. 511) .or. (nproc .eq. 512)
     &   .or. (nproc .eq. 503) .or. (nproc .eq. 513)
     &   .or. (nproc .eq. 506) .or. (nproc .eq. 516)) then
        case='qq_ttw'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(11)='pp'
        ndim=20
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 501) then 
c-- 501 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^+(nu(p9),mu^+(p10))'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='nl'
          plabel(10)='ea'
	  nqcdjets=2
        elseif   (nproc .eq. 502) then 
c-- 502 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^+(nu(p9),mu^+(p10))[rid]'
          case='ttwldk'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='nl'
          plabel(10)='ea'
	  nqcdjets=2
        elseif (nproc .eq. 503) then 
c-- 503 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+q(p7)+q~(p8))+W^+(nu(p9),mu^+(p10))'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='nl'
          plabel(10)='ea'
	  nqcdjets=4
c--- for CMS study: we can allow 2 jets to be lost
          notag=2	  
        elseif (nproc .eq. 506) then 
c-- 506 '  f(p1)+f(p2) --> t(-->q(p3)+q~(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^+(nu(p9),mu^+(p10))'
          nwz=+1
          plabel(3)='pp'
          plabel(4)='pp'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='nl'
          plabel(10)='ea'
	  nqcdjets=4
c--- for CMS study: we can allow 2 jets to be lost
          notag=2	  
        elseif (nproc .eq. 511) then 
c-- 511 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^-(mu^-(p9),nu~(p10))'
          nwz=-1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='el'
          plabel(10)='na'
	  nqcdjets=2
        elseif (nproc .eq. 512) then 
c-- 512 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^-(mu^-(p9),nu~(p10))[rid]'
          case='ttwldk'
          nwz=-1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='el'
          plabel(10)='na'
	  nqcdjets=2
        elseif (nproc .eq. 513) then 
c-- 513 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+q(p7)+q~(p8))+W^-(mu^-(p9),nu~(p10))'
          nwz=-1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='el'
          plabel(10)='na'
	  nqcdjets=4
c--- for CMS study: we can allow 2 jets to be lost
          notag=2	  
        elseif (nproc .eq. 516) then 
c-- 516 '  f(p1)+f(p2) --> t(-->q(p3)+q~(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^-(mu^-(p9),nu~(p10))'
          nwz=-1
          plabel(3)='pp'
          plabel(4)='pp'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='el'
          plabel(10)='na'
	  nqcdjets=4
c--- for CMS study: we can allow 2 jets to be lost
          notag=2	  
        endif

        if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brtop*brwen)**2*brwen
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            plabel(9)='ig'
            plabel(10)='ig'
	    nqcdjets=0
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 530) .or. (nproc .eq. 531)) then
        case='qq_ttz'
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        nwz=1

        ndim=22
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth

        if     (nproc .eq. 530) then 
c--  530 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->nu~(p7)+e^-(p8)+b~(p6))+Z(e(p9),e~(p10))'
          plabel(9)='el'
          plabel(10)='ea'
          q1=-1d0
          l1=le
          r1=re
	  mb=0d0
          if (removebr) then
            q1=0d0
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brtop*brwen)**2*brzee
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            plabel(9)='ig'
            plabel(10)='ig'
          endif
        elseif (nproc .eq. 531) then 
c--  531 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->nu~(p7)+e^-(p8)+b~(p6))+Z(b(p9),b~(p10))'
          plabel(9)='bq'
          plabel(10)='ba'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
	  mb=0d0
        endif

c-----------------------------------------------------------------------

      elseif (nproc/10 .ge. 90) then
        write(6,*) 'Setting part to lord and zerowidth to false'
        zerowidth=.false.
        part='lord'
        if     (nproc .eq. 902) then
          case='vlchk2'
          nwz=1
          ndim=4
          n3=1
          mass3=wmass
          width3=wwidth
        elseif (nproc .eq. 903) then
          case='vlchk3'
          nwz=1
          ndim=7
          n3=1
          mass3=wmass
          width3=wwidth
        elseif (nproc .eq. 904) then
          case='vlchk4'
          nwz=1
          ndim=10
          n2=1
          n3=1
          mass2=zmass
          width2=zwidth
          mass3=wmass
          width3=wwidth
        elseif (nproc .eq. 905) then
          case='vlchk5'
          nwz=1
          ndim=13
          n2=1
          n3=1
          mass2=hmass
          width2=hwidth
          mass3=wmass
          width3=wwidth
        elseif (nproc .eq. 906) then
          case='vlchk6'
          nwz=1
          ndim=16
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=mt
          width3=twidth
        elseif (nproc .eq. 908) then
          case='vlchk8'
          nwz=1
          ndim=22
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=mt
          width3=twidth
        elseif (nproc .eq. 909) then
          case='vlchkm'
          write(6,*) 'mb=',mb
          nwz=1
          ndim=10
          n2=1
          n3=1
          mass2=hmass
          width2=hwidth
          mass3=wmass
          width3=wwidth
        elseif (nproc .eq. 910) then
          case='vlchm3'
          write(6,*) 'mt=',mt
          nwz=1
          ndim=7
          n2=0
          n3=0
          mass2=mt
          width2=twidth
          mass3=mt
          width3=twidth
        elseif (nproc .eq. 911) then
          case='vlchwt'
          write(6,*) 'mt=',mt

          write(6,*) 'Setting zerowidth = .true.'
          zerowidth=.true.
             
          ndim=13
          mb=0
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=wmass
          width3=wwidth             
        elseif (nproc .eq. 912) then
          case='vlchwn'
          write(6,*) 'mt=',mt

          write(6,*) 'Setting zerowidth = .true.'
          zerowidth=.true.
             
          ndim=7
          mb=0
          n2=0
          n3=1
          mass3=wmass
          width3=wwidth
        elseif (nproc .eq. 913) then
          case='vlchwg'
          write(6,*) 'mt=',mt

          write(6,*) 'Setting zerowidth = .true.'
          zerowidth=.true.
             
          ndim=16
          mb=0
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=wmass
          width3=wwidth             
          
        elseif (nproc .eq. 914) then
          case='vlchwh'
          write(6,*) 'mt=',mt

          write(6,*) 'Setting zerowidth = .true.'
          zerowidth=.true.
             
          ndim=16
          mb=0
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=wmass
          width3=wwidth             
          
        endif
      else 
        call nprocinvalid()
      endif

c--- set up alpha-s again (in case nflav was changed)
      call coupling2

c--- remove 2 dimensions from integration if decay is not included
      if (nodecay) ndim=ndim-2

c--- report on the removed BR, if necessary
      if (removebr) then
        write(6,*)'****************************************************'
        write(6,98) BrnRat
        write(6,*)'****************************************************'
      endif

c--- check that calculation can be performed
      call checkorder(order)

c--- fill up CKM matrix
      call ckmfill(nwz)

c--- set flags to true unless we're doing W+2 jet or Z+2 jet
      if ( ((case .ne. 'W_2jet') .and. (case .ne. 'Z_2jet'))
     . .or. (part .eq. 'lord') ) then
        Qflag=.true.
        Gflag=.true.
      endif

      return

 43   write(6,*) 'problems opening process.DAT'
c      stop

 44   write(6,*) 'Unimplemented process number, nproc = ',nproc, 
     . ' mcfm halted'
c      stop
 
 98   format(' *             Brn.Rat. removed = ',  f11.7, '       *')
     
      end

      subroutine nprocinvalid()
      implicit none
      integer nproc
      common/nproc/nproc

      write(6,*) 'chooser: Unimplemented case'
      write(6,*) 'nproc=',nproc      
c      stop
      
      return 
      end
      
      subroutine checkminzmass(i)
c--- Checks that the minimum invariant mass specified in the options
c--- file is not zero for boson 34 (i=1) or boson 56 (i=2)
      implicit none
      include 'limits.f'
      include 'zerowidth.f'
      integer i

c--- if generating exactly on-shell, there's nothing to worry about      
      if (zerowidth) return
      
      if ((i .eq. 1) .and. (wsqmin .eq. 0d0)) then
        write(6,*)
        write(6,*) 'Please set m34min not equal to zero to'
        write(6,*) 'prevent the virtual photon from becoming real.'
c        stop
      endif

      if ((i .eq. 2) .and. (bbsqmin .eq. 0d0)) then
        write(6,*)
        write(6,*) 'Please set m56min not equal to zero to'
        write(6,*) 'prevent the virtual photon from becoming real.'
c        stop
      endif
      
      return
      end
      
      
