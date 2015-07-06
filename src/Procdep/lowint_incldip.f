      double precision function lowint(r,wgt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'noglue.f'
      include 'process.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'phasemin.f'
      include 'PDFerrors.f'
      include 'wts_bypart.f'
      include 'stopscales.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'outputoptions.f'
c --- DSW. To store flavour information :
      include 'nflav.f'
c      include 'b0.f'
c --- DSW.
      integer pflav,pbarflav
c --- To use VEGAS random number sequence :
      double precision ran2
      integer ih1,ih2,j,k,nvec,sgnj,sgnk,ii
      double precision r(mxdim),W,sqrts,xmsq,val,val2,ptmp,
     . fx1(-nf:nf),fx2(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     . pswt,rscalestart,fscalestart,
     . fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf),
     . fxb1(-nf:nf),fxb2(-nf:nf)
      double precision wgt,msq(-nf:nf,-nf:nf),m3,m4,m5,xmsqjk
     & ,msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
      double precision xx(2),flux,vol,vol_mass,vol3_mass,vol_wt,BrnRat
      double precision xmsq_bypart(-1:1,-1:1),lord_bypart(-1:1,-1:1)
      logical bin,first,includedipole
      logical creatent,dswhisto
c      double precision gx1(-nf:nf),gx2(-nf:nf)
c      integer idum
c      COMMON/ranno/idum
      character*30 runstring
      double precision b1scale,q2scale,q1scale,b2scale
      external qg_tbq,BSYqqb_QQbdk_gvec,qqb_QQbdk,qg_tbqdk,qg_tbqdk_gvec
      common/runstring/runstring
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/x1x2/xx
      common/BrnRat/BrnRat
      common/bypart/lord_bypart
      common/outputflags/creatent,dswhisto
      common/bqscale/b1scale,q2scale,q1scale,b2scale
      data p/48*0d0/
      data first/.true./
      save first,rscalestart,fscalestart
      
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif
      
      ntotshot=ntotshot+1
      lowint=0d0
c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0d0
      
      W=sqrts**2
c--- processes that use "gen2"     
      if     ( (case .eq. 'W_only')
     .    .or. (case .eq. 'Z_only')
     .    .or. (case .eq. 'Higaga')
     .    .or. (case .eq. 'vlchk2') ) then
        if (case .eq. 'vlchk2') then
          wsqmin=0d0
          wsqmax=sqrts**2
        endif
        npart=2
        if (new_pspace) then
          call gen2a(r,p,pswt,*999)
        else
          call gen2(r,p,pswt,*999)
        endif

c--- processes that use "gen2jet"     
      elseif ((case .eq. 'twojet') 
     .   .or. (case .eq. 'dirgam')
     .   .or. (case .eq. 'gamgam')) then
        npart=2
        call gen2jet(r,p,pswt,*999)

c--- processes that use "gen2m"     
      elseif ( (case .eq. 'ggfus0')
     .    .or. (case .eq. 'tt_tot')
     .    .or. (case .eq. 'bb_tot')
     .    .or. (case .eq. 'cc_tot') ) then
        npart=2
        call gen2m(r,p,pswt,*999)
          
c--- processes that use "gen3"     
      elseif ( (case .eq. 'W_cjet') 
     .   .or.  (case .eq. 'Wbfrmc')
     .   .or.  (case .eq. 'W_tndk')
     .   .or.  (case .eq. 'vlchwn')
     .   .or.  (case .eq. 'epem3j') ) then
        npart=3
        call gen3(r,p,pswt,*999)

c--- processes that use "gen3h"     
      elseif (case .eq. 'Hi_Zga') then
        npart=3
        call gen3h(r,p,pswt,*999)

c---  processes that use "gen3jet"     
      elseif ( (case .eq. 'Wgamma') 
     .   .or.  (case .eq. 'Zgamma')
     .   .or.  (case .eq. 'W_frag') 
     .   .or.  (case .eq. 'Z_frag')
     .   .or.  (case .eq. 'vlchk3') ) then
        if (case .eq. 'vlchk3') then
          wsqmin=0d0
          wsqmax=sqrts**2
        endif
        npart=3
        call gen_phots_jets(r,1,0,p,pswt,*999)
        
c--- processes that use "gen3jetgaga"     
      elseif (case .eq. 'gmgmjt') then
        npart=3
        call gen3jetgaga(r,p,pswt,*999)

c--- processes that use "gen_phots_jets"     
c--- special treatment for W+gamma+jet and Z+gamma+jet 
      elseif ( (case .eq. 'Wgajet') .or. (case .eq. 'Zgajet') ) then
        npart=4
	if (new_pspace) then
          call gen_vgamj(r,p,pswt,*999) ! New PS routine
	else
        if     (ipsgen .eq. 1) then
          call gen_phots_jets(r,1,1,p,pswt,*999)
        elseif (ipsgen .eq. 2) then
          call gen_phots_jets_dkrad(r,1,1,p,pswt,*999)
	else
	  write(6,*) 'Parameter ipsgen should be 1 or 2'
	  write(6,*) 'ipsgen = ',ipsgen
	  stop
        endif
	endif
	 
c--- special treatment for Z+gamma+gamma 
      elseif ( (case .eq. 'Z_2gam') ) then
        npart=4
        if  (ipsgen .eq. 1) then
            call gen_phots_jets(r,2,0,p,pswt,*999) !AA+AB
        elseif  (ipsgen .eq. 2) then
            call gen_phots_jets_dkrad2(r,2,0,p,pswt,*999) !BB+BC
        elseif  (ipsgen .eq. 3) then
           call gen_phots_jets_dkrad(r,2,0,p,pswt,*999) !CC+AC+CD
        elseif  (ipsgen .eq. 4) then
           call gen_phots_jets_dkrad(r,2,0,p,pswt,*999) !DD+AD+BD
           do ii=1,4                
              ptmp=p(5,ii)
              p(5,ii)=p(6,ii)
              p(6,ii)=ptmp
           enddo
        else
           write(6,*) 'Parameter ipsgen should be 1 or 2 or 3 or 4'
           write(6,*) 'ipsgen = ',ipsgen
           stop
        endif

c--- special treatment for Z+gamma+gamma+jet 
      elseif ( (case .eq. 'Z2gajt') ) then
        npart=5
        if (ipsgen .eq. 1) then
           call gen_phots_jets(r,2,1,p,pswt,*999)  !AA+AB
        elseif (ipsgen .eq. 2) then
           call gen_phots_jets_dkrad2(r,2,1,p,pswt,*999)  !BB+BC
        elseif (ipsgen .eq. 3) then
           call gen_phots_jets_dkrad(r,2,1,p,pswt,*999)  !CC+AC+CD 
        elseif (ipsgen .eq. 4) then
           call gen_phots_jets_dkrad(r,2,1,p,pswt,*999)  !DD+AD+BD
           do ii=1,4               
              ptmp=p(5,ii)
              p(5,ii)=p(6,ii)
              p(6,ii)=ptmp
           enddo
        else
           write(6,*) 'Parameter ipsgen should be 1 or 2 or 3 or 4'
           write(6,*) 'ipsgen = ',ipsgen
           stop
        endif

c--- special treatment for Z+gamma+jet+jet 
      elseif ( (case .eq. 'Zga2jt') ) then
        npart=5
        if (ipsgen .eq. 1) then
          call gen_phots_jets(r,1,2,p,pswt,*999)
        elseif (ipsgen .eq. 2) then
          call gen_phots_jets_dkrad(r,1,2,p,pswt,*999)
        else
          write(6,*) 'Parameter ipsgen should be 1 or 2'
          write(6,*) 'ipsgen = ',ipsgen
          stop
        endif

      elseif ((case .eq. 'thrjet') .or. (case .eq. 'gamjet')) then
!        call gen3jet(r,p,pswt,*999)
        call gen_photons_jets(r,1,2,p,pswt,*999)
        npart=3
c--- processes that use "gen3m"
      elseif ( (case .eq. 'tottth') ) then
        m3=mt
        m4=mt
        m5=hmass
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)
          
      elseif ( (case .eq. 'totttz') ) then
        m3=mt
        m4=mt
        m5=zmass
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)
	  
c--- processes that use "gen3m"
      elseif (case .eq. 'tt_glu') then
        m3=mt
        m4=mt
        m5=0d0
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)
	
c--- processes that use "gen3m"
      elseif ((case .eq. 'qg_tbq') .or. (case .eq. 'qq_tbg')) then
        m3=mt
        m4=mb
        m5=0d0
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)
	
c--- processes that use gen3mdk
      elseif ( (case .eq. '4ftwdk') ) then
        m3=mt
        m4=mb
        m5=0d0
        npart=5
        call gen3mdk(r,p,m3,m4,m5,pswt,*999)
        
c--- processes that use "gen3m_rap"     
      elseif ( (case .eq. 'vlchm3') ) then
        taumin=(2d0*mt/sqrts)**2
        m3=mt
        m4=mt
        npart=3
        call gen3m_rap(r,p,m3,m4,pswt,*999)
          
c--- processes that use "gen4"     
      elseif (case .eq. 'qgtbqq') then
        npart=4
        call gen4(r,p,pswt,*999)
                  
c--- processes that use "gen4h"     
      elseif ( (case .eq. 'HWW_4l')
     .    .or. (case .eq. 'HWW2lq')
     .    .or. (case .eq. 'HWW_tb')
     .    .or. (case .eq. 'HWWint')
     .    .or. (case .eq. 'HZZ_4l') ) then
        npart=4
        call gen4h(r,p,pswt,*999)
         
c--- processes that use "gen4mdk"     
      elseif (case .eq. '4ftjet') then
        npart=6
        call gen4mdk(r,p,pswt,*999)
                  
c--- processes that use "gen5" 
      elseif ( (case .eq. 'W_twdk') 
     . .or.    (case .eq. 'HWWjet') 
     . .or.    (case .eq. 'HZZjet') 
     . .or.    (case .eq. 'WW_jet') 
     . .or.    (case .eq. 'ZZ_jet') 
     . .or.    (case .eq. 'vlchwt') 
     . ) then 
        npart=5 
        call gen5(r,p,pswt,*999)

c--- processes that use "gen6"     
      elseif ( (case .eq. 'tt_bbl')
     .    .or. (case .eq. 'ttZbbl')
     .    .or. (case .eq. 'tt_bbh')
     .    .or. (case .eq. 'tt_bbu')
     .    .or. (case .eq. 'tautau')
     .    .or. (case .eq. 'Wtbwdk')
     .    .or. (case .eq. 'vlchk6')
     .    .or. (case .eq. 'vlchwg')
     .    .or. (case .eq. 'vlchwh')
     .    .or. (case .eq. 'qq_HWW')
     .    .or. (case .eq. 'qq_HZZ')
     .    .or. (case .eq. 'HWW2jt')
     .    .or. (case .eq. 'HZZ2jt')
     .    .or. (case .eq. 'WH__ZZ')
     .    .or. (case .eq. 'WH__WW')  
     .    .or. (case .eq. 'ZH__WW')
     .    .or. (case .eq. 'ZH__ZZ')
     .    .or. (case .eq. 'WpWp2j')) then
        npart=6
        call gen6(r,p,pswt,*999)

c--- processes that use "gen7"     
      elseif ( (case .eq. 'qq_ttg') ) then
        m3=mt
        m4=mt
        m5=0d0
        npart=7
        call gen7m(r,p,m3,m4,m5,pswt,*999)
c        npart=7
c        call gen7_rap(r,p,pswt,*999)
c	call writeout(p)

      elseif ( (case.eq.'WpWp3j')
     &     .or.(case.eq.'HWW3jt')
     &     .or.(case.eq.'HZZ3jt') )  then
        npart=7
	call gen7(r,p,pswt,*999)

c--- processes that use "gen8"     
      elseif ( (case .eq. 'qq_tth') 
     .    .or. (case .eq. 'qq_ttz') 
     .    .or. (case .eq. 'vlchk8') ) then
        npart=8
        call gen8(r,p,pswt,*999)
          
      elseif ( (case .eq. 'qq_ttw')  ) then
        npart=8
        taumin=(2d0*mt/sqrts)**2
        call gen8_rap(r,p,pswt,*999)
          
c--- processes that use "gen_njets" with an argument of "1"     
      elseif ( (case .eq. 'W_1jet')
     .    .or. (case .eq. 'Wcjet0')
     .    .or. (case .eq. 'Z_1jet')
     .    .or. (case .eq. 'H_1jet')
     .    .or. (case .eq. 'httjet')
     .    .or. (case .eq. 'attjet')
     .    .or. (case .eq. 'ggfus1')
     .    .or. (case .eq. 'Hgagaj')
     .    .or. (case .eq. 'gQ__ZQ') ) then
        npart=3
        call gen_njets(r,1,p,pswt,*999)
         
c--- processes that use "gen_njets" with an argument of "2"
      elseif ( (case .eq. 'Wbbbar')
     .    .or. (case .eq. 'W_2jet')
     .    .or. (case .eq. 'Z_2jet')
     .    .or. (case .eq. 'Zbbbar')
     .    .or. (case .eq. 'qq_Hqq')
     .    .or. (case .eq. 'qq_Hgg')
     .    .or. (case .eq. 'ggfus2')
     .    .or. (case .eq. 'gagajj')
     .    .or. (case .eq. 'W_bjet')
     .    .or. (case .eq. 'Wcjetg')
     .    .or. (case .eq. 'Z_bjet') ) then
        npart=4
        call gen_njets(r,2,p,pswt,*999)         
        
c--- processes that use "gen_njets" with an argument of "3"
      elseif ( (case .eq. 'W_3jet') 
     .    .or. (case .eq. 'Wbbjet') 
     .    .or. (case .eq. 'Z_3jet') 
     .    .or. (case .eq. 'Zbbjet') 
     .    .or. (case .eq. 'Wb2jet') 
     .    .or. (case .eq. 'qqHqqg')
     .    .or. (case .eq. 'ggfus3')
     .    .or. (case .eq. 'Zbjetg')
     .    .or. (case .eq. 'vlchk5') ) then
        npart=5
        call gen_njets(r,3,p,pswt,*999)      
c--- processes that use "gen_stop" with an argument of "1" (number of extra jets)
      elseif ( (case .eq. 'bq_tpq')
     .    .or. (case .eq. 'ttdkay')
     .    .or. (case .eq. 't_bbar')
     .    .or. (case .eq. 'tdecay') ) then
        npart=4
        call gen_stop(r,1,p,pswt,*999)
	
c--- DEFAULT: processes that use "gen4"
      else
        if ((case .eq. 'vlchk4') .or. (case .eq. 'vlchkm')) then
          wsqmin=0d0
          wsqmax=sqrts**2
        endif
        npart=4
        call gen4(r,p,pswt,*999)      
      endif

      nvec=npart+2
      call dotem(nvec,p,s)
      
c--- (moved to includedipole) impose cuts on final state
c      if (case .ne. 'vlchk6' .and. case .ne. 'tautau') then
c        call masscuts(p,*999)
c      endif

c---- reject event if any s(i,j) is too small
      call smalls(s,npart,*999)                                             

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif
      
      if (dynamicscale) call scaleset(rscalestart,fscalestart,p)
      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)      

c--- Calculate the required matrix elements      
      if     (case .eq. 'W_only') then
        call qqb_w(p,msq)
      elseif (case .eq. 'W_1jet') then
        call qqb_w_g(p,msq)
      elseif (case .eq. 'Wgamma') then
c        call compare_madgraph(p,qqb_wgam,qqb_wgam_mad) ! Checked: JC, 12/29/10
        call qqb_wgam(p,msq)
      elseif (case .eq. 'Wgajet') then
         call qqb_wgam_g(p,msq)
      elseif (case .eq. 'Wbfrmc') then
        call qqb_wbfromc(p,msq)
      elseif (case .eq. 'W_cjet') then
        call qqb_w_cjet(p,msq)
      elseif (case .eq. 'Wcjet0') then
        call qqb_w_cjet_massless(p,msq)
      elseif (case .eq. 'Wbbmas') then
        call qqb_wbbm(p,msq)
      elseif (case .eq. 'Wttmas') then
        call qqb_wbbm(p,msq)
      elseif (case .eq. 'Wbbbar') then
        call qqb_wbb(p,msq)
      elseif (case .eq. 'W_2jet') then
        call qqb_w2jet(p,msq)
      elseif (case .eq. 'W_3jet') then
        call qqb_w2jet_g(p,msq)
      elseif (case .eq. 'Wbbjet') then
        call qqb_wbb_g(p,msq)
      elseif (case .eq. 'Z_only') then
        call qqb_z(p,msq)
      elseif (case .eq. 'Z_1jet') then
        call qqb_z1jet(p,msq)
      elseif (case .eq. 'Z_2jet') then
        call qqb_z2jet(p,msq)
      elseif (case .eq. 'Z_3jet') then
        call qqb_z2jet_g(p,msq)
      elseif (case .eq. 'Zgamma') then
c        call compare_madgraph(p,qqb_zgam,qqb_zgam_mad) ! Checked: JC, 01/03/11
        call qqb_zgam(p,msq)
      elseif (case .eq. 'Z_2gam') then
        call qqb_zaa(p,msq)
      elseif (case .eq. 'Zgajet') then
        call qqb_zaj(p,msq)
c        call qqb_zgam_g(p,msq)	! Old routine
      elseif (case .eq. 'Z2gajt') then
        call qqb_zaa_g(p,msq)
      elseif (case .eq. 'Zga2jt') then
        call qqb_zaj_g(p,msq)
      elseif (case .eq. 'Zbbmas') then
        call qqb_zbbm(p,msq)
      elseif (case .eq. 'Zbbbar') then
        call qqb_zbb(p,msq)
      elseif (case .eq. 'Zbbjet') then
        call qqb_zbb_g(p,msq)
      elseif (case .eq. 'WWqqbr') then
        call qqb_ww(p,msq)
      elseif (case .eq. 'WWnpol') then
        call qqb_ww_unpol(p,msq)
      elseif (case .eq. 'WW_jet') then
        call qqb_ww_g(p,msq)
      elseif (case .eq. 'WpWp2j') then
        call qqb_wpwp_qqb(p(1:12,:),msq(-5:5,-5:5))
      elseif (case .eq. 'WpWp3j') then
        call qqb_wpwp_qqb_g(p(1:12,:),msq(-5:5,-5:5))
      elseif (case .eq. 'WZbbar') then
        call qqb_wz(p,msq)
      elseif (case .eq. 'ZZlept') then
        call qqb_zz(p,msq)
      elseif (case .eq. 'ZZ_jet') then
        call qqb_zz_g(p,msq)
      elseif (case .eq. 'WHbbar') then
        call qqb_wh(p,msq)
c      elseif (case .eq. 'twojet') then
c        call qqb_2jet(p,msq)
c      elseif (case .eq. 'thrjet') then
c        call qqb_3jet(p,msq)
      elseif (case .eq. 'dirgam') then
        call qqb_dirgam(p,msq)
      elseif (case .eq. 'gamgam') then
        call qqb_gamgam(p,msq)
c        write(6,*) 'gamgam',msq(1,-1),msq(2,-2),msq(-1,1),msq(-2,2) 
c        call qqb_gamgam_mad(p,msq)
c        write(6,*) 'gamgam',msq(1,-1),msq(2,-2),msq(-1,1),msq(-2,2) 
c        pause
      elseif (case .eq. 'gmgmjt') then
        call qqb_gamgam_g(p,msq)
      elseif (case .eq. 'gamjet') then
        call qqb_dirgam_g(p,msq)
      elseif (case .eq. 'WH__WW') then
        call qqb_wh_ww(p,msq)
      elseif (case .eq. 'WH__ZZ') then
        call qqb_wh_zz(p,msq)
      elseif (case .eq. 'WHgaga') then
        call qqb_wh_gaga(p,msq)
      elseif (case .eq. 'ZHbbar') then
        call qqb_zh(p,msq)
      elseif (case .eq. 'ZH__WW') then
        call qqb_zh_ww(p,msq)
      elseif (case .eq. 'ZH__ZZ') then
        call qqb_zh_zz(p,msq)
      elseif (case .eq. 'ZHgaga') then
        call qqb_zh_gaga(p,msq)
      elseif (case .eq. 'ggfus0') then
        call gg_h(p,msq)
      elseif (case .eq. 'Higaga') then
        call gg_hgamgam(p,msq)
      elseif (case .eq. 'Hi_Zga') then
        call gg_hzgam(p,msq)
      elseif (case .eq. 'HWW_4l') then
        call qqb_hww(p,msq)
      elseif (case .eq. 'HWW2lq') then
        call qqb_hww(p,msq)
      elseif (case .eq. 'HWW_tb') then
        call qqb_hww_tb(p,msq)
      elseif (case .eq. 'HWWint') then
        call gg_ww_int(p,msq)
      elseif (case .eq. 'HZZ_4l') then
        call qqb_hzz(p,msq)
      elseif (case .eq. 'H_1jet') then
        call qqb_hg(p,msq)
      elseif (case .eq. 'ttZbbl') then
        call qqbZtt(p,msq)
      elseif ((case .eq. 'tt_bbl') .or. (case .eq. 'tt_bbh')) then
        call qqb_QQbdk(p,msq)
      elseif (case .eq. 'tt_bbu') then
        call qqb_QQbdku(p,msq)
      elseif (case .eq. 'qq_ttg') then
        call qqb_QQbdk_g(p,msq)
      elseif (case .eq. 'tt_tot') then
        call qqb_QQb(p,msq)
      elseif (case .eq. 'bb_tot') then
        call qqb_QQb(p,msq)
      elseif (case .eq. 'cc_tot') then
        call qqb_QQb(p,msq)
      elseif (case .eq. 'tt_glu') then
        call qqb_QQb_g(p,msq)
      elseif (case .eq. 'bq_tpq') then
        call bq_tpq(p,msq)
      elseif (case .eq. 'ttdkay') then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (case .eq. 't_bbar') then
	call qqb_tbbdk(p,msq)
      elseif (case .eq. 'tdecay') then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (case .eq. 'W_tndk') then
        call qqb_w_tndk(p,msq)
      elseif (case .eq. 'W_twdk') then
        call qqb_w_twdk(p,msq)
      elseif (case .eq. 'Wtbwdk') then
        call qqb_wtbwdk(p,msq)
      elseif (case .eq. 'Wtbndk') then
        call qqb_wtbndk(p,msq)
      elseif (case .eq. 'tottth') then
        call qqb_tottth(p,msq)
      elseif (case .eq. 'qq_tth') then
        call qqb_tth(p,msq)
      elseif (case .eq. 'qq_ttz') then
        call qqb_ttz(p,msq)
      elseif (case .eq. 'qq_ttw') then
        call qqb_ttw(p,msq)
      elseif (case .eq. 'httjet') then
        call qqb_higgs(p,msq)
      elseif (case .eq. 'ggfus1') then
        call gg_hg(p,msq)
      elseif (case .eq. 'Hgagaj') then
        call gg_hgagag(p,msq)
      elseif (case .eq. 'HWWjet') then
        call gg_hWWg(p,msq)
      elseif (case .eq. 'HZZjet') then
        call gg_hZZg(p,msq)
      elseif (case .eq. 'HWW2jt') then
        call gg_hWWgg(p,msq)
      elseif (case .eq. 'HZZ2jt') then
        call gg_hZZgg(p,msq)
      elseif (case .eq. 'HWW3jt') then
        call gg_hWWggg(p,msq)
      elseif (case .eq. 'HZZ3jt') then
        call gg_hZZggg(p,msq)
      elseif (case .eq. 'attjet') then
        call qqb_higgs_odd(p,msq)
      elseif (case .eq. 'qq_Hqq') then
        call VV_hqq(p,msq)
      elseif (case .eq. 'qq_Hgg') then
        call VV_Hgaga(p,msq)
      elseif (case .eq. 'qqHqqg') then
        call VV_hqq_g(p,msq)
      elseif (case .eq. 'qq_HWW') then
        call VV_HWW(p,msq)
      elseif (case .eq. 'qq_HZZ') then
        call VV_HZZ(p,msq)
      elseif (case .eq. 'tautau') then
        call qqb_tautau(p,msq)
      elseif (case .eq. 'qg_tbq') then
        call qg_tbq(p,msq)
c--- Check of gvec routines
c      call checkgvec(+2,0,2,p,qg_tbq,qg_tbq_gvec)
c      call checkgvec(-1,0,2,p,qg_tbq,qg_tbq_gvec)
      elseif (case .eq. 'qgtbqq') then
        call qg_tbq_g(p,msq)
      elseif (case .eq. '4ftwdk') then
	call qg_tbqdk(p,msq)
      elseif (case .eq. '4ftjet') then
        call qg_tbqdk_g(p,msq)
      elseif (case .eq. 'qq_tbg') then
        call qq_tbg(p,msq)
c--- Check of gvec routines
c      call checkgvec(2,-1,5,p,qq_tbg,qq_tbg_gvec)
c      call checkgvec(-1,2,5,p,qq_tbg,qq_tbg_gvec)
      elseif (case .eq. 'qqtbgg') then
        call qq_tbg_g(p,msq)
      elseif (case .eq. 'epem3j') then
        call epem3j(p,msq)
      elseif (case .eq. 'gQ__ZQ') then
        call gQ_zQ(p,msq)
      elseif (case .eq. 'Zccmas') then
        call qqb_zccm(p,msq)
      elseif (case .eq. 'ggfus2') then
        call gg_hgg(p,msq)
      elseif (case .eq. 'gagajj') then
        call gg_hgg(p,msq)
      elseif (case .eq. 'ggfus3') then
        call gg_hggg(p,msq)
      elseif (case .eq. 'W_bjet') then
        call qqb_wbjet(p,msq)
      elseif (case .eq. 'Wcjetg') then
        call qqb_w_cjet_massless_g(p,msq)
      elseif (case .eq. 'Z_bjet') then
        call qqb_zbjet(p,msq)
      elseif (case .eq. 'Zbjetg') then
        call qqb_zbjet_g(p,msq)
      elseif (case .eq. 'vlchk2') then
        call qqb_vol(p,msq)
        flux=one/vol(W,2)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=1d0
        fx2(-1)=1d0
      elseif (case .eq. 'vlchk3') then
        call qqb_vol(p,msq)
        flux=one/vol(W,3)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=2d0
        fx2(-1)=2d0
      elseif (case .eq. 'vlchk4') then
        taumin=0.0001d0
        bbsqmax=W
        bbsqmin=0d0
        call qqb_vol(p,msq)
        flux=one/vol(W,4)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=4d0*xx(1)
        fx2(-1)=2d0/xx(2)
      elseif (case .eq. 'vlchk5') then
        taumin=0.0001d0
        bbsqmax=W
        bbsqmin=0d0
        call qqb_vol(p,msq)
        flux=one/vol(W,5)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=4d0
        fx2(-1)=4d0
      elseif (case .eq. 'vlchk6') then
        call qqb_vol(p,msq)
        flux=one/vol(W,6)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=6d0*xx(1)
        fx2(-1)=4d0/xx(2)
      elseif (case .eq. 'vlchk8') then
        call qqb_vol(p,msq)
        flux=one/vol(W,8)
        bbsqmax=W
        bbsqmin=0d0
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=6d0/xx(1)
        fx2(-1)=6d0/xx(2)
      elseif (case .eq. 'vlchkm') then
        taumin=0.0001d0
        bbsqmax=W
        bbsqmin=0d0
        call qqb_vol(p,msq)
        flux=one/vol_mass(mb,W)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=4d0*xx(1)
        fx2(-1)=2d0/xx(2)      
      elseif (case .eq. 'vlchm3') then
        taumin=(2d0*mt/sqrts)**2
        call qqb_vol(p,msq)
        flux=one/vol3_mass(mt,W)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=2d0
        fx2(-1)=2d0
      elseif ((case .eq. 'vlchwt') .or. (case .eq. 'vlchwn')
     .   .or. (case .eq. 'vlchwg') .or. (case .eq. 'vlchwh')) then
        taumin=0.0001d0
        call qqb_vol(p,msq)
        flux=one/vol_wt(W)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=1d0
        fx2(-1)=1d0
      else
        write(6,*) 'Unimplemented process in lowint : case=',case
        stop 
      endif
      
      
      do j=-1,1
      do k=-1,1
      xmsq_bypart(j,k)=0d0
      enddo
      enddo 


      currentPDF=0

c--- do not calculate the flux if we're only checking the volume      
      if (case(1:4) .ne. 'vlch') then      
        flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
c--- for mlm study, divide by (Ecm)**2=W
c	if (runstring(1:3) .eq. 'mlm') then
c	  flux=flux/W
c	endif
      endif
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0d0
      if (PDFerrors) then
        call InitPDF(currentPDF)
      endif
      
c--- calculate PDF's  
      if (((case .eq. 'qg_tbq') .or. (case .eq. '4ftwdk'))
     &     .and. (dynamicscale .eqv. .false.)) then
c--- for single top + b, make sure to use two different scales
        call fdist(ih1,xx(1),facscale_H,fx1_H)
        call fdist(ih2,xx(2),facscale_H,fx2_H)
        call fdist(ih1,xx(1),facscale_L,fx1_L)
        call fdist(ih2,xx(2),facscale_L,fx2_L)
	do j=-nf,nf
	  if (j .eq. 0) then   ! heavy quark line has gluon init. state
	    fx1(j)=fx1_H(j)
	    fx2(j)=fx2_H(j)
	  else
	    fx1(j)=fx1_L(j)
	    fx2(j)=fx2_L(j)
	  endif
	enddo
      elseif (case(1:4) .ne. 'vlch') then
cc--- for comparison with C. Oleari's e+e- --> QQbg calculation
c        if (runstring(1:5) .eq. 'carlo') then
c          flux=1d0/2d0/W/(as/twopi)
cc--- divide out by (ason2pi) and then the "LO" massless DY process
c	  flux=flux/(aveqq*xn*fourpi*(gwsq/fourpi)**2/3d0/sqrts**2)
c	  flux=flux/(xn/8d0)
c	  do j=-nf,nf
c	  fx1(j)=0d0
c	  fx2(j)=0d0
c	  enddo
c	  fx1(0)=1d0
c	  fx2(0)=1d0
c        endif
        if ((case .eq. 'bq_tpq') .or. (case .eq. 'qg_tbq')) then   
c--- single top: allow for different scales on each leg  
c---  (applies only if dynstring = 'DDIS')
          if (dynstring .eq. 'DDIS') then
            call fdist(ih1,xx(1),b1scale,fxb1)
            call fdist(ih2,xx(2),q2scale,fx2)
            call fdist(ih1,xx(1),q1scale,fx1)
            call fdist(ih2,xx(2),b2scale,fxb2)
	  else          
            call fdist(ih1,xx(1),facscale,fx1)
            call fdist(ih2,xx(2),facscale,fx2)
	  endif
        else   
c--- usual case            
          call fdist(ih1,xx(1),facscale,fx1)
          call fdist(ih2,xx(2),facscale,fx2)
        endif
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav    

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      if (omitgg) then 
      if ((j.eq.0) .and. (k.eq.0)) goto 20
      endif

      if     ((case .eq. 'bq_tpq') .and. (dynstring .eq. 'DDIS')) then
c--- special case for dynamic scale in t-channel single top
        if     (abs(j) .eq. 5) then
          xmsqjk=fxb1(j)*fx2(k)*msq(j,k)
	elseif (abs(k) .eq. 5) then
          xmsqjk=fx1(j)*fxb2(k)*msq(j,k)
	else
	  xmsqjk=0d0
	endif
      elseif ((case .eq. 'qg_tbq') .and. (dynstring .eq. 'DDIS')) then
c--- special case for dynamic scale in t-channel single top
        if     (j .eq. 0) then
          xmsqjk=fxb1(j)*fx2(k)*msq(j,k)
	elseif (k .eq. 0) then
          xmsqjk=fx1(j)*fxb2(k)*msq(j,k)
	else
	  xmsqjk=0d0
	endif
      else
c--- DEFAULT
        xmsqjk=fx1(j)*fx2(k)*msq(j,k)
      endif

      xmsq=xmsq+xmsqjk
      
      if     (j .gt. 0) then
        sgnj=+1
      elseif (j .lt. 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k .gt. 0) then
        sgnk=+1
      elseif (k .lt. 0) then
        sgnk=-1
      else
        sgnk=0
      endif

      if (currentPDF .eq. 0) then
        xmsq_bypart(sgnj,sgnk)=xmsq_bypart(sgnj,sgnk)+xmsqjk
      endif
      
 20   continue
      enddo
      enddo

      if (currentPDF .eq. 0) then
        lowint=flux*pswt*xmsq/BrnRat
      endif
            
c--- loop over all PDF error sets, if necessary
      if (PDFerrors) then
        PDFwgt(currentPDF)=flux*pswt*xmsq/BrnRat*wgt/itmx
        PDFxsec(currentPDF)=PDFxsec(currentPDF)
     .     +PDFwgt(currentPDF)
        currentPDF=currentPDF+1
        if (currentPDF .le. maxPDFsets) goto 777
      endif    

      if (creatent) then
        wt_gg=xmsq_bypart(0,0)*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_gq=(xmsq_bypart(+1,0)+xmsq_bypart(-1,0)
     .        +xmsq_bypart(0,+1)+xmsq_bypart(0,-1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_qq=(xmsq_bypart(+1,+1)+xmsq_bypart(-1,-1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_qqb=(xmsq_bypart(+1,-1)+xmsq_bypart(-1,+1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
      endif

      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

      do j=-1,1
      do k=-1,1
        lord_bypart(j,k)=lord_bypart(j,k)+
     .       wgt*flux*pswt*xmsq_bypart(j,k)/BrnRat
      enddo
      enddo

      val=lowint*wgt
      val2=lowint**2*wgt
c--- update the maximum weight so far, if necessary
c---  but not if we are already unweighting ...
      if ((.not.unweight) .and. (dabs(val) .gt. wtmax)) then
        wtmax=dabs(val)
      endif

c      if(rescale) then 
c         call rescale_pjet(pjet)
c      endif

      if (bin) then
c ---   DSW. If the user has not selected to generate
c ---   events, still call nplotter here in order to
c ---   fill histograms/ntuples with weighted events :
        if (.not.evtgen) then
          call nplotter(pjet,val,val2,0)
c--- POWHEG-style output if requested
	  if (writepwg) then
            call pwhgplotter(p,pjet,val,0)
	  endif
        endif
      endif

c --- Check weights :
      if (unweight) then
c       write(6,*) 'Test weight ',val,' against max ',wtmax
        wtabs = dabs(val)
        if (ran2() .lt. (wtabs/wtmax)) then
c         write(6,*) 'Keep event with weight',val
          if (wtabs.lt.wtmax) then
            newwt = 1d0
          else
            newwt = wtabs/wtmax
          endif
          if (newwt .gt. 1.0d0) then
            write(6,*) 'WARNING : lowint : event with |weight| > 1.',
     +            ' |weight| = ',newwt
          endif
c ---     just in case the weight was negative :
          newwt = newwt*dsign(1d0,val)
          call nplotter(pjet,newwt,newwt,0)
c ---     DSW. If I'm storing the event, I need to make a decision
c ---     about the flavours :
          call decide_flavour(pflav,pbarflav)
          call storeevent(pjet,newwt,pflav,pbarflav)
        endif
      endif

      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end


