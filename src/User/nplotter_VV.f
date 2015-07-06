      subroutine nplotter_VV(p,wt,wt2,switch,nd)
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c--- switch:  an integer equal to 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      implicit none
      include 'vegas_common.f'
      include 'constants.f'
      include 'histo.f'
      include 'jetlabel.f'
      include 'frag.f'
      include 'masses.f'
      double precision p(mxpart,4),wt,wt2
      double precision etarap,pt,r
      double precision y3,y4,y5,y6,pt3,pt4,pt5,pt6,re5,rea5,re6,rea6
      double precision ylep, yjet, ptlep, ptjet,ptphot,yWgam,pt34,pttwo
      double precision pt45,pt56,m45,mt45
      double precision yraptwo,m345,mtrans,yellgam
      double precision et_vec(4),etmiss,r2,delphi,met,m34,m56,m3456
      integer switch,n,nplotmax,nproc,nqcdjets,nqcdstart,nd 
      character*4 tag
      integer j,i
      logical first,creatent,dswhisto
      common/outputflags/creatent,dswhisto
      common/nplotmax/nplotmax
      common/nproc/nproc
      common/nqcdjets/nqcdjets,nqcdstart
      data first/.true./
      save first
  
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************


      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag='book'
	y3=1d3
        y4=1d3
        pt3=1d7
        pt4=1d7
        pt34=1d7
c---Intiailise photon 
        y5=1d3
        pt5=1d7
        mt45=1d7
        m45=1d7
        pt45=1d7
        pt34=1d7
        pt56=1d7
c----Initialise jet values will not pass cuts in there is an NLO jet
        y6=1d3
        pt6=1d7
c--- If re5 is not changed by the NLO value, it will be out of
c--- the plotting range
        re5=1d3
        rea5=1d3
        re6=1d3
        rea6=1d3
        jets=nqcdjets
c--- (Upper) limits for the plots
        ylep=5d0
        yjet=5d0
        ptlep=100d0
        ptjet=200d0
        ptphot=500d0
        goto 99
      else
c--- Add event in histograms
        tag='plot'
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
************************************************************************

!     121 '  f(p1)+f(p2) --> H(--> W^+(nu(p3)+e^+(p4)) + W^-(e^-(p5)+nu~(p6))) [top, bottom loops, exact]' 'L'
!     122 '  f(p1)+f(p2) --> H(--> W^+(nu(p3)+e^+(p4)) + W^-(e^-(p5)+nu~(p6))) [above + interf. with gg->WW]' 'L'


      
         pt3=pt(3,p)
         pt4=pt(4,p) 
         pt5=pt(5,p)
         pt6=pt(6,p) 
         y4=etarap(4,p)
         y5=etarap(5,p) 
         pt45=pttwo(4,5,p)
         pt34=pttwo(3,4,p) 
         pt56=pttwo(5,6,p)
         met=etmiss(p,et_vec) 
         r2= (p(4,1)*p(5,1)+p(4,2)*p(5,2))
     .        /dsqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
         if (r2 .gt. +0.9999999D0) r2=+1D0
         if (r2 .lt. -0.9999999D0) r2=-1D0
         delphi=dacos(r2)
!--- m_ll cut 
         m45=0d0
      do i=1,4
         if(i.ne.4) then 
            m45=m45-(p(4,i)+p(5,i))**2
         else
            m45=m45+(p(4,i)+p(5,i))**2 
         endif
      enddo
      m45=dsqrt(m45)
      m3456=0d0
      do i=1,4
         if(i.ne.4) then 
            m34=m34-(p(3,i)+p(4,i))**2
            m3456=m3456-(p(3,i)+p(4,i)+p(5,i)+p(6,i))**2
         else
            m34=m34+(p(4,i)+p(3,i))**2 
            m3456=m3456+(p(3,i)+p(4,i)+p(5,i)+p(6,i))**2
         endif
      enddo
      m34=dsqrt(m34)
      m3456=dsqrt(m3456)
       do i=1,4
         if(i.ne.4) then 
            m56=m56-(p(5,i)+p(6,i))**2
         else
            m56=m56+(p(5,i)+p(6,i))**2 
         endif
      enddo
      m56=dsqrt(m56)

      mt45=0d0 
      mt45=(dsqrt(dsqrt(pttwo(4,5,p)**2+m45**2)+etmiss(p,et_vec))**2)
 !    endif



************************************************************************
*                                                                      *
*     FILL HISTOGRAMS                                                  *
*                                                                      *
************************************************************************

c--- Call histogram routines
   99 continue

c--- Book and fill ntuple if that option is set, remembering to divide
c--- by # of iterations now that is handled at end for regular histograms
      if (creatent .eqv. .true.) then
        call bookfill(tag,p,wt/dfloat(itmx))  
      endif

c--- "n" will count the number of histograms
      n=1              

c--- Syntax of "bookplot" routine is:
c
c---   call bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot)
c
c---        n:  internal number of histogram
c---      tag:  "book" to initialize histogram, "plot" to fill
c---   titlex:  title of histogram
c---      var:  value of quantity being plotted
c---       wt:  weight of this event (passed in)
c---      wt2:  weight of this event (passed in)
c---     xmin:  lowest value to bin
c---     xmax:  highest value to bin
c---       dx:  bin width
c---   llplot:  equal to "lin"/"log" for linear/log scale
 
      call bookplot(n,tag,'pt_nu_1',pt3,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_e_1',pt4,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'y_e_1',y4,wt,wt2,-4d0,4d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_W_1',pt34,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_e_2',pt5,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'y_e_2',y5,wt,wt2,-4d0,4d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_nu_2',pt6,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_W_2',pt56,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_ll',pt45,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'m_ll',m45,wt,wt2,0d0,100d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'m_ll',m34,wt,wt2,0d0,250d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'m_ll_2',m56,wt,wt2,0d0,250d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'m_4l',m3456,wt,wt2,hmass-0.5d0,hmass+0.5d0,
     & 0.05d0,'lin')
      n=n+1
       call bookplot(n,tag,'m_4l',m3456,wt,wt2,0d0,2000d0,20d0
     &,'lin')
      n=n+1
      call bookplot(n,tag,'mt',mt45,wt,wt2,0d0,200d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'delphi',delphi,wt,wt2,0d0,3.14d0,0.1d0,'lin')
      n=n+1
  
************************************************************************
*                                                                      *
*     FINAL BOOKKEEPING                                                *
*                                                                      *
************************************************************************

c--- We have over-counted the number of histograms by 1 at this point
      n=n-1

c--- Ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(n)

c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
c
      return 
      end
