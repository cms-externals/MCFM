      subroutine phase5(r,p1,p2,p3,p4,p5,p6,p7,wt)
c----phase space for signal
      implicit none
      include 'constants.f'
      include 'heavyflav.f'
      include 'masses.f'
      include 'process.f'
      include 'mxdim.f'
      include 'debug.f'
c********* generate phase space for 2-->5 process
c********* r(mxdim),p1(4),p2(4) are inputs 
c--------- incoming p1 and p2 reversed in sign from physical values 
c---- i.e. phase space for -p1-p2 --> p3+p4+p5+p6+p7
c---- with all 2 pi's (ie 1/(2*pi)^11)

      double precision r(mxdim)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      double precision p127(4),p12(4),p56(4),p34(4),smin
      double precision wt,wt127,wt3456,wt34,wt56,wt0
      integer j

      parameter(wt0=1d0/twopi**3)

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo


      smin=mb**2

c--- In the case of HVV_4l, we should generate s127 according to
c--- a Breit-Wigner at mH, otherwise just linearly      
      if (  (case .eq. 'HWW_4l') 
     . .or. (case .eq. 'HWW2lq')
     . .or. (case .eq. 'HWW_tb')
     . .or. (case .eq. 'HWWint')
     . .or. (case .eq. 'HZZ_4l')
     . .or. (case .eq. 'HZZ_tb')
     . .or. (case .eq. 'HZZint')
     . .or. (case .eq. 'HWWjet')
     . .or. (case .eq. 'HZZjet')
     . ) then
        call phi1_2m_bw(zip,r(13),r(12),r(11),smin,p12,p7,p127,
     .   hmass,hwidth,wt127,*99)
      else
        call phi1_2m_nobw(zip,r(13),r(12),r(11),
     .   smin,p12,p7,p127,wt127,*99)
      endif
      
      call phi1_2(r(1),r(2),r(3),r(4),p127,p56,p34,wt3456,*99)

      if (   (case .eq. 'Wbbjet') .or. (case .eq. 'Wbbjem')
     .  .or. ((case .eq. 'Wbbmas') .and. (flav .eq. 5))  
     .  .or. (case .eq. 'WHbbar')  
     .  .or. (case .eq. 'ZHbbar')  
     .  .or. ((case .eq. 'W_bjet') .and. (mb .gt. 0d0)) ) then
        call phi3m(r(5),r(6),p56,p5,p6,mb,mb,wt56,*99)
      elseif ((case .eq. 'Wbbmas') .and. (flav .eq. 4)) then
        call phi3m(r(5),r(6),p56,p5,p6,mc,mc,wt56,*99)
      elseif ((case .eq. 'Wttmas') .or. (case .eq. 'qq_ttw')) then
        call phi3m(r(5),r(6),p56,p5,p6,mt,mt,wt56,*99)
      else
        call phi3m0(r(5),r(6),p56,p5,p6,wt56,*99)
      endif
      call phi3m0(r(7),r(8),p34,p3,p4,wt34,*99)
      wt=wt0*wt127*wt3456*wt56*wt34

  
      if (debug) write(6,*) 'wt127',wt127
      if (debug) write(6,*) 'wt3456',wt3456
      if (debug) write(6,*) 'wt34',wt34
      if (debug) write(6,*) 'wt56',wt56

      return
 99   continue
      wt=0d0
      return
      end

