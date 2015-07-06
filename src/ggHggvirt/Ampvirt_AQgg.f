      subroutine Ampvirt_AQgg(p1,p2,p3,p4,ab41,ba41,ab43,ba43)
      implicit none
      include 'constants.f'
      include 'zprods_com.f'
      integer p1,p2,p3,p4
      double complex ab41(2,2,2),ba41(2,2,2),ab43(2,2,2),ba43(2,2,2),
     . A41HAQggmppp,A41HAQggmpmm,A41HAQggmpmp,A41HAQggmppm,
     . A43HAQggmppp,A43HAQggmpmm,A43HAQggmpmp,A43HAQggmppm

c--- calculate all the A41 amplitudes      
      ab41(1,2,2)=A41HAQggmppp(p1,p2,p3,p4,za,zb)
      ab41(1,1,1)=A41HAQggmpmm(p1,p2,p3,p4,za,zb)
      ab41(1,1,2)=A41HAQggmpmp(p1,p2,p3,p4,za,zb)
      ab41(1,2,1)=A41HAQggmppm(p1,p2,p3,p4,za,zb)

      ba41(1,2,2)=A41HAQggmppp(p1,p2,p4,p3,za,zb)
      ba41(1,1,1)=A41HAQggmpmm(p1,p2,p4,p3,za,zb)
      ba41(1,1,2)=A41HAQggmppm(p1,p2,p4,p3,za,zb)
      ba41(1,2,1)=A41HAQggmpmp(p1,p2,p4,p3,za,zb)

c--- obtain opposite helicity amplitudes by parity
      ab41(2,1,1)=-ab41(1,2,2)
      ab41(2,2,2)=-ab41(1,1,1)
      ab41(2,2,1)=-ab41(1,1,2)
      ab41(2,1,2)=-ab41(1,2,1)

      ba41(2,1,1)=-ba41(1,2,2)
      ba41(2,2,2)=-ba41(1,1,1)
      ba41(2,2,1)=-ba41(1,1,2)
      ba41(2,1,2)=-ba41(1,2,1)
      

c--- calculate all the A43 amplitudes      
      ab43(1,2,2)=A43HAQggmppp(p1,p2,p3,p4,za,zb)
      ab43(1,1,1)=A43HAQggmpmm(p1,p2,p3,p4,za,zb)
      ab43(1,1,2)=A43HAQggmpmp(p1,p2,p3,p4,za,zb)
      ab43(1,2,1)=A43HAQggmppm(p1,p2,p3,p4,za,zb)

      ba43(1,2,2)=A43HAQggmppp(p1,p2,p4,p3,za,zb)
      ba43(1,1,1)=A43HAQggmpmm(p1,p2,p4,p3,za,zb)
      ba43(1,1,2)=A43HAQggmppm(p1,p2,p4,p3,za,zb)
      ba43(1,2,1)=A43HAQggmpmp(p1,p2,p4,p3,za,zb)

c--- obtain opposite helicity amplitudes by parity
      ab43(2,1,1)=-ab43(1,2,2)
      ab43(2,2,2)=-ab43(1,1,1)
      ab43(2,2,1)=-ab43(1,1,2)
      ab43(2,1,2)=-ab43(1,2,1)

      ba43(2,1,1)=-ba43(1,2,2)
      ba43(2,2,2)=-ba43(1,1,1)
      ba43(2,2,1)=-ba43(1,1,2)
      ba43(2,1,2)=-ba43(1,2,1)
      
      
      return
      end
