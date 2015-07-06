      subroutine Amplo_AQgg(p1,p2,p3,p4,ab,ba)
c--- this routine is a wrapper to the new versions of the leading
c--- order amplitudes that are calculated in the same way as the
c--- new virtual amplitudes
      implicit none
      include 'constants.f'
      include 'zprods_com.f'
      integer p1,p2,p3,p4
      double complex ab(2,2,2),ba(2,2,2),
     . A0HAQggmppp,A0HAQggmpmm,A0HAQggmpmp,A0HAQggmppm
      
      ab(1,2,2)=A0HAQggmppp(p1,p2,p3,p4,za,zb)
      ab(1,1,1)=A0HAQggmpmm(p1,p2,p3,p4,za,zb)
      ab(1,1,2)=A0HAQggmpmp(p1,p2,p3,p4,za,zb)
      ab(1,2,1)=A0HAQggmppm(p1,p2,p3,p4,za,zb)

      ba(1,2,2)=A0HAQggmppp(p1,p2,p4,p3,za,zb)
      ba(1,1,1)=A0HAQggmpmm(p1,p2,p4,p3,za,zb)
      ba(1,1,2)=A0HAQggmppm(p1,p2,p4,p3,za,zb)
      ba(1,2,1)=A0HAQggmpmp(p1,p2,p4,p3,za,zb)

c--- obtain opposite helicity amplitudes by parity
      ab(2,1,1)=-ab(1,2,2)
      ab(2,2,2)=-ab(1,1,1)
      ab(2,2,1)=-ab(1,1,2)
      ab(2,1,2)=-ab(1,2,1)

      ba(2,1,1)=-ba(1,2,2)
      ba(2,2,2)=-ba(1,1,1)
      ba(2,2,1)=-ba(1,1,2)
      ba(2,1,2)=-ba(1,2,1)
      
      return
      end
      
