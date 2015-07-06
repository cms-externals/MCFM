      double precision function etmiss(p,etvec)
      implicit none
      include 'constants.f'
      include 'plabel.f'
      integer j,k
      double precision etvec(4),p(mxpart,4)
      
      do k=1,4
        etvec(k)=0d0
      enddo
      
      do j=1,mxpart
        if ((plabel(j) .eq. 'nl') .or. (plabel(j) .eq. 'na')) then
          do k=1,4
            etvec(k)=etvec(k)+p(j,k)
          enddo
        endif
      enddo
      
      etmiss=dsqrt(etvec(1)**2+etvec(2)**2)
      
      return
      end
      
