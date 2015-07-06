      subroutine qqb_dirgam_frag(p,msq)
 !------- FRAGMENTATION CONTRIBUTION, Always fragment p3 (=j) , allow phase space integration to symmetrize. 

      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'frag.f'
      integer j,k,i
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     .  qqij_ij,aaij_ij,
     .  qaij_ij,aqij_ij,aqii_jj,qaii_jj,
     .  qqii_ii,aaii_ii,aqii_ii,qaii_ii,
     .  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_gg,gg_qa,qa_gg,ss,tt,uu,
     .  smalla,smallb,smallc,smalld,
     &  qg_gq,ag_ga,gq_qg,ga_ag,qqij_ji,aaij_ji,qaij_ji,aqij_ji
     &  ,aqii_ii_s,qaii_ii_s
      double precision D(0:5),d_sum_q,fsq
      common/D/D

      fsq=frag_scale**2
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=down 2=up ....
      do i=0,5
         D(i)=0d0
         if     (fragset .eq. 'BFGset_I') then
            call get_frag(z_frag,fsq,1,i,D(i))   
         elseif (fragset .eq. 'BFGsetII') then  
            call get_frag(z_frag,fsq,2,i,D(i))   
         elseif (fragset .eq. 'GdRG__LO') then 
            call GGdR_frag(z_frag,i,D(i),0) 
         else
            write(6,*) 'Unrecognized fragmentation set name: ',fragset
            stop        
         endif
      enddo

      

      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(tt,ss,uu)
      aqij_ij=fac*aveqq*smalla(uu,tt,ss)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aqii_ii=fac*aveqq*smallb(uu,tt,ss)
      qaii_ii=fac*aveqq*smallb(uu,tt,ss)

      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half

      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(tt,uu,ss)
      ga_ga=-fac*aveqg*smallc(tt,ss,uu)

      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)

      gg_gg=fac*avegg*smalld(ss,tt,uu)*half

c--- additional amplitudes
      qg_gq=-fac*aveqg*smallc(uu,ss,tt)
      ag_ga=-fac*aveqg*smallc(uu,tt,ss)
      gq_qg=-fac*aveqg*smallc(uu,tt,ss)
      ga_ag=-fac*aveqg*smallc(uu,ss,tt)

      aqii_ii_s=fac*aveqq*smallb(tt,uu,ss)
      qaii_ii_s=fac*aveqq*smallb(tt,uu,ss)
      qqij_ji=fac*aveqq*smalla(ss,uu,tt)
      aaij_ji=fac*aveqq*smalla(ss,uu,tt)     
      qaij_ji=fac*aveqq*smalla(ss,uu,tt)    
      aqij_ji=fac*aveqq*smalla(tt,uu,ss)
   
       do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qq      
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then

!---- p3 = j 
            msq(j,k)=two*qqii_ii*D(abs(j))
          else
!-----  p3 = j 
            msq(j,k)=qqij_ij*D(abs(j))+qqij_ji*D(k)
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
               
            msq(j,k)=qaii_ii*D(j)+qaii_ii_s*D(j)
     &            +two*qaii_jj*d_sum_q(j)+qa_gg*D(0)*two
          else
!-------  p3=j 
            msq(j,k)=qaij_ij*D(j)+qaij_ji*D(abs(k))
          endif

C--aa      
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
            msq(j,k)=two*aaii_ii*D(abs(j))
          else
            msq(j,k)=aaij_ij*D(abs(j))+aaij_ji*D(abs(k))
          endif

C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=(aqii_ii+aqii_ii_s)*D(abs(j))
     &            +two*(aqii_jj*d_sum_q(-j)+aq_gg*D(0))
          else
            msq(j,k)=aqij_ij*D(abs(j))+aqij_ji*D(k)
          endif

C--qg_qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qg_qg*D(j)+qg_gq*D(0)
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=ag_ag*D(abs(j))+ag_ga*D(0)
C--gq_gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=gq_qg*D(k)+gq_gq*D(0)
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=ga_ag*D(abs(k))+ga_ga*D(0)
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=gg_gg*D(0)+gg_qa*d_sum_q(0)
      endif

      enddo
      enddo


      return
      end


      double precision function d_sum_q(j)
!----- sum over i=1,5 i!=j   
      implicit none 
      double precision D(0:5) 
      common/D/D
      integer i,j
      d_sum_q=0d0
      do i=1,5 
         if(i.ne.j) then 
            d_sum_q=d_sum_q+D(i) 
         endif
      enddo 
      return 
      end


