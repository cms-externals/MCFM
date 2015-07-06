*********************************************************
*AUTHOR: FABIO MALTONI                                  *
*DATE  : 3/18/2002                                      *
*NOTES : PROGRAM GENERATED BY WRITE_VIRTUALS.M          *
*********************************************************

      DOUBLE PRECISION FUNCTION BBGHVIRT(I1,I2,I3)                                    
* ---------------------------------------------------------------------
* returns the interference between 1-loop amplitude and tree amplitude 
*                                                                      
* 2 Re(AMP_LO AMP_VIRT^* )                                             
*                                                                      
* for the process:                                                     
*          0 -> bbar1 b2 g3 h                                          
*                                                                      
* All momenta outgoing.                                                
* No averaging is performed for initial spins or colors.               
* Higgs-bottom and strong couplings are included.                      
* Virtual contributions calculated in DRED.                            
* UV divergences are subtracted here in the MSbar scheme.              
* IR divergences are present as epinv^2 and epinv terms.               
* If you want to use this subroutine for b bbar -> h g or h-> b bbar g 
* you only have to modifiy the sign of the born and of 'term'.         
* ---------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'constants.f'
      INCLUDE 'epinv.f'
      include 'epinv2.f'
      include 'sprods_com.f'
      INCLUDE 'qcdcouple.f'
      INCLUDE 'scale.f'
      INCLUDE 'b0.f'

*     ARGUMENTS 
      INTEGER I1,I2,I3

*     EXTERNAL 
      DOUBLE PRECISION ddilog

*     LOCAL 
      DOUBLE PRECISION SS,TT,UU,SN,TN,UN,MH2
      DOUBLE COMPLEX LNRAT
      DOUBLE COMPLEX xlns,xlnt,xlnu,xlnsm,xlntm,xlnum
      DOUBLE COMPLEX xlnmh
      DOUBLE COMPLEX xln1ms,xln1mt,xln1mu
      DOUBLE PRECISION dilogs,dilogt,dilogu
      DOUBLE PRECISION c0t,c0u,c0ms,c0mt,c0mu,d0st,d0su,d0tu
      DOUBLE PRECISION born,term

C     coupsq=susycoup**2*ghbb**2  => as in the Born function

      DOUBLE PRECISION coupsq
      coupsq=1d0   !change this



      mh2=s(i1,i2)+s(i1,i3)+s(i2,i3)
      ss =s(i1,i2)
      tt =s(i1,i3)
      uu =s(i2,i3)
      sn =s(i1,i2)/mh2
      tn =s(i1,i3)/mh2
      un =s(i2,i3)/mh2

* compute the logs and the dilogs
* the analytic continuation is performed by lnrat

               xlns =lnrat(-ss,-mh2)                           
               xlnt =lnrat(-tt,-mh2)                           

               xlnu =lnrat(-uu,-mh2)                           
               xlnmh=lnrat(-mh2,musq)                         
               xlnsm=lnrat(-ss ,musq)                         
               xlntm=lnrat(-tt ,musq)                         
               xlnum=lnrat(-uu ,musq)                         
               xln1ms=lnrat(-mh2+ss,-mh2)                      
               xln1mt=lnrat(-mh2+tt,-mh2)                      
               xln1mu=lnrat(-mh2+uu,-mh2)                      
               dilogs=ddilog(sn)                              
               dilogt=ddilog(tn)                              
               dilogu=ddilog(un)                              

* compute the scalar integrals 
* the following integrals are dimensionless,   
* and are defined as in the tex notes, 
* and have a 1/(16 pi^2)*c_\Gamma understood in front 
* where c_\Gamma= (4 pi)^eps Gamma(1+eps) Gamma(1-eps)^2/Gamma(1-2eps)

           c0t= dble(epinv2 - epinv*xlntm + xlntm**2/2d0)
           c0u= dble(epinv2 - epinv*xlnum + xlnum**2/2d0)
           c0ms=dble(epinv*xlns + (xlnmh**2 - xlnsm**2)/2d0)
           c0mt=dble(epinv*xlnt + (xlnmh**2 - xlntm**2)/2d0)
           c0mu=dble(epinv*xlnu + (xlnmh**2 - xlnum**2)/2d0)
           d0su=dble(
     .    2d0*(-dilogs - dilogu + epinv2 + pisq/6d0 - epinv*xlnmh + 
     -    xlnmh**2/2d0 - epinv*xlns - xln1ms*xlns + xlnmh*xlns + 
     -    xlns**2/2d0 - epinv*xlnu - xln1mu*xlnu + xlnmh*xlnu + 
     -    xlns*xlnu + xlnu**2/2d0))
           d0st=dble(
     .    2d0*(-dilogs - dilogt + epinv2 + pisq/6d0 - epinv*xlnmh + 
     -    xlnmh**2/2d0 - epinv*xlns - xln1ms*xlns + xlnmh*xlns + 
     -    xlns**2/2d0 - epinv*xlnt - xln1mt*xlnt + xlnmh*xlnt + 
     -    xlns*xlnt + xlnt**2/2d0))
           d0tu=dble(
     .    2d0*(-dilogt - dilogu + epinv2 + pisq/6d0 - epinv*xlnmh + 
     -    xlnmh**2/2d0 - epinv*xlnt - xln1mt*xlnt + xlnmh*xlnt + 
     -    xlnt**2/2d0 - epinv*xlnu - xln1mu*xlnu + xlnmh*xlnu + 
     -    xlnt*xlnu + xlnu**2/2d0))

* the sign of born and term is valid only for g b->b h 

           born  =  -(16D0*(1d0+sn**2)/tn/un)*gsq*coupsq
           term  =  -(16D0*(1d0-sn)   /tn/un)*gsq*coupsq

* final result 

           BBGHVIRT=(CA-CF)*term+
     .        BORN*(CF*(-2d0*c0ms - c0mt - c0mu + c0t + c0u - d0st - 
     -       d0su) + (CA*(2d0*c0ms - 2d0*c0t - 2d0*c0u + d0st + d0su - 
     -         d0tu))/2d0)

           BBGHVIRT=ason2pi*BBGHVIRT

* UV renormalization of the strong coupling in DRED-MSbar 
* and finite renormalization to standard MSbar

           BBGHVIRT=BBGHVIRT+ason2pi*born*(-(epinv)*b0+CA/6d0)

* UV renormalization of the Yukawa coupling in DRED-MSbar
* and finite renormalization to standard MSbar

           BBGHVIRT=BBGHVIRT+ason2pi*born*(-(epinv)*3d0*CF-CF)

         END
         
