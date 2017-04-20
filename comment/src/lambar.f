         subroutine lambar(i1,i2,srt,siglab)
*  srt    = DSQRT(s) in GeV                                               *
*  siglab = lambda-nuclar elastic cross section in mb 
*         = 12 + 0.43/p_lab**3.3 (mb)  
*                                                    
* (2) Calculate p(lab) from srt [GeV], since the formular in the 
* reference applies only to the case of a p_bar on a proton at rest
* Formula used: srt**2=2.*pmass*(pmass+sqrt(pmass**2+plab**2))
*****************************
        PARAMETER (MAXSTR=150001)
        COMMON /AA/ R(3,MAXSTR)
cc      SAVE /AA/
        COMMON /BB/ P(3,MAXSTR)
cc      SAVE /BB/
        COMMON /CC/ E(MAXSTR)
cc      SAVE /CC/
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      SAVE   
          siglab=1.e-06
        if( iabs(lb(i1)).ge.14.and.iabs(lb(i1)).le.17 )then
          eml = e(i1)
          emb = e(i2)
         else
          eml = e(i2)
          emb = e(i1)
        endif
       pthr = srt**2-eml**2-emb**2
        if(pthr .gt. 0.)then
       plab2=(pthr/2./emb)**2-eml**2
       if(plab2.gt.0)then
         plab=sqrt(plab2)
         siglab=12. + 0.43/(plab**3.3)
       if(siglab.gt.200.)siglab=200.
       endif
       endif
         return
      END
C------------------------------------------------------------------
clin-7/26/03 improve speed
***************************************
