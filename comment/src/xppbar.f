      real function xppbar(srt)
*  srt    = DSQRT(s) in GeV                                                   *
*  xppbar = pp_bar annihilation cross section in mb                           *
*                                                    
*  Reference: G.J. Wang, R. Bellwied, C. Pruneau and G. Welke
*             Proc. of the 14th Winter Workshop on Nuclear Dynamics, 
*             Snowbird, Utah 31, Eds. W. Bauer and H.G. Ritter 
*             (Plenum Publishing, 1998)                             *
*
******************************************
       Parameter (pmass=0.9383,xmax=400.)
      SAVE   
* Note:
* (1) we introduce a new parameter xmax=400 mb:
*     the maximum annihilation xsection 
* there are shadowing effects in pp_bar annihilation, with this parameter
* we can probably look at these effects  
* (2) Calculate p(lab) from srt [GeV], since the formular in the 
* reference applies only to the case of a p_bar on a proton at rest
* Formula used: srt**2=2.*pmass*(pmass+sqrt(pmass**2+plab**2))
       xppbar=1.e-06
       plab2=(srt**2/(2.*pmass)-pmass)**2-pmass**2
       if(plab2.gt.0)then
           plab=sqrt(plab2)
       xppbar=67./(plab**0.7)
       if(xppbar.gt.xmax)xppbar=xmax
       endif
         return
      END
cbali1/16/99 end
**********************************
cbali2/6/99
********************************************
* Purpose: To generate randomly the no. of pions in the final 
*          state of pp_bar annihilation according to a statistical 
*          model by using of the rejection method.  
cbz2/25/99
c      real*4 function pbarfs(srt,npion,iseed)
