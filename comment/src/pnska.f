      real function PNSKA(srt)
      SAVE   
***********************************
       if(srt.gt.3.0)then
       pnska=0
       return
       endif
      ala=1.116
      aka=0.498
      ana=0.939
      asa=1.197
      t1=asa+aka      
      if(srt.le.t1) THEN
      Pnska=0
       return
      Endif
      IF(SRT.LT.1.9)SBB1=(0.7/0.218)*(SRT-T1)
      IF(SRT.GE.1.9)SBB1=0.14/(SRT-1.7)
      sbb2=0.
       if(srt.gT.1.682)sbb2=0.5*(1.-0.75*(srt-1.682))
       pnska=0.25*(sbb1+sbb2)
* give the cross section in fm**2
       pnska=pnska/10.
      return
      end
********************************
*
*       Kaon momentum distribution in baryon-baryon-->N lamda K process
*
*       NOTE: dsima/dp is prototional to (1-p/p_max)(p/p_max)^2
*              we use rejection method to generate kaon momentum
*
*       Variables: Fkaon = F(p)/F_max
*                 srt   = cms energy of the colliding pair, 
*                          used to calculate the P_max
*       Date: Feb. 8, 1994
*
*       Reference: C. M. Ko et al.  
******************************** 
