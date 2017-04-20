      SUBROUTINE opioet(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
      PARAMETER (MAXSTR=150001)
      PARAMETER (AP1=0.13496,AP2=0.13957,ETAM=0.5475,aomega=0.782)
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
cc      SAVE /ppb1/
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
cc      SAVE /ppmm/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
      if((lb(i1).ge.3.and.lb(i1).le.5.and.lb(i2).eq.28).or.
     1     (lb(i2).ge.3.and.lb(i2).le.5.and.lb(i1).eq.28)) then
         iblock=1890
         ei1=aomega
         ei2=etam
c     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho)
c     thus the cross sections used are considered as the isospin-averaged ones.
         lbb1=28
         lbb2=0
      elseif((lb(i1).eq.28.and.lb(i2).eq.0).or.
     1        (lb(i1).eq.0.and.lb(i2).eq.28)) then
         iblock=1891
         lbb1=28
         lbb2=3+int(3*RANART(NSEED))
         ei1=aomega
         ei2=ap2
         if(lbb2.eq.4) ei2=ap1
      endif
      return
      END
*****************************************
* for rho rho <-> eta eta cross sections
