      SUBROUTINE ro2et2(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
      PARAMETER (MAXSTR=150001)
      parameter (ETAM=0.5475,arho=0.77)
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
cc      SAVE /ppb1/
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
cc      SAVE /ppmm/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
      if(lb(i1).ge.25.and.lb(i1).le.27.and.
     1     lb(i2).ge.25.and.lb(i2).le.27) then
         iblock=1895
         ei1=etam
         ei2=etam
c     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho)
c     thus the cross sections used are considered as the isospin-averaged ones.
         lbb1=0
         lbb2=0
      elseif(lb(i1).eq.0.and.lb(i2).eq.0) then
         iblock=1896
         lbb1=25+int(3*RANART(NSEED))
         lbb2=25+int(3*RANART(NSEED))
         ei1=arho
         ei2=arho
      endif
      return
      END
*****************************
* purpose: Xsection for K* Kbar or K*bar K to pi(eta) rho(omega)
