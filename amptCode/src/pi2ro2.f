      SUBROUTINE pi2ro2(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
      PARAMETER (MAXSTR=150001)
      PARAMETER (AP1=0.13496,AP2=0.13957)
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      COMMON/RNDF77/NSEED
      SAVE   
      if((lb(i1).ge.3.and.lb(i1).le.5)
     1     .and.(lb(i2).ge.3.and.lb(i2).le.5)) then
         iblock=1850
         ei1=0.77
         ei2=0.77
         lbb1=25+int(3*RANART(NSEED))
         lbb2=25+int(3*RANART(NSEED))
      elseif((lb(i1).ge.25.and.lb(i1).le.27)
     1     .and.(lb(i2).ge.25.and.lb(i2).le.27)) then
         iblock=1851
         lbb1=3+int(3*RANART(NSEED))
         lbb2=3+int(3*RANART(NSEED))
         ei1=ap2
         ei2=ap2
         if(lbb1.eq.4) ei1=ap1
         if(lbb2.eq.4) ei2=ap1
      endif
      return
      END
