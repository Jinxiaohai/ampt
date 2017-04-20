       subroutine Rmasdd(srt,am10,am20,
     &dmin1,dmin2,ISEED,ic,dm1,dm2)
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
       amn=0.94
       amp=0.14
* the maximum mass for resonance 1
         dmax1=srt-dmin2
* generate the mass for the first resonance
 5        NTRY1=0
         ntry2=0
         ntry=0
         ictrl=0
10        DM1 = RANART(NSEED) * (DMAX1-DMIN1) + DMIN1
          NTRY1=NTRY1+1
* the maximum mass for resonance 2 
         if(ictrl.eq.0)dmax2=srt-dm1
* generate the mass for the second resonance
20         dm2=RANART(NSEED)*(dmax2-dmin2)+dmin2
          NTRY2=NTRY2+1
* check the energy-momentum conservation with two masses
* q2 in the following is q**2*4*srt**2
         q2=((srt**2-dm1**2-dm2**2)**2-4.*dm1**2*dm2**2)
         if(q2.le.0)then
         dmax2=dm2-0.01
c         dmax1=dm1-0.01
         ictrl=1
         go to 20
         endif
* determine the weight of the mass pair         
          IF(DMAX1.LT.am10) THEN
          if(ic.eq.1)FM1=Fmassd(DMAX1)
          if(ic.eq.2)FM1=Fmassn(DMAX1)
          if(ic.eq.3)FM1=Fmassd(DMAX1)
          if(ic.eq.4)FM1=Fmassd(DMAX1)
          ELSE
          if(ic.eq.1)FM1=Fmassd(am10)
          if(ic.eq.2)FM1=Fmassn(am10)
          if(ic.eq.3)FM1=Fmassd(am10)
          if(ic.eq.4)FM1=Fmassd(am10)
          ENDIF
          IF(DMAX2.LT.am20) THEN
          if(ic.eq.1)FM2=Fmassd(DMAX2)
          if(ic.eq.2)FM2=Fmassn(DMAX2)
          if(ic.eq.3)FM2=Fmassn(DMAX2)
          if(ic.eq.4)FM2=Fmassr(DMAX2)
          ELSE
          if(ic.eq.1)FM2=Fmassd(am20)
          if(ic.eq.2)FM2=Fmassn(am20)
          if(ic.eq.3)FM2=Fmassn(am20)
          if(ic.eq.4)FM2=Fmassr(am20)
          ENDIF
          IF(FM1.EQ.0.)FM1=1.e-04
          IF(FM2.EQ.0.)FM2=1.e-04
         prob0=fm1*fm2
          if(ic.eq.1)prob=Fmassd(dm1)*fmassd(dm2)
          if(ic.eq.2)prob=Fmassn(dm1)*fmassn(dm2)
          if(ic.eq.3)prob=Fmassd(dm1)*fmassn(dm2)
          if(ic.eq.4)prob=Fmassd(dm1)*fmassr(dm2)
         if(prob.le.1.e-06)prob=1.e-06
         fff=prob/prob0
         ntry=ntry+1 
          IF(RANART(NSEED).GT.fff.AND.
     1    NTRY.LE.20) GO TO 10
clin-2/26/03 limit the mass of (rho,Delta,N*1440) below a certain value
c     (here taken as its central value + 2* B-W fullwidth):
          if((abs(am10-0.77).le.0.01.and.dm1.gt.1.07)
     1         .or.(abs(am10-1.232).le.0.01.and.dm1.gt.1.47)
     2         .or.(abs(am10-1.44).le.0.01.and.dm1.gt.2.14)) goto 5
          if((abs(am20-0.77).le.0.01.and.dm2.gt.1.07)
     1         .or.(abs(am20-1.232).le.0.01.and.dm2.gt.1.47)
     2         .or.(abs(am20-1.44).le.0.01.and.dm2.gt.2.14)) goto 5
       RETURN
       END
*FUNCTION Fmassd(DMASS) GIVES the delta MASS DISTRIBUTION 
