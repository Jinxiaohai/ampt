      SUBROUTINE CRPN(PX,PY,PZ,SRT,I1,I2,
     & IBLOCK,xkaon0,xkaon,Xphi,xphin)
*     PURPOSE:                                                         *
*           DEALING WITH PION+N-->L/S+KAON PROCESS AND PION PRODUCTION *
*     NOTE   :                                                         *
*          
*     QUANTITIES:                                                 *
*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
*           SRT      - SQRT OF S                                       *
*           IBLOCK   - THE INFORMATION BACK                            *
*                     7  PION+N-->L/S+KAON
*           iblock   - 77 pion+N-->Delta+pion
*           iblock   - 78 pion+N-->Delta+RHO
*           iblock   - 79 pion+N-->Delta+OMEGA
*           iblock   - 222 pion+N-->Phi 
**********************************
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,APHI=1.020,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974)
        parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
        COMMON /AA/ R(3,MAXSTR)
cc      SAVE /AA/
        COMMON /BB/ P(3,MAXSTR)
cc      SAVE /BB/
        COMMON /CC/ E(MAXSTR)
cc      SAVE /CC/
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
      PX0=PX
      PY0=PY
      PZ0=PZ
      iblock=1
      x1=RANART(NSEED)
      ianti=0
      if(lb(i1).lt.0 .or. lb(i2).lt.0) ianti=1
      if(xkaon0/(xkaon+Xphi).ge.x1)then
* kaon production
*-----------------------------------------------------------------------
        IBLOCK=7
        if(ianti .eq. 1)iblock=-7
        NTAG=0
* RELABLE PARTICLES FOR THE PROCESS PION+n-->LAMBDA K OR SIGMA k
* DECIDE LAMBDA OR SIGMA PRODUCTION, AND TO CALCULATE THE NEW
* MOMENTA FOR PARTICLES IN THE FINAL STATE.
       KAONC=0
       IF(PNLKA(SRT)/(PNLKA(SRT)
     &       +PNSKA(SRT)).GT.RANART(NSEED))KAONC=1
       IF(E(I1).LE.0.2)THEN
           LB(I1)=23
           E(I1)=AKA
           IF(KAONC.EQ.1)THEN
              LB(I2)=14
              E(I2)=ALA
           ELSE
              LB(I2) = 15 + int(3 * RANART(NSEED))
              E(I2)=ASA       
           ENDIF
           if(ianti .eq. 1)then
              lb(i1) = 21
              lb(i2) = -lb(i2)
           endif
       ELSE
           LB(I2)=23
           E(I2)=AKA
           IF(KAONC.EQ.1)THEN
              LB(I1)=14
              E(I1)=ALA
           ELSE
              LB(I1) = 15 + int(3 * RANART(NSEED))
              E(I1)=ASA       
           ENDIF
           if(ianti .eq. 1)then
              lb(i2) = 21
              lb(i1) = -lb(i1)
           endif
       ENDIF
        EM1=E(I1)
        EM2=E(I2)
        go to 50
* to gererate the momentum for the kaon and L/S
      elseif(Xphi/(xkaon+Xphi).ge.x1)then
          iblock=222
         if(xphin/Xphi .ge. RANART(NSEED))then
          LB(I1)= 1+int(2*RANART(NSEED))
           E(I1)=AMN
         else
          LB(I1)= 6+int(4*RANART(NSEED))
           E(I1)=AM0
         endif
c  !! at present only baryon
         if(ianti .eq. 1)lb(i1)=-lb(i1)
          LB(I2)= 29
           E(I2)=APHI
        EM1=E(I1)
        EM2=E(I2)
       go to 50
         else
* CHECK WHAT KIND OF PION PRODUCTION PROCESS HAS HAPPENED
       IF(RANART(NSEED).LE.TWOPI(SRT)/
     &  (TWOPI(SRT)+THREPI(SRT)+FOURPI(SRT)))THEN
       iblock=77
       ELSE 
        IF(THREPI(SRT)/(THREPI(SRT)+FOURPI(SRT)).
     &  GT.RANART(NSEED))THEN
       IBLOCK=78
       ELSE
       IBLOCK=79
       ENDIF
       endif
       ntag=0
* pion production (Delta+pion/rho/omega in the final state)
* generate the mass of the delta resonance
       X2=RANART(NSEED)
* relable the particles
       if(iblock.eq.77)then
* GENERATE THE DELTA MASS
       dmax=srt-ap1-0.02
       dm=rmass(dmax,iseed)
* pion+baryon-->pion+delta
* Relable particles, I1 is assigned to the Delta and I2 is assigned to the
* meson
*(1) for pi(+)+p-->D(+)+P(+) OR D(++)+p(0)
       if( ((lb(i1).eq.1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.1))
     &       .OR. ((lb(i1).eq.-1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=8
       e(i1)=dm
       lb(i2)=5
       e(i2)=ap1
       go to 40
       ELSE
       lb(i1)=9
       e(i1)=dm
       lb(i2)=4
        ipi = 4
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=8
       e(i2)=dm
       lb(i1)=5
       e(i1)=ap1
       go to 40
       ELSE
       lb(i2)=9
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
*(2) for pi(-)+p-->D(0)+P(0) OR D(+)+p(-),or D(-)+p(+)
       if( ((lb(i1).eq.1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.1))
     &        .OR. ((lb(i1).eq.-1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=6
       e(i1)=dm
       lb(i2)=5
       e(i2)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=7
       e(i1)=dm
       lb(i2)=4
       e(i2)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=8
       e(i1)=dm
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=6
       e(i2)=dm
       lb(i1)=5
       e(i1)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=7
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=8
       e(i2)=dm
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
*(3) for pi(+)+n-->D(+)+Pi(0) OR D(++)+p(-) or D(0)+pi(+)
       if( ((lb(i1).eq.2.and.lb(i2).eq.5).
     &   or.(lb(i1).eq.5.and.lb(i2).eq.2))
     & .OR. ((lb(i1).eq.-2.and.lb(i2).eq.3).
     &   or.(lb(i1).eq.3.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=8
       e(i1)=dm
       lb(i2)=4
       e(i2)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=7
       e(i1)=dm
       lb(i2)=5
       e(i2)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=9
       e(i1)=dm
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=8
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=7
       e(i2)=dm
       lb(i1)=5
       e(i1)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=9
       e(i2)=dm
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
*(4) for pi(0)+p-->D(+)+Pi(0) OR D(++)+p(-) or D(0)+pi(+)
       if((iabs(lb(i1)).eq.1.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.1))then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=8
       e(i1)=dm
       lb(i2)=4
       e(i2)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=7
       e(i1)=dm
       lb(i2)=5
       e(i2)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=9
       e(i1)=dm
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=8
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=7
       e(i2)=dm
       lb(i1)=5
       e(i1)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=9
       e(i2)=dm
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       endif 
*(5) for pi(-)+n-->D(-)+P(0) OR D(0)+p(-)
       if( ((lb(i1).eq.2.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.2))
     &         .OR. ((lb(i1).eq.-2.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=6
       e(i1)=dm
       lb(i2)=4
       e(i2)=ap1
       go to 40
       ELSE
       lb(i1)=7
       e(i1)=dm
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=6
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       ELSE
       lb(i2)=7
       e(i2)=dm
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       ENDIF
*(6) for pi(0)+n-->D(0)+P(0), D(-)+p(+) or D(+)+p(-)
       if((iabs(lb(i1)).eq.2.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.2))then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=7
       e(i1)=dm
       lb(i2)=4
       e(i2)=ap1
       go to 40
       Endif
       IF(X2.LE.0.67.AND.X2.GT.0.33)THEN       
       lb(i1)=6
       e(i1)=dm
       lb(i2)=5
       e(i2)=ap1
       go to 40
       endif
       IF(X2.GT.0.67)THEN
       LB(I1)=8
       E(I1)=DM
       LB(I2)=3
       E(I2)=AP1
       GO TO 40
       ENDIF
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=7
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       ENDIF
       IF(X2.LE.0.67.AND.X2.GT.0.33)THEN       
       lb(i2)=6
       e(i2)=dm
       lb(i1)=5
       e(i1)=ap1
       go to 40
       endif
       IF(X2.GT.0.67)THEN
       LB(I2)=8
       E(I2)=DM
       LB(I1)=3
       E(I1)=AP1
       GO TO 40
       ENDIF
              endif
       endif
                     ENDIF
       if(iblock.eq.78)then
       call Rmasdd(srt,1.232,0.77,1.08,
     &  0.28,ISEED,4,dm,ameson)
       arho=AMESON
* pion+baryon-->Rho+delta
*(1) for pi(+)+p-->D(+)+rho(+) OR D(++)+rho(0)
       if( ((lb(i1).eq.1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.1))
     &        .OR. ((lb(i1).eq.-1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=8
       e(i1)=dm
       lb(i2)=27
       e(i2)=arho
       go to 40
       ELSE
       lb(i1)=9
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=8
       e(i2)=dm
       lb(i1)=27
       e(i1)=arho
       go to 40
       ELSE
       lb(i2)=9
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       endif
              endif
       endif
*(2) for pi(-)+p-->D(+)+rho(-) OR D(0)+rho(0) or D(-)+rho(+)
       if( ((lb(i1).eq.1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.1))
     &        .OR. ((lb(i1).eq.-1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=6
       e(i1)=dm
       lb(i2)=27
       e(i2)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=7
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=8
       e(i1)=dm
       lb(i2)=25
       e(i2)=arho
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=6
       e(i2)=dm
       lb(i1)=27
       e(i1)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=7
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=8
       e(i2)=dm
       lb(i1)=25
       e(i1)=arho
       go to 40
       endif
              endif
       endif
*(3) for pi(+)+n-->D(+)+rho(0) OR D(++)+rho(-) or D(0)+rho(+)
       if( ((lb(i1).eq.2.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.2))
     &       .OR.((lb(i1).eq.-2.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=8
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=7
       e(i1)=dm
       lb(i2)=27
       e(i2)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=9
       e(i1)=dm
       lb(i2)=25
       e(i2)=arho
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=8
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=7
       e(i2)=dm
       lb(i1)=27
       e(i1)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=9
       e(i2)=dm
       lb(i1)=25
       e(i1)=arho
       go to 40
       endif
              endif
       endif
*(4) for pi(0)+p-->D(+)+rho(0) OR D(++)+rho(-) or D(0)+rho(+)
       if((iabs(lb(i1)).eq.1.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.1))then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=7
       e(i1)=dm
       lb(i2)=27
       e(i2)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=8
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=9
       e(i1)=dm
       lb(i2)=25
       e(i2)=arho
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=7
       e(i2)=dm
       lb(i1)=27
       e(i1)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=8
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=9
       e(i2)=dm
       lb(i1)=25
       e(i1)=arho
       go to 40
       endif
              endif
       endif 
*(5) for pi(-)+n-->D(-)+rho(0) OR D(0)+rho(-)
       if( ((lb(i1).eq.2.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.2))
     &        .OR. ((lb(i1).eq.-2.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=6
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       ELSE
       lb(i1)=7
       e(i1)=dm
       lb(i2)=25
       e(i2)=arho
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=6
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       ELSE
       lb(i2)=7
       e(i2)=dm
       lb(i1)=25
       e(i1)=arho
       go to 40
       endif
              endif
       ENDIF
*(6) for pi(0)+n-->D(0)+rho(0), D(-)+rho(+) and D(+)+rho(-)
       if((iabs(lb(i1)).eq.2.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.2))then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=7
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       endif
       if(x2.gt.0.33.and.x2.le.0.67)then       
       lb(i1)=6
       e(i1)=dm
       lb(i2)=27
       e(i2)=arho
       go to 40
       endif
       if(x2.gt.0.67)then
       lb(i1)=8
       e(i1)=dm
       lb(i2)=25
       e(i2)=arho
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=7
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       endif
       if(x2.le.0.67.and.x2.gt.0.33)then       
       lb(i2)=6
       e(i2)=dm
       lb(i1)=27
       e(i1)=arho
       go to 40
       endif
       if(x2.gt.0.67)then
       lb(i2)=8
       e(i2)=dm
       lb(i1)=25
       e(i1)=arho
       endif
              endif
       endif
                     Endif
       if(iblock.eq.79)then
       aomega=0.782
* GENERATE THE DELTA MASS
       dmax=srt-0.782-0.02
       dm=rmass(dmax,iseed)
* pion+baryon-->omega+delta
*(1) for pi(+)+p-->D(++)+omega(0)
       if( ((lb(i1).eq.1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.1))
     &  .OR.((lb(i1).eq.-1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       lb(i1)=9
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              else
        ii = i2
       lb(i2)=9
       e(i2)=dm
       lb(i1)=28
       e(i1)=aomega
       go to 40
              endif
       endif
*(2) for pi(-)+p-->D(0)+omega(0) 
       if( ((lb(i1).eq.1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.1))
     &        .OR. ((lb(i1).eq.-1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       lb(i1)=7
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              else
        ii = i2
       lb(i2)=7
       e(i2)=dm
       lb(i1)=28
       e(i1)=aomega
       go to 40
              endif
       endif
*(3) for pi(+)+n-->D(+)+omega(0) 
       if( ((lb(i1).eq.2.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.2))
     &       .OR. ((lb(i1).eq.-2.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       lb(i1)=8
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              else
        ii = i2
       lb(i2)=8
       e(i2)=dm
       lb(i1)=28
       e(i1)=aomega
       go to 40
              endif
       endif
*(4) for pi(0)+p-->D(+)+omega(0) 
       if((iabs(lb(i1)).eq.1.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.1))then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       lb(i1)=8
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              else
        ii = i2
       lb(i2)=8
       e(i2)=dm
       lb(i1)=28
       e(i1)=aomega
       go to 40
              endif
       endif 
*(5) for pi(-)+n-->D(-)+omega(0) 
       if( ((lb(i1).eq.2.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.2))
     &        .OR. ((lb(i1).eq.-2.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       lb(i1)=6
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              ELSE
        ii = i2
       lb(i2)=6
       e(i2)=dm
       lb(i1)=28
       e(i1)=aomega
              endif
       ENDIF
*(6) for pi(0)+n-->D(0)+omega(0) 
       if((iabs(lb(i1)).eq.2.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.2))then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       lb(i1)=7
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              else
        ii = i2
       lb(i2)=7
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
              endif
       endif
                     Endif
40       em1=e(i1)
       em2=e(i2)
       if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
         lb(ii) = -lb(ii)
           jj = i2
          if(ii .eq. i2)jj = i1
         if(iblock .eq. 77)then
          if(lb(jj).eq.3)then
           lb(jj) = 5
          elseif(lb(jj).eq.5)then
           lb(jj) = 3
          endif
         elseif(iblock .eq. 78)then
          if(lb(jj).eq.25)then
           lb(jj) = 27
          elseif(lb(jj).eq.27)then
           lb(jj) = 25
          endif
         endif
       endif
           endif
*-----------------------------------------------------------------------
* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
* ENERGY CONSERVATION
50          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=0.00000001
          PR=SQRT(PR2)/(2.*SRT)
* here we use the same transverse momentum distribution as for
* pp collisions, it might be necessary to use a different distribution
clin-10/25/02 get rid of argument usage mismatch in PTR():
          xptr=0.33*pr
c         cc1=ptr(0.33*pr,iseed)
         cc1=ptr(xptr,iseed)
clin-10/25/02-end
clin-9/2012: check argument in sqrt():
         scheck=pr**2-cc1**2
         if(scheck.lt.0) then
            write(99,*) 'scheck36: ', scheck
            scheck=0.
         endif
         c1=sqrt(scheck)/pr
c         c1=sqrt(pr**2-cc1**2)/pr
*          C1   = 1.0 - 2.0 * RANART(NSEED)
          T1   = 2.0 * PI * RANART(NSEED)
      S1   = SQRT( 1.0 - C1**2 )
      CT1  = COS(T1)
      ST1  = SIN(T1)
* THE MOMENTUM IN THE CMS IN THE FINAL STATE
      PZ   = PR * C1
      PX   = PR * S1*CT1 
      PY   = PR * S1*ST1
* ROTATE IT 
       CALL ROTATE(PX0,PY0,PZ0,PX,PY,PZ) 
      RETURN
      END
**********************************
*                                                                      *
*                                                                      *
