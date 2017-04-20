       SUBROUTINE CRPHIM(PX,PY,PZ,SRT,I1,I2,
     &  XSK1, XSK2, XSK3, XSK4, XSK5, XSK6, SIGPHI, IKKG, IKKL, IBLOCK)
*
*     QUANTITIES:                                                      *
*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
*           SRT      - SQRT OF S                                       *
*           IBLOCK   - THE INFORMATION BACK                            *
*                      20 -->  elastic
*                      223 --> phi + pi(rho,omega)
*                      224 --> phi + K -> K + pi(rho,omega)
*                      225 --> phi + K -> K* + pi(rho,omega)
*                      226 --> phi + K* -> K + pi(rho,omega)
*                      227 --> phi + K* -> K* + pi(rho,omega)
**********************************
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,ARHO=0.77,AOMEGA=0.7819,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        PARAMETER    (AKA=0.498,AKS=0.895)
        parameter   (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
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
c
       PX0=PX
       PY0=PY
       PZ0=PZ
         LB1 = LB(i1)
         LB2 = LB(i2)
        X1 = RANART(NSEED) * SIGPHI
        XSK2 = XSK1 + XSK2
        XSK3 = XSK2 + XSK3
        XSK4 = XSK3 + XSK4
        XSK5 = XSK4 + XSK5
        XSK6 = XSK5 + XSK6
        IF (X1 .LE. XSK1) THEN
c        !! elastic scatt
           IBLOCK=20
           GOTO 100
        ELSE
c
*phi + (K,K*)-bar
       if( lb1.eq.23.or.lb1.eq.21.or.iabs(lb1).eq.30 .OR.
     &     lb2.eq.23.or.lb2.eq.21.or.iabs(lb2).eq.30 )then
c
             if(lb1.eq.23.or.lb2.eq.23)then
               IKKL=1
               IBLOCK=224
               iad1 = 23
               iad2 = 30
              elseif(lb1.eq.30.or.lb2.eq.30)then
               IKKL=0
               IBLOCK=226
               iad1 = 23
               iad2 = 30
             elseif(lb1.eq.21.or.lb2.eq.21)then
               IKKL=1
               IBLOCK=124
               iad1 = 21
               iad2 = -30
c         !! -30
             else
               IKKL=0
               IBLOCK=126
               iad1 = 21
               iad2 = -30
              endif
         IF (X1 .LE. XSK2) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = iad1
           E(I1) = AP1
           E(I2) = AKA
           IKKG = 1
           GOTO 100
        ELSE IF (X1 .LE. XSK3) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = iad1
           E(I1) = ARHO
           E(I2) = AKA
           IKKG = 1
           GOTO 100
        ELSE IF (X1 .LE. XSK4) THEN
           LB(I1) = 28
           LB(I2) = iad1
           E(I1) = AOMEGA
           E(I2) = AKA
           IKKG = 1
           GOTO 100
        ELSE IF (X1 .LE. XSK5) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = iad2
           E(I1) = AP1
           E(I2) = AKS
           IKKG = 0
           IBLOCK=IBLOCK+1
           GOTO 100
        ELSE IF (X1 .LE. XSK6) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = iad2
           E(I1) = ARHO
           E(I2) = AKS
           IKKG = 0
           IBLOCK=IBLOCK+1
           GOTO 100
        ELSE 
           LB(I1) = 28
           LB(I2) = iad2
           E(I1) = AOMEGA
           E(I2) = AKS
           IKKG = 0
           IBLOCK=IBLOCK+1
           GOTO 100
         ENDIF
       else
c      !! phi destruction via (pi,rho,omega)
          IBLOCK=223
*phi + pi(rho,omega)
         IF (X1 .LE. XSK2) THEN
           LB(I1) = 23
           LB(I2) = 21
           E(I1) = AKA
           E(I2) = AKA
           IKKG = 2
           IKKL = 0
           GOTO 100
        ELSE IF (X1 .LE. XSK3) THEN
           LB(I1) = 23
c           LB(I2) = 30
           LB(I2) = -30
clin-2/10/03 currently take XSK3 to be the sum of KK*bar & KbarK*:
           if(RANART(NSEED).le.0.5) then
              LB(I1) = 21
              LB(I2) = 30
           endif
           E(I1) = AKA
           E(I2) = AKS
           IKKG = 1
           IKKL = 0
           GOTO 100
        ELSE IF (X1 .LE. XSK4) THEN
           LB(I1) = 30
c           LB(I2) = 30
           LB(I2) = -30
           E(I1) = AKS
           E(I2) = AKS
           IKKG = 0
           IKKL = 0
           GOTO 100
         ENDIF
       endif
         ENDIF
*
100    CONTINUE
       EM1=E(I1)
       EM2=E(I2)
*-----------------------------------------------------------------------
* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
* ENERGY CONSERVATION
          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.E-08
          PR=SQRT(PR2)/(2.*SRT)
* WE ASSUME AN ISOTROPIC ANGULAR DISTRIBUTION IN THE CMS 
          C1   = 1.0 - 2.0 * RANART(NSEED)
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
**********************************
cbz3/9/99 khyperon
*************************************
* purpose: Xsection for K+Y ->  piN                                       *
*          Xsection for K+Y-bar ->  piN-bar   !! sp03/29/01               *
*
