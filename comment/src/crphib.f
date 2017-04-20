        SUBROUTINE CRPHIB(PX,PY,PZ,SRT,I1,I2,
     &     XSK1, XSK2, XSK3, XSK4, XSK5, SIGP, IBLOCK)
*
*     PURPOSE:                                                         *
*             DEALING WITH PHI + N(D) --> pi+N(D), rho+N(D),  K+ + La
*     QUANTITIES:                                                      *
*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
*           SRT      - SQRT OF S                                       *
*           IBLOCK   - INFORMATION about the reaction channel          *
*                
*             iblock   - 20  elastic
*             iblock   - 221  K+ formation
*             iblock   - 223  others
**********************************
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AMRHO=0.769,AMOMGA=0.782,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974,ARHO=0.77)
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
c
       PX0=PX
       PY0=PY
       PZ0=PZ
       IBLOCK=223
c
        X1 = RANART(NSEED) * SIGP
        XSK2 = XSK1 + XSK2
        XSK3 = XSK2 + XSK3
        XSK4 = XSK3 + XSK4
        XSK5 = XSK4 + XSK5
c
c  !! elastic scatt.
        IF (X1 .LE. XSK1) THEN
           iblock=20
           GOTO 100
        ELSE IF (X1 .LE. XSK2) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 1 + int(2 * RANART(NSEED))
           E(I1) = AP1
           E(I2) = AMN
           GOTO 100
        ELSE IF (X1 .LE. XSK3) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 6 + int(4 * RANART(NSEED))
           E(I1) = AP1
           E(I2) = AM0
           GOTO 100
        ELSE IF (X1 .LE. XSK4) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = 1 + int(2 * RANART(NSEED))
           E(I1) = ARHO
           E(I2) = AMN
           GOTO 100
        ELSE IF (X1 .LE. XSK5) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = 6 + int(4 * RANART(NSEED))
           E(I1) = ARHO
           E(I2) = AM0
           GOTO 100
        ELSE 
           LB(I1) = 23
           LB(I2) = 14
           E(I1) = AKA
           E(I2) = ALA
          IBLOCK=221
         ENDIF
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
c
*****************************
* purpose: Xsection for Phi + B 
c!! in fm^2
