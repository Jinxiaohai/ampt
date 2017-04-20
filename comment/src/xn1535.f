        REAL FUNCTION XN1535(I1,I2,LA)
        PARAMETER (MAXSTR=150001,MAXR=1,
     1  AMN=0.939457,AMP=0.93828,ETAM=0.5475,
     2  AP1=0.13496,AP2=0.13957,AM0=1.232,PI=3.1415926)
clin-9/2012: improve precision for argument in sqrt():
        double precision e10,e20,scheck,p1,p2,p3
        COMMON /AA/ R(3,MAXSTR)
cc      SAVE /AA/
        COMMON /BB/ P(3,MAXSTR)
cc      SAVE /BB/
        COMMON /CC/ E(MAXSTR)
cc      SAVE /CC/
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
        COMMON   /RUN/NUM
cc      SAVE /RUN/
        COMMON   /PA/RPION(3,MAXSTR,MAXR)
cc      SAVE /PA/
        COMMON   /PB/PPION(3,MAXSTR,MAXR)
cc      SAVE /PB/
        COMMON   /PC/EPION(MAXSTR,MAXR)
cc      SAVE /PC/
        COMMON   /PD/LPION(MAXSTR,MAXR)
cc      SAVE /PD/
      SAVE   
        AVMASS=0.5*(AMN+AMP)
        AVPI=(2.*AP2+AP1)/3.
* 1. DETERMINE THE MOMENTUM COMPONENT OF N*(1535) IN THE LAB. FRAME
clin-9/2012: improve precision for argument in sqrt():
c        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2)
c        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2)
        E10=dSQRT(dble(E(I1))**2+dble(P(1,I1))**2
     1     +dble(P(2,I1))**2+dble(P(3,I1))**2)
        E20=dSQRT(dble(E(I2))**2+dble(P(1,I2))**2
     1       +dble(P(2,I2))**2+dble(P(3,I2))**2)
c        P1=P(1,I1)+P(1,I2)
c        P2=P(2,I1)+P(2,I2)
c        P3=P(3,I1)+P(3,I2)
        p1=dble(P(1,I1))+dble(P(1,I2))
        p2=dble(P(2,I1))+dble(P(2,I2))
        p3=dble(P(3,I1))+dble(P(3,I2))
* 2. DETERMINE THE MASS OF DELTA BY USING OF THE REACTION KINEMATICS
clin-9/2012: check argument in sqrt():
        scheck=(E10+E20)**2-p1**2-p2**2-p3**2
        if(scheck.lt.0) then
           write(99,*) 'scheck21: ', scheck
           scheck=0.d0
        endif
        DM=SQRT(sngl(scheck))
c        DM=SQRT((E10+E20)**2-P1**2-P2**2-P3**2)
        IF(DM.LE.1.1) THEN
        XN1535=1.E-06
        RETURN
        ENDIF
* 3. DETERMINE THE PION(ETA)+NUCLEON->N*(1535) CROSS SECTION ACCORDING TO THE
*    BREIT-WIGNER FORMULA IN UNIT OF FM**2
        GAM=W1535(DM)
       GAM0=0.15
        F1=0.25*GAM0**2/(0.25*GAM**2+(DM-1.535)**2)
        IF(LA.EQ.1)THEN
       XMAX=11.3
        ELSE
       XMAX=74.
        ENDIF
        XN1535=F1*XMAX/10.
        RETURN
        END
***************************8
*FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF
*KITAZOE'S FORMULA
