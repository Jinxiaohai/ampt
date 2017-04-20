        SUBROUTINE RHORES(I1,I2)
        PARAMETER (MAXSTR=150001,MAXR=1,
     1  AMN=0.939457,AMP=0.93828,
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
* 1. DETERMINE THE MOMENTUM COMPONENT OF THE RHO IN THE CMS OF NN FRAME
*    WE LET I1 TO BE THE RHO AND ABSORB I2
clin-9/2012: improve precision for argument in sqrt():
c        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2)
c        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2)
        E10=dSQRT(dble(E(I1))**2+dble(P(1,I1))**2
     1     +dble(P(2,I1))**2+dble(P(3,I1))**2)
        E20=dSQRT(dble(E(I2))**2+dble(P(1,I2))**2
     1       +dble(P(2,I2))**2+dble(P(3,I2))**2)
        p1=dble(P(1,I1))+dble(P(1,I2))
        p2=dble(P(2,I1))+dble(P(2,I2))
        p3=dble(P(3,I1))+dble(P(3,I2))
        P(1,I1)=P(1,I1)+P(1,I2)
        P(2,I1)=P(2,I1)+P(2,I2)
        P(3,I1)=P(3,I1)+P(3,I2)
* 2. DETERMINE THE MASS OF THE RHO BY USING THE REACTION KINEMATICS
clin-9/2012: check argument in sqrt():
        scheck=(E10+E20)**2-p1**2-p2**2-p3**2
        if(scheck.lt.0) then
           write(99,*) 'scheck18: ', scheck
           scheck=0.d0
        endif
        DM=SQRT(sngl(scheck))
c        DM=SQRT((E10+E20)**2-P(1,I1)**2-P(2,I1)**2-P(3,I1)**2)
        E(I1)=DM
       E(I2)=0
        RETURN
        END
*---------------------------------------------------------------------------
* PURPOSE : CALCULATE THE PION+NUCLEON CROSS SECTION ACCORDING TO THE
*           BREIT-WIGNER FORMULA/(p*)**2
* VARIABLE : LA = 1 FOR DELTA RESONANCE
*            LA = 0 FOR N*(1440) RESONANCE
*            LA = 2 FRO N*(1535) RESONANCE
* DATE    : JAN.29,1990
