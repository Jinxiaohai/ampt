        SUBROUTINE DRESON(I1,I2)
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
* 1. DETERMINE THE MOMENTUM COMPONENT OF DELTA/N* IN THE LAB. FRAME
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
        IF(iabs(LB(I2)) .EQ. 1 .OR. iabs(LB(I2)) .EQ. 2 .OR.
     &     (iabs(LB(I2)) .GE. 6 .AND. iabs(LB(I2)) .LE. 17)) THEN
        E(I1)=0.
        I=I2
        ELSE
        E(I2)=0.
        I=I1
        ENDIF
        P(1,I)=P(1,I1)+P(1,I2)
        P(2,I)=P(2,I1)+P(2,I2)
        P(3,I)=P(3,I1)+P(3,I2)
* 2. DETERMINE THE MASS OF DELTA/N* BY USING THE REACTION KINEMATICS
clin-9/2012: check argument in sqrt():
        scheck=(E10+E20)**2-p1**2-p2**2-p3**2
        if(scheck.lt.0) then
           write(99,*) 'scheck17: ', scheck
           write(99,*) 'scheck17', scheck,E10,E20,P(1,I),P(2,I),P(3,I)
           write(99,*) 'scheck17-1',E(I1),P(1,I1),P(2,I1),P(3,I1)
           write(99,*) 'scheck17-2',E(I2),P(1,I2),P(2,I2),P(3,I2)
           scheck=0.d0
        endif
        DM=SQRT(sngl(scheck))
c        DM=SQRT((E10+E20)**2-P(1,I)**2-P(2,I)**2-P(3,I)**2)
        E(I)=DM
        RETURN
        END
*---------------------------------------------------------------------------
* PURPOSE : CALCULATE THE MASS AND MOMENTUM OF RHO RESONANCE 
*           AFTER PION + PION COLLISION
* DATE    : NOV. 30,1994
