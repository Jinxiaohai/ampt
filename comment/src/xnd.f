      SUBROUTINE XND(px,py,pz,srt,I1,I2,xinel,
     &               sigk,xsk1,xsk2,xsk3,xsk4,xsk5)
*     PURPOSE:                                                         *
*             calculate NUCLEON-BARYON RESONANCE inelatic Xsection     *
*     NOTE   :                                                         *
*     QUANTITIES:                                                 *
*                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    *
*                      N12,                                            *
*                      M12=1 FOR p+n-->delta(+)+ n                     *
*                          2     p+n-->delta(0)+ p                     *
*                          3     p+p-->delta(++)+n                     *
*                          4     p+p-->delta(+)+p                      *
*                          5     n+n-->delta(0)+n                      *
*                          6     n+n-->delta(-)+p                      *
*                          7     n+p-->N*(0)(1440)+p                   *
*                          8     n+p-->N*(+)(1440)+n                   *
*                        9     p+p-->N*(+)(1535)+p                     *
*                        10    n+n-->N*(0)(1535)+n                     *
*                         11    n+p-->N*(+)(1535)+n                     *
*                        12    n+p-->N*(0)(1535)+p
*                        13    D(++)+D(-)-->N*(+)(1440)+n
*                         14    D(++)+D(-)-->N*(0)(1440)+p
*                        15    D(+)+D(0)--->N*(+)(1440)+n
*                        16    D(+)+D(0)--->N*(0)(1440)+p
*                        17    D(++)+D(0)-->N*(+)(1535)+p
*                        18    D(++)+D(-)-->N*(0)(1535)+p
*                        19    D(++)+D(-)-->N*(+)(1535)+n
*                        20    D(+)+D(+)-->N*(+)(1535)+p
*                        21    D(+)+D(0)-->N*(+)(1535)+n
*                        22    D(+)+D(0)-->N*(0)(1535)+p
*                        23    D(+)+D(-)-->N*(0)(1535)+n
*                        24    D(0)+D(0)-->N*(0)(1535)+n
*                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p
*                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n
*                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n
*                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p
*                        29    N*(+)(14)+D+-->N*(+)(15)+p
*                        30    N*(+)(14)+D0-->N*(+)(15)+n
*                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n
*                        32    N*(0)(14)+D++--->N*(+)(15)+p
*                        33    N*(0)(14)+D+--->N*(+)(15)+n
*                        34    N*(0)(14)+D+--->N*(0)(15)+p
*                        35    N*(0)(14)+D0-->N*(0)(15)+n
*                        36    N*(+)(14)+D0--->N*(0)(15)+p
*                            and more
***********************************
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AKA=0.498,APHI=1.020,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
        COMMON /AA/ R(3,MAXSTR)
cc      SAVE /AA/
        COMMON /BB/ P(3,MAXSTR)
cc      SAVE /BB/
        COMMON /CC/ E(MAXSTR)
cc      SAVE /CC/
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
        common /ff/f(-mx:mx,-my:my,-mz:mz,-mpx:mpx,-mpy:mpy,-mpz:mpzp)
cc      SAVE /ff/
        common /gg/ dx,dy,dz,dpx,dpy,dpz
cc      SAVE /gg/
        COMMON /INPUT/ NSTAR,NDIRCT,DIR
cc      SAVE /INPUT/
        COMMON /NN/NNN
cc      SAVE /NN/
        COMMON /BG/BETAX,BETAY,BETAZ,GAMMA
cc      SAVE /BG/
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
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      SAVE   
*-----------------------------------------------------------------------
       xinel=0.
       sigk=0
       xsk1=0
       xsk2=0
       xsk3=0
       xsk4=0
       xsk5=0
        EM1=E(I1)
        EM2=E(I2)
      PR  = SQRT( PX**2 + PY**2 + PZ**2 )
*     CAN HAPPEN ANY MORE ==> RETURN (2.04 = 2*AVMASS + PI-MASS+0.02)
        IF (SRT .LT. 2.04) RETURN
* Resonance absorption or Delta + N-->N*(1440), N*(1535)
* COM: TEST FOR DELTA OR N* ABSORPTION
*      IN THE PROCESS DELTA+N-->NN, N*+N-->NN
        PRF=SQRT(0.25*SRT**2-AVMASS**2)
        IF(EM1.GT.1.)THEN
        DELTAM=EM1
        ELSE
        DELTAM=EM2
        ENDIF
        RENOM=DELTAM*PRF**2/DENOM(SRT,1.)/PR
        RENOMN=DELTAM*PRF**2/DENOM(SRT,2.)/PR
        RENOM1=DELTAM*PRF**2/DENOM(SRT,-1.)/PR
* avoid the inelastic collisions between n+delta- -->N+N 
*       and p+delta++ -->N+N due to charge conservation,
*       but they can scatter to produce kaons 
       if((iabs(lb(i1)).eq.2).and.(iabs(lb(i2)).eq.6)) renom=0.
       if((iabs(lb(i2)).eq.2).and.(iabs(lb(i1)).eq.6)) renom=0.
       if((iabs(lb(i1)).eq.1).and.(iabs(lb(i2)).eq.9)) renom=0.
       if((iabs(lb(i2)).eq.1).and.(iabs(lb(i1)).eq.9)) renom=0.
       Call M1535(iabs(lb(i1)),iabs(lb(i2)),srt,x1535)
        X1440=(3./4.)*SIGMA(SRT,2,0,1)
* CROSS SECTION FOR KAON PRODUCTION from the four channels
* for NLK channel
       akp=0.498
       ak0=0.498
       ana=0.94
       ada=1.232
       al=1.1157
       as=1.1197
       xsk1=0
       xsk2=0
       xsk3=0
       xsk4=0
c      !! phi production
       xsk5=0
       t1nlk=ana+al+akp
       if(srt.le.t1nlk)go to 222
       XSK1=1.5*PPLPK(SRT)
* for DLK channel
       t1dlk=ada+al+akp
       t2dlk=ada+al-akp
       if(srt.le.t1dlk)go to 222
       es=srt
       pmdlk2=(es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
       pmdlk=sqrt(pmdlk2)
       XSK3=1.5*PPLPK(srt)
* for NSK channel
       t1nsk=ana+as+akp
       t2nsk=ana+as-akp
       if(srt.le.t1nsk)go to 222
       pmnsk2=(es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
       pmnsk=sqrt(pmnsk2)
       XSK2=1.5*(PPK1(srt)+PPK0(srt))
* for DSK channel
       t1DSk=aDa+aS+akp
       t2DSk=aDa+aS-akp
       if(srt.le.t1dsk)go to 222
       pmDSk2=(es**2-t1DSk**2)*(es**2-t2DSk**2)/(4.*es**2)
       pmDSk=sqrt(pmDSk2)
       XSK4=1.5*(PPK1(srt)+PPK0(srt))
csp11/21/01
c phi production
       if(srt.le.(2.*amn+aphi))go to 222
c  !! mb put the correct form
         xsk5 = 0.0001
csp11/21/01 end
* THE TOTAL KAON+ PRODUCTION CROSS SECTION IS THEN
222       SIGK=XSK1+XSK2+XSK3+XSK4
cbz3/7/99 neutralk
        XSK1 = 2.0 * XSK1
        XSK2 = 2.0 * XSK2
        XSK3 = 2.0 * XSK3
        XSK4 = 2.0 * XSK4
        SIGK = 2.0 * SIGK + xsk5
cbz3/7/99 neutralk end
* avoid the inelastic collisions between n+delta- -->N+N 
*       and p+delta++ -->N+N due to charge conservation,
*       but they can scatter to produce kaons 
       if(((iabs(lb(i1)).eq.2).and.(iabs(lb(i2)).eq.6)).OR. 
     &         ((iabs(lb(i2)).eq.2).and.(iabs(lb(i1)).eq.6)).OR.
     &         ((iabs(lb(i1)).eq.1).and.(iabs(lb(i2)).eq.9)).OR.
     &         ((iabs(lb(i2)).eq.1).and.(iabs(lb(i1)).eq.9)))THEN
       xinel=sigk
       return
       ENDIF
* WE DETERMINE THE REACTION CHANNELS IN THE FOLLOWING
* FOR n+delta(++)-->p+p or n+delta(++)-->n+N*(+)(1440),n+N*(+)(1535)
* REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN, 
        IF(LB(I1)*LB(I2).EQ.18.AND.
     &    (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))then
        SIGND=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        xinel=SIGDN+X1440+X1535+SIGK
       RETURN
       endif
* FOR p+delta(-)-->n+n or p+delta(-)-->n+N*(0)(1440),n+N*(0)(1535)
* REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN, 
        IF(LB(I1)*LB(I2).EQ.6.AND.
     &    (iabs(LB(I1)).EQ.1.OR.iabs(LB(I2)).EQ.1))THEN
        SIGND=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        xinel=SIGDN+X1440+X1535+SIGK
       RETURN
       endif
* FOR p+delta(+)-->p+p, N*(+)(144)+p, N*(+)(1535)+p
cbz11/25/98
        IF(LB(I1)*LB(I2).EQ.8.AND.
     &    (iabs(LB(I1)).EQ.1.OR.iabs(LB(I2)).EQ.1))THEN
        SIGND=1.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        xinel=SIGDN+x1440+x1535+SIGK
       RETURN
       endif
* FOR n+delta(0)-->n+n, N*(0)(144)+n, N*(0)(1535)+n
        IF(LB(I1)*LB(I2).EQ.14.AND.
     &   (iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2))THEN
        SIGND=1.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        xinel=SIGDN+x1440+x1535+SIGK
       RETURN
       endif
* FOR n+delta(+)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p,
*                       N*(+)(1535)+n,N*(0)(1535)+p
        IF(LB(I1)*LB(I2).EQ.16.AND.
     &     (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))THEN
        SIGND=0.5*SIGMA(SRT,1,1,1)+0.25*SIGMA(SRT,1,1,0)
        SIGDN=0.5*SIGND*RENOM
        xinel=SIGDN+2.*x1440+2.*x1535+SIGK
       RETURN
       endif
* FOR p+delta(0)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p,
*                       N*(+)(1535)+n,N*(0)(1535)+p
        IF(LB(I1)*LB(I2).EQ.7)THEN
        SIGND=0.5*SIGMA(SRT,1,1,1)+0.25*SIGMA(SRT,1,1,0)
        SIGDN=0.5*SIGND*RENOM
        xinel=SIGDN+2.*x1440+2.*x1535+SIGK
       RETURN
       endif
* FOR p+N*(0)(14)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p
* OR  P+N*(0)(14)-->D(+)+N, D(0)+P, 
        IF(LB(I1)*LB(I2).EQ.10.AND.
     &   (iabs(LB(I1)).EQ.1.OR.iabs(LB(I2)).EQ.1))then
        SIGND=(3./4.)*SIGMA(SRT,2,0,1)
        SIGDN=SIGND*RENOMN
        xinel=SIGDN+X1535+SIGK
       RETURN
       endif
* FOR n+N*(+)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p
        IF(LB(I1)*LB(I2).EQ.22.AND.
     &   (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))then
        SIGND=(3./4.)*SIGMA(SRT,2,0,1)
        SIGDN=SIGND*RENOMN
        xinel=SIGDN+X1535+SIGK
       RETURN
       endif
* FOR N*(1535)+N-->N+N COLLISIONS
        IF((iabs(LB(I1)).EQ.12).OR.(iabs(LB(I1)).EQ.13).OR.
     1  (iabs(LB(I2)).EQ.12).OR.(iabs(LB(I2)).EQ.13))THEN
        SIGND=X1535
        SIGDN=SIGND*RENOM1
        xinel=SIGDN+SIGK
       RETURN
       endif
        RETURN
       end
**********************************
*                                                                      *
*                                                                      *
