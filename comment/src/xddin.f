      SUBROUTINE XDDIN(PX,PY,PZ,SRT,I1,I2,
     &XINEL,SIGK,XSK1,XSK2,XSK3,XSK4,XSK5)
*     PURPOSE:                                                         *
*             DEALING WITH BARYON RESONANCE-BARYON RESONANCE COLLISIONS*
*     NOTE   :                                                         *
*           VALID ONLY FOR BARYON-BARYON-DISTANCES LESS THAN 1.32 FM   *
*           (1.32 = 2 * HARD-CORE-RADIUS [HRC] )                       *
*     QUANTITIES:                                                 *
*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
*           SRT      - SQRT OF S                                       *
*           NSTAR =1 INCLUDING N* RESORANCE,ELSE NOT                   *
*           NDIRCT=1 INCLUDING DIRECT PION PRODUCTION PROCESS         *
*           IBLOCK   - THE INFORMATION BACK                            *
*                      0-> COLLISION CANNOT HAPPEN                     *
*                      1-> N-N ELASTIC COLLISION                       *
*                      2-> N+N->N+DELTA,OR N+N->N+N* REACTION          *
*                      3-> N+DELTA->N+N OR N+N*->N+N REACTION          *
*                      4-> N+N->N+N+PION,DIRTCT PROCESS                *
*                     5-> DELTA(N*)+DELTA(N*)   TOTAL   COLLISIONS    *
*           N12       - IS USED TO SPECIFY BARYON-BARYON REACTION      *
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
*                        +++
*               AND MORE CHANNELS AS LISTED IN THE NOTE BOOK      
*
* NOTE ABOUT N*(1440) RESORANCE:                                       *
*     As it has been discussed in VerWest's paper,I= 1 (initial isospin)
*     channel can all be attributed to delta resorance while I= 0      *
*     channel can all be  attribured to N* resorance.Only in n+p       *
*     one can have I=0 channel so is the N*(1440) resorance            *
* REFERENCES:    J. CUGNON ET AL., NUCL. PHYS. A352, 505 (1981)        *
*                    Y. KITAZOE ET AL., PHYS. LETT. 166B, 35 (1986)    *
*                    B. VerWest el al., PHYS. PRV. C25 (1982)1979      *
*                    Gy. Wolf  et al, Nucl Phys A517 (1990) 615        *
*                    CUTOFF = 2 * AVMASS + 20 MEV                      *
*                                                                      *
*       for N*(1535) we use the parameterization by Gy. Wolf et al     *
*       Nucl phys A552 (1993) 349, added May 18, 1994                  *
**********************************
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
       XINEL=0
       SIGK=0
       XSK1=0
       XSK2=0
       XSK3=0
       XSK4=0
       XSK5=0
        EM1=E(I1)
        EM2=E(I2)
      PR  = SQRT( PX**2 + PY**2 + PZ**2 )
*     IF THERE WERE 2 N*(1535) AND THEY DIDN'T SCATT. ELAST., 
*     ALLOW THEM TO PRODUCE KAONS. NO OTHER INELASTIC CHANNELS
*     ARE KNOWN
C       if((lb(i1).ge.12).and.(lb(i2).ge.12))return
*     ALL the inelastic collisions between N*(1535) and Delta as well
*     as N*(1440) TO PRODUCE KAONS, NO OTHER CHANNELS ARE KNOWN
C       if((lb(i1).ge.12).and.(lb(i2).ge.3))return
C       if((lb(i2).ge.12).and.(lb(i1).ge.3))return
*     calculate the N*(1535) production cross section in I1+I2 collisions
       call N1535(iabs(lb(i1)),iabs(lb(i2)),srt,X1535)
c
* for Delta+Delta-->N*(1440 OR 1535)+N AND N*(1440)+N*(1440)-->N*(1535)+X 
*     AND DELTA+N*(1440)-->N*(1535)+X
* WE ASSUME THEY HAVE THE SAME CROSS SECTIONS as CORRESPONDING N+N COLLISION):
* FOR D++D0, D+D+,D+D-,D0D0,N*+N*+,N*0N*0,N*(+)D+,N*(+)D(-),N*(0)D(0)
* N*(1535) production, kaon production and reabsorption through 
* D(N*)+D(N*)-->NN are ALLOWED.
* CROSS SECTION FOR KAON PRODUCTION from the four channels are
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
        IDD=iabs(LB(I1)*LB(I2))
* The reabsorption cross section for the process
* D(N*)D(N*)-->NN is
       s2d=reab2d(i1,i2,srt)
cbz3/16/99 pion
        S2D = 0.
cbz3/16/99 pion end
*(1) N*(1535)+D(N*(1440)) reactions
*    we allow kaon production and reabsorption only
       if(((iabs(lb(i1)).ge.12).and.(iabs(lb(i2)).ge.12)).OR.
     &       ((iabs(lb(i1)).ge.12).and.(iabs(lb(i2)).ge.6)).OR.
     &       ((iabs(lb(i2)).ge.12).and.(iabs(lb(i1)).ge.6)))THEN
       XINEL=sigk+s2d
       RETURN
       ENDIF
* channels have the same charge as pp 
        IF((IDD.EQ.63).OR.(IDD.EQ.64).OR.(IDD.EQ.48).
     1  OR.(IDD.EQ.49).OR.(IDD.EQ.11*11).OR.(IDD.EQ.10*10).
     2  OR.(IDD.EQ.88).OR.(IDD.EQ.66).
     3  OR.(IDD.EQ.90).OR.(IDD.EQ.70))THEN
        XINEL=X1535+SIGK+s2d
       RETURN
        ENDIF
* IN DELTA+N*(1440) and N*(1440)+N*(1440) COLLISIONS, 
* N*(1535), kaon production and reabsorption are ALLOWED
* IN N*(1440)+N*(1440) COLLISIONS, ONLY N*(1535) IS ALLOWED
       IF((IDD.EQ.110).OR.(IDD.EQ.77).OR.(IDD.EQ.80))THEN
       XINEL=X1535+SIGK+s2d
       RETURN
       ENDIF       
       IF((IDD.EQ.54).OR.(IDD.EQ.56))THEN
* LIKE FOR N+P COLLISION, 
* IN DELTA+DELTA COLLISIONS BOTH N*(1440) AND N*(1535) CAN BE PRODUCED
        SIG2=(3./4.)*SIGMA(SRT,2,0,1)
        XINEL=2.*(SIG2+X1535)+SIGK+s2d
       RETURN
       ENDIF
       RETURN
       END
******************************************
