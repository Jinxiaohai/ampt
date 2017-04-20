      SUBROUTINE CRDD(IRUN,PX,PY,PZ,SRT,I1,I2,IBLOCK,
     1NTAG,SIGNN,SIG,NT,ipert1)
c     1NTAG,SIGNN,SIG)
*     PURPOSE:                                                         *
*             DEALING WITH BARYON RESONANCE-BARYON RESONANCE COLLISIONS*
*     NOTE   :                                                         *
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
        parameter (xmd=1.8756,npdmax=10000)
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
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1 px1n,py1n,pz1n,dp1n
cc      SAVE /leadng/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      common /dpi/em2,lb2
      common /para8/ idpert,npertd,idxsec
      dimension ppd(3,npdmax),lbpd(npdmax)
      SAVE   
*-----------------------------------------------------------------------
       n12=0
       m12=0
        IBLOCK=0
        NTAG=0
        EM1=E(I1)
        EM2=E(I2)
      PR  = SQRT( PX**2 + PY**2 + PZ**2 )
      C2  = PZ / PR
      IF(PX .EQ. 0.0 .AND. PY .EQ. 0.0) THEN
        T2 = 0.0
      ELSE
        T2=ATAN2(PY,PX)
      END IF
      X1  = RANART(NSEED)
      ianti=0
      if(lb(i1).lt.0 .and. lb(i2).lt.0)ianti=1
clin-6/2008 Production of perturbative deuterons for idpert=1:
      call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
      if(idpert.eq.1.and.ipert1.eq.1) then
         IF (SRT .LT. 2.012) RETURN
         if((iabs(lb(i1)).ge.6.and.iabs(lb(i1)).le.13)
     1        .and.(iabs(lb(i2)).ge.6.and.iabs(lb(i2)).le.13)) then
            goto 108
         else
            return
         endif
      endif
*-----------------------------------------------------------------------
*COM: TEST FOR ELASTIC SCATTERING (EITHER N-N OR DELTA-DELTA 0R
*      N-DELTA OR N*-N* or N*-Delta)
      IF (X1 .LE. SIGNN/SIG) THEN
*COM:  PARAMETRISATION IS TAKEN FROM THE CUGNON-PAPER
        AS  = ( 3.65 * (SRT - 1.8766) )**6
        A   = 6.0 * AS / (1.0 + AS)
        TA  = -2.0 * PR**2
        X   = RANART(NSEED)
clin-10/24/02        T1  = DLOG( (1-X) * DEXP(dble(A)*dble(TA)) + X )  /  A
        T1  = sngl(DLOG(dble(1.-X)*DEXP(dble(A)*dble(TA))+dble(X)))/  A
        C1  = 1.0 - T1/TA
        T1  = 2.0 * PI * RANART(NSEED)
        IBLOCK=20
       GO TO 107
      ELSE
*COM: TEST FOR INELASTIC SCATTERING
*     IF THE AVAILABLE ENERGY IS LESS THAN THE PION-MASS, NOTHING
*     CAN HAPPEN ANY MORE ==> RETURN (2.15 = 2*AVMASS +2*PI-MASS)
        IF (SRT .LT. 2.15) RETURN
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
       ana=0.938
       ada=1.232
       al=1.1157
       as=1.1197
       xsk1=0
       xsk2=0
       xsk3=0
       xsk4=0
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
       signd=sigk+s2d
clin-6/2008
       IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c       if(x1.gt.(signd+signn)/sig)return
       if(x1.gt.(signd+signn+sdprod)/sig)return
c
* if kaon production
clin-6/2008
c       IF(SIGK/SIG.GE.RANART(NSEED))GO TO 306
       IF((SIGK+sdprod)/SIG.GE.RANART(NSEED))GO TO 306
c
* if reabsorption
       go to 1012
       ENDIF
       IDD=iabs(LB(I1)*LB(I2))
* channels have the same charge as pp 
        IF((IDD.EQ.63).OR.(IDD.EQ.64).OR.(IDD.EQ.48).
     1  OR.(IDD.EQ.49).OR.(IDD.EQ.11*11).OR.(IDD.EQ.10*10).
     2  OR.(IDD.EQ.88).OR.(IDD.EQ.66).
     3  OR.(IDD.EQ.90).OR.(IDD.EQ.70))THEN
        SIGND=X1535+SIGK+s2d
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN
        IF (X1.GT.(SIGNN+SIGND+sdprod)/SIG)RETURN
c
* if kaon production
       IF(SIGK/SIGND.GT.RANART(NSEED))GO TO 306
* if reabsorption
       if(s2d/(x1535+s2d).gt.RANART(NSEED))go to 1012
* if N*(1535) production
       IF(IDD.EQ.63)N12=17
       IF(IDD.EQ.64)N12=20
       IF(IDD.EQ.48)N12=23
       IF(IDD.EQ.49)N12=24
       IF(IDD.EQ.121)N12=25
       IF(IDD.EQ.100)N12=26
       IF(IDD.EQ.88)N12=29
       IF(IDD.EQ.66)N12=31
       IF(IDD.EQ.90)N12=32
       IF(IDD.EQ.70)N12=35
       GO TO 1011
        ENDIF
* IN DELTA+N*(1440) and N*(1440)+N*(1440) COLLISIONS, 
* N*(1535), kaon production and reabsorption are ALLOWED
* IN N*(1440)+N*(1440) COLLISIONS, ONLY N*(1535) IS ALLOWED
       IF((IDD.EQ.110).OR.(IDD.EQ.77).OR.(IDD.EQ.80))THEN
clin-6/2008
          IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c       IF(X1.GT.(SIGNN+X1535+SIGK+s2d)/SIG)RETURN
          IF(X1.GT.(SIGNN+X1535+SIGK+s2d+sdprod)/SIG)RETURN
c
       IF(SIGK/(X1535+SIGK+s2d).GT.RANART(NSEED))GO TO 306
       if(s2d/(x1535+s2d).gt.RANART(NSEED))go to 1012
       IF(IDD.EQ.77)N12=30
       IF((IDD.EQ.77).AND.(RANART(NSEED).LE.0.5))N12=36
       IF(IDD.EQ.80)N12=34
       IF((IDD.EQ.80).AND.(RANART(NSEED).LE.0.5))N12=35
       IF(IDD.EQ.110)N12=27
       IF((IDD.EQ.110).AND.(RANART(NSEED).LE.0.5))N12=28
       GO TO 1011
        ENDIF
       IF((IDD.EQ.54).OR.(IDD.EQ.56))THEN
* LIKE FOR N+P COLLISION, 
* IN DELTA+DELTA COLLISIONS BOTH N*(1440) AND N*(1535) CAN BE PRODUCED
        SIG2=(3./4.)*SIGMA(SRT,2,0,1)
        SIGND=2.*(SIG2+X1535)+SIGK+s2d
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF(X1.GT.(SIGNN+SIGND)/SIG)RETURN
        IF(X1.GT.(SIGNN+SIGND+sdprod)/SIG)RETURN
c
       IF(SIGK/SIGND.GT.RANART(NSEED))GO TO 306
       if(s2d/(2.*(sig2+x1535)+s2d).gt.RANART(NSEED))go to 1012
       IF(RANART(NSEED).LT.X1535/(SIG2+X1535))THEN
* N*(1535) PRODUCTION
       IF(IDD.EQ.54)N12=18
       IF((IDD.EQ.54).AND.(RANART(NSEED).LE.0.5))N12=19
       IF(IDD.EQ.56)N12=21
       IF((IDD.EQ.56).AND.(RANART(NSEED).LE.0.5))N12=22
               ELSE 
* N*(144) PRODUCTION
       IF(IDD.EQ.54)N12=13
       IF((IDD.EQ.54).AND.(RANART(NSEED).LE.0.5))N12=14
       IF(IDD.EQ.56)N12=15
       IF((IDD.EQ.56).AND.(RANART(NSEED).LE.0.5))N12=16
              ENDIF
       ENDIF
1011       CONTINUE
       iblock=5
*PARAMETRIZATION OF THE SHAPE OF THE N*(1440) AND N*(1535) 
* RESONANCE ACCORDING
*     TO kitazoe's or J.D.JACKSON'S MASS FORMULA AND BREIT WIGNER
*     FORMULA FOR N* RESORANCE
*     DETERMINE DELTA MASS VIA REJECTION METHOD.
          DMAX = SRT - AVMASS-0.005
          DMIN = 1.078
          IF((n12.ge.13).and.(n12.le.16))then
* N*(1440) production
          IF(DMAX.LT.1.44) THEN
          FM=FNS(DMAX,SRT,0.)
          ELSE
clin-10/25/02 get rid of argument usage mismatch in FNS():
             xdmass=1.44
c          FM=FNS(1.44,SRT,1.)
          FM=FNS(xdmass,SRT,1.)
clin-10/25/02-end
          ENDIF
          IF(FM.EQ.0.)FM=1.E-09
          NTRY2=0
11        DM=RANART(NSEED)*(DMAX-DMIN)+DMIN
          NTRY2=NTRY2+1
          IF((RANART(NSEED).GT.FNS(DM,SRT,1.)/FM).AND.
     1    (NTRY2.LE.10)) GO TO 11
clin-2/26/03 limit the N* mass below a certain value 
c     (here taken as its central value + 2* B-W fullwidth):
          if(dm.gt.2.14) goto 11
              GO TO 13
              ENDIF
                    IF((n12.ge.17).AND.(N12.LE.36))then
* N*(1535) production
          IF(DMAX.LT.1.535) THEN
          FM=FD5(DMAX,SRT,0.)
          ELSE
clin-10/25/02 get rid of argument usage mismatch in FNS():
             xdmass=1.535
c          FM=FD5(1.535,SRT,1.)
          FM=FD5(xdmass,SRT,1.)
clin-10/25/02-end
          ENDIF
          IF(FM.EQ.0.)FM=1.E-09
          NTRY1=0
12        DM = RANART(NSEED) * (DMAX-DMIN) + DMIN
          NTRY1=NTRY1+1
          IF((RANART(NSEED) .GT. FD5(DM,SRT,1.)/FM).AND.
     1    (NTRY1.LE.10)) GOTO 12
clin-2/26/03 limit the N* mass below a certain value 
c     (here taken as its central value + 2* B-W fullwidth):
          if(dm.gt.1.84) goto 12
             ENDIF
13       CONTINUE
*-------------------------------------------------------
* RELABLE BARYON I1 AND I2
*13 D(++)+D(-)--> N*(+)(14)+n
          IF(N12.EQ.13)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=11
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=11
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*14 D(++)+D(-)--> N*(0)(14)+P
          IF(N12.EQ.14)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=10
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=10
          E(I1)=DM
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
*15 D(+)+D(0)--> N*(+)(14)+n
          IF(N12.EQ.15)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=11
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=11
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*16 D(+)+D(0)--> N*(0)(14)+P
          IF(N12.EQ.16)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=10
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=10
          E(I1)=DM
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
*17 D(++)+D(0)--> N*(+)(14)+P
          IF(N12.EQ.17)THEN
          LB(I2)=13
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
         go to 200
          ENDIF
*18 D(++)+D(-)--> N*(0)(15)+P
          IF(N12.EQ.18)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=12
          E(I1)=DM
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
*19 D(++)+D(-)--> N*(+)(15)+N
          IF(N12.EQ.19)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=13
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=13
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*20 D(+)+D(+)--> N*(+)(15)+P
          IF(N12.EQ.20)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=13
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=13
          E(I1)=DM
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
*21 D(+)+D(0)--> N*(+)(15)+N
          IF(N12.EQ.21)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=13
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=13
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*22 D(+)+D(0)--> N*(0)(15)+P
          IF(N12.EQ.22)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=12
          E(I1)=DM
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
*23 D(+)+D(-)--> N*(0)(15)+N
          IF(N12.EQ.23)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=12
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*24 D(0)+D(0)--> N*(0)(15)+N
          IF(N12.EQ.24)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
         go to 200
          ENDIF
*25 N*(+)+N*(+)--> N*(0)(15)+P
          IF(N12.EQ.25)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
         go to 200
          ENDIF
*26 N*(0)+N*(0)--> N*(0)(15)+N
          IF(N12.EQ.26)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
         go to 200
          ENDIF
*27 N*(+)+N*(0)--> N*(+)(15)+N
          IF(N12.EQ.27)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=13
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=13
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*28 N*(+)+N*(0)--> N*(0)(15)+P
          IF(N12.EQ.28)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=12
          E(I1)=DM
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
*27 N*(+)+N*(0)--> N*(+)(15)+N
          IF(N12.EQ.27)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=13
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=13
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*29 N*(+)+D(+)--> N*(+)(15)+P
          IF(N12.EQ.29)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=13
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=13
          E(I1)=DM
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
*30 N*(+)+D(0)--> N*(+)(15)+N
          IF(N12.EQ.30)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=13
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=13
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*31 N*(+)+D(-)--> N*(0)(15)+N
          IF(N12.EQ.31)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=12
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*32 N*(0)+D(++)--> N*(+)(15)+P
          IF(N12.EQ.32)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=13
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=13
          E(I1)=DM
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
*33 N*(0)+D(+)--> N*(+)(15)+N
          IF(N12.EQ.33)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=13
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=13
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*34 N*(0)+D(+)--> N*(0)(15)+P
          IF(N12.EQ.34)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=12
          E(I1)=DM
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
*35 N*(0)+D(0)--> N*(0)(15)+N
          IF(N12.EQ.35)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=12
          E(I1)=DM
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*36 N*(+)+D(0)--> N*(0)(15)+P
          IF(N12.EQ.36)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=12
          E(I2)=DM
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=12
          E(I1)=DM
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
1012         continue
         iblock=55
         lb1=lb(i1)
         lb2=lb(i2)
         ich=iabs(lb1*lb2)
*-------------------------------------------------------
* RELABLE BARYON I1 AND I2 in the reabsorption processes
*37 D(++)+D(-)--> n+p
          IF(ich.EQ.9*6)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=1
          E(I2)=amp
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=1
          E(I1)=amp
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*38 D(+)+D(0)--> n+p
          IF(ich.EQ.8*7)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=1
          E(I2)=amp
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=1
          E(I1)=amp
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*39 D(++)+D(0)--> p+p
          IF(ich.EQ.9*7)THEN
          LB(I2)=1
          E(I2)=amp
         LB(I1)=1
         E(I1)=AMP
         go to 200
          ENDIF
*40 D(+)+D(+)--> p+p
          IF(ich.EQ.8*8)THEN
          LB(I2)=1
          E(I2)=amp
         LB(I1)=1
         E(I1)=AMP
          go to 200
          ENDIF
*41 D(+)+D(-)--> n+n
          IF(ich.EQ.8*6)THEN
          LB(I2)=2
          E(I2)=amn
         LB(I1)=2
         E(I1)=AMN
          go to 200
          ENDIF
*42 D(0)+D(0)--> n+n
          IF(ich.EQ.6*6)THEN
          LB(I2)=2
          E(I2)=amn
         LB(I1)=2
         E(I1)=AMN
         go to 200
          ENDIF
*43 N*(+)+N*(+)--> p+p
          IF(ich.EQ.11*11.or.ich.eq.13*13.or.ich.eq.11*13)THEN
          LB(I2)=1
          E(I2)=amp
         LB(I1)=1
         E(I1)=AMP
         go to 200
          ENDIF
*44 N*(0)(1440)+N*(0)--> n+n
          IF(ich.EQ.10*10.or.ich.eq.12*12.or.ich.eq.10*12)THEN
          LB(I2)=2
          E(I2)=amn
         LB(I1)=2
         E(I1)=AMN
         go to 200
          ENDIF
*45 N*(+)+N*(0)--> n+p
          IF(ich.EQ.10*11.or.ich.eq.12*13.or.ich.
     &    eq.10*13.or.ich.eq.11*12)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=1
          E(I2)=amp
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=1
          E(I1)=amp
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*46 N*(+)+D(+)--> p+p
          IF(ich.eq.11*8.or.ich.eq.13*8)THEN
          LB(I2)=1
          E(I2)=amp
         LB(I1)=1
         E(I1)=AMP
          go to 200
          ENDIF
*47 N*(+)+D(0)--> n+p
          IF(ich.EQ.11*7.or.ich.eq.13*7)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=1
          E(I2)=amp
         LB(I1)=2
         E(I1)=AMN
          ELSE
          LB(I1)=1
          E(I1)=amp
         LB(I2)=2
         E(I2)=AMN
          ENDIF
         go to 200
          ENDIF
*48 N*(+)+D(-)--> n+n
          IF(ich.EQ.11*6.or.ich.eq.13*6)THEN
          LB(I2)=2
          E(I2)=amn
         LB(I1)=2
         E(I1)=AMN
          go to 200
          ENDIF
*49 N*(0)+D(++)--> p+p
          IF(ich.EQ.10*9.or.ich.eq.12*9)THEN
          LB(I2)=1
          E(I2)=amp
         LB(I1)=1
         E(I1)=AMP
         go to 200
          ENDIF
*50 N*(0)+D(0)--> n+n
          IF(ich.EQ.10*7.or.ich.eq.12*7)THEN
          LB(I2)=2
          E(I2)=amn
         LB(I1)=2
         E(I1)=AMN
          go to 200
          ENDIF
*51 N*(0)+D(+)--> n+p
          IF(ich.EQ.10*8.or.ich.eq.12*8)THEN
          IF(RANART(NSEED).LE.0.5)THEN
          LB(I2)=2
          E(I2)=amn
         LB(I1)=1
         E(I1)=AMP
          ELSE
          LB(I1)=2
          E(I1)=amn
         LB(I2)=1
         E(I2)=AMP
          ENDIF
         go to 200
          ENDIF
         lb(i1)=1
         e(i1)=amp
         lb(i2)=2
         e(i2)=amn
* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
* ENERGY CONSERVATION
* resonance production or absorption in resonance+resonance collisions is
* assumed to have the same pt distribution as pp
200       EM1=E(I1)
          EM2=E(I2)
          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.e-09
          PR=SQRT(PR2)/(2.*SRT)
             if(srt.le.2.14)C1= 1.0 - 2.0 * RANART(NSEED)
         if(srt.gt.2.14.and.srt.le.2.4)c1=ang(srt,iseed)       
         if(srt.gt.2.4)then
clin-10/25/02 get rid of argument usage mismatch in PTR():
             xptr=0.33*pr
c         cc1=ptr(0.33*pr,iseed)
         cc1=ptr(xptr,iseed)
clin-10/25/02-end
clin-9/2012: check argument in sqrt():
         scheck=pr**2-cc1**2
         if(scheck.lt.0) then
            write(99,*) 'scheck7: ', scheck
            scheck=0.
         endif
         c1=sqrt(scheck)/pr
c         c1=sqrt(pr**2-cc1**2)/pr
         endif
          T1   = 2.0 * PI * RANART(NSEED)
       if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
         lb(i1) = -lb(i1)
         lb(i2) = -lb(i2)
       endif
         ENDIF
*COM: SET THE NEW MOMENTUM COORDINATES
clin-9/2012: check argument in sqrt():
 107     scheck=1.0 - C1**2
         if(scheck.lt.0) then
            write(99,*) 'scheck8: ', scheck
            scheck=0.
         endif
         S1=SQRT(scheck)
c107   S1   = SQRT( 1.0 - C1**2 )
clin-9/2012: check argument in sqrt():
      scheck=1.0 - C2**2
      if(scheck.lt.0) then
         write(99,*) 'scheck9: ', scheck
         scheck=0.
      endif
      S2=SQRT(scheck)
c      S2  =  SQRT( 1.0 - C2**2 )
      CT1  = COS(T1)
      ST1  = SIN(T1)
      CT2  = COS(T2)
      ST2  = SIN(T2)
      PZ   = PR * ( C1*C2 - S1*S2*CT1 )
      SS   = C2 * S1 * CT1  +  S2 * C1
      PX   = PR * ( SS*CT2 - S1*ST1*ST2 )
      PY   = PR * ( SS*ST2 + S1*ST1*CT2 )
      RETURN
* FOR THE DD-->KAON+X PROCESS, FIND MOMENTUM OF THE FINAL PARTICLES IN 
* THE NUCLEUS-NUCLEUS CMS.
306     CONTINUE
csp11/21/01 phi production
              if(XSK5/sigK.gt.RANART(NSEED))then
              pz1=p(3,i1)
              pz2=p(3,i2)
                LB(I1) = 1 + int(2 * RANART(NSEED))
                LB(I2) = 1 + int(2 * RANART(NSEED))
              nnn=nnn+1
                LPION(NNN,IRUN)=29
                EPION(NNN,IRUN)=APHI
                iblock = 222
              GO TO 208
               ENDIF
              iblock=10
                if(ianti .eq. 1)iblock=-10
              pz1=p(3,i1)
              pz2=p(3,i2)
* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
              nnn=nnn+1
                LPION(NNN,IRUN)=23
                EPION(NNN,IRUN)=Aka
              if(srt.le.2.63)then
* only lambda production is possible
* (1.1)P+P-->p+L+kaon+
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              GO TO 208
                ENDIF
       if(srt.le.2.74.and.srt.gt.2.63)then
* both Lambda and sigma production are possible
              if(XSK1/(XSK1+XSK2).gt.RANART(NSEED))then
* lambda production
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              else
* sigma production
                LB(I1) = 1 + int(2 * RANART(NSEED))
                LB(I2) = 15 + int(3 * RANART(NSEED))
              ic=2
              endif
              GO TO 208
       endif
       if(srt.le.2.77.and.srt.gt.2.74)then
* then pp-->Delta lamda kaon can happen
              if(xsk1/(xsk1+xsk2+xsk3).gt.RANART(NSEED))then
* * (1.1)P+P-->p+L+kaon+
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              go to 208
              else
              if(xsk2/(xsk2+xsk3).gt.RANART(NSEED))then
* pp-->psk
              ic=2
                LB(I1) = 1 + int(2 * RANART(NSEED))
                LB(I2) = 15 + int(3 * RANART(NSEED))
              else
* pp-->D+l+k        
              ic=3
                LB(I1) = 6 + int(4 * RANART(NSEED))
              lb(i2)=14
              endif
              GO TO 208
              endif
       endif
       if(srt.gt.2.77)then
* all four channels are possible
              if(xsk1/(xsk1+xsk2+xsk3+xsk4).gt.RANART(NSEED))then
* p lambda k production
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              go to 208
       else
          if(xsk3/(xsk2+xsk3+xsk4).gt.RANART(NSEED))then
* delta l K production
              ic=3
                LB(I1) = 6 + int(4 * RANART(NSEED))
              lb(i2)=14
              go to 208
          else
              if(xsk2/(xsk2+xsk4).gt.RANART(NSEED))then
* n sigma k production
                LB(I1) = 1 + int(2 * RANART(NSEED))
                LB(I2) = 15 + int(3 * RANART(NSEED))
              ic=2
              else
* D sigma K
              ic=4
                LB(I1) = 6 + int(4 * RANART(NSEED))
                LB(I2) = 15 + int(3 * RANART(NSEED))
              endif
              go to 208
          endif
       endif
       endif
208             continue
         if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
          lb(i1) = - lb(i1)
          lb(i2) = - lb(i2)
          if(LPION(NNN,IRUN) .eq. 23)LPION(NNN,IRUN)=21
         endif
       lbi1=lb(i1)
       lbi2=lb(i2)
* KEEP ALL COORDINATES OF PARTICLE 2 FOR POSSIBLE PHASE SPACE CHANGE
           NTRY1=0
129        CALL BBKAON(ic,SRT,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 129
c       if(icou1.lt.0)return
* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
* NUCLEUS CMS. FRAME 
* (1) for the necleon/delta
*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
              E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
* (2) for the lambda/sigma
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
* GET the kaon'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
                EPCM=SQRT(aka**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
clin-5/2008:
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
clin-5/2008:
c2007        X01 = 1.0 - 2.0 * RANART(NSEED)
c            Y01 = 1.0 - 2.0 * RANART(NSEED)
c            Z01 = 1.0 - 2.0 * RANART(NSEED)
c        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2007
c                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
c                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
c                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
                    RPION(1,NNN,IRUN)=R(1,I1)
                    RPION(2,NNN,IRUN)=R(2,I1)
                    RPION(3,NNN,IRUN)=R(3,I1)
c
* assign the nucleon/delta and lambda/sigma to i1 or i2 to keep the 
* leadng particle behaviour
C              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lbi1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lbi2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
        LB1=LB(I1)
        LB2=LB(I2)
        AM1=EM1
       am2=em2
        E1= SQRT( EM1**2 + PX1**2 + PY1**2 + PZ1**2 )
       RETURN
clin-6/2008 D+D->Deuteron+pi:
*     FIND MOMENTUM OF THE FINAL PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
 108   CONTINUE
           if(idpert.eq.1.and.ipert1.eq.1.and.npertd.ge.1) then
c     For idpert=1: we produce npertd pert deuterons:
              ndloop=npertd
           elseif(idpert.eq.2.and.npertd.ge.1) then
c     For idpert=2: we first save information for npertd pert deuterons;
c     at the last ndloop we create the regular deuteron+pi 
c     and those pert deuterons:
              ndloop=npertd+1
           else
c     Just create the regular deuteron+pi:
              ndloop=1
           endif
c
           dprob1=sdprod/sig/float(npertd)
           do idloop=1,ndloop
              CALL bbdangle(pxd,pyd,pzd,nt,ipert1,ianti,idloop,pfinal,
     1 dprob1,lbm)
              CALL ROTATE(PX,PY,PZ,PXd,PYd,PZd)
*     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE 
*     FROM THE NN CMS FRAME INTO THE GLOBAL CMS FRAME:
*     For the Deuteron:
              xmass=xmd
              E1dCM=SQRT(xmass**2+PXd**2+PYd**2+PZd**2)
              P1dBETA=PXd*BETAX+PYd*BETAY+PZd*BETAZ
              TRANSF=GAMMA*(GAMMA*P1dBETA/(GAMMA+1.)+E1dCM)
              pxi1=BETAX*TRANSF+PXd
              pyi1=BETAY*TRANSF+PYd
              pzi1=BETAZ*TRANSF+PZd
              if(ianti.eq.0)then
                 lbd=42
              else
                 lbd=-42
              endif
              if(idpert.eq.1.and.ipert1.eq.1.and.npertd.ge.1) then
cccc  Perturbative production for idpert=1:
                 nnn=nnn+1
                 PPION(1,NNN,IRUN)=pxi1
                 PPION(2,NNN,IRUN)=pyi1
                 PPION(3,NNN,IRUN)=pzi1
                 EPION(NNN,IRUN)=xmd
                 LPION(NNN,IRUN)=lbd
                 RPION(1,NNN,IRUN)=R(1,I1)
                 RPION(2,NNN,IRUN)=R(2,I1)
                 RPION(3,NNN,IRUN)=R(3,I1)
clin-6/2008 assign the perturbative probability:
                 dppion(NNN,IRUN)=sdprod/sig/float(npertd)
              elseif(idpert.eq.2.and.idloop.le.npertd) then
clin-6/2008 For idpert=2, we produce NPERTD perturbative (anti)deuterons 
c     only when a regular (anti)deuteron+pi is produced in NN collisions.
c     First save the info for the perturbative deuterons:
                 ppd(1,idloop)=pxi1
                 ppd(2,idloop)=pyi1
                 ppd(3,idloop)=pzi1
                 lbpd(idloop)=lbd
              else
cccc  Regular production:
c     For the regular pion: do LORENTZ-TRANSFORMATION:
                 E(i1)=xmm
                 E2piCM=SQRT(xmm**2+PXd**2+PYd**2+PZd**2)
                 P2piBETA=-PXd*BETAX-PYd*BETAY-PZd*BETAZ
                 TRANSF=GAMMA*(GAMMA*P2piBETA/(GAMMA+1.)+E2piCM)
                 pxi2=BETAX*TRANSF-PXd
                 pyi2=BETAY*TRANSF-PYd
                 pzi2=BETAZ*TRANSF-PZd
                 p(1,i1)=pxi2
                 p(2,i1)=pyi2
                 p(3,i1)=pzi2
c     Remove regular pion to check the equivalence 
c     between the perturbative and regular deuteron results:
c                 E(i1)=0.
c
                 LB(I1)=lbm
                 PX1=P(1,I1)
                 PY1=P(2,I1)
                 PZ1=P(3,I1)
                 EM1=E(I1)
                 ID(I1)=2
                 ID1=ID(I1)
                 E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
                 lb1=lb(i1)
c     For the regular deuteron:
                 p(1,i2)=pxi1
                 p(2,i2)=pyi1
                 p(3,i2)=pzi1
                 lb(i2)=lbd
                 lb2=lb(i2)
                 E(i2)=xmd
                 EtI2=E(I2)
                 ID(I2)=2
c     For idpert=2: create the perturbative deuterons:
                 if(idpert.eq.2.and.idloop.eq.ndloop) then
                    do ipertd=1,npertd
                       nnn=nnn+1
                       PPION(1,NNN,IRUN)=ppd(1,ipertd)
                       PPION(2,NNN,IRUN)=ppd(2,ipertd)
                       PPION(3,NNN,IRUN)=ppd(3,ipertd)
                       EPION(NNN,IRUN)=xmd
                       LPION(NNN,IRUN)=lbpd(ipertd)
                       RPION(1,NNN,IRUN)=R(1,I1)
                       RPION(2,NNN,IRUN)=R(2,I1)
                       RPION(3,NNN,IRUN)=R(3,I1)
clin-6/2008 assign the perturbative probability:
                       dppion(NNN,IRUN)=1./float(npertd)
                    enddo
                 endif
              endif
           enddo
           IBLOCK=501
           return
clin-6/2008 D+D->Deuteron+pi over
        END
**********************************
**********************************
*                                                                      *
