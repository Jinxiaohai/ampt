      SUBROUTINE CRND(IRUN,PX,PY,PZ,SRT,I1,I2,IBLOCK,
     &SIGNN,SIG,sigk,xsk1,xsk2,xsk3,xsk4,xsk5,NT,ipert1)
*     PURPOSE:                                                         *
*             DEALING WITH NUCLEON-BARYON RESONANCE COLLISIONS         *
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
*                        ++    see the note book for more listing
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
        X1  = RANART(NSEED)
        ianti=0
        if(lb(i1).lt.0 .and. lb(i2).lt.0)ianti=1
clin-6/2008 Production of perturbative deuterons for idpert=1:
      call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
      if(idpert.eq.1.and.ipert1.eq.1) then
         IF (SRT .LT. 2.012) RETURN
         if((iabs(lb(i1)).eq.1.or.iabs(lb(i1)).eq.2)
     1        .and.(iabs(lb(i2)).ge.6.and.iabs(lb(i2)).le.13)) then
            goto 108
         elseif((iabs(lb(i2)).eq.1.or.iabs(lb(i2)).eq.2)
     1           .and.(iabs(lb(i1)).ge.6.and.iabs(lb(i1)).le.13)) then
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
clin-10/24/02        T1  = ALOG( (1-X) * DEXP(dble(A)*dble(TA)) + X )  /  A
        T1  = sngl(DLOG(dble(1.-X)*DEXP(dble(A)*dble(TA))+dble(X)))/  A
        C1  = 1.0 - T1/TA
        T1  = 2.0 * PI * RANART(NSEED)
        IBLOCK=1
       GO TO 107
      ELSE
*COM: TEST FOR INELASTIC SCATTERING
*     IF THE AVAILABLE ENERGY IS LESS THAN THE PION-MASS, NOTHING
*     CAN HAPPEN ANY MORE ==> RETURN (2.04 = 2*AVMASS + PI-MASS+0.02)
        IF (SRT .LT. 2.04) RETURN
clin-6/2008 add d+meson production for n*N*(0)(1440) and p*N*(+)(1440) channels
c     (they did not have any inelastic reactions before):
        if(((iabs(LB(I1)).EQ.2.or.iabs(LB(I2)).EQ.2).AND.
     1       (LB(I1)*LB(I2)).EQ.20).or.(LB(I1)*LB(I2)).EQ.13) then
           IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        ENDIF
c
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
* avoid the inelastic collisions between n+delta- -->N+N 
*       and p+delta++ -->N+N due to charge conservation,
*       but they can scatter to produce kaons 
       if(((iabs(lb(i1)).eq.2).and.(iabs(lb(i2)).eq.6)).OR. 
     &         ((iabs(lb(i2)).eq.2).and.(iabs(lb(i1)).eq.6)).OR.
     &         ((iabs(lb(i1)).eq.1).and.(iabs(lb(i2)).eq.9)).OR.
     &         ((iabs(lb(i2)).eq.1).and.(iabs(lb(i1)).eq.9)))THEN
clin-6/2008
          IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c          IF((SIGK+SIGNN)/SIG.GE.X1)GO TO 306
          IF((SIGK+SIGNN+sdprod)/SIG.GE.X1)GO TO 306
c
       ENDIF
* WE DETERMINE THE REACTION CHANNELS IN THE FOLLOWING
* FOR n+delta(++)-->p+p or n+delta(++)-->n+N*(+)(1440),n+N*(+)(1535)
* REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN, 
        IF(LB(I1)*LB(I2).EQ.18.AND.
     &  (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))then
        SIGND=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF(X1.GT.(SIGNN+SIGDN+X1440+X1535+SIGK)/SIG)RETURN
        IF(X1.GT.(SIGNN+SIGDN+X1440+X1535+SIGK+sdprod)/SIG)RETURN
c
       IF(SIGK/(SIGK+SIGDN+X1440+X1535).GT.RANART(NSEED))GO TO 306
* REABSORPTION:
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1440+X1535))THEN
        M12=3
       GO TO 206
       ELSE
* N* PRODUCTION
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
* N*(1440)
              M12=37
              ELSE
* N*(1535)       M12=38
clin-2/26/03 why is the above commented out? leads to M12=0 but 
c     particle mass is changed after 204 (causes energy violation).
c     replace by elastic process (return):
                   return
              ENDIF
       GO TO 204
       ENDIF
        ENDIF
* FOR p+delta(-)-->n+n or p+delta(-)-->n+N*(0)(1440),n+N*(0)(1535)
* REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN, 
        IF(LB(I1)*LB(I2).EQ.6.AND.
     &   ((iabs(LB(I1)).EQ.1).OR.(iabs(LB(I2)).EQ.1)))then
        SIGND=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF (X1.GT.(SIGNN+SIGDN+X1440+X1535+SIGK)/SIG)RETURN
        IF (X1.GT.(SIGNN+SIGDN+X1440+X1535+SIGK+sdprod)/SIG)RETURN
c
       IF(SIGK/(SIGK+SIGDN+X1440+X1535).GT.RANART(NSEED))GO TO 306
* REABSORPTION:
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1440+X1535))THEN
        M12=6
       GO TO 206
       ELSE
* N* PRODUCTION
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
* N*(1440)
              M12=47
              ELSE
* N*(1535)       M12=48
clin-2/26/03 causes energy violation, replace by elastic process (return):
                   return
              ENDIF
       GO TO 204
       ENDIF
        ENDIF
* FOR p+delta(+)-->p+p, N*(+)(144)+p, N*(+)(1535)+p
        IF(LB(I1)*LB(I2).EQ.8.AND.
     &   (iabs(LB(I1)).EQ.1.OR.iabs(LB(I2)).EQ.1))THEN
        SIGND=1.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF(X1.GT.(SIGNN+SIGDN+x1440+x1535+SIGK)/SIG)RETURN
        IF(X1.GT.(SIGNN+SIGDN+x1440+x1535+SIGK+sdprod)/SIG)RETURN
c
       IF(SIGK/(SIGK+SIGDN+X1440+X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1440+X1535))THEN
        M12=4
       GO TO 206
       ELSE
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
* N*(144)
              M12=39
              ELSE
              M12=40
              ENDIF
              GO TO 204
       ENDIF
        ENDIF
* FOR n+delta(0)-->n+n, N*(0)(144)+n, N*(0)(1535)+n
        IF(LB(I1)*LB(I2).EQ.14.AND.
     &   (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))THEN
        SIGND=1.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF(X1.GT.(SIGNN+SIGDN+x1440+x1535+SIGK)/SIG)RETURN
        IF(X1.GT.(SIGNN+SIGDN+x1440+x1535+SIGK+sdprod)/SIG)RETURN
c
       IF(SIGK/(SIGK+SIGDN+X1440+X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1440+X1535))THEN
        M12=5
       GO TO 206
       ELSE
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
* N*(144)
              M12=48
              ELSE
              M12=49
              ENDIF
              GO TO 204
       ENDIF
        ENDIF
* FOR n+delta(+)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p,
*                       N*(+)(1535)+n,N*(0)(1535)+p
        IF(LB(I1)*LB(I2).EQ.16.AND.
     &   (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))THEN
        SIGND=0.5*SIGMA(SRT,1,1,1)+0.25*SIGMA(SRT,1,1,0)
        SIGDN=0.5*SIGND*RENOM
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF(X1.GT.(SIGNN+SIGDN+2.*x1440+2.*x1535+SIGK)/SIG)RETURN
        IF(X1.GT.(SIGNN+SIGDN+2.*x1440+2.*x1535+SIGK+sdprod)/SIG)RETURN
c
       IF(SIGK/(SIGK+SIGDN+2*X1440+2*X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+2.*X1440+2.*X1535))THEN
        M12=1
       GO TO 206
       ELSE
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
              M12=41
              IF(RANART(NSEED).LE.0.5)M12=43
              ELSE
              M12=42
              IF(RANART(NSEED).LE.0.5)M12=44
              ENDIF
              GO TO 204
       ENDIF
        ENDIF
* FOR p+delta(0)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p,
*                       N*(+)(1535)+n,N*(0)(1535)+p
        IF(LB(I1)*LB(I2).EQ.7)THEN
        SIGND=0.5*SIGMA(SRT,1,1,1)+0.25*SIGMA(SRT,1,1,0)
        SIGDN=0.5*SIGND*RENOM
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF(X1.GT.(SIGNN+SIGDN+2.*x1440+2.*x1535+SIGK)/SIG)RETURN
        IF(X1.GT.(SIGNN+SIGDN+2.*x1440+2.*x1535+SIGK+sdprod)/SIG)RETURN
c
       IF(SIGK/(SIGK+SIGDN+2*X1440+2*X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+2.*X1440+2.*X1535))THEN
        M12=2
       GO TO 206
       ELSE
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
              M12=50
              IF(RANART(NSEED).LE.0.5)M12=51
              ELSE
              M12=52
              IF(RANART(NSEED).LE.0.5)M12=53
              ENDIF
              GO TO 204
       ENDIF
        ENDIF
* FOR p+N*(0)(14)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p
* OR  P+N*(0)(14)-->D(+)+N, D(0)+P, 
        IF(LB(I1)*LB(I2).EQ.10.AND.
     &  (iabs(LB(I1)).EQ.1.OR.iabs(LB(I2)).EQ.1))then
        SIGND=(3./4.)*SIGMA(SRT,2,0,1)
        SIGDN=SIGND*RENOMN
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF(X1.GT.(SIGNN+SIGDN+X1535+SIGK)/SIG)RETURN
        IF(X1.GT.(SIGNN+SIGDN+X1535+SIGK+sdprod)/SIG)RETURN
c
       IF(SIGK/(SIGK+SIGDN+X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1535))THEN
        M12=7
        GO TO 206
       ELSE
       M12=54
       IF(RANART(NSEED).LE.0.5)M12=55
       ENDIF
       GO TO 204
        ENDIF
* FOR n+N*(+)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p
        IF(LB(I1)*LB(I2).EQ.22.AND.
     &   (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))then
        SIGND=(3./4.)*SIGMA(SRT,2,0,1)
        SIGDN=SIGND*RENOMN
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF(X1.GT.(SIGNN+SIGDN+X1535+SIGK)/SIG)RETURN
        IF(X1.GT.(SIGNN+SIGDN+X1535+SIGK+sdprod)/SIG)RETURN
c
       IF(SIGK/(SIGK+SIGDN+X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1535))THEN
        M12=8
        GO TO 206
       ELSE
       M12=56
       IF(RANART(NSEED).LE.0.5)M12=57
       ENDIF
       GO TO 204
        ENDIF
* FOR N*(1535)+N-->N+N COLLISIONS
        IF((iabs(LB(I1)).EQ.12).OR.(iabs(LB(I1)).EQ.13).OR.
     1  (iabs(LB(I2)).EQ.12).OR.(iabs(LB(I2)).EQ.13))THEN
        SIGND=X1535
        SIGDN=SIGND*RENOM1
clin-6/2008
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
c        IF(X1.GT.(SIGNN+SIGDN+SIGK)/SIG)RETURN
        IF(X1.GT.(SIGNN+SIGDN+SIGK+sdprod)/SIG)RETURN
c
       IF(SIGK/(SIGK+SIGDN).GT.RANART(NSEED))GO TO 306
        IF(LB(I1)*LB(I2).EQ.24)M12=10
        IF(LB(I1)*LB(I2).EQ.12)M12=12
        IF(LB(I1)*LB(I2).EQ.26)M12=11
       IF(LB(I1)*LB(I2).EQ.13)M12=9
       GO TO 206
        ENDIF
204       CONTINUE
* (1) GENERATE THE MASS FOR THE N*(1440) AND N*(1535)
* (2) CALCULATE THE FINAL MOMENTUM OF THE n+N* SYSTEM
* (3) RELABLE THE FINAL STATE PARTICLES
*PARAMETRIZATION OF THE SHAPE OF THE N* RESONANCE ACCORDING
*     TO kitazoe's or J.D.JACKSON'S MASS FORMULA AND BREIT WIGNER
*     FORMULA FOR N* RESORANCE
*     DETERMINE DELTA MASS VIA REJECTION METHOD.
          DMAX = SRT - AVMASS-0.005
          DMIN = 1.078
          IF((M12.eq.37).or.(M12.eq.39).or.
     1    (M12.eQ.41).OR.(M12.eQ.43).OR.(M12.EQ.46).
     2     OR.(M12.EQ.48).OR.(M12.EQ.50).OR.(M12.EQ.51))then
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
              ELSE
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
* (2) DETERMINE THE FINAL MOMENTUM
       PRF=0.
       PF2=((SRT**2-DM**2+AVMASS**2)/(2.*SRT))**2-AVMASS**2
       IF(PF2.GT.0.)PRF=SQRT(PF2)
* (3) RELABLE FINAL STATE PARTICLES
* 37 D(++)+n-->N*(+)(14)+p
          IF(M12.EQ.37)THEN
          IF(iabs(LB(I1)).EQ.9)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=11
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=11
         E(I1)=DM
          ENDIF
         GO TO 207
          ENDIF
* 38 D(++)+n-->N*(+)(15)+p
          IF(M12.EQ.38)THEN
          IF(iabs(LB(I1)).EQ.9)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 39 D(+)+P-->N*(+)(14)+p
          IF(M12.EQ.39)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=11
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=11
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 40 D(+)+P-->N*(+)(15)+p
          IF(M12.EQ.40)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 41 D(+)+N-->N*(+)(14)+N
          IF(M12.EQ.41)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=11
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=11
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 42 D(+)+N-->N*(+)(15)+N
          IF(M12.EQ.42)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 43 D(+)+N-->N*(0)(14)+P
          IF(M12.EQ.43)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=10
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=10
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 44 D(+)+N-->N*(0)(15)+P
          IF(M12.EQ.44)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=12
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 46 D(-)+P-->N*(0)(14)+N
          IF(M12.EQ.46)THEN
          IF(iabs(LB(I1)).EQ.6)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=10
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=10
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 47 D(-)+P-->N*(0)(15)+N
          IF(M12.EQ.47)THEN
          IF(iabs(LB(I1)).EQ.6)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=12
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 48 D(0)+N-->N*(0)(14)+N
          IF(M12.EQ.48)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=11
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=11
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 49 D(0)+N-->N*(0)(15)+N
          IF(M12.EQ.49)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=12
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 50 D(0)+P-->N*(0)(14)+P
          IF(M12.EQ.50)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=10
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=10
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 51 D(0)+P-->N*(+)(14)+N
          IF(M12.EQ.51)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=11
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=11
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 52 D(0)+P-->N*(0)(15)+P
          IF(M12.EQ.52)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=12
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 53 D(0)+P-->N*(+)(15)+N
          IF(M12.EQ.53)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 54 N*(0)(14)+P-->N*(+)(15)+N
          IF(M12.EQ.54)THEN
          IF(iabs(LB(I1)).EQ.10)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 55 N*(0)(14)+P-->N*(0)(15)+P
          IF(M12.EQ.55)THEN
          IF(iabs(LB(I1)).EQ.10)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=12
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 56 N*(+)(14)+N-->N*(+)(15)+N
          IF(M12.EQ.56)THEN
          IF(iabs(LB(I1)).EQ.11)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
* 57 N*(+)(14)+N-->N*(0)(15)+P
          IF(M12.EQ.57)THEN
          IF(iabs(LB(I1)).EQ.11)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=12
         E(I1)=DM
          ENDIF
         ENDIF
          GO TO 207
*------------------------------------------------
* RELABLE NUCLEONS AFTER DELTA OR N* BEING ABSORBED
*(1) n+delta(+)-->n+p
206       IF(M12.EQ.1)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I2)=2
          LB(I1)=1
          E(I1)=AMP
          ELSE
          LB(I1)=2
          LB(I2)=1
          E(I2)=AMP
          ENDIF
         GO TO 207
          ENDIF
*(2) p+delta(0)-->p+n
          IF(M12.EQ.2)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I2)=1
          LB(I1)=2
          E(I1)=AMN
          ELSE
          LB(I1)=1
          LB(I2)=2
          E(I2)=AMN
          ENDIF
         GO TO 207
          ENDIF
*(3) n+delta(++)-->p+p
          IF(M12.EQ.3)THEN
          LB(I1)=1
          LB(I2)=1
          E(I1)=AMP
          E(I2)=AMP
         GO TO 207
          ENDIF
*(4) p+delta(+)-->p+p
          IF(M12.EQ.4)THEN
          LB(I1)=1
          LB(I2)=1
          E(I1)=AMP
          E(I2)=AMP
         GO TO 207
          ENDIF
*(5) n+delta(0)-->n+n
          IF(M12.EQ.5)THEN
          LB(I1)=2
          LB(I2)=2
          E(I1)=AMN
          E(I2)=AMN
         GO TO 207
          ENDIF
*(6) p+delta(-)-->n+n
          IF(M12.EQ.6)THEN
          LB(I1)=2
          LB(I2)=2
          E(I1)=AMN
          E(I2)=AMN
         GO TO 207
          ENDIF
*(7) p+N*(0)-->n+p
          IF(M12.EQ.7)THEN
          IF(iabs(LB(I1)).EQ.1)THEN
          LB(I1)=1
          LB(I2)=2
          E(I1)=AMP
          E(I2)=AMN
          ELSE
          LB(I1)=2
          LB(I2)=1
          E(I1)=AMN
          E(I2)=AMP
          ENDIF
         GO TO 207
          ENDIF
*(8) n+N*(+)-->n+p
          IF(M12.EQ.8)THEN
          IF(iabs(LB(I1)).EQ.2)THEN
          LB(I1)=2
          LB(I2)=1
          E(I1)=AMN
          E(I2)=AMP
          ELSE
          LB(I1)=1
          LB(I2)=2
          E(I1)=AMP
          E(I2)=AMN
          ENDIF
         GO TO 207
          ENDIF
clin-6/2008
c*(9) N*(+)p-->pp
*(9) N*(+)(1535) p-->pp
          IF(M12.EQ.9)THEN
          LB(I1)=1
          LB(I2)=1
          E(I1)=AMP
          E(I2)=AMP
         GO TO 207
         ENDIF
*(12) N*(0)P-->nP
          IF(M12.EQ.12)THEN
          LB(I1)=2
          LB(I2)=1
          E(I1)=AMN
          E(I2)=AMP
         GO TO 207
         ENDIF
*(11) N*(+)n-->nP
          IF(M12.EQ.11)THEN
          LB(I1)=2
          LB(I2)=1
          E(I1)=AMN
          E(I2)=AMP
         GO TO 207
         ENDIF
clin-6/2008
c*(12) N*(0)p-->Np
*(12) N*(0)(1535) p-->Np
          IF(M12.EQ.12)THEN
          LB(I1)=1
          LB(I2)=2
          E(I1)=AMP
          E(I2)=AMN
         ENDIF
*----------------------------------------------
207       PR   = PRF
          C1   = 1.0 - 2.0 * RANART(NSEED)
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
            write(99,*) 'scheck4: ', scheck
            scheck=0.
         endif
         c1=sqrt(scheck)/pr
c         c1=sqrt(pr**2-cc1**2)/pr
         endif
          T1   = 2.0 * PI * RANART(NSEED)
          IBLOCK=3
      ENDIF
      if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
         lb(i1) = -lb(i1)
         lb(i2) = -lb(i2)
      endif
*-----------------------------------------------------------------------
*COM: SET THE NEW MOMENTUM COORDINATES
 107  IF(PX .EQ. 0.0 .AND. PY .EQ. 0.0) THEN
         T2 = 0.0
      ELSE
         T2=ATAN2(PY,PX)
      END IF
clin-9/2012: check argument in sqrt():
      scheck=1.0 - C1**2
      if(scheck.lt.0) then
         write(99,*) 'scheck5: ', scheck
         scheck=0.
      endif
      S1=SQRT(scheck)
c      S1   = SQRT( 1.0 - C1**2 )
clin-9/2012: check argument in sqrt():
      scheck=1.0 - C2**2
      if(scheck.lt.0) then
         write(99,*) 'scheck6: ', scheck
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
* FOR THE NN-->KAON+X PROCESS, FIND MOMENTUM OF THE FINAL PARTICLES IN 
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
csp11/21/01 end
                IBLOCK=11
                if(ianti .eq. 1)iblock=-11
c
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
              if(xsk1/(xsk1+xsk2+xsk3).
     1          gt.RANART(NSEED))then
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
128        CALL BBKAON(ic,SRT,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 128
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
c2008        X01 = 1.0 - 2.0 * RANART(NSEED)
c            Y01 = 1.0 - 2.0 * RANART(NSEED)
c            Z01 = 1.0 - 2.0 * RANART(NSEED)
c        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2008
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
                if(LPION(NNN,IRUN) .ne. 29) IBLOCK=11
        LB1=LB(I1)
        LB2=LB(I2)
        AM1=EM1
       am2=em2
        E1= SQRT( EM1**2 + PX1**2 + PY1**2 + PZ1**2 )
       RETURN
clin-6/2008 N+D->Deuteron+pi:
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
clin-6/2008 N+D->Deuteron+pi over
      END
**********************************
*                                                                      *
*                                                                      *
