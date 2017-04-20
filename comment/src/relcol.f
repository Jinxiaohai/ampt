      SUBROUTINE RELCOL(LCOLL,LBLOC,LCNNE,LDD,LPP,lppk,
     &LPN,lpd,lrho,lomega,LKN,LNNK,LDDK,LNDK,LCNND,LCNDN,
     &LDIRT,LDECAY,LRES,LDOU,LDDRHO,LNNRHO,LNNOM,
     &NT,ntmax,sp,akaon,sk)
*                                                                      *
*       PURPOSE:    CHECK CONDITIONS AND CALCULATE THE KINEMATICS      * 
*                   FOR BINARY COLLISIONS AMONG PARTICLES              *
*                                 - RELATIVISTIC FORMULA USED          *
*                                                                      *
*       REFERENCES: HAGEDORN, RELATIVISTIC KINEMATICS (1963)           *
*                                                                      *
*       VARIABLES:                                                     *
*         MASSPR  - NUMBER OF NUCLEONS IN PROJECTILE   (INTEGER,INPUT) *
*         MASSTA  - NUMBER OF NUCLEONS IN TARGET       (INTEGER,INPUT) *
*         NUM     - NUMBER OF TESTPARTICLES PER NUCLEON(INTEGER,INPUT) *
*         ISEED   - SEED FOR RANDOM NUMBER GENERATOR   (INTEGER,INPUT) *
*         IAVOID  - (= 1 => AVOID FIRST CLLISIONS WITHIN THE SAME      *
*                   NUCLEUS, ELSE ALL COLLISIONS)      (INTEGER,INPUT) *
*         DELTAR  - MAXIMUM SPATIAL DISTANCE FOR WHICH A COLLISION     *
*                   STILL CAN OCCUR                       (REAL,INPUT) *
*         DT      - TIME STEP SIZE                        (REAL,INPUT) *
*         LCOLL   - NUMBER OF COLLISIONS              (INTEGER,OUTPUT) *
*         LBLOC   - NUMBER OF PULI-BLOCKED COLLISIONS (INTEGER,OUTPUT) *
*         LCNNE   - NUMBER OF ELASTIC COLLISION       (INTEGER,OUTPUT) *
*         LCNND   - NUMBER OF N+N->N+DELTA REACTION   (INTEGER,OUTPUT) *
*         LCNDN   - NUMBER OF N+DELTA->N+N REACTION   (INTEGER,OUTPUT) *
*         LDD     - NUMBER OF RESONANCE+RESONANCE COLLISIONS
*         LPP     - NUMBER OF PION+PION elastic COLIISIONS
*         lppk    - number of pion(RHO,OMEGA)+pion(RHO,OMEGA)
*                   -->K+K- collisions
*         LPN     - NUMBER OF PION+N-->KAON+X
*         lpd     - number of pion+n-->delta+pion
*         lrho    - number of pion+n-->Delta+rho
*         lomega  - number of pion+n-->Delta+omega
*         LKN     - NUMBER OF KAON RESCATTERINGS
*         LNNK    - NUMBER OF bb-->kAON PROCESS
*         LDDK    - NUMBER OF DD-->KAON PROCESS
*         LNDK    - NUMBER OF ND-->KAON PROCESS
*         LB(I) IS USED TO LABEL PARTICLE'S CHARGE STATE
*         LB(I)   = 
cbali2/7/99 
*                 -45 Omega baryon(bar)
*                 -41 cascade0(bar)
*                 -40 cascade-(bar)
clin-11/07/00:
*                 -30 K*-
*                 -17 sigma+(bar)
*                 -16 sigma0(bar)
*                 -15 sigma-(bar)
*                 -14 LAMBDA(bar)
clin-8/29/00
*                 -13 anti-N*(+1)(1535),s_11
*                 -12 anti-N*0(1535),s_11
*                 -11 anti-N*(+1)(1440),p_11
*                 -10 anti-N*0(1440), p_11
*                  -9 anti-DELTA+2
*                  -8 anti-DELTA+1
*                  -7 anti-DELTA0
*                  -6 anti-DELTA-1
*
*                  -2 antineutron 
*                  -1 antiproton
cbali2/7/99end 
*                   0 eta
*                   1 PROTON
*                   2 NUETRON
*                   3 PION-
*                   4 PION0
*                   5 PION+          
*                   6 DELTA-1
*                   7 DELTA0
*                   8 DELTA+1
*                   9 DELTA+2
*                   10 N*0(1440), p_11
*                   11 N*(+1)(1440),p_11
*                  12 N*0(1535),s_11
*                  13 N*(+1)(1535),s_11
*                  14 LAMBDA
*                   15 sigma-
*                   16 sigma0
*                   17 sigma+
*                   21 kaon-
clin-2/23/03        22 Kaon0Long (converted at the last timestep)
*                   23 KAON+
*                   24 Kaon0short (converted at the last timestep then decay)
*                   25 rho-
*                   26 rho0
*                   27 rho+
*                   28 omega meson
*                   29 phi
*                   30 K*+
* sp01/03/01
*                   31 eta-prime
*                   40 cascade-
*                   41 cascade0
*                   45 Omega baryon
* sp01/03/01 end
*                   
*                   ++  ------- SEE NOTE BOOK
*         NSTAR=1 INCLUDING N* RESORANCE
*         ELSE DELTA RESORANCE ONLY
*         NDIRCT=1 INCLUDING DIRECT PROCESS,ELSE NOT
*         DIR - PERCENTAGE OF DIRECT PION PRODUCTION PROCESS
**********************************
      PARAMETER      (MAXSTR=150001,MAXR=1,PI=3.1415926)
      parameter      (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
      PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974,aks=0.895)
      PARAMETER      (AA1=1.26,APHI=1.02,AP1=0.13496)
      parameter            (maxx=20,maxz=24)
      parameter            (rrkk=0.6,prkk=0.3,srhoks=5.,ESBIN=0.04)
      DIMENSION MASSRN(0:MAXR),RT(3,MAXSTR),PT(3,MAXSTR),ET(MAXSTR)
      DIMENSION LT(MAXSTR), PROT(MAXSTR)
      COMMON   /AA/  R(3,MAXSTR)
cc      SAVE /AA/
      COMMON   /BB/  P(3,MAXSTR)
cc      SAVE /BB/
      COMMON   /CC/  E(MAXSTR)
cc      SAVE /CC/
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
cc      SAVE /DD/
      COMMON   /EE/  ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      COMMON   /HH/  PROPER(MAXSTR)
cc      SAVE /HH/
      common /ff/f(-mx:mx,-my:my,-mz:mz,-mpx:mpx,-mpy:mpy,-mpz:mpzp)
cc      SAVE /ff/
      common   /gg/  dx,dy,dz,dpx,dpy,dpz
cc      SAVE /gg/
      COMMON   /INPUT/ NSTAR,NDIRCT,DIR
cc      SAVE /INPUT/
      COMMON   /NN/NNN
cc      SAVE /NN/
      COMMON   /RR/  MASSR(0:MAXR)
cc      SAVE /RR/
      common   /ss/  inout(20)
cc      SAVE /ss/
      COMMON   /BG/BETAX,BETAY,BETAZ,GAMMA
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
      COMMON   /PE/PROPI(MAXSTR,MAXR)
cc      SAVE /PE/
      COMMON   /KKK/TKAON(7),EKAON(7,0:2000)
cc      SAVE /KKK/
      COMMON  /KAON/    AK(3,50,36),SPECK(50,36,7),MF
cc      SAVE /KAON/
      COMMON/TABLE/ xarray(0:1000),earray(0:1000)
cc      SAVE /TABLE/
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1 px1n,py1n,pz1n,dp1n
cc      SAVE /leadng/
      COMMON/tdecay/tfdcy(MAXSTR),tfdpi(MAXSTR,MAXR),tft(MAXSTR)
cc      SAVE /tdecay/
      common /lastt/itimeh,bimp 
cc      SAVE /lastt/
c
      COMMON/ppbmas/niso(15),nstate,ppbm(15,2),thresh(15),weight(15)
cc      SAVE /ppbmas/
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
cc      SAVE /ppb1/
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
cc      SAVE /ppmm/
      COMMON/hbt/lblast(MAXSTR),xlast(4,MAXSTR),plast(4,MAXSTR),nlast
cc      SAVE /hbt/
      common/resdcy/NSAV,iksdcy
cc      SAVE /resdcy/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      COMMON/FTMAX/ftsv(MAXSTR),ftsvt(MAXSTR, MAXR)
      dimension ftpisv(MAXSTR,MAXR),fttemp(MAXSTR)
      common /dpi/em2,lb2
      common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
clin-5/2008:
      DIMENSION dptemp(MAXSTR)
      common /para8/ idpert,npertd,idxsec
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
c
      real zet(-45:45)
      SAVE   
      data zet /
     4     1.,0.,0.,0.,0.,
     3     1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     2     -1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     1     0.,0.,0.,-1.,0.,1.,0.,-1.,0.,-1.,
     s     0.,-2.,-1.,0.,1.,0.,0.,0.,0.,-1.,
     e     0.,
     s     1.,0.,-1.,0.,1.,-1.,0.,1.,2.,0.,
     1     1.,0.,1.,0.,-1.,0.,1.,0.,0.,0.,
     2     -1.,0.,1.,0.,-1.,0.,1.,0.,0.,1.,
     3     0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,
     4     0.,0.,0.,0.,-1./
clin-2/19/03 initialize n and nsav for resonance decay at each timestep
c     in order to prevent integer overflow:
      call inidcy
c OFF skip ART collisions to reproduce HJ:      
cc       if(nt.ne.ntmax) return
clin-11/07/00 rrkk is assumed to be 0.6mb(default) for mm->KKbar 
c     with m=rho or omega, estimated from Ko's paper:
c      rrkk=0.6
c prkk: cross section of pi (rho or omega) -> K* Kbar (AND) K*bar K:
c      prkk=0.3
c     cross section in mb for (rho or omega) K* -> pi K:
c      srhoks=5.
clin-11/07/00-end
c      ESBIN=0.04
      RESONA=5.
*-----------------------------------------------------------------------
*     INITIALIZATION OF COUNTING VARIABLES
      NODELT=0
      SUMSRT =0.
      LCOLL  = 0
      LBLOC  = 0
      LCNNE  = 0
      LDD  = 0
      LPP  = 0
      lpd  = 0
      lpdr=0
      lrho = 0
      lrhor=0
      lomega=0
      lomgar=0
      LPN  = 0
      LKN  = 0
      LNNK = 0
      LDDK = 0
      LNDK = 0
      lppk =0
      LCNND  = 0
      LCNDN  = 0
      LDIRT  = 0
      LDECAY = 0
      LRES   = 0
      Ldou   = 0
      LDDRHO = 0
      LNNRHO = 0
      LNNOM  = 0
      MSUM   = 0
      MASSRN(0)=0
* COM: MSUM IS USED TO COUNT THE TOTAL NO. OF PARTICLES 
*      IN PREVIOUS IRUN-1 RUNS
* KAON COUNTERS
      DO 1002 IL=1,5
         TKAON(IL)=0
         DO 1001 IS=1,2000
            EKAON(IL,IS)=0
 1001    CONTINUE
 1002 CONTINUE
c sp 12/19/00
      DO 1004 i =1,NUM
         DO 1003 j =1,MAXSTR
            PROPI(j,i) = 1.
 1003    CONTINUE
 1004 CONTINUE
      do 1102 i=1,maxstr
         fttemp(i)=0.
         do 1101 irun=1,maxr
            ftpisv(i,irun)=0.
 1101    continue
 1102 continue
c sp 12/19/00 end
      sp=0
* antikaon counters
      akaon=0
      sk=0
*-----------------------------------------------------------------------
*     LOOP OVER ALL PARALLEL RUNS
cbz11/17/98
c      MASS=MASSPR+MASSTA
      MASS = 0
cbz11/17/98end
      DO 1000 IRUN = 1,NUM
         NNN=0
         MSUM=MSUM+MASSR(IRUN-1)
*     LOOP OVER ALL PSEUDOPARTICLES 1 IN THE SAME RUN
         J10=2
         IF(NT.EQ.NTMAX)J10=1
c
ctest off skips the check of energy conservation after each timestep:
c         enetot=0.
c         do ip=1,MASSR(IRUN)
c            if(e(ip).ne.0.or.lb(ip).eq.10022) enetot=enetot
c     1           +sqrt(p(1,ip)**2+p(2,ip)**2+p(3,ip)**2+e(ip)**2)
c         enddo
c         write(91,*) 'A:',nt,enetot,massr(irun),bimp 
         DO 800 J1 = J10,MASSR(IRUN)
            I1  = J1 + MSUM
* E(I)=0 are for pions having been absorbed or photons which do not enter here:
clin-4/2012 option of pi0 decays:
c            IF(E(I1).EQ.0.)GO TO 800
            IF(E(I1).EQ.0.)GO TO 798
c     To include anti-(Delta,N*1440 and N*1535):
c          IF ((LB(I1) .LT. -13 .OR. LB(I1) .GT. 28)
c     1         .and.iabs(LB(I1)) .ne. 30 ) GOTO 800
clin-4/2012 option of pi0 decays:
c            IF (LB(I1) .LT. -45 .OR. LB(I1) .GT. 45) GOTO 800
            IF (LB(I1) .LT. -45 .OR. LB(I1) .GT. 45) GOTO 798
            X1  = R(1,I1)
            Y1  = R(2,I1)
            Z1  = R(3,I1)
            PX1 = P(1,I1)
            PY1 = P(2,I1)
            PZ1 = P(3,I1)
            EM1 = E(I1)
            am1= em1
            E1  = SQRT( EM1**2 + PX1**2 + PY1**2 + PZ1**2 )
            ID1 = ID(I1)
            LB1 = LB(I1)
c     generate k0short and k0long from K+ and K- at the last timestep:
            if(nt.eq.ntmax.and.(lb1.eq.21.or.lb1.eq.23)) then
               pk0=RANART(NSEED)
               if(pk0.lt.0.25) then
                  LB(I1)=22
               elseif(pk0.lt.0.50) then
                  LB(I1)=24
               endif
               LB1=LB(I1)
            endif
clin-8/07/02 these particles don't decay strongly, so skip decay routines:     
c            IF( (lb1.ge.-2.and.lb1.le.5) .OR. lb1.eq.31 .OR.
c     &           (iabs(lb1).ge.14.and.iabs(lb1).le.24) .OR.
c     &           (iabs(lb1).ge.40.and.iabs(lb1).le.45) .or. 
c     &           lb1.eq.31)GO TO 1 
c     only decay K0short when iksdcy=1:
            if(lb1.eq.0.or.lb1.eq.25.or.lb1.eq.26.or.lb1.eq.27
     &           .or.lb1.eq.28.or.lb1.eq.29.or.iabs(lb1).eq.30
     &           .or.(iabs(lb1).ge.6.and.iabs(lb1).le.13)
     &           .or.(iksdcy.eq.1.and.lb1.eq.24)
     &           .or.iabs(lb1).eq.16
     &           .or.(ipi0dcy.eq.1.and.nt.eq.ntmax.and.lb1.eq.4)) then
clin-4/2012-above for option of pi0 decay:
c     &           .or.iabs(lb1).eq.16) then
               continue
            else
               goto 1
            endif
* IF I1 IS A RESONANCE, CHECK WHETHER IT DECAYS DURING THIS TIME STEP
         IF(lb1.ge.25.and.lb1.le.27) then
             wid=0.151
         ELSEIF(lb1.eq.28) then
             wid=0.00841
         ELSEIF(lb1.eq.29) then
             wid=0.00443
          ELSEIF(iabs(LB1).eq.30) then
             WID=0.051
         ELSEIF(lb1.eq.0) then
             wid=1.18e-6
c     to give K0short ct0=2.676cm:
         ELSEIF(iksdcy.eq.1.and.lb1.eq.24) then
             wid=7.36e-15
clin-4/29/03 add Sigma0 decay to Lambda, ct0=2.22E-11m:
         ELSEIF(iabs(lb1).eq.16) then
             wid=8.87e-6
csp-07/25/01 test a1 resonance:
cc          ELSEIF(LB1.EQ.32) then
cc             WID=0.40
          ELSEIF(LB1.EQ.32) then
             call WIDA1(EM1,rhomp,WID,iseed)
          ELSEIF(iabs(LB1).ge.6.and.iabs(LB1).le.9) then
             WID=WIDTH(EM1)
          ELSEIF((iabs(LB1).EQ.10).OR.(iabs(LB1).EQ.11)) then
             WID=W1440(EM1)
          ELSEIF((iabs(LB1).EQ.12).OR.(iabs(LB1).EQ.13)) then
             WID=W1535(EM1)
clin-4/2012 for option of pi0 decay:
          ELSEIF(ipi0dcy.eq.1.and.nt.eq.ntmax.and.lb1.eq.4) then
             wid=7.85e-9
          ENDIF
* if it is the last time step, FORCE all resonance to strong-decay
* and go out of the loop
          if(nt.eq.ntmax)then
             pdecay=1.1
clin-5b/2008 forbid phi decay at the end of hadronic cascade:
             if(iphidcy.eq.0.and.iabs(LB1).eq.29) pdecay=0.
ctest off clin-9/2012 forbid long-time decays (eta,omega,K*,Sigma0)
c     at the end of hadronic cascade to analyze freezeout time:
c             if(LB1.eq.0.or.LB1.eq.28.or.iabs(LB1).eq.30
c     1            .or.iabs(LB1).eq.16) pdecay=0.
          else
             T0=0.19733/WID
             GFACTR=E1/EM1
             T0=T0*GFACTR
             IF(T0.GT.0.)THEN
                PDECAY=1.-EXP(-DT/T0)
             ELSE
                PDECAY=0.
             ENDIF
          endif
          XDECAY=RANART(NSEED)
cc dilepton production from rho0, omega, phi decay 
cc        if(lb1.eq.26 .or. lb1.eq.28 .or. lb1.eq.29)
cc     &   call dec_ceres(nt,ntmax,irun,i1)
cc
          IF(XDECAY.LT.PDECAY) THEN
clin-10/25/02 get rid of argument usage mismatch in rhocay():
             idecay=irun
             tfnl=nt*dt
clin-10/28/03 keep formation time of hadrons unformed at nt=ntmax-1:
             if(nt.eq.ntmax.and.ftsv(i1).gt.((ntmax-1)*dt)) 
     1            tfnl=ftsv(i1)
             xfnl=x1
             yfnl=y1
             zfnl=z1
* use PYTHIA to perform decays of eta,rho,omega,phi,K*,(K0s) and Delta:
             if(lb1.eq.0.or.lb1.eq.25.or.lb1.eq.26.or.lb1.eq.27
     &           .or.lb1.eq.28.or.lb1.eq.29.or.iabs(lb1).eq.30
     &           .or.(iabs(lb1).ge.6.and.iabs(lb1).le.9)
     &           .or.(iksdcy.eq.1.and.lb1.eq.24)
     &           .or.iabs(lb1).eq.16
     &           .or.(ipi0dcy.eq.1.and.nt.eq.ntmax.and.lb1.eq.4)) then
clin-4/2012 Above for option of pi0 decay:
c     &           .or.iabs(lb1).eq.16) then
c     previous rho decay performed in rhodecay():
c                nnn=nnn+1
c                call rhodecay(idecay,i1,nnn,iseed)
c
ctest off record decays of phi,K*,Lambda(1520) resonances:
c                if(lb1.eq.29.or.iabs(lb1).eq.30) 
c     1               write(18,112) 'decay',lb1,px1,py1,pz1,am1,nt
c
clin-4/2012 option of pi0 decays:
c                call resdec(i1,nt,nnn,wid,idecay)
                call resdec(i1,nt,nnn,wid,idecay,0)
                p(1,i1)=px1n
                p(2,i1)=py1n
                p(3,i1)=pz1n
clin-5/2008:
                dpertp(i1)=dp1n
c     add decay time to freezeout positions & time at the last timestep:
                if(nt.eq.ntmax) then
                   R(1,i1)=xfnl
                   R(2,i1)=yfnl
                   R(3,i1)=zfnl
                   tfdcy(i1)=tfnl
                endif
c
* decay number for baryon resonance or L/S decay
                if(iabs(lb1).ge.6.and.iabs(lb1).le.9) then
                   LDECAY=LDECAY+1
                endif
* for a1 decay 
c             elseif(lb1.eq.32)then
c                NNN=NNN+1
c                call a1decay(idecay,i1,nnn,iseed,rhomp)
* FOR N*(1440)
             elseif(iabs(LB1).EQ.10.OR.iabs(LB1).EQ.11) THEN
                NNN=NNN+1
                LDECAY=LDECAY+1
                PNSTAR=1.
                IF(E(I1).GT.1.22)PNSTAR=0.6
                IF(RANART(NSEED).LE.PNSTAR)THEN
* (1) DECAY TO SINGLE PION+NUCLEON
                   CALL DECAY(idecay,I1,NNN,ISEED,wid,nt)
                ELSE
* (2) DECAY TO TWO PIONS + NUCLEON
                   CALL DECAY2(idecay,I1,NNN,ISEED,wid,nt)
                   NNN=NNN+1
                ENDIF
c for N*(1535) decay
             elseif(iabs(LB1).eq.12.or.iabs(LB1).eq.13) then
                NNN=NNN+1
                CALL DECAY(idecay,I1,NNN,ISEED,wid,nt)
                LDECAY=LDECAY+1
             endif
c
*COM: AT HIGH ENERGIES WE USE VERY SHORT TIME STEPS,
*     IN ORDER TO TAKE INTO ACCOUNT THE FINITE FORMATIOM TIME, WE
*     DO NOT ALLOW PARTICLES FROM THE DECAY OF RESONANCE TO INTERACT 
*     WITH OTHERS IN THE SAME TIME STEP. CHANGE 9000 TO REVERSE THIS 
*     ASSUMPTION. EFFECTS OF THIS ASSUMPTION CAN BE STUDIED BY CHANGING 
*     THE STATEMENT OF 9000. See notebook for discussions on effects of
*     changing statement 9000.
c
c     kaons from K* decay are converted to k0short (and k0long), 
c     phi decay may produce rho, K0S or eta, N*(1535) decay may produce eta,
c     and these decay daughters need to decay again if at the last timestep:
c     (note: these daughters have been assigned to lb(i1) only, not to lpion)
c             if(nt.eq.ntmax.and.(lb1.eq.29.or.iabs(lb1).eq.30
c     1            .iabs(lb1).eq.12.or.iabs(lb1).eq.13)) then
             if(nt.eq.ntmax) then
                if(lb(i1).eq.25.or.lb(i1).eq.26.or.lb(i1).eq.27) then
                   wid=0.151
                elseif(lb(i1).eq.0) then
                   wid=1.18e-6
                elseif(lb(i1).eq.24.and.iksdcy.eq.1) then
clin-4/2012 corrected K0s decay width:
c                   wid=7.36e-17
                   wid=7.36e-15
clin-4/2012 option of pi0 decays:
                elseif(ipi0dcy.eq.1.and.lb(i1).eq.4) then
                   wid=7.85e-9
                else
                   goto 9000
                endif
                LB1=LB(I1)
                PX1=P(1,I1)
                PY1=P(2,I1)
                PZ1=P(3,I1)
                EM1=E(I1)
                E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
clin-4/2012 option of pi0 decays:
c                call resdec(i1,nt,nnn,wid,idecay)
                call resdec(i1,nt,nnn,wid,idecay,0)
                p(1,i1)=px1n
                p(2,i1)=py1n
                p(3,i1)=pz1n
                R(1,i1)=xfnl
                R(2,i1)=yfnl
                R(3,i1)=zfnl
                tfdcy(i1)=tfnl
clin-5/2008:
                dpertp(i1)=dp1n
             endif
c     Decay daughter of the above decay in lb(i1) may be a pi0:
             if(nt.eq.ntmax.and.ipi0dcy.eq.1.and.lb(i1).eq.4) then
                wid=7.85e-9
                LB1=LB(I1)
                PX1=P(1,I1)
                PY1=P(2,I1)
                PZ1=P(3,I1)
                EM1=E(I1)
                E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
                call resdec(i1,nt,nnn,wid,idecay,0)
                p(1,i1)=px1n
                p(2,i1)=py1n
                p(3,i1)=pz1n
                R(1,i1)=xfnl
                R(2,i1)=yfnl
                R(3,i1)=zfnl
                tfdcy(i1)=tfnl
                dpertp(i1)=dp1n
             endif
* negelecting the Pauli blocking at high energies
clin-4/2012 option of pi0 decays:
c 9000        go to 800
 9000        go to 798
          ENDIF
* LOOP OVER ALL PSEUDOPARTICLES 2 IN THE SAME RUN
* SAVE ALL THE COORDINATES FOR POSSIBLE CHANGE IN THE FOLLOWING COLLISION
clin-4/2012 option of pi0 decays:
c 1        if(nt.eq.ntmax)go to 800
 1        if(nt.eq.ntmax)go to 798
          X1 = R(1,I1)
          Y1 = R(2,I1)
          Z1 = R(3,I1)
c
           DO 600 J2 = 1,J1-1
            I2  = J2 + MSUM
* IF I2 IS A MESON BEING ABSORBED, THEN GO OUT OF THE LOOP
            IF(E(I2).EQ.0.) GO TO 600
clin-5/2008 in case the first particle is already destroyed:
            IF(E(I1).EQ.0.) GO TO 800
clin-4/2012 option of pi0 decays:
            IF (LB(I2) .LT. -45 .OR. LB(I2) .GT. 45) GOTO 600
clin-7/26/03 improve speed
            X2=R(1,I2)
            Y2=R(2,I2)
            Z2=R(3,I2)
            dr0max=5.
clin-9/2008 deuteron+nucleon elastic cross sections could reach ~2810mb:
            ilb1=iabs(LB(I1))
            ilb2=iabs(LB(I2))
            IF(ilb1.EQ.42.or.ilb2.EQ.42) THEN
               if((ILB1.GE.1.AND.ILB1.LE.2)
     1              .or.(ILB1.GE.6.AND.ILB1.LE.13)
     2              .or.(ILB2.GE.1.AND.ILB2.LE.2)
     3              .or.(ILB2.GE.6.AND.ILB2.LE.13)) then
                  if((lb(i1)*lb(i2)).gt.0) dr0max=10.
               endif
            ENDIF
c
            if(((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2).GT.dr0max**2)
     1           GO TO 600
            IF (ID(I1)*ID(I2).EQ.IAVOID) GOTO 400
            ID1=ID(I1)
            ID2 = ID(I2)
c
            ix1= nint(x1/dx)
            iy1= nint(y1/dy)
            iz1= nint(z1/dz)
            PX1=P(1,I1)
            PY1=P(2,I1)
            PZ1=P(3,I1)
            EM1=E(I1)
            AM1=EM1
            LB1=LB(I1)
            E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
            IPX1=NINT(PX1/DPX)
            IPY1=NINT(PY1/DPY)
            IPZ1=NINT(PZ1/DPZ)         
            LB2 = LB(I2)
            PX2 = P(1,I2)
            PY2 = P(2,I2)
            PZ2 = P(3,I2)
            EM2=E(I2)
            AM2=EM2
            lb1i=lb(i1)
            lb2i=lb(i2)
            px1i=P(1,I1)
            py1i=P(2,I1)
            pz1i=P(3,I1)
            em1i=E(I1)
            px2i=P(1,I2)
            py2i=P(2,I2)
            pz2i=P(3,I2)
            em2i=E(I2)
clin-2/26/03 ctest off check energy conservation after each binary search:
            eini=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2)
     1           +SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2)
            pxini=P(1,I1)+P(1,I2)
            pyini=P(2,I1)+P(2,I2)
            pzini=P(3,I1)+P(3,I2)
            nnnini=nnn
c
clin-4/30/03 initialize value:
            iblock=0
c
* TO SAVE COMPUTING TIME we do the following
* (1) make a ROUGH estimate to see whether particle i2 will collide with
* particle I1, and (2) skip the particle pairs for which collisions are 
* not modeled in the code.
* FOR MESON-BARYON AND MESON-MESON COLLISIONS, we use a maximum 
* interaction distance DELTR0=2.6
* for ppbar production from meson (pi rho omega) interactions:
c
            DELTR0=3.
        if( (iabs(lb1).ge.14.and.iabs(lb1).le.17) .or.
     &      (iabs(lb1).ge.30.and.iabs(lb1).le.45) ) DELTR0=5.0
        if( (iabs(lb2).ge.14.and.iabs(lb2).le.17) .or.
     &      (iabs(lb2).ge.30.and.iabs(lb2).le.45) ) DELTR0=5.0
            if(lb1.eq.28.and.lb2.eq.28) DELTR0=4.84
clin-10/08/00 to include pi pi -> rho rho:
            if((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.3.and.lb2.le.5)) then
               E2=SQRT(EM2**2+PX2**2+PY2**2+PZ2**2)
         spipi=(e1+e2)**2-(px1+px2)**2-(py1+py2)**2-(pz1+pz2)**2
               if(spipi.ge.(4*0.77**2)) DELTR0=3.5
            endif
c khyperon
        IF (LB1.EQ.23 .AND. (LB2.GE.14.AND.LB2.LE.17)) GOTO 3699
        IF (LB2.EQ.23 .AND. (LB1.GE.14.AND.LB1.LE.17)) GOTO 3699
* K(K*) + Kbar(K*bar) scattering including 
*     K(K*) + Kbar(K*bar) --> phi + pi(rho,omega) and pi pi(rho,omega)
       if(lb1.eq.21.and.lb2.eq.23)go to 3699
       if(lb2.eq.21.and.lb1.eq.23)go to 3699
       if(lb1.eq.30.and.lb2.eq.21)go to 3699
       if(lb2.eq.30.and.lb1.eq.21)go to 3699
       if(lb1.eq.-30.and.lb2.eq.23)go to 3699
       if(lb2.eq.-30.and.lb1.eq.23)go to 3699
       if(lb1.eq.-30.and.lb2.eq.30)go to 3699
       if(lb2.eq.-30.and.lb1.eq.30)go to 3699
c
clin-12/15/00
c     kaon+rho(omega,eta) collisions:
      if(lb1.eq.21.or.lb1.eq.23) then
         if(lb2.eq.0.or.(lb2.ge.25.and.lb2.le.28)) then
            go to 3699
         endif
      elseif(lb2.eq.21.or.lb2.eq.23) then
         if(lb1.eq.0.or.(lb1.ge.25.and.lb1.le.28)) then
            goto 3699
         endif
      endif
clin-8/14/02 K* (pi, rho, omega, eta) collisions:
      if(iabs(lb1).eq.30 .and.
     1     (lb2.eq.0.or.(lb2.ge.25.and.lb2.le.28)
     2     .or.(lb2.ge.3.and.lb2.le.5))) then
         go to 3699
      elseif(iabs(lb2).eq.30 .and.
     1        (lb1.eq.0.or.(lb1.ge.25.and.lb1.le.28)
     2        .or.(lb1.ge.3.and.lb1.le.5))) then
         goto 3699
clin-8/14/02-end
c K*/K*-bar + baryon/antibaryon collisions:
        elseif( iabs(lb1).eq.30 .and.
     1     (iabs(lb2).eq.1.or.iabs(lb2).eq.2.or.
     2     (iabs(lb2).ge.6.and.iabs(lb2).le.13)) )then
              go to 3699
           endif
         if( iabs(lb2).eq.30 .and.
     1         (iabs(lb1).eq.1.or.iabs(lb1).eq.2.or.
     2         (iabs(lb1).ge.6.and.iabs(lb1).le.13)) )then
                go to 3699
        endif                                                              
* K^+ baryons and antibaryons:
c** K+ + B-bar  --> La(Si)-bar + pi
* K^- and antibaryons, note K^- and baryons are included in newka():
* note that we fail to satisfy charge conjugation for these cross sections:
        if((lb1.eq.23.or.lb1.eq.21).and.
     1       (iabs(lb2).eq.1.or.iabs(lb2).eq.2.or.
     2       (iabs(lb2).ge.6.and.iabs(lb2).le.13))) then
           go to 3699
        elseif((lb2.eq.23.or.lb2.eq.21).and.
     1       (iabs(lb1).eq.1.or.iabs(lb1).eq.2.or.
     2       (iabs(lb1).ge.6.and.iabs(lb1).le.13))) then
           go to 3699
        endif
*
* For anti-nucleons annihilations:
* Assumptions: 
* (1) for collisions involving a p_bar or n_bar,
* we allow only collisions between a p_bar and a baryon or a baryon 
* resonance (as well as a n_bar and a baryon or a baryon resonance),
* we skip all other reactions involving a p_bar or n_bar, 
* such as collisions between p_bar (n_bar) and mesons, 
* and collisions between two p_bar's (n_bar's). 
* (2) we introduce a new parameter rppmax: the maximum interaction 
* distance to make the quick collision check,rppmax=3.57 fm 
* corresponding to a cutoff of annihilation xsection= 400mb which is
* also used consistently in the actual annihilation xsection to be 
* used in the following as given in the subroutine xppbar(srt)
        rppmax=3.57   
* anti-baryon on baryons
        if((lb1.eq.-1.or.lb1.eq.-2.or.(lb1.ge.-13.and.lb1.le.-6))
     1 .and.(lb2.eq.1.or.lb2.eq.2.or.(lb2.ge.6.and.lb2.le.13))) then
            DELTR0 = RPPMAX
            GOTO 2699
       else if((lb2.eq.-1.or.lb2.eq.-2.or.(lb2.ge.-13.and.lb2.le.-6))
     1 .and.(lb1.eq.1.or.lb1.eq.2.or.(lb1.ge.6.and.lb1.le.13))) then
            DELTR0 = RPPMAX
            GOTO 2699
         END IF
c*  ((anti) lambda, cascade, omega  should not be rejected)
        if( (iabs(lb1).ge.14.and.iabs(lb1).le.17) .or.
     &      (iabs(lb2).ge.14.and.iabs(lb2).le.17) )go to 3699
c
clin-9/2008 maximum sigma~2810mb for deuteron+nucleon elastic collisions:
         IF (iabs(LB1).EQ.42.or.iabs(LB2).EQ.42) THEN
            ilb1=iabs(LB1)
            ilb2=iabs(LB2)
            if((ILB1.GE.1.AND.ILB1.LE.2)
     1           .or.(ILB1.GE.6.AND.ILB1.LE.13)
     2           .or.(ILB2.GE.1.AND.ILB2.LE.2)
     3           .or.(ILB2.GE.6.AND.ILB2.LE.13)) then
               if((lb1*lb2).gt.0) deltr0=9.5
            endif
         ENDIF
c
        if( (iabs(lb1).ge.40.and.iabs(lb1).le.45) .or. 
     &      (iabs(lb2).ge.40.and.iabs(lb2).le.45) )go to 3699
c
c* phi channel --> elastic + inelastic scatt.  
         IF( (lb1.eq.29 .and.((lb2.ge.1.and.lb2.le.13).or.  
     &       (lb2.ge.21.and.lb2.le.28).or.iabs(lb2).eq.30)) .OR.
     &     (lb2.eq.29 .and.((lb1.ge.1.and.lb1.le.13).or.
     &       (lb1.ge.21.and.lb1.le.28).or.iabs(lb1).eq.30)) )THEN
             DELTR0=3.0
             go to 3699
        endif
c
c  La/Si, Cas, Om (bar)-meson elastic colln
* pion vs. La & Ca (bar) coll. are treated in resp. subroutines
* SKIP all other K* RESCATTERINGS
        If(iabs(lb1).eq.30.or.iabs(lb2).eq.30) go to 400
* SKIP KAON(+) RESCATTERINGS WITH particles other than pions and baryons 
         If(lb1.eq.23.and.(lb2.lt.1.or.lb2.gt.17))go to 400
         If(lb2.eq.23.and.(lb1.lt.1.or.lb1.gt.17))go to 400
c
c anti-baryon proccess: B-bar+M, N-bar+R-bar, N-bar+N-bar, R-bar+R-bar
c  R = (D,N*)
         if( ((lb1.le.-1.and.lb1.ge.-13)
     &        .and.(lb2.eq.0.or.(lb2.ge.3.and.lb2.le.5)
     &            .or.(lb2.ge.25.and.lb2.le.28))) 
     &      .OR.((lb2.le.-1.and.lb2.ge.-13)
     &         .and.(lb1.eq.0.or.(lb1.ge.3.and.lb1.le.5)
     &              .or.(lb1.ge.25.and.lb1.le.28))) ) then
         elseIF( ((LB1.eq.-1.or.lb1.eq.-2).
     &             and.(LB2.LT.-5.and.lb2.ge.-13))
     &      .OR. ((LB2.eq.-1.or.lb2.eq.-2).
     &             and.(LB1.LT.-5.and.lb1.ge.-13)) )then
         elseIF((LB1.eq.-1.or.lb1.eq.-2)
     &     .AND.(LB2.eq.-1.or.lb2.eq.-2))then
         elseIF((LB1.LT.-5.and.lb1.ge.-13).AND.
     &          (LB2.LT.-5.and.lb2.ge.-13)) then
c        elseif((lb1.lt.0).or.(lb2.lt.0)) then
c         go to 400
       endif               
 2699    CONTINUE
* for baryon-baryon collisions
         IF (LB1 .EQ. 1 .OR. LB1 .EQ. 2 .OR. (LB1 .GE. 6 .AND.
     &        LB1 .LE. 17)) THEN
            IF (LB2 .EQ. 1 .OR. LB2 .EQ. 2 .OR. (LB2 .GE. 6 .AND.
     &           LB2 .LE. 17)) THEN
               DELTR0 = 2.
            END IF
         END IF
c
 3699   RSQARE = (X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2
        IF (RSQARE .GT. DELTR0**2) GO TO 400
*NOW PARTICLES ARE CLOSE ENOUGH TO EACH OTHER !
* KEEP ALL COORDINATES FOR POSSIBLE PHASE SPACE CHANGE
            ix2 = nint(x2/dx)
            iy2 = nint(y2/dy)
            iz2 = nint(z2/dz)
            ipx2 = nint(px2/dpx)
            ipy2 = nint(py2/dpy)
            ipz2 = nint(pz2/dpz)
* FIND MOMENTA OF PARTICLES IN THE CMS OF THE TWO COLLIDING PARTICLES
* AND THE CMS ENERGY SRT
            CALL CMS(I1,I2,PCX,PCY,PCZ,SRT)
clin-7/26/03 improve speed
          drmax=dr0max
          call distc0(drmax,deltr0,DT,
     1         Ifirst,PCX,PCY,PCZ,
     2         x1,y1,z1,px1,py1,pz1,em1,x2,y2,z2,px2,py2,pz2,em2)
          if(Ifirst.eq.-1) goto 400
         ISS=NINT(SRT/ESBIN)
clin-4/2008 use last bin if ISS is out of EKAON's upper bound of 2000:
         if(ISS.gt.2000) ISS=2000
*Sort collisions
c
clin-8/2008 Deuteron+Meson->B+B; 
c     meson=(pi,rho,omega,eta), B=(n,p,Delta,N*1440,N*1535):
         IF (iabs(LB1).EQ.42.or.iabs(LB2).EQ.42) THEN
            ilb1=iabs(LB1)
            ilb2=iabs(LB2)
            if(LB1.eq.0.or.(LB1.GE.3.AND.LB1.LE.5)
     1           .or.(LB1.GE.25.AND.LB1.LE.28)
     2           .or.
     3           LB2.eq.0.or.(LB2.GE.3.AND.LB2.LE.5)
     4           .or.(LB2.GE.25.AND.LB2.LE.28)) then
               GOTO 505
clin-9/2008 Deuteron+Baryon or antiDeuteron+antiBaryon elastic collisions:
            elseif(((ILB1.GE.1.AND.ILB1.LE.2)
     1              .or.(ILB1.GE.6.AND.ILB1.LE.13)
     2              .or.(ILB2.GE.1.AND.ILB2.LE.2)
     3              .or.(ILB2.GE.6.AND.ILB2.LE.13))
     4              .and.(lb1*lb2).gt.0) then
               GOTO 506
            else
               GOTO 400
            endif
         ENDIF
c
* K+ + (N,N*,D)-bar --> L/S-bar + pi
          if( ((lb1.eq.23.or.lb1.eq.30).and.
     &         (lb2.eq.-1.or.lb2.eq.-2.or.(lb2.ge.-13.and.lb2.le.-6))) 
     &         .OR.((lb2.eq.23.or.lb2.eq.30).and.
     &         (lb1.eq.-1.or.lb1.eq.-2.or.(lb1.ge.-13.and.lb1.le.-6))) )
     &         then
             bmass=0.938
             if(srt.le.(bmass+aka)) then
                pkaon=0.
             else
                pkaon=sqrt(((srt**2-(aka**2+bmass**2))
     1               /2./bmass)**2-aka**2)
             endif
clin-10/31/02 cross sections are isospin-averaged, same as those in newka
c     for K- + (N,N*,D) --> L/S + pi:
             sigela = 0.5 * (AKPEL(PKAON) + AKNEL(PKAON))
             SIGSGM = 1.5 * AKPSGM(PKAON) + AKNSGM(PKAON)
             SIG = sigela + SIGSGM + AKPLAM(PKAON)
             if(sig.gt.1.e-7) then
c     ! K+ + N-bar reactions
                icase=3
                brel=sigela/sig
                brsgm=sigsgm/sig
                brsig = sig
                nchrg = 1
                go to 3555
             endif
             go to 400
          endif
c
c
c  meson + hyperon-bar -> K+ + N-bar
          if(((lb1.ge.-17.and.lb1.le.-14).and.(lb2.ge.3.and.lb2.le.5)) 
     &         .OR.((lb2.ge.-17.and.lb2.le.-14)
     &         .and.(lb1.ge.3.and.lb1.le.5)))then
             nchrg=-100
C*       first classify the reactions due to total charge.
             if((lb1.eq.-15.and.(lb2.eq.5.or.lb2.eq.27)).OR.
     &            (lb2.eq.-15.and.(lb1.eq.5.or.lb1.eq.27))) then
                nchrg=-2
c     ! D-(bar)
                bmass=1.232
                go to 110
             endif
             if( (lb1.eq.-15.and.(lb2.eq.0.or.lb2.eq.4.or.lb2.eq.26.or.
     &            lb2.eq.28)).OR.(lb2.eq.-15.and.(lb1.eq.0.or.
     &            lb1.eq.4.or.lb1.eq.26.or.lb1.eq.28)).OR.
     &   ((lb1.eq.-14.or.lb1.eq.-16).and.(lb2.eq.5.or.lb2.eq.27)).OR.
     &   ((lb2.eq.-14.or.lb2.eq.-16).and.(lb1.eq.5.or.lb1.eq.27)) )then
                nchrg=-1
c     ! n-bar
                bmass=0.938
                go to 110
             endif
             if(  (lb1.eq.-15.and.(lb2.eq.3.or.lb2.eq.25)).OR.
     &            (lb2.eq.-15.and.(lb1.eq.3.or.lb1.eq.25)).OR.
     &            (lb1.eq.-17.and.(lb2.eq.5.or.lb2.eq.27)).OR.
     &            (lb2.eq.-17.and.(lb1.eq.5.or.lb1.eq.27)).OR.
     &            ((lb1.eq.-14.or.lb1.eq.-16).and.(lb2.eq.0.or.lb2.eq.4
     &            .or.lb2.eq.26.or.lb2.eq.28)).OR.
     &            ((lb2.eq.-14.or.lb2.eq.-16).and.(lb1.eq.0.or.lb1.eq.4
     &            .or.lb1.eq.26.or.lb1.eq.28)) )then
               nchrg=0
c     ! p-bar
                bmass=0.938
                go to 110
             endif
             if( (lb1.eq.-17.and.(lb2.eq.0.or.lb2.eq.4.or.lb2.eq.26.or.
     &            lb2.eq.28)).OR.(lb2.eq.-17.and.(lb1.eq.0.or.
     &            lb1.eq.4.or.lb1.eq.26.or.lb1.eq.28)).OR.
     &  ((lb1.eq.-14.or.lb1.eq.-16).and.(lb2.eq.3.or.lb2.eq.25)).OR.
     &  ((lb2.eq.-14.or.lb2.eq.-16).and.(lb1.eq.3.or.lb1.eq.25)))then
               nchrg=1
c     ! D++(bar)
                bmass=1.232
             endif
c
c 110     if(nchrg.ne.-100.and.srt.ge.(aka+bmass))then !! for elastic
 110         sig = 0.
c !! for elastic
         if(nchrg.ne.-100.and.srt.ge.(aka+bmass))then
cc110        if(nchrg.eq.-100.or.srt.lt.(aka+bmass)) go to 400
c             ! PI + La(Si)-bar => K+ + N-bar reactions
            icase=4
cc       pkaon=sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
            pkaon=sqrt(((srt**2-(aka**2+0.938**2))/2./0.938)**2-aka**2)
c ! lambda-bar + Pi
            if(lb1.eq.-14.or.lb2.eq.-14) then
               if(nchrg.ge.0) sigma0=akPlam(pkaon)
               if(nchrg.lt.0) sigma0=akNlam(pkaon)
c                ! sigma-bar + pi
            else
c !K-p or K-D++
               if(nchrg.ge.0) sigma0=akPsgm(pkaon)
c !K-n or K-D-
               if(nchrg.lt.0) sigma0=akNsgm(pkaon)
               SIGMA0 = 1.5 * AKPSGM(PKAON) + AKNSGM(PKAON)
            endif
            sig=(srt**2-(aka+bmass)**2)*(srt**2-(aka-bmass)**2)/
     &           (srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)*sigma0
c ! K0barD++, K-D-
            if(nchrg.eq.-2.or.nchrg.eq.2) sig=2.*sig
C*     the factor 2 comes from spin of delta, which is 3/2
C*     detailed balance. copy from Page 423 of N.P. A614 1997
            IF (LB1 .EQ. -14 .OR. LB2 .EQ. -14) THEN
               SIG = 4.0 / 3.0 * SIG
            ELSE IF (NCHRG .EQ. -2 .OR. NCHRG .EQ. 2) THEN
               SIG = 8.0 / 9.0 * SIG
            ELSE
               SIG = 4.0 / 9.0 * SIG
            END IF
cc        brel=0.
cc        brsgm=0.
cc        brsig = sig
cc          if(sig.lt.1.e-7) go to 400
*-
         endif
c                ! PI + La(Si)-bar => elastic included
         icase=4
         sigela = 10.
         sig = sig + sigela
         brel= sigela/sig
         brsgm=0.
         brsig = sig
*-
         go to 3555
      endif
** MULTISTRANGE PARTICLE (Cas,Omega -bar) PRODUCTION - (NON)PERTURBATIVE
* K-/K*0bar + La/Si --> cascade + pi/eta
      if( ((lb1.eq.21.or.lb1.eq.-30).and.(lb2.ge.14.and.lb2.le.17)).OR.
     &  ((lb2.eq.21.or.lb2.eq.-30).and.(lb1.ge.14.and.lb1.le.17)) )then
          kp = 0
          go to 3455
        endif
c K+/K*0 + La/Si(bar) --> cascade-bar + pi/eta
      if( ((lb1.eq.23.or.lb1.eq.30).and.(lb2.le.-14.and.lb2.ge.-17)).OR.
     &  ((lb2.eq.23.or.lb2.eq.30).and.(lb1.le.-14.and.lb1.ge.-17)) )then
          kp = 1
          go to 3455
        endif
* K-/K*0bar + cascade --> omega + pi
       if( ((lb1.eq.21.or.lb1.eq.-30).and.(lb2.eq.40.or.lb2.eq.41)).OR.
     & ((lb2.eq.21.or.lb2.eq.-30).and.(lb1.eq.40.or.lb1.eq.41)) )then
          kp = 0
          go to 3455
        endif
* K+/K*0 + cascade-bar --> omega-bar + pi
       if( ((lb1.eq.23.or.lb1.eq.30).and.(lb2.eq.-40.or.lb2.eq.-41)).OR.
     &  ((lb2.eq.23.or.lb2.eq.30).and.(lb1.eq.-40.or.lb1.eq.-41)) )then
          kp = 1
          go to 3455
        endif
* Omega + Omega --> Di-Omega + photon(eta)
cc        if( lb1.eq.45.and.lb2.eq.45 ) go to 3455
c annhilation of cascade(bar), omega(bar)
         kp = 3
* K- + L/S <-- cascade(bar) + pi/eta
       if( (((lb1.ge.3.and.lb1.le.5).or.lb1.eq.0) 
     &       .and.(iabs(lb2).eq.40.or.iabs(lb2).eq.41))
     & .OR. (((lb2.ge.3.and.lb2.le.5).or.lb2.eq.0) 
     &       .and.(iabs(lb1).eq.40.or.iabs(lb1).eq.41)) )go to 3455
* K- + cascade(bar) <-- omega(bar) + pi
*         if(  (lb1.eq.0.and.iabs(lb2).eq.45)
*    &       .OR. (lb2.eq.0.and.iabs(lb1).eq.45) )go to 3455
        if( ((lb1.ge.3.and.lb1.le.5).and.iabs(lb2).eq.45)
     &  .OR.((lb2.ge.3.and.lb2.le.5).and.iabs(lb1).eq.45) )go to 3455
c
***  MULTISTRANGE PARTICLE PRODUCTION  (END)
c* K+ + La(Si) --> Meson + B
        IF (LB1.EQ.23 .AND. (LB2.GE.14.AND.LB2.LE.17)) GOTO 5699
        IF (LB2.EQ.23 .AND. (LB1.GE.14.AND.LB1.LE.17)) GOTO 5699
c* K- + La(Si)-bar --> Meson + B-bar
       IF (LB1.EQ.21 .AND. (LB2.GE.-17.AND.LB2.LE.-14)) GOTO 5699
       IF (LB2.EQ.21 .AND. (LB1.GE.-17.AND.LB1.LE.-14)) GOTO 5699
c La/Si-bar + B --> pi + K+
       IF( (((LB1.eq.1.or.LB1.eq.2).or.(LB1.ge.6.and.LB1.le.13))
     &       .AND.(LB2.GE.-17.AND.LB2.LE.-14)) .OR.
     &     (((LB2.eq.1.or.LB2.eq.2).or.(LB2.ge.6.and.LB2.le.13))
     &      .AND.(LB1.GE.-17.AND.LB1.LE.-14)) )go to 5999
c La/Si + B-bar --> pi + K-
       IF( (((LB1.eq.-1.or.LB1.eq.-2).or.(LB1.le.-6.and.LB1.ge.-13))
     &       .AND.(LB2.GE.14.AND.LB2.LE.17)) .OR.
     &     (((LB2.eq.-1.or.LB2.eq.-2).or.(LB2.le.-6.and.LB2.ge.-13))
     &       .AND.(LB1.GE.14.AND.LB1.LE.17)) )go to 5999 
*
*
* K(K*) + Kbar(K*bar) --> phi + pi(rho,omega), M + M (M=pi,rho,omega,eta)
       if(lb1.eq.21.and.lb2.eq.23) go to 8699
       if(lb2.eq.21.and.lb1.eq.23) go to 8699
       if(lb1.eq.30.and.lb2.eq.21) go to 8699
       if(lb2.eq.30.and.lb1.eq.21) go to 8699
       if(lb1.eq.-30.and.lb2.eq.23) go to 8699
       if(lb2.eq.-30.and.lb1.eq.23) go to 8699
       if(lb1.eq.-30.and.lb2.eq.30) go to 8699
       if(lb2.eq.-30.and.lb1.eq.30) go to 8699
c* (K,K*)-bar + rho(omega) --> phi +(K,K*)-bar, piK and elastic
       IF( ((lb1.eq.23.or.lb1.eq.21.or.iabs(lb1).eq.30) .and.
     &      (lb2.ge.25.and.lb2.le.28)) .OR.
     &     ((lb2.eq.23.or.lb2.eq.21.or.iabs(lb2).eq.30) .and.
     &      (lb1.ge.25.and.lb1.le.28)) ) go to 8799
c
c* K*(-bar) + pi --> phi + (K,K*)-bar
       IF( (iabs(lb1).eq.30.and.(lb2.ge.3.and.lb2.le.5)) .OR.
     &     (iabs(lb2).eq.30.and.(lb1.ge.3.and.lb1.le.5)) )go to 8799
*
c
c* phi + N --> pi+N(D),  rho+N(D),  K+ +La
c* phi + D --> pi+N(D),  rho+N(D)
       IF( (lb1.eq.29 .and.(lb2.eq.1.or.lb2.eq.2.or.
     &       (lb2.ge.6.and.lb2.le.9))) .OR.
     &     (lb2.eq.29 .and.(lb1.eq.1.or.lb1.eq.2.or.
     &       (lb1.ge.6.and.lb1.le.9))) )go to 7222
c
c* phi + (pi,rho,ome,K,K*-bar) --> K+K, K+K*, K*+K*, (pi,rho,omega)+(K,K*-bar)
       IF( (lb1.eq.29 .and.((lb2.ge.3.and.lb2.le.5).or.
     &      (lb2.ge.21.and.lb2.le.28).or.iabs(lb2).eq.30)) .OR.
     &     (lb2.eq.29 .and.((lb1.ge.3.and.lb1.le.5).or.
     &      (lb1.ge.21.and.lb1.le.28).or.iabs(lb1).eq.30)) )THEN
             go to 7444
      endif
*
c
* La/Si, Cas, Om (bar)-(rho,omega,phi) elastic colln
* pion vs. La, Ca, Omega-(bar) elastic coll. treated in resp. subroutines
      if( ((iabs(lb1).ge.14.and.iabs(lb1).le.17).or.iabs(lb1).ge.40)
     &    .and.((lb2.ge.25.and.lb2.le.29).or.lb2.eq.0) )go to 888
      if( ((iabs(lb2).ge.14.and.iabs(lb2).le.17).or.iabs(lb2).ge.40)
     &    .and.((lb1.ge.25.and.lb1.le.29).or.lb1.eq.0) )go to 888
c
c K+/K* (N,R)  OR   K-/K*- (N,R)-bar  elastic scatt
        if( ((lb1.eq.23.or.lb1.eq.30).and.(lb2.eq.1.or.lb2.eq.2.or.
     &         (lb2.ge.6.and.lb2.le.13))) .OR.
     &      ((lb2.eq.23.or.lb2.eq.30).and.(lb1.eq.1.or.lb1.eq.2.or.
     &         (lb1.ge.6.and.lb1.le.13))) ) go to 888
        if( ((lb1.eq.21.or.lb1.eq.-30).and.(lb2.eq.-1.or.lb2.eq.-2.or.
     &       (lb2.ge.-13.and.lb2.le.-6))) .OR. 
     &      ((lb2.eq.21.or.lb2.eq.-30).and.(lb1.eq.-1.or.lb1.eq.-2.or.
     &       (lb1.ge.-13.and.lb1.le.-6))) ) go to 888
c
* L/S-baryon elastic collision 
       If( ((lb1.ge.14.and.lb1.le.17).and.(lb2.ge.6.and.lb2.le.13))
     & .OR.((lb2.ge.14.and.lb2.le.17).and.(lb1.ge.6.and.lb1.le.13)) )
     &   go to 7799
       If(((lb1.le.-14.and.lb1.ge.-17).and.(lb2.le.-6.and.lb2.ge.-13))
     &.OR.((lb2.le.-14.and.lb2.ge.-17).and.(lb1.le.-6.and.lb1.ge.-13)))
     &   go to 7799
c
c skip other collns with perturbative particles or hyperon-bar
       if( iabs(lb1).ge.40 .or. iabs(lb2).ge.40
     &    .or. (lb1.le.-14.and.lb1.ge.-17) 
     &    .or. (lb2.le.-14.and.lb2.ge.-17) )go to 400
c
c
* anti-baryon on baryon resonaces 
        if((lb1.eq.-1.or.lb1.eq.-2.or.(lb1.ge.-13.and.lb1.le.-6))
     1 .and.(lb2.eq.1.or.lb2.eq.2.or.(lb2.ge.6.and.lb2.le.13))) then
            GOTO 2799
       else if((lb2.eq.-1.or.lb2.eq.-2.or.(lb2.ge.-13.and.lb2.le.-6))
     1 .and.(lb1.eq.1.or.lb1.eq.2.or.(lb1.ge.6.and.lb1.le.13))) then
            GOTO 2799
         END IF
c
clin-10/25/02 get rid of argument usage mismatch in newka():
         inewka=irun
c        call newka(icase,irun,iseed,dt,nt,
clin-5/01/03 set iblock value in art1f.f, necessary for resonance studies:
c        call newka(icase,inewka,iseed,dt,nt,
c     &                  ictrl,i1,i2,srt,pcx,pcy,pcz)
        call newka(icase,inewka,iseed,dt,nt,
     &                  ictrl,i1,i2,srt,pcx,pcy,pcz,iblock)
clin-10/25/02-end
        IF (ICTRL .EQ. 1) GOTO 400
c
* SEPARATE NUCLEON+NUCLEON( BARYON RESONANCE+ BARYON RESONANCE ELASTIC
* COLLISION), BARYON RESONANCE+NUCLEON AND BARYON-PION
* COLLISIONS INTO THREE PARTS TO CHECK IF THEY ARE GOING TO SCATTER,
* WE only allow L/S to COLLIDE elastically with a nucleon and meson
       if((iabs(lb1).ge.14.and.iabs(lb1).le.17).
     &  or.(iabs(lb2).ge.14.and.iabs(lb2).le.17))go to 400
* IF PION+PION COLLISIONS GO TO 777
* if pion+eta, eta+eta to create kaons go to 777 
       IF((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.3.and.lb2.le.5))GO TO 777
       if(lb1.eq.0.and.(lb2.ge.3.and.lb2.le.5)) go to 777
       if(lb2.eq.0.and.(lb1.ge.3.and.lb1.le.5)) go to 777
       if(lb1.eq.0.and.lb2.eq.0)go to 777
* we assume that rho and omega behave the same way as pions in 
* kaon production
* (1) rho(omega)+rho(omega)
       if( (lb1.ge.25.and.lb1.le.28).and.
     &     (lb2.ge.25.and.lb2.le.28) )goto 777
* (2) rho(omega)+pion
      If((lb1.ge.25.and.lb1.le.28).and.(lb2.ge.3.and.lb2.le.5))go to 777
      If((lb2.ge.25.and.lb2.le.28).and.(lb1.ge.3.and.lb1.le.5))go to 777
* (3) rho(omega)+eta
       if((lb1.ge.25.and.lb1.le.28).and.lb2.eq.0)go to 777
       if((lb2.ge.25.and.lb2.le.28).and.lb1.eq.0)go to 777
c
* if kaon+pion collisions go to 889
       if((lb1.eq.23.or.lb1.eq.21).and.(lb2.ge.3.and.lb2.le.5))go to 889
       if((lb2.eq.23.or.lb2.eq.21).and.(lb1.ge.3.and.lb1.le.5))go to 889
c
clin-2/06/03 skip all other (K K* Kbar K*bar) channels:
* SKIP all other K and K* RESCATTERINGS
        If(iabs(lb1).eq.30.or.iabs(lb2).eq.30) go to 400
        If(lb1.eq.21.or.lb2.eq.21) go to 400
        If(lb1.eq.23.or.lb2.eq.23) go to 400
c
* IF PION+baryon COLLISION GO TO 3
           IF( (LB1.ge.3.and.LB1.le.5) .and. 
     &         (iabs(LB2).eq.1.or.iabs(LB2).eq.2.or.
     &          (iabs(LB2).ge.6.and.iabs(LB2).le.13)) )GO TO 3
           IF( (LB2.ge.3.and.LB2.le.5) .and. 
     &         (iabs(LB1).eq.1.or.iabs(LB1).eq.2.or.
     &          (iabs(LB1).ge.6.and.iabs(LB1).le.13)) )GO TO 3
c
* IF rho(omega)+NUCLEON (baryon resonance) COLLISION GO TO 33
           IF( (LB1.ge.25.and.LB1.le.28) .and. 
     &         (iabs(LB2).eq.1.or.iabs(LB2).eq.2.or.
     &          (iabs(LB2).ge.6.and.iabs(LB2).le.13)) )GO TO 33
           IF( (LB2.ge.25.and.LB2.le.28) .and. 
     &         (iabs(LB1).eq.1.or.iabs(LB1).eq.2.or.
     &          (iabs(LB1).ge.6.and.iabs(LB1).le.13)) )GO TO 33
c
* IF ETA+NUCLEON (baryon resonance) COLLISIONS GO TO 547
           IF( LB1.eq.0 .and. 
     &         (iabs(LB2).eq.1.or.iabs(LB2).eq.2.or.
     &          (iabs(LB2).ge.6.and.iabs(LB2).le.13)) )GO TO 547
           IF( LB2.eq.0 .and. 
     &         (iabs(LB1).eq.1.or.iabs(LB1).eq.2.or.
     &          (iabs(LB1).ge.6.and.iabs(LB1).le.13)) )GO TO 547
c
* IF NUCLEON+BARYON RESONANCE COLLISION GO TO 44
            IF((LB1.eq.1.or.lb1.eq.2).
     &        AND.(LB2.GT.5.and.lb2.le.13))GOTO 44
            IF((LB2.eq.1.or.lb2.eq.2).
     &        AND.(LB1.GT.5.and.lb1.le.13))GOTO 44
            IF((LB1.eq.-1.or.lb1.eq.-2).
     &        AND.(LB2.LT.-5.and.lb2.ge.-13))GOTO 44
            IF((LB2.eq.-1.or.lb2.eq.-2).
     &        AND.(LB1.LT.-5.and.lb1.ge.-13))GOTO 44
c
* IF NUCLEON+NUCLEON COLLISION GO TO 4
       IF((LB1.eq.1.or.lb1.eq.2).AND.(LB2.eq.1.or.lb2.eq.2))GOTO 4
       IF((LB1.eq.-1.or.lb1.eq.-2).AND.(LB2.eq.-1.or.lb2.eq.-2))GOTO 4
c
* IF BARYON RESONANCE+BARYON RESONANCE COLLISION GO TO 444
            IF((LB1.GT.5.and.lb1.le.13).AND.
     &         (LB2.GT.5.and.lb2.le.13)) GOTO 444
            IF((LB1.LT.-5.and.lb1.ge.-13).AND.
     &         (LB2.LT.-5.and.lb2.ge.-13)) GOTO 444
c
* if L/S+L/S or L/s+nucleon go to 400
* otherwise, develop a model for their collisions
       if((lb1.lt.3).and.(lb2.ge.14.and.lb2.le.17))goto 400
       if((lb2.lt.3).and.(lb1.ge.14.and.lb1.le.17))goto 400
       if((lb1.ge.14.and.lb1.le.17).and.
     &  (lb2.ge.14.and.lb2.le.17))goto 400
c
* otherwise, go out of the loop
              go to 400
*
*
547           IF(LB1*LB2.EQ.0)THEN
* (1) FOR ETA+NUCLEON SYSTEM, we allow both elastic collision, 
*     i.e. N*(1535) formation and kaon production
*     the total kaon production cross section is
*     ASSUMED to be THE SAME AS PION+NUCLEON COLLISIONS
* (2) for eta+baryon resonance we only allow kaon production
           ece=(em1+em2+0.02)**2
           xkaon0=0.
           if(srt.ge.1.63.AND.SRT.LE.1.7)xkaon0=pnlka(srt)
           IF(SRT.GT.1.7)XKAON0=PNLKA(SRT)+pnska(srt)
cbz3/7/99 neutralk
            XKAON0 = 2.0 * XKAON0
cbz3/7/99 neutralk end
* Here we negelect eta+n inelastic collisions other than the 
* kaon production, therefore the total inelastic cross section
* xkaon equals to the xkaon0 (kaon production cross section)
           xkaon=xkaon0
* note here the xkaon is in unit of fm**2
            XETA=XN1535(I1,I2,0)
        If((iabs(LB(I1)).ge.6.and.iabs(LB(I1)).le.13).or.
     &     (iabs(LB(I2)).ge.6.and.iabs(LB(I2)).le.13)) xeta=0.      
            IF((XETA+xkaon).LE.1.e-06)GO TO 400
            DSE=SQRT((XETA+XKAON)/PI)
           DELTRE=DSE+0.1
        px1cm=pcx
        py1cm=pcy
        pz1cm=pcz
* CHECK IF N*(1535) resonance CAN BE FORMED
         CALL DISTCE(I1,I2,DELTRE,DSE,DT,ECE,SRT,IC,
     1   PCX,PCY,PCZ)
         IF(IC.EQ.-1) GO TO 400
         ekaon(4,iss)=ekaon(4,iss)+1
        IF(XKAON0/(XKAON+XETA).GT.RANART(NSEED))then
* kaon production, USE CREN TO CALCULATE THE MOMENTUM OF L/S K+
        CALL CREN(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)
* kaon production
       IF(IBLOCK.EQ.7) then
          LPN=LPN+1
       elseIF(IBLOCK.EQ.-7) then
       endif
c
       em1=e(i1)
       em2=e(i2)
       GO TO 440
       endif
* N*(1535) FORMATION
        resona=1.
         GO TO 98
         ENDIF
*IF PION+NUCLEON (baryon resonance) COLLISION THEN
3           CONTINUE
           px1cm=pcx
           py1cm=pcy
           pz1cm=pcz
* the total kaon production cross section for pion+baryon (resonance) is
* assumed to be the same as in pion+nucleon
           xkaon0=0.
           if(srt.ge.1.63.AND.SRT.LE.1.7)xkaon0=pnlka(srt)
           IF(SRT.GT.1.7)XKAON0=PNLKA(SRT)+pnska(srt)
            XKAON0 = 2.0 * XKAON0
c
c sp11/21/01  phi production: pi +N(D) -> phi + N(D)
         Xphi = 0.
       if( ( ((lb1.ge.1.and.lb1.le.2).or.
     &        (lb1.ge.6.and.lb1.le.9))
     &   .OR.((lb2.ge.1.and.lb2.le.2).or.
     &        (lb2.ge.6.and.lb2.le.9)) )
     &       .AND. srt.gt.1.958)
     &        call pibphi(srt,lb1,lb2,em1,em2,Xphi,xphin)
c !! in fm^2 above
* if a pion collide with a baryon resonance, 
* we only allow kaon production AND the reabsorption 
* processes: Delta+pion-->N+pion, N*+pion-->N+pion
* Later put in pion+baryon resonance elastic
* cross through forming higher resonances implicitly.
c          If(em1.gt.1.or.em2.gt.1.)go to 31
         If((iabs(LB(I1)).ge.6.and.iabs(LB(I1)).le.13).or.
     &      (iabs(LB(I2)).ge.6.and.iabs(LB(I2)).le.13)) go to 31
* For pion+nucleon collisions: 
* using the experimental pion+nucleon inelastic cross section, we assume it
* is exhausted by the Delta+pion, Delta+rho and Delta+omega production 
* and kaon production. In the following we first check whether 
* inelastic pion+n collision can happen or not, then determine in 
* crpn whether it is through pion production or through kaon production
* note that the xkaon0 is the kaon production cross section
* Note in particular that: 
* xkaon in the following is the total pion+nucleon inelastic cross section
* note here the xkaon is in unit of fm**2, xnpi is also in unit of fm**2
* FOR PION+NUCLEON SYSTEM, THE MINIMUM S IS 1.2056 the minimum srt for 
* elastic scattering, and it is 1.60 for pion production, 1.63 for LAMBDA+kaon 
* production and 1.7 FOR SIGMA+KAON
* (EC = PION MASS+NUCLEON MASS+20MEV)**2
            EC=(em1+em2+0.02)**2
           xkaon=0.
           if(srt.gt.1.23)xkaon=(pionpp(srt)+PIPP1(SRT))/2.
* pion+nucleon elastic cross section is divided into two parts:
* (1) forming D(1232)+N*(1440) +N*(1535)
* (2) cross sections forming higher resonances are calculated as
*     the difference between the total elastic and (1), this part is 
*     treated as direct process since we do not explicitLY include
*     higher resonances.
* the following is the resonance formation cross sections.
*1. PION(+)+PROTON-->DELTA++,PION(-)+NEUTRON-->DELTA(-)
           IF( (LB1*LB2.EQ.5.OR.((LB1*LB2.EQ.6).AND.
     &         (LB1.EQ.3.OR.LB2.EQ.3)))
     &    .OR. (LB1*LB2.EQ.-3.OR.((LB1*LB2.EQ.-10).AND.
     &         (LB1.EQ.5.OR.LB2.EQ.5))) )then    
              XMAX=190.
              xmaxn=0
              xmaxn1=0
              xdirct=dirct1(srt)
               go to 678
           endif
*2. PION(-)+PROTON-->DELTA0,PION(+)+NEUTRON-->DELTA+ 
*   or N*(+)(1440) or N*(+)(1535)
* note the factor 2/3 is from the isospin consideration and
* the factor 0.6 or 0.5 is the branching ratio for the resonance to decay
* into pion+nucleon
            IF( (LB1*LB2.EQ.3.OR.((LB1*LB2.EQ.10).AND.
     &          (LB1.EQ.5.OR.LB2.EQ.5)))
     &     .OR. (LB1*LB2.EQ.-5.OR.((LB1*LB2.EQ.-6).AND.
     &          (LB1.EQ.3.OR.LB2.EQ.3))) )then      
              XMAX=27.
              xmaxn=2./3.*25.*0.6
               xmaxn1=2./3.*40.*0.5
              xdirct=dirct2(srt)
               go to 678
              endif
*3. PION0+PROTON-->DELTA+,PION0+NEUTRON-->DELTA0, or N*(0)(1440) or N*(0)(1535)
            IF((LB1.EQ.4.OR.LB2.EQ.4).AND.
     &         (iabs(LB1*LB2).EQ.4.OR.iabs(LB1*LB2).EQ.8))then
              XMAX=50.
              xmaxn=1./3.*25*0.6
              xmaxn1=1/3.*40.*0.5
              xdirct=dirct3(srt)
                go to 678
              endif
678           xnpin1=0
           xnpin=0
            XNPID=XNPI(I1,I2,1,XMAX)
           if(xmaxn1.ne.0)xnpin1=XNPI(i1,i2,2,XMAXN1)
            if(xmaxn.ne.0)XNPIN=XNPI(I1,I2,0,XMAXN)
* the following 
           xres=xnpid+xnpin+xnpin1
           xnelas=xres+xdirct 
           icheck=1
           go to 34
* For pion + baryon resonance the reabsorption 
* cross section is calculated from the detailed balance
* using reab(i1,i2,srt,ictrl), ictrl=1, 2 and 3
* for pion, rho and omega + baryon resonance
31           ec=(em1+em2+0.02)**2
           xreab=reab(i1,i2,srt,1)
clin-12/02/00 to satisfy detailed balance, forbid N* absorptions:
          if((iabs(lb1).ge.10.and.iabs(lb1).le.13)
     1         .or.(iabs(lb2).ge.10.and.iabs(lb2).le.13)) XREAB = 0.
           xkaon=xkaon0+xreab
* a constant of 10 mb IS USED FOR PION + N* RESONANCE, 
        IF((iabs(LB1).GT.9.AND.iabs(LB1).LE.13) .OR.
     &      (iabs(LB2).GT.9.AND.iabs(LB2).LE.13))THEN
           Xnelas=1.0
        ELSE
           XNELAS=DPION(EM1,EM2,LB1,LB2,SRT)
        ENDIF
           icheck=2
34          IF((Xnelas+xkaon+Xphi).LE.0.000001)GO TO 400
            DS=SQRT((Xnelas+xkaon+Xphi)/PI)
csp09/20/01
c           totcr = xnelas+xkaon
c           if(srt .gt. 3.5)totcr = max1(totcr,3.)
c           DS=SQRT(totcr/PI)
csp09/20/01 end
           deltar=ds+0.1
         CALL DISTCE(I1,I2,DELTAR,DS,DT,EC,SRT,IC,
     1   PCX,PCY,PCZ)
         IF(IC.EQ.-1) GO TO 400
       ekaon(4,iss)=ekaon(4,iss)+1
c***
* check what kind of collision has happened
* (1) pion+baryon resonance
* if direct elastic process
        if(icheck.eq.2)then
c  !!sp11/21/01
      if(xnelas/(xnelas+xkaon+Xphi).ge.RANART(NSEED))then
c               call Crdir(PX1CM,PY1CM,PZ1CM,SRT,I1,I2)
               call Crdir(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)
              go to 440
              else
* for inelastic process, go to 96 to check
* kaon production and pion reabsorption : pion+D(N*)-->pion+N
               go to 96
                endif
              endif
*(2) pion+n
* CHECK IF inELASTIC COLLISION IS POSSIBLE FOR PION+N COLLISIONS
clin-8/17/00 typo corrected, many other occurences:
c        IF(XKAON/(XKAON+Xnelas).GT.RANART(NSEED))GO TO 95
       IF((XKAON+Xphi)/(XKAON+Xphi+Xnelas).GT.RANART(NSEED))GO TO 95
* direct process
        if(xdirct/xnelas.ge.RANART(NSEED))then
c               call Crdir(PX1CM,PY1CM,PZ1CM,SRT,I1,I2)
               call Crdir(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)
              go to 440
              endif
* now resonance formation or direct process (higher resonances)
           IF( (LB1*LB2.EQ.5.OR.((LB1*LB2.EQ.6).AND.
     &         (LB1.EQ.3.OR.LB2.EQ.3)))
     &    .OR. (LB1*LB2.EQ.-3.OR.((LB1*LB2.EQ.-10).AND.
     &         (LB1.EQ.5.OR.LB2.EQ.5))) )then    
c
* ONLY DELTA RESONANCE IS POSSIBLE, go to 99
        GO TO 99
       else
* NOW BOTH DELTA AND N* RESORANCE ARE POSSIBLE
* DETERMINE THE RESORANT STATE BY USING THE MONTRE CARLO METHOD
            XX=(XNPIN+xnpin1)/xres
            IF(RANART(NSEED).LT.XX)THEN
* N* RESONANCE IS SELECTED
* decide N*(1440) or N*(1535) formation
        xx0=xnpin/(xnpin+xnpin1)
        if(RANART(NSEED).lt.xx0)then
         RESONA=0.
* N*(1440) formation
         GO TO 97
        else
* N*(1535) formation
        resona=1.
         GO TO 98
        endif
         ELSE
* DELTA RESONANCE IS SELECTED
         GO TO 99
         ENDIF
         ENDIF
97       CONTINUE
            IF(RESONA.EQ.0.)THEN
*N*(1440) IS PRODUCED,WE DETERMINE THE CHARGE STATE OF THE PRODUCED N*
            I=I1
            IF(EM1.LT.0.6)I=I2
* (0.1) n+pion(+)-->N*(+)
           IF( (LB1*LB2.EQ.10.AND.(LB1.EQ.5.OR.LB2.EQ.5))
     &      .OR.(LB1*LB2.EQ.-6.AND.(LB1.EQ.3.OR.LB2.EQ.3)) )THEN
            LB(I)=11
           go to 303
            ENDIF
* (0.2) p+pion(0)-->N*(+)
c            IF(LB(I1)*LB(I2).EQ.4.AND.(LB(I1).EQ.1.OR.LB(I2).EQ.1))THEN
            IF(iabs(LB(I1)*LB(I2)).EQ.4.AND.
     &         (LB(I1).EQ.4.OR.LB(I2).EQ.4))THEN    
            LB(I)=11
           go to 303
            ENDIF
* (0.3) n+pion(0)-->N*(0)
c            IF(LB(I1)*LB(I2).EQ.8.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN
            IF(iabs(LB(I1)*LB(I2)).EQ.8.AND.
     &        (LB(I1).EQ.4.OR.LB(I2).EQ.4))THEN    
            LB(I)=10
           go to 303
            ENDIF
* (0.4) p+pion(-)-->N*(0)
c            IF(LB(I1)*LB(I2).EQ.3)THEN
            IF( (LB(I1)*LB(I2).EQ.3)
     &      .OR.(LB(I1)*LB(I2).EQ.-5) )THEN
            LB(I)=10
            ENDIF
303         CALL DRESON(I1,I2)
            if(LB1.lt.0.or.LB2.lt.0) LB(I)=-LB(I)
            lres=lres+1
            GO TO 101
*COM: GO TO 101 TO CHANGE THE PHASE SPACE DENSITY OF THE NUCLEON
            ENDIF
98          IF(RESONA.EQ.1.)THEN
*N*(1535) IS PRODUCED, WE DETERMINE THE CHARGE STATE OF THE PRODUCED N*
            I=I1
            IF(EM1.LT.0.6)I=I2
* note: this condition applies to both eta and pion
* (0.1) n+pion(+)-->N*(+)
c            IF(LB1*LB2.EQ.10.AND.(LB1.EQ.2.OR.LB2.EQ.2))THEN
            IF( (LB1*LB2.EQ.10.AND.(LB1.EQ.5.OR.LB2.EQ.5))
     &      .OR.(LB1*LB2.EQ.-6.AND.(LB1.EQ.3.OR.LB2.EQ.3)) )THEN
            LB(I)=13
           go to 304
            ENDIF
* (0.2) p+pion(0)-->N*(+)
c            IF(LB(I1)*LB(I2).EQ.4.AND.(LB(I1).EQ.1.OR.LB(I2).EQ.1))THEN
            IF(iabs(LB(I1)*LB(I2)).EQ.4.AND.
     &           (LB(I1).EQ.4.OR.LB(I2).EQ.4))THEN 
            LB(I)=13
           go to 304
            ENDIF
* (0.3) n+pion(0)-->N*(0)
c            IF(LB(I1)*LB(I2).EQ.8.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN
            IF(iabs(LB(I1)*LB(I2)).EQ.8.AND.
     &           (LB(I1).EQ.4.OR.LB(I2).EQ.4))THEN      
            LB(I)=12
           go to 304
            ENDIF
* (0.4) p+pion(-)-->N*(0)
c            IF(LB(I1)*LB(I2).EQ.3)THEN
            IF( (LB(I1)*LB(I2).EQ.3)
     &      .OR.(LB(I1)*LB(I2).EQ.-5) )THEN
            LB(I)=12
           go to 304
           endif
* (0.5) p+eta-->N*(+)(1535),n+eta-->N*(0)(1535)
           if(lb(i1)*lb(i2).eq.0)then
c            if((lb(i1).eq.1).or.(lb(i2).eq.1))then
            if(iabs(lb(i1)).eq.1.or.iabs(lb(i2)).eq.1)then
           LB(I)=13
           go to 304
           ELSE
           LB(I)=12
           ENDIF
           endif
304         CALL DRESON(I1,I2)
            if(LB1.lt.0.or.LB2.lt.0) LB(I)=-LB(I) 
            lres=lres+1
            GO TO 101
*COM: GO TO 101 TO CHANGE THE PHASE SPACE DENSITY OF THE NUCLEON
            ENDIF
*DELTA IS PRODUCED,IN THE FOLLOWING WE DETERMINE THE
*CHARGE STATE OF THE PRODUCED DELTA
99      LRES=LRES+1
        I=I1
        IF(EM1.LE.0.6)I=I2
* (1) p+pion(+)-->DELTA(++)
c        IF(LB(I1)*LB(I2).EQ.5)THEN
            IF( (LB(I1)*LB(I2).EQ.5)
     &      .OR.(LB(I1)*LB(I2).EQ.-3) )THEN
        LB(I)=9
       go to 305
        ENDIF
* (2) p+pion(0)-->delta(+)
c        IF(LB(I1)*LB(I2).EQ.4.AND.(LB(I1).EQ.1.OR.LB(I2).EQ.1))then
       IF(iabs(LB(I1)*LB(I2)).EQ.4.AND.(LB(I1).EQ.4.OR.LB(I2).EQ.4))then
        LB(I)=8
       go to 305
        ENDIF
* (3) n+pion(+)-->delta(+)
c        IF(LB(I1)*LB(I2).EQ.10.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN
       IF( (LB(I1)*LB(I2).EQ.10.AND.(LB(I1).EQ.5.OR.LB(I2).EQ.5))
     & .OR.(LB(I1)*LB(I2).EQ.-6.AND.(LB(I1).EQ.3.OR.LB(I2).EQ.3)) )THEN
        LB(I)=8
       go to 305
        ENDIF
* (4) n+pion(0)-->delta(0)
c        IF(LB(I1)*LB(I2).EQ.8.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN
       IF(iabs(LB(I1)*LB(I2)).EQ.8.AND.(LB(I1).EQ.4.OR.LB(I2).EQ.4))THEN
        LB(I)=7
       go to 305
        ENDIF
* (5) p+pion(-)-->delta(0)
c        IF(LB(I1)*LB(I2).EQ.3)THEN
            IF( (LB(I1)*LB(I2).EQ.3)
     &      .OR.(LB(I1)*LB(I2).EQ.-5) )THEN
        LB(I)=7
       go to 305
        ENDIF
* (6) n+pion(-)-->delta(-)
c        IF(LB(I1)*LB(I2).EQ.6.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN
       IF( (LB(I1)*LB(I2).EQ.6.AND.(LB(I1).EQ.3.OR.LB(I2).EQ.3))
     & .OR.(LB(I1)*LB(I2).EQ.-10.AND.(LB(I1).EQ.5.OR.LB(I2).EQ.5)) )THEN 
        LB(I)=6
        ENDIF
305     CALL DRESON(I1,I2)
        if(LB1.lt.0.or.LB2.lt.0) LB(I)=-LB(I) 
       GO TO 101
csp-11/08/01 K*
* FOR kaON+pion COLLISIONS, form K* (bar) or
c La/Si-bar + N <-- pi + K+
c La/Si + N-bar <-- pi + K-                                             
c phi + K <-- pi + K                                             
clin (rho,omega) + K* <-- pi + K
889       CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
* the cross section is from C.M. Ko, PRC 23, 2760 (1981).
       spika=60./(1.+4.*(srt-0.895)**2/(0.05)**2)
c
cc       if(lb(i1).eq.23.or.lb(i2).eq.23)then   !! block  K- + pi->La + B-bar
        call Crkpla(PX1CM,PY1CM,PZ1CM,EC,SRT,spika,
     &                  emm1,emm2,lbp1,lbp2,I1,I2,icase,srhoks)
cc
c* only K* or K*bar formation
c       else 
c      DSkn=SQRT(spika/PI/10.)
c      dsknr=dskn+0.1
c      CALL DISTCE(I1,I2,dsknr,DSkn,DT,EC,SRT,IC,
c    1     PX1CM,PY1CM,PZ1CM)
c        IF(IC.EQ.-1) GO TO 400
c       icase = 1
c      endif
c
         if(icase .eq. 0) then
            iblock=0
            go to 400
         endif
       if(icase .eq. 1)then
             call KSRESO(I1,I2)
clin-4/30/03 give non-zero iblock for resonance selections:
             iblock = 171
ctest off for resonance (phi, K*) studies:
c             if(iabs(lb(i1)).eq.30) then
c             write(17,112) 'ks',lb(i1),p(1,i1),p(2,i1),p(3,i1),e(i1),nt
c             elseif(iabs(lb(i2)).eq.30) then
c             write(17,112) 'ks',lb(i2),p(1,i2),p(2,i2),p(3,i2),e(i2),nt
c             endif
c
              lres=lres+1
              go to 101
       elseif(icase .eq. 2)then
             iblock = 71
c
* La/Si (bar) formation                                                   
       elseif(iabs(icase).eq.5)then
             iblock = 88
       else
*
* phi formation
             iblock = 222
       endif
             LB(I1) = lbp1
             LB(I2) = lbp2
             E(I1) = emm1
             E(I2) = emm2
             em1=e(i1)
             em2=e(i2)
             ntag = 0
             go to 440
c             
33       continue
       em1=e(i1)
       em2=e(i2)
* (1) if rho or omega collide with a nucleon we allow both elastic 
*     scattering and kaon production to happen if collision conditions 
*     are satisfied.
* (2) if rho or omega collide with a baryon resonance we allow
*     kaon production, pion reabsorption: rho(omega)+D(N*)-->pion+N
*     and NO elastic scattering to happen
           xelstc=0
            if((lb1.ge.25.and.lb1.le.28).and.
     &    (iabs(lb2).eq.1.or.iabs(lb2).eq.2))
     &      xelstc=ERHON(EM1,EM2,LB1,LB2,SRT)
            if((lb2.ge.25.and.lb2.le.28).and.
     &   (iabs(lb1).eq.1.or.iabs(lb1).eq.2))
     &      xelstc=ERHON(EM1,EM2,LB1,LB2,SRT)
            ec=(em1+em2+0.02)**2
* the kaon production cross section is
           xkaon0=0
           if(srt.ge.1.63.AND.SRT.LE.1.7)xkaon0=pnlka(srt)
           IF(SRT.GT.1.7)XKAON0=PNLKA(SRT)+pnska(srt)
           if(xkaon0.lt.0)xkaon0=0
cbz3/7/99 neutralk
            XKAON0 = 2.0 * XKAON0
cbz3/7/99 neutralk end
* the total inelastic cross section for rho(omega)+N is
           xkaon=xkaon0
           ichann=0
* the total inelastic cross section for rho (omega)+D(N*) is 
* xkaon=xkaon0+reab(**) 
c sp11/21/01  phi production: rho + N(D) -> phi + N(D)
         Xphi = 0.
       if( ( (((lb1.ge.1.and.lb1.le.2).or.
     &         (lb1.ge.6.and.lb1.le.9))
     &         .and.(lb2.ge.25.and.lb2.le.27))
     &   .OR.(((lb2.ge.1.and.lb2.le.2).or.
     &         (lb2.ge.6.and.lb2.le.9))
     &        .and.(lb1.ge.25.and.lb1.le.27)) ).AND. srt.gt.1.958)
     &    call pibphi(srt,lb1,lb2,em1,em2,Xphi,xphin)
c !! in fm^2 above
c
        if((iabs(lb1).ge.6.and.lb2.ge.25).or.
     &    (lb1.ge.25.and.iabs(lb2).ge.6))then
           ichann=1
           ictrl=2
           if(lb1.eq.28.or.lb2.eq.28)ictrl=3
            xreab=reab(i1,i2,srt,ictrl)
clin-12/02/00 to satisfy detailed balance, forbid N* absorptions:
            if((iabs(lb1).ge.10.and.iabs(lb1).le.13)
     1           .or.(iabs(lb2).ge.10.and.iabs(lb2).le.13)) XREAB = 0.
        if(xreab.lt.0)xreab=1.E-06
            xkaon=xkaon0+xreab
          XELSTC=1.0
           endif
            DS=SQRT((XKAON+Xphi+xelstc)/PI)
c
csp09/20/01
c           totcr = xelstc+xkaon
c           if(srt .gt. 3.5)totcr = max1(totcr,3.)
c           DS=SQRT(totcr/PI)
csp09/20/01 end
c
        DELTAR=DS+0.1
       px1cm=pcx
       py1cm=pcy
       pz1cm=pcz
* CHECK IF the collision can happen
         CALL DISTCE(I1,I2,DELTAR,DS,DT,EC,SRT,IC,
     1   PCX,PCY,PCZ)
         IF(IC.EQ.-1) GO TO 400
        ekaon(4,iss)=ekaon(4,iss)+1
c*
* NOW rho(omega)+N or D(N*) COLLISION IS POSSIBLE
* (1) check elastic collision
       if(xelstc/(xelstc+xkaon+Xphi).gt.RANART(NSEED))then
c       call crdir(px1CM,py1CM,pz1CM,srt,I1,i2)
       call crdir(px1CM,py1CM,pz1CM,srt,I1,i2,IBLOCK)
       go to 440
       endif
* (2) check pion absorption or kaon production
        CALL CRRD(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK,xkaon0,xkaon,Xphi,xphin)
* kaon production
csp05/16/01
       IF(IBLOCK.EQ.7) then
          LPN=LPN+1
       elseIF(IBLOCK.EQ.-7) then
       endif
csp05/16/01 end
* rho obsorption
       if(iblock.eq.81) lrhor=lrhor+1
* omega obsorption
       if(iblock.eq.82) lomgar=lomgar+1
       em1=e(i1)
       em2=e(i2)
       GO TO 440
* for pion+n now using the subroutine crpn to change 
* the particle label and set the new momentum of L/S+K final state
95       continue
* NOW PION+N INELASTIC COLLISION IS POSSIBLE
* check pion production or kaon production
        CALL CRPN(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK,xkaon0,xkaon,Xphi,xphin)
* kaon production
csp05/16/01
       IF(IBLOCK.EQ.7) then
          LPN=LPN+1
       elseIF(IBLOCK.EQ.-7) then
       endif
csp05/16/01 end
* pion production
       if(iblock.eq.77) lpd=lpd+1
* rho production
       if(iblock.eq.78) lrho=lrho+1
* omega production
       if(iblock.eq.79) lomega=lomega+1
       em1=e(i1)
       em2=e(i2)
       GO TO 440
* for pion+D(N*) now using the subroutine crpd to 
* (1) check kaon production or pion reabsorption 
* (2) change the particle label and set the new 
*     momentum of L/S+K final state
96       continue
        CALL CRPD(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK,xkaon0,xkaon,Xphi,xphin)
* kaon production
csp05/16/01
       IF(IBLOCK.EQ.7) then
          LPN=LPN+1
       elseIF(IBLOCK.EQ.-7) then
       endif
csp05/16/01 end
* pion obserption
       if(iblock.eq.80) lpdr=lpdr+1
       em1=e(i1)
       em2=e(i2)
       GO TO 440
* CALCULATE KAON PRODUCTION PROBABILITY FROM PION + N COLLISIONS
C        IF(SRT.GT.1.615)THEN
C        CALL PKAON(SRT,XXp,PK)
C        TKAON(7)=TKAON(7)+PK 
C        EKAON(7,ISS)=EKAON(7,ISS)+1
c        CALL KSPEC1(SRT,PK)
C        call LK(3,srt,iseed,pk)
C        ENDIF
* negelecting the pauli blocking at high energies
101       continue
        IF(E(I2).EQ.0.)GO TO 600
        IF(E(I1).EQ.0.)GO TO 800
* IF NUCLEON+BARYON RESONANCE COLLISIONS
44      CONTINUE
* CALCULATE THE TOTAL CROSS SECTION OF NUCLEON+ BARYON RESONANCE COLLISION
* WE ASSUME THAT THE ELASTIC CROSS SECTION IS THE SAME AS NUCLEON+NUCLEON
* COM: WE USE THE PARAMETERISATION BY CUGNON FOR LOW ENERGIES
*      AND THE PARAMETERIZATIONS FROM CERN DATA BOOK FOR HIGHER 
*      ENERGIES. THE CUTOFF FOR THE TOTAL CROSS SECTION IS 55 MB 
       cutoff=em1+em2+0.02
       IF(SRT.LE.CUTOFF)GO TO 400
        IF(SRT.GT.2.245)THEN
       SIGNN=PP2(SRT)
       ELSE
        SIGNN = 35.0 / (1. + (SRT - CUTOFF) * 100.0)  +  20.0
       ENDIF 
        call XND(pcx,pcy,pcz,srt,I1,I2,xinel,
     &               sigk,xsk1,xsk2,xsk3,xsk4,xsk5)
       sig=signn+xinel
* For nucleon+baryon resonance collision, the minimum cms**2 energy is
        EC=(EM1+EM2+0.02)**2
* CHECK THE DISTENCE BETWEEN THE TWO PARTICLES
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
clin-6/2008 Deuteron production:
        ianti=0
        if(lb(i1).lt.0 .and. lb(i2).lt.0) ianti=1
        call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
        sig=sig+sdprod
clin-6/2008 perturbative treatment of deuterons:
        ipdflag=0
        if(idpert.eq.1) then
           ipert1=1
           sigr0=sig
           dspert=sqrt(sigr0/pi/10.)
           dsrpert=dspert+0.1
           CALL DISTCE(I1,I2,dsrpert,dspert,DT,EC,SRT,IC,
     1          PX1CM,PY1CM,PZ1CM)
           IF(IC.EQ.-1) GO TO 363
           signn0=0.
           CALL CRND(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     &  IBLOCK,SIGNN0,SIGr0,sigk,xsk1,xsk2,xsk3,xsk4,xsk5,NT,ipert1)
c     &  IBLOCK,SIGNN,SIG,sigk,xsk1,xsk2,xsk3,xsk4,xsk5)
           ipdflag=1
 363       continue
           ipert1=0
        endif
        if(idpert.eq.2) ipert1=1
c
        DS=SQRT(SIG/(10.*PI))
        DELTAR=DS+0.1
        CALL DISTCE(I1,I2,DELTAR,DS,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
c        IF(IC.EQ.-1)GO TO 400
        IF(IC.EQ.-1) then
           if(ipdflag.eq.1) iblock=501
           GO TO 400
        endif
        ekaon(3,iss)=ekaon(3,iss)+1
* CALCULATE KAON PRODUCTION PROBABILITY FROM NUCLEON + BARYON RESONANCE 
* COLLISIONS
        go to 361
* CHECK WHAT KIND OF COLLISION HAS HAPPENED
 361    continue 
        CALL CRND(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     &     IBLOCK,SIGNN,SIG,sigk,xsk1,xsk2,xsk3,xsk4,xsk5,NT,ipert1)
c     &  IBLOCK,SIGNN,SIG,sigk,xsk1,xsk2,xsk3,xsk4,xsk5)
        IF(iblock.eq.0.and.ipdflag.eq.1) iblock=501
        IF(IBLOCK.EQ.11)THEN
           LNDK=LNDK+1
           GO TO 400
c        elseIF(IBLOCK.EQ.-11) then
        elseIF(IBLOCK.EQ.-11.or.iblock.eq.501) then
           GO TO 400
        ENDIF
        if(iblock .eq. 222)then
c    !! sp12/17/01 
           GO TO 400
        ENDIF
        em1=e(i1)
        em2=e(i2)
        GO TO 440
* IF NUCLEON+NUCLEON OR BARYON RESONANCE+BARYON RESONANCE COLLISIONS
4       CONTINUE
* PREPARE THE EALSTIC CROSS SECTION FOR BARYON+BARYON COLLISIONS
* COM: WE USE THE PARAMETERISATION BY CUGNON FOR SRT LEQ 2.0 GEV
*      AND THE PARAMETERIZATIONS FROM CERN DATA BOOK FOR HIGHER 
*      ENERGIES. THE CUTOFF FOR THE TOTAL CROSS SECTION IS 55 MB 
*      WITH LOW-ENERGY-CUTOFF
        CUTOFF=em1+em2+0.14
* AT HIGH ENERGIES THE ISOSPIN DEPENDENCE IS NEGLIGIBLE
* THE TOTAL CROSS SECTION IS TAKEN AS THAT OF THE PP 
* ABOVE E_KIN=800 MEV, WE USE THE ISOSPIN INDEPENDNET XSECTION
        IF(SRT.GT.2.245)THEN
           SIG=ppt(srt)
           SIGNN=SIG-PP1(SRT)
        ELSE
* AT LOW ENERGIES THE ISOSPIN DEPENDENCE FOR NN COLLISION IS STRONG
           SIG=XPP(SRT)
           IF(ZET(LB(I1))*ZET(LB(I2)).LE.0)SIG=XNP(SRT)
           IF(ZET(LB(I1))*ZET(LB(I2)).GT.0)SIG=XPP(SRT)
           IF(ZET(LB(I1)).EQ.0.
     &          AND.ZET(LB(I2)).EQ.0)SIG=XPP(SRT)
           if((lb(i1).eq.-1.and.lb(i2).eq.-2) .or.
     &          (lb(i2).eq.-1.and.lb(i1).eq.-2))sig=xnp(srt)
*     WITH LOW-ENERGY-CUTOFF
           IF (SRT .LT. 1.897) THEN
              SIGNN = SIG
           ELSE 
              SIGNN = 35.0 / (1. + (SRT - 1.897) * 100.0)  +  20.0
           ENDIF
        ENDIF 
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
clin-5/2008 Deuteron production cross sections were not included 
c     in the previous parameterized inelastic cross section of NN collisions  
c     (SIGinel=SIG-SIGNN), so they are added here:
        ianti=0
        if(lb(i1).lt.0 .and. lb(i2).lt.0) ianti=1
        call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
        sig=sig+sdprod
c
clin-5/2008 perturbative treatment of deuterons:
        ipdflag=0
        if(idpert.eq.1) then
c     For idpert=1: ipert1=1 means we will first treat deuteron perturbatively,
c     then we set ipert1=0 to treat regular NN or NbarNbar collisions including
c     the regular deuteron productions.
c     ipdflag=1 means perturbative deuterons are produced here:
           ipert1=1
           EC=2.012**2
c     Use the same cross section for NN/NNBAR collisions 
c     to trigger perturbative production
           sigr0=sig
c     One can also trigger with X*sbbdm() so the weight will not be too small;
c     but make sure to limit the maximum trigger Xsec:
c           sigr0=sdprod*25.
c           if(sigr0.ge.100.) sigr0=100.
           dspert=sqrt(sigr0/pi/10.)
           dsrpert=dspert+0.1
           CALL DISTCE(I1,I2,dsrpert,dspert,DT,EC,SRT,IC,
     1          PX1CM,PY1CM,PZ1CM)
           IF(IC.EQ.-1) GO TO 365
           signn0=0.
           CALL CRNN(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK,
     1          NTAG,signn0,sigr0,NT,ipert1)
           ipdflag=1
 365       continue
           ipert1=0
        endif
        if(idpert.eq.2) ipert1=1
c
clin-5/2008 in case perturbative deuterons are produced for idpert=1:
c        IF(SIGNN.LE.0)GO TO 400
        IF(SIGNN.LE.0) then
           if(ipdflag.eq.1) iblock=501
           GO TO 400
        endif
c
        EC=3.59709
        ds=sqrt(sig/pi/10.)
        dsr=ds+0.1
        IF((E(I1).GE.1.).AND.(e(I2).GE.1.))EC=4.75
        CALL DISTCE(I1,I2,dsr,ds,DT,EC,SRT,IC,
     1       PX1CM,PY1CM,PZ1CM)
clin-5/2008 in case perturbative deuterons are produced above:
c        IF(IC.EQ.-1) GO TO 400
        IF(IC.EQ.-1) then
           if(ipdflag.eq.1) iblock=501
           GO TO 400
        endif
c
* CALCULATE KAON PRODUCTION PROBABILITY FROM NUCLEON+NUCLEON OR 
* RESONANCE+RESONANCE COLLISIONS
        go to 362
C CHECK WHAT KIND OF COLLISION HAS HAPPENED 
 362    ekaon(1,iss)=ekaon(1,iss)+1
        CALL CRNN(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK,
     1       NTAG,SIGNN,SIG,NT,ipert1)
clin-5/2008 give iblock # in case pert deuterons are produced for idpert=1:
        IF(iblock.eq.0.and.ipdflag.eq.1) iblock=501
clin-5/2008 add iblock # for deuteron formation:
c        IF(IBLOCK.EQ.4.OR.IBLOCK.Eq.9.or.iblock.ge.44.OR.IBLOCK.EQ.-9
c     &       .or.iblock.eq.222)THEN
        IF(IBLOCK.EQ.4.OR.IBLOCK.Eq.9.or.iblock.ge.44.OR.IBLOCK.EQ.-9
     &       .or.iblock.eq.222.or.iblock.eq.501)THEN
c
c     !! sp12/17/01 above
* momentum of the three particles in the final state have been calculated
* in the crnn, go out of the loop
           LCOLL=LCOLL+1
           if(iblock.eq.4)then
              LDIRT=LDIRT+1
           elseif(iblock.eq.44)then
              LDdrho=LDdrho+1
           elseif(iblock.eq.45)then
              Lnnrho=Lnnrho+1
           elseif(iblock.eq.46)then
              Lnnom=Lnnom+1
           elseif(iblock .eq. 222)then
           elseIF(IBLOCK.EQ.9) then
              LNNK=LNNK+1
           elseIF(IBLOCK.EQ.-9) then
           endif
           GO TO 400
        ENDIF
        em1=e(i1)
        em2=e(i2)
        GO TO 440
clin-8/2008 B+B->Deuteron+Meson over
c
clin-8/2008 Deuteron+Meson->B+B collisions:
 505    continue
        ianti=0
        if(lb(i1).lt.0 .or. lb(i2).lt.0) ianti=1
        call sdmbb(SRT,sdm,ianti)
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
c     minimum srt**2, note a 2.012GeV lower cutoff is used in N+N->Deuteron+pi:
        EC=2.012**2
        ds=sqrt(sdm/31.4)
        dsr=ds+0.1
        CALL DISTCE(I1,I2,dsr,ds,DT,EC,SRT,IC,PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL crdmbb(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK,
     1       NTAG,sdm,NT,ianti)
        LCOLL=LCOLL+1
        GO TO 400
clin-8/2008 Deuteron+Meson->B+B collisions over
c
clin-9/2008 Deuteron+Baryon elastic collisions:
 506    continue
        ianti=0
        if(lb(i1).lt.0 .or. lb(i2).lt.0) ianti=1
        call sdbelastic(SRT,sdb)
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
c     minimum srt**2, note a 2.012GeV lower cutoff is used in N+N->Deuteron+pi:
        EC=2.012**2
        ds=sqrt(sdb/31.4)
        dsr=ds+0.1
        CALL DISTCE(I1,I2,dsr,ds,DT,EC,SRT,IC,PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL crdbel(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK,
     1       NTAG,sdb,NT,ianti)
        LCOLL=LCOLL+1
        GO TO 400
clin-9/2008 Deuteron+Baryon elastic collisions over
c
* IF BARYON RESONANCE+BARYON RESONANCE COLLISIONS
 444    CONTINUE
* PREPARE THE EALSTIC CROSS SECTION FOR BARYON+BARYON COLLISIONS
       CUTOFF=em1+em2+0.02
* AT HIGH ENERGIES THE ISOSPIN DEPENDENCE IS NEGLIGIBLE
* THE TOTAL CROSS SECTION IS TAKEN AS THAT OF THE PP 
       IF(SRT.LE.CUTOFF)GO TO 400
        IF(SRT.GT.2.245)THEN
       SIGNN=PP2(SRT)
       ELSE
        SIGNN = 35.0 / (1. + (SRT - CUTOFF) * 100.0)  +  20.0
       ENDIF 
       IF(SIGNN.LE.0)GO TO 400
      CALL XDDIN(PCX,PCY,PCZ,SRT,I1,I2,
     &XINEL,SIGK,XSK1,XSK2,XSK3,XSK4,XSK5)
       SIG=SIGNN+XINEL
       EC=(EM1+EM2+0.02)**2
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
clin-6/2008 Deuteron production:
        ianti=0
        if(lb(i1).lt.0 .and. lb(i2).lt.0) ianti=1
        call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
        sig=sig+sdprod
clin-6/2008 perturbative treatment of deuterons:
        ipdflag=0
        if(idpert.eq.1) then
           ipert1=1
           sigr0=sig
           dspert=sqrt(sigr0/pi/10.)
           dsrpert=dspert+0.1
           CALL DISTCE(I1,I2,dsrpert,dspert,DT,EC,SRT,IC,
     1          PX1CM,PY1CM,PZ1CM)
           IF(IC.EQ.-1) GO TO 367
           signn0=0.
           CALL CRDD(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1          IBLOCK,NTAG,SIGNN0,SIGr0,NT,ipert1)
c     1          IBLOCK,NTAG,SIGNN,SIG)
           ipdflag=1
 367       continue
           ipert1=0
        endif
        if(idpert.eq.2) ipert1=1
c
        ds=sqrt(sig/31.4)
        dsr=ds+0.1
        CALL DISTCE(I1,I2,dsr,ds,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
c        IF(IC.EQ.-1) GO TO 400
        IF(IC.EQ.-1) then
           if(ipdflag.eq.1) iblock=501
           GO TO 400
        endif
* CALCULATE KAON PRODUCTION PROBABILITY FROM NUCLEON+NUCLEON OR 
* RESONANCE+RESONANCE COLLISIONS
       go to 364
C CHECK WHAT KIND OF COLLISION HAS HAPPENED 
364       ekaon(2,iss)=ekaon(2,iss)+1
* for resonance+resonance
clin-6/2008:
        CALL CRDD(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK,NTAG,SIGNN,SIG,NT,ipert1)
c     1  IBLOCK,NTAG,SIGNN,SIG)
        IF(iblock.eq.0.and.ipdflag.eq.1) iblock=501
c
        IF(iabs(IBLOCK).EQ.10)THEN
* momentum of the three particles in the final state have been calculated
* in the crnn, go out of the loop
           LCOLL=LCOLL+1
           IF(IBLOCK.EQ.10)THEN
              LDDK=LDDK+1
           elseIF(IBLOCK.EQ.-10) then
           endif
           GO TO 400
        ENDIF
clin-6/2008
c        if(iblock .eq. 222)then
        if(iblock .eq. 222.or.iblock.eq.501)then
c    !! sp12/17/01 
           GO TO 400
        ENDIF
        em1=e(i1)
        em2=e(i2)
        GO TO 440
* FOR PION+PION,pion+eta, eta+eta and rho(omega)+pion(rho,omega) or eta 
777       CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
* energy thresh for collisions
       ec0=em1+em2+0.02
       IF(SRT.LE.ec0)GO TO 400
       ec=(em1+em2+0.02)**2
* we negelect the elastic collision between mesons except that betwen
* two pions because of the lack of information about these collisions
* However, we do let them to collide inelastically to produce kaons
clin-8/15/02       ppel=1.e-09
       ppel=20.
        ipp=1
       if(lb1.lt.3.or.lb1.gt.5.or.lb2.lt.3.or.lb2.gt.5)go to 778       
       CALL PPXS(LB1,LB2,SRT,PPSIG,spprho,IPP)
       ppel=ppsig
778       ppink=pipik(srt)
* pi+eta and eta+eta are assumed to be the same as pipik( for pi+pi -> K+K-) 
* estimated from Ko's paper:
        ppink = 2.0 * ppink
       if(lb1.ge.25.and.lb2.ge.25) ppink=rrkk
clin-2/13/03 include omega the same as rho, eta the same as pi:
c        if(((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.25.and.lb2.le.27))
c     1  .or.((lb2.ge.3.and.lb2.le.5).and.(lb1.ge.25.and.lb1.le.27)))
        if( ( (lb1.eq.0.or.(lb1.ge.3.and.lb1.le.5))
     1       .and.(lb2.ge.25.and.lb2.le.28))
     2       .or. ( (lb2.eq.0.or.(lb2.ge.3.and.lb2.le.5))
     3       .and.(lb1.ge.25.and.lb1.le.28))) then
           ppink=0.
           if(srt.ge.(aka+aks)) ppink = prkk
        endif
c pi pi <-> rho rho:
        call spprr(lb1,lb2,srt)
clin-4/03/02 pi pi <-> eta eta:
        call sppee(lb1,lb2,srt)
clin-4/03/02 pi pi <-> pi eta:
        call spppe(lb1,lb2,srt)
clin-4/03/02 rho pi <-> rho eta:
        call srpre(lb1,lb2,srt)
clin-4/03/02 omega pi <-> omega eta:
        call sopoe(lb1,lb2,srt)
clin-4/03/02 rho rho <-> eta eta:
        call srree(lb1,lb2,srt)
        ppinnb=0.
        if(srt.gt.thresh(1)) then
           call getnst(srt)
           if(lb1.ge.3.and.lb1.le.5.and.lb2.ge.3.and.lb2.le.5) then
              ppinnb=ppbbar(srt)
           elseif((lb1.ge.3.and.lb1.le.5.and.lb2.ge.25.and.lb2.le.27)
     1 .or.(lb2.ge.3.and.lb2.le.5.and.lb1.ge.25.and.lb1.le.27)) then
              ppinnb=prbbar(srt)
           elseif(lb1.ge.25.and.lb1.le.27
     1             .and.lb2.ge.25.and.lb2.le.27) then
              ppinnb=rrbbar(srt)
           elseif((lb1.ge.3.and.lb1.le.5.and.lb2.eq.28)
     1             .or.(lb2.ge.3.and.lb2.le.5.and.lb1.eq.28)) then
              ppinnb=pobbar(srt)
           elseif((lb1.ge.25.and.lb1.le.27.and.lb2.eq.28)
     1             .or.(lb2.ge.25.and.lb2.le.27.and.lb1.eq.28)) then
              ppinnb=robbar(srt)
           elseif(lb1.eq.28.and.lb2.eq.28) then
              ppinnb=oobbar(srt)
           else
              if(lb1.ne.0.and.lb2.ne.0) 
     1             write(6,*) 'missed MM lb1,lb2=',lb1,lb2
           endif
        endif
        ppin=ppink+ppinnb+pprr+ppee+pppe+rpre+xopoe+rree
* check if a collision can happen
       if((ppel+ppin).le.0.01)go to 400
       DSPP=SQRT((ppel+ppin)/31.4)
       dsppr=dspp+0.1
        CALL DISTCE(I1,I2,dsppr,DSPP,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
       if(ppel.eq.0)go to 400
* the collision can happen
* check what kind collision has happened
       ekaon(5,iss)=ekaon(5,iss)+1
        CALL CRPP(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK,ppel,ppin,spprho,ipp)
* rho formation, go to 400
c       if(iblock.eq.666)go to 600
       if(iblock.eq.666)go to 555
       if(iblock.eq.6)LPP=LPP+1
       if(iblock.eq.66)then
          LPPk=LPPk+1
       elseif(iblock.eq.366)then
          LPPk=LPPk+1
       elseif(iblock.eq.367)then
          LPPk=LPPk+1
       endif
       em1=e(i1)
       em2=e(i2)
       go to 440
* In this block we treat annihilations of
clin-9/28/00* an anti-nucleon and a baryon or baryon resonance  
* an anti-baryon and a baryon (including resonances)
2799        CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
clin assume the same cross section (as a function of sqrt s) as for PPbar:
clin-ctest annih maximum
c        DSppb=SQRT(amin1(xppbar(srt),30.)/PI/10.)
       DSppb=SQRT(xppbar(srt)/PI/10.)
       dsppbr=dsppb+0.1
        CALL DISTCE(I1,I2,dsppbr,DSppb,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL Crppba(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK)
       em1=e(i1)
       em2=e(i2)
       go to 440
c
3555    PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
       DSkk=SQRT(SIG/PI/10.)
       dskk0=dskk+0.1
        CALL DISTCE(I1,I2,dskk0,DSkk,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL Crlaba(PX1CM,PY1CM,PZ1CM,SRT,brel,brsgm,
     &                  I1,I2,nt,IBLOCK,nchrg,icase)
       em1=e(i1)
       em2=e(i2)
       go to 440
*
c perturbative production of cascade and omega
3455    PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        call pertur(PX1CM,PY1CM,PZ1CM,SRT,IRUN,I1,I2,nt,kp,icontp)
        if(icontp .eq. 0)then
c     inelastic collisions:
         em1 = e(i1)
         em2 = e(i2)
         iblock = 727
          go to 440
        endif
c     elastic collisions:
        if (e(i1) .eq. 0.) go to 800
        if (e(i2) .eq. 0.) go to 600
        go to 400
*
c* phi + N --> pi+N(D),  N(D,N*)+N(D,N*),  K+ +La
c* phi + D --> pi+N(D)
7222        CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
        CALL XphiB(LB1, LB2, EM1, EM2, SRT,
     &             XSK1, XSK2, XSK3, XSK4, XSK5, SIGP)
       DSkk=SQRT(SIGP/PI/10.)
       dskk0=dskk+0.1
        CALL DISTCE(I1,I2,dskk0,DSkk,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL CRPHIB(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     &     XSK1, XSK2, XSK3, XSK4, XSK5, SIGP, IBLOCK)
       em1=e(i1)
       em2=e(i2)
       go to 440
*
c* phi + M --> K+ + K* .....
7444        CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
        CALL PHIMES(I1, I2, SRT, XSK1, XSK2, XSK3, XSK4, XSK5,
     1     XSK6, XSK7, SIGPHI)
       DSkk=SQRT(SIGPHI/PI/10.)
       dskk0=dskk+0.1
        CALL DISTCE(I1,I2,dskk0,DSkk,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
c*---
        PZRT = p(3,i1)+p(3,i2)
        ER1 = sqrt( p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+E(i1)**2 )
        ER2 = sqrt( p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+E(i2)**2 )
        ERT = ER1+ER2
        yy = 0.5*log( (ERT+PZRT)/(ERT-PZRT) )
c*------
        CALL CRPHIM(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     &  XSK1, XSK2, XSK3, XSK4, XSK5, XSK6, SIGPHI, IKKG, IKKL, IBLOCK)
       em1=e(i1)
       em2=e(i2)
       go to 440
c
c lambda-N elastic xsection, Li & Ko, PRC 54(1996)1897.
 7799    CONTINUE
         PX1CM=PCX
         PY1CM=PCY
         PZ1CM=PCZ
         EC=(em1+em2+0.02)**2
         call lambar(i1,i2,srt,siglab)
        DShn=SQRT(siglab/PI/10.)
        dshnr=dshn+0.1
         CALL DISTCE(I1,I2,dshnr,DShn,DT,EC,SRT,IC,
     1    PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
         CALL Crhb(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)
        em1=e(i1)
        em2=e(i2)
        go to 440
c
c* K+ + La(Si) --> Meson + B
c* K- + La(Si)-bar --> Meson + B-bar
5699        CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
        CALL XKHYPE(I1, I2, SRT, XKY1, XKY2, XKY3, XKY4, XKY5,
     &     XKY6, XKY7, XKY8, XKY9, XKY10, XKY11, XKY12, XKY13,
     &     XKY14, XKY15, XKY16, XKY17, SIGK)
       DSkk=SQRT(sigk/PI)
       dskk0=dskk+0.1
        CALL DISTCE(I1,I2,dskk0,DSkk,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
c
       if(lb(i1).eq.23 .or. lb(i2).eq.23)then
             IKMP = 1
        else
             IKMP = -1
        endif
        CALL Crkhyp(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     &     XKY1, XKY2, XKY3, XKY4, XKY5,
     &     XKY6, XKY7, XKY8, XKY9, XKY10, XKY11, XKY12, XKY13,
     &     XKY14, XKY15, XKY16, XKY17, SIGK, IKMP,
     1  IBLOCK)
       em1=e(i1)
       em2=e(i2)
       go to 440
c khyperon end
*
csp11/03/01 La/Si-bar + N --> pi + K+
c  La/Si + N-bar --> pi + K-
5999     CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
        sigkp = 15.
c      if((lb1.ge.14.and.lb1.le.17)
c     &    .or.(lb2.ge.14.and.lb2.le.17))sigkp=10.
        DSkk=SQRT(SIGKP/PI/10.)
        dskk0=dskk+0.1
        CALL DISTCE(I1,I2,dskk0,DSkk,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
c
        CALL CRLAN(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)
        em1=e(i1)
        em2=e(i2)
        go to 440
c
c*
* K(K*) + K(K*) --> phi + pi(rho,omega)
8699     CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
*  CALL CROSSKKPHI(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)  used for KK*->phi+rho
         CALL Crkphi(PX1CM,PY1CM,PZ1CM,EC,SRT,IBLOCK,
     &                  emm1,emm2,lbp1,lbp2,I1,I2,ikk,icase,rrkk,prkk)
         if(icase .eq. 0) then
            iblock=0
            go to 400
         endif
c*---
         if(lbp1.eq.29.or.lbp2.eq.29) then
        PZRT = p(3,i1)+p(3,i2)
        ER1 = sqrt( p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+E(i1)**2 )
        ER2 = sqrt( p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+E(i2)**2 )
        ERT = ER1+ER2
        yy = 0.5*log( (ERT+PZRT)/(ERT-PZRT) )
c*------
             iblock = 222
             ntag = 0
          endif
             LB(I1) = lbp1
             LB(I2) = lbp2
             E(I1) = emm1
             E(I2) = emm2
             em1=e(i1)
             em2=e(i2)
             go to 440
c*
* rho(omega) + K(K*)  --> phi + K(K*)
8799     CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
*  CALL CROSSKKPHI(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)  used for KK*->phi+rho
         CALL Crksph(PX1CM,PY1CM,PZ1CM,EC,SRT,
     &       emm1,emm2,lbp1,lbp2,I1,I2,ikkg,ikkl,iblock,icase,srhoks)
         if(icase .eq. 0) then
            iblock=0
            go to 400
         endif
c
         if(lbp1.eq.29.or.lbp2.eq.20) then
c*---
        PZRT = p(3,i1)+p(3,i2)
        ER1 = sqrt( p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+E(i1)**2 )
        ER2 = sqrt( p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+E(i2)**2 )
        ERT = ER1+ER2
        yy = 0.5*log( (ERT+PZRT)/(ERT-PZRT) )
          endif
             LB(I1) = lbp1
             LB(I2) = lbp2
             E(I1) = emm1
             E(I2) = emm2
             em1=e(i1)
             em2=e(i2)
             go to 440
* for kaon+baryon scattering, using a constant xsection of 10 mb.
888       CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
         sig = 10.
         if(iabs(lb1).eq.14.or.iabs(lb2).eq.14 .or.
     &      iabs(lb1).eq.30.or.iabs(lb2).eq.30)sig=20.
         if(lb1.eq.29.or.lb2.eq.29)sig=5.0
       DSkn=SQRT(sig/PI/10.)
       dsknr=dskn+0.1
        CALL DISTCE(I1,I2,dsknr,DSkn,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL Crkn(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK)
       em1=e(i1)
       em2=e(i2)
       go to 440
***
 440    CONTINUE
*                IBLOCK = 0 ; NOTHING HAS HAPPENED
*                IBLOCK = 1 ; ELASTIC N-N COLLISION
*                IBLOCK = 2 ; N + N -> N + DELTA
*                IBLOCK = 3 ; N + DELTA -> N + N
*                IBLOCK = 4 ; N + N -> d + d + PION,DIRECT PROCESS
*               IBLOCK = 5 ; D(N*)+D(N*) COLLISIONS
*                IBLOCK = 6 ; PION+PION COLLISIONS
*                iblock = 7 ; pion+nucleon-->l/s+kaon
*               iblock =77;  pion+nucleon-->delta+pion
*               iblock = 8 ; kaon+baryon rescattering
*                IBLOCK = 9 ; NN-->KAON+X
*                IBLOCK = 10; DD-->KAON+X
*               IBLOCK = 11; ND-->KAON+X
cbali2/1/99
*                
*           iblock   - 1902 annihilation-->pion(+)+pion(-)   (2 pion)
*           iblock   - 1903 annihilation-->pion(+)+rho(-)    (3 pion)
*           iblock   - 1904 annihilation-->rho(+)+rho(-)     (4 pion)
*           iblock   - 1905 annihilation-->rho(0)+omega      (5 pion)
*           iblock   - 1906 annihilation-->omega+omega       (6 pion)
cbali3/5/99
*           iblock   - 1907 K+K- to pi+pi-
cbali3/5/99 end
cbz3/9/99 khyperon
*           iblock   - 1908 K+Y -> piN
cbz3/9/99 khyperon end
cbali2/1/99end
clin-9/28/00 Processes: m(pi rho omega)+m(pi rho omega)
c     to anti-(p n D N*1 N*2)+(p n D N*1 N*2):
*           iblock   - 1801  mm -->pbar p 
*           iblock   - 18021 mm -->pbar n 
*           iblock   - 18022 mm -->nbar p 
*           iblock   - 1803  mm -->nbar n 
*           iblock   - 18041 mm -->pbar Delta 
*           iblock   - 18042 mm -->anti-Delta p
*           iblock   - 18051 mm -->nbar Delta 
*           iblock   - 18052 mm -->anti-Delta n
*           iblock   - 18061 mm -->pbar N*(1400) 
*           iblock   - 18062 mm -->anti-N*(1400) p
*           iblock   - 18071 mm -->nbar N*(1400)
*           iblock   - 18072 mm -->anti-N*(1400) n
*           iblock   - 1808  mm -->anti-Delta Delta 
*           iblock   - 18091 mm -->pbar N*(1535)
*           iblock   - 18092 mm -->anti-N*(1535) p
*           iblock   - 18101 mm -->nbar N*(1535)
*           iblock   - 18102 mm -->anti-N*(1535) n
*           iblock   - 18111 mm -->anti-Delta N*(1440)
*           iblock   - 18112 mm -->anti-N*(1440) Delta
*           iblock   - 18121 mm -->anti-Delta N*(1535)
*           iblock   - 18122 mm -->anti-N*(1535) Delta
*           iblock   - 1813  mm -->anti-N*(1440) N*(1440)
*           iblock   - 18141 mm -->anti-N*(1440) N*(1535)
*           iblock   - 18142 mm -->anti-N*(1535) N*(1440)
*           iblock   - 1815  mm -->anti-N*(1535) N*(1535)
clin-9/28/00-end
clin-10/08/00 Processes: pi pi <-> rho rho
*           iblock   - 1850  pi pi -> rho rho
*           iblock   - 1851  rho rho -> pi pi
clin-10/08/00-end
clin-08/14/02 Processes: pi pi <-> eta eta
*           iblock   - 1860  pi pi -> eta eta
*           iblock   - 1861  eta eta -> pi pi
* Processes: pi pi <-> pi eta
*           iblock   - 1870  pi pi -> pi eta
*           iblock   - 1871  pi eta -> pi pi
* Processes: rho pi <-> rho eta
*           iblock   - 1880  pi pi -> pi eta
*           iblock   - 1881  pi eta -> pi pi
* Processes: omega pi <-> omega eta
*           iblock   - 1890  pi pi -> pi eta
*           iblock   - 1891  pi eta -> pi pi
* Processes: rho rho <-> eta eta
*           iblock   - 1895  rho rho -> eta eta
*           iblock   - 1896  eta eta -> rho rho
clin-08/14/02-end
clin-11/07/00 Processes: 
*           iblock   - 366  pi rho -> K* Kbar or K*bar K
*           iblock   - 466  pi rho <- K* Kbar or K*bar K
clin-9/2008 Deuteron:
*           iblock   - 501  B+B -> Deuteron+Meson
*           iblock   - 502  Deuteron+Meson -> B+B
*           iblock   - 503  Deuteron+Baryon elastic
*           iblock   - 504  Deuteron+Meson elastic
c
                 IF(IBLOCK.EQ.0)        GOTO 400
*COM: FOR DIRECT PROCESS WE HAVE TREATED THE PAULI BLOCKING AND FIND
*     THE MOMENTUM OF PARTICLES IN THE ''LAB'' FRAME. SO GO TO 400
* A COLLISION HAS TAKEN PLACE !!
              LCOLL = LCOLL +1
* WAS COLLISION PAULI-FORBIDEN? IF YES, NTAG = -1
              NTAG = 0
*
*             LORENTZ-TRANSFORMATION INTO CMS FRAME
              E1CM    = SQRT (EM1**2 + PX1CM**2 + PY1CM**2 + PZ1CM**2)
              P1BETA  = PX1CM*BETAX + PY1CM*BETAY + PZ1CM*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1I1 = BETAX * TRANSF + PX1CM
              Pt2I1 = BETAY * TRANSF + PY1CM
              Pt3I1 = BETAZ * TRANSF + PZ1CM
* negelect the pauli blocking at high energies
              go to 90002
clin-10/25/02-comment out following, since there is no path to it:
c*CHECK IF PARTICLE #1 IS PAULI BLOCKED
c              CALL PAULat(I1,occup)
c              if (RANART(NSEED) .lt. occup) then
c                ntag = -1
c              else
c                ntag = 0
c              end if
clin-10/25/02-end
90002              continue
*IF PARTICLE #1 IS NOT PAULI BLOCKED
c              IF (NTAG .NE. -1) THEN
                E2CM    = SQRT (EM2**2 + PX1CM**2 + PY1CM**2 + PZ1CM**2)
                TRANSF  = GAMMA * (-GAMMA*P1BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF - PX1CM
                Pt2I2 = BETAY * TRANSF - PY1CM
                Pt3I2 = BETAZ * TRANSF - PZ1CM
              go to 90003
clin-10/25/02-comment out following, since there is no path to it:
c*CHECK IF PARTICLE #2 IS PAULI BLOCKED
c                CALL PAULat(I2,occup)
c                if (RANART(NSEED) .lt. occup) then
c                  ntag = -1
c                else
c                  ntag = 0
c                end if
cc              END IF
c* IF COLLISION IS BLOCKED,RESTORE THE MOMENTUM,MASSES
c* AND LABELS OF I1 AND I2
cc             IF (NTAG .EQ. -1) THEN
c                LBLOC  = LBLOC + 1
c                P(1,I1) = PX1
c                P(2,I1) = PY1
c                P(3,I1) = PZ1
c                P(1,I2) = PX2
c                P(2,I2) = PY2
c                P(3,I2) = PZ2
c                E(I1)   = EM1
c                E(I2)   = EM2
c                LB(I1)  = LB1
c                LB(I2)  = LB2
cc              ELSE
clin-10/25/02-end
90003           IF(IBLOCK.EQ.1) LCNNE=LCNNE+1
              IF(IBLOCK.EQ.5) LDD=LDD+1
                if(iblock.eq.2) LCNND=LCNND+1
              IF(IBLOCK.EQ.8) LKN=LKN+1
                   if(iblock.eq.43) Ldou=Ldou+1
c                IF(IBLOCK.EQ.2) THEN
* CALCULATE THE AVERAGE SRT FOR N + N---> N + DELTA PROCESS
C                NODELT=NODELT+1
C                SUMSRT=SUMSRT+SRT
c                ENDIF
                IF(IBLOCK.EQ.3) LCNDN=LCNDN+1
* assign final momenta to particles while keep the leadng particle
* behaviour
C              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
C              else
C              p(1,i1)=pt1i2
C              p(2,i1)=pt2i2
C              p(3,i1)=pt3i2
C              p(1,i2)=pt1i1
C              p(2,i2)=pt2i1
C              p(3,i2)=pt3i1
C              endif
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
                EM1     = E(I1)
                EM2     = E(I2)
                LB1     = LB(I1)
                LB2     = LB(I2)
                ID(I1)  = 2
                ID(I2)  = 2
                E1      = SQRT( EM1**2 + PX1**2 + PY1**2 + PZ1**2 )
                ID1     = ID(I1)
              go to 90004
clin-10/25/02-comment out following, since there is no path to it:
c* change phase space density FOR NUCLEONS INVOLVED :
c* NOTE THAT f is the phase space distribution function for nucleons only
c                if ((abs(ix1).le.mx) .and. (abs(iy1).le.my) .and.
c     &              (abs(iz1).le.mz)) then
c                  ipx1p = nint(p(1,i1)/dpx)
c                  ipy1p = nint(p(2,i1)/dpy)
c                  ipz1p = nint(p(3,i1)/dpz)
c                  if ((ipx1p.ne.ipx1) .or. (ipy1p.ne.ipy1) .or.
c     &                (ipz1p.ne.ipz1)) then
c                    if ((abs(ipx1).le.mpx) .and. (abs(ipy1).le.my)
c     &                .and. (ipz1.ge.-mpz) .and. (ipz1.le.mpzp)
c     &                .AND. (AM1.LT.1.))
c     &                f(ix1,iy1,iz1,ipx1,ipy1,ipz1) =
c     &                f(ix1,iy1,iz1,ipx1,ipy1,ipz1) - 1.
c                    if ((abs(ipx1p).le.mpx) .and. (abs(ipy1p).le.my)
c     &                .and. (ipz1p.ge.-mpz).and. (ipz1p.le.mpzp)
c     &                .AND. (EM1.LT.1.))
c     &                f(ix1,iy1,iz1,ipx1p,ipy1p,ipz1p) =
c     &                f(ix1,iy1,iz1,ipx1p,ipy1p,ipz1p) + 1.
c                  end if
c                end if
c                if ((abs(ix2).le.mx) .and. (abs(iy2).le.my) .and.
c     &              (abs(iz2).le.mz)) then
c                  ipx2p = nint(p(1,i2)/dpx)
c                  ipy2p = nint(p(2,i2)/dpy)
c                  ipz2p = nint(p(3,i2)/dpz)
c                  if ((ipx2p.ne.ipx2) .or. (ipy2p.ne.ipy2) .or.
c     &                (ipz2p.ne.ipz2)) then
c                    if ((abs(ipx2).le.mpx) .and. (abs(ipy2).le.my)
c     &                .and. (ipz2.ge.-mpz) .and. (ipz2.le.mpzp)
c     &                .AND. (AM2.LT.1.))
c     &                f(ix2,iy2,iz2,ipx2,ipy2,ipz2) =
c     &                f(ix2,iy2,iz2,ipx2,ipy2,ipz2) - 1.
c                    if ((abs(ipx2p).le.mpx) .and. (abs(ipy2p).le.my)
c     &                .and. (ipz2p.ge.-mpz) .and. (ipz2p.le.mpzp)
c     &                .AND. (EM2.LT.1.))
c     &                f(ix2,iy2,iz2,ipx2p,ipy2p,ipz2p) =
c     &                f(ix2,iy2,iz2,ipx2p,ipy2p,ipz2p) + 1.
c                  end if
c                end if
clin-10/25/02-end
90004              continue
            AM1=EM1
            AM2=EM2
c            END IF
  400       CONTINUE
c
clin-6/10/03 skips the info output on resonance creations:
c            goto 550
cclin-4/30/03 study phi,K*,Lambda(1520) resonances at creation:
cc     note that no decays give these particles, so don't need to consider nnn:
c            if(iblock.ne.0.and.(lb(i1).eq.29.or.iabs(lb(i1)).eq.30
c     1           .or.lb(i2).eq.29.or.iabs(lb(i2)).eq.30
c     2           .or.lb1i.eq.29.or.iabs(lb1i).eq.30
c     3           .or.lb2i.eq.29.or.iabs(lb2i).eq.30)) then
c               lb1now=lb(i1)
c               lb2now=lb(i2)
cc
c               nphi0=0
c               nksp0=0
c               nksm0=0
cc               nlar0=0
cc               nlarbar0=0
c               if(lb1i.eq.29) then
c                  nphi0=nphi0+1
c               elseif(lb1i.eq.30) then
c                  nksp0=nksp0+1
c               elseif(lb1i.eq.-30) then
c                  nksm0=nksm0+1
c               endif
c               if(lb2i.eq.29) then
c                  nphi0=nphi0+1
c               elseif(lb2i.eq.30) then
c                  nksp0=nksp0+1
c               elseif(lb2i.eq.-30) then
c                  nksm0=nksm0+1
c               endif
cc
c               nphi=0
c               nksp=0
c               nksm=0
c               nlar=0
c               nlarbar=0
c               if(lb1now.eq.29) then
c                  nphi=nphi+1
c               elseif(lb1now.eq.30) then
c                  nksp=nksp+1
c               elseif(lb1now.eq.-30) then
c                  nksm=nksm+1
c               endif
c               if(lb2now.eq.29) then
c                  nphi=nphi+1
c               elseif(lb2now.eq.30) then
c                  nksp=nksp+1
c               elseif(lb2now.eq.-30) then
c                  nksm=nksm+1
c               endif
cc     
c               if(nphi.eq.2.or.nksp.eq.2.or.nksm.eq.2) then
c                  write(91,*) '2 same resonances in one reaction!'
c                  write(91,*) nphi,nksp,nksm,iblock
c               endif
c
cc     All reactions create or destroy no more than 1 these resonance,
cc     otherwise file "fort.91" warns us:
c               do 222 ires=1,3
c                  if(ires.eq.1.and.nphi.ne.nphi0) then
c                     idr=29
c                  elseif(ires.eq.2.and.nksp.ne.nksp0) then
c                     idr=30
c                  elseif(ires.eq.3.and.nksm.ne.nksm0) then
c                     idr=-30
c                  else
c                     goto 222
c                  endif
cctest off for resonance (phi, K*) studies:
cc               if(lb1now.eq.idr) then
cc       write(17,112) 'collision',lb1now,P(1,I1),P(2,I1),P(3,I1),e(I1),nt
cc               elseif(lb2now.eq.idr) then
cc       write(17,112) 'collision',lb2now,P(1,I2),P(2,I2),P(3,I2),e(I2),nt
cc               elseif(lb1i.eq.idr) then
cc       write(18,112) 'collision',lb1i,px1i,py1i,pz1i,em1i,nt
cc               elseif(lb2i.eq.idr) then
cc       write(18,112) 'collision',lb2i,px2i,py2i,pz2i,em2i,nt
cc               endif
c 222           continue
c
c            else
c            endif
cc 112        format(a10,I4,4(1x,f9.3),1x,I4)
c
clin-2/26/03 skips the check of energy conservation after each binary search:
c 550        goto 555
c            pxfin=0
c            pyfin=0
c            pzfin=0
c            efin=0
c            if(e(i1).ne.0.or.lb(i1).eq.10022) then
c               efin=efin+SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2)
c               pxfin=pxfin+P(1,I1)
c               pyfin=pyfin+P(2,I1)
c               pzfin=pzfin+P(3,I1)
c            endif
c            if(e(i2).ne.0.or.lb(i2).eq.10022) then
c               efin=efin+SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2)
c               pxfin=pxfin+P(1,I2)
c               pyfin=pyfin+P(2,I2)
c               pzfin=pzfin+P(3,I2)
c            endif
c            if((nnn-nnnini).ge.1) then
c               do imore=nnnini+1,nnn
c                  if(EPION(imore,IRUN).ne.0) then
c                     efin=efin+SQRT(EPION(imore,IRUN)**2
c     1                    +PPION(1,imore,IRUN)**2+PPION(2,imore,IRUN)**2
c     2                    +PPION(3,imore,IRUN)**2)
c                     pxfin=pxfin+PPION(1,imore,IRUN)
c                     pyfin=pyfin+PPION(2,imore,IRUN)
c                     pzfin=pzfin+PPION(3,imore,IRUN)
c                  endif
c               enddo
c            endif
c            devio=sqrt((pxfin-pxini)**2+(pyfin-pyini)**2
c     1           +(pzfin-pzini)**2+(efin-eini)**2)
cc
c            if(devio.ge.0.1) then
c               write(92,'a20,5(1x,i6),2(1x,f8.3)') 'iblock,lb,npi=',
c     1              iblock,lb1i,lb2i,lb(i1),lb(i2),e(i1),e(i2)
c               do imore=nnnini+1,nnn
c                  if(EPION(imore,IRUN).ne.0) then
c                     write(92,'a10,2(1x,i6)') 'ipi,lbm=',
c     1                    imore,LPION(imore,IRUN)
c                  endif
c               enddo
c               write(92,'a3,4(1x,f8.3)') 'I:',eini,pxini,pyini,pzini
c               write(92,'a3,5(1x,f8.3)') 
c     1              'F:',efin,pxfin,pyfin,pzfin,devio
c            endif
c
 555        continue
ctest off only one collision for the same 2 particles in the same timestep:
c            if(iblock.ne.0) then
c               goto 800
c            endif
ctest off collisions history:
c            if(iblock.ne.0) then 
c               write(10,*) nt,i1,i2,iblock,x1,z1,x2,z2
c            endif
  600     CONTINUE
clin-4/2012 option of pi0 decays:
c     particles in lpion() may be a pi0, and when ipi0dcy=1 
c     we need to decay them at nt=ntmax after all lb(i1) decays are done:
 798      if(nt.eq.ntmax.and.ipi0dcy.eq.1
     1         .and.i1.eq.(MASSR(IRUN)+MSUM)) then
             do ipion=1,NNN
                if(LPION(ipion,IRUN).eq.4) then
                   wid=7.85e-9
                   call resdec(i1,nt,nnn,wid,idecay,ipion)
                endif
             enddo
          endif
ctest off
c          if(nt.eq.ntmax.and.i1.eq.(MASSR(IRUN)+MSUM)) then
c             do ip=1,i1
c                write(98,*) lb(ip),e(ip),ip
c             enddo
c          endif
clin-4/2012 option of pi0 decays-end
  800   CONTINUE
* RELABLE MESONS LEFT IN THIS RUN EXCLUDING THOSE BEING CREATED DURING
* THIS TIME STEP AND COUNT THE TOTAL NO. OF PARTICLES IN THIS RUN
* note that the first mass=mta+mpr particles are baryons
c        write(*,*)'I: NNN,massr ', nnn,massr(irun)
        N0=MASS+MSUM
        DO 1005 N=N0+1,MASSR(IRUN)+MSUM
cbz11/25/98
clin-2/19/03 lb>5000: keep particles with no LB codes in ART(photon,lepton,..):
c        IF(E(N).GT.0.)THEN
        IF(E(N) .GT. 0. .OR. LB(N) .GT. 5000)THEN
cbz11/25/98end
        NNN=NNN+1
        RPION(1,NNN,IRUN)=R(1,N)
        RPION(2,NNN,IRUN)=R(2,N)
        RPION(3,NNN,IRUN)=R(3,N)
clin-10/28/03:
        if(nt.eq.ntmax) then
           ftpisv(NNN,IRUN)=ftsv(N)
           tfdpi(NNN,IRUN)=tfdcy(N)
        endif
c
        PPION(1,NNN,IRUN)=P(1,N)
        PPION(2,NNN,IRUN)=P(2,N)
        PPION(3,NNN,IRUN)=P(3,N)
        EPION(NNN,IRUN)=E(N)
        LPION(NNN,IRUN)=LB(N)
c       !! sp 12/19/00
        PROPI(NNN,IRUN)=PROPER(N)
clin-5/2008:
        dppion(NNN,IRUN)=dpertp(N)
c        if(lb(n) .eq. 45)
c    &   write(*,*)'IN-1  NT,NNN,LB,P ',nt,NNN,lb(n),proper(n)
        ENDIF
 1005 CONTINUE
        MASSRN(IRUN)=NNN+MASS
c        write(*,*)'F: NNN,massrn ', nnn,massrn(irun)
1000   CONTINUE
* CALCULATE THE AVERAGE SRT FOR N + N--->N +DELTA PROCESSES
C        IF(NODELT.NE.0)THEN
C        AVSRT=SUMSRT/FLOAT(NODELT)
C        ELSE
C        AVSRT=0.
C        ENDIF
C        WRITE(1097,'(F8.2,2X,E10.3)')FLOAT(NT)*DT,AVSRT
* RELABLE ALL THE PARTICLES EXISTING AFTER THIS TIME STEP
        IA=0
        IB=0
        DO 10001 IRUN=1,NUM
        IA=IA+MASSR(IRUN-1)
        IB=IB+MASSRN(IRUN-1)
        DO 10001 IC=1,MASSRN(IRUN)
        IE=IA+IC
        IG=IB+IC
        IF(IC.LE.MASS)THEN
        RT(1,IG)=R(1,IE)
        RT(2,IG)=R(2,IE)
        RT(3,IG)=R(3,IE)
clin-10/28/03:
        if(nt.eq.ntmax) then
           fttemp(IG)=ftsv(IE)
           tft(IG)=tfdcy(IE)
        endif
c
        PT(1,IG)=P(1,IE)
        PT(2,IG)=P(2,IE)
        PT(3,IG)=P(3,IE)
        ET(IG)=E(IE)
        LT(IG)=LB(IE)
        PROT(IG)=PROPER(IE)
clin-5/2008:
        dptemp(IG)=dpertp(IE)
        ELSE
        I0=IC-MASS
        RT(1,IG)=RPION(1,I0,IRUN)
        RT(2,IG)=RPION(2,I0,IRUN)
        RT(3,IG)=RPION(3,I0,IRUN)
clin-10/28/03:
        if(nt.eq.ntmax) then
           fttemp(IG)=ftpisv(I0,IRUN)
           tft(IG)=tfdpi(I0,IRUN)
        endif
c
        PT(1,IG)=PPION(1,I0,IRUN)
        PT(2,IG)=PPION(2,I0,IRUN)
        PT(3,IG)=PPION(3,I0,IRUN)
        ET(IG)=EPION(I0,IRUN)
        LT(IG)=LPION(I0,IRUN)
        PROT(IG)=PROPI(I0,IRUN)
clin-5/2008:
        dptemp(IG)=dppion(I0,IRUN)
        ENDIF
10001   CONTINUE
c
        IL=0
clin-10/26/01-hbt:
c        DO 10002 IRUN=1,NUM
        DO 10003 IRUN=1,NUM
        MASSR(IRUN)=MASSRN(IRUN)
        IL=IL+MASSR(IRUN-1)
        DO 10002 IM=1,MASSR(IRUN)
        IN=IL+IM
        R(1,IN)=RT(1,IN)
        R(2,IN)=RT(2,IN)
        R(3,IN)=RT(3,IN)
clin-10/28/03:
        if(nt.eq.ntmax) then
           ftsv(IN)=fttemp(IN)
           tfdcy(IN)=tft(IN)
        endif
        P(1,IN)=PT(1,IN)
        P(2,IN)=PT(2,IN)
        P(3,IN)=PT(3,IN)
        E(IN)=ET(IN)
        LB(IN)=LT(IN)
        PROPER(IN)=PROT(IN)
clin-5/2008:
        dpertp(IN)=dptemp(IN)
       IF(LB(IN).LT.1.OR.LB(IN).GT.2)ID(IN)=0
10002   CONTINUE
clin-ctest off check energy conservation after each timestep
c         enetot=0.
c         do ip=1,MASSR(IRUN)
c            if(e(ip).ne.0.or.lb(ip).eq.10022) enetot=enetot
c     1           +sqrt(p(1,ip)**2+p(2,ip)**2+p(3,ip)**2+e(ip)**2)
c         enddo
c         write(91,*) 'B:',nt,enetot,massr(irun),bimp 
clin-3/2009 move to the end of a timestep to take care of freezeout spacetime:
c        call hbtout(MASSR(IRUN),nt,ntmax)
10003 CONTINUE
c
      RETURN
      END
clin-9/2012: use double precision for S in CMS(): to avoid crash 
c     (segmentation fault due to s<0, which happened at high energies 
c     such as LHC with large NTMAX for two almost-comoving hadrons
c     that have small Pt but large |Pz|):
****************************************
c            SUBROUTINE CMS(I1,I2,PX1CM,PY1CM,PZ1CM,SRT)
* PURPOSE : FIND THE MOMENTA OF PARTICLES IN THE CMS OF THE
*          TWO COLLIDING PARTICLES
* VARIABLES :
*****************************************
c            PARAMETER (MAXSTR=150001)
c            COMMON   /AA/  R(3,MAXSTR)
ccc      SAVE /AA/
c            COMMON   /BB/  P(3,MAXSTR)
ccc      SAVE /BB/
c            COMMON   /CC/  E(MAXSTR)
ccc      SAVE /CC/
c            COMMON   /BG/  BETAX,BETAY,BETAZ,GAMMA
ccc      SAVE /BG/
c            SAVE   
c            PX1=P(1,I1)
c            PY1=P(2,I1)
c            PZ1=P(3,I1)
c            PX2=P(1,I2)
c            PY2=P(2,I2)
c            PZ2=P(3,I2)
c            EM1=E(I1)
c            EM2=E(I2)
c            E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
c            E2=SQRT(EM2**2 + PX2**2 + PY2**2 + PZ2**2 )
c            S=(E1+E2)**2-(PX1+PX2)**2-(PY1+PY2)**2-(PZ1+PZ2)**2
c            SRT=SQRT(S)
c*LORENTZ-TRANSFORMATION IN I1-I2-C.M. SYSTEM
c              ETOTAL = E1 + E2
c              BETAX  = (PX1+PX2) / ETOTAL
c              BETAY  = (PY1+PY2) / ETOTAL
c              BETAZ  = (PZ1+PZ2) / ETOTAL
c              GAMMA  = 1.0 / SQRT(1.0-BETAX**2-BETAY**2-BETAZ**2)
c*TRANSFORMATION OF MOMENTA (PX1CM = - PX2CM)
c              P1BETA = PX1*BETAX + PY1*BETAY + PZ1 * BETAZ
c              TRANSF = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) - E1 )
c              PX1CM  = BETAX * TRANSF + PX1
c              PY1CM  = BETAY * TRANSF + PY1
c              PZ1CM  = BETAZ * TRANSF + PZ1
c              RETURN
c              END
****************************************
