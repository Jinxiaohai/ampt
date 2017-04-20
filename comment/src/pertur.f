      SUBROUTINE pertur(PX,PY,PZ,SRT,IRUN,I1,I2,nt,kp,icont)
*                                                                         *
*       PURPOSE:   TO PRODUCE CASCADE AND OMEGA PERTURBATIVELY            *
c sp 01/03/01
*                   40 cascade-
*                  -40 cascade-(bar)
*                   41 cascade0
*                  -41 cascade0(bar)
*                   45 Omega baryon
*                  -45 Omega baryon(bar)
*                   44 Di-Omega
**********************************
      PARAMETER      (MAXSTR=150001,MAXR=1,PI=3.1415926)
      parameter      (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
      PARAMETER (AMN=0.939457,AMP=0.93828,AP1=0.13496,AP2=0.13957)
      PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974,aks=0.895)
      PARAMETER      (ACAS=1.3213,AOME=1.6724,AMRHO=0.769,AMOMGA=0.782)
      PARAMETER      (AETA=0.548,ADIOMG=3.2288)
      parameter            (maxx=20,maxz=24)
      COMMON   /AA/  R(3,MAXSTR)
cc      SAVE /AA/
      COMMON   /BB/  P(3,MAXSTR)
cc      SAVE /BB/
      COMMON   /CC/  E(MAXSTR)
cc      SAVE /CC/
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
      COMMON   /RR/  MASSR(0:MAXR)
cc      SAVE /RR/
      COMMON   /BG/BETAX,BETAY,BETAZ,GAMMA
cc      SAVE /BG/
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
c     perturbative method is disabled:
c      common /imulst/ iperts
c
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      SAVE   
      px0 = px
      py0 = py
      pz0 = pz
      LB1 = LB(I1)
      EM1 = E(I1)
      X1  = R(1,I1)
      Y1  = R(2,I1)
      Z1  = R(3,I1)
      prob1 = PROPER(I1)
c     
      LB2 = LB(I2)
      EM2 = E(I2)
      X2  = R(1,I2)
      Y2  = R(2,I2)
      Z2  = R(3,I2)
      prob2 = PROPER(I2)
c
c                 !! flag for real 2-body process (1/0=no/yes)
      icont = 1
c                !! flag for elastic scatt only (-1=no)
      icsbel = -1
* K-/K*0bar + La/Si --> cascade + pi
* K+/K*0 + La/Si (bar) --> cascade-bar + pi
       if( (lb1.eq.21.or.lb1.eq.23.or.iabs(lb1).eq.30) .and.
     &     (iabs(lb2).ge.14.and.iabs(lb2).le.17) )go to 60
       if( (lb2.eq.21.or.lb2.eq.23.or.iabs(lb2).eq.30) .and.
     &     (iabs(lb1).ge.14.and.iabs(lb1).le.17) )go to 60
* K-/K*0bar + cascade --> omega + pi
* K+/K*0 + cascade-bar --> omega-bar + pi
        if( (lb1.eq.21.or.lb1.eq.23.or.iabs(lb1).eq.30) .and.
     &      (iabs(lb2).eq.40.or.iabs(lb2).eq.41) )go to 70
        if( (lb2.eq.21.or.lb2.eq.23.or.iabs(lb2).eq.30) .and.
     &      (iabs(lb1).eq.40.or.iabs(lb1).eq.41) )go to 70
c
c annhilation of cascade,cascade-bar, omega,omega-bar
c
* K- + La/Si <-- cascade + pi(eta,rho,omega)
* K+ + La/Si(bar) <-- cascade-bar + pi(eta,rho,omega)
       if( (((lb1.ge.3.and.lb1.le.5).or.lb1.eq.0) 
     &        .and.(iabs(lb2).eq.40.or.iabs(lb2).eq.41))
     & .OR. (((lb2.ge.3.and.lb2.le.5).or.lb2.eq.0) 
     &        .and.(iabs(lb1).eq.40.or.iabs(lb1).eq.41)) )go to 90
* K- + cascade <-- omega + pi
* K+ + cascade-bar <-- omega-bar + pi
c         if( (lb1.eq.0.and.iabs(lb2).eq.45)
c    &    .OR. (lb2.eq.0.and.iabs(lb1).eq.45) ) go to 110
       if( ((lb1.ge.3.and.lb1.le.5).and.iabs(lb2).eq.45)
     & .OR.((lb2.ge.3.and.lb2.le.5).and.iabs(lb1).eq.45) )go to 110
c
c----------------------------------------------------
*  for process:  K-bar + L(S) --> Ca + pi 
*
60         if(iabs(lb1).ge.14 .and. iabs(lb1).le.17)then 
             asap = e(i1)
             akap = e(i2)
             idp = i1
           else
             asap = e(i2)
             akap = e(i1)
             idp = i2
           endif
          app = 0.138
         if(srt .lt. (acas+app))return
          srrt = srt - (acas+app) + (amn+akap)
          pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2 - akap**2)
          sigca = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
clin pii & pff should be each divided by (4*srt**2), 
c     but these two factors cancel out in the ratio pii/pff:
          pii = sqrt((srt**2-(amn+akap)**2)*(srt**2-(amn-akap)**2))
          pff = sqrt((srt**2-(asap+app)**2)*(srt**2-(asap-app)**2))
         cmat = sigca*pii/pff
         sigpi = cmat*
     &            sqrt((srt**2-(acas+app)**2)*(srt**2-(acas-app)**2))/
     &            sqrt((srt**2-(asap+akap)**2)*(srt**2-(asap-akap)**2))
c 
         sigeta = 0.
        if(srt .gt. (acas+aeta))then
           srrt = srt - (acas+aeta) + (amn+akap)
         pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2 - akap**2)
            sigca = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
         cmat = sigca*pii/pff
         sigeta = cmat*
     &            sqrt((srt**2-(acas+aeta)**2)*(srt**2-(acas-aeta)**2))/
     &            sqrt((srt**2-(asap+akap)**2)*(srt**2-(asap-akap)**2))
        endif
c
         sigca = sigpi + sigeta
         sigpe = 0.
clin-2/25/03 disable the perturb option:
c        if(iperts .eq. 1) sigpe = 40.   !! perturbative xsecn
           sig = amax1(sigpe,sigca)     
         ds = sqrt(sig/31.4)
         dsr = ds + 0.1
         ec = (em1+em2+0.02)**2
         call distce(i1,i2,dsr,ds,dt,ec,srt,ic,px,py,pz)
           if(ic .eq. -1)return
          brpp = sigca/sig
c
c else particle production
          if( (lb1.ge.14.and.lb1.le.17) .or.
     &          (lb2.ge.14.and.lb2.le.17) )then
c   !! cascade- or cascde0
            lbpp1 = 40 + int(2*RANART(NSEED))
          else
* elseif(lb1 .eq. -14 .or. lb2 .eq. -14)
c     !! cascade-bar- or cascde0 -bar
            lbpp1 = -40 - int(2*RANART(NSEED))
          endif
              empp1 = acas
           if(RANART(NSEED) .lt. sigpi/sigca)then
c    !! pion
            lbpp2 = 3 + int(3*RANART(NSEED))
            empp2 = 0.138
           else
c    !! eta
            lbpp2 = 0
            empp2 = aeta
           endif        
c* check real process of cascade(bar) and pion formation
          if(RANART(NSEED) .lt. brpp)then
c       !! real process flag
            icont = 0
            lb(i1) = lbpp1
            e(i1) = empp1
c  !! cascade formed with prob Gam
            proper(i1) = brpp
            lb(i2) = lbpp2
            e(i2) = empp2
c         !! pion/eta formed with prob 1.
            proper(i2) = 1.
           endif
c else only cascade(bar) formed perturbatively
             go to 700
c----------------------------------------------------
*  for process:  Cas(bar) + K_bar(K) --> Om(bar) + pi  !! eta
*
70         if(iabs(lb1).eq.40 .or. iabs(lb1).eq.41)then 
             acap = e(i1)
             akap = e(i2)
             idp = i1
           else
             acap = e(i2)
             akap = e(i1)
             idp = i2
           endif
           app = 0.138
*         ames = aeta
c  !! only pion
           ames = 0.138
         if(srt .lt. (aome+ames))return 
          srrt = srt - (aome+ames) + (amn+akap)
         pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2 - akap**2)
c use K(bar) + Ca --> Om + eta  xsecn same as  K(bar) + N --> Si + Pi
*  as Omega have no resonances
c** using same matrix elements as K-bar + N -> Si + pi
         sigomm = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
         cmat = sigomm*
     &          sqrt((srt**2-(amn+akap)**2)*(srt**2-(amn-akap)**2))/
     &          sqrt((srt**2-(asa+app)**2)*(srt**2-(asa-app)**2))
        sigom = cmat*
     &           sqrt((srt**2-(aome+ames)**2)*(srt**2-(aome-ames)**2))/
     &           sqrt((srt**2-(acap+akap)**2)*(srt**2-(acap-akap)**2))
          sigpe = 0.
clin-2/25/03 disable the perturb option:
c         if(iperts .eq. 1) sigpe = 40.   !! perturbative xsecn
          sig = amax1(sigpe,sigom)     
         ds = sqrt(sig/31.4)
         dsr = ds + 0.1
         ec = (em1+em2+0.02)**2
         call distce(i1,i2,dsr,ds,dt,ec,srt,ic,px,py,pz)
           if(ic .eq. -1)return
           brpp = sigom/sig
c
c else particle production
           if( (lb1.ge.40.and.lb1.le.41) .or.
     &           (lb2.ge.40.and.lb2.le.41) )then
c    !! omega
            lbpp1 = 45
           else
* elseif(lb1 .eq. -40 .or. lb2 .eq. -40)
c    !! omega-bar
            lbpp1 = -45
           endif
           empp1 = aome
*           lbpp2 = 0    !! eta
c    !! pion
           lbpp2 = 3 + int(3*RANART(NSEED))
           empp2 = ames
c
c* check real process of omega(bar) and pion formation
           xrand=RANART(NSEED)
         if(xrand .lt. (proper(idp)*brpp))then
c       !! real process flag
            icont = 0
            lb(i1) = lbpp1
            e(i1) = empp1
c  !! P_Om = P_Cas*Gam
            proper(i1) = proper(idp)*brpp
            lb(i2) = lbpp2
            e(i2) = empp2
c   !! pion formed with prob 1.
            proper(i2) = 1.
          elseif(xrand.lt.brpp) then
c else omega(bar) formed perturbatively and cascade destroyed
             e(idp) = 0.
          endif
             go to 700
c-----------------------------------------------------------
*  for process:  Ca + pi/eta --> K-bar + L(S)
*
90         if(iabs(lb1).eq.40 .or. iabs(lb1).eq.41)then 
             acap = e(i1)
             app = e(i2)
             idp = i1
             idn = i2
           else
             acap = e(i2)
             app = e(i1)
             idp = i2
             idn = i1
           endif
c            akal = (aka+aks)/2.  !! average of K and K* taken
c  !! using K only
            akal = aka
c
         alas = ala
       if(srt .le. (alas+aka))return
           srrt = srt - (acap+app) + (amn+aka)
         pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2 - aka**2)
c** using same matrix elements as K-bar + N -> La/Si + pi
         sigca = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
         cmat = sigca*
     &          sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/
     &          sqrt((srt**2-(alas+0.138)**2)*(srt**2-(alas-0.138)**2))
         sigca = cmat*
     &            sqrt((srt**2-(acap+app)**2)*(srt**2-(acap-app)**2))/
     &            sqrt((srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2))
c    !! pi
            dfr = 1./3.
c       !! eta
           if(lb(idn).eq.0)dfr = 1.
        sigcal = sigca*dfr*(srt**2-(alas+aka)**2)*
     &           (srt**2-(alas-aka)**2)/(srt**2-(acap+app)**2)/
     &           (srt**2-(acap-app)**2)
c
          alas = ASA
       if(srt .le. (alas+aka))then
         sigcas = 0.
       else
           srrt = srt - (acap+app) + (amn+aka)
        pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2 - aka**2)
c use K(bar) + La/Si --> Ca + Pi  xsecn same as  K(bar) + N --> Si + Pi
c** using same matrix elements as K-bar + N -> La/Si + pi
          sigca = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
         cmat = sigca*
     &          sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/
     &          sqrt((srt**2-(alas+0.138)**2)*(srt**2-(alas-0.138)**2))
         sigca = cmat*
     &            sqrt((srt**2-(acap+app)**2)*(srt**2-(acap-app)**2))/
     &            sqrt((srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2))
c    !! pi
            dfr = 1.
c    !! eta
           if(lb(idn).eq.0)dfr = 3.
        sigcas = sigca*dfr*(srt**2-(alas+aka)**2)*
     &           (srt**2-(alas-aka)**2)/(srt**2-(acap+app)**2)/
     &           (srt**2-(acap-app)**2)
       endif
c
         sig = sigcal + sigcas
         brpp = 1.                                                   
         ds = sqrt(sig/31.4)
         dsr = ds + 0.1
         ec = (em1+em2+0.02)**2
         call distce(i1,i2,dsr,ds,dt,ec,srt,ic,px,py,pz)
c
clin-2/25/03: checking elastic scatt after failure of inelastic scatt gives 
c     conditional probability (in general incorrect), tell Pal to correct:
       if(ic .eq. -1)then
c check for elastic scatt, no particle annhilation
c  !! elastic cross section of 20 mb
         ds = sqrt(20.0/31.4)
         dsr = ds + 0.1
         call distce(i1,i2,dsr,ds,dt,ec,srt,icsbel,px,py,pz)
           if(icsbel .eq. -1)return
            empp1 = EM1
            empp2 = EM2
             go to 700
       endif
c
c else pert. produced cascade(bar) is annhilated OR real process
c
* DECIDE LAMBDA OR SIGMA PRODUCTION
c
       IF(sigcal/sig .GT. RANART(NSEED))THEN  
          if(lb1.eq.40.or.lb1.eq.41.or.lb2.eq.40.or.lb2.eq.41)then
          lbpp1 = 21
           lbpp2 = 14
          else
           lbpp1 = 23
           lbpp2 = -14
          endif
         alas = ala
       ELSE
          if(lb1.eq.40.or.lb1.eq.41.or.lb2.eq.40.or.lb2.eq.41)then
           lbpp1 = 21
            lbpp2 = 15 + int(3 * RANART(NSEED))
          else
            lbpp1 = 23
            lbpp2 = -15 - int(3 * RANART(NSEED))
          endif
         alas = ASA       
        ENDIF
             empp1 = aka  
             empp2 = alas 
c
c check for real process for L/S(bar) and K(bar) formation
          if(RANART(NSEED) .lt. proper(idp))then
* real process
c       !! real process flag
            icont = 0
            lb(i1) = lbpp1
            e(i1) = empp1
c   !! K(bar) formed with prob 1.
            proper(i1) = 1.
            lb(i2) = lbpp2
            e(i2) = empp2
c   !! L/S(bar) formed with prob 1.
            proper(i2) = 1.
             go to 700
           else
c else only cascade(bar) annhilation & go out
            e(idp) = 0.
           endif
          return
c
c----------------------------------------------------
*  for process:  Om(bar) + pi --> Cas(bar) + K_bar(K)
*
110         if(lb1 .eq. 45 .or. lb1 .eq. -45)then 
             aomp = e(i1)
             app = e(i2)
             idp = i1
             idn = i2
           else
             aomp = e(i2)
             app = e(i1)
             idp = i2
             idn = i1
           endif
c            akal = (aka+aks)/2.  !! average of K and K* taken 
c  !! using K only
            akal = aka
       if(srt .le. (acas+aka))return
           srrt = srt - (aome+app) + (amn+aka)
         pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2 - aka**2)
c use K(bar) + Ca --> Om + eta  xsecn same as  K(bar) + N --> Si + Pi
c** using same matrix elements as K-bar + N -> La/Si + pi
           sigca = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
         cmat = sigca*
     &          sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/
     &          sqrt((srt**2-(asa+0.138)**2)*(srt**2-(asa-0.138)**2))
         sigom = cmat*
     &            sqrt((srt**2-(aomp+app)**2)*(srt**2-(aomp-app)**2))/
     &            sqrt((srt**2-(acas+aka)**2)*(srt**2-(acas-aka)**2))
c            dfr = 2.    !! eta
c    !! pion
           dfr = 2./3.
        sigom = sigom*dfr*(srt**2-(acas+aka)**2)*
     &           (srt**2-(acas-aka)**2)/(srt**2-(aomp+app)**2)/
     &           (srt**2-(aomp-app)**2)
c                                                                         
         brpp = 1.
         ds = sqrt(sigom/31.4)
         dsr = ds + 0.1
         ec = (em1+em2+0.02)**2
         call distce(i1,i2,dsr,ds,dt,ec,srt,ic,px,py,pz)
c
clin-2/25/03: checking elastic scatt after failure of inelastic scatt gives 
c     conditional probability (in general incorrect), tell Pal to correct:
       if(ic .eq. -1)then
c check for elastic scatt, no particle annhilation
c  !! elastic cross section of 20 mb
         ds = sqrt(20.0/31.4)
         dsr = ds + 0.1
         call distce(i1,i2,dsr,ds,dt,ec,srt,icsbel,px,py,pz)
           if(icsbel .eq. -1)return
            empp1 = EM1
            empp2 = EM2
             go to 700
       endif
c
c else pert. produced omega(bar) annhilated  OR real process
c annhilate only pert. omega, rest from hijing go out WITHOUT annhil.
           if(lb1.eq.45 .or. lb2.eq.45)then
c  !! Ca
             lbpp1 = 40 + int(2*RANART(NSEED))
c   !! K-
             lbpp2 = 21
            else
* elseif(lb1 .eq. -45 .or. lb2 .eq. -45)
c    !! Ca-bar
            lbpp1 = -40 - int(2*RANART(NSEED))
c      !! K+
            lbpp2 = 23
           endif
             empp1 = acas
             empp2 = aka  
c
c check for real process for Cas(bar) and K(bar) formation
          if(RANART(NSEED) .lt. proper(idp))then
c       !! real process flag
            icont = 0
            lb(i1) = lbpp1
            e(i1) = empp1
c   !! P_Cas(bar) = P_Om(bar)
            proper(i1) = proper(idp)
            lb(i2) = lbpp2
            e(i2) = empp2
c   !! K(bar) formed with prob 1.
            proper(i2) = 1.
c
           else
c else Cascade(bar)  produced and Omega(bar) annhilated
            e(idp) = 0.
           endif
c   !! for produced particles
             go to 700
c
c-----------------------------------------------------------
700    continue
* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
* ENERGY CONSERVATION
          PR2   = (SRT**2 - EMpp1**2 - EMpp2**2)**2
     &                - 4.0 * (EMpp1*EMpp2)**2
          IF(PR2.LE.0.)PR2=0.00000001
          PR=SQRT(PR2)/(2.*SRT)
* using isotropic
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
       if(icont .eq. 0)return
c
* LORENTZ-TRANSFORMATION INTO CMS FRAME
              E1CM    = SQRT (EMpp1**2 + PX**2 + PY**2 + PZ**2)
              P1BETA  = PX*BETAX + PY*BETAY + PZ*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Ppt11 = BETAX * TRANSF + PX
              Ppt12 = BETAY * TRANSF + PY
              Ppt13 = BETAZ * TRANSF + PZ
c
cc** for elastic scattering update the momentum of pertb particles
         if(icsbel .ne. -1)then
c            if(EMpp1 .gt. 0.9)then
              p(1,i1) = Ppt11
              p(2,i1) = Ppt12
              p(3,i1) = Ppt13
c            else
              E2CM    = SQRT (EMpp2**2 + PX**2 + PY**2 + PZ**2)
              TRANSF  = GAMMA * ( -GAMMA * P1BETA / (GAMMA + 1) + E2CM )
              Ppt21 = BETAX * TRANSF - PX
              Ppt22 = BETAY * TRANSF - PY
              Ppt23 = BETAZ * TRANSF - PZ
              p(1,i2) = Ppt21
              p(2,i2) = Ppt22
              p(3,i2) = Ppt23
c            endif
             return
          endif
clin-5/2008:
c2008        X01 = 1.0 - 2.0 * RANART(NSEED)
c            Y01 = 1.0 - 2.0 * RANART(NSEED)
c            Z01 = 1.0 - 2.0 * RANART(NSEED)
c        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2008
c                Xpt=X1+0.5*x01
c                Ypt=Y1+0.5*y01
c                Zpt=Z1+0.5*z01
                Xpt=X1
                Ypt=Y1
                Zpt=Z1
c
c
c          if(lbpp1 .eq. 45)then
c           write(*,*)'II lb1,lb2,lbpp1,empp1,proper(idp),brpp'
c           write(*,*)lb1,lb2,lbpp1,empp1,proper(idp),brpp
c          endif
c
               NNN=NNN+1
               PROPI(NNN,IRUN)= proper(idp)*brpp
               LPION(NNN,IRUN)= lbpp1
               EPION(NNN,IRUN)= empp1
                RPION(1,NNN,IRUN)=Xpt
                RPION(2,NNN,IRUN)=Ypt
                RPION(3,NNN,IRUN)=Zpt
               PPION(1,NNN,IRUN)=Ppt11
               PPION(2,NNN,IRUN)=Ppt12
               PPION(3,NNN,IRUN)=Ppt13
clin-5/2008:
               dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
            RETURN
            END
**********************************
*  sp 12/08/00                                                         *
