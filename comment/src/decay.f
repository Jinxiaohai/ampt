        SUBROUTINE DECAY(IRUN,I,NNN,ISEED,wid,nt)
        PARAMETER (MAXSTR=150001,MAXR=1,
     1  AMN=0.939457,ETAM=0.5475,AMP=0.93828,AP1=0.13496,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926)
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
        COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &       IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
cc      SAVE /INPUT2/
      COMMON/RNDF77/NSEED
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
cc      SAVE /RNDF77/
      SAVE   
        lbanti=LB(I)
c
        DM=E(I)
*1. FOR N*+(1440) DECAY
        IF(iabs(LB(I)).EQ.11)THEN
           X3=RANART(NSEED)
           IF(X3.GT.(1./3.))THEN
              LB(I)=2
              NLAB=2
              LPION(NNN,IRUN)=5
              EPION(NNN,IRUN)=AP2
           ELSE
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=4
              EPION(NNN,IRUN)=AP1
           ENDIF
*2. FOR N*0(1440) DECAY
        ELSEIF(iabs(LB(I)).EQ.10)THEN
           X4=RANART(NSEED)
           IF(X4.GT.(1./3.))THEN
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=3
              EPION(NNN,IRUN)=AP2
           ELSE
              LB(I)=2
              NALB=2
              LPION(NNN,IRUN)=4
              EPION(NNN,IRUN)=AP1
           ENDIF
* N*(1535) CAN DECAY TO A PION OR AN ETA IF DM > 1.49 GeV
*3 N*(0)(1535) DECAY
        ELSEIF(iabs(LB(I)).EQ.12)THEN
           CTRL=0.65
           IF(DM.lE.1.49)ctrl=-1.
           X5=RANART(NSEED)
           IF(X5.GE.ctrl)THEN
* DECAY TO PION+NUCLEON
              X6=RANART(NSEED)
              IF(X6.GT.(1./3.))THEN
                 LB(I)=1
                 NLAB=1
                 LPION(NNN,IRUN)=3
                 EPION(NNN,IRUN)=AP2
              ELSE
                 LB(I)=2
                 NALB=2
                 LPION(NNN,IRUN)=4
                 EPION(NNN,IRUN)=AP1
              ENDIF
           ELSE
* DECAY TO ETA+NEUTRON
              LB(I)=2
              NLAB=2
              LPION(NNN,IRUN)=0
              EPION(NNN,IRUN)=ETAM
           ENDIF
*4. FOR N*+(1535) DECAY
        ELSEIF(iabs(LB(I)).EQ.13)THEN
           CTRL=0.65
           IF(DM.lE.1.49)ctrl=-1.
           X5=RANART(NSEED)
           IF(X5.GE.ctrl)THEN
* DECAY TO PION+NUCLEON
              X8=RANART(NSEED)
              IF(X8.GT.(1./3.))THEN
                 LB(I)=2
                 NLAB=2
                 LPION(NNN,IRUN)=5
                 EPION(NNN,IRUN)=AP2
              ELSE
                 LB(I)=1
                 NLAB=1
                 LPION(NNN,IRUN)=4
                 EPION(NNN,IRUN)=AP1
              ENDIF
           ELSE
* DECAY TO ETA+NUCLEON
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=0
              EPION(NNN,IRUN)=ETAM
           ENDIF
        ENDIF
c
        CALL DKINE(IRUN,I,NNN,NLAB,ISEED,wid,nt)
c
c     anti-particle ID for anti-N* decays:
        if(lbanti.lt.0) then
           lbi=LB(I)
           if(lbi.eq.1.or.lbi.eq.2) then
              lbi=-lbi
           elseif(lbi.eq.3) then
              lbi=5
           elseif(lbi.eq.5) then
              lbi=3
           endif
           LB(I)=lbi
c
           lbi=LPION(NNN,IRUN)
           if(lbi.eq.3) then
              lbi=5
           elseif(lbi.eq.5) then
              lbi=3
           elseif(lbi.eq.1.or.lbi.eq.2) then
              lbi=-lbi
           endif
           LPION(NNN,IRUN)=lbi
        endif
c
        if(nt.eq.ntmax) then
c     at the last timestep, assign rho or eta (decay daughter) 
c     to lb(i1) only (not to lpion) in order to decay them again:
           lbm=LPION(NNN,IRUN)
           if(lbm.eq.0.or.lbm.eq.25
     1          .or.lbm.eq.26.or.lbm.eq.27) then
c     switch rho or eta with baryon, positions are the same (no change needed):
              lbsave=lbm
              xmsave=EPION(NNN,IRUN)
              pxsave=PPION(1,NNN,IRUN)
              pysave=PPION(2,NNN,IRUN)
              pzsave=PPION(3,NNN,IRUN)
clin-5/2008:
              dpsave=dppion(NNN,IRUN)
              LPION(NNN,IRUN)=LB(I)
              EPION(NNN,IRUN)=E(I)
              PPION(1,NNN,IRUN)=P(1,I)
              PPION(2,NNN,IRUN)=P(2,I)
              PPION(3,NNN,IRUN)=P(3,I)
clin-5/2008:
              dppion(NNN,IRUN)=dpertp(I)
              LB(I)=lbsave
              E(I)=xmsave
              P(1,I)=pxsave
              P(2,I)=pysave
              P(3,I)=pzsave
clin-5/2008:
              dpertp(I)=dpsave
           endif
        endif
       RETURN
       END
*-------------------------------------------------------------------
*-------------------------------------------------------------------
* PURPOSE:
*         CALCULATE THE MOMENTUM OF NUCLEON AND PION (OR ETA) 
*         IN THE LAB. FRAME AFTER DELTA OR N* DECAY
* DATE   : JAN. 24,1990, MODIFIED ON MAY 17, 1994 TO INCLUDE ETA PRODUCTION
