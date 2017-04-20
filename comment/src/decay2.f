        SUBROUTINE DECAY2(IRUN,I,NNN,ISEED,wid,nt)
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
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
        lbanti=LB(I)
c
        DM=E(I)
* DETERMINE THE DECAY PRODUCTS
* FOR N*+(1440) DECAY
        IF(iabs(LB(I)).EQ.11)THEN
           X3=RANART(NSEED)
           IF(X3.LT.(1./3))THEN
              LB(I)=2
              NLAB=2
              LPION(NNN,IRUN)=5
              EPION(NNN,IRUN)=AP2
              LPION(NNN+1,IRUN)=4
              EPION(NNN+1,IRUN)=AP1
           ELSEIF(X3.LT.2./3.AND.X3.GT.1./3.)THEN
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=5
              EPION(NNN,IRUN)=AP2
              LPION(NNN+1,IRUN)=3
              EPION(NNN+1,IRUN)=AP2
           ELSE
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=4
              EPION(NNN,IRUN)=AP1
              LPION(NNN+1,IRUN)=4
              EPION(NNN+1,IRUN)=AP1
           ENDIF
* FOR N*0(1440) DECAY
        ELSEIF(iabs(LB(I)).EQ.10)THEN
           X3=RANART(NSEED)
           IF(X3.LT.(1./3))THEN
              LB(I)=2
              NLAB=2
              LPION(NNN,IRUN)=4
              EPION(NNN,IRUN)=AP1
              LPION(NNN+1,IRUN)=4
              EPION(NNN+1,IRUN)=AP1
           ELSEIF(X3.LT.2./3.AND.X3.GT.1./3.)THEN
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=3
              EPION(NNN,IRUN)=AP2
              LPION(NNN+1,IRUN)=4
              EPION(NNN+1,IRUN)=AP1
           ELSE
              LB(I)=2
              NLAB=2
              LPION(NNN,IRUN)=5
              EPION(NNN,IRUN)=AP2
              LPION(NNN+1,IRUN)=3
              EPION(NNN+1,IRUN)=AP2
           ENDIF
        ENDIF
        CALL DKINE2(IRUN,I,NNN,NLAB,ISEED,wid,nt)
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
c
           lbi=LPION(NNN+1,IRUN)
           if(lbi.eq.3) then
              lbi=5
           elseif(lbi.eq.5) then
              lbi=3
           elseif(lbi.eq.1.or.lbi.eq.2) then
              lbi=-lbi
           endif
           LPION(NNN+1,IRUN)=lbi
        endif
c
       RETURN
       END
*-------------------------------------------------------------------
*--------------------------------------------------------------------------
*         CALCULATE THE MOMENTUM OF NUCLEON AND PION (OR ETA) 
*         IN THE LAB. FRAME AFTER DELTA OR N* DECAY
* DATE   : JAN. 24,1990, MODIFIED ON MAY 17, 1994 TO INCLUDE ETA PRODUCTION
*--------------------------------------------------------------------------
