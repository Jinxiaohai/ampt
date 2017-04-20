        SUBROUTINE DKINE(IRUN,I,NNN,NLAB,ISEED,wid,nt)
        PARAMETER (hbarc=0.19733)
        PARAMETER (MAXSTR=150001,MAXR=1,
     1  AMN=0.939457,AMP=0.93828,ETAM=0.5475,
     2  AP1=0.13496,AP2=0.13957,AM0=1.232,PI=3.1415926)
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
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1 px1n,py1n,pz1n,dp1n
cc      SAVE /leadng/
        COMMON/tdecay/tfdcy(MAXSTR),tfdpi(MAXSTR,MAXR),tft(MAXSTR)
cc      SAVE /tdecay/
        COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &       IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
cc      SAVE /INPUT2/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
        EXTERNAL IARFLV, INVFLV
      SAVE   
* READ IN THE COORDINATES OF DELTA OR N* UNDERGOING DECAY
        PX=P(1,I)
        PY=P(2,I)
        PZ=P(3,I)
        RX=R(1,I)
        RY=R(2,I)
        RZ=R(3,I)
        DM=E(I)
        EDELTA=SQRT(DM**2+PX**2+PY**2+PZ**2)
        PM=EPION(NNN,IRUN)
        AM=AMP
        IF(NLAB.EQ.2)AM=AMN
* FIND OUT THE MOMENTUM AND ENERGY OF PION AND NUCLEON IN DELTA REST FRAME
* THE MAGNITUDE OF MOMENTUM IS DETERMINED BY ENERGY CONSERVATION ,THE FORMULA
* CAN BE FOUND ON PAGE 716,W BAUER P.R.C40,1989
* THE DIRECTION OF THE MOMENTUM IS ASSUMED ISOTROPIC. NOTE THAT P(PION)=-P(N)
        Q2=((DM**2-AM**2+PM**2)/(2.*DM))**2-PM**2
        IF(Q2.LE.0.)Q2=1.e-09
        Q=SQRT(Q2)
11      QX=1.-2.*RANART(NSEED)
        QY=1.-2.*RANART(NSEED)
        QZ=1.-2.*RANART(NSEED)
        QS=QX**2+QY**2+QZ**2
        IF(QS.GT.1.) GO TO 11
        PXP=Q*QX/SQRT(QS)
        PYP=Q*QY/SQRT(QS)
        PZP=Q*QZ/SQRT(QS)
        EP=SQRT(Q**2+PM**2)
        PXN=-PXP
        PYN=-PYP
        PZN=-PZP
        EN=SQRT(Q**2+AM**2)
* TRANSFORM INTO THE LAB. FRAME. THE GENERAL LORENTZ TRANSFORMATION CAN
* BE FOUND ON PAGE 34 OF R. HAGEDORN " RELATIVISTIC KINEMATICS"
        GD=EDELTA/DM
        FGD=GD/(1.+GD)
        BDX=PX/EDELTA
        BDY=PY/EDELTA
        BDZ=PZ/EDELTA
        BPP=BDX*PXP+BDY*PYP+BDZ*PZP
        BPN=BDX*PXN+BDY*PYN+BDZ*PZN
        P(1,I)=PXN+BDX*GD*(FGD*BPN+EN)
        P(2,I)=PYN+BDY*GD*(FGD*BPN+EN)
        P(3,I)=PZN+BDZ*GD*(FGD*BPN+EN)
        E(I)=AM
* WE ASSUME THAT THE SPACIAL COORDINATE OF THE NUCLEON
* IS THAT OF THE DELTA
        PPION(1,NNN,IRUN)=PXP+BDX*GD*(FGD*BPP+EP)
        PPION(2,NNN,IRUN)=PYP+BDY*GD*(FGD*BPP+EP)
        PPION(3,NNN,IRUN)=PZP+BDZ*GD*(FGD*BPP+EP)
clin-5/2008:
        dppion(NNN,IRUN)=dpertp(I)
* WE ASSUME THE PION OR ETA COMING FROM DELTA DECAY IS LOCATED ON THE SPHERE
* OF RADIUS 0.5FM AROUND DELTA, THIS POINT NEED TO BE CHECKED 
* AND OTHER CRIERTION MAY BE TRIED
clin-2/20/03 no additional smearing for position of decay daughters:
c200         X0 = 1.0 - 2.0 * RANART(NSEED)
c            Y0 = 1.0 - 2.0 * RANART(NSEED)
c            Z0 = 1.0 - 2.0 * RANART(NSEED)
c        IF ((X0*X0+Y0*Y0+Z0*Z0) .GT. 1.0) GOTO 200
c        RPION(1,NNN,IRUN)=R(1,I)+0.5*x0
c        RPION(2,NNN,IRUN)=R(2,I)+0.5*y0
c        RPION(3,NNN,IRUN)=R(3,I)+0.5*z0
        RPION(1,NNN,IRUN)=R(1,I)
        RPION(2,NNN,IRUN)=R(2,I)
        RPION(3,NNN,IRUN)=R(3,I)
c
        devio=SQRT(EPION(NNN,IRUN)**2+PPION(1,NNN,IRUN)**2
     1       +PPION(2,NNN,IRUN)**2+PPION(3,NNN,IRUN)**2)
     2       +SQRT(E(I)**2+P(1,I)**2+P(2,I)**2+P(3,I)**2)-e1
c        if(abs(devio).gt.0.02) write(93,*) 'decay(): nt=',nt,devio,lb1
c     add decay time to daughter's formation time at the last timestep:
        if(nt.eq.ntmax) then
           tau0=hbarc/wid
           taudcy=tau0*(-1.)*alog(1.-RANART(NSEED))
c     lorentz boost:
           taudcy=taudcy*e1/em1
           tfnl=tfnl+taudcy
           xfnl=xfnl+px1/e1*taudcy
           yfnl=yfnl+py1/e1*taudcy
           zfnl=zfnl+pz1/e1*taudcy
           R(1,I)=xfnl
           R(2,I)=yfnl
           R(3,I)=zfnl
           tfdcy(I)=tfnl
           RPION(1,NNN,IRUN)=xfnl
           RPION(2,NNN,IRUN)=yfnl
           RPION(3,NNN,IRUN)=zfnl
           tfdpi(NNN,IRUN)=tfnl
        endif
 200    format(a30,2(1x,e10.4))
 210    format(i6,5(1x,f8.3))
 220    format(a2,i5,5(1x,f8.3))
        RETURN
        END
*-----------------------------------------------------------------------------
*-----------------------------------------------------------------------------
* PURPOSE:1. N*-->N+PION+PION  DECAY PRODUCTS
*         2. DETERMINE THE MOMENTUM AND COORDINATES OF NUCLEON AND PION
*            AFTER THE DELTA OR N* DECAYING
* DATE   : NOV.7,1994
*----------------------------------------------------------------------------
