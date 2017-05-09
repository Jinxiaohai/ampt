       SUBROUTINE ARTMN
cbz11/16/98end
**********************************
* PARAMETERS:                                                           *
*  MAXPAR     - MAXIMUM NUMBER OF PARTICLES      PROGRAM CAN HANDLE     *
*  MAXP       - MAXIMUM NUMBER OF CREATED MESONS PROGRAM CAN HANDLE     *
*  MAXR       - MAXIMUM NUMBER OF EVENTS AT EACH IMPACT PARAMETER       *
*  MAXX       - NUMBER OF MESHPOINTS IN X AND Y DIRECTION = 2 MAXX + 1  *
*  MAXZ       - NUMBER OF MESHPOINTS IN Z DIRECTION       = 2 MAXZ + 1  *
*  AMU        - 1 ATOMIC MASS UNIT "GEV/C**2"                           *
*  MX,MY,MZ   - MESH SIZES IN COORDINATE SPACE [FM] FOR PAULI LATTICE   *
*  MPX,MPY,MPZ- MESH SIZES IN MOMENTUM SPACE [GEV/C] FOR PAULI LATTICE  *
*---------------------------------------------------------------------- *
clin      PARAMETER     (maxpar=200000,MAXR=50,AMU= 0.9383,
      PARAMETER     (MAXSTR=150001,MAXR=1,AMU= 0.9383,
     1               AKA=0.498,etaM=0.5475)
      PARAMETER     (MAXX   =   20,  MAXZ  =    24)
      PARAMETER     (ISUM   =   1001,  IGAM  =    1100)
      parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
clin      PARAMETER (MAXP = 14000)
*----------------------------------------------------------------------*
      INTEGER   OUTPAR, zta,zpr
      COMMON  /AA/      R(3,MAXSTR)
cc      SAVE /AA/
      COMMON  /BB/      P(3,MAXSTR)
cc      SAVE /BB/
      COMMON  /CC/      E(MAXSTR)
cc      SAVE /CC/
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
cc      SAVE /DD/
      COMMON  /EE/      ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      COMMON  /HH/  PROPER(MAXSTR)
cc      SAVE /HH/
      common  /ff/f(-mx:mx,-my:my,-mz:mz,-mpx:mpx,-mpy:mpy,-mpz:mpzp)
cc      SAVE /ff/
      common  /gg/      dx,dy,dz,dpx,dpy,dpz
cc      SAVE /gg/
      COMMON  /INPUT/ NSTAR,NDIRCT,DIR
cc      SAVE /INPUT/
      COMMON  /PP/      PRHO(-20:20,-24:24)
      COMMON  /QQ/      PHRHO(-MAXZ:MAXZ,-24:24)
      COMMON  /RR/      MASSR(0:MAXR)
cc      SAVE /RR/
      common  /ss/      inout(20)
cc      SAVE /ss/
      common  /zz/      zta,zpr
cc      SAVE /zz/
      COMMON  /RUN/     NUM
cc      SAVE /RUN/
clin-4/2008:
c      COMMON  /KKK/     TKAON(7),EKAON(7,0:200)
      COMMON  /KKK/     TKAON(7),EKAON(7,0:2000)
cc      SAVE /KKK/
      COMMON  /KAON/    AK(3,50,36),SPECK(50,36,7),MF
cc      SAVE /KAON/
      COMMON/TABLE/ xarray(0:1000),earray(0:1000)
cc      SAVE /TABLE/
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      COMMON  /DDpi/    piRHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
cc      SAVE /DDpi/
      common  /tt/  PEL(-maxx:maxx,-maxx:maxx,-maxz:maxz)
     &,rxy(-maxx:maxx,-maxx:maxx,-maxz:maxz)
cc      SAVE /tt/
clin-4/2008:
c      DIMENSION TEMP(3,MAXSTR),SKAON(7),SEKAON(7,0:200)
      DIMENSION TEMP(3,MAXSTR),SKAON(7),SEKAON(7,0:2000)
cbz12/2/98
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
cc      SAVE /INPUT2/
      COMMON /INPUT3/ PLAB, ELAB, ZEROPT, B0, BI, BM, DENCUT, CYCBOX
cc      SAVE /INPUT3/
cbz12/2/98end
cbz11/16/98
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
cc      SAVE /ARPRNT/
c.....note in the below, since a common block in ART is called EE,
c.....the variable EE in /ARPRC/is changed to PEAR.
clin-9/29/03 changed name in order to distinguish from /prec2/
c        COMMON /ARPRC/ ITYPAR(MAXSTR),
c     &       GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
c     &       PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
c     &       XMAR(MAXSTR)
cc      SAVE /ARPRC/
clin-9/29/03-end
      COMMON /ARERCP/PRO1(MAXSTR, MAXR)
cc      SAVE /ARERCP/
      COMMON /ARERC1/MULTI1(MAXR)
cc      SAVE /ARERC1/
      COMMON /ARPRC1/ITYP1(MAXSTR, MAXR),
     &     GX1(MAXSTR, MAXR), GY1(MAXSTR, MAXR), GZ1(MAXSTR, MAXR), 
     &     FT1(MAXSTR, MAXR),
     &     PX1(MAXSTR, MAXR), PY1(MAXSTR, MAXR), PZ1(MAXSTR, MAXR),
     &     EE1(MAXSTR, MAXR), XM1(MAXSTR, MAXR)
cc      SAVE /ARPRC1/
c
      DIMENSION NPI(MAXR)
      DIMENSION RT(3, MAXSTR, MAXR), PT(3, MAXSTR, MAXR)
     &     , ET(MAXSTR, MAXR), LT(MAXSTR, MAXR), PROT(MAXSTR, MAXR)
      EXTERNAL IARFLV, INVFLV
cbz11/16/98end
      common /lastt/itimeh,bimp 
cc      SAVE /lastt/
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
cc      SAVE /snn/
      COMMON/hbt/lblast(MAXSTR),xlast(4,MAXSTR),plast(4,MAXSTR),nlast
cc      SAVE /hbt/
      common/resdcy/NSAV,iksdcy
cc      SAVE /resdcy/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      COMMON/FTMAX/ftsv(MAXSTR),ftsvt(MAXSTR, MAXR)
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
clin-4/2008 zet() expanded to avoid out-of-bound errors:
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
      nlast=0



      do 1002 i=1,MAXSTR
         ftsv(i)=0.
         do 1101 irun=1,maxr
            ftsvt(i,irun)=0.
 1101    continue
         lblast(i)=999
         do 1001 j=1,4
clin-4/2008 bugs pointed out by Vander Molen & Westfall:
c            xlast(i,j)=0.
c            plast(i,j)=0.
            xlast(j,i)=0.
            plast(j,i)=0.
 1001    continue
 1002 continue


      
*-------------------------------------------------------------------*
* Input information about the reaction system and contral parameters* 
*-------------------------------------------------------------------*
*              input section starts here                           *
*-------------------------------------------------------------------*
cbz12/2/98
c.....input section is moved to subroutine ARTSET
cbz12/2/98end
*-----------------------------------------------------------------------*
*                   input section ends here                            *
*-----------------------------------------------------------------------*
* read in the table for gengrating the transverse momentum
* IN THE NN-->DDP PROCESS
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  tablem
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
       call tablem
* several control parameters, keep them fixed in this code. 


       ikaon=1
       nstar=1
       ndirct=0
       dir=0.02
       asy=0.032
       ESBIN=0.04
       MF=36
*----------------------------------------------------------------------*
c      CALL FRONT(12,MASSTA,MASSPR,ELAB)
*----------------------------------------------------------------------*
      RADTA  = 1.124 * FLOAT(MASSTA)**(1./3.)
      RADPR  = 1.124 * FLOAT(MASSPR)**(1./3.)
      ZDIST  = RADTA + RADPR
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      write(9897,*)" RADTA     RADPR     ZDIST"
      write(9897,*)RADTA, RADPR, ZDIST
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c      if ( cycbox.ne.0 ) zdist=0
      BMAX   = RADTA + RADPR
      MASS   = MASSTA + MASSPR
      NTOTAL = NUM * MASS
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      write(9897,*)"BMAX      MASS       NTOTAL"
      write(9897,*)bmax, mass, ntotal
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      
*
      IF (NTOTAL .GT. MAXSTR) THEN
        WRITE(12,'(//10X,''**** FATAL ERROR: TOO MANY TEST PART. ****'//
     & ' '')')
        STOP
      END IF
*
*-----------------------------------------------------------------------
*       RELATIVISTIC KINEMATICS
*
*       1) LABSYSTEM
*
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    实验室系的处理.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ETA    = FLOAT(MASSTA) * AMU
      PZTA   = 0.0
      BETATA = 0.0
      GAMMTA = 1.0      
*
      EPR    = FLOAT(MASSPR) * (AMU + 0.001 * ELAB)
      PZPR   = SQRT( EPR**2 - (AMU * FLOAT(MASSPR))**2 )
      BETAPR = PZPR / EPR
      GAMMPR = 1.0 / SQRT( 1.0 - BETAPR**2 )
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      write(9896,*)"ETA     PZTA    BETATA    GAMMTA"
      write(9896,*)ETA, PZTA, BETATA, GAMMTA
      write(9896,*)"EPR     PZPR    BETAPR    GAMMPR"
      write(9896,*)EPR, PZPR, BETAPR, GAMMPR
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      
*
* BETAC AND GAMMAC OF THE C.M. OBSERVED IN THE LAB. FRAME
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    质心系的速度和gamma
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        BETAC=(PZPR+PZTA)/(EPR+ETA)
        GAMMC=1.0 / SQRT(1.-BETAC**2)
*
c      WRITE(12,'(/10x,''****    KINEMATICAL PARAMETERS    ****''/)')
c      WRITE(12,'(10x,''1) LAB-FRAME:        TARGET PROJECTILE'')')
c      WRITE(12,'(10x,''   ETOTAL "GEV" '',2F11.4)') ETA, EPR
c      WRITE(12,'(10x,''   P "GEV/C"    '',2F11.4)') PZTA, PZPR
c      WRITE(12,'(10x,''   BETA         '',2F11.4)') BETATA, BETAPR
c      WRITE(12,'(10x,''   GAMMA        '',2F11.4)') GAMMTA, GAMMPR
      IF (INSYS .NE. 0) THEN
*
*       2) C.M. SYSTEM
*
        S      = (EPR+ETA)**2 - PZPR**2
        xx1=4.*alog(float(massta))
        xx2=4.*alog(float(masspr))
        xx1=exp(xx1)
        xx2=exp(xx2)
        PSQARE = (S**2 + (xx1+ xx2) * AMU**4
     &             - 2.0 * S * AMU**2 * FLOAT(MASSTA**2 + MASSPR**2)
     &             - 2.0 * FLOAT(MASSTA**2 * MASSPR**2) * AMU**4)
     &           / (4.0 * S)
*
        ETA    = SQRT ( PSQARE + (FLOAT(MASSTA) * AMU)**2 )
        PZTA   = - SQRT(PSQARE)
        BETATA = PZTA / ETA
        GAMMTA = 1.0 / SQRT( 1.0 - BETATA**2 )
*
        EPR    = SQRT ( PSQARE + (FLOAT(MASSPR) * AMU)**2 )
        PZPR   = SQRT(PSQARE)
        BETAPR = PZPR/ EPR
        GAMMPR = 1.0 / SQRT( 1.0 - BETAPR**2 )
*
c        WRITE(12,'(10x,''2) C.M.-FRAME:  '')')
c        WRITE(12,'(10x,''   ETOTAL "GEV" '',2F11.4)') ETA, EPR
c        WRITE(12,'(10x,''   P "GEV/C"    '',2F11.4)') PZTA, PZPR
c        WRITE(12,'(10x,''   BETA         '',2F11.4)') BETATA, BETAPR
c        WRITE(12,'(10x,''   GAMMA        '',2F11.4)') GAMMTA, GAMMPR
c        WRITE(12,'(10x,''S "GEV**2"      '',F11.4)')  S
c        WRITE(12,'(10x,''PSQARE "GEV/C"2 '',E14.3)')  PSQARE
c        WRITE(12,'(/10x,''*** CALCULATION DONE IN CM-FRAME ***''/)')
      ELSE
c        WRITE(12,'(/10x,''*** CALCULATION DONE IN LAB-FRAME ***''/)')
      END IF
* MOMENTUM PER PARTICLE
      PZTA = PZTA / FLOAT(MASSTA)
      PZPR = PZPR / FLOAT(MASSPR)
* total initial energy in the N-N cms frame
      ECMS0=ETA+EPR
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      write(9896,*)"PZTA    PZPR     ECMS0"
      write(9896,*)pzta, pzpr, ecms0
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<

*-----------------------------------------------------------------------
*
* Start loop over many runs of different impact parameters
* IF MANYB=1, RUN AT A FIXED IMPACT PARAMETER B0, OTHERWISE GENERATE 
* MINIMUM BIAS EVENTS WITHIN THE IMPACT PARAMETER RANGE OF B_MIN AND B_MAX
       DO 50000 IMANY=1,MANYB
*------------------------------------------------------------------------
* Initialize the impact parameter B
       if (manyb. gt.1) then
 111      BX=1.0-2.0*RANART(NSEED)
          BY=1.0-2.0*RANART(NSEED)
          B2=BX*BX+BY*BY
          IF(B2.GT.1.0) GO TO 111       
          B=SQRT(B2)*(BM-BI)+BI
       ELSE
          B=B0
       ENDIF
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
       write(9895,*) "B2     B"
       write(9895,*)B2, B
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
       
c      WRITE(12,'(///10X,''RUN NUMBER:'',I6)') IMANY       
c      WRITE(12,'(//10X,''IMPACT PARAMETER B FOR THIS RUN:'',
c     &             F9.3,'' FM''/10X,49(''*'')/)') B
*
*-----------------------------------------------------------------------
*       INITIALIZATION
*1 INITIALIZATION IN ISOSPIN SPACE FOR BOTH THE PROJECTILE AND TARGET
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  这个初始化是干嘛的?????
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      call coulin(masspr,massta,NUM)
*2 INITIALIZATION IN PHASE SPACE FOR THE TARGET
      CALL INIT(1       ,MASSTA   ,NUM     ,RADTA,
     &          B/2.    ,ZEROPT+ZDIST/2.   ,PZTA,
     &          GAMMTA  ,ISEED    ,MASS    ,IMOMEN)
*3.1 INITIALIZATION IN PHASE SPACE FOR THE PROJECTILE
      CALL INIT(1+MASSTA,MASS     ,NUM     ,RADPR,
     &          -B/2.   ,ZEROPT-ZDIST/2.   ,PZPR,
     &          GAMMPR  ,ISEED    ,MASS    ,IMOMEN)
*3.2 OUTPAR IS THE NO. OF ESCAPED PARTICLES
      OUTPAR = 0
*3.3 INITIALIZATION FOR THE NO. OF PARTICLES IN EACH SAMPLE
*    THIS IS NEEDED DUE TO THE FACT THAT PIONS CAN BE PRODUCED OR ABSORBED
      MASSR(0)=0
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    num的数值为1.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DO 1003 IR =1,NUM
         MASSR(IR)=MASS
 1003 CONTINUE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      write(9894,*)"num    ", num
      do 464 IR = 1, NUM
         write(9894,*)massr(IR)
 464     continue
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      
*3.4 INITIALIZation FOR THE KAON SPECTRUM
*      CALL KSPEC0(BETAC,GAMMC)
*     calculate the local baryon density matrix
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  计算局部的重子数密度.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      CALL DENS(IPOT,MASS,NUM,OUTPAR)
*
*-----------------------------------------------------------------------
*       CONTROL PRINTOUT OF INITIAL CONFIGURATION
*
*      WRITE(12,'(''**********  INITIAL CONFIGURATION  **********''/)')
*
c print out the INITIAL density matrix in the reaction plane
c       do ix=-10,10
c       do iz=-10,10
c       write(1053,992)ix,iz,rho(ix,0,iz)/0.168
c       end do
c       end do
*-----------------------------------------------------------------------
*       CALCULATE MOMENTA FOR T = 0.5 * DT 
*       (TO OBTAIN 2ND DEGREE ACCURACY!)
*       "Reference: J. AICHELIN ET AL., PHYS. REV. C31, 1730 (1985)"
*
      IF (ICOLL .NE. -1) THEN
        DO 700 I = 1,NTOTAL
          IX = NINT( R(1,I) )
          IY = NINT( R(2,I) )
          IZ = NINT( R(3,I) )
clin-4/2008 check bounds:
          IF(IX.GE.MAXX.OR.IY.GE.MAXX.OR.IZ.GE.MAXZ
     1         .OR.IX.LE.-MAXX.OR.IY.LE.-MAXX.OR.IZ.LE.-MAXZ) goto 700
          CALL GRADU(IPOT,IX,IY,IZ,GRADX,GRADY,GRADZ)
          P(1,I) = P(1,I) - (0.5 * DT) * GRADX
          P(2,I) = P(2,I) - (0.5 * DT) * GRADY
          P(3,I) = P(3,I) - (0.5 * DT) * GRADZ
  700   CONTINUE
      END IF
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*4 INITIALIZATION OF TIME-LOOP VARIABLES
*4.1 COLLISION NUMBER COUNTERS
clin 51      RCNNE  = 0
        RCNNE  = 0
       RDD  = 0
       RPP  = 0
       rppk = 0
       RPN  = 0
       rpd  = 0
       RKN  = 0
       RNNK = 0
       RDDK = 0
       RNDK = 0
      RCNND  = 0
      RCNDN  = 0
      RCOLL  = 0
      RBLOC  = 0
      RDIRT  = 0
      RDECAY = 0
      RRES   = 0
*4.11 KAON PRODUCTION PROBABILITY COUNTER FOR PERTURBATIVE CALCULATIONS ONLY
      DO 1005 KKK=1,5
         SKAON(KKK)  = 0
         DO 1004 IS=1,2000
            SEKAON(KKK,IS)=0
 1004    CONTINUE
 1005 CONTINUE
*4.12 anti-proton and anti-kaon counters
       pr0=0.
       pr1=0.
       ska0=0.
       ska1=0.
*       ============== LOOP OVER ALL TIME STEPS ================       *
*                             STARTS HERE                              *
*       ========================================================       *
cbz11/16/98
      IF (IAPAR2(1) .NE. 1) THEN
         DO 1016 I = 1, MAXSTR
            DO 1015 J = 1, 3
               R(J, I) = 0.
               P(J, I) = 0.
 1015       CONTINUE
            E(I) = 0.
            LB(I) = 0
cbz3/25/00
            ID(I)=0
c     sp 12/19/00
           PROPER(I) = 1.
 1016   CONTINUE
         MASS = 0
cbz12/22/98
c         MASSR(1) = 0
c         NP = 0
c         NPI = 1
         NP = 0
         DO 1017 J = 1, NUM
            MASSR(J) = 0
            NPI(J) = 1
 1017    CONTINUE
         DO 1019 I = 1, MAXR
            DO 1018 J = 1, MAXSTR
               RT(1, J, I) = 0.
               RT(2, J, I) = 0.
               RT(3, J, I) = 0.
               PT(1, J, I) = 0.
               PT(2, J, I) = 0.
               PT(3, J, I) = 0.
               ET(J, I) = 0.
               LT(J, I) = 0
c     sp 12/19/00
               PROT(J, I) = 1.
 1018       CONTINUE
 1019    CONTINUE
cbz12/22/98end
      END IF
cbz11/16/98end
      DO 10000 NT = 1,NTMAX
*TEMPORARY PARTICLE COUNTERS
*4.2 PION COUNTERS : LP1,LP2 AND LP3 ARE THE NO. OF P+,P0 AND P-
      LP1=0
      LP2=0
      LP3=0
*4.3 DELTA COUNTERS : LD1,LD2,LD3 AND LD4 ARE THE NO. OF D++,D+,D0 AND D-
      LD1=0
      LD2=0
      LD3=0
      LD4=0
*4.4 N*(1440) COUNTERS : LN1 AND LN2 ARE THE NO. OF N*+ AND N*0
      LN1=0
      LN2=0
*4.5 N*(1535) counters
      LN5=0
*4.6 ETA COUNTERS
      LE=0
*4.7 KAON COUNTERS
      LKAON=0
clin-11/09/00:
* KAON* COUNTERS
      LKAONS=0
*-----------------------------------------------------------------------
        IF (ICOLL .NE. 1) THEN
* STUDYING BINARY COLLISIONS AMONG PARTICLES DURING THIS TIME INTERVAL *
clin-10/25/02 get rid of argument usage mismatch in relcol(.nt.):
           numnt=nt
          CALL RELCOL(LCOLL,LBLOC,LCNNE,LDD,LPP,lppk,
     &    LPN,lpd,LRHO,LOMEGA,LKN,LNNK,LDDK,LNDK,LCNND,
     &    LCNDN,LDIRT,LDECAY,LRES,LDOU,LDDRHO,LNNRHO,
     &    LNNOM,numnt,ntmax,sp,akaon,sk)
c     &    LNNOM,NT,ntmax,sp,akaon,sk)
clin-10/25/02-end
*-----------------------------------------------------------------------
c dilepton production from Dalitz decay
c of pi0 at final time
*      if(nt .eq. ntmax) call dalitz_pi(nt,ntmax)
*                                                                      *
**********************************
*                Lables of collision channels                             *
**********************************
*         LCOLL   - NUMBER OF COLLISIONS              (INTEGER,OUTPUT) *
*         LBLOC   - NUMBER OF PULI-BLOCKED COLLISIONS (INTEGER,OUTPUT) *
*         LCNNE   - NUMBER OF ELASTIC COLLISION       (INTEGER,OUTPUT) *
*         LCNND   - NUMBER OF N+N->N+DELTA REACTION   (INTEGER,OUTPUT) *
*         LCNDN   - NUMBER OF N+DELTA->N+N REACTION   (INTEGER,OUTPUT) *
*         LDD     - NUMBER OF RESONANCE+RESONANCE COLLISIONS
*         LPP     - NUMBER OF PION+PION elastic COLIISIONS
*         lppk    - number of pion(RHO,OMEGA)+pion(RHO,OMEGA)
*                 -->K+K- collisions
*         LPN     - NUMBER OF PION+N-->KAON+X
*         lpd     - number of pion+n-->delta+pion
*         lrho    - number of pion+n-->Delta+rho
*         lomega  - number of pion+n-->Delta+omega
*         LKN     - NUMBER OF KAON RESCATTERINGS
*         LNNK    - NUMBER OF bb-->kAON PROCESS
*         LDDK    - NUMBER OF DD-->KAON PROCESS
*         LNDK    - NUMBER OF ND-->KAON PROCESS
***********************************
* TIME-INTEGRATED COLLISIONS NUMBERS OF VARIOUS PROCESSES
          RCOLL = RCOLL + FLOAT(LCOLL)/num
          RBLOC = RBLOC + FLOAT(LBLOC)/num
          RCNNE = RCNNE + FLOAT(LCNNE)/num
          RDD   = RDD   + FLOAT(LDD)/num
          RPP   = RPP   + FLOAT(LPP)/NUM
          rppk  =rppk   + float(lppk)/num
          RPN   = RPN   + FLOAT(LPN)/NUM
          rpd   =rpd    + float(lpd)/num
          RKN   = RKN   + FLOAT(LKN)/NUM
          RNNK  =RNNK   + FLOAT(LNNK)/NUM
          RDDK  =RDDK   + FLOAT(LDDK)/NUM
          RNDK  =RNDK   + FLOAT(LNDK)/NUM
          RCNND = RCNND + FLOAT(LCNND)/num
          RCNDN = RCNDN + FLOAT(LCNDN)/num
          RDIRT = RDIRT + FLOAT(LDIRT)/num
          RDECAY= RDECAY+ FLOAT(LDECAY)/num
          RRES  = RRES  + FLOAT(LRES)/num
* AVERAGE RATES OF VARIOUS COLLISIONS IN THE CURRENT TIME STEP
          ADIRT=LDIRT/DT/num
          ACOLL=(LCOLL-LBLOC)/DT/num
          ACNND=LCNND/DT/num
          ACNDN=LCNDN/DT/num
          ADECAY=LDECAY/DT/num
          ARES=LRES/DT/num
          ADOU=LDOU/DT/NUM
          ADDRHO=LDDRHO/DT/NUM
          ANNRHO=LNNRHO/DT/NUM
          ANNOM=LNNOM/DT/NUM
          ADD=LDD/DT/num
          APP=LPP/DT/num
          appk=lppk/dt/num
          APN=LPN/DT/num
          apd=lpd/dt/num
          arh=lrho/dt/num
          aom=lomega/dt/num
          AKN=LKN/DT/num
          ANNK=LNNK/DT/num
          ADDK=LDDK/DT/num
          ANDK=LNDK/DT/num
* PRINT OUT THE VARIOUS COLLISION RATES
* (1)N-N COLLISIONS 
c       WRITE(1010,9991)NT*DT,ACNND,ADOU,ADIRT,ADDRHO,ANNRHO+ANNOM
c9991       FORMAT(6(E10.3,2X))
* (2)PION-N COLLISIONS
c       WRITE(1011,'(5(E10.3,2X))')NT*DT,apd,ARH,AOM,APN
* (3)KAON PRODUCTION CHANNELS
c        WRITE(1012,9993)NT*DT,ANNK,ADDK,ANDK,APN,Appk
* (4)D(N*)+D(N*) COLLISION
c       WRITE(1013,'(4(E10.3,2X))')NT*DT,ADDK,ADD,ADD+ADDK
* (5)MESON+MESON
c       WRITE(1014,'(4(E10.3,2X))')NT*DT,APPK,APP,APP+APPK
* (6)DECAY AND RESONANCE
c       WRITE(1016,'(3(E10.3,2X))')NT*DT,ARES,ADECAY
* (7)N+D(N*)
c       WRITE(1017,'(4(E10.3,2X))')NT*DT,ACNDN,ANDK,ACNDN+ANDK
c9992    FORMAT(5(E10.3,2X))
c9993    FORMAT(6(E10.3,2X))
* PRINT OUT TIME-INTEGRATED COLLISION INFORMATION
cbz12/28/98
c        write(1018,'(5(e10.3,2x),/, 4(e10.3,2x))')
c     &           RCNNE,RCNND,RCNDN,RDIRT,rpd,
c     &           RDECAY,RRES,RDD,RPP
c        write(1018,'(6(e10.3,2x),/, 5(e10.3,2x))')
c     &           NT*DT,RCNNE,RCNND,RCNDN,RDIRT,rpd,
c     &           NT*DT,RDECAY,RRES,RDD,RPP
cbz12/18/98end
* PRINT OUT TIME-INTEGRATED KAON MULTIPLICITIES FROM DIFFERENT CHANNELS
c       WRITE(1019,'(7(E10.3,2X))')NT*DT,RNNK,RDDK,RNDK,RPN,Rppk,
c     &                           RNNK+RDDK+RNDK+RPN+Rppk
*                                                                      *
        END IF
*
*       UPDATE BARYON DENSITY
*
        CALL DENS(IPOT,MASS,NUM,OUTPAR)
*
*       UPDATE POSITIONS FOR ALL THE PARTICLES PRESENT AT THIS TIME
*
       sumene=0
        ISO=0
        DO 201 MRUN=1,NUM
        ISO=ISO+MASSR(MRUN-1)
        DO 201 I0=1,MASSR(MRUN)
        I =I0+ISO
        ETOTAL = SQRT( E(I)**2 + P(1,I)**2 + P(2,I)**2 +P(3,I)**2 )
       sumene=sumene+etotal
C for kaons, if there is a potential
C CALCULATE THE ENERGY OF THE KAON ACCORDING TO THE IMPULSE APPROXIMATION
C REFERENCE: B.A. LI AND C.M. KO, PHYS. REV. C 54 (1996) 3283. 
         if(kpoten.ne.0.and.lb(i).eq.23)then
             den=0.
              IX = NINT( R(1,I) )
              IY = NINT( R(2,I) )
              IZ = NINT( R(3,I) )
clin-4/2008:
c       IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
c     & ABS(IZ) .LT. MAXZ) den=rho(ix,iy,iz)
              IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1             .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ)
     2             den=rho(ix,iy,iz)
c         ecor=0.1973**2*0.255*kmul*4*3.14159*(1.+0.4396/0.938)
c         etotal=sqrt(etotal**2+ecor*den)
c** G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV, m^*=m
c     GeV^2 fm^3
          akg = 0.1727
c     GeV fm^3
          bkg = 0.333
         rnsg = den
         ecor = - akg*rnsg + (bkg*den)**2
         etotal = sqrt(etotal**2 + ecor)
         endif
c
         if(kpoten.ne.0.and.lb(i).eq.21)then
             den=0.
              IX = NINT( R(1,I) )
              IY = NINT( R(2,I) )
              IZ = NINT( R(3,I) )
clin-4/2008:
c       IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
c     & ABS(IZ) .LT. MAXZ) den=rho(ix,iy,iz)
              IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1             .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ)
     2             den=rho(ix,iy,iz)
c* for song potential no effect on position
c** G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV, m^*=m
c     GeV^2 fm^3
          akg = 0.1727
c     GeV fm^3
          bkg = 0.333
         rnsg = den
         ecor = - akg*rnsg + (bkg*den)**2
         etotal = sqrt(etotal**2 + ecor)
          endif
c
C UPDATE POSITIONS
          R(1,I) = R(1,I) + DT*P(1,I)/ETOTAL
          R(2,I) = R(2,I) + DT*P(2,I)/ETOTAL
          R(3,I) = R(3,I) + DT*P(3,I)/ETOTAL
c use cyclic boundary conitions
            if ( cycbox.ne.0 ) then
              if ( r(1,i).gt. cycbox/2 ) r(1,i)=r(1,i)-cycbox
              if ( r(1,i).le.-cycbox/2 ) r(1,i)=r(1,i)+cycbox
              if ( r(2,i).gt. cycbox/2 ) r(2,i)=r(2,i)-cycbox
              if ( r(2,i).le.-cycbox/2 ) r(2,i)=r(2,i)+cycbox
              if ( r(3,i).gt. cycbox/2 ) r(3,i)=r(3,i)-cycbox
              if ( r(3,i).le.-cycbox/2 ) r(3,i)=r(3,i)+cycbox
            end if
* UPDATE THE DELTA, N* AND PION COUNTERS
          LB1=LB(I)
* 1. FOR DELTA++
        IF(LB1.EQ.9)LD1=LD1+1
* 2. FOR DELTA+
        IF(LB1.EQ.8)LD2=LD2+1
* 3. FOR DELTA0
        IF(LB1.EQ.7)LD3=LD3+1
* 4. FOR DELTA-
        IF(LB1.EQ.6)LD4=LD4+1
* 5. FOR N*+(1440)
        IF(LB1.EQ.11)LN1=LN1+1
* 6. FOR N*0(1440)
        IF(LB1.EQ.10)LN2=LN2+1
* 6.1 FOR N*(1535)
       IF((LB1.EQ.13).OR.(LB1.EQ.12))LN5=LN5+1
* 6.2 FOR ETA
       IF(LB1.EQ.0)LE=LE+1
* 6.3 FOR KAONS
       IF(LB1.EQ.23)LKAON=LKAON+1
clin-11/09/00: FOR KAON*
       IF(LB1.EQ.30)LKAONS=LKAONS+1
* UPDATE PION COUNTER
* 7. FOR PION+
        IF(LB1.EQ.5)LP1=LP1+1
* 8. FOR PION0
        IF(LB1.EQ.4)LP2=LP2+1
* 9. FOR PION-
        IF(LB1.EQ.3)LP3=LP3+1
201     CONTINUE
        LP=LP1+LP2+LP3
        LD=LD1+LD2+LD3+LD4
        LN=LN1+LN2
        ALP=FLOAT(LP)/FLOAT(NUM)
        ALD=FLOAT(LD)/FLOAT(NUM)
        ALN=FLOAT(LN)/FLOAT(NUM)
       ALN5=FLOAT(LN5)/FLOAT(NUM)
        ATOTAL=ALP+ALD+ALN+0.5*ALN5
       ALE=FLOAT(LE)/FLOAT(NUM)
       ALKAON=FLOAT(LKAON)/FLOAT(NUM)
* UPDATE MOMENTUM DUE TO COULOMB INTERACTION 
        if (icou .eq. 1) then
*       with Coulomb interaction
          iso=0
          do 1026 irun = 1,num
            iso=iso+massr(irun-1)
            do 1021 il = 1,massr(irun)
               temp(1,il) = 0.
               temp(2,il) = 0.
               temp(3,il) = 0.
 1021       continue
            do 1023 il = 1, massr(irun)
              i=iso+il
              if (zet(lb(i)).ne.0) then
                do 1022 jl = 1,il-1
                  j=iso+jl
                  if (zet(lb(j)).ne.0) then
                    ddx=r(1,i)-r(1,j)
                    ddy=r(2,i)-r(2,j)
                    ddz=r(3,i)-r(3,j)
                    rdiff = sqrt(ddx**2+ddy**2+ddz**2)
                    if (rdiff .le. 1.) rdiff = 1.
                    grp=zet(lb(i))*zet(lb(j))/rdiff**3
                    ddx=ddx*grp
                    ddy=ddy*grp
                    ddz=ddz*grp
                    temp(1,il)=temp(1,il)+ddx
                    temp(2,il)=temp(2,il)+ddy
                    temp(3,il)=temp(3,il)+ddz
                    temp(1,jl)=temp(1,jl)-ddx
                    temp(2,jl)=temp(2,jl)-ddy
                    temp(3,jl)=temp(3,jl)-ddz
                  end if
 1022          continue
              end if
 1023      continue
            do 1025 il = 1,massr(irun)
              i= iso+il
              if (zet(lb(i)).ne.0) then
                do 1024 idir = 1,3
                  p(idir,i) = p(idir,i) + temp(idir,il)
     &                                    * dt * 0.00144
 1024          continue
              end if
 1025      continue
 1026   continue
        end if
*       In the following, we shall:  
*       (1) UPDATE MOMENTA DUE TO THE MEAN FIELD FOR BARYONS AND KAONS,
*       (2) calculate the thermalization, temperature in a sphere of 
*           radius 2.0 fm AROUND THE CM
*       (3) AND CALCULATE THE NUMBER OF PARTICLES IN THE HIGH DENSITY REGION 
       spt=0
       spz=0
       ncen=0
       ekin=0
          NLOST = 0
          MEAN=0
         nquark=0
         nbaryn=0
csp06/18/01
           rads = 2.
           zras = 0.1
           denst = 0.
           edenst = 0.
csp06/18/01 end
          DO 6000 IRUN = 1,NUM
          MEAN=MEAN+MASSR(IRUN-1)
          DO 5800 J = 1,MASSR(irun)
          I=J+MEAN
c
csp06/18/01
           radut = sqrt(r(1,i)**2+r(2,i)**2)
       if( radut .le. rads )then
        if( abs(r(3,i)) .le. zras*nt*dt )then
c         vols = 3.14159*radut**2*abs(r(3,i))      ! cylinder pi*r^2*l
c     cylinder pi*r^2*l
         vols = 3.14159*rads**2*zras
         engs=sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)
         gammas=1.
         if(e(i).ne.0.)gammas=engs/e(i)
c     rho
         denst = denst + 1./gammas/vols
c     energy density
         edenst = edenst + engs/gammas/gammas/vols
        endif
       endif
csp06/18/01 end
c
         drr=sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
         if(drr.le.2.0)then
         spt=spt+p(1,i)**2+p(2,i)**2
         spz=spz+p(3,i)**2
         ncen=ncen+1
         ekin=ekin+sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)-e(i)
         endif
              IX = NINT( R(1,I) )
              IY = NINT( R(2,I) )
              IZ = NINT( R(3,I) )
C calculate the No. of particles in the high density region
clin-4/2008:
c              IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
c     & ABS(IZ) .LT. MAXZ) THEN
              IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1          .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ) THEN
       if(rho(ix,iy,iz)/0.168.gt.dencut)go to 5800
       if((rho(ix,iy,iz)/0.168.gt.5.).and.(e(i).gt.0.9))
     &  nbaryn=nbaryn+1
       if(pel(ix,iy,iz).gt.2.0)nquark=nquark+1
       endif
c*
c If there is a kaon potential, propogating kaons 
        if(kpoten.ne.0.and.lb(i).eq.23)then
        den=0.
clin-4/2008:
c       IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
c     & ABS(IZ) .LT. MAXZ)then
        IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1       .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ) THEN
           den=rho(ix,iy,iz)
c        ecor=0.1973**2*0.255*kmul*4*3.14159*(1.+0.4396/0.938)
c       etotal=sqrt(P(1,i)**2+p(2,I)**2+p(3,i)**2+e(i)**2+ecor*den)
c** for G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV
c     !! GeV^2 fm^3
            akg = 0.1727
c     !! GeV fm^3
            bkg = 0.333
          rnsg = den
          ecor = - akg*rnsg + (bkg*den)**2
          etotal = sqrt(P(1,i)**2+p(2,I)**2+p(3,i)**2+e(i)**2 + ecor)
          ecor = - akg + 2.*bkg**2*den + 2.*bkg*etotal
c** G.Q. Li potential (END)           
        CALL GRADUK(IX,IY,IZ,GRADXk,GRADYk,GRADZk)
        P(1,I) = P(1,I) - DT * GRADXk*ecor/(2.*etotal)
        P(2,I) = P(2,I) - DT * GRADYk*ecor/(2.*etotal)
        P(3,I) = P(3,I) - DT * GRADZk*ecor/(2.*etotal)
        endif
         endif
c
        if(kpoten.ne.0.and.lb(i).eq.21)then
         den=0.
clin-4/2008:
c           IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
c     &        ABS(IZ) .LT. MAXZ)then
         IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1        .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ) THEN
               den=rho(ix,iy,iz)
        CALL GRADUK(IX,IY,IZ,GRADXk,GRADYk,GRADZk)
c        P(1,I) = P(1,I) - DT * GRADXk*(-0.12/0.168)    !! song potential
c        P(2,I) = P(2,I) - DT * GRADYk*(-0.12/0.168)
c        P(3,I) = P(3,I) - DT * GRADZk*(-0.12/0.168)
c** for G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV
c    !! GeV^2 fm^3
            akg = 0.1727
c     !! GeV fm^3
            bkg = 0.333
          rnsg = den
          ecor = - akg*rnsg + (bkg*den)**2
          etotal = sqrt(P(1,i)**2+p(2,I)**2+p(3,i)**2+e(i)**2 + ecor)
          ecor = - akg + 2.*bkg**2*den - 2.*bkg*etotal
        P(1,I) = P(1,I) - DT * GRADXk*ecor/(2.*etotal)
        P(2,I) = P(2,I) - DT * GRADYk*ecor/(2.*etotal)
        P(3,I) = P(3,I) - DT * GRADZk*ecor/(2.*etotal)
c** G.Q. Li potential (END)           
        endif
         endif
c
c for other mesons, there is no potential
       if(j.gt.mass)go to 5800         
c  with mean field interaction for baryons   (open endif below) !!sp05
**      if( (iabs(lb(i)).eq.1.or.iabs(lb(i)).eq.2) .or.
**    &     (iabs(lb(i)).ge.6.and.iabs(lb(i)).le.17) .or.
**    &      iabs(lb(i)).eq.40.or.iabs(lb(i)).eq.41 )then  
        IF (ICOLL .NE. -1) THEN
* check if the baryon has run off the lattice
*             IX0=NINT(R(1,I)/DX)
*             IY0=NINT(R(2,I)/DY)
*             IZ0=NINT(R(3,I)/DZ)
*             IPX0=NINT(P(1,I)/DPX)
*             IPY0=NINT(P(2,I)/DPY)
*             IPZ0=NINT(P(3,I)/DPZ)
*      if ( (abs(ix0).gt.mx) .or. (abs(iy0).gt.my) .or. (abs(iz0).gt.mz)
*     &  .or. (abs(ipx0).gt.mpx) .or. (abs(ipy0) 
*     &  .or. (ipz0.lt.-mpz) .or. (ipz0.gt.mpzp)) NLOST=NLOST+1
clin-4/2008:
c              IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
c     &                                    ABS(IZ) .LT. MAXZ     ) THEN
           IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1          .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ) THEN
                CALL GRADU(IPOT,IX,IY,IZ,GRADX,GRADY,GRADZ)
              TZ=0.
              GRADXN=0
              GRADYN=0
              GRADZN=0
              GRADXP=0
              GRADYP=0
              GRADZP=0
             IF(ICOU.EQ.1)THEN
                CALL GRADUP(IX,IY,IZ,GRADXP,GRADYP,GRADZP)
                CALL GRADUN(IX,IY,IZ,GRADXN,GRADYN,GRADZN)
               IF(ZET(LB(I)).NE.0)TZ=-1
               IF(ZET(LB(I)).EQ.0)TZ= 1
             END IF
           if(iabs(lb(i)).ge.14.and.iabs(lb(i)).le.17)then
              facl = 2./3.
            elseif(iabs(lb(i)).eq.40.or.iabs(lb(i)).eq.41)then
              facl = 1./3.
            else
              facl = 1.
            endif
        P(1,I) = P(1,I) - facl*DT * (GRADX+asy*(GRADXN-GRADXP)*TZ)
        P(2,I) = P(2,I) - facl*DT * (GRADY+asy*(GRADYN-GRADYP)*TZ)
        P(3,I) = P(3,I) - facl*DT * (GRADZ+asy*(GRADZN-GRADZP)*TZ)
                end if                                                       
              ENDIF
**          endif          !!sp05     
 5800       CONTINUE
 6000       CONTINUE
c print out the average no. of particles in regions where the local 
c baryon density is higher than 5*rho0 
c       write(1072,'(e10.3,2x,e10.3)')nt*dt,float(nbaryn)/float(num)
C print out the average no. of particles in regions where the local 
c energy density is higher than 2 GeV/fm^3. 
c       write(1073,'(e10.3,2x,e10.3)')nt*dt,float(nquark)/float(num)
c print out the no. of particles that have run off the lattice
*          IF (NLOST .NE. 0 .AND. (NT/NFREQ)*NFREQ .EQ. NT) THEN
*            WRITE(12,'(5X,''***'',I7,'' TESTPARTICLES LOST AFTER '',
*     &                   ''TIME STEP NUMBER'',I4)') NLOST, NT
*         END IF
*
*       update phase space density
*        call platin(mode,mass,num,dx,dy,dz,dpx,dpy,dpz,fnorm)
*
*       CONTROL-PRINTOUT OF CONFIGURATION (IF REQUIRED)
*
*        if (inout(5) .eq. 2) CALL ENERGY(NT,IPOT,NUM,MASS,EMIN,EMAX)
*
* 
* print out central baryon density as a function of time
       CDEN=RHO(0,0,0)/0.168
cc        WRITE(1002,990)FLOAT(NT)*DT,CDEN
c        WRITE(1002,1990)FLOAT(NT)*DT,CDEN,denst/real(num)
* print out the central energy density as a function of time
cc        WRITE(1003,990)FLOAT(NT)*DT,PEL(0,0,0)
c        WRITE(1003,1990)FLOAT(NT)*DT,PEL(0,0,0),edenst/real(num)
* print out the no. of pion-like particles as a function of time 
c        WRITE(1004,9999)FLOAT(NT)*DT,ALD,ALN,ALP,ALN5,
c     &               ALD+ALN+ALP+0.5*ALN5
* print out the no. of eta-like particles as a function of time
c        WRITE(1005,991)FLOAT(NT)*DT,ALN5,ALE,ALE+0.5*ALN5
c990       FORMAT(E10.3,2X,E10.3)
c1990       FORMAT(E10.3,2X,E10.3,2X,E10.3)
c991       FORMAT(E10.3,2X,E10.3,2X,E10.3,2X,E10.3)
c9999    FORMAT(e10.3,2X,e10.3,2X,E10.3,2X,E10.3,2X,
c     1  E10.3,2X,E10.3)
C THE FOLLOWING OUTPUTS CAN BE TURNED ON/OFF by setting icflow and icrho=0  
c print out the baryon and meson density matrix in the reaction plane
        IF ((NT/NFREQ)*NFREQ .EQ. NT ) THEN
       if(icflow.eq.1)call flow(nt)
cbz11/18/98
c       if(icrho.ne.1)go to 10000 
c       if (icrho .eq. 1) then 
cbz11/18/98end
c       do ix=-10,10
c       do iz=-10,10
c       write(1053,992)ix,iz,rho(ix,0,iz)/0.168
c       write(1054,992)ix,iz,pirho(ix,0,iz)/0.168
c       write(1055,992)ix,iz,pel(ix,0,iz)
c       end do
c       end do
cbz11/18/98
c        end if
cbz11/18/98end
c992       format(i3,i3,e11.4)
       endif
c print out the ENERGY density matrix in the reaction plane
C CHECK LOCAL MOMENTUM EQUILIBRIUM IN EACH CELL, 
C AND PERFORM ON-LINE FLOW ANALYSIS AT A FREQUENCY OF NFREQ
c        IF ((NT/NFREQ)*NFREQ .EQ. NT ) THEN
c       call flow(nt)
c       call equ(ipot,mass,num,outpar)
c       do ix=-10,10
c       do iz=-10,10
c       write(1055,992)ix,iz,pel(ix,0,iz)
c       write(1056,992)ix,iz,rxy(ix,0,iz)
c       end do
c       end do
c       endif
C calculate the volume of high BARYON AND ENERGY density 
C matter as a function of time
c       vbrho=0.
c       verho=0.
c       do ix=-20,20
c       do iy=-20,20
c       do iz=-20,20
c       if(rho(ix,iy,iz)/0.168.gt.5.)vbrho=vbrho+1.
c       if(pel(ix,iy,iz).gt.2.)verho=verho+1.
c       end do
c       end do
c       end do
c       write(1081,993)dt*nt,vbrho
c       write(1082,993)dt*nt,verho
c993       format(e11.4,2x,e11.4)
*-----------------------------------------------------------------------
cbz11/16/98
c.....for read-in initial conditions produce particles from read-in 
c.....common block.
c.....note that this part is only for cascade with number of test particles
c.....NUM = 1.
      IF (IAPAR2(1) .NE. 1) THEN
         CT = NT * DT
cbz12/22/98
c         NP = MASSR(1)
c         DO WHILE (FTAR(NPI) .GT. CT - DT .AND. FTAR(NPI) .LE. CT)
c            NP = NP + 1
c            R(1, NP) = GXAR(NPI) + PXAR(NPI) / PEAR(NPI) * (CT - FTAR(NPI))
c            R(2, NP) = GYAR(NPI) + PYAR(NPI) / PEAR(NPI) * (CT - FTAR(NPI))
c            R(3, NP) = GZAR(NPI) + PZAR(NPI) / PEAR(NPI) * (CT - FTAR(NPI))
c            P(1, NP) = PXAR(NPI)
c            P(2, NP) = PYAR(NPI)
c            P(3, NP) = PZAR(NPI)
c            E(NP) = XMAR(NPI)
c            LB(NP) = IARFLV(ITYPAR(NPI))
c            NPI = NPI + 1
c         END DO
c         MASSR(1) = NP
         IA = 0
         DO 1028 IRUN = 1, NUM
            DO 1027 IC = 1, MASSR(IRUN)
               IE = IA + IC
               RT(1, IC, IRUN) = R(1, IE)
               RT(2, IC, IRUN) = R(2, IE)
               RT(3, IC, IRUN) = R(3, IE)
               PT(1, IC, IRUN) = P(1, IE)
               PT(2, IC, IRUN) = P(2, IE)
               PT(3, IC, IRUN) = P(3, IE)
               ET(IC, IRUN) = E(IE)
               LT(IC, IRUN) = LB(IE)
c         !! sp 12/19/00
               PROT(IC, IRUN) = PROPER(IE)
clin-5/2008:
               dpertt(IC, IRUN)=dpertp(IE)
 1027       CONTINUE
            NP = MASSR(IRUN)
            NP1 = NPI(IRUN)
cbz10/05/99
c            DO WHILE (FT1(NP1, IRUN) .GT. CT - DT .AND. 
c     &           FT1(NP1, IRUN) .LE. CT)
cbz10/06/99
c            DO WHILE (NPI(IRUN).LE.MULTI1(IRUN).AND.
cbz10/06/99 end
clin-11/13/00 finally read in all unformed particles and do the decays in ART:
c           DO WHILE (NP1.LE.MULTI1(IRUN).AND.
c    &           FT1(NP1, IRUN) .GT. CT - DT .AND. 
c    &           FT1(NP1, IRUN) .LE. CT)
c
               ctlong = ct
             if(nt .eq. (ntmax-1))then
               ctlong = 1.E30
             elseif(nt .eq. ntmax)then
               go to 1111
             endif
c
            DO WHILE (NP1.LE.MULTI1(IRUN).AND.
     &           FT1(NP1, IRUN) .GT. ((NT-1) * DT) .AND. 
     &           FT1(NP1, IRUN) .LE. ctlong)
clin-ma-5/2016 changed the following to 2nd line above to avoid bug 
c     that leads to loss of hadrons inside ART due to finite accuracy 
c     [which results in (ct-dt) + dt != ct exactly]:
c     &           FT1(NP1, IRUN) .GT. (CT - DT) .AND. 
               NP = NP + 1
               UDT = (CT - FT1(NP1, IRUN)) / EE1(NP1, IRUN)
clin-10/28/03 since all unformed hadrons at time ct are read in at nt=ntmax-1, 
c     their positions should not be propagated to time ct:
               if(nt.eq.(ntmax-1)) then
                  ftsvt(NP,IRUN)=FT1(NP1, IRUN)
                  if(FT1(NP1, IRUN).gt.ct) UDT=0.
               endif
               RT(1, NP, IRUN) = GX1(NP1, IRUN) + 
     &              PX1(NP1, IRUN) * UDT
               RT(2, NP, IRUN) = GY1(NP1, IRUN) + 
     &              PY1(NP1, IRUN) * UDT
               RT(3, NP, IRUN) = GZ1(NP1, IRUN) + 
     &              PZ1(NP1, IRUN) * UDT
               PT(1, NP, IRUN) = PX1(NP1, IRUN)
               PT(2, NP, IRUN) = PY1(NP1, IRUN)
               PT(3, NP, IRUN) = PZ1(NP1, IRUN)
               ET(NP, IRUN) = XM1(NP1, IRUN)
               LT(NP, IRUN) = IARFLV(ITYP1(NP1, IRUN))
clin-5/2008:
               dpertt(NP,IRUN)=dpp1(NP1,IRUN)
clin-4/30/03 ctest off 
c     record initial phi,K*,Lambda(1520) resonances formed during the timestep:
c               if(LT(NP, IRUN).eq.29.or.iabs(LT(NP, IRUN)).eq.30)
c     1              write(17,112) 'formed',LT(NP, IRUN),PX1(NP1, IRUN),
c     2 PY1(NP1, IRUN),PZ1(NP1, IRUN),XM1(NP1, IRUN),nt
c 112           format(a10,1x,I4,4(1x,f9.3),1x,I4)
c
               NP1 = NP1 + 1
c     !! sp 12/19/00
               PROT(NP, IRUN) = 1.
            END DO
*
 1111      continue
            NPI(IRUN) = NP1
            IA = IA + MASSR(IRUN)
            MASSR(IRUN) = NP
 1028    CONTINUE
         IA = 0
         DO 1030 IRUN = 1, NUM
            IA = IA + MASSR(IRUN - 1)
            DO 1029 IC = 1, MASSR(IRUN)
               IE = IA + IC
               R(1, IE) = RT(1, IC, IRUN)
               R(2, IE) = RT(2, IC, IRUN)
               R(3, IE) = RT(3, IC, IRUN)
               P(1, IE) = PT(1, IC, IRUN)
               P(2, IE) = PT(2, IC, IRUN)
               P(3, IE) = PT(3, IC, IRUN)
               E(IE) = ET(IC, IRUN)
               LB(IE) = LT(IC, IRUN)
c     !! sp 12/19/00
               PROPER(IE) = PROT(IC, IRUN)
               if(nt.eq.(ntmax-1)) ftsv(IE)=ftsvt(IC,IRUN)
clin-5/2008:
               dpertp(IE)=dpertt(IC, IRUN)
 1029       CONTINUE
clin-3/2009 Moved here to better take care of freezeout spacetime:
            call hbtout(MASSR(IRUN),nt,ntmax)
 1030    CONTINUE
cbz12/22/98end
      END IF
cbz11/16/98end
clin-5/2009 ctest off:
c      call flowh(ct) 
10000       continue
*                                                                      *
*       ==============  END OF TIME STEP LOOP   ================       *
************************************
*     WRITE OUT particle's MOMENTA ,and/OR COORDINATES ,
*     label and/or their local baryon density in the final state
        iss=0
        do 1032 lrun=1,num
           iss=iss+massr(lrun-1)
           do 1031 l0=1,massr(lrun)
              ipart=iss+l0
 1031      continue
 1032   continue
cbz11/16/98
      IF (IAPAR2(1) .NE. 1) THEN
cbz12/22/98
c        NSH = MASSR(1) - NPI + 1
c        IAINT2(1) = IAINT2(1) + NSH
c.....to shift the unformed particles to the end of the common block
c        IF (NSH .GT. 0) THEN
c           IB = IAINT2(1)
c           IE = MASSR(1) + 1
c           II = -1
c        ELSE IF (NSH .LT. 0) THEN
c           IB = MASSR(1) + 1
c           IE = IAINT2(1)
c           II = 1
c        END IF
c        IF (NSH .NE. 0) THEN
c           DO I = IB, IE, II
c              J = I - NSH
c              ITYPAR(I) = ITYPAR(J)
c              GXAR(I) = GXAR(J)
c              GYAR(I) = GYAR(J)
c              GZAR(I) = GZAR(J)
c              FTAR(I) = FTAR(J)
c              PXAR(I) = PXAR(J)
c              PYAR(I) = PYAR(J)
c              PZAR(I) = PZAR(J)
c              PEAR(I) = PEAR(J)
c              XMAR(I) = XMAR(J)
c           END DO
c        END IF
c.....to copy ART particle info to COMMON /ARPRC/
c        DO I = 1, MASSR(1)
c           ITYPAR(I) = INVFLV(LB(I))
c           GXAR(I) = R(1, I)
c           GYAR(I) = R(2, I)
c           GZAR(I) = R(3, I)
c           FTAR(I) = CT
c           PXAR(I) = P(1, I)
c           PYAR(I) = P(2, I)
c           PZAR(I) = P(3, I)
c           XMAR(I) = E(I)
c           PEAR(I) = SQRT(PXAR(I) ** 2 + PYAR(I) ** 2 + PZAR(I) ** 2
c     &        + XMAR(I) ** 2)
c        END DO
        IA = 0
        DO 1035 IRUN = 1, NUM
           IA = IA + MASSR(IRUN - 1)
           NP1 = NPI(IRUN)
           NSH = MASSR(IRUN) - NP1 + 1
           MULTI1(IRUN) = MULTI1(IRUN) + NSH
c.....to shift the unformed particles to the end of the common block
           IF (NSH .GT. 0) THEN
              IB = MULTI1(IRUN)
              IE = MASSR(IRUN) + 1
              II = -1
           ELSE IF (NSH .LT. 0) THEN
              IB = MASSR(IRUN) + 1
              IE = MULTI1(IRUN)
              II = 1
           END IF
           IF (NSH .NE. 0) THEN
              DO 1033 I = IB, IE, II
                 J = I - NSH
                 ITYP1(I, IRUN) = ITYP1(J, IRUN)
                 GX1(I, IRUN) = GX1(J, IRUN)
                 GY1(I, IRUN) = GY1(J, IRUN)
                 GZ1(I, IRUN) = GZ1(J, IRUN)
                 FT1(I, IRUN) = FT1(J, IRUN)
                 PX1(I, IRUN) = PX1(J, IRUN)
                 PY1(I, IRUN) = PY1(J, IRUN)
                 PZ1(I, IRUN) = PZ1(J, IRUN)
                 EE1(I, IRUN) = EE1(J, IRUN)
                 XM1(I, IRUN) = XM1(J, IRUN)
c     !! sp 12/19/00
                 PRO1(I, IRUN) = PRO1(J, IRUN)
clin-5/2008:
                 dpp1(I,IRUN)=dpp1(J,IRUN)
 1033         CONTINUE
           END IF
c.....to copy ART particle info to COMMON /ARPRC1/
           DO 1034 I = 1, MASSR(IRUN)
              IB = IA + I
              ITYP1(I, IRUN) = INVFLV(LB(IB))
              GX1(I, IRUN) = R(1, IB)
              GY1(I, IRUN) = R(2, IB)
              GZ1(I, IRUN) = R(3, IB)
clin-10/28/03:
c since all unformed hadrons at time ct are read in at nt=ntmax-1, 
c their formation time ft1 should be kept to determine their freezeout(x,t):
c              FT1(I, IRUN) = CT
              if(FT1(I, IRUN).lt.CT) FT1(I, IRUN) = CT
              PX1(I, IRUN) = P(1, IB)
              PY1(I, IRUN) = P(2, IB)
              PZ1(I, IRUN) = P(3, IB)
              XM1(I, IRUN) = E(IB)
              EE1(I, IRUN) = SQRT(PX1(I, IRUN) ** 2 + 
     &             PY1(I, IRUN) ** 2 +
     &             PZ1(I, IRUN) ** 2 + 
     &             XM1(I, IRUN) ** 2)
c     !! sp 12/19/00
              PRO1(I, IRUN) = PROPER(IB)
 1034      CONTINUE
 1035   CONTINUE
cbz12/22/98end
      END IF
cbz11/16/98end
c
**********************************
*                                                                      *
*       ======= END OF MANY LOOPS OVER IMPACT PARAMETERS ==========    *
*                                                               *
**********************************
50000   CONTINUE
*
*-----------------------------------------------------------------------
*                       ==== ART COMPLETED ====
*-----------------------------------------------------------------------
cbz11/16/98
c      STOP
      RETURN
cbz11/16/98end
      END
**********************************
