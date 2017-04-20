        SUBROUTINE JETINI(JP,JT,itrig)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  itrig = 0:对于正常的过程。
c$$$      itrig = 1:对于triggered过程
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
C*******Initialize PYTHIA for jet production**********************
C        itrig=0: for normal processes
C        itrig=1: for triggered processes
C       JP: sequence number of the projectile
C       JT: sequence number of the target
C     For A+A collisions, one has to initilize pythia
C     separately for each type of collisions, pp, pn,np and nn,
C     or hp and hn for hA collisions. In this subroutine we use the following
C     catalogue for different type of collisions:
C     h+h: h+h (itype=1)
C     h+A: h+p (itype=1), h+n (itype=2)
C     A+h: p+h (itype=1), n+h (itype=2)
C     A+A: p+p (itype=1), p+n (itype=2), n+p (itype=3), n+n (itype=4)
C*****************************************************************
        CHARACTER BEAM*16,TARG*16
        DIMENSION XSEC0(8,0:200),COEF0(8,200,20),INI(8),
     &                MINT44(8),MINT45(8)
        COMMON/hjcrdn/YP(3,300),YT(3,300)
cc      SAVE /hjcrdn/
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
        COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
cc      SAVE /HSTRNG/
        COMMON/HPINT/MINT4,MINT5,ATCO(200,20),ATXS(0:200)
cc      SAVE /HPINT/
C
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
cc      SAVE /LUDAT1/
        COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
cc      SAVE /LUDAT3/
        COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
cc      SAVE /PYSUBS/
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
cc      SAVE /PYPARS/
        COMMON/PYINT1/MINT(400),VINT(400)
cc      SAVE /PYINT1/
        COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
cc      SAVE /PYINT2/
        COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3)
        common /xiaohai/xiaohaiflag
cc      SAVE /PYINT5/
        SAVE
clin        DATA INI/8*0/ilast/-1/
        DATA INI/8*0/,ilast/-1/
C
        IHNT2(11)=JP
        IHNT2(12)=JT
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9991, *)"IHNT2(5) = ", IHNT2(5),
     &      "  IHNT2(6) = ", IHNT2(6)
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    IHNT2(5)和IHNT2(6):弹核和靶核的flavor 代码。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IF(IHNT2(5).NE.0 .AND. IHNT2(6).NE.0) THEN
           itype=1
        ELSE IF(IHNT2(5).NE.0 .AND. IHNT2(6).EQ.0) THEN
           itype=1
           IF(NFT(JT,4).EQ.2112) itype=2
        ELSE IF(IHNT2(5).EQ.0 .AND. IHNT2(6).NE.0) THEN
           itype=1
           IF(NFP(JP,4).EQ.2112) itype=2
        ELSE
           IF(NFP(JP,4).EQ.2212 .AND. NFT(JT,4).EQ.2212) THEN
              itype=1
           ELSE IF(NFP(JP,4).EQ.2212 .AND. NFT(JT,4).EQ.2112) THEN
              itype=2
           ELSE IF(NFP(JP,4).EQ.2112 .AND. NFT(JT,4).EQ.2212) THEN
              itype=3
           ELSE
              itype=4
           ENDIF
        ENDIF
clin-12/2012 correct NN differential cross section in HIJING:
c        write(94,*) 'In JETINI: ',jp,jt,NFP(JP,4),NFT(JT,4),itype
c
        IF(itrig.NE.0) GO TO 160
        IF(itrig.EQ.ilast) GO TO 150
        MSTP(2)=2
c                        ********second order running alpha_strong
        MSTP(33)=1
        PARP(31)=HIPR1(17)
C                        ********inclusion of K factor
        MSTP(51)=3
C                        ********Duke-Owens set 1 structure functions
        MSTP(61)=1
C                        ********INITIAL STATE RADIATION
        MSTP(71)=1
C                        ********FINAL STATE RADIATION
        IF(IHPR2(2).EQ.0.OR.IHPR2(2).EQ.2) MSTP(61)=0
        IF(IHPR2(2).EQ.0.OR.IHPR2(2).EQ.1) MSTP(71)=0
c
        MSTP(81)=0
C                        ******** NO MULTIPLE INTERACTION
        MSTP(82)=1
C                        *******STRUCTURE OF MUTLIPLE INTERACTION
        MSTP(111)=0
C                ********frag off(have to be done by local call)
        IF(IHPR2(10).EQ.0) MSTP(122)=0
C                ********No printout of initialization information
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(8): (D=2.0 GeV/c) minimum P T transfer in hard or
c$$$        semihard scatterings.硬和半硬的散射。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PARP(81)=HIPR1(8)
        CKIN(5)=HIPR1(8)
        CKIN(3)=HIPR1(8)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(9): (D=−1.0 GeV/c) maximum P T transfer in hard
c$$$        or semihard scatterings. If
c$$$        negative, the limit is set by the colliding energy.最大的横
c$$$      动量。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        CKIN(4)=HIPR1(9)
        IF(HIPR1(9).LE.HIPR1(8)) CKIN(4)=-1.0
        CKIN(9)=-10.0
        CKIN(10)=10.0
        MSEL=0
        DO 100 ISUB=1,200
           MSUB(ISUB)=0
 100    CONTINUE
        MSUB(11)=1
        MSUB(12)=1
        MSUB(13)=1
        MSUB(28)=1
        MSUB(53)=1
        MSUB(68)=1
        MSUB(81)=1
        MSUB(82)=1
        DO 110 J=1,MIN(8,MDCY(21,3))
 110    MDME(MDCY(21,2)+J-1,1)=0
        ISEL=4
        IF(HINT1(1).GE.20.0 .and. IHPR2(18).EQ.1) ISEL=5
        MDME(MDCY(21,2)+ISEL-1,1)=1
C                        ********QCD subprocesses
        MSUB(14)=1
        MSUB(18)=1
        MSUB(29)=1
C                       ******* direct photon production
 150    IF(INI(itype).NE.0) GO TO 800
        GO TO 400
C
C        *****triggered subprocesses, jet, photon, heavy quark and DY
C
 160    itype=4+itype
        IF(itrig.EQ.ilast) GO TO 260
        PARP(81)=ABS(HIPR1(10))-0.25
        CKIN(5)=ABS(HIPR1(10))-0.25
        CKIN(3)=ABS(HIPR1(10))-0.25
        CKIN(4)=ABS(HIPR1(10))+0.25
        IF(HIPR1(10).LT.HIPR1(8)) CKIN(4)=-1.0
c
        MSEL=0
        DO 101 ISUB=1,200
           MSUB(ISUB)=0
 101    CONTINUE
        IF(IHPR2(3).EQ.1) THEN
           MSUB(11)=1
           MSUB(12)=1
           MSUB(13)=1
           MSUB(28)=1
           MSUB(53)=1
           MSUB(68)=1
           MSUB(81)=1
           MSUB(82)=1
           MSUB(14)=1
           MSUB(18)=1
           MSUB(29)=1
           DO 102 J=1,MIN(8,MDCY(21,3))
 102           MDME(MDCY(21,2)+J-1,1)=0
           ISEL=4
           IF(HINT1(1).GE.20.0 .and. IHPR2(18).EQ.1) ISEL=5
           MDME(MDCY(21,2)+ISEL-1,1)=1
C                        ********QCD subprocesses
        ELSE IF(IHPR2(3).EQ.2) THEN
           MSUB(14)=1
           MSUB(18)=1
           MSUB(29)=1
C                ********Direct photon production
c                q+qbar->g+gamma,q+qbar->gamma+gamma, q+g->q+gamma
        ELSE IF(IHPR2(3).EQ.3) THEN
           CKIN(3)=MAX(0.0,HIPR1(10))
           CKIN(5)=HIPR1(8)
           PARP(81)=HIPR1(8)
           MSUB(81)=1
           MSUB(82)=1
           DO 105 J=1,MIN(8,MDCY(21,3))
 105           MDME(MDCY(21,2)+J-1,1)=0
           ISEL=4
           IF(HINT1(1).GE.20.0 .and. IHPR2(18).EQ.1) ISEL=5
           MDME(MDCY(21,2)+ISEL-1,1)=1
C             **********Heavy quark production
        ENDIF
260        IF(INI(itype).NE.0) GO TO 800
C
C
400        INI(itype)=1
        IF(IHPR2(10).EQ.0) MSTP(122)=0
        IF(NFP(JP,4).EQ.2212) THEN
                BEAM='P'
        ELSE IF(NFP(JP,4).EQ.-2212) THEN
                BEAM='P~'
        ELSE IF(NFP(JP,4).EQ.2112) THEN
                BEAM='N'
        ELSE IF(NFP(JP,4).EQ.-2112) THEN
                BEAM='N~'
        ELSE IF(NFP(JP,4).EQ.211) THEN
                BEAM='PI+'
        ELSE IF(NFP(JP,4).EQ.-211) THEN
                BEAM='PI-'
        ELSE IF(NFP(JP,4).EQ.321) THEN
                BEAM='PI+'
        ELSE IF(NFP(JP,4).EQ.-321) THEN
                BEAM='PI-'
        ELSE
                WRITE(6,*) 'unavailable beam type', NFP(JP,4)
        ENDIF
        IF(NFT(JT,4).EQ.2212) THEN
                TARG='P'
        ELSE IF(NFT(JT,4).EQ.-2212) THEN
                TARG='P~'
        ELSE IF(NFT(JT,4).EQ.2112) THEN
                TARG='N'
        ELSE IF(NFT(JT,4).EQ.-2112) THEN
                TARG='N~'
        ELSE IF(NFT(JT,4).EQ.211) THEN
                TARG='PI+'
        ELSE IF(NFT(JT,4).EQ.-211) THEN
                TARG='PI-'
        ELSE IF(NFT(JT,4).EQ.321) THEN
                TARG='PI+'
        ELSE IF(NFT(JT,4).EQ.-321) THEN
                TARG='PI-'
        ELSE
                WRITE(6,*) 'unavailable target type', NFT(JT,4)
        ENDIF
C
        IHNT2(16)=1
C       ******************indicate for initialization use when
C                         structure functions are called in PYTHIA
C
        CALL PYINIT('CMS',BEAM,TARG,HINT1(1))
        MINT4=MINT(44)
        MINT5=MINT(45)
        MINT44(itype)=MINT(44)
        MINT45(itype)=MINT(45)
        ATXS(0)=XSEC(0,1)
        XSEC0(itype,0)=XSEC(0,1)
        DO 500 I=1,200
                ATXS(I)=XSEC(I,1)
                XSEC0(itype,I)=XSEC(I,1)
                DO 500 J=1,20
                        ATCO(I,J)=COEF(I,J)
                        COEF0(itype,I,J)=COEF(I,J)
500        CONTINUE
C
        IHNT2(16)=0
C
        RETURN
C                ********Store the initialization information for
C                                late use
C
C
800        MINT(44)=MINT44(itype)
        MINT(45)=MINT45(itype)
        MINT4=MINT(44)
        MINT5=MINT(45)
        XSEC(0,1)=XSEC0(itype,0)
        ATXS(0)=XSEC(0,1)
        DO 900 I=1,200
                XSEC(I,1)=XSEC0(itype,I)
                ATXS(I)=XSEC(I,1)
        DO 900 J=1,20
                COEF(I,J)=COEF0(itype,I,J)
                ATCO(I,J)=COEF(I,J)
900        CONTINUE
        ilast=itrig
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  NFP(I,4):弹核I的原始的original flavor代码。
c$$$        NFT(I,4):靶核I的原始的orininal flavor代码。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        MINT(11)=NFP(JP,4)
        MINT(12)=NFT(JT,4)
        RETURN
        END
C            
C
C
