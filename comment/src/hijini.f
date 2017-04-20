c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        此程序的这些数据参看HIJING的开头那些数据.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        SUBROUTINE HIJINI
        PARAMETER (MAXSTR=150001)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
        COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
cc      SAVE /HSTRNG/
        COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &                PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &                PJPM(300,500),NTJ(300),KFTJ(300,500),
     &                PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &                PJTE(300,500),PJTM(300,500)
cc      SAVE /HJJET1/
        COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &       K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &       PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
cc      SAVE /HJJET2/
c        COMMON/HJJET4/NDR,IADR(900,2),KFDR(900),PDR(900,5)
        COMMON/HJJET4/NDR,IADR(MAXSTR,2),KFDR(MAXSTR),PDR(MAXSTR,5)
cc      SAVE /HJJET4/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      common /xiaohai/xiaohaiflag
      SAVE   
C****************Reset the momentum of initial particles************
C             and assign flavors to the proj and targ string       *
C*******************************************************************
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NSG: the total number of such string systems.(string的个数.)
c$$$        NDR:total number of directly produced particles.(直接产生的粒子
c$$$        总数)
c$$$        将IPP和IPT都赋值为质子.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        NSG=0
        NDR=0
        IPP=2212
        IPT=2212
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IPP and IPT应该靶核和弹核核子的flavor code.
c$$$        程序默认先从质子进行的处理,即IPP和IPT都为2212.
c$$$        IHNT2(5):the flavor code of the projectile hadron.
c$$$        IHNT2(6):the flavor code of the target hadron.
c$$$        理论上都应该为0,即下面的条件语句不执行.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9993, *) "IHNT2(5) = ",IHNT2(5),"   IHNT2(6)",IHNT2(6)
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        IF(IHNT2(5).NE.0) IPP=IHNT2(5)
        IF(IHNT2(6).NE.0) IPT=IHNT2(6)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9993, *) "IPP = ", IPP,"   ITT = ",ITT
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
C                ********in case the proj or targ is a hadron.
C
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        对于弹核的每个核子。
c$$$        PP(I,1),PP(I,2),PP(I,3),PP(I,4),PP(I,5):弹核的第I个核子的四动量
c$$$                        和不变质量。
c$$$        PP(I,6),PP(I,7):弹核的第I个核子的价夸克的横向动量(Px, Py).
c$$$        PP(I,8),PP(I,9):弹核的第I个核子的反价夸克(anti-quark)的
c$$$                       横向动量(Px, Py).
c$$$       PP(I,10),PP(I,11),PP(I,12):弹核的第I个核子在最后一次硬散射
c$$$                       夸克或者反夸克所转移的三动量(Px,Py,Pz).
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        遍历所有的核子IHNT2(1)为质量数。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<

        DO 100 I=1,IHNT2(1)
        PP(I,1)=0.0
        PP(I,2)=0.0
        PP(I,3)=SQRT(HINT1(1)**2/4.0-HINT1(8)**2)
        PP(I,4)=HINT1(1)/2
        PP(I,5)=HINT1(8)
        PP(I,6)=0.0
        PP(I,7)=0.0
        PP(I,8)=0.0
        PP(I,9)=0.0
        PP(I,10)=0.0
cbzdbg2/22/99
ctest OFF
        PP(I, 11) = 0.0
        PP(I, 12) = 0.0
cbzdbg2/22/99end
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,3):present flavor code of the projectile nucleon
c$$$        (hadron) I(a necleon of meson can be excired to its vector
c$$$        resonance).

c$$$        NFP(I,4):original flavor code of projectile
c$$$        nucleon(hadron)I.

c$$$        NFP(I,5):collision status of projectile nucleon(hadron) I.
c$$$        = 0 : suffered no collision.
c$$$        = 1 : suffered an elastic collision.
c$$$        = 2 : being the diffractive one in a single-diffractive
c$$$        collision.
c$$$        = 3 : became an excited string after an inelastic collision. 

c$$$        NFP(I,6): the total number of hard scattering associated
c$$$        with projectile nucleon(hadron) I. if NEP(I,6)<0, it can not
c$$$        produce jets any more due to energy conservation.

c$$$        NFP(I,10): to indicate whether the valence quarks of
c$$$        diquarks (anti_quarks) in projectile nucleon (hadron) I suffered
c$$$        a hard scattering.
c$$$        = 0 : has not suffered a hard scattering.
c$$$        = 1 : suffered one or more hard scatterings in current binary
c$$$        nucleon-nucleon collision.
c$$$        = -1: suffered one or more hard scatterings in previous binary
c$$$        nucleon-nucleon collisions.

c$$$         NFP(I,11): total number of interactions projectile nucleon
c$$$        (hadron) I has suffered so far.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        将弹核的原始味道代码和现在的味道代码(NFP(I,3)和NFP(I,4))
c$$$        都赋为质子的2212.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,3)和NFP(I,4)貌似和npart-xy的后两行有直接的关系。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        NFP(I,3)=IPP
        NFP(I,4)=IPP
        NFP(I,5)=0
        NFP(I,6)=0
        NFP(I,7)=0
        NFP(I,8)=0
        NFP(I,9)=0
        NFP(I,10)=0
        NFP(I,11)=0
        NPJ(I)=0
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        对于剩下的核子全部初始为中子。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IF(I.GT.ABS(IHNT2(2))) NFP(I,3)=2112
clin-12/2012 correct NN differential cross section in HIJING:
        IF(I.GT.ABS(IHNT2(2))) NFP(I,4)=2112
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9993, *) "IDQ = ", IDQ, " IDQQ = ", IDQQ
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        CALL ATTFLV(NFP(I,3),IDQ,IDQQ)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9993, *) "IDQ = ", IDQ, " IDQQ = ", IDQQ
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<

        NFP(I,1)=IDQ
        NFP(I,2)=IDQQ
        NFP(I,15)=-1
        IF(ABS(IDQ).GT.1000.OR.(ABS(IDQ*IDQQ).LT.100.AND.
     &                RANART(NSEED).LT.0.5)) NFP(I,15)=1
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      弹核中夸克和反夸克的质量。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PP(I,14)=ULMASS(IDQ)
        PP(I,15)=ULMASS(IDQQ)
100        CONTINUE
C
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    靶核的处理，和上面的处理方法类似。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 200 I=1,IHNT2(3)
        PT(I,1)=0.0
        PT(I,2)=0.0
        PT(I,3)=-SQRT(HINT1(1)**2/4.0-HINT1(9)**2)
        PT(I,4)=HINT1(1)/2.0
        PT(I,5)=HINT1(9)
        PT(I,6)=0.0
        PT(I,7)=0.0
        PT(I,8)=0.0
        PT(I,9)=0.0
        PT(I,10)=0.0
ctest OFF
cbzdbg2/22/99
        PT(I, 11) = 0.0
        PT(I, 12) = 0.0
cbzdbg2/22/99end
        NFT(I,3)=IPT
        NFT(I,4)=IPT
        NFT(I,5)=0
        NFT(I,6)=0
        NFT(I,7)=0
        NFT(I,8)=0
        NFT(I,9)=0
        NFT(I,10)=0
        NFT(I,11)=0
        NTJ(I)=0
        IF(I.GT.ABS(IHNT2(4))) NFT(I,3)=2112
clin-12/2012 correct NN differential cross section in HIJING:
        IF(I.GT.ABS(IHNT2(4))) NFT(I,4)=2112
        CALL ATTFLV(NFT(I,3),IDQ,IDQQ)
        NFT(I,1)=IDQ
        NFT(I,2)=IDQQ
        NFT(I,15)=1
        IF(ABS(IDQ).GT.1000.OR.(ABS(IDQ*IDQQ).LT.100.AND.
     &       RANART(NSEED).LT.0.5)) NFT(I,15)=-1
        PT(I,14)=ULMASS(IDQ)
        PT(I,15)=ULMASS(IDQQ)
200        CONTINUE
        RETURN
        END
C
C
C
