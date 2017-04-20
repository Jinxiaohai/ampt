c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        碰撞参数最大值和最小值以及坐标系
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        SUBROUTINE HIJING(FRAME,BMIN0,BMAX0)
cbz1/25/99
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        MAXPTN: CONSTANT NUMBER = 400001
c$$$        MAXSTR: CONSTANT NUMBER = 150001
c$$$        MAXSTR: CONSTANT NUMBER = 4001
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PARAMETER (MAXPTN=400001)
clin-4/20/01        PARAMETER (MAXSTR = 1600)
        PARAMETER (MAXSTR=150001)
cbz1/25/99end
clin-4/26/01:
        PARAMETER (MAXIDL=4001)
cbz1/31/99
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        声明一些双精度的物理量(三坐标,时间,三动量,能量和质量)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DOUBLE PRECISION  GX0, GY0, GZ0, FT0, PX0, PY0, PZ0, E0, XMASS0
        DOUBLE PRECISION  GX5, GY5, GZ5, FT5, PX5, PY5, PZ5, E5, XMASS5
        DOUBLE PRECISION  ATAUI, ZT1, ZT2, ZT3
        DOUBLE PRECISION  xnprod,etprod,xnfrz,etfrz,
     & dnprod,detpro,dnfrz,detfrz
clin-8/2015:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        这是三速度吗？
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DOUBLE PRECISION vxp0,vyp0,vzp0,xstrg0,ystrg0,xstrg,ystrg
cbz1/31/99end
        CHARACTER FRAME*8
        DIMENSION SCIP(300,300),RNIP(300,300),SJIP(300,300),JTP(3),
     &                        IPCOL(90000),ITCOL(90000)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
C
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        分别为弹核和靶核的三坐标YP(3,300), YT(3,300)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/hjcrdn/YP(3,300),YT(3,300)
cc      SAVE /hjcrdn/
clin-7/16/03 NINT is a intrinsic fortran function, rename it to NINTHJ
c        COMMON/HJGLBR/NELT,NINT,NELP,NINP
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        弹核和靶核核子碰撞的状态(弹性和非弹性)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HJGLBR/NELT,NINTHJ,NELP,NINP
cc      SAVE /HJGLBR/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        EATT:总的能量,用来检验能量守恒。
c$$$        JATT:the total number of jet-pairs in the current event.
c$$$        NATT:当前事件产生的稳定的粒子的总数。
c$$$        NT:靶核受伤的核子数。
c$$$        NP:弹核受伤的核子数。
c$$$        N0:       number of N-N, N-N_wounded, N_wounded-N, and 
c$$$        N01:      N_wounded-N_wounded collisions in the current
c$$$        N10:      events.(这几个变量似乎是受伤的核子)。
c$$$        N11:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
cc      SAVE /HMAIN1/
clin-4/26/01
c        COMMON/HMAIN2/KATT(130000,4),PATT(130000,4)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KATT两维数组的介绍！！！！！
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KATT(I,1): particle ID of the I_th produced particle. Users
c$$$        have to refer to JETSET7.2 for identifying particles with
c$$$        their ID's.(产生的第I个粒子的ID).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KATT(I,2): This is a code to identify the sources from which
c$$$        the particle comes.
c$$$        = 0 : projectile which has not interacted at all.
c$$$        = 1 : projectile nucleon (or hadron) which only suffers an
c$$$              elastic collision.
c$$$        = 2 : from a diffractive projectile nucleon (or hadron) in a
c$$$              single diffractive interaction.
c$$$        = 3 : from the fragmentation of a projectile string
c$$$              system(including gluon jets).
c$$$        = 10: target nucleon (or hadron) which has not interacted at
c$$$              all.
c$$$        = 11: target nucleon (or hadron) which only suffers an elastic
c$$$              collision.
c$$$        = 12: from a diffractive target nucleon (or hadron) in a single
c$$$              diffractive interaction.
c$$$        = 13: from the fragmentation of a target string system
c$$$              (including gluon jets).
c$$$        = 20: from scattered partons which form string systems
c$$$              themeselves.
c$$$        = 40: from direct production in the hard process (currently,
c$$$              only direct photons are included).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KATT(I,3): (I=1,...,NATT)line number of the parent particle.
c$$$        For finally produced or directly produced (not from the decay of
c$$$        another particle)particles, it is set to 0(The option to keep
c$$$        the information of all particles including the decayed ones is
c$$$        IHPR2(21)=1).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KATT(I,4): (I=1,...,NATT)status number of the particle.
c$$$        = 1 : finally of directly produced particles.
c$$$        = 11: particles which has already decayed.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PATT(I,1-4): (I=1,...,NATT) four-momentum(px, py, pz, E) of
c$$$        the produced particles.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
cc      SAVE /HMAIN2/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,1):(I=1,...,IAP)flavor code of the valence quark in
c$$$        projectile  nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,2):flavor code of diquark in projectile nucleon
c$$$       (anti_quark in projectile meson) I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,3):present flavor code of the projectile nucleon
c$$$        (hadron) I(a necleon of meson can be excired to its vector
c$$$        resonance).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,4):original flavor code of projectile
c$$$        nucleon(hadron)I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,5):collision status of projectile nucleon(hadron) I.
c$$$        = 0 : suffered no collision.
c$$$        = 1 : suffered an elastic collision.
c$$$        = 2 : being the diffractive one in a single-diffractive
c$$$        collision.
c$$$        = 3 : became an excited string after an inelastic collision. 
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,6): the total number of hard scattering associated
c$$$        with projectile nucleon(hadron) I. if NEP(I,6)<0, it can not
c$$$        produce jets any more due to energy conservation.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,10): to indicate whether the valence quarks of
c$$$        diquarks (anti_quarks) in projectile nucleon (hadron) I suffered
c$$$        a hard scattering.
c$$$        = 0 : has not suffered a hard scattering.
c$$$        = 1 : suffered one or more hard scatterings in current binary
c$$$        nucleon-nucleon collision.
c$$$        = -1: suffered one or more hard scatterings in previous binary
c$$$        nucleon-nucleon collisions.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,11): total number of interactions projectile nucleon
c$$$        (hadron) I has suffered so far.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,1),PP(I,2),PP(I,3),PP(I,4),PP(I,5):four momentum and
c$$$        the invariant mass (px,py,pz,E,M) of projectile nucleon(hadron)I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,6), PP(I,7): transverse momentum (px,py)of the valence
c$$$        quark in projectile nucleon(hadron)I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,8), PP(I,9):transverse momentum (px,py)of the diquark
c$$$        (anti_quark) in projectile nucleon(hadron)I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,10), PP(I,11), PP(I,12):three momentum (px,py,pz)
c$$$        transferred to the quark or diquark (anti_quark) in projectile
c$$$        nucleon (hadron) I from the last hard scattering.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,14):mass of the quark in projectile nucleon(hadron) I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,15):mass of the diquark (anti_quark) in projectile
c$$$        nucleon (hadron) I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
cc      SAVE /HSTRNG/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NPJ(I):(I=1,...,IAP)number of partons associated with
c$$$        projectile nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KFPJ(I,J):(I=1,...,IAP,J=1,...,NPJ(I))parton flavor code of
c$$$        the parton J associated with projectile nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PJPX(I,J),PJPY(I,J),PJPZ(I,J),PJPM(I,J):the four momentum
c$$$        and mass (px,py,pz,E,M) of parton J associated with the
c$$$        projectile nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NTJ(I):(I=1,...,IAT)number of partons associated with
c$$$        target nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KFTJ(I,J):(I=1,...,IAT,J=1,...,NTJ(I)):parton flavor code
c$$$        of the parton J associated with target nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PJTX(I,J),PJTY(I,J),PJTZ(I,J),PJTE(I,J),PJTM(I,J):
c$$$        the four momentum and mass (px,py,pz,E,M) of parton J associated
c$$$        with the target nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &                PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &                PJPM(300,500),NTJ(300),KFTJ(300,500),
     &                PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &                PJTE(300,500),PJTM(300,500)
cc      SAVE /HJJET1/
clin-4/2008
c        COMMON/HJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100),
c     &       K2SG(900,100),PXSG(900,100),PYSG(900,100),
c     &       PZSG(900,100),PESG(900,100),PMSG(900,100)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NSG:the total number of such string systems.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NJSG(I):(I=1,...,NSG)number of partons in the string system I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IASG(I,1),IASG(I,2): to specify which projectile and target
c$$$        nucleons produce string system I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        K1SG(I,J):(J=1,...,NJSG(I))color flow information of parton
c$$$        J in string system I (see JETSET 7.2 for detailed explanation).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        K2SG(I,J):(J=1,...,NJSG(I))flavor code of parton J in
c$$$        string system I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PXSG(I,J),PYSG(I,J),PZSG(I,J),PESG(I,J),PMSG(I,J):four
c$$$        momentum and mass (px,py,pz,E,M) of parton J in string system I
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &       K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &       PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
cc      SAVE /HJJET2/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NDR:total number of directly produced particles.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IADR(I,1),IADR(I,2):the sequence numbers of projectile and
c$$$        target nucleons which produce particle I during their
c$$$        interaction.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KFDR(I):the flavor code of directly produced particle I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PDR(I,1,...,5):four momentum and mass (px,py,pz,E,M)of
c$$$        particle I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HJJET4/NDR,IADR(MAXSTR,2),KFDR(MAXSTR),PDR(MAXSTR,5)
clin-4/2008:
c        common/xydr/rtdr(900,2)
        common/xydr/rtdr(MAXSTR,2)
cc      SAVE /HJJET4/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
C
        COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)   
cc      SAVE /LUJETS/
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
cc      SAVE /LUDAT1/
clin-9/29/03 changed name in order to distinguish from /prec2/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        changed name in order to distinguish from /prec2/(进行区分)
c$$$        参看下面下面下面的数据块进行区分。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON /ARPRC/ ITYPAR(MAXSTR),
     &       GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &       PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &       XMAR(MAXSTR)
ccbz11/11/98
c        COMMON /ARPRC/ ITYP(MAXSTR),
c     &     GX(MAXSTR), GY(MAXSTR), GZ(MAXSTR), FT(MAXSTR),
c     &     PX(MAXSTR), PY(MAXSTR), PZ(MAXSTR), EE(MAXSTR),
c     &     XM(MAXSTR)
cc      SAVE /ARPRC/
ccbz11/11/98end
cbz1/25/99
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        下面是储存什么信息的?????????????????
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON /PARA1/ MUL
cc      SAVE /PARA1/
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &     PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &     XMASS0(MAXPTN), ITYP0(MAXPTN)
cc      SAVE /prec1/
        COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &       PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &       XMASS5(MAXPTN), ITYP5(MAXPTN)
cc      SAVE /prec2/
        COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
cc      SAVE /ilist7/
        COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)
cc      SAVE /ilist8/
        COMMON /SREC1/ NSP, NST, NSI
cc      SAVE /SREC1/
        COMMON /SREC2/ATAUI(MAXSTR),ZT1(MAXSTR),ZT2(MAXSTR),ZT3(MAXSTR)
cc      SAVE /SREC2/
cbz1/25/99end
clin-2/25/00
        COMMON /frzout/ xnprod(30),etprod(30),xnfrz(30),etfrz(30),
     & dnprod(30),detpro(30),dnfrz(30),detfrz(30)
cc      SAVE /frzout/ 
clin-4/11/01 soft:
      common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
clin-4/25/01 soft3:
      DOUBLE PRECISION PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
cc      SAVE /SOFT/
clin-4/26/01 lepton and photon info:
        COMMON /NOPREC/ NNOZPC, ITYPN(MAXIDL),
     &       GXN(MAXIDL), GYN(MAXIDL), GZN(MAXIDL), FTN(MAXIDL),
     &       PXN(MAXIDL), PYN(MAXIDL), PZN(MAXIDL), EEN(MAXIDL),
     &       XMN(MAXIDL)
cc      SAVE /NOPREC/
clin-6/22/01:
        common /lastt/itimeh,bimp
cc      SAVE /lastt/
        COMMON /AREVT/ IAEVT, IARUN, MISS
        common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
clin-7/2011 ioscar value is needed:
        common /para7/ ioscar,nsmbbbar,nsmmeson
clin-2/2012 allow random orientation of reaction plane:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        考虑反应平面的随机方向问题,下面的两个变量都是以PHI角开头。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        common /phiHJ/iphirp,phiRP
clin-8/2015:
        common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
     1       xstrg0(MAXPTN),ystrg0(MAXPTN),
     2       xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)
        common /xiaohai/xiaohaiflag
        SAVE   
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(34): maximum radial coordinate for projectile nucleons
c$$$        to be given by the initialization program HIJSET.(弹核半径)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(35): maximum radial coordiante for target nucleons
c$$$        to be given by the initialization program HIJSET.(靶核半径)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(1):The mass number of the projectile nucleus.(弹核
c$$$        核子的质量数)(1 for a hadron).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(3):The mass number of the target nucleus.(靶核
c$$$        核子的质量数)(1 for a hadron).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        当弹核或者靶核的核子数都小于1时,设定BMAX的数值.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(31): the cross section sigma which characterizes the
c$$$        geometrical size of a nucleon(此处看文献上的等式).The default
c$$$        value is only for high-energy limit > 200GeV. At lower energies,
c$$$        a slight decrease which depends on energy is parametrized in the
c$$$        program. The default values of the two parameters HIPR1(30),
c$$$        HIPR1(31) are only for NN type interactions.for other kinds of
c$$$        projectile or target hadrons, users should change these values
c$$$        so that correct inelastic and total cross sections are obtained
c$$$        by the program.(散射截面和核子大小的关系).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(40): value of PION.(pion 的具体数值).
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9983,*)"HIPR1(34) = ",HIPR1(34),
     &      "HIPR1(35) = ", HIPR1(35)
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        BMAX=MIN(BMAX0,HIPR1(34)+HIPR1(35))
        BMIN=MIN(BMIN0,BMAX)
        IF(IHNT2(1).LE.1 .AND. IHNT2(3).LE.1) THEN
                BMIN=0.0
                BMAX=2.5*SQRT(HIPR1(31)*0.1/HIPR1(40))
        ENDIF
C                        ********HIPR1(31) is in mb =0.1fm**2
C*******THE FOLLOWING IS TO SELECT THE COORDINATIONS OF NUCLEONS 
C       BOTH IN PROJECTILE AND TARGET NUCLEAR( in fm)
C
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        弹核和靶核每个核子坐标的初始化问题.
c$$$        YP(1,KP):第KP个核子的X坐标.
c$$$        YP(2,KP):第KP个核子的Y坐标.
c$$$        YP(3,KP):第KP个核子的Z坐标.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        YP(1,1)=0.0
        YP(2,1)=0.0
        YP(3,1)=0.0
        IF(IHNT2(1).LE.1) GO TO 14
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        下面的循环给出初始弹核核子的坐标信息.
c$$$        此循环共计22行。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 10 KP=1,IHNT2(1)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     HIRND(1):这是一个随机的函数。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
5        R=HIRND(1)
        X=RANART(NSEED)
        CX=2.0*X-1.0
        SX=SQRT(1.0-CX*CX)
C                ********choose theta from uniform cos(theta) distr
        PHI=RANART(NSEED)*2.0*HIPR1(40)
C                ********choose phi form uniform phi distr 0 to 2*pi
        YP(1,KP)=R*SX*COS(PHI)
        YP(2,KP)=R*SX*SIN(PHI)
        YP(3,KP)=R*CX
        IF(HIPR1(29).EQ.0.0) GO TO 10
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        两个核子之间的距离不能小于0.4
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 8  KP2=1,KP-1
                DNBP1=(YP(1,KP)-YP(1,KP2))**2
                DNBP2=(YP(2,KP)-YP(2,KP2))**2
                DNBP3=(YP(3,KP)-YP(3,KP2))**2
                DNBP=DNBP1+DNBP2+DNBP3
                IF(DNBP.LT.HIPR1(29)*HIPR1(29)) GO TO 5
C                        ********two neighbors cannot be closer than 
C                                HIPR1(29)
8        CONTINUE
10        CONTINUE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      输出弹核的坐标信息
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
          if(xiaohaiflag .eq. 6) then
             DO 521 KP = 1, IHNT2(1)
                write(9987, *)YP(1, KP), "    ", YP(2, KP), "    ",
     &    YP(3, KP)
 521         CONTINUE
          endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
          
clin-1/27/03 Hulthen wavefn for deuteron borrowed from hijing1.382.f, 
c     but modified [divide by 2, & x(p)=-x(n)]: 
c     (Note: hijing1.383.f has corrected this bug in hijing1.382.f)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        氘核的初始化
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(IHNT2(1).EQ.2) then
           rnd1=max(RANART(NSEED),1.0e-20)
           rnd2=max(RANART(NSEED),1.0e-20)
           rnd3=max(RANART(NSEED),1.0e-20)
           R=-(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0
     &          +4.38*0.85*log(rnd3)/(4.38+0.85))
           X=RANART(NSEED)
           CX=2.0*X-1.0
           SX=SQRT(1.0-CX*CX)
           PHI=RANART(NSEED)*2.0*HIPR1(40)
c     R above is the relative distance between p & n in a deuteron:
           R=R/2.
           YP(1,1)=R*SX*COS(PHI)
           YP(2,1)=R*SX*SIN(PHI)
           YP(3,1)=R*CX
c     p & n has opposite coordinates in the deuteron frame:
           YP(1,2)=-YP(1,1)
           YP(2,2)=-YP(2,1)
           YP(3,2)=-YP(3,1)
        endif
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        下面的貌似还是给出弹核的坐标信息.好像弹核的核子坐标必须满足
c$$$        YP(3,I) <= YP(3,J),即纵向的坐标条件。但为什么按着z的数值从小
c$$$        到大进行一个纵向的排列
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 12 I=1,IHNT2(1)-1
        DO 12 J=I+1,IHNT2(1)
        IF(YP(3,I).GT.YP(3,J)) GO TO 12
        Y1=YP(1,I)
        Y2=YP(2,I)
        Y3=YP(3,I)
        YP(1,I)=YP(1,J)
        YP(2,I)=YP(2,J)
        YP(3,I)=YP(3,J)
        YP(1,J)=Y1
        YP(2,J)=Y2
        YP(3,J)=Y3
12        CONTINUE
C
C******************************
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        弹核和靶核每个核子坐标的初始化问题.
c$$$        YT(1,KT):第KT个核子的X坐标.
c$$$        YT(2,KT):第KT个核子的Y坐标.
c$$$        YT(3,KT):第KT个核子的Z坐标.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
14        YT(1,1)=0.0
        YT(2,1)=0.0
        YT(3,1)=0.0
        IF(IHNT2(3).LE.1) GO TO 24
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        弹核的核子坐标初始化,此段代码共计22行。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 20 KT=1,IHNT2(3)
15        R=HIRND(2)
        X=RANART(NSEED)
        CX=2.0*X-1.0
        SX=SQRT(1.0-CX*CX)
C                ********choose theta from uniform cos(theta) distr
        PHI=RANART(NSEED)*2.0*HIPR1(40)
C                ********chose phi form uniform phi distr 0 to 2*pi
        YT(1,KT)=R*SX*COS(PHI)
        YT(2,KT)=R*SX*SIN(PHI)
        YT(3,KT)=R*CX
        IF(HIPR1(29).EQ.0.0) GO TO 20
        DO 18  KT2=1,KT-1
                DNBT1=(YT(1,KT)-YT(1,KT2))**2
                DNBT2=(YT(2,KT)-YT(2,KT2))**2
                DNBT3=(YT(3,KT)-YT(3,KT2))**2
                DNBT=DNBT1+DNBT2+DNBT3
                IF(DNBT.LT.HIPR1(29)*HIPR1(29)) GO TO 15
C                        ********two neighbors cannot be closer than 
C                                HIPR1(29)
18        CONTINUE
20        CONTINUE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
          if(xiaohaiflag .eq. 6) then
             DO 520 KT = 1, IHNT2(3)
                write(9986, *)YT(1, KT), "    ", YT(2, KT), "    ",
     &    YT(3, KT)
 520         CONTINUE
          endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c
clin-1/27/03 Hulthen wavefn for deuteron borrowed from hijing1.382.f, 
c     but modified [divide by 2, & x(p)=-x(n)]:
        if(IHNT2(3).EQ.2) then
           rnd1=max(RANART(NSEED),1.0e-20)
           rnd2=max(RANART(NSEED),1.0e-20)
           rnd3=max(RANART(NSEED),1.0e-20)
           R=-(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0
     &          +4.38*0.85*log(rnd3)/(4.38+0.85))
           X=RANART(NSEED)
           CX=2.0*X-1.0
           SX=SQRT(1.0-CX*CX)
           PHI=RANART(NSEED)*2.0*HIPR1(40)
           R=R/2.
           YT(1,1)=R*SX*COS(PHI)
           YT(2,1)=R*SX*SIN(PHI)
           YT(3,1)=R*CX
           YT(1,2)=-YT(1,1)
           YT(2,2)=-YT(2,1)
           YT(3,2)=-YT(3,1)
        endif
c
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        下面的貌似还是给出靶核的坐标信息.好像靶核的核子坐标必须满足
c$$$        YT(3,I) <= YT(3,J),即纵向的坐标条件。但为什么按照z的纵向的数
c$$$        值从小到大进行排列
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 22 I=1,IHNT2(3)-1
        DO 22 J=I+1,IHNT2(3)
        IF(YT(3,I).LT.YT(3,J)) GO TO 22
        Y1=YT(1,I)
        Y2=YT(2,I)
        Y3=YT(3,I)
        YT(1,I)=YT(1,J)
        YT(2,I)=YT(2,J)
        YT(3,I)=YT(3,J)
        YT(1,J)=Y1
        YT(2,J)=Y2
        YT(3,J)=Y3
22        CONTINUE
C********************
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$          对MISS进行初始化MISS = 0, MISS 是初始化失败的次数。缺省值
c$$$          为1000.
c$$$          如果MISS大于maxmiss则进入死循环.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
24        MISS=-1
50        MISS=MISS+1
clin-6/2009
c        IF(MISS.GT.50) THEN
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        MISS应该是hijing程序此Subroutine重复的次数,如果重复的次数
c$$$        大于maxmiss时，就显示发生了无限的循环。其中maxmiss时从
c$$$        ampt.input读进去的，缺省值时1000.
c$$$        此条件语句共计四行。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IF(MISS.GT.maxmiss) THEN
           WRITE(6,*) 'infinite loop happened in  HIJING'
           STOP
        ENDIF
clin-4/30/01:
        itest=0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        EATT:总的能量,用来检验能量守恒。
c$$$        JATT:the total number of jet-pairs in the current event.
c$$$        NATT:当前事件产生的稳定的粒子的总数。
c$$$        NT:靶核受伤的核子数。
c$$$        NP:弹核受伤的核子数。
c$$$        N0:       number of N-N, N-N_wounded, N_wounded-N, and 
c$$$        N01:      N_wounded-N_wounded collisions in the current
c$$$        N10:      events.(这几个变量似乎是受伤的核子)。
c$$$        N11:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        NATT=0
        JATT=0
        EATT=0.0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        初始化核子的每个夸克和反夸克的详细的信息。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        CALL HIJINI
        NLOP=0
C                        ********Initialize for a new event
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     下面的这一堆初始化有重要的意义
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
60        NT=0
        NP=0
        N0=0
        N01=0
        N10=0
        N11=0
        NELT=0
        NINTHJ=0
        NELP=0
        NINP=0
        NSG=0
        NCOLT=0
C****        BB IS THE ABSOLUTE VALUE OF IMPACT PARAMETER,BB**2 IS 
C       RANDOMLY GENERATED AND ITS ORIENTATION IS RANDOMLY SET 
C       BY THE ANGLE PHI  FOR EACH COLLISION.******************
C
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        给出每个事件的碰撞参数和碰撞的方向。包含下面的四行。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        BB=SQRT(BMIN**2+RANART(NSEED)*(BMAX**2-BMIN**2))
cbz6/28/99 flow1
clin-2/2012:
        PHI=0.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    iphirp : 数值大小为0， 也就是说PHI=0.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(iphirp.eq.1) PHI=2.0*HIPR1(40)*RANART(NSEED)
        phiRP=phi
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9985, *)"PHIRP =   ", phiRP
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
cbz6/28/99 flow1 end
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$       BBX = BB, BBY = 0
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        BBX=BB*COS(PHI)
        BBY=BB*SIN(PHI)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(19):碰撞参数。
c$$$        HINT1(20):碰撞的方向。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9994, *) "BB = ", BB, "  PHI = ", PHI
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        HINT1(19)=BB
        HINT1(20)=PHI
C
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        下面的关于DO 70的循环完成弹性散射.
c$$$        此两层循环共计47行
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 70 JP=1,IHNT2(1)
        DO 70 JT=1,IHNT2(3)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        有空自己把下面代码推导下。!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
           SCIP(JP,JT)=-1.0
           B2=(YP(1,JP)+BBX-YT(1,JT))**2+(YP(2,JP)+BBY-YT(2,JT))**2
           R2=B2*HIPR1(40)/HIPR1(31)/0.1
C                ********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
           RRB1=MIN((YP(1,JP)**2+YP(2,JP)**2)
     &          /1.2**2/REAL(IHNT2(1))**0.6666667,1.0)
           RRB2=MIN((YT(1,JT)**2+YT(2,JT)**2)
     &          /1.2**2/REAL(IHNT2(3))**0.6666667,1.0)
           APHX1=HIPR1(6)*4.0/3.0*(IHNT2(1)**0.3333333-1.0)
     &           *SQRT(1.0-RRB1)
           APHX2=HIPR1(6)*4.0/3.0*(IHNT2(3)**0.3333333-1.0)
     &           *SQRT(1.0-RRB2)
           HINT1(18)=HINT1(14)-APHX1*HINT1(15)
     &                        -APHX2*HINT1(16)+APHX1*APHX2*HINT1(17)
           IF(IHPR2(14).EQ.0.OR.
     &          (IHNT2(1).EQ.1.AND.IHNT2(3).EQ.1)) THEN
              GS=1.0-EXP(-(HIPR1(30)+HINT1(18))*ROMG(R2)/HIPR1(31))
              RANTOT=RANART(NSEED)
              IF(RANTOT.GT.GS) GO TO 70
              GO TO 65
           ENDIF
           GSTOT0=2.0*(1.0-EXP(-(HIPR1(30)+HINT1(18))
     &             /HIPR1(31)/2.0*ROMG(0.0)))
           R2=R2/GSTOT0
           GS=1.0-EXP(-(HIPR1(30)+HINT1(18))/HIPR1(31)*ROMG(R2))
           GSTOT=2.0*(1.0-SQRT(1.0-GS))
           RANTOT=RANART(NSEED)*GSTOT0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           if(xiaohaiflag .eq. 6) then
              write(9994, *) "RANTOT : ", RANTOT, "  GSTOT = ", GSTOT
           endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           IF(RANTOT.GT.GSTOT) GO TO 70
           IF(RANTOT.GT.GS) THEN
              CALL HIJCSC(JP,JT)
              GO TO 70
C                        ********perform elastic collisions
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           if(xiaohaiflag .eq. 6) then
              write(9994, *) "RANTOT : ", RANTOT, "  GSTOT = ", GSTOT
           endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           ENDIF
 65           SCIP(JP,JT)=R2
           RNIP(JP,JT)=RANTOT
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  HINT1(18):有效的散射的截面。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
           SJIP(JP,JT)=HINT1(18)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  相互作用的碰撞次数,ipcol,itcol记录的是每个核子参加的是第几次碰撞
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
           NCOLT=NCOLT+1
           IPCOL(NCOLT)=JP
           ITCOL(NCOLT)=JT
70        CONTINUE
C                ********total number interactions proj and targ has
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
          do 517 i=1, NCOLT
             write(9962,*)IPCOL(I),"    ", ITCOL(I)
 517         continue
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
C                                suffered
clin-5/22/01 write impact parameter:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
             write(9994, *)"bimp = ",bimp,
     &        " nlop = ",nlop," ncolt = ",ncolt
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        bimp=bb
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    如果碰撞的次数小于0,那就重来呗。goto 60
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(6,*) '#impact parameter,nlop,ncolt=',bimp,nlop,ncolt
        IF(NCOLT.EQ.0) THEN
           NLOP=NLOP+1
           IF(NLOP.LE.20.OR.
     &           (IHNT2(1).EQ.1.AND.IHNT2(3).EQ.1)) GO TO 60
           RETURN
        ENDIF
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     对于较大的碰撞参数，选择重复事件，值到碰撞的发生。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
C               ********At large impact parameter, there maybe no
C                       interaction at all. For NN collision
C                       repeat the event until interaction happens
C
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  IHPR2(3):swich for triggered hard scattering with specified
c$$$        Pt >= HIPR1(10),记录硬散射在哪次碰撞中创建的。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
cxiaohai-17/4/20  下面的硬散射机制不发生。
c_______________________________________________________________________
        IF(IHPR2(3).NE.0) THEN
           NHARD=1+INT(RANART(NSEED)*(NCOLT-1)+0.5)
           NHARD=MIN(NHARD,NCOLT)
           JPHARD=IPCOL(NHARD)
           JTHARD=ITCOL(NHARD)
clin-6/2009 ctest off:
c           write(99,*) IAEVT,NHARD,NCOLT,JPHARD,JTHARD
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           write(9961, *) IAEVT,"  ",NHARD,"  ",
     &          NCOLT,"  ",JPHARD,"  ",JTHARD
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        ENDIF
c_______________________________________________________________________
C     
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    ihpr2(9):确保至少一对minijet产生。确定哪次碰撞产生minijet。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IF(IHPR2(9).EQ.1) THEN
                NMINI=1+INT(RANART(NSEED)*(NCOLT-1)+0.5)
                NMINI=MIN(NMINI,NCOLT)
                JPMINI=IPCOL(NMINI)
                JTMINI=ITCOL(NMINI)
        ENDIF
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           write(9994, *) NMINI,"  ",JPMINI,"  ",JTMINI
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
C                ********Specifying the location of the hard and
C                        minijet if they are enforced by user
C
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,1):(I=1,...,IAP)flavor code of the valence quark in
c$$$        projectile  nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,2):flavor code of diquark in projectile nucleon
c$$$       (anti_quark in projectile meson) I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,3):present flavor code of the projectile nucleon
c$$$        (hadron) I(a necleon of meson can be excired to its vector
c$$$        resonance).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,4):original flavor code of projectile
c$$$        nucleon(hadron)I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,5):collision status of projectile nucleon(hadron) I.
c$$$        = 0 : suffered no collision.
c$$$        = 1 : suffered an elastic collision.
c$$$        = 2 : being the diffractive one in a single-diffractive
c$$$        collision.
c$$$        = 3 : became an excited string after an inelastic collision. 
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,6): the total number of hard scattering associated
c$$$        with projectile nucleon(hadron) I. if NEP(I,6)<0, it can not
c$$$        produce jets any more due to energy conservation.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,10): to indicate whether the valence quarks of
c$$$        diquarks (anti_quarks) in projectile nucleon (hadron) I suffered
c$$$        a hard scattering.
c$$$        = 0 : has not suffered a hard scattering.
c$$$        = 1 : suffered one or more hard scatterings in current binary
c$$$        nucleon-nucleon collision.
c$$$        = -1: suffered one or more hard scatterings in previous binary
c$$$        nucleon-nucleon collisions.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,11): total number of interactions projectile nucleon
c$$$        (hadron) I has suffered so far.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    遍历初始的两个反应核子
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 200 JP=1,IHNT2(1)
        DO 200 JT=1,IHNT2(3)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
              write(9994, *)"SCIP(JP, JT) = ", SCIP(JP, JT)
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    NFP(JP, 5)：核子的碰撞状态。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IF(SCIP(JP,JT).EQ.-1.0) GO
     $          TO 200
                NFP(JP,11)=NFP(JP,11)+1
                NFT(JT,11)=NFT(JT,11)+1
        IF(NFP(JP,5).LE.1 .AND. NFT(JT,5).GT.1) THEN
                NP=NP+1
                N01=N01+1
        ELSE IF(NFP(JP,5).GT.1 .AND. NFT(JT,5).LE.1) THEN
                NT=NT+1
                N10=N10+1
        ELSE IF(NFP(JP,5).LE.1 .AND. NFT(JT,5).LE.1) THEN
                NP=NP+1
                NT=NT+1
                N0=N0+1
        ELSE IF(NFP(JP,5).GT.1 .AND. NFT(JT,5).GT.1) THEN
                N11=N11+1
        ENDIF
        JOUT=0
        NFP(JP,10)=0
        NFT(JT,10)=0
C*****************************************************************
        IF(IHPR2(8).EQ.0 .AND. IHPR2(3).EQ.0) GO TO 160
C                ********When IHPR2(8)=0 no jets are produced
        IF(NFP(JP,6).LT.0 .OR. NFT(JT,6).LT.0) GO TO 160
C                ********jets can not be produced for (JP,JT)
C                        because not enough energy avaible for 
C                                JP or JT 
        R2=SCIP(JP,JT)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    HINT1(18):有效的散射截面。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        HINT1(18)=SJIP(JP,JT)
        TT=ROMG(R2)*HINT1(18)/HIPR1(31)
        TTS=HIPR1(30)*ROMG(R2)/HIPR1(31)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9994, *)"TT = ", TT, "TTS = ", TTS
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        NJET=0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9994, *)"IHPR2(3) = ", IHPR2(3), "  JPHARD = ", JPHARD,
     &      "   JTHARD = ", JTHARD
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  貌似满足jet的初始化条件。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
cxiaohai17/4/20 下面的硬散射不发生。        
c_______________________________________________________________________
        IF(IHPR2(3).NE.0 .AND.
     $       JP.EQ.JPHARD .AND. JT.EQ.JTHARD) THEN
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           write(9960,*)"JET WAS HAPPENED !"
           write(9960,*)"IHPR2(3) = ", IHPR2(3),
     &          "JP = ", JP, "JPHARD = ", JPHARD,
     &          "JT = ", JT, "JTHARD = ", JTHARD
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           CALL JETINI(JP,JT,1)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           if(xiaohaiflag .eq. 6) then
              write(9994, *) "JFLG : ", JFLG
           endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    HIJHRD(JP,JT,0,JFLG,0)硬散射的处理过程.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
           CALL HIJHRD(JP,JT,0,JFLG,0)
           HINT1(26)=HINT1(47)
           HINT1(27)=HINT1(48)
           HINT1(28)=HINT1(49)
           HINT1(29)=HINT1(50)
           HINT1(36)=HINT1(67)
           HINT1(37)=HINT1(68)
           HINT1(38)=HINT1(69)
           HINT1(39)=HINT1(70)
C
           IF(ABS(HINT1(46)).GT.HIPR1(11).AND.JFLG.EQ.2) NFP(JP,7)=1
           IF(ABS(HINT1(56)).GT.HIPR1(11).AND.JFLG.EQ.2) NFT(JT,7)=1
           IF(MAX(ABS(HINT1(46)),ABS(HINT1(56))).GT.HIPR1(11).AND.
     &                                JFLG.GE.3) IASG(NSG,3)=1
           IHNT2(9)=IHNT2(14)
           IHNT2(10)=IHNT2(15)
           DO 105 I05=1,5
              HINT1(20+I05)=HINT1(40+I05)
              HINT1(30+I05)=HINT1(50+I05)
 105           CONTINUE
clin-6/2009 ctest off:
c           write(99,*) jp,jt,IHPR2(3),HIPR1(10),njet,
c     1          ihnt2(9),hint1(21),hint1(22),hint1(23),
c     2          ihnt2(10),hint1(31),hint1(32),hint1(33)
c           write(99,*) ' '
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
                  write(9990,*) jp, jt, IHPR2(3), HIPR1(10), njet,
     &             ihnt2(9), hint1(21), hint1(22), hint1(23),
     &             ihnt2(10), hint1(31), hint1(32), hint1(33)
                  write(9990,*)
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           JOUT=1
           IF(IHPR2(8).EQ.0) GO TO 160
           RRB1=MIN((YP(1,JP)**2+YP(2,JP)**2)/1.2**2
     &                /REAL(IHNT2(1))**0.6666667,1.0)
           RRB2=MIN((YT(1,JT)**2+YT(2,JT)**2)/1.2**2
     &                /REAL(IHNT2(3))**0.6666667,1.0)
           APHX1=HIPR1(6)*4.0/3.0*(IHNT2(1)**0.3333333-1.0)
     &           *SQRT(1.0-RRB1)
           APHX2=HIPR1(6)*4.0/3.0*(IHNT2(3)**0.3333333-1.0)
     &           *SQRT(1.0-RRB2)
           HINT1(65)=HINT1(61)-APHX1*HINT1(62)
     &                        -APHX2*HINT1(63)+APHX1*APHX2*HINT1(64)
           TTRIG=ROMG(R2)*HINT1(65)/HIPR1(31)
           NJET=-1
C                ********subtract the trigger jet from total number
C                        of jet production  to be done since it has
C                                already been produced here
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     统计产生的jet的数量信息。并从中间去激发的jet数量。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
           XR1=-ALOG(EXP(-TTRIG)+RANART(NSEED)*(1.0-EXP(-TTRIG)))
 106           NJET=NJET+1
           XR1=XR1-ALOG(max(RANART(NSEED),1.0e-20))
           IF(XR1.LT.TTRIG) GO TO 106
           XR=0.0
 107           NJET=NJET+1
           XR=XR-ALOG(max(RANART(NSEED),1.0e-20))
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
              write(9989,*)"njet", njet
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           IF(XR.LT.TT-TTRIG) GO TO 107
           NJET=NJET-1
           GO TO 112
        ENDIF
c_______________________________________________________________________        
C                ********create a hard interaction with specified P_T
c                                 when IHPR2(3)>0
        IF(IHPR2(9).EQ.1.AND.JP.EQ.JPMINI.AND.JT.EQ.JTMINI) GO TO 110
C                ********create at least one pair of mini jets 
C                        when IHPR2(9)=1
C
clin-4/15/2010 changed .LT. to .LE. to avoid problem when two sides are equal; 
c     this problem may lead to a jet production when there should be none and 
c     crash the run; crashes at low energies were reported by P. Bhaduri.
c        IF(IHPR2(8).GT.0 .AND.RNIP(JP,JT).LT.EXP(-TT)*
c     &                (1.0-EXP(-TTS))) GO TO 160
        IF(IHPR2(8).GT.0 .AND.RNIP(JP,JT).LE.EXP(-TT)*
     &                 (1.0-EXP(-TTS))) GO TO 160
c
C                ********this is the probability for no jet production
110        XR=-ALOG(EXP(-TT)+RANART(NSEED)*(1.0-EXP(-TT)))
111        NJET=NJET+1
        XR=XR-ALOG(max(RANART(NSEED),1.0e-20))
        IF(XR.LT.TT) GO TO 111
112        NJET=MIN(NJET,IHPR2(8))
        IF(IHPR2(8).LT.0)  NJET=ABS(IHPR2(8))
C                ******** Determine number of mini jet production
C
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           write(9988, *) "JP = ", jp, "jt = ",jt, "NJET = ", NJET
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     对其中的两个核子碰撞产生的jet进行初始化。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 150 ijet=1,NJET
           CALL JETINI(JP,JT,0)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    HIJHRD(JP,JT,0,JFLG,0)硬散射的处理过程.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
           CALL HIJHRD(JP,JT,JOUT,JFLG,1)
C                ********JFLG=1 jets valence quarks, JFLG=2 with 
C                        gluon jet, JFLG=3 with q-qbar prod for
C                        (JP,JT). If JFLG=0 jets can not be produced 
C                        this time. If JFLG=-1, error occured abandon
C                        this event. JOUT is the total hard scat for
C                        (JP,JT) up to now.
           IF(JFLG.EQ.0) GO TO 160
           IF(JFLG.LT.0) THEN
              IF(IHPR2(10).NE.0) WRITE(6,*) 'error occured in HIJHRD'
              GO TO 50
           ENDIF
           JOUT=JOUT+1
           IF(ABS(HINT1(46)).GT.HIPR1(11).AND.JFLG.EQ.2) NFP(JP,7)=1
           IF(ABS(HINT1(56)).GT.HIPR1(11).AND.JFLG.EQ.2) NFT(JT,7)=1
           IF(MAX(ABS(HINT1(46)),ABS(HINT1(56))).GT.HIPR1(11).AND.
     &                        JFLG.GE.3) IASG(NSG,3)=1
C                ******** jet with PT>HIPR1(11) will be quenched
 150        CONTINUE
 160        CONTINUE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      完成JP和JT之间的软散射过程。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        CALL HIJSFT(JP,JT,JOUT,IERROR)
        IF(IERROR.NE.0) THEN
           IF(IHPR2(10).NE.0) WRITE(6,*) 'error occured in HIJSFT'
           GO TO 50
        ENDIF
C
C                ********conduct soft scattering between JP and JT
        JATT=JATT+JOUT
200        CONTINUE
c
c**************************
c
clin-6/2009 write out initial minijet information:
clin-2/2012:
c           call minijet_out(BB)
           call minijet_out(BB,phiRP)
           if(pttrig.gt.0.and.ntrig.eq.0) goto 50
clin-4/2012 
clin-6/2009 write out initial transverse positions of initial nucleons:
c           write(94,*) IAEVT,MISS,IHNT2(1),IHNT2(3)
        DO 201 JP=1,IHNT2(1)
clin-6/2009:
c           write(94,203) YP(1,JP)+0.5*BB, YP(2,JP), JP, NFP(JP,5)
clin-2/2012:
c       write(94,203) YP(1,JP)+0.5*BB, YP(2,JP), JP, NFP(JP,5),yp(3,jp)
clin-4/2012:
c           write(94,203) YP(1,JP)+0.5*BB*cos(phiRP), 
c     1 YP(2,JP)+0.5*BB*sin(phiRP), JP, NFP(JP,5),yp(3,jp)
           IF(NFP(JP,5).GT.2) THEN
              NINP=NINP+1
           ELSE IF(NFP(JP,5).EQ.2.OR.NFP(JP,5).EQ.1) THEN
              NELP=NELP+1
           ENDIF
 201    continue
        DO 202 JT=1,IHNT2(3)
clin-6/2009 target nucleon # has a minus sign for distinction from projectile:
c           write(94,203) YT(1,JT)-0.5*BB, YT(2,JT), -JT, NFT(JT,5)
clin-2/2012:
c       write(94,203) YT(1,JT)-0.5*BB, YT(2,JT), -JT, NFT(JT,5),yt(3,jt)
clin-4/2012:
c           write(94,203) YT(1,JT)-0.5*BB*cos(phiRP), 
c     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt)
           IF(NFT(JT,5).GT.2) THEN
              NINTHJ=NINTHJ+1
           ELSE IF(NFT(JT,5).EQ.2.OR.NFT(JT,5).EQ.1) THEN
              NELT=NELT+1
           ENDIF
 202    continue
c 203    format(f10.3,1x,f10.3,2(1x,I5))
c 203    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3)
c     
c*******************************
C********perform jet quenching for jets with PT>HIPR1(11)**********
        IF((IHPR2(8).NE.0.OR.IHPR2(3).NE.0).AND.IHPR2(4).GT.0.AND.
     &                        IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
                DO 271 I=1,IHNT2(1)
                        IF(NFP(I,7).EQ.1) CALL QUENCH(I,1)
271                CONTINUE
                DO 272 I=1,IHNT2(3)
                        IF(NFT(I,7).EQ.1) CALL QUENCH(I,2)
272                CONTINUE
                DO 273 ISG=1,NSG
                        IF(IASG(ISG,3).EQ.1) CALL QUENCH(ISG,3)
273                CONTINUE
        ENDIF
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    处理strings的方式的不同。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$clin*****4/09/01-soft1, default way of treating strings:
        if(isoft.eq.1) then
c$$$clin-4/16/01 allow fragmentation:
c$$$           isflag=1
c$$$cbz1/25/99
c$$$c.....transfer data from HIJING to ZPC
c$$$        NSP = IHNT2(1)
c$$$        NST = IHNT2(3)
c$$$        NSI = NSG
c$$$        ISTR = 0
c$$$        NPAR = 0
c$$$        DO 1008 I = 1, IHNT2(1)
c$$$           ISTR = ISTR + 1
c$$$           DO 1007 J = 1, NPJ(I)
c$$$cbz1/27/99
c$$$c.....for now only consider gluon cascade
c$$$              IF (KFPJ(I, J) .EQ. 21) THEN
c$$$cbz1/27/99end
c$$$              NPAR = NPAR + 1
c$$$              LSTRG0(NPAR) = ISTR
c$$$              LPART0(NPAR) = J
c$$$              ITYP0(NPAR) = KFPJ(I, J)
c$$$cbz6/28/99 flow1
c$$$clin-7/20/01 add dble or sngl to make precisions consistent
c$$$c              GX0(NPAR) = YP(1, I)
c$$$clin-2/2012:
c$$$c              GX0(NPAR) = dble(YP(1, I) + 0.5 * BB)
c$$$              GX0(NPAR) = dble(YP(1, I)+0.5*BB*cos(phiRP))
c$$$cbz6/28/99 flow1 end
c$$$c              GY0(NPAR) = dble(YP(2, I))
c$$$              GY0(NPAR) = dble(YP(2, I)+0.5*BB*sin(phiRP))
c$$$              GZ0(NPAR) = 0d0
c$$$              FT0(NPAR) = 0d0
c$$$              PX0(NPAR) = dble(PJPX(I, J))
c$$$              PY0(NPAR) = dble(PJPY(I, J))
c$$$              PZ0(NPAR) = dble(PJPZ(I, J))
c$$$              XMASS0(NPAR) = dble(PJPM(I, J))
c$$$c              E0(NPAR) = dble(PJPE(I, J))
c$$$              E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
c$$$     1             +PZ0(NPAR)**2+XMASS0(NPAR)**2)
c$$$clin-7/20/01-end
c$$$cbz1/27/99
c$$$c.....end gluon selection
c$$$              END IF
c$$$cbz1/27/99end
c$$$ 1007      CONTINUE
c$$$ 1008   CONTINUE
c$$$        DO 1010 I = 1, IHNT2(3)
c$$$           ISTR = ISTR + 1
c$$$           DO 1009 J = 1, NTJ(I)
c$$$cbz1/27/99
c$$$c.....for now only consider gluon cascade
c$$$              IF (KFTJ(I, J) .EQ. 21) THEN
c$$$cbz1/27/99end
c$$$              NPAR = NPAR + 1
c$$$              LSTRG0(NPAR) = ISTR
c$$$              LPART0(NPAR) = J
c$$$              ITYP0(NPAR) = KFTJ(I, J)
c$$$cbz6/28/99 flow1
c$$$clin-7/20/01 add dble or sngl to make precisions consistent
c$$$c              GX0(NPAR) = YT(1, I)
c$$$clin-2/2012:
c$$$c              GX0(NPAR) = dble(YT(1, I) - 0.5 * BB)
c$$$              GX0(NPAR) = dble(YT(1, I)-0.5*BB*cos(phiRP))
c$$$cbz6/28/99 flow1 end
c$$$c              GY0(NPAR) = dble(YT(2, I))
c$$$              GY0(NPAR) = dble(YT(2, I)-0.5*BB*sin(phiRP))
c$$$              GZ0(NPAR) = 0d0
c$$$              FT0(NPAR) = 0d0
c$$$              PX0(NPAR) = dble(PJTX(I, J))
c$$$              PY0(NPAR) = dble(PJTY(I, J))
c$$$              PZ0(NPAR) = dble(PJTZ(I, J))
c$$$              XMASS0(NPAR) = dble(PJTM(I, J))
c$$$c              E0(NPAR) = dble(PJTE(I, J))
c$$$              E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
c$$$     1             +PZ0(NPAR)**2+XMASS0(NPAR)**2)
c$$$cbz1/27/99
c$$$c.....end gluon selection
c$$$              END IF
c$$$cbz1/27/99end
c$$$ 1009      CONTINUE
c$$$ 1010   CONTINUE
c$$$        DO 1012 I = 1, NSG
c$$$           ISTR = ISTR + 1
c$$$           DO 1011 J = 1, NJSG(I)
c$$$cbz1/27/99
c$$$c.....for now only consider gluon cascade
c$$$              IF (K2SG(I, J) .EQ. 21) THEN
c$$$cbz1/27/99end
c$$$              NPAR = NPAR + 1
c$$$              LSTRG0(NPAR) = ISTR
c$$$              LPART0(NPAR) = J
c$$$              ITYP0(NPAR) = K2SG(I, J)
c$$$clin-7/20/01 add dble or sngl to make precisions consistent:
c$$$              GX0(NPAR) = 0.5d0 * 
c$$$     1             dble(YP(1, IASG(I, 1)) + YT(1, IASG(I, 2)))
c$$$              GY0(NPAR) = 0.5d0 * 
c$$$     2             dble(YP(2, IASG(I, 1)) + YT(2, IASG(I, 2)))
c$$$              GZ0(NPAR) = 0d0
c$$$              FT0(NPAR) = 0d0
c$$$              PX0(NPAR) = dble(PXSG(I, J))
c$$$              PY0(NPAR) = dble(PYSG(I, J))
c$$$              PZ0(NPAR) = dble(PZSG(I, J))
c$$$              XMASS0(NPAR) = dble(PMSG(I, J))
c$$$c              E0(NPAR) = dble(PESG(I, J))
c$$$              E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
c$$$     1             +PZ0(NPAR)**2+XMASS0(NPAR)**2)
c$$$cbz1/27/99
c$$$c.....end gluon selection
c$$$              END IF
c$$$cbz1/27/99end
c$$$ 1011      CONTINUE
c$$$ 1012   CONTINUE
c$$$        MUL = NPAR
c$$$cbz2/4/99
c$$$        CALL HJANA1
c$$$cbz2/4/99end
c$$$clin-6/2009:
c$$$        if(ioscar.eq.3) WRITE (95, *) IAEVT, mul
c$$$c.....call ZPC for parton cascade
c$$$        CALL ZPCMN
c$$$c     write out parton and wounded nucleon information to ana/zpc1.mom:
c$$$clin-6/2009:
c$$$c        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
c$$$        WRITE (14, 395) IAEVT, MISS, MUL, bimp, NELP,NINP,NELT,NINTHJ
c$$$        DO 1013 I = 1, MUL
c$$$cc           WRITE (14, 411) PX5(I), PY5(I), PZ5(I), ITYP5(I),
c$$$c     &        XMASS5(I), E5(I)
c$$$           if(dmax1(abs(GX5(I)),abs(GY5(I)),abs(GZ5(I)),abs(FT5(I)))
c$$$     1          .lt.9999) then
c$$$              write(14,210) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
c$$$     1             GX5(I), GY5(I), GZ5(I), FT5(I)
c$$$           else
c$$$c     change format for large numbers:
c$$$              write(14,211) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
c$$$     1             GX5(I), GY5(I), GZ5(I), FT5(I)
c$$$           endif
c$$$ 1013   CONTINUE
 210    format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
 211    format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
 395    format(3I8,f10.4,4I5)
c$$$clin-4/09/01:
c$$$        itest=itest+1
c$$$c 411    FORMAT(1X, 3F10.3, I6, 2F10.3)
c$$$cbz3/19/99 end
c$$$clin-5/2009 ctest off:
c$$$c        call frztm(1,1)
c$$$c.....transfer data back from ZPC to HIJING
c$$$        DO 1014 I = 1, MUL
c$$$           IF (LSTRG1(I) .LE. NSP) THEN
c$$$              NSTRG = LSTRG1(I)
c$$$              NPART = LPART1(I)
c$$$              KFPJ(NSTRG, NPART) = ITYP5(I)
c$$$clin-7/20/01 add dble or sngl to make precisions consistent
c$$$              PJPX(NSTRG, NPART) = sngl(PX5(I))
c$$$              PJPY(NSTRG, NPART) = sngl(PY5(I))
c$$$              PJPZ(NSTRG, NPART) = sngl(PZ5(I))
c$$$              PJPE(NSTRG, NPART) = sngl(E5(I))
c$$$              PJPM(NSTRG, NPART) = sngl(XMASS5(I))
c$$$           ELSE IF (LSTRG1(I) .LE. NSP + NST) THEN
c$$$              NSTRG = LSTRG1(I) - NSP
c$$$              NPART = LPART1(I)
c$$$              KFTJ(NSTRG, NPART) = ITYP5(I)
c$$$              PJTX(NSTRG, NPART) = sngl(PX5(I))
c$$$              PJTY(NSTRG, NPART) = sngl(PY5(I))
c$$$              PJTZ(NSTRG, NPART) = sngl(PZ5(I))
c$$$              PJTE(NSTRG, NPART) = sngl(E5(I))
c$$$              PJTM(NSTRG, NPART) = sngl(XMASS5(I))
c$$$           ELSE
c$$$              NSTRG = LSTRG1(I) - NSP - NST
c$$$              NPART = LPART1(I)
c$$$              K2SG(NSTRG, NPART) = ITYP5(I)
c$$$              PXSG(NSTRG, NPART) = sngl(PX5(I))
c$$$              PYSG(NSTRG, NPART) = sngl(PY5(I))
c$$$              PZSG(NSTRG, NPART) = sngl(PZ5(I))
c$$$              PESG(NSTRG, NPART) = sngl(E5(I))
c$$$              PMSG(NSTRG, NPART) = sngl(XMASS5(I))
c$$$           END IF
c$$$ 1014   CONTINUE
c$$$cbz1/25/99end
c$$$cbz2/4/99
c$$$        CALL HJANA2
c$$$cbz2/4/99end
c$$$clin*****4/09/01-soft2, put q+dq+X in strings into ZPC:
c$$$        elseif(isoft.eq.2) then
c$$$        NSP = IHNT2(1)
c$$$        NST = IHNT2(3)
c$$$clin-4/27/01:
c$$$        NSI = NSG
c$$$        NPAR=0
c$$$        ISTR=0
c$$$C
c$$$clin  No fragmentation to hadrons, only on parton level, 
c$$$c     and transfer minijet and string data from HIJING to ZPC:
c$$$        MSTJ(1)=0
c$$$clin-4/12/01 forbid soft radiation before ZPC to avoid small-mass strings,
c$$$c     and forbid jet order reversal before ZPC to avoid unphysical flavors:
c$$$        IHPR2(1)=0
c$$$        isflag=0
c$$$        IF(IHPR2(20).NE.0) THEN
c$$$           DO 320 NTP=1,2
c$$$              DO 310 jjtp=1,IHNT2(2*NTP-1)
c$$$                 ISTR = ISTR + 1
c$$$c change: do gluon kink only once: either here or in fragmentation.
c$$$                 CALL HIJFRG(jjtp,NTP,IERROR)
c$$$c                 call lulist(1)
c$$$                 if(NTP.eq.1) then
c$$$c 354                continue
c$$$                    NPJ(jjtp)=MAX0(N-2,0)
c$$$clin-4/12/01:                    NPJ(jjtp)=MAX0(ipartn-2,0)
c$$$                 else
c$$$c 355                continue
c$$$                    NTJ(jjtp)=MAX0(N-2,0)
c$$$clin-4/12/01:                    NTJ(jjtp)=MAX0(ipartn-2,0)
c$$$                 endif
c$$$                 do 300 ii=1,N
c$$$                 NPAR = NPAR + 1
c$$$                 LSTRG0(NPAR) = ISTR
c$$$                 LPART0(NPAR) = II
c$$$                 ITYP0(NPAR) = K(II,2)
c$$$                 GZ0(NPAR) = 0d0
c$$$                 FT0(NPAR) = 0d0
c$$$clin-7/20/01 add dble or sngl to make precisions consistent
c$$$                 PX0(NPAR) = dble(P(II,1))
c$$$                 PY0(NPAR) = dble(P(II,2))
c$$$                 PZ0(NPAR) = dble(P(II,3))
c$$$                 XMASS0(NPAR) = dble(P(II,5))
c$$$c                 E0(NPAR) = dble(P(II,4))
c$$$                 E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
c$$$     1                +PZ0(NPAR)**2+XMASS0(NPAR)**2)
c$$$                 IF (NTP .EQ. 1) THEN
c$$$clin-7/20/01 add dble or sngl to make precisions consistent
c$$$clin-2/2012:
c$$$c                    GX0(NPAR) = dble(YP(1, jjtp)+0.5 * BB)
c$$$c                    GY0(NPAR) = dble(YP(2, jjtp))
c$$$                    GX0(NPAR) = dble(YP(1, jjtp)+0.5*BB*cos(phiRP))
c$$$                    GY0(NPAR) = dble(YP(2, jjtp)+0.5*BB*sin(phiRP))
c$$$                    IITYP=ITYP0(NPAR)
c$$$                    nstrg=LSTRG0(NPAR)
c$$$                    if(IITYP.eq.2112.or.IITYP.eq.2212) then
c$$$                    elseif((IITYP.eq.1.or.IITYP.eq.2).and.
c$$$     1 (II.eq.1.or.II.eq.N)) then
c$$$                       PP(nstrg,6)=sngl(PX0(NPAR))
c$$$                       PP(nstrg,7)=sngl(PY0(NPAR))
c$$$                       PP(nstrg,14)=sngl(XMASS0(NPAR))
c$$$                    elseif((IITYP.eq.1103.or.IITYP.eq.2101
c$$$     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
c$$$     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
c$$$     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
c$$$     4 .and.(II.eq.1.or.II.eq.N)) then
c$$$                       PP(nstrg,8)=sngl(PX0(NPAR))
c$$$                       PP(nstrg,9)=sngl(PY0(NPAR))
c$$$                       PP(nstrg,15)=sngl(XMASS0(NPAR))
c$$$                    else
c$$$                       NPART = LPART0(NPAR)-1
c$$$                       KFPJ(NSTRG, NPART) = ITYP0(NPAR)
c$$$                       PJPX(NSTRG, NPART) = sngl(PX0(NPAR))
c$$$                       PJPY(NSTRG, NPART) = sngl(PY0(NPAR))
c$$$                       PJPZ(NSTRG, NPART) = sngl(PZ0(NPAR))
c$$$                       PJPE(NSTRG, NPART) = sngl(E0(NPAR))
c$$$                       PJPM(NSTRG, NPART) = sngl(XMASS0(NPAR))
c$$$                    endif
c$$$                 ELSE
c$$$clin-2/2012:
c$$$c                    GX0(NPAR) = dble(YT(1, jjtp)-0.5 * BB)
c$$$c                    GY0(NPAR) = dble(YT(2, jjtp)) 
c$$$                    GX0(NPAR) = dble(YT(1, jjtp)-0.5*BB*cos(phiRP))
c$$$                    GY0(NPAR) = dble(YT(2, jjtp)-0.5*BB*sin(phiRP))
c$$$                    IITYP=ITYP0(NPAR)
c$$$                    nstrg=LSTRG0(NPAR)-NSP
c$$$                    if(IITYP.eq.2112.or.IITYP.eq.2212) then
c$$$                    elseif((IITYP.eq.1.or.IITYP.eq.2).and.
c$$$     1 (II.eq.1.or.II.eq.N)) then
c$$$                       PT(nstrg,6)=sngl(PX0(NPAR))
c$$$                       PT(nstrg,7)=sngl(PY0(NPAR))
c$$$                       PT(nstrg,14)=sngl(XMASS0(NPAR))
c$$$                    elseif((IITYP.eq.1103.or.IITYP.eq.2101
c$$$     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
c$$$     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
c$$$     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
c$$$     4 .and.(II.eq.1.or.II.eq.N)) then
c$$$                       PT(nstrg,8)=sngl(PX0(NPAR))
c$$$                       PT(nstrg,9)=sngl(PY0(NPAR))
c$$$                       PT(nstrg,15)=sngl(XMASS0(NPAR))
c$$$                    else
c$$$                       NPART = LPART0(NPAR)-1
c$$$                       KFTJ(NSTRG, NPART) = ITYP0(NPAR)
c$$$                       PJTX(NSTRG, NPART) = sngl(PX0(NPAR))
c$$$                       PJTY(NSTRG, NPART) = sngl(PY0(NPAR))
c$$$                       PJTZ(NSTRG, NPART) = sngl(PZ0(NPAR))
c$$$                       PJTE(NSTRG, NPART) = sngl(E0(NPAR))
c$$$                       PJTM(NSTRG, NPART) = sngl(XMASS0(NPAR))
c$$$                    endif
c$$$                 END IF
c$$$ 300          continue
c$$$ 310          continue
c$$$ 320       continue
c$$$           DO 330 ISG=1,NSG
c$$$              ISTR = ISTR + 1
c$$$              CALL HIJFRG(ISG,3,IERROR)
c$$$c              call lulist(2)
c$$$c
c$$$              NJSG(ISG)=N
c$$$c
c$$$              do 1001 ii=1,N
c$$$                 NPAR = NPAR + 1
c$$$                 LSTRG0(NPAR) = ISTR
c$$$                 LPART0(NPAR) = II
c$$$                 ITYP0(NPAR) = K(II,2)
c$$$                 GX0(NPAR)=0.5d0*
c$$$     1                dble(YP(1,IASG(ISG,1))+YT(1,IASG(ISG,2)))
c$$$                 GY0(NPAR)=0.5d0*
c$$$     2                dble(YP(2,IASG(ISG,1))+YT(2,IASG(ISG,2)))
c$$$                 GZ0(NPAR) = 0d0
c$$$                 FT0(NPAR) = 0d0
c$$$                 PX0(NPAR) = dble(P(II,1))
c$$$                 PY0(NPAR) = dble(P(II,2))
c$$$                 PZ0(NPAR) = dble(P(II,3))
c$$$                 XMASS0(NPAR) = dble(P(II,5))
c$$$c                 E0(NPAR) = dble(P(II,4))
c$$$                 E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
c$$$     1                +PZ0(NPAR)**2+XMASS0(NPAR)**2)
c$$$ 1001         continue
c$$$ 330       continue
c$$$        endif
c$$$        MUL = NPAR
c$$$cbz2/4/99
c$$$        CALL HJANA1
c$$$cbz2/4/99end
c$$$clin-6/2009:
c$$$        if(ioscar.eq.3) WRITE (95, *) IAEVT, mul
c$$$c.....call ZPC for parton cascade
c$$$        CALL ZPCMN
c$$$cbz3/19/99
c$$$clin-6/2009:
c$$$c        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
c$$$        WRITE (14, 395) IAEVT, MISS, MUL, bimp, NELP,NINP,NELT,NINTHJ
c$$$        itest=itest+1
c$$$        DO 1015 I = 1, MUL
c$$$c           WRITE (14, 311) PX5(I), PY5(I), PZ5(I), ITYP5(I),
c$$$c     &        XMASS5(I), E5(I)
c$$$clin-4/2012 write parton freeze-out position in zpc.dat for this test scenario:
c$$$c           WRITE (14, 312) PX5(I), PY5(I), PZ5(I), ITYP5(I),
c$$$c     &        XMASS5(I), E5(I),LSTRG1(I), LPART1(I)
c$$$           if(dmax1(abs(GX5(I)),abs(GY5(I)),abs(GZ5(I)),abs(FT5(I)))
c$$$     1          .lt.9999) then
c$$$              write(14,210) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
c$$$     1             GX5(I), GY5(I), GZ5(I), FT5(I)
c$$$           else
c$$$              write(14,211) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
c$$$     1             GX5(I), GY5(I), GZ5(I), FT5(I)
c$$$           endif
c$$$c
c$$$ 1015   CONTINUE
c$$$c 311    FORMAT(1X, 3F10.4, I6, 2F10.4)
c$$$c 312    FORMAT(1X, 3F10.3, I6, 2F10.3,1X,I6,1X,I3)
c$$$cbz3/19/99 end
c$$$clin-5/2009 ctest off:
c$$$c        call frztm(1,1)
c$$$clin-4/13/01 initialize four momenta and invariant mass of strings after ZPC:
c$$$        do 1004 nmom=1,5
c$$$           do 1002 nstrg=1,nsp
c$$$              PP(nstrg,nmom)=0.
c$$$ 1002      continue
c$$$           do 1003 nstrg=1,nst
c$$$              PT(nstrg,nmom)=0.
c$$$ 1003      continue
c$$$ 1004   continue
c$$$clin-4/13/01-end
c$$$        DO 1005 I = 1, MUL
c$$$           IITYP=ITYP5(I)
c$$$           IF (LSTRG1(I) .LE. NSP) THEN
c$$$              NSTRG = LSTRG1(I)
c$$$c     nucleons without interactions:
c$$$              if(IITYP.eq.2112.or.IITYP.eq.2212) then
c$$$clin-7/20/01 add dble or sngl to make precisions consistent
c$$$                 PP(nstrg,1)=sngl(PX5(I))
c$$$                 PP(nstrg,2)=sngl(PY5(I))
c$$$                 PP(nstrg,3)=sngl(PZ5(I))
c$$$                 PP(nstrg,4)=sngl(E5(I))
c$$$                 PP(nstrg,5)=sngl(XMASS5(I))
c$$$c     valence quark:
c$$$              elseif((IITYP.eq.1.or.IITYP.eq.2).and.
c$$$     1 (LPART1(I).eq.1.or.LPART1(I).eq.(NPJ(NSTRG)+2))) then
c$$$                 PP(nstrg,6)=sngl(PX5(I))
c$$$                 PP(nstrg,7)=sngl(PY5(I))
c$$$                 PP(nstrg,14)=sngl(XMASS5(I))
c$$$                 PP(nstrg,1)=PP(nstrg,1)+sngl(PX5(I))
c$$$                 PP(nstrg,2)=PP(nstrg,2)+sngl(PY5(I))
c$$$                 PP(nstrg,3)=PP(nstrg,3)+sngl(PZ5(I))
c$$$                 PP(nstrg,4)=PP(nstrg,4)+sngl(E5(I))
c$$$                 PP(nstrg,5)=sqrt(PP(nstrg,4)**2-PP(nstrg,1)**2
c$$$     1                -PP(nstrg,2)**2-PP(nstrg,3)**2)
c$$$c     diquark:
c$$$              elseif((IITYP.eq.1103.or.IITYP.eq.2101
c$$$     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
c$$$     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
c$$$     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
c$$$     4 .and.(LPART1(I).eq.1.or.LPART1(I).eq.(NPJ(NSTRG)+2))) then
c$$$                 PP(nstrg,8)=sngl(PX5(I))
c$$$                 PP(nstrg,9)=sngl(PY5(I))
c$$$                 PP(nstrg,15)=sngl(XMASS5(I))
c$$$                 PP(nstrg,1)=PP(nstrg,1)+sngl(PX5(I))
c$$$                 PP(nstrg,2)=PP(nstrg,2)+sngl(PY5(I))
c$$$                 PP(nstrg,3)=PP(nstrg,3)+sngl(PZ5(I))
c$$$                 PP(nstrg,4)=PP(nstrg,4)+sngl(E5(I))
c$$$                 PP(nstrg,5)=sqrt(PP(nstrg,4)**2-PP(nstrg,1)**2
c$$$     1                -PP(nstrg,2)**2-PP(nstrg,3)**2)
c$$$c     partons in projectile or target strings:
c$$$              else
c$$$                 NPART = LPART1(I)-1
c$$$                 KFPJ(NSTRG, NPART) = ITYP5(I)
c$$$                 PJPX(NSTRG, NPART) = sngl(PX5(I))
c$$$                 PJPY(NSTRG, NPART) = sngl(PY5(I))
c$$$                 PJPZ(NSTRG, NPART) = sngl(PZ5(I))
c$$$                 PJPE(NSTRG, NPART) = sngl(E5(I))
c$$$                 PJPM(NSTRG, NPART) = sngl(XMASS5(I))
c$$$              endif
c$$$           ELSE IF (LSTRG1(I) .LE. NSP + NST) THEN
c$$$              NSTRG = LSTRG1(I) - NSP
c$$$              if(IITYP.eq.2112.or.IITYP.eq.2212) then
c$$$                 PT(nstrg,1)=sngl(PX5(I))
c$$$                 PT(nstrg,2)=sngl(PY5(I))
c$$$                 PT(nstrg,3)=sngl(PZ5(I))
c$$$                 PT(nstrg,4)=sngl(E5(I))
c$$$                 PT(nstrg,5)=sngl(XMASS5(I))
c$$$              elseif((IITYP.eq.1.or.IITYP.eq.2).and.
c$$$     1 (LPART1(I).eq.1.or.LPART1(I).eq.(NTJ(NSTRG)+2))) then
c$$$                 PT(nstrg,6)=sngl(PX5(I))
c$$$                 PT(nstrg,7)=sngl(PY5(I))
c$$$                 PT(nstrg,14)=sngl(XMASS5(I))
c$$$                 PT(nstrg,1)=PT(nstrg,1)+sngl(PX5(I))
c$$$                 PT(nstrg,2)=PT(nstrg,2)+sngl(PY5(I))
c$$$                 PT(nstrg,3)=PT(nstrg,3)+sngl(PZ5(I))
c$$$                 PT(nstrg,4)=PT(nstrg,4)+sngl(E5(I))
c$$$                 PT(nstrg,5)=sqrt(PT(nstrg,4)**2-PT(nstrg,1)**2
c$$$     1                -PT(nstrg,2)**2-PT(nstrg,3)**2)
c$$$              elseif((IITYP.eq.1103.or.IITYP.eq.2101
c$$$     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
c$$$     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
c$$$     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
c$$$     4 .and.(LPART1(I).eq.1.or.LPART1(I).eq.(NTJ(NSTRG)+2))) then
c$$$                 PT(nstrg,8)=sngl(PX5(I))
c$$$                 PT(nstrg,9)=sngl(PY5(I))
c$$$                 PT(nstrg,15)=sngl(XMASS5(I))
c$$$                 PT(nstrg,1)=PT(nstrg,1)+sngl(PX5(I))
c$$$                 PT(nstrg,2)=PT(nstrg,2)+sngl(PY5(I))
c$$$                 PT(nstrg,3)=PT(nstrg,3)+sngl(PZ5(I))
c$$$                 PT(nstrg,4)=PT(nstrg,4)+sngl(E5(I))
c$$$                 PT(nstrg,5)=sqrt(PT(nstrg,4)**2-PT(nstrg,1)**2
c$$$     1                -PT(nstrg,2)**2-PT(nstrg,3)**2)
c$$$              else
c$$$                 NPART = LPART1(I)-1
c$$$                 KFTJ(NSTRG, NPART) = ITYP5(I)
c$$$                 PJTX(NSTRG, NPART) = sngl(PX5(I))
c$$$                 PJTY(NSTRG, NPART) = sngl(PY5(I))
c$$$                 PJTZ(NSTRG, NPART) = sngl(PZ5(I))
c$$$                 PJTE(NSTRG, NPART) = sngl(E5(I))
c$$$                 PJTM(NSTRG, NPART) = sngl(XMASS5(I))
c$$$              endif
c$$$           ELSE
c$$$              NSTRG = LSTRG1(I) - NSP - NST
c$$$              NPART = LPART1(I)
c$$$              K2SG(NSTRG, NPART) = ITYP5(I)
c$$$              PXSG(NSTRG, NPART) = sngl(PX5(I))
c$$$              PYSG(NSTRG, NPART) = sngl(PY5(I))
c$$$              PZSG(NSTRG, NPART) = sngl(PZ5(I))
c$$$              PESG(NSTRG, NPART) = sngl(E5(I))
c$$$              PMSG(NSTRG, NPART) = sngl(XMASS5(I))
c$$$           END IF
c$$$ 1005   CONTINUE
c$$$cbz1/25/99end
c$$$clin-4/09/01  turn on fragmentation with soft radiation 
c$$$c     and jet order reversal to form hadrons after ZPC:
c$$$        MSTJ(1)=1
c$$$        IHPR2(1)=1
c$$$        isflag=1
c$$$clin-4/13/01 allow small mass strings (D=1.5GeV):
c$$$        HIPR1(1)=0.94
c$$$cbz2/4/99
c$$$        CALL HJANA2
c$$$cbz2/4/99end
c$$$clin-4/19/01-soft3, fragment strings, then convert hadrons to partons 
c$$$c     and input to ZPC:
        elseif(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
clin-4/24/01 normal fragmentation first:
        isflag=0
c        write(99,*) 'IAEVT,NSG,NDR=',IAEVT,NSG,NDR
        IF(IHPR2(20).NE.0) THEN
           DO 560 ISG=1,NSG
                CALL HIJFRG(ISG,3,IERROR)
C
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
 551                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  551
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
C                ******** boost back to lab frame(if it was in)
C
                nsbstR=0
                DO 560 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 560
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=20
                   KATT(NATT,4)=K(I,1)
c     from Yasushi, to avoid violation of array limits:
c                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
clin-4/2008 to avoid out-of-bound error in K():
c                   IF(K(I,3).EQ.0 .OR. 
c     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
c                      KATT(NATT,3)=0
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
clin-4/2008-end
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
C       ****** identify the mother particle
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
                   GXAR(NATT) = 0.5 * (YP(1, IASG(ISG, 1)) +
     &                YT(1, IASG(ISG, 2)))
                   GYAR(NATT) = 0.5 * (YP(2, IASG(ISG, 1)) +
     &                YT(2, IASG(ISG, 2)))
                   GZAR(NATT) = 0.
                   FTAR(NATT) = 0.
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
clin-8/2015: record hadron information, to be used for its constituent partons:
                   xstrg0(NATT)=dble(GXAR(NATT))
                   ystrg0(NATT)=dble(GYAR(NATT))
                   istrg0(NATT)=ISG
c                   write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT),
c     1                  K(I,2),P(I, 1),P(I, 2),P(I, 3)
cbz11/11/98end
 560            CONTINUE
C                ********Fragment the q-qbar jets systems *****
C
           JTP(1)=IHNT2(1)
           JTP(2)=IHNT2(3)
           DO 600 NTP=1,2
           DO 600 jjtp=1,JTP(NTP)
                CALL HIJFRG(jjtp,NTP,IERROR)
C
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
 581                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  581
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
C                ******** boost back to lab frame(if it was in)
C
                NFTP=NFP(jjtp,5)
                IF(NTP.EQ.2) NFTP=10+NFT(jjtp,5)
                nsbstR=0
                DO 590 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 590
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=NFTP
                   KATT(NATT,4)=K(I,1)
c                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
clin-4/2008
c                   IF(K(I,3).EQ.0 .OR.
c     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
c                      KATT(NATT,3)=0
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
clin-4/2008-end
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
C       ****** identify the mother particle
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
                   IF (NTP .EQ. 1) THEN
clin-2/2012:
c                      GXAR(NATT) = YP(1, jjtp)+0.5 * BB
c                      GYAR(NATT) = YP(2, jjtp)
                      GXAR(NATT) = YP(1, jjtp)+0.5*BB*cos(phiRP)
                      GYAR(NATT) = YP(2, jjtp)+0.5*BB*sin(phiRP)
                   ELSE
clin-2/2012:
c                      GXAR(NATT) = YT(1, jjtp)-0.5 * BB
c                      GYAR(NATT) = YT(2, jjtp)
                      GXAR(NATT) = YT(1, jjtp)-0.5*BB*cos(phiRP)
                      GYAR(NATT) = YT(2, jjtp)-0.5*BB*sin(phiRP)
                   END IF
                   GZAR(NATT) = 0.
                   FTAR(NATT) = 0.
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
clin-8/2015: record hadron information, to be used for its constituent partons:
                   xstrg0(NATT)=dble(GXAR(NATT))
                   ystrg0(NATT)=dble(GYAR(NATT))
c     String ID is separated for projectile/target strings:
                   istrg0(NATT)=NTP*10000+jjtp
c              if(N.eq.nsbst.and.(K(I,2).eq.2112.or.K(I,2).eq.2212)) then
c                      write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT)
c     1                  ,K(I,2),P(I, 1),P(I, 2),P(I, 3),'spectator'
c                   else
c                      write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT)
c     1                  ,K(I,2),P(I, 1),P(I, 2),P(I, 3)
c                   endif
cbz11/11/98end
 590                CONTINUE 
 600           CONTINUE
C     ********Fragment the q-qq related string systems
        ENDIF
clin-4/2008 check for zero NDR value:
        if(NDR.ge.1) then
c
        DO 650 I=1,NDR
                NATT=NATT+1
                KATT(NATT,1)=KFDR(I)
                KATT(NATT,2)=40
                KATT(NATT,3)=0
                PATT(NATT,1)=PDR(I,1)
                PATT(NATT,2)=PDR(I,2)
                PATT(NATT,3)=PDR(I,3)
                PATT(NATT,4)=PDR(I,4)
                EATT=EATT+PDR(I,4)
clin-11/11/03     set direct photons positions and time at formation:
                GXAR(NATT) = rtdr(I,1)
                GYAR(NATT) = rtdr(I,2)
                GZAR(NATT) = 0.
                FTAR(NATT) = 0.
                ITYPAR(NATT) =KATT(NATT,1) 
                PXAR(NATT) = PATT(NATT,1)
                PYAR(NATT) = PATT(NATT,2)
                PZAR(NATT) = PATT(NATT,3)
                PEAR(NATT) = PATT(NATT,4)
                XMAR(NATT) = PDR(I,5)
 650        CONTINUE
clin-4/2008:
         endif
clin-6/2009
         call embedHighPt
c
        CALL HJANA1
clin-4/19/01 convert hadrons to partons for ZPC (with GX0 given):
        call htop
clin-7/03/01 move up, used in zpstrg (otherwise not set and incorrect):
        nsp=0
        nst=0
        nsg=natt
        NSI=NSG
clin-7/03/01-end
clin-6/2009:
        if(ioscar.eq.3) WRITE (95, *) IAEVT, mul
c.....call ZPC for parton cascade
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9980,*)"call zpcmn...."
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        CALL ZPCMN
clin-6/2009:
c        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
        WRITE (14, 395) IAEVT, MISS, MUL, bimp, NELP,NINP,NELT,NINTHJ
        itest=itest+1
        DO 1016 I = 1, MUL
c           WRITE (14, 511) PX5(I), PY5(I), PZ5(I), ITYP5(I),
c     &        XMASS5(I), E5(I)
clin-4/2012 write parton freeze-out position in zpc.dat 
c     for string melting version:
c           WRITE (14, 512) ITYP5(I), PX5(I), PY5(I), PZ5(I), 
c     &        XMASS5(I), LSTRG1(I), LPART1(I), FT5(I)
           if(dmax1(abs(GX5(I)),abs(GY5(I)),abs(GZ5(I)),abs(FT5(I)))
     1          .lt.9999) then
              write(14,210) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           else
              write(14,211) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           endif
c
 1016   CONTINUE
c 511    FORMAT(1X, 3F10.4, I6, 2F10.4)
c 512    FORMAT(I6,4(1X,F10.3),1X,I6,1X,I3,1X,F10.3)
c 513    FORMAT(1X, 4F10.4)
clin-5/2009 ctest off:
c        call frztm(1,1)
clin  save data after ZPC for fragmentation purpose:
c.....transfer data back from ZPC to HIJING
        DO 1018 I = 1, MAXSTR
           DO 1017 J = 1, 3
              K1SGS(I, J) = 0
              K2SGS(I, J) = 0
              PXSGS(I, J) = 0d0
              PYSGS(I, J) = 0d0
              PZSGS(I, J) = 0d0
              PESGS(I, J) = 0d0
              PMSGS(I, J) = 0d0
              GXSGS(I, J) = 0d0
              GYSGS(I, J) = 0d0
              GZSGS(I, J) = 0d0
              FTSGS(I, J) = 0d0
 1017      CONTINUE
 1018   CONTINUE
        DO 1019 I = 1, MUL
           IITYP=ITYP5(I)
           NSTRG = LSTRG1(I)
           NPART = LPART1(I)
           K2SGS(NSTRG, NPART) = ITYP5(I)
           PXSGS(NSTRG, NPART) = PX5(I)
           PYSGS(NSTRG, NPART) = PY5(I)
           PZSGS(NSTRG, NPART) = PZ5(I)
           PMSGS(NSTRG, NPART) = XMASS5(I)
clin-7/20/01 E5(I) does no include the finite parton mass XMASS5(I), 
c     so define it anew:
c           PESGS(NSTRG, NPART) = E5(I)
c           if(abs(PZ5(i)/E5(i)).gt.0.9999999d0) 
c     1          write(91,*) 'a',PX5(i),PY5(i),XMASS5(i),PZ5(i),E5(i)
           E5(I)=dsqrt(PX5(I)**2+PY5(I)**2+PZ5(I)**2+XMASS5(I)**2)
           PESGS(NSTRG, NPART) = E5(I)
c           if(abs(PZ5(i)/E5(i)).gt.0.9999999d0) 
c     1          write(91,*) 'b: new E5(I)=',E5(i)
clin-7/20/01-end
           GXSGS(NSTRG, NPART) = GX5(I)
           GYSGS(NSTRG, NPART) = GY5(I)
           GZSGS(NSTRG, NPART) = GZ5(I)
           FTSGS(NSTRG, NPART) = FT5(I)
 1019   CONTINUE
        CALL HJANA2
clin-4/19/01-end
        endif
clin-4/09/01-end
C
C**************fragment all the string systems in the following*****
C
C********nsbst is where particle information starts
C********nsbstR+1 is the number of strings in fragmentation
C********the number of strings before a line is stored in K(I,4)
C********IDSTR is id number of the string system (91,92 or 93)
C
clin-4/30/01 convert partons to hadrons after ZPC:
        if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
           NATT=0
           EATT=0.
           call ptoh
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      nnozpc是什么意思????????????
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           write(9972,*)"the value of nnozpc is ===>  ", nnozpc
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           do 1006 I=1,nnozpc
              NATT=NATT+1
              KATT(NATT,1)=ITYPN(I)
              PATT(NATT,1)=PXN(I)
              PATT(NATT,2)=PYN(I)
              PATT(NATT,3)=PZN(I)
              PATT(NATT,4)=EEN(I)
              EATT=EATT+EEN(I)
              GXAR(NATT)=GXN(I)
              GYAR(NATT)=GYN(I)
              GZAR(NATT)=GZN(I)
              FTAR(NATT)=FTN(I)
              ITYPAR(NATT)=ITYPN(I)
              PXAR(NATT)=PXN(I)
              PYAR(NATT)=PYN(I)
              PZAR(NATT)=PZN(I)
              PEAR(NATT)=EEN(I)
              XMAR(NATT)=XMN(I)
 1006      continue
           goto 565
        endif
clin-4/30/01-end        
        IF(IHPR2(20).NE.0) THEN
           DO 360 ISG=1,NSG
                CALL HIJFRG(ISG,3,IERROR)
                IF(MSTU(24).NE.0 .OR.IERROR.GT.0) THEN
                   MSTU(24)=0
                   MSTU(28)=0
                   IF(IHPR2(10).NE.0) THEN
c                      call lulist(2)
                      WRITE(6,*) 'error occured ISG, repeat the event'
                  write(6,*) ISG
                   ENDIF
                   GO TO 50
                ENDIF
C                        ********Check errors
C
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
351                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  351
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
C
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
C                ******** boost back to lab frame(if it was in)
C
                nsbstR=0
                DO 360 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 360
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=20
                   KATT(NATT,4)=K(I,1)
c                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
clin-4/2008:
c                   IF(K(I,3).EQ.0 .OR. 
c     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
c                      KATT(NATT,3)=0
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
clin-4/2008-end
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
C       ****** identify the mother particle
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
cbz11/11/98
cbz1/25/99
c                   GXAR(NATT) = 0.5 * (YP(1, IASG(ISG, 1)) +
c     &                YT(1, IASG(ISG, 2)))
c                   GYAR(NATT) = 0.5 * (YP(2, IASG(ISG, 1)) +
c     &                YT(2, IASG(ISG, 2)))
                   LSG = NSP + NST + ISG
                   GXAR(NATT) = sngl(ZT1(LSG))
                   GYAR(NATT) = sngl(ZT2(LSG))
                   GZAR(NATT) = sngl(ZT3(LSG))
                   FTAR(NATT) = sngl(ATAUI(LSG))
cbz1/25/99end
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
cbz11/11/98end
360           CONTINUE
C                ********Fragment the q-qbar jets systems *****
C
           JTP(1)=IHNT2(1)
           JTP(2)=IHNT2(3)
           DO 400 NTP=1,2
           DO 400 jjtp=1,JTP(NTP)
                CALL HIJFRG(jjtp,NTP,IERROR)
                IF(MSTU(24).NE.0 .OR. IERROR.GT.0) THEN
                   MSTU(24)=0
                   MSTU(28)=0
                   IF(IHPR2(10).NE.0) THEN
c                  call lulist(2)
                  WRITE(6,*) 'error occured P&T, repeat the event'
                  WRITE(6,*) NTP,jjtp
clin-6/2009 when this happens, the event will be repeated, 
c     and another record for the same event number will be written into
c     zpc.dat, zpc.res, minijet-initial-beforePropagation.dat,
c     parton-initial-afterPropagation.dat, parton-after-coalescence.dat, 
c     and parton-collisionsHistory.dat. 
                   ENDIF
                   GO TO 50
                ENDIF
C                        ********check errors
C
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
381                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  381
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
C                ******** boost back to lab frame(if it was in)
C
                NFTP=NFP(jjtp,5)
                IF(NTP.EQ.2) NFTP=10+NFT(jjtp,5)
                nsbstR=0
                DO 390 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 390
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=NFTP
                   KATT(NATT,4)=K(I,1)
c                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
clin-4/2008:
c                   IF(K(I,3).EQ.0 .OR. 
c     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
c                      KATT(NATT,3)=0
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
clin-4/2008-end
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
C       ****** identify the mother particle
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
cbz11/11/98
cbz1/25/99
c                   IF (NTP .EQ. 1) THEN
c                      GXAR(NATT) = YP(1, jjtp)
c                   ELSE
c                      GXAR(NATT) = YT(1, jjtp)
c                   END IF
c                   IF (NTP .EQ. 1) THEN
c                      GYAR(NATT) = YP(2, jjtp)
c                   ELSE
c                      GYAR(NATT) = YT(2, jjtp)
c                   END IF
                   IF (NTP .EQ. 1) THEN
                      LSG = jjtp
                   ELSE
                      LSG = jjtp + NSP
                   END IF
                   GXAR(NATT) = sngl(ZT1(LSG))
                   GYAR(NATT) = sngl(ZT2(LSG))
                   GZAR(NATT) = sngl(ZT3(LSG))
                   FTAR(NATT) = sngl(ATAUI(LSG))
cbz1/25/99end
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
cbz11/11/98end
390                CONTINUE 
400           CONTINUE
C     ********Fragment the q-qq related string systems
        ENDIF
        DO 450 I=1,NDR
           NATT=NATT+1
           KATT(NATT,1)=KFDR(I)
           KATT(NATT,2)=40
           KATT(NATT,3)=0
           PATT(NATT,1)=PDR(I,1)
           PATT(NATT,2)=PDR(I,2)
           PATT(NATT,3)=PDR(I,3)
           PATT(NATT,4)=PDR(I,4)
           EATT=EATT+PDR(I,4)
clin-11/11/03     set direct photons positions and time at formation:
           GXAR(NATT) = rtdr(I,1)
           GYAR(NATT) = rtdr(I,2)
           GZAR(NATT) = 0.
           FTAR(NATT) = 0.
           ITYPAR(NATT) =KATT(NATT,1) 
           PXAR(NATT) = PATT(NATT,1)
           PYAR(NATT) = PATT(NATT,2)
           PZAR(NATT) = PATT(NATT,3)
           PEAR(NATT) = PATT(NATT,4)
           XMAR(NATT) = PDR(I,5)
 450    CONTINUE
C                        ********store the direct-produced particles
C
clin-4/19/01 soft3:
 565    continue
        DENGY=EATT/(IHNT2(1)*HINT1(6)+IHNT2(3)*HINT1(7))-1.0
        IF(ABS(DENGY).GT.HIPR1(43).AND.IHPR2(20).NE.0
     &     .AND.IHPR2(21).EQ.0) THEN
         IF(IHPR2(10).NE.0) 
     &        WRITE(6,*) 'Energy not conserved, repeat the event'
c                call lulist(1)
         write(6,*) 'violated:EATT(GeV),NATT,B(fm)=',EATT,NATT,bimp
         GO TO 50
        ENDIF
        write(6,*) 'satisfied:EATT(GeV),NATT,B(fm)=',EATT,NATT,bimp
        write(6,*) ' '
c





        
clin-4/2012 write out initial transverse positions of initial nucleons:

c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  输出初始的核子的坐标的信息
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(xiaohaiflag .eq. 6) then
           write(9984,*)"PHIRP =  ", phiRP
        endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     因为采用的是固定反应平面的方式，所以只对X进行平移。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(94,*) IAEVT,MISS,IHNT2(1),IHNT2(3),bimp
        DO JP=1,IHNT2(1)
clin-12/2012 write out present and original flavor code of nucleons:
c           write(94,243) YP(1,JP)+0.5*BB*cos(phiRP), 
c     1 YP(2,JP)+0.5*BB*sin(phiRP), JP, NFP(JP,5),yp(3,jp)
           write(94,243) YP(1,JP)+0.5*BB*cos(phiRP), 
     1 YP(2,JP)+0.5*BB*sin(phiRP),JP, NFP(JP,5),yp(3,jp),
     2 NFP(JP,3),NFP(JP,4)
        ENDDO
        DO JT=1,IHNT2(3)
c target nucleon # has a minus sign for distinction from projectile:
clin-12/2012 write out present and original flavor code of nucleons:
c           write(94,243) YT(1,JT)-0.5*BB*cos(phiRP), 
c     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt)
           write(94,243) YT(1,JT)-0.5*BB*cos(phiRP), 
     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt),
     2 NFT(JT,3),NFT(JT,4)
        ENDDO
clin-12/2012 write out present and original flavor code of nucleons:
c 243    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3)
 243    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3,2(1x,I5))
clin-4/2012-end
        RETURN
        END
C
C
C
