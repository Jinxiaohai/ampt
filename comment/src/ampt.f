c.....driver program for A Multi-Phase Transport model
      PROGRAM AMPT
c
      double precision xmp, xmu, alpha, rscut2, cutof2, dshadow
      double precision smearp,smearh,dpcoal,drcoal,ecritl
      CHARACTER FRAME*8, PROJ*8, TARG*8
      character*25 amptvn
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON /HPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /AROUT/ IOUT
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /smearz/smearp,smearh
      COMMON/RNDF77/NSEED
      common/anim/nevent,isoft,isflag,izpc
c     parton coalescence radii in case of string melting:
      common /coal/dpcoal,drcoal,ecritl
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
c     initialization value for parton cascade:
      common /para2/ xmp, xmu, alpha, rscut2, cutof2
      common /para7/ ioscar,nsmbbbar,nsmmeson
      common /para8/ idpert,npertd,idxsec
      common /rndm3/ iseedp
c     initialization value for hadron cascade:
      COMMON /RUN/ NUM
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      common/oscar1/iap,izp,iat,izt
      common/oscar2/FRAME,amptvn
      common/resdcy/NSAV,iksdcy
clin-4/2012-6/2009:
c      common/phidcy/iphidcy
      common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
      common/embed/iembed,nsembd,pxqembd,pyqembd,xembd,yembd,
     1     psembd,tmaxembd,phidecomp
clin-7/2009:
      common/cmsflag/dshadow,ishadow
clin-2/2012 allow random orientation of reaction plane:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        初始的反映平面的方位角
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      common /phiHJ/iphirp,phiRP
      EXTERNAL HIDATA, PYDATA, LUDATA, ARDATA, PPBDAT, zpcbdt
      common /xiaohai/xiaohaiflag
      SAVE   
c****************
      OPEN (24, FILE = 'input.ampt', STATUS = 'UNKNOWN')
      OPEN (12, FILE = 'ana/version', STATUS = 'UNKNOWN')
      READ (24, *) EFRM
c     format-read characters (for ALPHA compilers):
      READ (24, 111) FRAME
      READ (24, 111) PROJ
      READ (24, 111) TARG
      READ (24, *) IAP
      READ (24, *) IZP
      READ (24, *) IAT
      READ (24, *) IZT
      READ (24, *) NEVNT
      READ (24, *) BMIN
      READ (24, *) BMAX
c     flag to select default AMPT or string melting:
      READ (24, *) isoft
c     read initialization value for hadron cascade:
      READ (24, *) NTMAX
      READ (24, *) DT
c     parj(41) and (42) are a and b parameters in Lund string fragmentation:
      READ (24, *) PARJ(41)
      READ (24, *) PARJ(42)
c     IHPR2(11)=3 (or 2) allows the popcorn mechanism in PYTHIA and 
c     increase the net-baryon stopping in rapidity (value HIJING is 1):
      READ (24, *) ipop
      if(ipop.eq.1) IHPR2(11)=3
c     PARJ(5) controls the fraction of BMBbar vs BBbar in popcorn:
      READ (24, *) PARJ(5)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        注意下面的两个flag
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c     shadowing flag in HIJING:
      READ (24, *) IHPR2(6)
c     quenching flag in HIJING:
      READ (24, *) IHPR2(4)
c     quenching rate when quenching flag is on (=1.0 GeV/fm):
      READ (24, *) HIPR1(14)
c     Minimum pt of hard or semihard scatterings in HIJING: D=2.0 GeV. 
      READ (24, *) HIPR1(8)
c     read initialization value for parton cascade:
      READ (24, *) xmu
      READ (24, *) izpc
      READ (24, *) alpha
c     quark coalescence radii in momentum and space for string melting:
      READ (24, *) dpcoal
      READ (24, *) drcoal
c     flag: read in HIJING random # seed at runtime(1) or from input.ampt(D=0):
      READ (24, *) ihjsed
c     2 seeds for random number generators in HIJING/hadron cascade and ZPC:
      READ (24, *) nseed
      READ (24, *) iseedp
      READ (24, *) iksdcy
      READ (24, *) iphidcy
      READ (24, *) ipi0dcy
c     flag for OSCAR output for final partons and hadrons:
      READ (24, *) ioscar
clin-5/2008     flag for perturbative treatment of deuterons:
      READ (24, *) idpert
      READ (24, *) npertd
      READ (24, *) idxsec
clin-6/2009 To select events that have at least 1 high-Pt minijet parton:
      READ (24, *) pttrig
      READ (24, *) maxmiss
      READ (24, *) IHPR2(2)
      READ (24, *) IHPR2(5)
clin-6/2009 To embed a back-to-back q/qbar pair into each event:
      READ (24, *) iembed
      READ (24, *) pxqembd, pyqembd
      READ (24, *) xembd, yembd
      READ (24, *) nsembd,psembd,tmaxembd
clin-7/2009 Allow modification of nuclear shadowing:
      READ (24, *) ishadow
      READ (24, *) dshadow
      READ (24, *) iphirp
c$$$  xiaohai flag
      READ (24, *) xiaohaiflag
c
      CLOSE (24)
 111  format(a8)
clin-6/2009 ctest off turn on jet triggering:
c      IHPR2(3)=1
c     Trigger Pt of high-pt jets in HIJING:
c      HIPR1(10)=7.
c
      if(isoft.eq.1) then
         amptvn = '1.26t7 (Default)'
      elseif(isoft.eq.4) then
         amptvn = '2.26t7 (StringMelting)'
      else
         amptvn = 'Test-Only'
      endif
      WRITE(6,50) amptvn
      WRITE(12,50) amptvn
 50   FORMAT(' '/
     &11X,'##################################################'/1X,
     &10X,'#      AMPT (A Multi-Phase Transport) model      #'/1X,
     &10X,'#          Version ',a25,                  '     #'/1X,
     &10X,'#               10/28/2016                       #'/1X,
     &10X,'##################################################'/1X,
     &10X,' ')
c     when ihjsed=11: use environment variable at run time for HIJING nseed:
      if(ihjsed.eq.11) then
         PRINT *,
     1 '# Read in NSEED in HIJING at run time (e.g. 20030819):'
      endif
      READ (*, *) nseedr
      if(ihjsed.eq.11) then
         nseed=nseedr
      endif
      if(ihjsed.eq.11) then      
         PRINT *, '#   read in: ', nseed
         WRITE(12,*) '# Read in NSEED in HIJING at run time:',nseed
      endif
      CLOSE(12)
clin-5/2015 an odd number is needed for the random number generator:
c      if(mod(NSEED,2).eq.0) NSEED=NSEED+1
      NSEED=2*NSEED+1
c     9/26/03 random number generator for f77 compiler:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        随机数生成器的初始化函数
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      CALL SRAND(NSEED)
c
c.....turn on warning messages in nohup.out when an event is repeated:
      IHPR2(10) = 1
c     string formation time:
      ARPAR1(1) = 0.7
c     smearp is the smearing halfwidth on parton z0, 
c     set to 0 for now to avoid overflow in eta.
c     smearh is the smearing halfwidth on string production point z0.
      smearp=0d0
      IAmax=max(iap,iat)
      smearh=1.2d0*IAmax**0.3333d0/(dble(EFRM)/2/0.938d0)
      nevent=NEVNT
c
c     AMPT momentum and space info at freezeout:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        冻结后的相空间的信息
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     文件号的递减编号
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      OPEN (16, FILE = 'ana/ampt.dat', STATUS = 'UNKNOWN')
      OPEN (14, FILE = 'ana/zpc.dat', STATUS = 'UNKNOWN')
      open(10000, file="testOut/subroutinexiaohai.dat")
      open(9999, file = 'testOut/hijsetxiaohai.dat', status = "unknown")
      open(9998, file = 'testOut/hijwdsxiaohai.dat', status = 'unknown')
c$$$  9997 和 9996给了parton.oscar和hadron.oscar
      open(9995, file = 'testOut/ampt.dat', status = 'unknown')
      open(9994, file = "testOut/amptHijing.dat", status = 'unknown')
      open(9993, file = "testOut/hijinixiaohai.dat", status = "unknown")
      open(9992, file = "testOut/hijcscxiaohai.dat", status = 'unknown')
      open(9991, file = "testOut/jetinixiaohai.dat",status='unknown')
      open(9990, file = "testOut/jetinfoxiaohai.dat", status='unknown')
      open(9989, file = "testOut/jetnumxiaohai.dat", status='unknown')
      open(9988, file = "testOut/shijijetxiaohai.dat", status='unknown')
      open(9987,file="testOut/projectilexiaohai.dat",status="unknown")
      open(9986,file="testOut/targetxiaohai.dat",status="unknown")
      open(9985,file="testOut/phixiaohai.dat", status="unknown")
      open(9984,file="testOut/phi2xiaohai.dat",status="unknown")
      open(9983,file="testOut/radiusxiaohai.dat",status="unknown")
      open(9982,file="testOut/inirecxiaohai.dat",status="unknown")
      open(9981,file="testOut/zpcmnxiaohai.dat",status="unknown")
      open(9980,file="testOut/callzpcxiaohai.dat")
      open(9979,file="testOut/nsbrunzpcmnxiaohai.dat")
      open(9978,file="testOut/neventzpcmnxiaohai.dat")
      open(9977,file="testOut/htopxiaohai.dat")
      open(9976,file="testOut/coalesxiaohai.dat")
      open(9975,file="testOut/coordinatexiaohai.dat")
      open(9974,file="testOut/partoncoordinatexiaohai.dat")
      open(9973,file="testOut/hadroncoordinatexiaohai.dat")
      open(9972,file="testOut/nnozpcxiaohai.dat")
      open(9971,file="testOut/getnpxiaohai.dat")
      open(9970,file="testOut/arinixiaohai.dat")
      open(9969,file="testOut/arini1xiaohai.dat")
      open(9968,file="testOut/addhadxiaohai.dat")
      open(9967,file="testOut/artordxiaohai.dat")
      open(9966,file="testOut/inizpcxiaohai.dat")
      open(9965,file="testOut/iniparxiaohai.dat")
      open(9964,file="testOut/inian1xiaohai.dat")
      open(9963,file="testOut/readpaxiaohai.dat")
      open(9962,file="testOut/nclotxiaohai.dat")
      open(9961,file="testOut/nhardxiaohai.dat")
      open(9960,file="testOut/jethappenxiaohai.dat")
      open(9959,file="testOut/initypexiaohai.dat")
      open(9958,file="testOut/typejpjtxiaohai.dat")
      open(9957,file="testOut/hijhrdxiaohai.dat")
      open(9956,file="testOut/firstscatterxiaohai.dat")
      open(9955,file="testOut/quarkantiquarkxiaohai.dat")
      open(9954,file="testOut/minijetoutxiaohai.dat")
      open(9953,file="testOut/minijetxiaohai.dat")
      open(9952,file="testOut/hboostxiaohai.dat")
      open(9951,file="testOut/nptjxiaohai.dat")
      open(9950,file="testOut/quenchxiaohai.dat")
      open(9949,file="testOut/nsgxiaohai.dat")
      open(9948,file="testOut/Nxiaohai.dat")
      open(9947,file="testOut/k2flavorxiaohai.dat")
      open(9946,file="testOut/nattxiaohai.dat")
      open(9945,file="testOut/natt2xiaohai.dat")
      open(9944,file="testOut/natt3xiaohai.dat")
      open(9943,file="testOut/embedhighptxiaohai.dat")
      open(9942,file="testOut/natt4xiaohai.dat")
      open(9941,file="testOut/hjana1xiaohai.dat")
      open(9940,file="testOut/pargluxiaohai.dat")
      open(9939,file="testOut/mulxiaohai.dat")
      open(9938,file="testOut/countnumxiaohai.dat")
      open(9937,file="testOut/natt5xiaohai.dat")
      open(9936,file="testOut/nattcountxiaohai.dat")
      open(9935,file="testOut/jinruzpcxiaohai.dat")
      open(9934,file="testOut/bujinruzpcxiaohai.dat")
      open(9933,file="testOut/sancanshuxiaohai.dat")
      open(9932,file="testOut/readixiaohai.dat")
      open(9931,file="testOut/inievtxiaohai.dat")
      open(9930,file="testOut/inirunxiaohai.dat")
      open(9929,file="testOut/ftimexiaohai.dat")
      open(9928,file="testOut/indxxiaohai.dat")
      open(9927,file="testOut/inixiaohai.dat")
      open(9926,file="testOut/equaldistancexiaohai.dat")
      open(9925,file="testOut/iilistxiaohai.dat")
      open(9924,file="testOut/zpcrunxiaohai.dat")
      open(9923,file="testOut/getictxiaohai.dat")
      open(9922,file="testOut/collisionxiaohai.dat")
      open(9921,file="testOut/newposxiaohai.dat")
      open(9920,file="testOut/iaevtpxiaohai.dat")
      open(9919,file="testOut/t2timexiaohai.dat")
      open(9918,file="testOut/xncollxiaohai.dat")
      open(9917,file="testOut/iloopxiaohai.dat")
      open(9916,file="testOut/v2partonxiaohai.dat")
      open(9915,file="testOut/filenum49xiaohai.dat")
      open(9914,file="testOut/ntmaxdtxiaohai.dat")
      open(9913,file="testOut/p-animatexiaohai.dat")
      open(9912,file="testOut/parton-txiaohai.dat")
      open(9911,file="testOut/zpca2xiaohai.dat")
      open(9910,file="testOut/zpcouxiaohai.dat")
      open(9909,file="testOut/zpcou1xiaohai.dat")
      open(9908,file="testOut/zpstrgxiaohai.dat")
      open(9907,file="testOut/istrgxiaohai.dat")
      open(9906,file="testOut/nstrxiaohai.dat")
      open(9905,file="testOut/afterzpcxiaohai.dat")
      open(9904,file="testOut/hjana2xiaohai.dat")
      open(9903,file="testOut/hjana2partonxiaohai.dat")
      open(9902,file="testOut/nisgxiaohai.dat")
      open(9901,file="testOut/ihpr2_20xiaohai.dat")
      open(9900,file="testOut/ihpr2_20newxiaohai.dat")
      open(9899,file="testOut/sorthadronxiaohai.dat")
      open(9898,file="testOut/arini2xiaohai.dat")
      open(9897,file="testOut/centralityxiaohai.dat")
      open(9896,file="testOut/labkinematicsxiaohai.dat")
      open(9895,file="testOut/artimpactparxiaohai.dat")
      open(9894,file="testOut/artmassxiaohai.dat")
      open(9893,file="testOut/traceparton1xiaohai.dat")
      open(9892,file="testOut/traceparton2xiaohai.dat")
      open(9891,file="testOut/traceparton3xiaohai.dat")
      open(9890,file="testOut/timevaluexiaohai.dat")
      open(9889,file="testOut/zpcrunreturnxiaohai.dat")
      open(9888,file="testOut/updatalistxiaohai.dat")
      open(9887,file="testOut/tracetimexiaohai.dat")
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  行号的递减编号
c$$$  521
c$$$  520  
c$$$  519
c$$$  518
c$$$  517
c$$$  516
c$$$  515
c$$$  514
c$$$  513
c$$$  512
c$$$  511
c$$$  510
c$$$  509
c$$$  508
c$$$  507
c$$$  506
c$$$  505
c$$$  504
c$$$  503
c$$$  502
c$$$  501
c$$$  500
c$$$  499
c$$$  498
c$$$  497
c$$$  496
c$$$  495
c$$$  494
c$$$  493
c$$$  492
c$$$  491  
c$$$  490
c$$$  489
c$$$  488
c$$$  487
c$$$  486
c$$$  485
c$$$  484
c$$$  483
c$$$  482
c$$$  481
c$$$  480
c$$$  479
c$$$  478
c$$$  477
c$$$  476
c$$$  475
c$$$  474
c$$$  473
c$$$  472
c$$$  471
c$$$  470
c$$$  469
c$$$  468
c$$$  467
c$$$  466
c$$$  465
c$$$  464
c$$$  463
c$$$  462
c$$$  461
c$$$  460
c$$$  459
c$$$  458
c$$$  457
c$$$  456
c$$$ 455  
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
ctest off for resonance (phi, K*) studies:
c      OPEN (17, FILE = 'ana/res-gain.dat', STATUS = 'UNKNOWN')
c     OPEN (18, FILE = 'ana/res-loss.dat', STATUS = 'UNKNOWN')
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        能量，坐标系，弹核和靶核碰撞的类型，剩下的四个参数为核子数
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      CALL HIJSET(EFRM, FRAME, PROJ, TARG, IAP, IZP, IAT, IZT)
      CALL ARTSET
      CALL INIZPC
clin-5/2009 ctest off:
c      call flowp(0)
c      call flowh0(NEVNT,0)
c      call iniflw(NEVNT,0)
c      call frztm(NEVNT,0)
c
       DO 2000 J = 1, NEVNT
          IAEVT = J
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$         WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$          WW     WW WW     WW  RR  RR   II     TT     EE
c$$$           WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$            WW WW     WW WW    RR RR    II     TT     EE
c$$$             WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
          write(9995, *) "NUM is ", NUM
c$$$          NUM的数值永远为1.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$         WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$          WW     WW WW     WW  RR  RR   II     TT     EE
c$$$           WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$            WW WW     WW WW    RR RR    II     TT     EE
c$$$             WW        WW      RR  RR   II     TT     EEEEEE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
          DO 1000 K = 1, NUM
             IARUN = K
             IF (IAEVT .EQ. NEVNT .AND. IARUN .EQ. NUM) THEN
                IOUT = 1
             END IF
             PRINT *, ' EVENT ', J, ', RUN ', K
             imiss=0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        传入坐标系，碰撞参数的上下限。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 100         CALL HIJING(FRAME, BMIN, BMAX)
             IAINT2(1) = NATT             
clin-6/2009 ctest off
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$         WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$          WW     WW WW     WW  RR  RR   II     TT     EE
c$$$           WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$            WW WW     WW WW    RR RR    II     TT     EE
c$$$             WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
                write(9994,*) HIPR1
                write(9994,*) ' '
                write(9994,*) IHPR2
                write(9994,*) ' '
                write(9994,*) (HINT1(i),i=1,20)
                write(9994,*) ' '
                write(9994,*) (HINT1(i),i=21,40)
                write(9994,*) ' '
                write(9994,*) (HINT1(i),i=41,60)
                write(9994,*) ' '
                write(9994,*) (HINT1(i),i=61,80)
                write(9994,*) ' '
                write(9994,*) (HINT1(i),i=81,100)
                write(9994,*) ' '
                write(9994,*) IHNT2
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$         WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$          WW     WW WW     WW  RR  RR   II     TT     EE
c$$$           WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$            WW WW     WW WW    RR RR    II     TT     EE
c$$$             WW        WW      RR  RR   II     TT     EEEEEE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           if(J.eq.-2) then 
              write(98,*) HIPR1
              write(98,*) ' '
              write(98,*) IHPR2
              write(98,*) ' '
              write(98,*) (HINT1(i),i=1,20)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=21,40)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=41,60)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=61,80)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=81,100)
              write(98,*) ' '
              write(98,*) IHNT2
           endif
c     evaluate Npart (from primary NN collisions) for both proj and targ:
             call getnp
c     switch for final parton fragmentation:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$         IHPR2(20) = 1
         write(9900,*)"IHPR2(20) ==>  ", IHPR2(20)
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
             IF (IHPR2(20) .EQ. 0) GOTO 2000
c     In the unlikely case of no interaction (even after loop of 20 in HIJING),
c     still repeat the event to get an interaction 
c     (this may have an additional "trigger" effect):
             if(NATT.eq.0) then
                imiss=imiss+1
                if(imiss.le.20) then
                   write(6,*) 'repeated event: natt=0,j,imiss=',j,imiss
                   goto 100
                else
                   write(6,*) 'missed event: natt=0,j=',j
                   goto 2000
                endif
             endif
c.....ART initialization and run
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  arini:初始化粒子的信息。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
             CALL ARINI
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  arini2():初始粒子的末态的时间到ntmax*dt.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
             CALL ARINI2(K)
 1000     CONTINUE
c

c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  初始化那一堆数组.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
          CALL ARTAN1
clin-9/2012 Analysis is not used:
c          CALL HJANA3


          CALL ARTMN
clin-9/2012 Analysis is not used:
c          CALL HJANA4
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
          CALL ARTAN2
 2000  CONTINUE
c
       CALL ARTOUT(NEVNT)
clin-5/2009 ctest off:
c       call flowh0(NEVNT,2)
c       call flowp(2)
c       call iniflw(NEVNT,2)
c       call frztm(NEVNT,2)
c
       STOP
       END
c     FYI: taken file unit numbers are 10-88, 91-93; 
c     so free file unit numbers are 1-4,7-9,89,97-99.
c....................amptsub.f
c.....this file contains 4 sections:
c.....1. ART subroutines;
c.....2. ART functions;
c.....3. ART block data;
c.....4. subprocesses borrowed from other codes.
c.....5. the previous artana.f
c.....6. the previous zpcsub.f
c.....7. subroutine getnp
c.....Note that Parts1-4 are the previous artsub.f
c
c=======================================================================
c.....subroutine to set up ART parameters and analysis files
c.....before looping different events
