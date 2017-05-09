        subroutine minijet_out(BB,phiRP)
        PARAMETER (MAXSTR=150001)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
        COMMON/hjcrdn/YP(3,300),YT(3,300)
        COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &                PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &                PJPM(300,500),NTJ(300),KFTJ(300,500),
     &                PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &                PJTE(300,500),PJTM(300,500)
        COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &       K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &       PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
        COMMON /AREVT/ IAEVT, IARUN, MISS
        common /para7/ ioscar,nsmbbbar,nsmmeson
        common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
        SAVE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9954,*)"pttrig = ", pttrig
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  其中NPJ(I)是第I个核子产生的partons.而下面的程序是想统计产生的
c$$$        minijet的数量。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ntrig=0
        do I = 1, IHNT2(1)
           do J = 1, NPJ(I)
              pt=sqrt(PJPX(I,J)**2+PJPY(I,J)**2)
              if(pt.ge.pttrig) ntrig=ntrig+1
           enddo
        enddo
        do I = 1, IHNT2(3)
           do J = 1, NTJ(I)
              pt=sqrt(PJTX(I,J)**2+PJTY(I,J)**2)
              if(pt.ge.pttrig) ntrig=ntrig+1
           enddo
        enddo
        do I = 1, NSG
           do J = 1, NJSG(I)
              pt=sqrt(PXSG(I,J)**2+PYSG(I,J)**2)
              if(pt.ge.pttrig) ntrig=ntrig+1
           enddo
        enddo
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    从最后的标号来区分parton的来源。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9953,*) IAEVT,MISS,IHNT2(1),IHNT2(3)
        DO 516 I = 1, IHNT2(1)
           DO 515 J = 1, NPJ(I)
              xiaohaiityp = KFPJ(I,J)
              xiaohaigx = YP(1,I)+0.5*BB*cos(phiRP)
              xiaohaigy = YP(2,I)+0.5*BB*sin(phiRP)
              xiaohaigz = 0
              xiaohaift = 0
              xiaohaipx = PJPX(I,J)
              xiaohaipy = PJPY(I,J)
              xiaohaipz = PJPZ(I,J)
              xiaohaixmass = PJPM(I,J)
              write(9953,*)xiaohaiityp, xiaohaipx, xiaohaipy,
     &             xiaohaipz, xiaohaixmass, xiaohaigx,
     &             xiaohaigy, xiaohaigz, xiaohaift, 1
 515       continue
           if (NPJ(I).NE.0) write(9953,*)
 516    continue
        write(9953,*)"new line."
        DO 514 I = 1, IHNT2(3)
           DO 513 J = 1, NTJ(I)
              xiaohaiityp = KFTJ(I,J)
              xiaohaigx = YT(1,I)+0.5*BB*cos(phiRP)
              xiaohaigy = YT(2,I)+0.5*BB*sin(phiRP)
              xiaohaigz = 0
              xiaohaift = 0
              xiaohaipx = PJTX(I,J)
              xiaohaipy = PJTY(I,J)
              xiaohaipz = PJTZ(I,J)
              xiaohaixmass = PJTM(I,J)
              write(9953,*)xiaohaiityp, xiaohaipx, xiaohaipy,
     &             xiaohaipz, xiaohaixmass, xiaohaigx,
     &             xiaohaigy, xiaohaigz, xiaohaift, 2
 513       continue 
           if (NTJ(I).NE.0) write(9953,*)
 514    continue
        write(9953,*)"new line."
      DO 512 I = 1, NSG
         DO 511 J = 1, NJSG(I)
            xiaohaiityp = K2SG(I,J)
            xiaohaigx = 0.5*(YP(1,IASG(I,1))+YT(1,IASG(I,2)))
            xiaohaigy = 0.5*(YP(2,IASG(I,1))+YT(2,IASG(I,2)))
            xiaohaigz = 0
            xiaohaift = 0
            xiaohaipx = PXSG(I,J)
            xiaohaipy = PYSG(I,J)
            xiaohaipz = PZSG(I,J)
            xiaohaixmass = PMSG(I,J)
            write(9953,*)xiaohaiityp, xiaohaipx, xiaohaipy,
     &           xiaohaipz, xiaohaixmass, xiaohaigx,
     &           xiaohaigy, xiaohaigz, xiaohaift, 3
 511     continue
           if (NJSG(I).NE.0) write(9953,*)
 512  continue
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        
c     Require at least 1 initial minijet parton above the trigger Pt value:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    没有minijet产生。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(ntrig.eq.0) return
c.....transfer data from HIJING to ZPC
        if(ioscar.eq.3) write(96,*) IAEVT,MISS,IHNT2(1),IHNT2(3)
        DO 1008 I = 1, IHNT2(1)
           DO 1007 J = 1, NPJ(I)
              ityp=KFPJ(I,J)
c     write out not only gluons:
c              if(ityp.ne.21) goto 1007
clin-2/2012:
c              gx=YP(1,I)+0.5*BB
c              gy=YP(2,I)
              gx=YP(1,I)+0.5*BB*cos(phiRP)
              gy=YP(2,I)+0.5*BB*sin(phiRP)
              gz=0.
              ft=0.
              px=PJPX(I,J)
              py=PJPY(I,J)
              pz=PJPZ(I,J)
              xmass=PJPM(I,J)
              if(ioscar.eq.3) then
                 if(amax1(abs(gx),abs(gy),
     1                abs(gz),abs(ft)).lt.9999) then
                    write(96,200) ityp,px,py,pz,xmass,gx,gy,gz,ft,1
                 else
                    write(96,201) ityp,px,py,pz,xmass,gx,gy,gz,ft,1
                 endif
              endif
 1007      CONTINUE
 1008   CONTINUE
        DO 1010 I = 1, IHNT2(3)
           DO 1009 J = 1, NTJ(I)
              ityp=KFTJ(I,J)
c              if(ityp.ne.21) goto 1009
clin-2/2012:
c              gx=YT(1,I)-0.5*BB
c              gy=YT(2,I)
              gx=YT(1,I)-0.5*BB*cos(phiRP)
              gy=YT(2,I)-0.5*BB*sin(phiRP)
              gz=0.
              ft=0.
              px=PJTX(I,J)
              py=PJTY(I,J)
              pz=PJTZ(I,J)
              xmass=PJTM(I,J)
              if(ioscar.eq.3) then
                 if(amax1(abs(gx),abs(gy),
     1                abs(gz),abs(ft)).lt.9999) then
                    write(96,200) ityp,px,py,pz,xmass,gx,gy,gz,ft,2
                 else
                    write(96,201) ityp,px,py,pz,xmass,gx,gy,gz,ft,2
                 endif
              endif
 1009      CONTINUE
 1010   CONTINUE
        DO 1012 I = 1, NSG
           DO 1011 J = 1, NJSG(I)
              ityp=K2SG(I,J)
c              if(ityp.ne.21) goto 1011
              gx=0.5*(YP(1,IASG(I,1))+YT(1,IASG(I,2)))
              gy=0.5*(YP(2,IASG(I,1))+YT(2,IASG(I,2)))
              gz=0.
              ft=0.
              px=PXSG(I,J)
              py=PYSG(I,J)
              pz=PZSG(I,J)
              xmass=PMSG(I,J)
              if(ioscar.eq.3) then
                 if(amax1(abs(gx),abs(gy),
     1                abs(gz),abs(ft)).lt.9999) then
                    write(96,200) ityp,px,py,pz,xmass,gx,gy,gz,ft,3
                 else
                    write(96,201) ityp,px,py,pz,xmass,gx,gy,gz,ft,3
                 endif
              endif
 1011      CONTINUE
 1012   CONTINUE
 200  format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,2(1x,f8.2),2(2x,f2.0),2x,I2)
 201  format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,2(1x,e8.2),2(2x,f2.0),2x,I2)
 510  format(1x,f8.0, 7(1x,f8.3), 1x,f8.0)
c
        return
        end
c=======================================================================
clin-6/2009 embed back-to-back high-Pt quark/antiquark pair
c     via embedding back-to-back high-Pt pion pair then melting the pion pair
c     by generating the internal quark and antiquark momentum parallel to 
c      the pion momentum (in order to produce a high-Pt and a low Pt parton):
