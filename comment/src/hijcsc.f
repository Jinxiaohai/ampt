c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  接受的两个参数是弹核和靶核的单个核子
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        SUBROUTINE HIJCSC(JP,JT)
        DIMENSION PSC1(5),PSC2(5)
        COMMON/hjcrdn/YP(3,300),YT(3,300)
cc      SAVE /hjcrdn/
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
        COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
cc      SAVE /HSTRNG/
        common /xiaohai/xiaohaiflag
        SAVE   
        IF(JP.EQ.0 .OR. JT.EQ.0) GO TO 25
        DO 10 I=1,5
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  两个核子的四动量和不变质量
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSC1(I)=PP(JP,I)
        PSC2(I)=PT(JT,I)
10        CONTINUE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  输出弹散之前的两个核子和弹散之后的两个核子的四动量。处理的
c$$$          四动量是在质心系下进行处理的
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
          if(xiaohaiflag .eq. 6) then
             write(9992, *), PSC1(1),"  ", PSC1(2), "  ", PSC1(3), "  ",
     &       PSC1(4)
             write(9992, *), PSC2(1),"  ", PSC2(2), "  ", PSC2(3), "  ",
     &        PSC2(4)
          endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        CALL HIJELS(PSC1,PSC2)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
          if(xiaohaiflag .eq. 6) then
             write(9992, *), PSC1(1),"  ", PSC1(2), "  ", PSC1(3), "  ",
     &       PSC1(4)
             write(9992, *), PSC2(1),"  ", PSC2(2), "  ", PSC2(3), "  ",
     &       PSC2(4)
          endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     这次碰撞在横向方向交换的动量
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DPP1=PSC1(1)-PP(JP,1)
        DPP2=PSC1(2)-PP(JP,2)
        DPT1=PSC2(1)-PT(JT,1)
        DPT2=PSC2(2)-PT(JT,2)
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
        PP(JP,6)=PP(JP,6)+DPP1/2.0
        PP(JP,7)=PP(JP,7)+DPP2/2.0
        PP(JP,8)=PP(JP,8)+DPP1/2.0
        PP(JP,9)=PP(JP,9)+DPP2/2.0
        PT(JT,6)=PT(JT,6)+DPT1/2.0
        PT(JT,7)=PT(JT,7)+DPT2/2.0
        PT(JT,8)=PT(JT,8)+DPT1/2.0
        PT(JT,9)=PT(JT,9)+DPT2/2.0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     末态四动量
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 20 I=1,4
        PP(JP,I)=PSC1(I)
        PT(JT,I)=PSC2(I)
20        CONTINUE
        NFP(JP,5)=MAX(1,NFP(JP,5))
        NFT(JT,5)=MAX(1,NFT(JT,5))
C                ********Perform elastic scattering between JP and JT
        RETURN
C                ********The following is for possible elastic cascade
c
25        IF(JP.EQ.0) GO TO 45
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     多次弹散的发生的处理
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PABS=SQRT(PP(JP,1)**2+PP(JP,2)**2+PP(JP,3)**2)
        BX=PP(JP,1)/PABS
        BY=PP(JP,2)/PABS
        BZ=PP(JP,3)/PABS
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    对于弹核中的每个核子遍历弹核的每个核子
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 40 I=1,IHNT2(1)
                IF(I.EQ.JP) GO TO 40
                DX=YP(1,I)-YP(1,JP)
                DY=YP(2,I)-YP(2,JP)
                DZ=YP(3,I)-YP(3,JP)
                DIS=DX*BX+DY*BY+DZ*BZ
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
                if(xiaohaiflag .eq. 6) then
                   write(9992, *)"elastic cascade happened : ", DIS
                endif
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
                IF(DIS.LE.0) GO TO 40
                BB=DX**2+DY**2+DZ**2-DIS**2
                R2=BB*HIPR1(40)/HIPR1(31)/0.1
C                ********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
                GS=1.0-EXP(-(HIPR1(30)+HINT1(11))/HIPR1(31)/2.0
     &                        *ROMG(R2))**2
                GS0=1.0-EXP(-(HIPR1(30)+HINT1(11))/HIPR1(31)/2.0
     &                        *ROMG(0.0))**2
                IF(RANART(NSEED).GT.GS/GS0) GO TO 40
                DO 30 K=1,5
                        PSC1(K)=PP(JP,K)
                        PSC2(K)=PP(I,K)
30                CONTINUE
                CALL HIJELS(PSC1,PSC2)
                DPP1=PSC1(1)-PP(JP,1)
                DPP2=PSC1(2)-PP(JP,2)
                DPT1=PSC2(1)-PP(I,1)
                DPT2=PSC2(2)-PP(I,2)
                PP(JP,6)=PP(JP,6)+DPP1/2.0
                PP(JP,7)=PP(JP,7)+DPP2/2.0
                PP(JP,8)=PP(JP,8)+DPP1/2.0
                PP(JP,9)=PP(JP,9)+DPP2/2.0
                PP(I,6)=PP(I,6)+DPT1/2.0
                PP(I,7)=PP(I,7)+DPT2/2.0
                PP(I,8)=PP(I,8)+DPT1/2.0
                PP(I,9)=PP(I,9)+DPT2/2.0
                DO 35 K=1,5
                        PP(JP,K)=PSC1(K)
                        PP(I,K)=PSC2(K)
35                CONTINUE
                NFP(I,5)=MAX(1,NFP(I,5))
                GO TO 45
40        CONTINUE
45        IF(JT.EQ.0) GO TO 80
clin 50        PABS=SQRT(PT(JT,1)**2+PT(JT,2)**2+PT(JT,3)**2)
        PABS=SQRT(PT(JT,1)**2+PT(JT,2)**2+PT(JT,3)**2)
        BX=PT(JT,1)/PABS
        BY=PT(JT,2)/PABS
        BZ=PT(JT,3)/PABS
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     对于弹核中的每个核子遍历靶核中的每个核子。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 70 I=1,IHNT2(3)
                IF(I.EQ.JT) GO TO 70
                DX=YT(1,I)-YT(1,JT)
                DY=YT(2,I)-YT(2,JT)
                DZ=YT(3,I)-YT(3,JT)
                DIS=DX*BX+DY*BY+DZ*BZ
                IF(DIS.LE.0) GO TO 70
                BB=DX**2+DY**2+DZ**2-DIS**2
                R2=BB*HIPR1(40)/HIPR1(31)/0.1
C                ********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
                GS=(1.0-EXP(-(HIPR1(30)+HINT1(11))/HIPR1(31)/2.0
     &                        *ROMG(R2)))**2
                GS0=(1.0-EXP(-(HIPR1(30)+HINT1(11))/HIPR1(31)/2.0
     &                        *ROMG(0.0)))**2
                IF(RANART(NSEED).GT.GS/GS0) GO TO 70
                DO 60 K=1,5
                        PSC1(K)=PT(JT,K)
                        PSC2(K)=PT(I,K)
60                CONTINUE
                CALL HIJELS(PSC1,PSC2)
                DPP1=PSC1(1)-PT(JT,1)
                DPP2=PSC1(2)-PT(JT,2)
                DPT1=PSC2(1)-PT(I,1)
                DPT2=PSC2(2)-PT(I,2)
                PT(JT,6)=PT(JT,6)+DPP1/2.0
                PT(JT,7)=PT(JT,7)+DPP2/2.0
                PT(JT,8)=PT(JT,8)+DPP1/2.0
                PT(JT,9)=PT(JT,9)+DPP2/2.0
                PT(I,6)=PT(I,6)+DPT1/2.0
                PT(I,7)=PT(I,7)+DPT2/2.0
                PT(I,8)=PT(I,8)+DPT1/2.0
                PT(I,9)=PT(I,9)+DPT2/2.0
                DO 65 K=1,5
                        PT(JT,K)=PSC1(K)
                        PT(I,K)=PSC2(K)
65                CONTINUE
                NFT(I,5)=MAX(1,NFT(I,5))
                GO TO 80
70        CONTINUE
80        RETURN
        END
C
C
C*******************************************************************
CThis subroutine performs elastic scattering between two nucleons
C
C*******************************************************************
