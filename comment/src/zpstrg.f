      SUBROUTINE zpstrg
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXPTN=400001)
      PARAMETER (MAXSTR=150001)
c      REAL*4 YP, YT, PXSG, PYSG, PZSG, PESG, PMSG, HIPR1, HINT1, BB
      REAL YP, YT, PXSG, PYSG, PZSG, PESG, PMSG, HIPR1, HINT1, BB
      COMMON /PARA1/ MUL
cc      SAVE /PARA1/
      COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &   PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &   XMASS5(MAXPTN), ITYP5(MAXPTN)
cc      SAVE /prec2/
      COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)
cc      SAVE /ilist8/
      COMMON /SREC1/ NSP, NST, NSI
cc      SAVE /SREC1/
      COMMON /SREC2/ATAUI(MAXSTR),ZT1(MAXSTR),ZT2(MAXSTR),ZT3(MAXSTR)
cc      SAVE /SREC2/
      COMMON/hjcrdn/YP(3,300),YT(3,300)
cc      SAVE /hjcrdn/
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &   K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &   PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
cc      SAVE /HJJET2/
cbz6/28/99 flow1
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
cbz6/28/99 flow1 end
      common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
      common/strg/np(maxstr)
cc      SAVE /strg/
clin-6/06/02 test local freezeout:
      common /frzprc/ 
     &     gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &     pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &     xmfrz(MAXPTN), 
     &     tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
cc      SAVE /frzprc/
      SAVE   
clin-6/06/02 test local freezeout for string melting,
c     use space-time values at local freezeout saved in /frzprc/:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      if(isoft .eq. 4) then
         do 472 ihai = 1, MUL
            write(9908,*)idfrz(ihai), pxfrz(ihia), pyfrz(ihai),
     &           pzfrz(ihai), xmfrz(ihai), gxfrz(ihai), gyfrz(ihai),
     &           gzfrz(ihai), ftfrz(ihai)
 472     continue      
      endif
      
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      if(isoft.eq.5) then
         do 1001 I = 1, MUL
            ITYP5(i)=idfrz(i)
            GX5(i)=gxfrz(i)
            GY5(i)=gyfrz(i)
            GZ5(i)=gzfrz(i)
            FT5(i)=ftfrz(i)
            PX5(i)=pxfrz(i)
            PY5(i)=pyfrz(i)
            PZ5(i)=pzfrz(i)
            E5(i)=efrz(i)
            XMASS5(i)=xmfrz(i)
 1001    continue
      endif
clin-6/06/02-end
      DO 1002 I = 1, MAXSTR
         ATAUI(I) = 0d0
         ZT1(I) = 0d0
         ZT2(I) = 0d0
clin-4/25/03 add zt3(I) to track longitudinal positions of partons/strings:
         ZT3(I) = 0d0
         NP(I) = 0
 1002 CONTINUE
      DO 1003 I = 1, MUL
         ISTRG = LSTRG1(I)
clin-9/2015 to avoid Floating-Point Exception:
c         TAU7 = SQRT(FT5(I) ** 2 - GZ5(I) ** 2)
         diff2=FT5(I)**2 - GZ5(I)**2
         if(diff2.lt.0d0) then
            write(6,*) '2:I,ft5,gz5,diff2=',I,ft5(i),gz5(i),diff2
            TAU7=1d-6
         else
            TAU7 = dSQRT(diff2)
         endif
c
         ATAUI(ISTRG) = ATAUI(ISTRG) + TAU7
         ZT1(ISTRG) = ZT1(ISTRG) + GX5(I)
         ZT2(ISTRG) = ZT2(ISTRG) + GY5(I)
         ZT3(ISTRG) = ZT3(ISTRG) + GZ5(I)
         NP(ISTRG) = NP(ISTRG) + 1
 1003 CONTINUE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      do 471 I = 1, MUL
         ISTRG = LSTRG1(I)
         write(9907,*)"ISTRG = ", istrg
 471     continue
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<

         NSTR = NSP + NST + NSI
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
         write(9906,*)"NSTR    NSP    NST    NSI"
         write(9906,*)nstr, nsp, nst, nsi
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
clin-7/03/01 correct averaging on transverse coordinates, no shift needed:
      if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
         DO 1004 I = 1, NSTR
            IF (NP(I) .NE. 0) THEN
               ATAUI(I) = ATAUI(I) / NP(I)
               ZT1(I) = ZT1(I) / NP(I)
               ZT2(I) = ZT2(I) / NP(I)
               ZT3(I) = ZT3(I) / NP(I)
            ENDIF
 1004    CONTINUE
         return
      endif
clin-7/03/01-end
      DO 1005 I = 1, NSTR
         IF (NP(I) .NE. 0) THEN
            ATAUI(I) = ATAUI(I) / NP(I)
            ZT1(I) = ZT1(I) / NP(I)
            ZT2(I) = ZT2(I) / NP(I)
            ZT3(I) = ZT3(I) / NP(I)
         ELSE
            IF (I .LE. NSP) THEN
               J = I
               ZT1(I) = dble(YP(1, J))
               ZT2(I) = dble(YP(2, J))
               ZT3(I) = 0d0
            ELSE IF (I .GT. NSP .AND. I .LE. NSP + NST) THEN
               J = I - NSP
               ZT1(I) = dble(YT(1, J))
               ZT2(I) = dble(YT(2, J))
               ZT3(I) = 0d0
            ELSE
               J = I - NSP - NST
               ZT1(I) = 0.5d0*
     1              dble((YP(1, IASG(J, 1)) + YT(1, IASG(J, 2))))
               ZT2(I) = 0.5d0 *
     1              dble((YP(2, IASG(J, 1)) + YT(2, IASG(J, 2))))
               ZT3(I) = 0d0
            END IF
         END IF
 1005 CONTINUE
cbz6/28/99 flow1
      BB = HINT1(19)
      DO 1006 I = 1, NSTR
         IF (NP(I).NE.0) THEN
            SHIFT=0d0
         ELSE
            SHIFT=0.5d0*dble(BB)
         END IF
         IF (I .LE. NSP) THEN
            ZT1(I) = ZT1(I) + SHIFT
         ELSE IF (I .GT. NSP .AND. I .LE. NSP + NST) THEN
            ZT1(I) = ZT1(I) - SHIFT
         END IF
 1006 CONTINUE
cbz6/28/99 flow1 end
c
      RETURN
      END
clin-10/01/03 random number generator for f77:
