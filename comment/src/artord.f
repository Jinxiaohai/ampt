      SUBROUTINE ARTORD
c.....before invoking ARTORD:
c.....IAINT2(1) must be set:
      PARAMETER (MAXSTR=150001,MAXR=1)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
cc      SAVE /ARPRNT/
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
cc      SAVE /ARPRC/
clin-3/2009 Take care of particle weights when user inserts initial hadrons:
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      DIMENSION dptemp(MAXSTR)
c
      DIMENSION ITYP0(MAXSTR), 
     &   GX0(MAXSTR), GY0(MAXSTR), GZ0(MAXSTR), FT0(MAXSTR),
     &   PX0(MAXSTR), PY0(MAXSTR), PZ0(MAXSTR), EE0(MAXSTR),
     &   XM0(MAXSTR)
      DIMENSION INDX(MAXSTR)
      EXTERNAL ARINDX
      SAVE   
c
      NPAR = 0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  输出初始的信息.
      do 518 I = 1, IAINT2(1)
         write(9967,*)ITYPAR(I),"    ",PXAR(I),"    ", PYAR(I),"    ",
     &        PZAR(I),"    ",XMAR(I),"    ",GXAR(I),"    ",
     &        GYAR(I),"    ",GZAR(I),"    ",FTAR(I),"    ",PEAR(I)
 518     CONTINUE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      NP = IAINT2(1)
      DO 1001 I = 1, NP
         ITYP0(I) = ITYPAR(I)
         GX0(I) = GXAR(I)
         GY0(I) = GYAR(I)
         GZ0(I) = GZAR(I)
         FT0(I) = FTAR(I)
         PX0(I) = PXAR(I)
         PY0(I) = PYAR(I)
         PZ0(I) = PZAR(I)
         EE0(I) = PEAR(I)
         XM0(I) = XMAR(I)
clin-3/2009:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    给出粒子的权重dpertp(i)
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
         dptemp(I) = dpertp(I)
 1001 CONTINUE
      CALL ARINDX(MAXSTR, NP, FT0, INDX)
      DO 1002 I = 1, NP
cbz12/3/98
c         IF (ITYP0(INDX(I)) .EQ. 211) THEN
c         IF (ITYP0(INDX(I)) .EQ. 211 .OR. ITYP0(INDX(I)) .EQ. 321) THEN
c         IF (ITYP0(INDX(I)) .EQ. 211 .OR. ITYP0(INDX(I)) .EQ. 2212 .OR.
c     &      ITYP0(INDX(I)) .EQ. 2112 .OR. ITYP0(INDX(I)) .EQ. -211 .OR.
c     &      ITYP0(INDX(I)) .EQ. 111) THEN
c         IF (ITYP0(INDX(I)) .EQ. 211 .OR. ITYP0(INDX(I)) .EQ. 2212 .OR.
c     &      ITYP0(INDX(I)) .EQ. 2112) THEN
         NPAR = NPAR + 1
c         ITYPAR(I) = ITYP0(INDX(I))
c         GXAR(I) = GX0(INDX(I))
c         GYAR(I) = GY0(INDX(I))
c         GZAR(I) = GZ0(INDX(I))
c         FTAR(I) = FT0(INDX(I))
c         PXAR(I) = PX0(INDX(I))
c         PYAR(I) = PY0(INDX(I))
c         PZAR(I) = PZ0(INDX(I))
c         PEAR(I) = EE0(INDX(I))
c         XMAR(I) = XM0(INDX(I))
         ITYPAR(NPAR) = ITYP0(INDX(I))
         GXAR(NPAR) = GX0(INDX(I))
         GYAR(NPAR) = GY0(INDX(I))
         GZAR(NPAR) = GZ0(INDX(I))
         FTAR(NPAR) = FT0(INDX(I))
         PXAR(NPAR) = PX0(INDX(I))
         PYAR(NPAR) = PY0(INDX(I))
         PZAR(NPAR) = PZ0(INDX(I))
         PEAR(NPAR) = EE0(INDX(I))
         XMAR(NPAR) = XM0(INDX(I))
clin-3/2009:
         dpertp(NPAR)=dptemp(INDX(I))
c         END IF
cbz12/3/98end
 1002 CONTINUE
      IAINT2(1) = NPAR
c
      RETURN
      END
c-----------------------------------------------------------------------
c.....subroutine to copy individually generated particle record into
c.....particle record for many test particle runs.
