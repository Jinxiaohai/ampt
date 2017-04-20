      SUBROUTINE ARTAN2
      PARAMETER (MAXSTR=150001, MAXR=1)
c.....y cut for mt spectrum
cbz3/17/99
c      PARAMETER (YMT1 = -0.4, YMT2 = 0.4)
      PARAMETER (YMT1 = -1.0, YMT2 = 1.0)
cbz3/17/99 end
c.....bin width for mt spectrum and y spectrum
c      PARAMETER (BMT = 0.05, BY = 0.2)
      PARAMETER (BMT = 0.05, BY = 0.4)
      COMMON /RUN/ NUM
cc      SAVE /RUN/
      COMMON /ARERC1/MULTI1(MAXR)
cc      SAVE /ARERC1/
      COMMON /ARPRC1/ITYP1(MAXSTR, MAXR),
     &     GX1(MAXSTR, MAXR), GY1(MAXSTR, MAXR), GZ1(MAXSTR, MAXR), 
     &     FT1(MAXSTR, MAXR),
     &     PX1(MAXSTR, MAXR), PY1(MAXSTR, MAXR), PZ1(MAXSTR, MAXR),
     &     EE1(MAXSTR, MAXR), XM1(MAXSTR, MAXR)
cbz3/17/99
c     &     dm2k0s(50), DMT2LA(50), DMT2LB(50)
cc      SAVE /ARPRC1/
      COMMON /ARANA2/
     &     dy2ntb(50), dy2ntp(50), DY2HM(50), 
     &     DY2KP(50), DY2KM(50), DY2K0S(50),
     &     DY2LA(50), DY2LB(50), DY2PHI(50),
     &     dm2pip(50), dm2pim(50), DMT2PR(50),
     &     DMT2PB(50), DMT2KP(50), dm2km(50),
     &     dm2k0s(50), DMT2LA(50), DMT2LB(50),
     &     dy2msn(50), DY2PIP(50), DY2PIM(50), 
     &     DY2PI0(50), DY2PR(50), DY2PB(50)
     &     ,DY2NEG(50), DY2CH(50), DE2NEG(50), DE2CH(50)
cbz3/17/99 end
cc      SAVE /ARANA2/
      SAVE   
      DO 1002 J = 1, NUM
         DO 1001 I = 1, MULTI1(J)
            ITYP = ITYP1(I, J)
            PX = PX1(I, J)
            PY = PY1(I, J)
            PZ = PZ1(I, J)
            EE = EE1(I, J)
            XM = XM1(I, J)
            XMT = SQRT(PX ** 2 + PY ** 2 + XM ** 2)
c     2/24/03 leptons and photons:
            if(xm.lt.0.01) goto 200
            DXMT = XMT - XM
            ptot = sqrt(PX ** 2 + PY ** 2 + pz ** 2)
clin-9/2012 determine pseudo-rapidity more generally:
c            eta = 0.5*alog((Ptot+pz+1e-5)/(ptot-pz+1e-5))
            if((PX**2+PY**2).gt.0.) then
               eta=asinh(PZ/sqrt(PX**2+PY**2))
            else
               eta = 1000000.0*sign(1.,PZ)
clin-2/2013 for spectator target nucleons in LAB frame:
               if(abs(pz).le.1e-3) eta=0.
            endif
clin-9/2012 determine rapidity more generally:
c            IF (ABS(PZ) .GE. EE) THEN
c               PRINT *, 'IN ARTAN2'
c               PRINT *, 'PARTICLE ', I, ' RUN ', J, 'PREC ERR'
ccbzdbg2/16/99
c               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
ccbzdbg2/16/99
ccbzdbg2/15/99
c               PRINT *, ' PZ = ', PZ, ' EE = ', EE
ccbzdbg2/16/99
c               PRINT *, ' XM = ', XM
ccbzdbg2/16/99end
c               GOTO 200
cc               STOP
ccbzdbg2/15/99end
c            END IF
c            Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ + 1e-5))
            if(XMT.gt.0.) then
               Y=asinh(PZ/XMT)
            else
               PRINT *, ' IN ARTAN2 mt=0'
               Y = 1000000.0*sign(1.,PZ)
            endif
c.....rapidity cut for the rapidity distribution
            IF (ABS(Y) .GE. 10.0) GOTO 100
c            IY = 1 + int(ABS(Y) / BY)
c            Ieta = 1 + int(ABS(eta) / BY)
            IF (ABS(eta) .GE. 10.0) GOTO 100
            IY = 1 + int((Y+10.) / BY)
            Ieta = 1 + int((eta+10.) / BY)
            IF (ITYP .LT. -1000) THEN
               dy2ntb(IY) = dy2ntb(IY) - 1.0
            END IF
            IF (ITYP .GT. 1000) THEN
               dy2ntb(IY) = dy2ntb(IY) + 1.0
            END IF
            IF (ITYP .EQ. -2212) THEN
               dy2ntp(IY) = dy2ntp(IY) - 1.0
            END IF
            IF (ITYP .EQ. 2212) THEN
               dy2ntp(IY) = dy2ntp(IY) + 1.0
            END IF
            IF (ITYP .EQ. -2112) THEN
               DY2HM(IY) = DY2HM(IY) + 1.0
            END IF
            IF (LUCHGE(ITYP).ne.0) THEN
               DY2CH(IY) = DY2CH(IY) + 1.0
               DE2CH(Ieta) = DE2CH(Ieta) + 1.0
               IF (LUCHGE(ITYP).lt.0) THEN
                  DY2NEG(IY) = DY2NEG(IY) + 1.0
                  DE2NEG(Ieta) = DE2NEG(Ieta) + 1.0
               endif
            END IF
cbz3/17/99
            IF ((ITYP .GE. 100 .AND. ITYP .LT. 1000) .OR. 
     &         (ITYP .GT. -1000 .AND. ITYP .LE. -100)) THEN
               dy2msn(IY) = dy2msn(IY) + 1.0
            END IF
            IF (ITYP .EQ. 211) THEN
               DY2PIP(IY) = DY2PIP(IY) + 1.0
            END IF
            IF (ITYP .EQ. -211) THEN
               DY2PIM(IY) = DY2PIM(IY) + 1.0
            END IF
            IF (ITYP .EQ. 111) THEN
               DY2PI0(IY) = DY2PI0(IY) + 1.0
            END IF
            IF (ITYP .EQ. 2212) THEN
               DY2PR(IY) = DY2PR(IY) + 1.0
            END IF
            IF (ITYP .EQ. -2212) THEN
               DY2PB(IY) = DY2PB(IY) + 1.0
            END IF
cbz3/17/99 end
            IF (ITYP .EQ. 321) THEN
               DY2KP(IY) = DY2KP(IY) + 1.0
            END IF
            IF (ITYP .EQ. -321) THEN
               DY2KM(IY) = DY2KM(IY) + 1.0
            END IF
clin-4/24/03:
c            IF (ITYP .EQ. 310) THEN
            IF (ITYP .EQ. 130) THEN
               DY2K0S(IY) = DY2K0S(IY) + 1.0
            END IF
            IF (ITYP .EQ. 3122) THEN
               DY2LA(IY) = DY2LA(IY) + 1.0
            END IF
            IF (ITYP .EQ. -3122) THEN
               DY2LB(IY) = DY2LB(IY) + 1.0
            END IF
            IF (ITYP .EQ. 333) THEN
               DY2PHI(IY) = DY2PHI(IY) + 1.0
            END IF
c.....insert rapidity cut for mt spectrum here
 100        IF (Y .LT. YMT1 .OR. Y .GT. YMT2) GOTO 200
            IF (DXMT .GE. 50.0 * BMT .OR. DXMT .EQ. 0) GOTO 200
            IMT = 1 + int(DXMT / BMT)
            IF (ITYP .EQ. 211) THEN
               dm2pip(IMT) = dm2pip(IMT) + 1.0 / XMT
            END IF
            IF (ITYP .EQ. -211) THEN
               dm2pim(IMT) = dm2pim(IMT) + 
     &            1.0 / XMT
            END IF
            IF (ITYP .EQ. 2212) THEN
               DMT2PR(IMT) = DMT2PR(IMT) + 1.0 / XMT
            END IF
            IF (ITYP .EQ. -2212) THEN
               DMT2PB(IMT) = DMT2PB(IMT) + 1.0 / XMT
            END IF
            IF (ITYP .EQ. 321) THEN
               DMT2KP(IMT) = DMT2KP(IMT) + 1.0 / XMT
            END IF
            IF (ITYP .EQ. -321) THEN
               dm2km(IMT) = dm2km(IMT) + 1.0 / XMT
            END IF
clin-4/24/03:
c            IF (ITYP .EQ. 310) THEN
            IF (ITYP .EQ. 130) THEN
               dm2k0s(IMT) = dm2k0s(IMT) + 1.0 / XMT
            END IF
            IF (ITYP .EQ. 3122) THEN
               DMT2LA(IMT) = DMT2LA(IMT) + 1.0 / XMT
            END IF
            IF (ITYP .EQ. -3122) THEN
               DMT2LB(IMT) = DMT2LB(IMT) + 1.0 / XMT
            END IF
 200        CONTINUE
 1001    CONTINUE
 1002 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
c.....output analysis results at the end of the simulation
