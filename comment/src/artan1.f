      SUBROUTINE ARTAN1
      PARAMETER (MAXSTR=150001, MAXR=1)
c.....y cut for mt spectrum
cbz3/17/99
c      PARAMETER (YMT1 = -0.4, YMT2 = 0.4)
      PARAMETER (YMT1 = -1.0, YMT2 = 1.0)
cbz3/17/99 end
c.....bin width for mt spectrum and y spectrum
clin-9/26/03 no symmetrization in y (or eta) for ana/*.dat:
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
c     &     dm1k0s(50), DMT1LA(50), DMT1LB(50)
cc      SAVE /ARPRC1/
      COMMON /ARANA1/
     &     dy1ntb(50), dy1ntp(50), DY1HM(50), 
     &     DY1KP(50), DY1KM(50), DY1K0S(50),
     &     DY1LA(50), DY1LB(50), DY1PHI(50),
     &     dm1pip(50), dm1pim(50), DMT1PR(50),
     &     DMT1PB(50), DMT1KP(50), dm1km(50),
     &     dm1k0s(50), DMT1LA(50), DMT1LB(50),
     &     dy1msn(50), DY1PIP(50), DY1PIM(50), 
     &     DY1PI0(50), DY1PR(50), DY1PB(50)
     &     ,DY1NEG(50), DY1CH(50), DE1NEG(50), DE1CH(50)
cc      SAVE /ARANA1/
      SAVE   
cbz3/17/99 end
      DO 1002 J = 1, NUM
         DO 1001 I = 1, MULTI1(J)
            ITYP = ITYP1(I, J)
            PX = PX1(I, J)
            PY = PY1(I, J)
            PZ = PZ1(I, J)
            EE = EE1(I, J)
            XM = XM1(I, J)
c     2/24/03 leptons and photons:
            if(xm.lt.0.01) goto 200
            ptot = sqrt(PX ** 2 + PY ** 2 + pz ** 2)
clin-9/2012 determine pseudo-rapidity more generally:
c            eta = 0.5*alog((Ptot+pz+1e-5)/(ptot-pz+1e-5))


            if((PX**2+PY**2).gt.0.) then
               eta=asinh(PZ/sqrt(PX**2+PY**2))
            else
               eta = 1000000.0*sign(1.,PZ)
clin-2/2013 for spectator target nucleons in LAB frame,
c     note that finite precision of HBOOST 
c     would give spectator target nucleons a small but non-zero pz:
               if(abs(pz).le.1e-3) eta=0.
            endif

            
            XMT = SQRT(PX ** 2 + PY ** 2 + XM ** 2)
            DXMT = XMT - XM
clin-9/2012 determine rapidity more generally:
c            IF (ABS(PZ) .GE. EE) THEN
c               PRINT *, 'IN ARTAN1'
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
c            else
cc            Y = 0.5 * LOG((EE + PZ) / (EE - PZ))
c               Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ +1e-5))
cc               STOP
ccbzdbg2/15/99end
c            END IF


            if(XMT.gt.0.) then
               Y=asinh(PZ/XMT)
            else
               PRINT *, ' IN ARTAN1 mt=0'
               Y = 1000000.0*sign(1.,PZ)
            endif
c.....rapidity cut for the rapidity distribution


            IF (ABS(Y) .GE. 10.0) GOTO 100
clin-9/26/03 no symmetrization in y (or eta) for ana/*.dat:
c            IY = 1 + int(ABS(Y) / BY)
c            Ieta = 1 + int(ABS(eta) / BY)
            IF (ABS(eta) .GE. 10.0) GOTO 100
            IY = 1 + int((Y+10.) / BY)
            Ieta = 1 + int((eta+10.) / BY)


            IF (ITYP .LT. -1000) THEN
               dy1ntb(IY) = dy1ntb(IY) - 1.0
            END IF


            IF (ITYP .GT. 1000) THEN
               dy1ntb(IY) = dy1ntb(IY) + 1.0
            END IF


            IF (ITYP .EQ. -2212) THEN
               dy1ntp(IY) = dy1ntp(IY) - 1.0
            END IF


            IF (ITYP .EQ. 2212) THEN
               dy1ntp(IY) = dy1ntp(IY) + 1.0
            END IF
c            IF (ITYP .EQ. -211 .OR. ITYP .EQ. -321 .OR.
c     &         ITYP .EQ. -2212) THEN


            IF (ITYP .EQ. -2112) THEN
               DY1HM(IY) = DY1HM(IY) + 1.0
            END IF
c


            IF (LUCHGE(ITYP).ne.0) THEN
               DY1CH(IY) = DY1CH(IY) + 1.0
               DE1CH(Ieta) = DE1CH(Ieta) + 1.0
               IF (LUCHGE(ITYP).lt.0) THEN
                  DY1NEG(IY) = DY1NEG(IY) + 1.0
                  DE1NEG(Ieta) = DE1NEG(Ieta) + 1.0
               endif
            END IF
cbz3/17/99


            IF ((ITYP .GE. 100 .AND. ITYP .LT. 1000) .OR. 
     &         (ITYP .GT. -1000 .AND. ITYP .LE. -100)) THEN
               dy1msn(IY) = dy1msn(IY) + 1.0
            END IF


            IF (ITYP .EQ. 211) THEN
               DY1PIP(IY) = DY1PIP(IY) + 1.0
            END IF


            IF (ITYP .EQ. -211) THEN
               DY1PIM(IY) = DY1PIM(IY) + 1.0
            END IF


            IF (ITYP .EQ. 111) THEN
               DY1PI0(IY) = DY1PI0(IY) + 1.0
            END IF


            IF (ITYP .EQ. 2212) THEN
               DY1PR(IY) = DY1PR(IY) + 1.0
            END IF


            IF (ITYP .EQ. -2212) THEN
               DY1PB(IY) = DY1PB(IY) + 1.0
            END IF
cbz3/17/99 end


            IF (ITYP .EQ. 321) THEN
               DY1KP(IY) = DY1KP(IY) + 1.0
            END IF


            IF (ITYP .EQ. -321) THEN
               DY1KM(IY) = DY1KM(IY) + 1.0
            END IF
clin-4/24/03 evaluate K0L instead of K0S, since sometimes we may decay K0S:
c            IF (ITYP .EQ. 310) THEN


            IF (ITYP .EQ. 130) THEN
               DY1K0S(IY) = DY1K0S(IY) + 1.0
            END IF


            IF (ITYP .EQ. 3122) THEN
               DY1LA(IY) = DY1LA(IY) + 1.0
            END IF


            IF (ITYP .EQ. -3122) THEN
               DY1LB(IY) = DY1LB(IY) + 1.0
            END IF


            IF (ITYP .EQ. 333) THEN
               DY1PHI(IY) = DY1PHI(IY) + 1.0
            END IF
c.....insert rapidity cut for mt spectrum here


 100        IF (Y .LT. YMT1 .OR. Y .GT. YMT2) GOTO 200
            IF (DXMT .GE. 50.0 * BMT .OR. DXMT .EQ. 0) GOTO 200
            IMT = 1 + int(DXMT / BMT)


            IF (ITYP .EQ. 211) THEN
               dm1pip(IMT) = dm1pip(IMT) + 1.0 / XMT
            END IF


            IF (ITYP .EQ. -211) THEN
               dm1pim(IMT) = dm1pim(IMT) + 
     &            1.0 / XMT
            END IF


            IF (ITYP .EQ. 2212) THEN
               DMT1PR(IMT) = DMT1PR(IMT) + 1.0 / XMT
            END IF



            IF (ITYP .EQ. -2212) THEN
               DMT1PB(IMT) = DMT1PB(IMT) + 1.0 / XMT
            END IF


            IF (ITYP .EQ. 321) THEN
               DMT1KP(IMT) = DMT1KP(IMT) + 1.0 / XMT
            END IF

            

            IF (ITYP .EQ. -321) THEN
               dm1km(IMT) = dm1km(IMT) + 1.0 / XMT
            END IF
clin-4/24/03:
c            IF (ITYP .EQ. 310) THEN


            IF (ITYP .EQ. 130) THEN
               dm1k0s(IMT) = dm1k0s(IMT) + 1.0 / XMT
            END IF


            
            IF (ITYP .EQ. 3122) THEN
               DMT1LA(IMT) = DMT1LA(IMT) + 1.0 / XMT
            END IF


            
            IF (ITYP .EQ. -3122) THEN
               DMT1LB(IMT) = DMT1LB(IMT) + 1.0 / XMT
            END IF
 200        CONTINUE
 1001    CONTINUE
 1002 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
c.....analysis subroutine after the hadronic space-time evolution
