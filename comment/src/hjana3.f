      SUBROUTINE HJANA3
c
      PARAMETER (MAXSTR=150001, MAXR=1)
c.....y cut for mt spectrum
      PARAMETER (YMIN = -1.0, YMAX = 1.0)
cbz11/7/99 end
c.....bin width for mt spectrum and y spectrum
      PARAMETER (DMT = 0.05, DY = 0.2)
      DOUBLE PRECISION v2i,eti,xmulti,v2mi,s2mi,xmmult,
     1     v2bi,s2bi,xbmult
      DIMENSION dndyh3(50), DMYH3(50), DEYH3(50)
      COMMON /RUN/ NUM
cc      SAVE /RUN/
      COMMON /ARERC1/MULTI1(MAXR)
cc      SAVE /ARERC1/
      COMMON /ARPRC1/ITYP1(MAXSTR, MAXR),
     &     GX1(MAXSTR, MAXR), GY1(MAXSTR, MAXR), GZ1(MAXSTR, MAXR), 
     &     FT1(MAXSTR, MAXR),
     &     PX1(MAXSTR, MAXR), PY1(MAXSTR, MAXR), PZ1(MAXSTR, MAXR),
     &     EE1(MAXSTR, MAXR), XM1(MAXSTR, MAXR)
cc      SAVE /ARPRC1/
      COMMON /AROUT/ IOUT
cc      SAVE /AROUT/
      COMMON/iflow/v2i,eti,xmulti,v2mi,s2mi,xmmult,v2bi,s2bi,xbmult
cc      SAVE /iflow/
      SAVE   
      DATA IW/0/
      IW = IW + 1
      DO 1002 J = 1, NUM
         DO 1001 I = 1, MULTI1(J)
            ITYP = ITYP1(I, J)
            IF (ITYP .GT. -100 .AND. ITYP .LT. 100) GOTO 200
            PX = PX1(I, J)
            PY = PY1(I, J)
            PZ = PZ1(I, J)
            EE = EE1(I, J)
            XM = XM1(I, J)
            XMT = SQRT(PX ** 2 + PY ** 2 + XM ** 2)
            DXMT = XMT - XM
clin-9/2012 determine rapidity more generally:
c            IF (ABS(PZ) .GE. EE) THEN
c               PRINT *, 'IN HJANA3'
c               PRINT *, ' PARTICLE ', I, ' RUN ', J, 'PREC ERR'
c               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
c               PRINT *, ' PZ = ', PZ, ' EE = ', EE
c               PRINT *, ' XM = ', XM
c               GOTO 200
c            END IF
c            Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ + 1e-5))
            if(XMT.gt.0.) then
               Y=asinh(PZ/XMT)
            else
               PRINT *, ' IN HJANA3 mt=0'
               Y = 1000000.0*sign(1.,PZ)
            endif
c.....rapidity cut for the rapidity distribution
c            IY = 1 + int(ABS(Y) / DY)
clin-8/2014 no rapidity shift here:
c            IY = 1 + int((Y+10.) / DY)
            IY = 1 + int(Y/DY)
clin-9/2012 prevent possible segmentation fault (due to IY<=0):
c            IF (IY .GT. 50) GOTO 100
            IF (IY.lt.1 .or.IY .GT. 50) GOTO 100
            dndyh3(IY) = dndyh3(IY) + 1.0
            DEYH3(IY) = DEYH3(IY) + XMT
 100        CONTINUE
c.....insert rapidity cut for mt spectrum here
            IF (Y. LT. YMIN .OR. Y .GE. YMAX) GOTO 200
            IMT = 1 + int(DXMT / DMT)
            IF (IMT .GT. 50) GOTO 200
            DMYH3(IMT) = DMYH3(IMT) + 1.0 / XMT
 200        CONTINUE
 1001    CONTINUE
 1002 CONTINUE
c
      RETURN
      END
c-----------------------------------------------------------------------
c.....analysis subroutine after ARTMN
