      SUBROUTINE ARINI1
c.....before invoking ARINI1:
c.....ARPAR1(1), IAINT2(1) must be set:
      PARAMETER (MAXSTR=150001)
      double precision  smearp,smearh
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
cc      SAVE /ARPRNT/
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
cc      SAVE /ARPRC/
      COMMON /smearz/smearp,smearh
cc      SAVE /smearz/
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
      common /nzpc/nattzp
cc      SAVE /nzpc/
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      common /para8/ idpert,npertd,idxsec
      SAVE   
clin-5/2008 for perturbatively-produced hadrons (currently only deuterons):
      OPEN (91, FILE = 'ana/deuteron_processes.dat', 
     1     STATUS = 'UNKNOWN')
      if(idpert.eq.1.or.idpert.eq.2) then
         OPEN (90, FILE = 'ana/ampt_pert.dat', STATUS = 'UNKNOWN')
      endif
c.....generate formation time and position at formation time.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    利用"信息时间"生成信息时间和位置
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
         write(9969,*)"TAU0 = ",ARPAR1(1), "NP = ", IAINT2(1)
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      TAU0 = ARPAR1(1)
      NP = IAINT2(1)
clin-7/10/01     initial positions already given for hadrons 
c     formed from partons inside ZPC (from string melting):
      if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
clin-8/2015 fixed a bug that may skip "dpertp(I)=1." in addhad and
c     cause the first few events to be missing in ampt.dat 
c     (mostly for low-multiplicity events such as PP collisions):
c         if(NP.le.nattzp) return
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
         write(9969, *)"NP = ", NP, "nattzp =", nattzp
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
         if(NP.gt.nattzp) then
         do 1001 I = nattzp+1, NP
clin-9/2012 determine rapidity more generally 
c     to prevent overflow when Pt~=0 and E=|Pz|:
c            IF (ABS(PZAR(I)) .GE. PEAR(I)) THEN
c               PRINT *, ' IN ARINI1'
c               PRINT *, 'ABS(PZ) .GE. EE for particle ', I
c               PRINT *, ' FLAV = ', ITYPAR(I), ' PX = ', PXAR(I), 
c     &              ' PY = ', PYAR(I)
c               PRINT *, ' PZ = ', PZAR(I), ' EE = ', PEAR(I)
c               PRINT *, ' XM = ', XMAR(I)
c               RAP = 1000000.0
c               GOTO 50
c            END IF
cc            RAP=0.5*LOG((PEAR(I)+PZAR(I))/(PEAR(I)-PZAR(I)))
c            RAP=0.5*LOG((PEAR(I)+PZAR(I)+1e-5)/(PEAR(I)-PZAR(I)+1e-5))
c 50         CONTINUE
            if((XMAR(I)**2+PXAR(I)**2+PYAR(I)**2).gt.0.) then
               RAP=asinh(PZAR(I)/sqrt(XMAR(I)**2+PXAR(I)**2+PYAR(I)**2))
            else
               PRINT *, ' IN ARINI1 mt=0'
               RAP = 1000000.0*sign(1.,PZAR(I))
            endif
            VX = PXAR(I) / PEAR(I)
            VY = PYAR(I) / PEAR(I)
            FTAR(I) = TAU0 * COSH(RAP)
            GXAR(I) = GXAR(I) + VX * FTAR(I)
            GYAR(I) = GYAR(I) + VY * FTAR(I)
            GZAR(I) = TAU0 * SINH(RAP)
clin-5/2009 No formation time for spectator projectile or target nucleons:
            if(PXAR(I).eq.0.and.PYAR(I).eq.0
     2           .and.(ITYPAR(I).eq.2112.or.ITYPAR(I).eq.2212)) then
clin-2/2013 for spectator target nucleons in LAB frame:
c     1           .and.(PEAR(I)*2/HINT1(1)).gt.0.99
              if((PEAR(I)/HINT1(6).gt.0.99.and.PEAR(I)/HINT1(6).lt.1.01)
     1 .or.(PEAR(I)/HINT1(7).gt.0.99.and.PEAR(I)/HINT1(7).lt.1.01)) then
c
               TAUI=1.E-20
               FTAR(I)=TAUI*COSH(RAP)
               GZAR(I)=TAUI*SINH(RAP)
               endif
            endif
 1001    continue
clin-8/2015:
         endif
clin-7/10/01-end
clin-3/2009 cleanup of program flow:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    下面的else是在default版本中使用的.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      else
c$$$         DO 1002 I = 1, NP
c$$$clin-9/2012 determine rapidity more generally:
c$$$c            IF (ABS(PZAR(I)) .GE. PEAR(I)) THEN
c$$$c               PRINT *, ' IN ARINI1'
c$$$c               PRINT *, 'ABS(PZ) .GE. EE for particle ', I
c$$$c               PRINT *, ' FLAV = ', ITYPAR(I), ' PX = ', PXAR(I), 
c$$$c     &              ' PY = ', PYAR(I)
c$$$c               PRINT *, ' PZ = ', PZAR(I), ' EE = ', PEAR(I)
c$$$c               PRINT *, ' XM = ', XMAR(I)
c$$$c               RAP = 1000000.0
c$$$c               GOTO 100
c$$$cc               STOP
c$$$c            END IF
c$$$c 100        CONTINUE
c$$$c            RAP=0.5*LOG((PEAR(I)+PZAR(I)+1e-5)/(PEAR(I)-PZAR(I)+1e-5))
c$$$            if((XMAR(I)**2+PXAR(I)**2+PYAR(I)**2).gt.0.) then
c$$$               RAP=asinh(PZAR(I)/sqrt(XMAR(I)**2+PXAR(I)**2+PYAR(I)**2))
c$$$            else
c$$$               PRINT *, ' IN ARINI1 mt=0'
c$$$               RAP = 1000000.0*sign(1.,PZAR(I))
c$$$            endif
c$$$            VX = PXAR(I) / PEAR(I)
c$$$            VY = PYAR(I) / PEAR(I)
c$$$c.....give initial formation time shift
c$$$            TAUI = FTAR(I) + TAU0
c$$$            FTAR(I) = TAUI * COSH(RAP)
c$$$            GXAR(I) = GXAR(I) + VX * TAU0 * COSH(RAP)
c$$$            GYAR(I) = GYAR(I) + VY * TAU0 * COSH(RAP)
c$$$c     4/25/03: hadron z-position upon formation determined the same way as x,y:
c$$$            GZAR(I) = TAUI * SINH(RAP)
c$$$c     the old prescription:
c$$$c            GZAR(I) = GZAR(I) + TAU0 * SINH(RAP)
c$$$            zsmear=sngl(smearh)*(2.*RANART(NSEED)-1.)
c$$$            GZAR(I)=GZAR(I)+zsmear
c$$$cbz1/28/99end
c$$$c     10/05/01 no formation time for spectator projectile or target nucleons:
c$$$            if(PXAR(I).eq.0.and.PYAR(I).eq.0
c$$$     2           .and.(ITYPAR(I).eq.2112.or.ITYPAR(I).eq.2212)) then
c$$$clin-2/2013 for spectator target nucleons in LAB frame:
c$$$c     1           .and.(PEAR(I)*2/HINT1(1)).gt.0.99
c$$$              if((PEAR(I)/HINT1(6).gt.0.99.and.PEAR(I)/HINT1(6).lt.1.01)
c$$$     1 .or.(PEAR(I)/HINT1(7).gt.0.99.and.PEAR(I)/HINT1(7).lt.1.01)) then
c$$$c
c$$$clin-5/2008:
c$$$c               TAUI=0.00001
c$$$               TAUI=1.E-20
c$$$               FTAR(I)=TAUI*COSH(RAP)
c$$$               GZAR(I)=TAUI*SINH(RAP)+zsmear
c$$$               endif
c$$$            endif
c$$$ 1002    CONTINUE
clin-3/2009 cleanup of program flow:
      endif
clin-3/2009 Add initial hadrons before the hadron cascade starts:
      call addhad
      RETURN
      END
c-----------------------------------------------------------------------
c.....subroutine to order particle labels according to increasing 
c.....formation time
