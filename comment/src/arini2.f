      SUBROUTINE ARINI2(K)
      PARAMETER (MAXSTR=150001,MAXR=1)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
cc      SAVE /ARPRNT/
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
cc      SAVE /ARPRC/
      COMMON /ARERC1/MULTI1(MAXR)
cc      SAVE /ARERC1/
      COMMON /ARPRC1/ITYP1(MAXSTR, MAXR),
     &     GX1(MAXSTR, MAXR), GY1(MAXSTR, MAXR), GZ1(MAXSTR, MAXR), 
     &     FT1(MAXSTR, MAXR),
     &     PX1(MAXSTR, MAXR), PY1(MAXSTR, MAXR), PZ1(MAXSTR, MAXR),
     &     EE1(MAXSTR, MAXR), XM1(MAXSTR, MAXR)
cc      SAVE /ARPRC1/
      COMMON/tdecay/tfdcy(MAXSTR),tfdpi(MAXSTR,MAXR),tft(MAXSTR)
cc      SAVE /tdecay/
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &     IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
cc      SAVE /INPUT2/
      COMMON/RNDF77/NSEED
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
cc      SAVE /RNDF77/
      SAVE   
      MULTI1(K) = IAINT2(1)
      DO 1001 I = 1, MULTI1(K)
         ITYP1(I, K) = ITYPAR(I)
         GX1(I, K) = GXAR(I)
         GY1(I, K) = GYAR(I)
         GZ1(I, K) = GZAR(I)
         FT1(I, K) = FTAR(I)
         PX1(I, K) = PXAR(I)
         PY1(I, K) = PYAR(I)
         PZ1(I, K) = PZAR(I)
         EE1(I, K) = PEAR(I)
         XM1(I, K) = XMAR(I)
clin-3/2009 hadron weights are initialized in addhad():
clin-5/2008 all hadrons not perturbatively-produced have the weight of 1:
c         dpp1(I,K)=1.
         dpp1(I,K)=dpertp(I)
 1001 CONTINUE
c     initialize final time of each particle to ntmax*dt except for 
c     decay daughters, which have values given by tfdcy() and >(ntmax*dt):
      do 1002 ip=1,MAXSTR
         tfdcy(ip)=NTMAX*DT
         tft(ip)=NTMAX*DT
 1002 continue
c
      do 1004 irun=1,MAXR
         do 1003 ip=1,MAXSTR
            tfdpi(ip,irun)=NTMAX*DT
 1003    continue
 1004 continue
      RETURN
      END
c=======================================================================
c.....function to convert PDG flavor code into ART flavor code.
