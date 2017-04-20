      subroutine addhad
      PARAMETER (MAXSTR=150001,MAXR=1,xmd=1.8756)
      double precision  smearp,smearh
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      COMMON /smearz/smearp,smearh
      COMMON/RNDF77/NSEED
      common /para8/ idpert,npertd,idxsec
      SAVE   
c     All hadrons at the start of hadron cascade have the weight of 1
c     except those inserted by the user in this subroutine:
      np0=IAINT2(1)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      初始所有的强子的cascade权重为1,dpertp(I)貌似就应该是权重了.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DO i=1,np0
         dpertp(I)=1.
      ENDDO
c     Specify number, species, weight, initial x,p,m for inserted hadrons here:
      nadd=0
      tau0=ARPAR1(1)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    对于插入的强子,给出ITYPAR(42)是什么鬼???????
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     下面的循环永远都不进行.所以注释掉了2017/4/19
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      DO 100 i=np0+1,np0+nadd
c$$$         ITYPAR(I)=42
c$$$clin-5/2012 fix type mismatch:
c$$$c         dpertp(I)=1d0/dble(nadd)
c$$$         dpertp(I)=1./float(nadd)
c$$$         GXAR(I)=5.*(1.-2.*RANART(NSEED))
c$$$         GYAR(I)=5.*(1.-2.*RANART(NSEED))
c$$$         GZAR(I)=2.*(1.-2.*RANART(NSEED))
c$$$         FTAR(I)=0.
c$$$         PXAR(I)=1.
c$$$         PYAR(I)=0.
c$$$         PZAR(I)=1.
c$$$         XMAR(I)=xmd
c$$$c
c$$$         PEAR(I)=sqrt(PXAR(I)**2+PYAR(I)**2+PZAR(I)**2+XMAR(I)**2)
c$$$clin-9/2012 determine rapidity more generally:
c$$$c         RAP=0.5*alog((PEAR(I)+PZAR(I)+1e-5)/(PEAR(I)-PZAR(I)+1e-5))
c$$$         RAP=asinh(PZAR(I)/sqrt(XMAR(I)**2+PXAR(I)**2+PYAR(I)**2))
c$$$c
c$$$         VX=PXAR(I)/PEAR(I)
c$$$         VY=PYAR(I)/PEAR(I)
c$$$c.....give initial formation time shift and boost according to rapidity:
c$$$         TAUI=FTAR(I)+TAU0
c$$$         FTAR(I)=TAUI*COSH(RAP)
c$$$         GXAR(I)=GXAR(I)+VX*TAU0*COSH(RAP)
c$$$         GYAR(I)=GYAR(I)+VY*TAU0*COSH(RAP)
c$$$c     Allow the intial z-position to be different from the Bjorken picture:
c$$$         GZAR(I)=TAUI*SINH(RAP)+GZAR(I)
c$$$c         GZAR(I)=TAUI*SINH(RAP)
c$$$         zsmear=sngl(smearh)*(2.*RANART(NSEED)-1.)
c$$$         GZAR(I)=GZAR(I)+zsmear
c$$$ 100  CONTINUE
      
      IAINT2(1)=IAINT2(1)+nadd
c
      if(nadd.ge.1.and.idpert.ne.1.and.idpert.ne.2) then
         write(16,*) 'IDPERT must be 1 or 2 to add initial hadrons,
     1 set NPERTD to 0 if you do not need perturbative deuterons'
         stop
      endif
      if(IAINT2(1).gt.MAXSTR) then
         write(16,*) 'Too many initial hadrons, array size is exceeded!'
         stop
      endif
c
      return
      end
clin-8/2014 define function asinh():
