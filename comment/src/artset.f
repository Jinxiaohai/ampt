c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        设定了一些参数和打开了oscar的几个标准的文件进行后面的输出。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      SUBROUTINE ARTSET
c
      PARAMETER (AMU= 0.9383,nxymax=10001)
      double precision dpcoal,drcoal,ecritl
      INTEGER ZTA, ZPR
      common  /gg/      dx,dy,dz,dpx,dpy,dpz
clin-10/03/03 
c     "SAVE   " (without argument) is used for most subroutines and functions,
c     this is important for the success when using "f77" to compile:
cc      SAVE /gg/
      common  /zz/      zta,zpr
cc      SAVE /zz/
      COMMON  /RUN/     NUM
cc      SAVE /RUN/
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
cc      SAVE /INPUT2/
      COMMON /INPUT3/ PLAB, ELAB, ZEROPT, B0, BI, BM, DENCUT, CYCBOX
cc      SAVE /INPUT3/
      common /imulst/ iperts
cc      SAVE /imulst/
      common /coal/dpcoal,drcoal,ecritl
      common/anim/nevent,isoft,isflag,izpc
      common /para7/ ioscar,nsmbbbar,nsmmeson
      common/embed/iembed,nsembd,pxqembd,pyqembd,xembd,yembd,
     1     psembd,tmaxembd,phidecomp
      common/xyembed/nxyjet,xyjet(nxymax,2)
      SAVE   
clin-10/03/03  ecritl: local energy density below which a parton 
c     will freeze out (in GeV/fm^3), for improvements on string melting, 
c     not used in this version of AMPT:
clin-4/2008
c      data ecritl/1.d0/
      ecritl=1.d0
c
c     combine ART initialization into ampt.ini:
c     (Note that the following values are relics from the old ART structure)
c.....input parameter file
c      OPEN(13, FILE = 'art1.ini', STATUS = 'UNKNOWN')
c      READ (13, *) MASSTA, ZTA
      MASSTA=1
      ZTA=1
c      write(12,*) massta, zta, ' massta, zta'
c      READ (13, *) MASSPR, ZPR
      MASSPR=1
      ZPR=1
c      write(12,*) masspr, zpr, ' masspr, zpr'
c      READ (13, *) PLAB, IPLAB
      PLAB=14.6 
      IPLAB=2
c      write(12,*) plab, iplab, ' plab, iplab'
      if(iplab.eq.2)then
         elab=sqrt(plab**2+amu**2)-amu
      else
         elab=plab
      endif
      elab=elab*1000.
c      READ (13, *) ZEROPT
      ZEROPT=0.
c      write(12,*) zeropt, ' zeropt'
clin-10/03/03 ISEED was used as a seed for random number inside ART, 
c     not used in AMPT:
      ISEED=700721
c     0/1: (Normal or Perturbative) multistrange partice production.
c     Perturbative option is disabled for now:
      iperts=0
c      READ (13, *) MANYB, B0, BI, BM
c     2/04/00 MANYB MUST BE SET TO 1 !
c     in order to skip impact parameter setting by ART, then B0 has no effect.
      MANYB=1
      B0=1
      BI=0
      BM=0
c      write(12,*) manyb, b0, bi, bm, ' manyb, b0, bi, bm'
c      READ (13, *) ISEED
c      write(12,*) iseed, ' iseed'
c      READ (13, *) DT
c      write(12,*) dt, ' dt'
c      READ (13, *) NTMAX
c      write(12,*) ntmax, ' ntmax'
c      READ (13, *) ICOLL
      ICOLL=-1
c      write(12,*) icoll, ' icoll'
c      READ (13, *) NUM
c     2/11/03 run events without test particles for now:
      NUM=1
c      write(12,*) num, ' num'
c      READ (13, *) INSYS
      INSYS=1
c      write(12,*) insys, ' insys'
c      READ (13, *) IPOT
      IPOT=3
c      write(12,*) ipot, ' ipot'
c      READ (13, *) MODE
      MODE=0
      IF(ICOLL.EQ.-1)IPOT=0
c      write(12,*) mode, ' mode'
c      READ (13, *) DX, DY, DZ
      DX=2.73
      DY=2.73
      DZ=2.73
c      write(12,*) dx,dy,dz,' dx,dy,dz'
c      READ (13, *) DPX, DPY, DPZ
      DPX=0.6
      DPY=0.6
      DPZ=0.6
c      write(12,*) dpx,dpy,dpz,' dpx,dpy,dpz'
c      READ (13, *) IAVOID
      IAVOID=1
c      write(12,*) iavoid, ' iavoid'
c      READ (13, *) IMOMEN
      IMOMEN=1
c      write(12,*) imomen, ' imomen'
      if(icoll.eq.-1)imomen=3
c      READ (13, *) NFREQ
      NFREQ=10
c      write(12,*) nfreq, ' nfreq'
c      READ (13, *) ICFLOW
      ICFLOW=0
c      write(12,*) ICFLOW, ' ICFLOW'
c      READ (13, *) ICRHO
      ICRHO=0
c      write(12,*) ICRHO, ' ICRHO'
c      READ (13, *) ICOU
      ICOU=0
c      write(12,*)icou, ' icou'
* kaon potential control parameter
* KMUL IS A MULTIPLIER TO THE STANDARD K-N SCATTERING LENGTH
c      READ (13, *) KPOTEN, KMUL
      KPOTEN=0
      KMUL=1
c      write(12,*)kpoten,kmul, ' kpoten, kmul'
* mean field control parameter FOR BARYONS
* no mean filed is used for baryons if their 
* local density is higher than dencut. 
c      READ (13, *) DENCUT
      DENCUT=15
c      write(12,*)dencut, ' dencut'
* test reactions in a box of side-length cycbox
* input cycbox
c      READ (13, *) CYCBOX
      CYCBOX=0
c      write(12,*) cycbox, ' cycbox'
c
clin-5b/2008
c      if(ioscar.eq.2) then
      if(ioscar.eq.2.or.ioscar.eq.3) then
         OPEN (92,FILE='ana/parton-initial-afterPropagation.dat',
     1        STATUS = 'UNKNOWN')
      endif
      if(ioscar.eq.3) then
clin-6/2009 write out full parton collision history:
         OPEN (95,FILE='ana/parton-collisionsHistory.dat',
     1        STATUS='UNKNOWN')
clin-6/2009 write out initial minijet information:
         OPEN (96,FILE='ana/minijet-initial-beforePropagation.dat',
     1        STATUS='UNKNOWN')
clin-6/2009 write out parton info after coalescence:
         if(isoft.eq.4.or.isoft.eq.5) then
            OPEN (85,FILE='ana/parton-after-coalescence.dat',
     1           STATUS='UNKNOWN')
         endif
      endif
clin-6/2009 write out initial transverse positions of initial nucleons:
      OPEN (94,FILE='ana/npart-xy.dat',STATUS='UNKNOWN')
c
clin-8/2009 In case that random positions are used to embed high-Pt jets:
      if(iembed.eq.3.or.iembed.eq.4) then
         OPEN (97,FILE='embed-jet-xy.txt',STATUS = 'UNKNOWN')
         read(97,*) nxyjet
c     Save positions in array to reuse when embedding more jet pairs 
c     than the number of entries in the position file:
         if(nevent.gt.nxyjet) then
            if(nxyjet.gt.nxymax) then
               print *, 'Too many lines in embed-jet-xy.txt: 
     1 increase value of the parameter nxymax'
               stop
            elseif(nxyjet.le.0) then
               print *, 'Check number of entries in embed-jet-xy.txt' 
               stop
            endif
            do ixy=1,nxyjet
               read(97,*) xyjet(ixy,1),xyjet(ixy,2)
            enddo
         endif
      endif
      RETURN
      END
c-----------------------------------------------------------------------
c.....subroutine to initialize cascade.
