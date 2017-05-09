      subroutine embedHighPt
      PARAMETER (MAXSTR=150001,MAXR=1,pichmass=0.140,pi0mass=0.135,
     1     pi=3.1415926,nxymax=10001)
      common/embed/iembed,nsembd,pxqembd,pyqembd,xembd,yembd,
     1     psembd,tmaxembd,phidecomp
      COMMON/RNDF77/NSEED
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      common/anim/nevent,isoft,isflag,izpc
      COMMON /AREVT/ IAEVT, IARUN, MISS
      common/xyembed/nxyjet,xyjet(nxymax,2)
      SAVE
c
      if(iembed.eq.1.or.iembed.eq.2) then
         xjet=xembd
         yjet=yembd
      elseif(iembed.eq.3.or.iembed.eq.4) then
         if(nevent.le.nxyjet) then
            read(97,*) xjet,yjet
         else
            ixy=mod(IAEVT,nxyjet)
            if(ixy.eq.0) ixy=nxyjet
            xjet=xyjet(ixy,1)
            yjet=xyjet(ixy,2)
         endif
      else
         return
      endif
c
      ptq=sqrt(pxqembd**2+pyqembd**2)
      if(ptq.lt.(pichmass/2.)) then
         print *, 'Embedded quark transverse momentum is too small'
         stop
      endif
c     Randomly embed u/ubar or d/dbar at high Pt:
      idqembd=1+int(2*RANART(NSEED))
c     Flavor content for the charged pion that contains the leading quark:
      if(idqembd.eq.1) then 
         idqsoft=-2
         idpi1=-211
      elseif(idqembd.eq.2) then 
         idqsoft=-1
         idpi1=211
      else
         print *, 'Wrong quark flavor embedded'
         stop
      endif
c     Caculate transverse momentum of the parent charged pion:
      xmq=ulmass(idqembd)
      xmqsoft=ulmass(idqsoft)
      ptpi=((pichmass**2+xmq**2-xmqsoft**2)*ptq
     1     -sqrt((xmq**2+ptq**2)*(pichmass**4
     2     -2.*pichmass**2*(xmq**2+xmqsoft**2)+(xmq**2-xmqsoft**2)**2)))
     3     /(2.*xmq**2)
      if(iembed.eq.1.or.iembed.eq.3) then
         pxpi1=ptpi*pxqembd/ptq
         pypi1=ptpi*pyqembd/ptq
         phidecomp=acos(pxqembd/ptq)
         if(pyqembd.lt.0) phidecomp=2.*pi-phidecomp
      else
         phidecomp=2.*pi*RANART(NSEED)
         pxpi1=ptpi*cos(phidecomp)
         pypi1=ptpi*sin(phidecomp)
      endif
c     Embedded quark/antiquark are assumed to have pz=0:
      pzpi1=0.
c     Insert the two parent charged pions, 
c     ipion=1 for the pion containing the leading quark, 
c     ipion=2 for the pion containing the leading antiquark of the same flavor:
      do ipion=1,2
         if(ipion.eq.1) then
            idpi=idpi1
            pxpi=pxpi1
            pypi=pypi1
            pzpi=pzpi1
         elseif(ipion.eq.2) then
            idpi=-idpi1
            pxpi=-pxpi1
            pypi=-pypi1
            pzpi=-pzpi1
         endif
         NATT=NATT+1
         KATT(NATT,1)=idpi
         KATT(NATT,2)=40
         KATT(NATT,3)=0
         PATT(NATT,1)=pxpi
         PATT(NATT,2)=pypi
         PATT(NATT,3)=pzpi
         PATT(NATT,4)=sqrt(pxpi**2+pypi**2+pzpi**2+pichmass**2)
         EATT=EATT+PATT(NATT,4)
         GXAR(NATT)=xjet
         GYAR(NATT)=yjet
         GZAR(NATT)=0.
         FTAR(NATT)=0.
         ITYPAR(NATT)=KATT(NATT,1) 
         PXAR(NATT)=PATT(NATT,1)
         PYAR(NATT)=PATT(NATT,2)
         PZAR(NATT)=PATT(NATT,3)
         PEAR(NATT)=PATT(NATT,4)
         XMAR(NATT)=pichmass
      enddo
c
clin-8/2009
c     Randomly embed a number of soft pions around each high-Pt quark in pair:
      if(nsembd.gt.0) then
         do ipion=1,2
            do ispion=1,nsembd
               idsart=3+int(3*RANART(NSEED))
               if(idsart.eq.3) then 
                  pimass=pichmass
                  idpis=-211
               elseif(idsart.eq.4) then 
                  pimass=pi0mass
                  idpis=111
               else
                  pimass=pichmass
                  idpis=211
               endif
               NATT=NATT+1
               KATT(NATT,1)=idpis
               KATT(NATT,2)=40
               KATT(NATT,3)=0
c     theta: relative angle between soft pion & associated high-Pt q or qbar,
c     generate theta and phi uniformly:
c     Note: it is not generated uniformly in solid angle because that gives 
c     a valley at theta=0, unlike the jet-like correlation (a peak at theta=0).
               theta=tmaxembd*RANART(NSEED)
               phi=2.*pi*RANART(NSEED)
               pxspi=psembd*sin(theta)*cos(phi)
               pyspi=psembd*sin(theta)*sin(phi)
               pzspi=psembd*cos(theta)
               if(ipion.eq.1) then
                  call rotate(pxpi1,pypi1,pzpi1,pxspi,pyspi,pzspi)
               else
                  call rotate(-pxpi1,-pypi1,-pzpi1,pxspi,pyspi,pzspi)
               endif
ctest off
c               write(99,*) "2  ",pxspi,pyspi,pzspi
               PATT(NATT,1)=pxspi
               PATT(NATT,2)=pyspi
               PATT(NATT,3)=pzspi
               PATT(NATT,4)=sqrt(psembd**2+pimass**2)
               EATT=EATT+PATT(NATT,4)
               GXAR(NATT)=xjet
               GYAR(NATT)=yjet
               GZAR(NATT)=0.
               FTAR(NATT)=0.
               ITYPAR(NATT)=KATT(NATT,1) 
               PXAR(NATT)=PATT(NATT,1)
               PYAR(NATT)=PATT(NATT,2)
               PZAR(NATT)=PATT(NATT,3)
               PEAR(NATT)=PATT(NATT,4)
               XMAR(NATT)=pimass
            enddo
         enddo
      endif
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
                DO 502 iy=1, NATT
                   write(9943,*)itypar(iy),pxar(iy),pyar(iy),pzar(iy),
     &                  xmar(iy),gxar(iy),gyar(iy),gzar(iy),ftar(iy)
 502            continue
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      
clin-8/2009-end
c
      return
      end
c.................... zpc.f
c	PROGRAM ZPC
