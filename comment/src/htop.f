      SUBROUTINE HTOP
c
      PARAMETER (MAXSTR=150001)
      PARAMETER (MAXPTN=400001)
      PARAMETER (MAXIDL=4001)
      DOUBLE PRECISION  GX0, GY0, GZ0, FT0, PX0, PY0, PZ0, E0, XMASS0
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS, ptwo, xmdq, ptwox, ptwoy, ptwoz
      dimension it(4)
      COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
cc      SAVE /HMAIN2/
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
cc      SAVE /HMAIN1/
      COMMON /PARA1/ MUL
cc      SAVE /PARA1/
      COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &     PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &     XMASS0(MAXPTN), ITYP0(MAXPTN)
cc      SAVE /prec1/
      COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
cc      SAVE /ilist7/
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
cc      SAVE /ARPRC/
      common /decom/ptwo(2,5)
cc      SAVE /decom/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  不进入ZPC的粒子。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      COMMON /NOPREC/ NNOZPC, ITYPN(MAXIDL),
     &     GXN(MAXIDL), GYN(MAXIDL), GZN(MAXIDL), FTN(MAXIDL),
     &     PXN(MAXIDL), PYN(MAXIDL), PZN(MAXIDL), EEN(MAXIDL),
     &     XMN(MAXIDL)
cc      SAVE /NOPREC/
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
c     7/20/01: use double precision
c     otherwise sometimes beta>1 and gamma diverge in lorenz():
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
cc      SAVE /SOFT/
      common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
clin-8/2015:
      DOUBLE PRECISION vxp0,vyp0,vzp0,xstrg0,ystrg0,xstrg,ystrg
      common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
     1     xstrg0(MAXPTN),ystrg0(MAXPTN),
     2     xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)
c      DOUBLE PRECISION  vxp0,vyp0,vzp0
c      common /precpa/ vxp0(MAXPTN), vyp0(MAXPTN), vzp0(MAXPTN)
cc      SAVE /precpa/
      common /para7/ ioscar,nsmbbbar,nsmmeson
      COMMON /AREVT/ IAEVT, IARUN, MISS
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
      SAVE
c
        npar=0
        nnozpc=0
clin-5b/2008 calculate the number of hadrons to be converted to q/qbar:
        if((isoft.eq.4.or.isoft.eq.5).and.(ioscar.eq.2.or.ioscar.eq.3))
     1       then
           nsmbbbar=0
           nsmmeson=0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     遍历所有的生成的部分子(包括胶子，光子，和中微子)注意ID号。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           write(9936,*)"NATT = ", natt
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<

           do i=1,natt
              id=ITYPAR(i)
              idabs=iabs(id)
              i2=MOD(idabs/10,10)
clin-9/2011 determine spectator nucleons consistently
c              if(PXAR(i).eq.0.and.PYAR(i).eq.0.and.PEAR(i)
c     1             .ge.(HINT1(1)/2*0.99).and.
c     2             .and.(id.eq.2112.or.id.eq.2212)) then
              if(abs(PXAR(i)).le.epsiPt.and.abs(PYAR(i)).le.epsiPt
     1             .and.(PZAR(i).gt.amax1(0.,PZPROJ-epsiPz)
     2                .or.PZAR(i).lt.(-PZTARG+epsiPz))
     3             .and.(id.eq.2112.or.id.eq.2212)) then
c     spectator proj or targ nucleons without interactions, do not enter ZPC:
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$       重子
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
              elseif(idabs.gt.1000.and.i2.ne.0) then
c     baryons to be converted to q/qbar:
                 nsmbbbar=nsmbbbar+1
              elseif((idabs.gt.100.and.idabs.lt.1000)
     1                .or.idabs.gt.10000) then
c     mesons to be converted to q/qbar:
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$       介子
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
                 nsmmeson=nsmmeson+1
              endif
           enddo
clin-6/2009:
           if(ioscar.eq.2.or.ioscar.eq.3) then
              write(92,*) iaevt,miss,3*nsmbbbar+2*nsmmeson,
     1             nsmbbbar,nsmmeson,natt,natt-nsmbbbar-nsmmeson
           endif
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           write(9977,*) iaevt, 3*nsmbbbar+2*nsmmeson
           write(9977,*) 'event#, total # of initial partons after
     &     string melting'
           write(9977,*) 'String melting converts ',nsmbbbar,
     &      ' baryons ', nsmmeson, 'mesons'
           write(9977,*) 'Total # of initial particles= ',natt
           write(9977,*) 'Total # of initial particles
     &    (gamma,e,muon,...) not entering ZPC= ',natt-nsmbbbar-nsmmeson
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c           write(92,*) iaevt, 3*nsmbbbar+2*nsmmeson
c           write(92,*) ' event#, total # of initial partons after string
c     1 melting'
c           write(92,*) 'String melting converts ',nsmbbbar, ' baryons &'
c     1, nsmmeson, ' mesons'
c           write(92,*) 'Total # of initial particles= ',natt
c           write(92,*) 'Total # of initial particles (gamma,e,muon,...)
c     1 not entering ZPC= ',natt-nsmbbbar-nsmmeson
        endif


c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      开始新的循环和条件
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
clin-5b/2008-over
        do 100 i=1,natt
           id=ITYPAR(i)
           idabs=iabs(id)
           i4=MOD(idabs/1000,10)
           i3=MOD(idabs/100,10)
           i2=MOD(idabs/10,10)
           i1=MOD(idabs,10)
           rnum=RANART(NSEED)
           ftime=0.197*PEAR(i)/(PXAR(i)**2+PYAR(i)**2+XMAR(i)**2)
           inozpc=0
           it(1)=0
           it(2)=0
           it(3)=0
           it(4)=0
c
clin-9/2011 determine spectator nucleons consistently
c           if(PXAR(i).eq.0.and.PYAR(i).eq.0.and.PEAR(i)
c     1 .ge.(HINT1(1)/2*0.99).and.((id.eq.2112).or.(id.eq.2212))) then
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    各种生成的类型的核子的判断，以及是否进入ZPC给出判定。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
              if(abs(PXAR(i)).le.epsiPt.and.abs(PYAR(i)).le.epsiPt
     1             .and.(PZAR(i).gt.amax1(0.,PZPROJ-epsiPz)
     2                .or.PZAR(i).lt.(-PZTARG+epsiPz))
     3             .and.(id.eq.2112.or.id.eq.2212)) then
c     spectator proj or targ nucleons without interactions, do not enter ZPC:
              inozpc=1
           elseif(idabs.gt.1000.and.i2.ne.0) then
c     baryons:
              if(((i4.eq.1.or.i4.eq.2).and.i4.eq.i3)
     1 .or.(i4.eq.3.and.i3.eq.3)) then
                 if(i1.eq.2) then
                    if(rnum.le.(1./2.)) then
                       it(1)=i4
                       it(2)=i3*1000+i2*100+1
                    elseif(rnum.le.(2./3.)) then
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    else
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    endif
                 elseif(i1.eq.4) then
                    if(rnum.le.(2./3.)) then
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    else
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    endif
                 endif
              elseif(i4.eq.1.or.i4.eq.2) then
                 if(i1.eq.2) then
                    if(rnum.le.(1./2.)) then
                       it(1)=i2
                       it(2)=i4*1000+i3*100+1
                    elseif(rnum.le.(2./3.)) then
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    else
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    endif
                 elseif(i1.eq.4) then
                    if(rnum.le.(2./3.)) then
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    else
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    endif
                 endif
              elseif(i4.ge.3) then
                 it(1)=i4
                 if(i3.lt.i2) then
                    it(2)=i2*1000+i3*100+1
                 else
                    it(2)=i3*1000+i2*100+3
                 endif
              endif
c       antibaryons:
              if(id.lt.0) then
                 it(1)=-it(1)
                 it(2)=-it(2)
              endif
c     isoft=4or5 decompose diquark flavor it(2) to two quarks it(3)&(4):

              if(isoft.eq.4.or.isoft.eq.5) then
                 it(3)=MOD(it(2)/1000,10)
                 it(4)=MOD(it(2)/100,10)
              endif
           elseif((idabs.gt.100.and.idabs.lt.1000)
     1 .or.idabs.gt.10000) then
c     mesons:
              if(i3.eq.i2) then
                 if(i3.eq.1.or.i3.eq.2) then
                    if(rnum.le.0.5) then
                       it(1)=1
                       it(2)=-1
                    else
                       it(1)=2
                       it(2)=-2
                    endif
                 else
                    it(1)=i3
                    it(2)=-i3
                 endif
              else
                 if((isign(1,id)*(-1)**i3).eq.1) then
                    it(1)=i3
                    it(2)=-i2
                 else
                    it(1)=i2
                    it(2)=-i3
                 endif
              endif
           else
c     save other particles (leptons and photons) outside of ZPC:
              inozpc=1
           endif
c






           if(inozpc.eq.1) then
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  非强子的处理。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
              NJSGS(i)=0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  不进入zpc的粒子的个数的统计。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
              nnozpc=nnozpc+1
              itypn(nnozpc)=ITYPAR(i)
              pxn(nnozpc)=PXAR(i)
              pyn(nnozpc)=PYAR(i)
              pzn(nnozpc)=PZAR(i)
              een(nnozpc)=PEAR(i)
              xmn(nnozpc)=XMAR(i)
              gxn(nnozpc)=GXAR(i)
              gyn(nnozpc)=GYAR(i)
              gzn(nnozpc)=GZAR(i)
              ftn(nnozpc)=FTAR(i)
           else
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    对于介子的处理。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
              NJSGS(i)=2
              ptwo(1,5)=dble(ulmass(it(1)))
              ptwo(2,5)=dble(ulmass(it(2)))
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  更新的ptwo数组的数值。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
              call decomp(dble(patt(i,1)),dble(patt(i,2)),
     1             dble(patt(i,3)),dble(XMAR(i)),i,it(1))
              ipamax=2
              if((isoft.eq.4.or.isoft.eq.5)
     1 .and.iabs(it(2)).gt.1000) ipamax=1
              do 1001 ipar=1,ipamax
                 npar=npar+1
                 ityp0(npar)=it(ipar)
                 px0(npar)=ptwo(ipar,1)
                 py0(npar)=ptwo(ipar,2)
                 pz0(npar)=ptwo(ipar,3)
                 e0(npar)=ptwo(ipar,4)
                 xmass0(npar)=ptwo(ipar,5)
                 gx0(npar)=dble(GXAR(i))
                 gy0(npar)=dble(GYAR(i))
                 gz0(npar)=dble(GZAR(i))
                 ft0(npar)=dble(ftime)
                 lstrg0(npar)=i
                 lpart0(npar)=ipar
                 vxp0(npar)=dble(patt(i,1)/patt(i,4))
                 vyp0(npar)=dble(patt(i,2)/patt(i,4))
                 vzp0(npar)=dble(patt(i,3)/patt(i,4))
clin-8/2015: set parent string information for this parton:
                 xstrg(npar)=xstrg0(i)
                 ystrg(npar)=ystrg0(i)
                 istrg(npar)=istrg0(i)
 1001     continue
 200      format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
 201      format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
c
              if((isoft.eq.4.or.isoft.eq.5)
     1 .and.iabs(it(2)).gt.1000) then
                 NJSGS(i)=3
                 xmdq=ptwo(2,5)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  重子的处理。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
                 ptwo(1,5)=dble(ulmass(it(3)))
                 ptwo(2,5)=dble(ulmass(it(4)))
c     8/19/02 avoid actual argument in common blocks of DECOMP:
c                 call decomp(ptwo(2,1),ptwo(2,2),ptwo(2,3),xmdq)
                 ptwox=ptwo(2,1)
                 ptwoy=ptwo(2,2)
                 ptwoz=ptwo(2,3)
                 call decomp(ptwox,ptwoy,ptwoz,xmdq,i,it(1))
c
                 do 1002 ipar=1,2
                    npar=npar+1
                    ityp0(npar)=it(ipar+2)
                    px0(npar)=ptwo(ipar,1)
                    py0(npar)=ptwo(ipar,2)
                    pz0(npar)=ptwo(ipar,3)
                    e0(npar)=ptwo(ipar,4)
                    xmass0(npar)=ptwo(ipar,5)
                    gx0(npar)=dble(GXAR(i))
                    gy0(npar)=dble(GYAR(i))
                    gz0(npar)=dble(GZAR(i))
                    ft0(npar)=dble(ftime)
                    lstrg0(npar)=i
                    lpart0(npar)=ipar+1
                    vxp0(npar)=dble(patt(i,1)/patt(i,4))
                    vyp0(npar)=dble(patt(i,2)/patt(i,4))
                    vzp0(npar)=dble(patt(i,3)/patt(i,4))
clin-8/2015: set parent string information for this parton:
                    xstrg(npar)=xstrg0(i)
                    ystrg(npar)=ystrg0(i)
                    istrg(npar)=istrg0(i)
 1002        continue
              endif
c
           endif
 100        continue
      MUL=NPAR




c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
      write(9935,*)"id, px, py, pz, e, mass, x, y, z, t,
     & lstrg, lpart, vx, vy, vz, xstrg, ystrg, istrg"
      DO 489 i = 1, NPAR
         write(9935,487)ityp0(i), px0(i), py0(i), pz0(i), e0(i),
     &        xmass0(i), gx0(i), gy0(i), gz0(i), ft0(i),
     &        lstrg0(i), lpart0(i), vxp0(i), vyp0(i), vzp0(i),
     &        xstrg(i), ystrg(i), istrg(i)
 489  continue

      DO 488 i = 1, nnozpc
         write(9934,486)itypn(i), pxn(i), pyn(i), pzn(i), een(i),
     &        xmn(i), gxn(i), gyn(i), gzn(i), ftn(i)
 488  continue

 487  format (i6, 2x,9(f12.4, 2x), 2x,i6, 2x,i6, 2x,5(f12.4, 2x)
     &     2x,i10)
 486  format (i10, 2x,9(f12.4, 2x))

c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<




c
clin-5b/2008:
      if((isoft.eq.4.or.isoft.eq.5).and.(ioscar.eq.2.or.ioscar.eq.3))
     1     then
         if((natt-nsmbbbar-nsmmeson).ne.nnozpc)
     1        write(92,*) 'Problem with the total # of initial particles
     2 (gamma,e,muon,...) not entering ZPC'
         if((3*nsmbbbar+2*nsmmeson).ne.npar)
     1        write(92,*) 'Problem with the total # of initial partons
     2 after string melting'
      endif
c
      RETURN
      END
c=======================================================================
