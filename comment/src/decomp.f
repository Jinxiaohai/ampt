        SUBROUTINE decomp(px0,py0,pz0,xm0,i,itq1)
c
        IMPLICIT DOUBLE PRECISION(D)  
        DOUBLE PRECISION  enenew, pxnew, pynew, pznew
clin-8/2015 changed ptwo(2,5) and related variables to double precision
c     to avoid IEEE_DIVIDE_BY_ZERO or IEEE_INVALID or IEEE_OVERFLOW_FLAG:
        DOUBLE PRECISION  de0, beta2, gam, ptwo, px0, py0, pz0, xm0
        common /lor/ enenew, pxnew, pynew, pznew
cc      SAVE /lor/
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
        common /decom/ptwo(2,5)
cc      SAVE /decom/
        COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
        COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
        common/embed/iembed,nsembd,pxqembd,pyqembd,xembd,yembd,
     1       psembd,tmaxembd,phidecomp
        SAVE   
c
        dcth=dble(RANART(NSEED))*2.d0-1.d0
        dPHI=dble(RANART(NSEED)*HIPR1(40))*2.d0
clin-6/2009 Added if embedding a high-Pt quark pair after string melting:
        if(iembed.ge.1.and.iembed.le.4) then
c     Decompose the parent high-Pt pion to q and qbar with an internal momentum
c     parallel to the pion direction so that one parton has ~the same hight Pt
c     and the other parton has a very soft Pt:
c     Note: htop() decomposes a meson to q as it(1) followed by qbar as it(2):
           if(i.eq.(natt-2*nsembd).or.i.eq.(natt-2*nsembd-1)) then
              dcth=0.d0
              dphi=dble(phidecomp)
           endif
        endif
c
        ds=xm0**2
        dpcm=dsqrt((ds-(ptwo(1,5)+ptwo(2,5))**2)
     1 *(ds-(ptwo(1,5)-ptwo(2,5))**2)/ds/4d0)
        dpz=dpcm*dcth
        dpx=dpcm*dsqrt(1.d0-dcth**2)*dcos(dphi)
        dpy=dpcm*dsqrt(1.d0-dcth**2)*dsin(dphi)
        de1=dsqrt(ptwo(1,5)**2+dpcm**2)
        de2=dsqrt(ptwo(2,5)**2+dpcm**2)
c
      de0=dsqrt(px0**2+py0**2+pz0**2+xm0**2)
        dbex=px0/de0
        dbey=py0/de0
        dbez=pz0/de0
c     boost the reference frame up by beta (pznew=gam(pz+beta e)):
      beta2 = dbex ** 2 + dbey ** 2 + dbez ** 2
      gam = 1.d0 / dsqrt(1.d0 - beta2)
      if(beta2.ge.0.9999999999999d0) then
         write(6,*) '1',dbex,dbey,dbez,beta2,gam
      endif
c
      call lorenz(de1,dpx,dpy,dpz,-dbex,-dbey,-dbez)
        ptwo(1,1)=pxnew
        ptwo(1,2)=pynew
        ptwo(1,3)=pznew
        ptwo(1,4)=enenew
      call lorenz(de2,-dpx,-dpy,-dpz,-dbex,-dbey,-dbez)
        ptwo(2,1)=pxnew
        ptwo(2,2)=pynew
        ptwo(2,3)=pznew
        ptwo(2,4)=enenew
c
      RETURN
      END
c=======================================================================
