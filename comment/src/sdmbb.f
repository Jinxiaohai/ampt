      subroutine sdmbb(SRT,sdm,ianti)
      PARAMETER (AMN=0.939457,AMP=0.93828,
     1     AM0=1.232,AM1440=1.44,AM1535=1.535,srt0=2.012)
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1     px1n,py1n,pz1n,dp1n
      common /dpi/em2,lb2
      common /dpifsl/lbnn1,lbnn2,lbnd1,lbnd2,lbns1,lbns2,lbnp1,lbnp2,
     1     lbdd1,lbdd2,lbds1,lbds2,lbdp1,lbdp2,lbss1,lbss2,
     2     lbsp1,lbsp2,lbpp1,lbpp2
      common /dpifsm/xmnn1,xmnn2,xmnd1,xmnd2,xmns1,xmns2,xmnp1,xmnp2,
     1     xmdd1,xmdd2,xmds1,xmds2,xmdp1,xmdp2,xmss1,xmss2,
     2     xmsp1,xmsp2,xmpp1,xmpp2
      common /dpisig/sdmel,sdmnn,sdmnd,sdmns,sdmnp,sdmdd,sdmds,sdmdp,
     1     sdmss,sdmsp,sdmpp
      common /para8/ idpert,npertd,idxsec
      COMMON/RNDF77/NSEED
      SAVE   
c
      sdm=0.
      sdmel=0.
      sdmnn=0.
      sdmnd=0.
      sdmns=0.
      sdmnp=0.
      sdmdd=0.
      sdmds=0.
      sdmdp=0.
      sdmss=0.
      sdmsp=0.
      sdmpp=0.
ctest off check Xsec using fixed mass for resonances:
c      if(lb1.ge.25.and.lb1.le.27) then
c         em1=0.776
c      elseif(lb1.eq.28) then
c         em1=0.783
c      elseif(lb1.eq.0) then
c         em1=0.548
c      endif
c      if(lb2.ge.25.and.lb2.le.27) then
c         em2=0.776
c      elseif(lb2.eq.28) then
c         em2=0.783
c      elseif(lb2.eq.0) then
c         em2=0.548
c      endif
c
      if(srt.le.(em1+em2)) return
      s=srt**2
      pinitial=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt
      fs=fnndpi(s)
c     Determine isospin and spin factors for the ratio between 
c     Deuteron+Meson->BB and BB->Deuteron+Meson cross sections:
      if(idxsec.eq.1.or.idxsec.eq.2) then
c     Assume B+B -> d+Meson has the same cross sections as N+N -> d+pi, 
c     then determine d+Meson -> B+B cross sections:
         if((lb1.ge.3.and.lb1.le.5).or.
     1        (lb2.ge.3.and.lb2.le.5)) then
            xnnfactor=8./9.
         elseif((lb1.ge.25.and.lb1.le.27).or.
     1           (lb2.ge.25.and.lb2.le.27)) then
            xnnfactor=8./27.
         elseif(lb1.eq.28.or.lb2.eq.28) then
            xnnfactor=8./9.
         elseif(lb1.eq.0.or.lb2.eq.0) then
            xnnfactor=8./3.
         endif
      else
c     Assume d+Meson -> B+B has the same cross sections as d+pi -> N+N:
      endif
clin-9/2008 For elastic collisions:
      if(idxsec.eq.1.or.idxsec.eq.3) then
c     1/3: assume the same |matrix element|**2/s (after averaging over initial 
c     spins and isospins) for d+Meson elastic at the same sqrt(s);
         sdmel=fdpiel(s)
      elseif(idxsec.eq.2.or.idxsec.eq.4) then
c     2/4: assume the same |matrix element|**2/s (after averaging over initial 
c     spins and isospins) for d+Meson elastic at the same sqrt(s)-threshold:
         threshold=em1+em2
         snew=(srt-threshold+srt0)**2
         sdmel=fdpiel(snew)
      endif
c
*     NN: DETERMINE THE CHARGE STATES OF PARTICLESIN THE FINAL STATE
      IF(((lb1.eq.5.or.lb2.eq.5.or.lb1.eq.27.or.lb2.eq.27)
     1     .and.ianti.eq.0).or.
     2     ((lb1.eq.3.or.lb2.eq.3.or.lb1.eq.25.or.lb2.eq.25)
     3     .and.ianti.eq.1))THEN
*     (1) FOR Deuteron+(pi+,rho+) -> P+P or DeuteronBar+(pi-,rho-)-> PBar+PBar:
         lbnn1=1
         lbnn2=1
         xmnn1=amp
         xmnn2=amp
      ELSEIF(lb1.eq.3.or.lb2.eq.3.or.lb1.eq.26.or.lb2.eq.26
     1        .or.lb1.eq.28.or.lb2.eq.28.or.lb1.eq.0.or.lb2.eq.0)THEN
*     (2) FOR Deuteron+(pi0,rho0,omega,eta) -> N+P 
*     or DeuteronBar+(pi0,rho0,omega,eta) ->NBar+PBar:
         lbnn1=2
         lbnn2=1
         xmnn1=amn
         xmnn2=amp
      ELSE
*     (3) FOR Deuteron+(pi-,rho-) -> N+N or DeuteronBar+(pi+,rho+)-> NBar+NBar:
         lbnn1=2
         lbnn2=2
         xmnn1=amn
         xmnn2=amn
      ENDIF
      if(srt.gt.(xmnn1+xmnn2)) then
         pfinal=sqrt((s-(xmnn1+xmnn2)**2)*(s-(xmnn1-xmnn2)**2))/2./srt
         if(idxsec.eq.1) then
c     1: assume the same |matrix element|**2/s (after averaging over initial 
c     spins and isospins) for B+B -> deuteron+meson at the same sqrt(s);
            sdmnn=fs*pfinal/pinitial*3./16.*xnnfactor
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmnn1+xmnn2,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
c     2: assume the same |matrix element|**2/s for B+B -> deuteron+meson 
c     at the same sqrt(s)-threshold:
               sdmnn=fnndpi(snew)*pfinal/pinitial*3./16.*xnnfactor
            elseif(idxsec.eq.4) then
c     4: assume the same |matrix element|**2/s for B+B <- deuteron+meson 
c     at the same sqrt(s)-threshold:
               sdmnn=fnndpi(snew)*pfinal/pinitial/6.
            endif
         elseif(idxsec.eq.3) then
c     3: assume the same |matrix element|**2/s for B+B <- deuteron+meson 
c     at the same sqrt(s):
            sdmnn=fs*pfinal/pinitial/6.
         endif
      endif
c     
*     ND: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
      lbnd1=1+int(2*RANART(NSEED))
      lbnd2=6+int(4*RANART(NSEED))
      if(lbnd1.eq.1) then
         xmnd1=amp
      elseif(lbnd1.eq.2) then
         xmnd1=amn
      endif
      xmnd2=am0
      if(srt.gt.(xmnd1+xmnd2)) then
         pfinal=sqrt((s-(xmnd1+xmnd2)**2)*(s-(xmnd1-xmnd2)**2))/2./srt
         if(idxsec.eq.1) then
c     The spin- and isospin-averaged factor is 8-times larger for ND:
            sdmnd=fs*pfinal/pinitial*3./16.*(xnnfactor*8.)
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmnd1+xmnd2,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sdmnd=fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*8.)
            elseif(idxsec.eq.4) then
               sdmnd=fnndpi(snew)*pfinal/pinitial/6.
            endif
         elseif(idxsec.eq.3) then
            sdmnd=fs*pfinal/pinitial/6.
         endif
      endif
c
*     NS: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
      lbns1=1+int(2*RANART(NSEED))
      lbns2=10+int(2*RANART(NSEED))
      if(lbns1.eq.1) then
         xmns1=amp
      elseif(lbns1.eq.2) then
         xmns1=amn
      endif
      xmns2=am1440
      if(srt.gt.(xmns1+xmns2)) then
         pfinal=sqrt((s-(xmns1+xmns2)**2)*(s-(xmns1-xmns2)**2))/2./srt
         if(idxsec.eq.1) then
            sdmns=fs*pfinal/pinitial*3./16.*(xnnfactor*2.)
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmns1+xmns2,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sdmns=fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*2.)
            elseif(idxsec.eq.4) then
               sdmns=fnndpi(snew)*pfinal/pinitial/6.
            endif
         elseif(idxsec.eq.3) then
            sdmns=fs*pfinal/pinitial/6.
         endif
      endif
c
*     NP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
      lbnp1=1+int(2*RANART(NSEED))
      lbnp2=12+int(2*RANART(NSEED))
      if(lbnp1.eq.1) then
         xmnp1=amp
      elseif(lbnp1.eq.2) then
         xmnp1=amn
      endif
      xmnp2=am1535
      if(srt.gt.(xmnp1+xmnp2)) then
         pfinal=sqrt((s-(xmnp1+xmnp2)**2)*(s-(xmnp1-xmnp2)**2))/2./srt
         if(idxsec.eq.1) then
            sdmnp=fs*pfinal/pinitial*3./16.*(xnnfactor*2.)
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmnp1+xmnp2,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sdmnp=fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*2.)
            elseif(idxsec.eq.4) then
               sdmnp=fnndpi(snew)*pfinal/pinitial/6.
            endif
         elseif(idxsec.eq.3) then
            sdmnp=fs*pfinal/pinitial/6.
         endif
      endif
c
*     DD: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
      lbdd1=6+int(4*RANART(NSEED))
      lbdd2=6+int(4*RANART(NSEED))
      xmdd1=am0
      xmdd2=am0
      if(srt.gt.(xmdd1+xmdd2)) then
         pfinal=sqrt((s-(xmdd1+xmdd2)**2)*(s-(xmdd1-xmdd2)**2))/2./srt
         if(idxsec.eq.1) then
            sdmdd=fs*pfinal/pinitial*3./16.*(xnnfactor*16.)
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmdd1+xmdd2,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sdmdd=fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*16.)
            elseif(idxsec.eq.4) then
               sdmdd=fnndpi(snew)*pfinal/pinitial/6.
            endif
         elseif(idxsec.eq.3) then
            sdmdd=fs*pfinal/pinitial/6.
         endif
      endif
c
*     DS: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
      lbds1=6+int(4*RANART(NSEED))
      lbds2=10+int(2*RANART(NSEED))
      xmds1=am0
      xmds2=am1440
      if(srt.gt.(xmds1+xmds2)) then
         pfinal=sqrt((s-(xmds1+xmds2)**2)*(s-(xmds1-xmds2)**2))/2./srt
         if(idxsec.eq.1) then
            sdmds=fs*pfinal/pinitial*3./16.*(xnnfactor*8.)
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmds1+xmds2,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sdmds=fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*8.)
            elseif(idxsec.eq.4) then
               sdmds=fnndpi(snew)*pfinal/pinitial/6.
            endif
         elseif(idxsec.eq.3) then
            sdmds=fs*pfinal/pinitial/6.
         endif
      endif
c
*     DP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
      lbdp1=6+int(4*RANART(NSEED))
      lbdp2=12+int(2*RANART(NSEED))
      xmdp1=am0
      xmdp2=am1535
      if(srt.gt.(xmdp1+xmdp2)) then
         pfinal=sqrt((s-(xmdp1+xmdp2)**2)*(s-(xmdp1-xmdp2)**2))/2./srt
         if(idxsec.eq.1) then
            sdmdp=fs*pfinal/pinitial*3./16.*(xnnfactor*8.)
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmdp1+xmdp2,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sdmdp=fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*8.)
            elseif(idxsec.eq.4) then
               sdmdp=fnndpi(snew)*pfinal/pinitial/6.
            endif
         elseif(idxsec.eq.3) then
            sdmdp=fs*pfinal/pinitial/6.
         endif
      endif
c
*     SS: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
      lbss1=10+int(2*RANART(NSEED))
      lbss2=10+int(2*RANART(NSEED))
      xmss1=am1440
      xmss2=am1440
      if(srt.gt.(xmss1+xmss2)) then
         pfinal=sqrt((s-(xmss1+xmss2)**2)*(s-(xmss1-xmss2)**2))/2./srt
         if(idxsec.eq.1) then
            sdmss=fs*pfinal/pinitial*3./16.*xnnfactor
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmss1+xmss2,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sdmss=fnndpi(snew)*pfinal/pinitial*3./16.*xnnfactor
            elseif(idxsec.eq.4) then
               sdmss=fnndpi(snew)*pfinal/pinitial/6.
            endif
         elseif(idxsec.eq.3) then
            sdmns=fs*pfinal/pinitial/6.
         endif
      endif
c
*     SP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
      lbsp1=10+int(2*RANART(NSEED))
      lbsp2=12+int(2*RANART(NSEED))
      xmsp1=am1440
      xmsp2=am1535
      if(srt.gt.(xmsp1+xmsp2)) then
         pfinal=sqrt((s-(xmsp1+xmsp2)**2)*(s-(xmsp1-xmsp2)**2))/2./srt
         if(idxsec.eq.1) then
            sdmsp=fs*pfinal/pinitial*3./16.*(xnnfactor*2.)
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmsp1+xmsp2,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sdmsp=fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*2.)
            elseif(idxsec.eq.4) then
               sdmsp=fnndpi(snew)*pfinal/pinitial/6.
            endif
         elseif(idxsec.eq.3) then
            sdmsp=fs*pfinal/pinitial/6.
         endif
      endif
c
*     PP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
      lbpp1=12+int(2*RANART(NSEED))
      lbpp2=12+int(2*RANART(NSEED))
      xmpp1=am1535
      xmpp2=am1535
      if(srt.gt.(xmpp1+xmpp2)) then
         pfinal=sqrt((s-(xmpp1+xmpp2)**2)*(s-(xmpp1-xmpp2)**2))/2./srt
         if(idxsec.eq.1) then
            sdmpp=fs*pfinal/pinitial*3./16.*xnnfactor
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmpp1+xmpp2,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sdmpp=fnndpi(snew)*pfinal/pinitial*3./16.*xnnfactor
            elseif(idxsec.eq.4) then
               sdmpp=fnndpi(snew)*pfinal/pinitial/6.
            endif
         elseif(idxsec.eq.3) then
            sdmpp=fs*pfinal/pinitial/6.
         endif
      endif
c
      sdm=sdmel+sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp
     1     +sdmss+sdmsp+sdmpp
      if(ianti.eq.1) then
         lbnn1=-lbnn1
         lbnn2=-lbnn2
         lbnd1=-lbnd1
         lbnd2=-lbnd2
         lbns1=-lbns1
         lbns2=-lbns2
         lbnp1=-lbnp1
         lbnp2=-lbnp2
         lbdd1=-lbdd1
         lbdd2=-lbdd2
         lbds1=-lbds1
         lbds2=-lbds2
         lbdp1=-lbdp1
         lbdp2=-lbdp2
         lbss1=-lbss1
         lbss2=-lbss2
         lbsp1=-lbsp1
         lbsp2=-lbsp2
         lbpp1=-lbpp1
         lbpp2=-lbpp2
      endif
ctest off
c      write(98,100) srt,sdmnn,sdmnd,sdmns,sdmnp,sdmdd,sdmds,sdmdp,
c     1     sdmss,sdmsp,sdmpp,sdm
c 100  format(f5.2,11(1x,f5.1))
c
      return
      end
c
clin-9/2008 Deuteron+Meson ->B+B and elastic collisions
