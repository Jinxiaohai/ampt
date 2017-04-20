      subroutine sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
      PARAMETER (xmd=1.8756,AP1=0.13496,AP2=0.13957,
     1     xmrho=0.770,xmomega=0.782,xmeta=0.548,srt0=2.012)
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1     px1n,py1n,pz1n,dp1n
      common /dpi/em2,lb2
      common /para8/ idpert,npertd,idxsec
      COMMON/RNDF77/NSEED
      SAVE   
c
      sdprod=0.
      sbbdpi=0.
      sbbdrho=0.
      sbbdomega=0.
      sbbdeta=0.
      if(srt.le.(em1+em2)) return
c
      ilb1=iabs(lb1)
      ilb2=iabs(lb2)
ctest off check Xsec using fixed mass for resonances:
c      if(ilb1.ge.6.and.ilb1.le.9) then
c         em1=1.232
c      elseif(ilb1.ge.10.and.ilb1.le.11) then
c         em1=1.44
c      elseif(ilb1.ge.12.and.ilb1.le.13) then
c         em1=1.535
c      endif
c      if(ilb2.ge.6.and.ilb2.le.9) then
c         em2=1.232
c      elseif(ilb2.ge.10.and.ilb2.le.11) then
c         em2=1.44
c      elseif(ilb2.ge.12.and.ilb2.le.13) then
c         em2=1.535
c      endif
c
      s=srt**2
clin-9/2012: check argument in sqrt():
      scheck=(s-(em1+em2)**2)*(s-(em1-em2)**2)
      if(scheck.le.0) then
         write(99,*) 'scheck50: ', scheck
         stop
      endif
      pinitial=sqrt(scheck)/2./srt
c      pinitial=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt
      fs=fnndpi(s)
c     Determine isospin and spin factors for the ratio between 
c     BB->Deuteron+Meson and Deuteron+Meson->BB cross sections:
      if(idxsec.eq.1.or.idxsec.eq.2) then
c     Assume B+B -> d+Meson has the same cross sections as N+N -> d+pi:
      else
c     Assume d+Meson -> B+B has the same cross sections as d+pi -> N+N, 
c     then determine B+B -> d+Meson cross sections:
         if(ilb1.ge.1.and.ilb1.le.2.and.
     1        ilb2.ge.1.and.ilb2.le.2) then
            pifactor=9./8.
         elseif((ilb1.ge.1.and.ilb1.le.2.and.
     1           ilb2.ge.6.and.ilb2.le.9).or.
     2           (ilb2.ge.1.and.ilb2.le.2.and.
     1           ilb1.ge.6.and.ilb1.le.9)) then
            pifactor=9./64.
         elseif((ilb1.ge.1.and.ilb1.le.2.and.
     1           ilb2.ge.10.and.ilb2.le.13).or.
     2           (ilb2.ge.1.and.ilb2.le.2.and.
     1           ilb1.ge.10.and.ilb1.le.13)) then
            pifactor=9./16.
         elseif(ilb1.ge.6.and.ilb1.le.9.and.
     1           ilb2.ge.6.and.ilb2.le.9) then
            pifactor=9./128.
         elseif((ilb1.ge.6.and.ilb1.le.9.and.
     1           ilb2.ge.10.and.ilb2.le.13).or.
     2           (ilb2.ge.6.and.ilb2.le.9.and.
     1           ilb1.ge.10.and.ilb1.le.13)) then
            pifactor=9./64.
         elseif((ilb1.ge.10.and.ilb1.le.11.and.
     1           ilb2.ge.10.and.ilb2.le.11).or.
     2           (ilb2.ge.12.and.ilb2.le.13.and.
     1           ilb1.ge.12.and.ilb1.le.13)) then
            pifactor=9./8.
         elseif((ilb1.ge.10.and.ilb1.le.11.and.
     1           ilb2.ge.12.and.ilb2.le.13).or.
     2           (ilb2.ge.10.and.ilb2.le.11.and.
     1           ilb1.ge.12.and.ilb1.le.13)) then
            pifactor=9./16.
         endif
      endif
c     d pi: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
*     (1) FOR P+P->Deuteron+pi+:
      IF((ilb1*ilb2).EQ.1)THEN
         lbm=5
         if(ianti.eq.1) lbm=3
         xmm=ap2
*     (2)FOR N+N->Deuteron+pi-:
      ELSEIF(ilb1.EQ.2.AND.ilb2.EQ.2)THEN
         lbm=3
         if(ianti.eq.1) lbm=5
         xmm=ap2
*     (3)FOR N+P->Deuteron+pi0:
      ELSEIF((ilb1*ilb2).EQ.2)THEN
         lbm=4
         xmm=ap1
      ELSE
c     For baryon resonances, use isospin-averaged cross sections:
         lbm=3+int(3 * RANART(NSEED))
         if(lbm.eq.4) then
            xmm=ap1
         else
            xmm=ap2
         endif
      ENDIF
c
      if(srt.ge.(xmd+xmm)) then
         pfinal=sqrt((s-(xmd+xmm)**2)*(s-(xmd-xmm)**2))/2./srt
         if((ilb1.eq.1.and.ilb2.eq.1).or.
     1        (ilb1.eq.2.and.ilb2.eq.2)) then
c     for pp or nn initial states:
            sbbdpi=fs*pfinal/pinitial/4.
         elseif((ilb1.eq.1.and.ilb2.eq.2).or.
     1           (ilb1.eq.2.and.ilb2.eq.1)) then
c     factor of 1/2 for pn or np initial states:
            sbbdpi=fs*pfinal/pinitial/4./2.
         else
c     for other BB initial states (spin- and isospin averaged):
            if(idxsec.eq.1) then
c     1: assume the same |matrix element|**2/s (after averaging over initial 
c     spins and isospins) for B+B -> deuteron+meson at the same sqrt(s);
               sbbdpi=fs*pfinal/pinitial*3./16.
            elseif(idxsec.eq.2.or.idxsec.eq.4) then
               threshold=amax1(xmd+xmm,em1+em2)
               snew=(srt-threshold+srt0)**2
               if(idxsec.eq.2) then
c     2: assume the same |matrix element|**2/s for B+B -> deuteron+meson 
c     at the same sqrt(s)-threshold:
                  sbbdpi=fnndpi(snew)*pfinal/pinitial*3./16.
               elseif(idxsec.eq.4) then
c     4: assume the same |matrix element|**2/s for B+B <- deuteron+meson 
c     at the same sqrt(s)-threshold:
                  sbbdpi=fnndpi(snew)*pfinal/pinitial/6.*pifactor
               endif
            elseif(idxsec.eq.3) then
c     3: assume the same |matrix element|**2/s for B+B <- deuteron+meson 
c     at the same sqrt(s):
               sbbdpi=fs*pfinal/pinitial/6.*pifactor
            endif
c
         endif
      endif
c     
*     d rho: DETERMINE THE CROSS SECTION TO THIS FINAL STATE:
      if(srt.gt.(xmd+xmrho)) then
         pfinal=sqrt((s-(xmd+xmrho)**2)*(s-(xmd-xmrho)**2))/2./srt
         if(idxsec.eq.1) then
            sbbdrho=fs*pfinal/pinitial*3./16.
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmd+xmrho,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sbbdrho=fnndpi(snew)*pfinal/pinitial*3./16.
            elseif(idxsec.eq.4) then
c     The spin- and isospin-averaged factor is 3-times larger for rho:
               sbbdrho=fnndpi(snew)*pfinal/pinitial/6.*(pifactor*3.)
            endif
         elseif(idxsec.eq.3) then
            sbbdrho=fs*pfinal/pinitial/6.*(pifactor*3.)
         endif
      endif
c
*     d omega: DETERMINE THE CROSS SECTION TO THIS FINAL STATE:
      if(srt.gt.(xmd+xmomega)) then
         pfinal=sqrt((s-(xmd+xmomega)**2)*(s-(xmd-xmomega)**2))/2./srt
         if(idxsec.eq.1) then
            sbbdomega=fs*pfinal/pinitial*3./16.
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmd+xmomega,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sbbdomega=fnndpi(snew)*pfinal/pinitial*3./16.
            elseif(idxsec.eq.4) then
               sbbdomega=fnndpi(snew)*pfinal/pinitial/6.*pifactor
            endif
         elseif(idxsec.eq.3) then
            sbbdomega=fs*pfinal/pinitial/6.*pifactor
         endif
      endif
c
*     d eta: DETERMINE THE CROSS SECTION TO THIS FINAL STATE:
      if(srt.gt.(xmd+xmeta)) then
         pfinal=sqrt((s-(xmd+xmeta)**2)*(s-(xmd-xmeta)**2))/2./srt
         if(idxsec.eq.1) then
            sbbdeta=fs*pfinal/pinitial*3./16.
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmd+xmeta,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sbbdeta=fnndpi(snew)*pfinal/pinitial*3./16.
            elseif(idxsec.eq.4) then
               sbbdeta=fnndpi(snew)*pfinal/pinitial/6.*(pifactor/3.)
            endif
         elseif(idxsec.eq.3) then
            sbbdeta=fs*pfinal/pinitial/6.*(pifactor/3.)
         endif
      endif
c
      sdprod=sbbdpi+sbbdrho+sbbdomega+sbbdeta
ctest off
c      write(99,111) srt,sbbdpi,sbbdrho,sbbdomega,sbbdeta,sdprod
c 111  format(6(f8.2,1x))
c
      if(sdprod.le.0) return
c
c     choose final state and assign masses here:
      x1=RANART(NSEED)
      if(x1.le.sbbdpi/sdprod) then
c     use the above-determined lbm and xmm.
      elseif(x1.le.(sbbdpi+sbbdrho)/sdprod) then
         lbm=25+int(3*RANART(NSEED))
         xmm=xmrho
      elseif(x1.le.(sbbdpi+sbbdrho+sbbdomega)/sdprod) then
         lbm=28
         xmm=xmomega
      else
         lbm=0
         xmm=xmeta
      endif
c
      return
      end
c
c     Generate angular distribution of Deuteron in the CMS frame:
