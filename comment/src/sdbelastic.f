      subroutine sdbelastic(SRT,sdb)
      PARAMETER (srt0=2.012)
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1     px1n,py1n,pz1n,dp1n
      common /dpi/em2,lb2
      common /para8/ idpert,npertd,idxsec
      SAVE   
c
      sdb=0.
      sdbel=0.
      if(srt.le.(em1+em2)) return
      s=srt**2
c     For elastic collisions:
      if(idxsec.eq.1.or.idxsec.eq.3) then
c     1/3: assume the same |matrix element|**2/s (after averaging over initial 
c     spins and isospins) for d+Baryon elastic at the same sqrt(s);
         sdbel=fdbel(s)
      elseif(idxsec.eq.2.or.idxsec.eq.4) then
c     2/4: assume the same |matrix element|**2/s (after averaging over initial 
c     spins and isospins) for d+Baryon elastic at the same sqrt(s)-threshold:
         threshold=em1+em2
         snew=(srt-threshold+srt0)**2
         sdbel=fdbel(snew)
      endif
      sdb=sdbel
      return
      end
clin-9/2008 Deuteron+Baryon elastic collisions
