      SUBROUTINE bbarfs(lbb1,lbb2,ei1,ei2,iblock,iseed)
*****************************************
      COMMON/ppbmas/niso(15),nstate,ppbm(15,2),thresh(15),weight(15)
cc      SAVE /ppbmas/
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
cc      SAVE /ppb1/
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
cc      SAVE /ppmm/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
c     determine which final BbarB channel occurs:
      rd=RANART(NSEED)
      wsum=0.
      do 1001 i=1,nstate
         wsum=wsum+weight(i)
         if(rd.le.(wsum/wtot)) then
            ifs=i
            ei1=ppbm(i,1)
            ei2=ppbm(i,2)
            goto 10
         endif
 1001 continue
 10   continue
c1    pbar p
      if(ifs.eq.1) then
         iblock=1801
         lbb1=-1
         lbb2=1
      elseif(ifs.eq.2) then
c2    pbar n
         if(RANART(NSEED).le.0.5) then
            iblock=18021
            lbb1=-1
            lbb2=2
c2    nbar p
         else
            iblock=18022
            lbb1=1
            lbb2=-2
         endif
c3    nbar n
      elseif(ifs.eq.3) then
         iblock=1803
         lbb1=-2
         lbb2=2
c4&5  (pbar nbar) Delta, (p n) anti-Delta
      elseif(ifs.eq.4.or.ifs.eq.5) then
         rd=RANART(NSEED)
         if(rd.le.0.5) then
c     (pbar nbar) Delta
            if(ifs.eq.4) then
               iblock=18041
               lbb1=-1
            else
               iblock=18051
               lbb1=-2
            endif
            rd2=RANART(NSEED)
            if(rd2.le.0.25) then
               lbb2=6
            elseif(rd2.le.0.5) then
               lbb2=7
            elseif(rd2.le.0.75) then
               lbb2=8
            else
               lbb2=9
            endif
         else
c     (p n) anti-Delta
            if(ifs.eq.4) then
               iblock=18042
               lbb1=1
            else
               iblock=18052
               lbb1=2
            endif
            rd2=RANART(NSEED)
            if(rd2.le.0.25) then
               lbb2=-6
            elseif(rd2.le.0.5) then
               lbb2=-7
            elseif(rd2.le.0.75) then
               lbb2=-8
            else
               lbb2=-9
            endif
         endif
c6&7  (pbar nbar) N*(1440), (p n) anti-N*(1440)
      elseif(ifs.eq.6.or.ifs.eq.7) then
         rd=RANART(NSEED)
         if(rd.le.0.5) then
c     (pbar nbar) N*(1440)
            if(ifs.eq.6) then
               iblock=18061
               lbb1=-1
            else
               iblock=18071
               lbb1=-2
            endif
            rd2=RANART(NSEED)
            if(rd2.le.0.5) then
               lbb2=10
            else
               lbb2=11
            endif
         else
c     (p n) anti-N*(1440)
            if(ifs.eq.6) then
               iblock=18062
               lbb1=1
            else
               iblock=18072
               lbb1=2
            endif
            rd2=RANART(NSEED)
            if(rd2.le.0.5) then
               lbb2=-10
            else
               lbb2=-11
            endif
         endif
c8    Delta anti-Delta
      elseif(ifs.eq.8) then
         iblock=1808
         rd1=RANART(NSEED)
         if(rd1.le.0.25) then
            lbb1=6
         elseif(rd1.le.0.5) then
            lbb1=7
         elseif(rd1.le.0.75) then
            lbb1=8
         else
            lbb1=9
         endif
         rd2=RANART(NSEED)
         if(rd2.le.0.25) then
            lbb2=-6
         elseif(rd2.le.0.5) then
            lbb2=-7
         elseif(rd2.le.0.75) then
            lbb2=-8
         else
            lbb2=-9
         endif
c9&10 (pbar nbar) N*(1535), (p n) anti-N*(1535)
      elseif(ifs.eq.9.or.ifs.eq.10) then
         rd=RANART(NSEED)
         if(rd.le.0.5) then
c     (pbar nbar) N*(1440)
            if(ifs.eq.9) then
               iblock=18091
               lbb1=-1
            else
               iblock=18101
               lbb1=-2
            endif
            rd2=RANART(NSEED)
            if(rd2.le.0.5) then
               lbb2=12
            else
               lbb2=13
            endif
         else
c     (p n) anti-N*(1535)
            if(ifs.eq.9) then
               iblock=18092
               lbb1=1
            else
               iblock=18102
               lbb1=2
            endif
            rd2=RANART(NSEED)
            if(rd2.le.0.5) then
               lbb2=-12
            else
               lbb2=-13
            endif
         endif
c11&12 anti-Delta N*, Delta anti-N*
      elseif(ifs.eq.11.or.ifs.eq.12) then
         rd=RANART(NSEED)
         if(rd.le.0.5) then
c     anti-Delta N*
            rd1=RANART(NSEED)
            if(rd1.le.0.25) then
               lbb1=-6
            elseif(rd1.le.0.5) then
               lbb1=-7
            elseif(rd1.le.0.75) then
               lbb1=-8
            else
               lbb1=-9
            endif
            if(ifs.eq.11) then
               iblock=18111
               rd2=RANART(NSEED)
               if(rd2.le.0.5) then
                  lbb2=10
               else
                  lbb2=11
               endif
            else
               iblock=18121
               rd2=RANART(NSEED)
               if(rd2.le.0.5) then
                  lbb2=12
               else
                  lbb2=13
               endif
            endif
         else
c     Delta anti-N*
            rd1=RANART(NSEED)
            if(rd1.le.0.25) then
               lbb1=6
            elseif(rd1.le.0.5) then
               lbb1=7
            elseif(rd1.le.0.75) then
               lbb1=8
            else
               lbb1=9
            endif
            if(ifs.eq.11) then
               iblock=18112
               rd2=RANART(NSEED)
               if(rd2.le.0.5) then
                  lbb2=-10
               else
                  lbb2=-11
               endif
            else
               iblock=18122
               rd2=RANART(NSEED)
               if(rd2.le.0.5) then
                  lbb2=-12
               else
                  lbb2=-13
               endif
            endif
         endif
c13   N*(1440) anti-N*(1440)
      elseif(ifs.eq.13) then
         iblock=1813
         rd1=RANART(NSEED)
         if(rd1.le.0.5) then
            lbb1=10
         else
            lbb1=11
         endif
         rd2=RANART(NSEED)
         if(rd2.le.0.5) then
            lbb2=-10
         else
            lbb2=-11
         endif
c14   anti-N*(1440) N*(1535), N*(1440) anti-N*(1535)
      elseif(ifs.eq.14) then
         rd=RANART(NSEED)
         if(rd.le.0.5) then
c     anti-N*(1440) N*(1535)
            iblock=18141
            rd1=RANART(NSEED)
            if(rd1.le.0.5) then
               lbb1=-10
            else
               lbb1=-11
            endif
            rd2=RANART(NSEED)
            if(rd2.le.0.5) then
               lbb2=12
            else
               lbb2=13
            endif
         else
c     N*(1440) anti-N*(1535)
            iblock=18142
            rd1=RANART(NSEED)
            if(rd1.le.0.5) then
               lbb1=10
            else
               lbb1=11
            endif
            rd2=RANART(NSEED)
            if(rd2.le.0.5) then
               lbb2=-12
            else
               lbb2=-13
            endif
         endif
c15   N*(1535) anti-N*(1535)
      elseif(ifs.eq.15) then
         iblock=1815
         rd1=RANART(NSEED)
         if(rd1.le.0.5) then
            lbb1=12
         else
            lbb1=13
         endif
         rd2=RANART(NSEED)
         if(rd2.le.0.5) then
            lbb2=-12
         else
            lbb2=-13
         endif
      else
      endif
      RETURN
      END
*****************************************
* for pi pi <-> rho rho cross sections
