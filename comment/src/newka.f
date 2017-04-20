        subroutine newka(icase,irun,iseed,dt,nt,ictrl,i1,i2,
     &                                   srt,pcx,pcy,pcz,iblock)
      PARAMETER      (MAXSTR=150001,MAXR=1)
      PARAMETER      (AKA=0.498)
      COMMON   /AA/  R(3,MAXSTR)
cc      SAVE /AA/
      COMMON   /BB/  P(3,MAXSTR)
cc      SAVE /BB/
      COMMON   /CC/  E(MAXSTR)
cc      SAVE /CC/
      COMMON   /EE/  ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      COMMON   /BG/BETAX,BETAY,BETAZ,GAMMA
cc      SAVE /BG/
      COMMON   /NN/NNN
cc      SAVE /NN/
      COMMON   /RUN/NUM
cc      SAVE /RUN/
      COMMON   /PA/RPION(3,MAXSTR,MAXR)
cc      SAVE /PA/
      COMMON   /PB/PPION(3,MAXSTR,MAXR)
cc      SAVE /PB/
      COMMON   /PC/EPION(MAXSTR,MAXR)
cc      SAVE /PC/
      COMMON   /PD/LPION(MAXSTR,MAXR)
cc      SAVE /PD/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
c
        logical lb1bn, lb2bn,lb1mn,lb2mn
cbz3/7/99 neutralk
c        logical lb1bn1, lb2bayon1, lb1bn0, lb2bn0
        logical lb1bn1, lb2bn1, lb1bn0, lb2bn0
cbz3/7/99 neutralk end
        logical lb1mn0, lb2mn0, lb1mn1, lb2mn1
        logical lb1mn2, lb2mn2
        icase=-1
c        icase: flag for the type of reaction that is going to happen.
c        icase=-1,  no desired reaction, return to main program.
c              1,  NN,ND,DD
c              2,  PI+N, PI+D
c              3,  K(-) absorption.
        nchrg=-100
c        nchrg: Net charges of the two incoming particles.
        ictrl = 1
        lb1=lb(i1)
        lb2=lb(i2)
        em1=e(i1)
        em2=e(i2)
        lb1bn=lb1.eq.1.or.lb1.eq.2.or.(lb1.gt.5.and.lb1.le.13)
        lb2bn=lb2.eq.1.or.lb2.eq.2.or.(lb2.gt.5.and.lb2.le.13)
        lb1bn0=lb1.eq.2.or.lb1.eq.7.or.lb1.eq.10.or.lb1.eq.12
        lb2bn0=lb2.eq.2.or.lb2.eq.7.or.lb2.eq.10.or.lb2.eq.12
        lb1bn1=lb1.eq.1.or.lb1.eq.8.or.lb1.eq.11.or.lb1.eq.13
        lb2bn1=lb2.eq.1.or.lb2.eq.8.or.lb2.eq.11.or.lb2.eq.13
        lb1mn=em1.lt.0.2.or.lb1.eq.0.or.(lb1.ge.25.and.lb1.le.29)
        lb2mn=em2.lt.0.2.or.lb2.eq.0.or.(lb2.ge.25.and.lb2.le.29)
        lb1mn0=lb1.eq.0.or.lb1.eq.4.or.lb1.eq.26.or.
     &                        lb1.eq.28.or.lb1.eq.29
        lb2mn0=lb2.eq.0.or.lb2.eq.4.or.lb2.eq.26.or.
     &                        lb2.eq.28.or.lb2.eq.29
        lb1mn1= lb1.eq.5.or.lb1.eq.27
        lb2mn1= lb2.eq.5.or.lb2.eq.27
        lb1mn2=lb1.eq.3.or.lb1.eq.25
        lb2mn2=lb2.eq.3.or.lb2.eq.25
c        1. consider N+N, N+Resonance, R + R reactions
        if(lb1bn.and.lb2bn) then
c     NN,ND,DD:
           icase=1
c     total cross section
           sig=40.
           if(lb1.eq.9.and.lb2.eq.9) then
                nchrg=4
           endif   
           if((lb1bn1.and.lb2.eq.9)
     &        .or.(lb2bn1.and.lb1.eq.9))then
                nchrg=3
           endif
           if((lb1bn0.and.lb2.eq.9)
     &        .or.(lb2bn0.and.lb1.eq.9)
     &        .or.(lb1bn1.and.lb2bn1)) then
                   nchrg=2
           endif
           if((lb1bn1.and.lb2bn0).or.(lb1.eq.6.and.lb2.eq.9)
     &        .or.(lb2bn1.and.lb1bn0)
     &        .or.(lb2.eq.6.and.lb1.eq.9))then
                   nchrg=1
           endif
           if((lb1bn0.and.lb2bn0).or.(lb1bn1.and.lb2.eq.6)
     &              .or.(lb2bn1.and.lb1.eq.6)) then
                   nchrg=0
           endif
           if((lb1bn0.and.lb2.eq.6)
     &        .or.(lb2bn0.and.lb1.eq.6))then
                nchrg=-1
           endif
           if(lb1.eq.6.and.lb2.eq.6) then
                nchrg=-2
           endif
c     brsig = x2kaon_no_isospin(srt)
           if(nchrg.ge.-1.and.nchrg.le.2) then
c     K,Kbar prduction x sect.
                   brsig = x2kaon(srt)
           else
                   brsig=0.0
c                if(nchrg.eq.-2.or.nchrg.eq.3) then
c                   brsig = x2kaon(srt+0.938-1.232)
c                else
c     nchrg=4
c                   brsig = x2kaon(srt+2.*(0.938-1.232))
c                endif
           endif
cbz3/7/99 neutralk
           BRSIG = 2.0 * BRSIG
cbz3/7/99 neutralk end
        endif
c        2. consider PI(meson:eta,omega,rho,phi) + N(N*,D)
        if((lb1bn.and.lb2mn).OR.(lb2bn.and.lb1mn)) then
c     PN,PD
          icase=2
          sig=20.
          sigma0 = piNsg0(srt)
          brsig=0.0
          if((lb1bn1.and.lb2mn0)
     &       .or.(lb2bn1.and.lb1mn0).
     & or.(lb1bn0.and.lb2mn1).or.(lb2bn0.and.lb1mn1).
     & or.(lb1.eq.9.and.lb2mn2).or.(lb2.eq.9.and.lb1mn2))then
                nchrg=1
cbz3/2/99/song
c                if(lb1bn1.or.lb2bn1) brsig=2.0*sigma0
c                if(lb1bn0.or.lb2bn0) brsig=0.5*sigma0
                if(lb1bn1.or.lb2bn1) brsig=0.5*sigma0
                if(lb1bn0.or.lb2bn0) brsig=2.0*sigma0
cbz3/2/99/song end
c                if(lb1.eq.9.or.lb2.eq.9) brsig=1.5*sigma0
          endif
          if( (lb1bn0.and.lb2mn0 )
     &       .or.(lb2bn0.and.lb1mn0)
     &  .or.(lb1bn1.and.lb2mn2).or.(lb2bn1.and.lb1mn2)
     &  .or.(lb1.eq.6.and.lb2mn1).or.(lb2.eq.6.and.lb1mn1)) then
                nchrg=0
                if(lb1bn1.or.lb2bn1) then
cbz3/2/99/song
c                  brsig=1.5*sigma0
                  brsig=3.0*sigma0
cbz3/2/99/song end
cbz3/11/99/song
c                  ratiok = 1./3.
                  ratiok = 2./3.
cbz3/11/99/song end
c                  ratiok: the ratio of channels: ->nK+k- vs. -> pK0K-
                endif
                if(lb1bn0.or.lb2bn0) then
                  brsig=2.5*sigma0
cbz3/2/99/song
c                  ratiok = 0.8
                  ratiok = 0.2
cbz3/2/99/song end
                endif
c                if(lb1.eq.6.or.lb2.eq.6) then
c     lb=6 : D-
c                  brsig=1.5*sigma0
c                  ratiok = 0.5
c                endif
          endif
          if( (lb1bn0.and.lb2mn2)
     &       .or.(lb2bn0.and.lb1mn2)
     & .or.(lb1.eq.6.and.lb2mn0).or.(lb2.eq.6.and.lb1mn0)) then
                nchrg=-1
                if(lb1bn0.or.lb2bn0) brsig=sigma0
c                if(lb1.eq.6.or.lb2.eq.6) brsig=sigma0
          endif
c          if((lb1.eq.6.and.lb2mn2).or.(lb2.eq.6.and.lb1mn2))then
c                nchrg=-2
c          endif
c          if((lb1bn1.and.lb2mn1).or.(lb2bn1.and.lb1mn1)
c    &           .or.(lb1.eq.9.and.lb2mn0).or.(lb2.eq.9.and.lb1mn0)) then
c                nchrg=2
c          endif
cbz3/11/99 neutralk
          if((lb1.eq.6.and.lb2mn2)
     &       .or.(lb2.eq.6.and.lb1mn2))then
                nchrg=-2
          endif
cbz3/11/99 neutralk
cbz3/8/99 neutralk
          if((lb1bn1.and.lb2mn1)
     &       .or.(lb2bn1.and.lb1mn1)
     & .or.(lb1.eq.9.and.lb2mn0).or.(lb2.eq.9.and.lb1mn0)) then
                nchrg=2
          endif
cbz3/8/99 neutralk end
cbz3/7/99 neutralk
          IF (NCHRG .GE. -2 .AND. NCHRG .LE. 2) THEN
             BRSIG = 3.0 * SIGMA0
          END IF
cbz3/7/99 neutralk end
        endif
c        3. consider K- + N(N*,D) absorption.
c        if((lb1bn.and.lb2.eq.21).OR.(lb2bn.and.lb1.eq.21)) then
        if( (lb1bn.and.(lb2.eq.21.or.lb2.eq.-30)).OR.
     &     (lb2bn.and.(lb1.eq.21.or.lb1.eq.-30)) )then 
c          bmass=em1+em2-aka
          bmass=0.938
          if(srt.le.(bmass+aka)) then
cbz3/2/99
c                write(100,*)'--lb1,lb2,em1,em2,srt',lb1,lb2,em1,em2,srt
cbz3/2/99end
                pkaon=0.
          else
            pkaon=sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
          endif
          sig=0.
          if(lb1.eq.1.or.lb2.eq.1.or.lb1.eq.8.or.lb2.eq.8.or.
     &    lb1.eq.11.or.lb2.eq.11.or.lb1.eq.13.or.lb2.eq.13) then
c          K- + (D+,N*+)p ->
              nchrg=0
              sigela=akPel(pkaon)
              sigsgm=3.*akPsgm(pkaon)
              sig=sigela+sigsgm+akPlam(pkaon)
          endif
          if(lb1.eq.2.or.lb2.eq.2.or.lb1.eq.7.or.lb2.eq.7.or.
     &    lb1.eq.10.or.lb2.eq.10.or.lb1.eq.12.or.lb2.eq.12) then
c          K- + (D0, N*0)n ->
              nchrg=-1
              sigela=akNel(pkaon)
              sigsgm=2.*akNsgm(pkaon)
              sig=sigela+sigsgm+akNlam(pkaon)
          endif
          if(lb1.eq.6.or.lb2.eq.6) then
c     K- + D-
              nchrg=-2
              sigela=akNel(pkaon)
              sigsgm=akNsgm(pkaon)
              sig=sigela+sigsgm
          endif
          if(lb1.eq.9.or.lb2.eq.9) then
c     K- + D++
              nchrg=1
              sigela=akPel(pkaon)
              sigsgm=2.*akPsgm(pkaon)
              sig=sigela+sigsgm+akPlam(pkaon)
          endif
cbz3/8/99 neutralk
          sigela = 0.5 * (AKPEL(PKAON) + AKNEL(PKAON))
          SIGSGM = 1.5 * AKPSGM(PKAON) + AKNSGM(PKAON)
          SIG = sigela + SIGSGM + AKPLAM(PKAON)
cbz3/8/99 neutralk end
          if(sig.gt.1.e-7) then
c     K(-) + N reactions
              icase=3
              brel=sigela/sig
              brsgm=sigsgm/sig
c              branch_lambda=akNlam(pkaon)/sig
              brsig = sig
          endif
        endif
c        4. meson + hyperon -> K- + N
c        if(((lb1.ge.14.and.lb1.le.17).and.lb2mn).OR.
c     &     ((lb2.ge.14.and.lb2.le.17).and.lb1mn)) then
        if(((lb1.ge.14.and.lb1.le.17).and.(lb2.ge.3.and.lb2.le.5)).OR.
     &     ((lb2.ge.14.and.lb2.le.17).and.(lb1.ge.3.and.lb1.le.5)))then
c        first classify the reactions due to total charge.
           nchrg=-100
           if((lb1.eq.15.and.(lb2.eq.3.or.lb2.eq.25)).OR.
     &              (lb2.eq.15.and.(lb1.eq.3.or.lb1.eq.25))) then
                nchrg=-2
c     D-
                  bmass=1.232
           endif
           if((lb1.eq.15.and.lb2mn0).or.(lb2.eq.15.and.lb1mn0).OR.
     &       ((lb1.eq.14.or.lb1.eq.16).and.(lb2.eq.3.or.lb2.eq.25)).OR.
     &       ((lb2.eq.14.or.lb2.eq.16).and.(lb1.eq.3.or.lb1.eq.25)))then
                nchrg=-1
c     n
                 bmass=0.938
           endif
           if((lb1.eq.15.and.(lb2.eq.5.or.lb2.eq.27)).OR.
     &              (lb2.eq.15.and.(lb1.eq.5.or.lb1.eq.27)).or.
     &        (lb1.eq.17.and.(lb2.eq.3.or.lb2.eq.25)).OR.
     &              (lb2.eq.17.and.(lb1.eq.3.or.lb1.eq.25)).or.
     &       ((lb1.eq.14.or.lb1.eq.16).and.lb2mn0).OR.
     &       ((lb2.eq.14.or.lb2.eq.16).and.lb1mn0)) then
                nchrg=0
c     p
                 bmass=0.938
           endif
           if((lb1.eq.17.and.lb2mn0).or.(lb2.eq.17.and.lb1mn0).OR.
     &       ((lb1.eq.14.or.lb1.eq.16).and.(lb2.eq.5.or.lb2.eq.27)).OR.
     &       ((lb2.eq.14.or.lb2.eq.16).and.(lb1.eq.5.or.lb1.eq.27)))then
                nchrg=1
c     D++
                 bmass=1.232
           endif
           sig = 0.
           if(nchrg.ne.-100.and.srt.gt.(aka+bmass)) then
c     PI+sigma or PI + Lambda => Kbar + N reactions
             icase=4
c             pkaon=sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
             pkaon=sqrt(((srt**2-(aka**2+0.938**2))/2./0.938)**2-aka**2)
c     lambda + Pi
             if(lb1.eq.14.or.lb2.eq.14) then
                if(nchrg.ge.0) sigma0=akPlam(pkaon)
                if(nchrg.lt.0) sigma0=akNlam(pkaon)
c     sigma + pi
             else
c     K-p or K-D++
                if(nchrg.ge.0) sigma0=akPsgm(pkaon)
c     K-n or K-D-
                if(nchrg.lt.0) sigma0=akNsgm(pkaon)
cbz3/8/99 neutralk
                SIGMA0 = 1.5 * AKPSGM(PKAON) + AKNSGM(PKAON)
cbz3/8/99 neutralk end
             endif
             sig=(srt**2-(aka+bmass)**2)*(srt**2-(aka-bmass)**2)/
     &         (srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)*sigma0
cbz3/8/99 neutralk
c     if(nchrg.eq.-2.or.nchrg.eq.1) sig=2.*sig K-D++, K-D-
c     K0barD++, K-D-
             if(nchrg.eq.-2.or.nchrg.eq.2) sig=2.*sig
cbz3/8/99 neutralk end
c             the factor 2 comes from spin of delta, which is 3/2
c             detailed balance. copy from Page 423 of N.P. A614 1997
cbz3/8/99 neutralk
             IF (LB1 .EQ. 14 .OR. LB2 .EQ. 14) THEN
                SIG = 4.0 / 3.0 * SIG
             ELSE IF (NCHRG .EQ. -2 .OR. NCHRG .EQ. 2) THEN
                SIG = 8.0 / 9.0 * SIG
             ELSE
                SIG = 4.0 / 9.0 * SIG
             END IF
cbz3/8/99 neutralk end
             brsig = sig
             if(sig.lt.1.e-7) sig = 1.e-7
           endif
csp05/07/01
* comment icase=4 statement below if only inelastic
c     PI+L/Si => Kbar + N  OR ELASTIC SCATTERING
           icase=4
           brsig = sig
c     elastic xsecn of 10mb
           sigela = 10.
           sig = sig + sigela
           brel = sigela/sig
cc          brsig = sig
csp05/07/01 end   
        endif
c
c        if(em2.lt.0.2.and.em1.lt.0.2) then
c     PI + PI 
c             icase=5
c     assumed PI PI total x section.
c              sig=50.
c     Mk + Mkbar
c              s0=aka+aka
c              brsig = 0.
c              if(srt.gt.s0) brsig = 2.7*(1.-s0**2/srt**2)**0.76
c              x section for PIPI->KKbar   PRC43 (1991) 1881
c        endif
        if(icase.eq.-1) then
           ictrl = -1
           return
        endif
        px1cm=pcx
        py1cm=pcy
        pz1cm=pcz
        ds=sqrt(sig/31.4)
        dsr=ds+0.1
        ec=(em1+em2+0.02)**2
c        ec=3.59709
c        if((e(i1).ge.1.).and.(e(i2).ge.1.)) ec = 4.75
        call distce(i1,i2,dsr,ds,dt,ec,srt,ic,px1cm,py1cm,pz1cm)
        if(ic.eq.-1) then
c     no anti-kaon production
           ictrl = -1
c           in=in+1
c           write(60,*)'--------------distance-----',in
           return
        endif
clin-10/24/02 set to 0: ik,ik0-3,il,im,im3-4,in,inpion,ipipi, 
c     sgsum,sgsum1,sgsum3:
        ik=0
        ik0=0
        ik1=0
        ik2=0
        ik3=0
        il=0
        im=0
        im3=0
        im4=0
        in=0
        inpion=0
        ipipi=0
        sgsum=0.
        sgsum1=0.
        sgsum3=0.
        if(icase.eq.1) then
           ik=ik+1
           if(srt.gt.2.8639) then
                ik0=ik0+1
                if(em1.lt.1.0.and.em2.lt.1.0) then
                        ik1=ik1+1
                        sgsum1=sgsum1+brsig
c                        ratio_1=sgsum1/ik1/40.
                endif
                if(em1.gt.1.0.and.em2.gt.1.0) then
                        ik3=ik3+1
                        sgsum3=sgsum3+brsig
c                        ratio_3=sgsum3/ik3/40.
                endif
                if(em1.gt.1.0.and.em2.lt.1.0) ik2=ik2+1
                if(em1.lt.1.0.and.em2.gt.1.0) ik2=ik2+1
                sgsum=sgsum+brsig
c                ratio=sgsum/ik0/40.
           endif
        endif
        if(icase.eq.2) inpion=inpion+1
        if(icase.eq.5) ipipi=ipipi+1
c        write(62,*)'ik1,ik2,ik3',ik1,ik2,ik3,ratio_1,ratio_3,ratio
c        write(62,*)'inpion,ipipi',inpion,ipipi
        if(RANART(NSEED).gt.(brsig/sig)) then
c     no anti-kaon production
           ictrl = -1
           return
        endif
        il=il+1
c        kaons could be created now.
        if(icase.eq.1) then
          in=in+1
c          write(60,*)'------in,s2kaon,sig=',in,brsig,sig,lb1,lb2
          call nnkaon(irun,iseed,
     &          ictrl,i1,i2,iblock,srt,pcx,pcy,pcz,nchrg)
        endif
        if(icase.eq.2) then
          im=im+1
c          call npik(irun,iseed,dt,nt,ictrl,i1,i2,srt,
c     &              pcx,pcy,pcz,nchrg,ratiok)
          call npik(irun,iseed,dt,nt,ictrl,i1,i2,srt,
     &              pcx,pcy,pcz,nchrg,ratiok,iblock)
        endif
c
        if(icase.eq.3) then
          im3=im3+1
c          write(63,*)'im3,lb1,lb2,pkaon',im3,lb1,lb2,pkaon
c          write(63,*)'sig,el,sigma',sig,brel,brsgm
c          write(63,*)'srt,pcx,pcy,pcz,em1,em2',srt,pcx,pcy,pcz,em1,em2
          call kaonN(brel,brsgm,irun,iseed,dt,nt,ictrl,
     &                i1,i2,iblock,srt,pcx,pcy,pcz,nchrg)
c         this subroutine format is diff. since three final states are possible
        endif
c
        if(icase.eq.4) then
          im4=im4+1
c          write(64,*)'im4,sigma0,branch,sig=',im4,sigma0,brsig,sig
c          write(64,*)'lb1,lb2,em1,em2,pkaon=',lb1,lb2,em1,em2,pkaon
csp06/07/01
      if(RANART(NSEED).lt.brel) then
         ielstc = 1
      else
         ielstc = 0
      endif                  
c          call Pihypn(ielstc,irun,iseed,dt,nt,ictrl,i1,i2,srt,
c     &                   pcx,pcy,pcz,nchrg)
          call Pihypn(ielstc,irun,iseed,dt,nt,ictrl,i1,i2,srt,
     &                   pcx,pcy,pcz,nchrg,iblock)
csp06/07/01 end
        endif
c        if(icase.eq.5) then
c          im5=im5+1
c          write(65,*)'---im5,s2kaon,sig=',im5,brsig,sig
c          call pipikaon(irun,iseed,dt,nt,ictrl,i1,i2,srt,pcx,pcy,pcz)
c        endif
cbz3/2/99
c        write(101,*)lb1,lb2,lb(i1),lb(i2)
c        write(101,*)em1,em2,e(i1),e(i2),srt
cbz3/2/99end
        return
        end
******************************************
* for pp-->pp + kaon + anti-kaon
c      real*4 function X2kaon(srt)
