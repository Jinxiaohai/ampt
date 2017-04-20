        subroutine npik(irun,iseed,dt,nt,ictrl,i1,i2,srt,
     &                  pcx,pcy,pcz,nchrg,ratiok,iblock)
*
* Process: PI + N -> K(-) + ANYTHING
* 1.  PI- + P -> P + K0 + K-
* 2.  PI+ + N -> P + K+ + K- 
* 3.  PI0 + P -> P + K+ + K-
* 4.  PI0 + N -> P + K0 + K-
* 5.  PI0 + N -> N + K+ + K-
* 6.  PI- + P -> N + K+ + K-
* 7.  PI- + N -> N + K0 + K-
* NOTE: the mass of K is assumed to be same as K0. ie. 0.498 NOT 0.494
****************************************
      PARAMETER      (MAXSTR=150001,MAXR=1,PI=3.1415926)
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
      dimension bb(3),p1(4),p2(4),p3(4),px(4),py(4),pz(4)
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      SAVE   
        px1cm=pcx
        py1cm=pcy
        pz1cm=pcz
        ictrl = 1
        lb1=lb(i1)
        lb2=lb(i2)
        k1=i1
        k2=i2
c        k1 must be bayron. k2 be meson. If not, exchange.
        if(lb2.eq.1.or.lb2.eq.2.or.(lb2.ge.6.and.lb2.le.13)) then
            k1=i2
            k2=i1
        endif
cbz3/8/99 neutralk
cbz10/12/99
c        LB(I1) = 1 + 2 * RANART(NSEED)
c        LB(I2) = 23
        LB(k1) = 1 + int(2*RANART(NSEED))
        LB(k2) = 23
c       pkmax=sqrt((srt**2-(aka+0.938+aka)**2)*(srt**2-(aka+0.938-aka)**2))
c     &           /2./srt
        pkmax=sqrt((srt**2-(aka+0.938+aka)**2)
     &           *(srt**2-(aka+0.938-aka)**2))/2./srt
        pk = RANART(NSEED)*pkmax
c-----------------------------------------------------
        css=1.-2.*RANART(NSEED)
        sss=sqrt(1.-css**2)
        fai=2*3.1415926*RANART(NSEED)
        p3(1)=pk*sss*cos(fai)
        p3(2)=pk*sss*sin(fai)
        p3(3)=pk*css
        eip = srt - sqrt(aka**2 + pk**2)
        rmnp=sqrt(eip**2-pk**2)
        do 1001 i= 1, 3
           bb(i) = -1.*p3(i)/eip
 1001   continue
c        bb: velocity of the other two particles as a whole.
        pznp=sqrt((rmnp**2-(aka+0.938)**2)
     c  *(rmnp**2-(0.938-aka)**2))/2./rmnp    
c-----------------------------------------------------
        css=1.-2.*RANART(NSEED)
        sss=sqrt(1.-css**2)
        fai=2*3.1415926*RANART(NSEED)
        p1(1)=pznp*sss*cos(fai)
        p1(2)=pznp*sss*sin(fai)
        p1(3)=pznp*css
        p1(4)=sqrt(0.938**2+pznp**2)
        p2(4)=sqrt(aka**2+pznp**2)
        do 1002 i=1,3
           p2(i)=-1.*p1(i)
 1002   continue
c        p1,p2: the momenta of the two particles in their cms
c        p1: momentum of N or P
c        p2: momentum of anti_kaon
c        p3: momentum of K0 or K+
        ilo=1
c        write(61,*)'--------p1,p2',p1,p2
c        write(61,*)'--------bb',bb
        call lorntz(ilo,bb,p1,p2)
c******* Checking *************
c        pxsum = p1(1)+p2(1)+p3(1)
c        pysum = p1(2)+p2(2)+p3(2)
c        pzsum = p1(3)+p2(3)+p3(3)
c        pesum = p1(4)+p2(4)+sqrt(p3(1)**2+p3(2)**2+p3(3)**2+aka**2)-srt
c        write(61,*)'---p1,pxsum',p1,pxsum
c        write(61,*)'---p2,pysum',p2,pysum
c        write(61,*)'---p3,pzsum',p3,pzsum
c        write(61,*)'---pesum',pesum
c***********************************
* Rotate the momenta of particles in the cms of I1 & I2
* px(1), py(1), pz(1): momentum of I1
* px(2), py(2), pz(2): momentum of I2
* px(3), py(3), pz(3): momentum of anti-kaon
c     10/28/02 get rid of argument usage mismatch in rotate():
        pxrota=p1(1)
        pyrota=p1(2)
        pzrota=p1(3)
c        call rotate(pcx,pcy,pcz,p1(1),p1(2),p1(3))
        call rotate(pcx,pcy,pcz,pxrota,pyrota,pzrota)
        p1(1)=pxrota
        p1(2)=pyrota
        p1(3)=pzrota
c
        pxrota=p2(1)
        pyrota=p2(2)
        pzrota=p2(3)
c        call rotate(pcx,pcy,pcz,p2(1),p2(2),p2(3))
        call rotate(pcx,pcy,pcz,pxrota,pyrota,pzrota)
        p2(1)=pxrota
        p2(2)=pyrota
        p2(3)=pzrota
c
        pxrota=p3(1)
        pyrota=p3(2)
        pzrota=p3(3)
c        call rotate(pcx,pcy,pcz,p3(1),p3(2),p3(3))
        call rotate(pcx,pcy,pcz,pxrota,pyrota,pzrota)
        p3(1)=pxrota
        p3(2)=pyrota
        p3(3)=pzrota
        nnn=nnn+1
c     K(-)
        lpion(nnn,irun)=21
c     aka: rest mass of K
        epion(nnn,irun)=aka
* Find the momenta of particles in the final state in the nucleus_nucleus
* cms frame.   Lorentz transformation into lab frame.
        e1cm   = sqrt(0.938**2 + p1(1)**2 + p1(2)**2 + p1(3)**2)
        p1beta = p1(1)*betax + p1(2)*betay + p1(3)*betaz
        transf = gamma * ( gamma*p1beta / (gamma+1) + e1cm)
        pt1i1 = betax*transf + p1(1)
        pt2i1 = betay*transf + p1(2)
        pt3i1 = betaz*transf + p1(3)
        eti1  = 0.938
        lb1   = lb(k1)
* For second nulceon, same
        e2cm   = sqrt(aka**2 + p3(1)**2 + p3(2)**2 + p3(3)**2)
        p2beta = p3(1)*betax + p3(2)*betay + p3(3)*betaz
        transf = gamma * ( gamma*p2beta / (gamma+1) + e2cm)
        pt1i2 = betax*transf + p3(1)
        pt2i2 = betay*transf + p3(2)
        pt3i2 = betaz*transf + p3(3)
        eti2  = aka
        lb2   = lb(k2)
c        if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0.)then
*       k1 stand for nucleon, k2 stand for kaon. lpion stand for Kbar.
                p(1,k1)=pt1i1
                p(2,k1)=pt2i1
                p(3,k1)=pt3i1
                e(k1)=eti1
                lb(k1)=lb1
                p(1,k2)=pt1i2
                p(2,k2)=pt2i2
                p(3,k2)=pt3i2
                e(k2)=eti2
                lb(k2)=lb2
c                px1 = p(1,i1)
c                py1 = p(2,i1)
c                pz1 = p(3,i1)
c                em1 = e(i1)
c                id(i1) = 2
c                id(i2) = 2
c                id1 = id(i1)
c     K(+)K(-) production
                iblock = 101
* Get Kaons' momenta and coordinates in nucleus-nucleus cms. frame.
c  p2:  momentum of anti-kaon.
c        epcmk = sqrt(epion(nnn,irun)**2 + p2(1)**2 + p2(2)**2 + p2(3)**2)
        epcmk = sqrt(epion(nnn,irun)**2 + p2(1)**2+p2(2)**2+p2(3)**2)
        betak = p2(1)*betax + p2(2)*betay + p2(3)*betaz
        transf= gamma*(gamma*betak/(gamma+1.) + epcmk)
        ppion(1,nnn,irun)=betax*transf + p2(1)
        ppion(2,nnn,irun)=betay*transf + p2(2)
        ppion(3,nnn,irun)=betaz*transf + p2(3)
clin-5/2008:
        dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
cbz3/2/99
c        write(400,*)'2 ', ppion(1,nnn,irun), ppion(2,nnn,irun),
c     &                    ppion(3,nnn,irun), dt*nt, srt
cbz3/2/99end
c        write(420,*)ppion(1,nnn,irun), ppion(2,nnn,irun),
c     &                    ppion(3,nnn,irun), dt*nt, srt
        k=i2
        if(lb(i1).eq.1.or.lb(i1).eq.2) k=i1
        rpion(1,nnn,irun)=r(1,k)
        rpion(2,nnn,irun)=r(2,k)
        rpion(3,nnn,irun)=r(3,k)
        return
        end
c-----------------------------------------------------------------------
c.....extracted from G. Song's ART expasion including K- interactions
c.....file `PIHYPN.FOR'
******************************************
