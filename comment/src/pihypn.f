        subroutine pihypn(ielstc,irun,iseed,dt,nt,ictrl,i1,i2,
     &     srt,pcx,pcy,pcz,nchrg,iblock)
*
* Process: PI + sigma(or Lambda) -> Kbar + N
* NOTE: the mass of K is assumed to be same as K0. ie. 0.498 NOT 0.494
******************************************
c NOTE: for PI + Hyperon: the produced kaons have mass 0.498
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
      dimension p1(4),p2(4)
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
        px1cm=pcx
        py1cm=pcy
        pz1cm=pcz
        ictrl = 1
csp06/07/01
        if(ielstc .eq. 1) then
*    L/Si + meson -> L/Si + meson
             k1=i1
             k2=i2
           dm3=e(k1)
           dm4=e(k2)
           iblock = 10
        else
           iblock = 12
csp06/07/01 end  
c        PI + Sigma(or Lambda) -> Kbar + N
        k1=i1
        k2=i2
c        k1 must be bayron! So if I1 is PI, exchange k1 & k2.
        if(lb(i1).lt.14.or.lb(i1).gt.17) then
           k1=i2
           k2=i1
        endif
cbz3/8/99 neutralk
        LB(K1) = 1 + int(2*RANART(NSEED))
        if(nchrg.eq.-2) lb(k1)=6
c     if(nchrg.eq.-1) lb(k1)=2
c     if(nchrg.eq. 0) lb(k1)=1
c     if(nchrg.eq. 1) lb(k1)=9
        IF (NCHRG .EQ. 2) LB(K1) = 9
cbz3/8/99 neutralk end
c     K-
        lb(k2)=21
        dm3=0.938
        if(nchrg.eq.-2.or.nchrg.eq.1) dm3=1.232
        dm4=aka
c        dm3,dm4: the mass of final state particles.
         endif
********Now, antikaon will be created.
c        call antikaon_fstate(iseed,srt,dm1,dm2,dm3,dm4,px,py,pz,icou1)
c        pkmax: the maximum momentum of anti-kaon
        pkmax=sqrt((srt**2-(dm3+dm4)**2)*(srt**2-(dm3-dm4)**2))
     &         /2./srt
        pk=pkmax
c-----------------------------------------------------
        css=1.-2.*RANART(NSEED)
        sss=sqrt(1.-css**2)
        fai=2*3.1415926*RANART(NSEED)
        p1(1)=pk*sss*cos(fai)
        p1(2)=pk*sss*sin(fai)
        p1(3)=pk*css
        do 1001 i=1,3
           p2(i)=-1.*p1(i)
 1001   continue
c        p1,p2: the momenta of the two particles in their cms
c        p1: momentum of kaon
c        p2: momentum of Kbar
* Rotate the momenta of particles in the cms of I1 & I2
clin-10/28/02 get rid of argument usage mismatch in rotate():
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
clin-10/28/02-end
* Find the momenta of particles in the final state in the nucleus_nucleus
* cms frame.   Lorentz transformation into lab frame.
        e1cm   = sqrt(dm3**2 + p1(1)**2 + p1(2)**2 + p1(3)**2)
        p1beta = p1(1)*betax + p1(2)*betay + p1(3)*betaz
        transf = gamma * ( gamma*p1beta / (gamma+1) + e1cm)
        pt1i1 = betax*transf + p1(1)
        pt2i1 = betay*transf + p1(2)
        pt3i1 = betaz*transf + p1(3)
        eti1  = dm3
        lb1   = lb(k1)
* For second kaon, same
        e2cm   = sqrt(dm4**2 + p2(1)**2 + p2(2)**2 + p2(3)**2)
        p2beta = p2(1)*betax + p2(2)*betay + p2(3)*betaz
        transf = gamma * ( gamma*p2beta / (gamma+1) + e2cm)
        pt1i2 = betax*transf + p2(1)
        pt2i2 = betay*transf + p2(2)
        pt3i2 = betaz*transf + p2(3)
cbz3/2/99
c        write(400,*)'3 ', pt1i2, pt2i2, pt3i2, dt*nt, srt
cbz3/2/99end
c        write(430,*)pt1i2, pt2i2, pt3i2, dt*nt, srt
        eti2  = dm4
        lb2   = lb(k2)
c        if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0.)then
c        k1=i1
c        k2=i2
*       k1 stand for nucleon, k2 stand for kaon.
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
cc                iblock = 101  ! K(+)K(-) production
* Get Kaons' momenta and coordinates in nucleus-nucleus cms. frame.
        return
        end
c-----------------------------------------------------------------------
c.....extracted from G. Song's ART expasion including K- interactions
c.....file `KAONN.FOR'
****************************************
