        subroutine pihypn(ielstc,irun,iseed,dt,nt,ictrl,i1,i2,
     &     srt,pcx,pcy,pcz,nchrg,iblock)
      PARAMETER      (MAXSTR=150001,MAXR=1,PI=3.1415926)
      PARAMETER      (AKA=0.498)
      COMMON   /AA/  R(3,MAXSTR)
      COMMON   /BB/  P(3,MAXSTR)
      COMMON   /CC/  E(MAXSTR)
      COMMON   /EE/  ID(MAXSTR),LB(MAXSTR)
      COMMON   /BG/BETAX,BETAY,BETAZ,GAMMA
      COMMON   /NN/NNN
      COMMON   /RUN/NUM
      COMMON   /PA/RPION(3,MAXSTR,MAXR)
      COMMON   /PB/PPION(3,MAXSTR,MAXR)
      COMMON   /PC/EPION(MAXSTR,MAXR)
      COMMON   /PD/LPION(MAXSTR,MAXR)
      dimension p1(4),p2(4)
      COMMON/RNDF77/NSEED
      SAVE   
        px1cm=pcx
        py1cm=pcy
        pz1cm=pcz
        ictrl = 1
        if(ielstc .eq. 1) then
             k1=i1
             k2=i2
           dm3=e(k1)
           dm4=e(k2)
           iblock = 10
        else
           iblock = 12
        k1=i1
        k2=i2
        if(lb(i1).lt.14.or.lb(i1).gt.17) then
           k1=i2
           k2=i1
        endif
        LB(K1) = 1 + int(2*RANART(NSEED))
        if(nchrg.eq.-2) lb(k1)=6
        IF (NCHRG .EQ. 2) LB(K1) = 9
        lb(k2)=21
        dm3=0.938
        if(nchrg.eq.-2.or.nchrg.eq.1) dm3=1.232
        dm4=aka
         endif
        pkmax=sqrt((srt**2-(dm3+dm4)**2)*(srt**2-(dm3-dm4)**2))
     &         /2./srt
        pk=pkmax
        css=1.-2.*RANART(NSEED)
        sss=sqrt(1.-css**2)
        fai=2*3.1415926*RANART(NSEED)
        p1(1)=pk*sss*cos(fai)
        p1(2)=pk*sss*sin(fai)
        p1(3)=pk*css
        do 1001 i=1,3
           p2(i)=-1.*p1(i)
 1001   continue
        pxrota=p1(1)
        pyrota=p1(2)
        pzrota=p1(3)
        call rotate(pcx,pcy,pcz,pxrota,pyrota,pzrota)
        p1(1)=pxrota
        p1(2)=pyrota
        p1(3)=pzrota
        pxrota=p2(1)
        pyrota=p2(2)
        pzrota=p2(3)
        call rotate(pcx,pcy,pcz,pxrota,pyrota,pzrota)
        p2(1)=pxrota
        p2(2)=pyrota
        p2(3)=pzrota
        e1cm   = sqrt(dm3**2 + p1(1)**2 + p1(2)**2 + p1(3)**2)
        p1beta = p1(1)*betax + p1(2)*betay + p1(3)*betaz
        transf = gamma * ( gamma*p1beta / (gamma+1) + e1cm)
        pt1i1 = betax*transf + p1(1)
        pt2i1 = betay*transf + p1(2)
        pt3i1 = betaz*transf + p1(3)
        eti1  = dm3
        lb1   = lb(k1)
        e2cm   = sqrt(dm4**2 + p2(1)**2 + p2(2)**2 + p2(3)**2)
        p2beta = p2(1)*betax + p2(2)*betay + p2(3)*betaz
        transf = gamma * ( gamma*p2beta / (gamma+1) + e2cm)
        pt1i2 = betax*transf + p2(1)
        pt2i2 = betay*transf + p2(2)
        pt3i2 = betaz*transf + p2(3)
        eti2  = dm4
        lb2   = lb(k2)
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
        return
        end
