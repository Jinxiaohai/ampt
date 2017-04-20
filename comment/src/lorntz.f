        subroutine lorntz(ilo,b,pi,pj)
c       It uses to perform Lorentz (or inverse Lorentz) transformation
        dimension pi(4),pj(4),b(3)
      SAVE   
c       dimension db(3)
        bb=b(1)*b(1)+b(2)*b(2)+b(3)*b(3)
        deno3=sqrt(1.-bb)
        if(deno3.eq.0.)deno3=1.e-10
        gam=1./deno3
        ga=gam*gam/(gam+1.)
        if(ilo.eq.1) goto 100
c       Lorentz transformation
        pib=pi(1)*b(1)+pi(2)*b(2)+pi(3)*b(3)
        pjb=pj(1)*b(1)+pj(2)*b(2)+pj(3)*b(3)
c       drb=drd(1)*b(1)+drd(2)*b(2)+drd(3)*b(3)
c       drdb=db(1)*b(1)+db(2)*b(2)+db(3)*b(3)
        do 1001 i=1,3
           pi(i)=pi(i)+b(i)*(ga*pib-gam*pi(4))
           pj(i)=pj(i)+b(i)*(ga*pjb-gam*pj(4))
c       drd(i)=drd(i)+b(i)*ga*drb
c       db(i)=db(i)+b(i)*ga*drdb
 1001   continue
        pi(4)=gam*(pi(4)-pib)
        pj(4)=gam*(pj(4)-pjb)
        return
100     continue
c       inverse Lorentz transformation
        pib=pi(1)*b(1)+pi(2)*b(2)+pi(3)*b(3)
        pjb=pj(1)*b(1)+pj(2)*b(2)+pj(3)*b(3)
        do 1002 i=1,3
           pi(i)=pi(i)+b(i)*(ga*pib+gam*pi(4))
           pj(i)=pj(i)+b(i)*(ga*pjb+gam*pj(4))
 1002   continue
        pi(4)=gam*(pi(4)+pib)
        pj(4)=gam*(pj(4)+pjb)
        return
        end
