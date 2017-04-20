c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    icall是传进去的常数,其数值为2和数值3.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      SUBROUTINE locldr(icall,drlocl)
c
      implicit double precision (a-h, o-z)
      dimension ftp0(3),pxp0(3),pyp0(3),pzp0(3),pep0(3)
      common /loclco/gxp(3),gyp(3),gzp(3),ftp(3),
     1     pxp(3),pyp(3),pzp(3),pep(3),pmp(3)
cc      SAVE /loclco/
      common /prtn23/ gxp0(3),gyp0(3),gzp0(3),ft0fom
cc      SAVE /prtn23/
      common /lor/ enenew, pxnew, pynew, pznew
cc      SAVE /lor/
      SAVE   
c     for 2-body kinematics:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    对于两体的动力学
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      if(icall.eq.2) then
         etot=pep(1)+pep(2)
         bex=(pxp(1)+pxp(2))/etot
         bey=(pyp(1)+pyp(2))/etot
         bez=(pzp(1)+pzp(2))/etot
c     boost the reference frame down by beta to get to the pair rest frame:
         do 1001 j=1,2
            beta2 = bex ** 2 + bey ** 2 + bez ** 2
            gam = 1.d0 / dsqrt(1.d0 - beta2)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    对于合成的粒子的beta过大的处理.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if(beta2.ge.0.9999999999999d0) then
               write(6,*) '4',pxp(1),pxp(2),pyp(1),pyp(2),
     1              pzp(1),pzp(2),pep(1),pep(2),pmp(1),pmp(2),
     2          dsqrt(pxp(1)**2+pyp(1)**2+pzp(1)**2+pmp(1)**2)/pep(1),
     3          dsqrt(pxp(1)**2+pyp(1)**2+pzp(1)**2)/pep(1)
               write(6,*) '4a',pxp(1)+pxp(2),pyp(1)+pyp(2),
     1              pzp(1)+pzp(2),etot
               write(6,*) '4b',bex,bey,bez,beta2,gam
            endif
c
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    对两个粒子的两个"四动量"进行lorentz变换
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
            call lorenz(ftp(j),gxp(j),gyp(j),gzp(j),bex,bey,bez)
            gxp0(j)=pxnew
            gyp0(j)=pynew
            gzp0(j)=pznew
            ftp0(j)=enenew
            call lorenz(pep(j),pxp(j),pyp(j),pzp(j),bex,bey,bez)
            pxp0(j)=pxnew
            pyp0(j)=pynew
            pzp0(j)=pznew
            pep0(j)=enenew
 1001    continue
c     
         if(ftp0(1).ge.ftp0(2)) then
            ilate=1
            iearly=2
         else
            ilate=2
            iearly=1
         endif
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     转换到相同的时间进行处理
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
         ft0fom=ftp0(ilate)
c     
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    对较早的那个粒子进行坐标的"运动"
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
         dt0=ftp0(ilate)-ftp0(iearly)
         gxp0(iearly)=gxp0(iearly)+pxp0(iearly)/pep0(iearly)*dt0
         gyp0(iearly)=gyp0(iearly)+pyp0(iearly)/pep0(iearly)*dt0
         gzp0(iearly)=gzp0(iearly)+pzp0(iearly)/pep0(iearly)*dt0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    计算两个粒子之间的距离
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
         drlocl=dsqrt((gxp0(ilate)-gxp0(iearly))**2
     1        +(gyp0(ilate)-gyp0(iearly))**2
     2        +(gzp0(ilate)-gzp0(iearly))**2)
c     for 3-body kinematics, used for baryons formation:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    下面是对baryon的处理,但是未啥不计算三个粒子之间的距离
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      elseif(icall.eq.3) then
         etot=pep(1)+pep(2)+pep(3)
         bex=(pxp(1)+pxp(2)+pxp(3))/etot
         bey=(pyp(1)+pyp(2)+pyp(3))/etot
         bez=(pzp(1)+pzp(2)+pzp(3))/etot
         beta2 = bex ** 2 + bey ** 2 + bez ** 2
         gam = 1.d0 / dsqrt(1.d0 - beta2)
         if(beta2.ge.0.9999999999999d0) then
            write(6,*) '5',bex,bey,bez,beta2,gam
         endif
c     boost the reference frame down by beta to get to the 3-parton rest frame:
         do 1002 j=1,3
            call lorenz(ftp(j),gxp(j),gyp(j),gzp(j),bex,bey,bez)
            gxp0(j)=pxnew
            gyp0(j)=pynew
            gzp0(j)=pznew
            ftp0(j)=enenew
            call lorenz(pep(j),pxp(j),pyp(j),pzp(j),bex,bey,bez)
            pxp0(j)=pxnew
            pyp0(j)=pynew
            pzp0(j)=pznew
            pep0(j)=enenew
 1002    continue
c     
         if(ftp0(1).gt.ftp0(2)) then
            ilate=1
            if(ftp0(3).gt.ftp0(1)) ilate=3
         else
            ilate=2
            if(ftp0(3).ge.ftp0(2)) ilate=3
         endif
         ft0fom=ftp0(ilate)
c     
         if(ilate.eq.1) then
            imin=2
            imax=3
            istep=1
         elseif(ilate.eq.2) then
            imin=1
            imax=3
            istep=2
         elseif(ilate.eq.3) then
            imin=1
            imax=2
            istep=1
         endif
c     
         do 1003 iearly=imin,imax,istep
            dt0=ftp0(ilate)-ftp0(iearly)
            gxp0(iearly)=gxp0(iearly)+pxp0(iearly)/pep0(iearly)*dt0
            gyp0(iearly)=gyp0(iearly)+pyp0(iearly)/pep0(iearly)*dt0
            gzp0(iearly)=gzp0(iearly)+pzp0(iearly)/pep0(iearly)*dt0
 1003    continue
      endif
c
      RETURN
      END
c=======================================================================
