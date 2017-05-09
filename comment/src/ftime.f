        subroutine ftime
c       this subroutine generates formation time for the particles
c       indexing ft(i)
c       input e(i)
c       output ft(i), indx(i)
        implicit double precision (a-h, o-z)
        external ftime1
        parameter (MAXPTN=400001)
        common /para1/ mul
cc      SAVE /para1/
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
cc      SAVE /para2/
        common /para4/ iftflg, ireflg, igeflg, ibstfg
cc      SAVE /para4/
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
cc      SAVE /prec1/
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
cc      SAVE /prec4/
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
cc      SAVE /ilist4/
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        common /par1/ formt
cc      SAVE /par1/
        common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
        common /rndm3/ iseedp
cc      SAVE /rndm3/
        SAVE   
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9929,*)"call ftime subroutine."
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        iseed=iseedp
clin-6/07/02 initialize here to expedite compiling, instead in zpcbdt:
        do 1001 i = 1, MAXPTN
           ct(i)=0d0
           ot(i)=0d0
 1001   continue
        tlarge=1000000.d0
clin-6/07/02-end
        if (iftflg .eq. 0) then
c     5/01/01 different prescription for parton initial formation time:
           if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
              do 1002 i = 1, mul
                 if (ft0(i) .gt. tlarge) ft0(i) = tlarge
 1002         continue
              goto 150
           else
c     5/01/01-end
           do 1003 i = 1, MAXPTN
              ft0(i) = tlarge
 1003      continue
           do 1004 i = 1, mul
              xmt2 = px0(i) ** 2 + py0(i) ** 2 + xmp ** 2
              formt = xmt2 / e0(i)           
              ft0(i) = ftime1(iseed)
              if (ft0(i) .gt. tlarge) ft0(i) = tlarge
 1004      continue
c     5/01/01:
        endif
        end if
c     5/01/01:
 150        continue
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           DO 485 IXIAOHAI = 1, MUL
              write(9929,*)ft0(i)
 485          CONTINUE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c        call index1(MAXPTN, mul, ft0, indx)
        if (mul .gt. 1) then
           call index1(MAXPTN, mul, ft0, indx)
        else
clin-7/09/03: need to set value for mul=1:
           indx(1)=1
        end if
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 483 ihaihai=1, MUL
           write(9928,*)"indx(i)",indx(ihaihai)
 483       continue
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        
c
        return
        end
