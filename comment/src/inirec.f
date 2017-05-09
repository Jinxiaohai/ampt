c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      得到parton的信息。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        subroutine inirec
        implicit double precision (a-h, o-z)
        external ran1
        parameter (MAXPTN=400001)
        common /para1/ mul
cc      SAVE /para1/
        common /para4/ iftflg, ireflg, igeflg, ibstfg
cc      SAVE /para4/
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
cc      SAVE /prec1/
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        common /prec3/gxs(MAXPTN),gys(MAXPTN),gzs(MAXPTN),fts(MAXPTN),
     &     pxs(MAXPTN), pys(MAXPTN), pzs(MAXPTN), es(MAXPTN),
     &     xmasss(MAXPTN), ityps(MAXPTN)
cc      SAVE /prec3/
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
cc      SAVE /prec4/
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
cc      SAVE /prec5/
        common /prec6/ etas(MAXPTN), raps(MAXPTN), taus(MAXPTN)
cc      SAVE /prec6/
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
cc      SAVE /ilist4/
cbz1/25/99
        COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
cc      SAVE /ilist7/
        COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)
cc      SAVE /ilist8/
cbz1/25/99end
        COMMON /smearz/smearp,smearh
cc      SAVE /smearz/
clin-8/2015:
c        dimension vxp(MAXPTN), vyp(MAXPTN), vzp(MAXPTN)
        common /precpb/vxp(MAXPTN),vyp(MAXPTN),vzp(MAXPTN)
clin-8/2015:
        common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
     1       xstrg0(MAXPTN),ystrg0(MAXPTN),
     2       xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)
c        common /precpa/ vxp0(MAXPTN), vyp0(MAXPTN), vzp0(MAXPTN)
cc      SAVE /precpa/
        common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
clin-6/06/02 local parton freezeout:
        common /frzprc/ 
     &       gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &       pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &       xmfrz(MAXPTN), 
     &       tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
cc      SAVE /frzprc/
        common /rndm3/ iseedp
cc      SAVE /rndm3/
        common /para7/ ioscar,nsmbbbar,nsmmeson
        COMMON /AREVT/ IAEVT, IARUN, MISS
        common /xiaohai/xiaohaiflag
        SAVE   
        iseed=iseedp
clin-6/06/02 local freezeout initialization:
c$$$        if(isoft.eq.5) then
c$$$           itlast=0
c$$$           call inifrz
c$$$        endif
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           do 519 i = 1, mul
              write(9982,*)"indx = ", indx(i)
 519          continue
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$   1001循环重新的规整粒子的索引号问题
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do 1001 i = 1, mul
clin-7/09/01 define indx(i) to save time:
c           ityp(i) = ityp0(indx(i))
c           gx(i) = gx0(indx(i))
c           gy(i) = gy0(indx(i))
c           gz(i) = gz0(indx(i))
c           ft(i) = ft0(indx(i))
c           px(i) = px0(indx(i))
c           py(i) = py0(indx(i))
c           pz(i) = pz0(indx(i))
c           e(i) = e0(indx(i))
c           xmass(i) = xmass0(indx(i))
ccbz1/25/99
c           LSTRG1(I) = LSTRG0(INDX(I))
c           LPART1(I) = LPART0(INDX(I))
ccbz1/25/99end
           indxi=indx(i)
           ityp(i) = ityp0(indxi)
           gx(i) = gx0(indxi)
           gy(i) = gy0(indxi)
           gz(i) = gz0(indxi)
           ft(i) = ft0(indxi)
           px(i) = px0(indxi)
           py(i) = py0(indxi)
           pz(i) = pz0(indxi)
           e(i) = e0(indxi)
           xmass(i) = xmass0(indxi)
           LSTRG1(I) = LSTRG0(INDXI)
           LPART1(I) = LPART0(INDXI)
           vxp(I) = vxp0(INDXI)
           vyp(I) = vyp0(INDXI)
           vzp(I) = vzp0(INDXI)
clin-8/2015:
           xstrg0(I) = xstrg(INDXI)
           ystrg0(I) = ystrg(INDXI)
           istrg0(I) = istrg(INDXI)
clin-7/09/01-end
c
clin-6/06/02 local freezeout initialization:
c$$$         if(isoft.eq.5) then
c$$$            idfrz(i)=ityp(i)
c$$$            gxfrz(i)=gx(i)
c$$$            gyfrz(i)=gy(i)
c$$$            gzfrz(i)=gz(i)
c$$$            ftfrz(i)=ft(i)
c$$$            pxfrz(i)=px(i)
c$$$            pyfrz(i)=py(i)
c$$$            pzfrz(i)=pz(i)
c$$$            efrz(i)=e(i)
c$$$            xmfrz(i)=xmass(i)
c$$$            ifrz(i)=0
c$$$         endif
clin-6/06/02-end
 1001   continue
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9927,*)"indx, ityp, x, y, z, t, px, py, pz, e, mass
     &, vx, vy, vz"
        DO 482 ihai=1, MUL
           write(9927,481)indx(ihai),ityp(ihai), gx(ihai), gy(ihai),
     &          gz(ihai), ft(ihai), px(ihai), py(ihai), pz(ihai),
     &          e(ihai), xmass(ihai),
     &          vxp(ihai), vyp(ihai), vzp(ihai)
 482    continue
 481    format (2x,i8, 2x,i8, 12(2x,f8.4))
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        
c       save particle info for fixed time analysis
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  储存粒子的信息为了固定时间的分析，也就是说gxs(i)是存诸的信息，
c$$$      gxs中的s是save的意思。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do 1002 i = 1, mul
           ityps(i) = ityp(i)
           gxs(i) = gx(i)
           gys(i) = gy(i)
           gzs(i) = gz(i)
           fts(i) = ft(i)
           pxs(i) = px(i)
           pys(i) = py(i)
           pzs(i) = pz(i)
           es(i) = e(i)
           xmasss(i) = xmass(i)
 1002   continue
clin-6/2009
        if(isoft.eq.1.and.(ioscar.eq.2.or.ioscar.eq.3))
     1       write(92,*) iaevt,miss,mul
        do 1003 i = 1, mul
           energy = e(i)
           vx(i) = px(i) / energy
           vy(i) = py(i) / energy
           vz(i) = pz(i) / energy
           if (iftflg .eq. 0) then
              formt = ft(i)
c     7/09/01 propagate partons with parent velocity till formation
c     so that partons in same hadron have 0 distance:
c            gx(i) = gx(i) + vx(i) * formt
c            gy(i) = gy(i) + vy(i) * formt
c            gz(i) = gz(i) + vz(i) * formt
            if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
               gx(i) = gx(i) + vxp(i) * formt
               gy(i) = gy(i) + vyp(i) * formt
               gz(i) = gz(i) + vzp(i) * formt
            else
               gx(i) = gx(i) + vx(i) * formt
               gy(i) = gy(i) + vy(i) * formt
               gz(i) = gz(i) + vz(i) * formt
            endif
c     7/09/01-end
c
c     3/27/00-ctest off no smear z on partons to avoid eta overflow:
c              gz(i) = gz(i)+smearp*(2d0 * ran1(iseed) - 1d0)
c     to give eta=y +- smearp*random:
c              smeary=smearp*(2d0 * ran1(iseed) - 1d0)
c              smearf=dexp(2*smeary)*(1+vz(i))/(1-vz(i)+1.d-8)
c              gz(i) = gz(i)+formt*(smearf-1)/(smearf+1)
c     3/27/00-end
           end if
clin-6/2009 write out initial parton information after string melting
c     and after propagating to its format time:
           if(ioscar.eq.2.or.ioscar.eq.3) then
              if(dmax1(abs(gx(i)),abs(gy(i)),
     1             abs(gz(i)),abs(ft(i))).lt.9999) then
clin-8/2015:
                 write(92,200) ityp(i),px(i),py(i),pz(i),xmass(i),
     1           gx(i),gy(i),gz(i),ft(i),istrg0(i),xstrg0(i),ystrg0(i)
              else
clin-8/2015:
                 write(92,201) ityp(i),px(i),py(i),pz(i),xmass(i),
     1           gx(i),gy(i),gz(i),ft(i),istrg0(i),xstrg0(i),ystrg0(i)
              endif
           endif
clin-8/2015:
c 200       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
c 201       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
c     reduce file size:
c 200       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f9.3),
c     1          1x,I6,2(1x,f8.3))
c 201       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e9.3),
c     1          1x,I6,2(1x,f8.3))
 200       format(I3,2(1x,f7.2),1x,f8.2,1x,f6.3,4(1x,f8.2),
     1          1x,I5,2(1x,f7.2))
 201       format(I3,2(1x,f7.2),1x,f8.2,1x,f6.3,4(1x,e8.2),
     1          1x,I5,2(1x,f7.2))
c
 1003   continue
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9926,*)"indx, ityp, x, y, z, t, px, py, pz, e, mass
     &, vx, vy, vz"
        DO 480 ihai=1, MUL
           write(9926,479)indx(ihai),ityp(ihai), gx(ihai), gy(ihai),
     &          gz(ihai), ft(ihai), px(ihai), py(ihai), pz(ihai),
     &          e(ihai), xmass(ihai),
     &          vxp(ihai), vyp(ihai), vzp(ihai)
 480    continue
 479    format (2x,i8, 2x,i8, 12(2x,f8.4))
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /PARA5/ i_config, iord_sch
c$$$        i_config : choice of geometric configuration and space cell
c$$$                   dividion optimization.
c$$$          = 1 : the system is undergoing 3-d expansion.
c$$$          = 2 : the system is undergoing 1-d expansion with space cell
c$$$                division.
c$$$          = 3 : the system is undergoing 1-d expansion without space cell
c$$$                division.
c$$$          = 4 : the system is confined in a box with space cell division.
c$$$          = 5 : the system is confined in a box without space cell division.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if (iconfg .le. 3) then
           do 1004 i = 1, mul
              if (ft(i) .le. abs(gz(i))) then
                 eta(i) = 1000000.d0
              else
                 eta(i) = 0.5d0 * log((ft(i) + gz(i)) / (ft(i) - gz(i)))
              end if
              if (e(i) .le. abs(pz(i))) then
                 rap(i) = 1000000.d0
              else
                 rap(i) = 0.5d0 * log((e(i) + pz(i)) / (e(i) - pz(i)))
              end if
clin-8/2015 to avoid IEEE_OVERFLOW_FLAG:
c              tau(i) = ft(i) / cosh(eta(i))
              if(eta(i).lt.1000000.d0) then
                 tau(i) = ft(i) / cosh(eta(i))
              else
                 tau(i) = 1d-10
              endif
c
 1004      continue
           do 1005 i = 1, mul
              etas(i) = eta(i)
              raps(i) = rap(i)
              taus(i) = tau(i)
 1005      continue
        end if
        return
        end
