        subroutine zpca1a(i)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
cc      SAVE /para2/
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        common /prec3/gxs(MAXPTN),gys(MAXPTN),gzs(MAXPTN),fts(MAXPTN),
     &     pxs(MAXPTN), pys(MAXPTN), pzs(MAXPTN), es(MAXPTN),
     &     xmasss(MAXPTN), ityps(MAXPTN)
cc      SAVE /prec3/
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
cc      SAVE /prec5/
        common /prec6/ etas(MAXPTN), raps(MAXPTN), taus(MAXPTN)
cc      SAVE /prec6/
        common /ana1/ ts(12)
cc      SAVE /ana1/
        SAVE   
        if (iconfg .eq. 1) then
           t1 = fts(i)
           t2 = ft(i)
           ipic = 11
        else if (iconfg .eq. 2 .or.
     &     iconfg .eq. 3) then
cd           t1 = fts(i)
cd           t2 = ft(i)
           t1 = taus(i)
           t2 = tau(i)
           ipic = 12
        else if (iconfg .eq. 4 .or.
     &     iconfg .eq. 5) then
           t1 = fts(i)
           t2 = ft(i)
           ipic = 12
        end if
        if (iconfg .le. 3) then
           do 1002 ian = 1, ipic
              if (t1 .le. ts(ian) .and.
     &           t2 .gt. ts(ian)) then
                 rapi = raps(i)
c     7/20/01:
c                 et = sqrt(pxs(i) ** 2 + pys(i) ** 2 + xmp ** 2)
                 et = dsqrt(pxs(i) ** 2 + pys(i) ** 2 + xmp ** 2)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$   subroutine zpca1b()更新ana2数据块里面的数值。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
                 call zpca1b(rapi, et, ian)
              end if
 1002      continue
        else
           do 1003 ian = 1, ipic
              if (t1 .le. ts(ian) .and.
     &           t2 .gt. ts(ian)) then
                 p0 = es(i)
                 p1 = pxs(i)
                 p2 = pys(i)
                 p3 = pzs(i)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  subroutine zpca1c()更新ana3数据块里面的数值。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
                 call zpca1c(p0, p1, p2, p3, ian)
              end if
 1003      continue
        end if
        return
        end
