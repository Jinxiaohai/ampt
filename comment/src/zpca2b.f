        subroutine zpca2b
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
cc      SAVE /ilist4/
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        common /ana1/ ts(12)
cc      SAVE /ana1/
        SAVE   
        do 1002 i = 1, ichkpt
           t1 = ft(i)
           t2 = tlarge
           ipic = 12
           do 1001 ian = 1, ipic
              if (t1 .le. ts(ian) .and.
     &           t2 .gt. ts(ian)) then
                 p0 = e(i)
                 p1 = px(i)
                 p2 = py(i)
                 p3 = pz(i)
                 call zpca1c(p0, p1, p2, p3, ian)
              end if
 1001      continue
 1002   continue
        return
        end
