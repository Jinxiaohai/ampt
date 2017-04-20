        subroutine getict(t1)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para1/ mul
cc      SAVE /para1/
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
cc      SAVE /ilist4/
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        SAVE   
c       neglect possibility of 2 collisions at the same time
c0       set initial conditions
        t1 = tlarge
        iscat = 0
        jscat = 0
c1      get next collision between particles
        do 1001 i = 1, ichkpt
           if (ot(i) .lt. t1) then
              t1 = ot(i)
              iscat = i
           end if
 1001   continue
        if (iscat .ne. 0) jscat = next(iscat)
c2      get ictype
c     10/30/02 ictype=0:collision; 1:parton formation
        if (iscat .ne. 0 .and. jscat .ne. 0) then
           if (icsta(iscat) .eq. 0 .and. icsta(jscat) .eq. 0) then
              ictype = 0
           else
              ictype = 4
           end if
        else if (iscat .ne. 0 .or. jscat .ne. 0) then
           ictype = 3
        end if
c
        if (ifmpt .le. mul) then
           if (ft(ifmpt) .lt. t1) then
              ictype = 1
              t1 = ft(ifmpt)
           else if (ft(ifmpt) .eq. t1) then
              if (ictype .eq. 0) ictype = 2
              if (ictype .eq. 3) ictype = 5
              if (ictype .eq. 4) ictype = 6
           end if
        end if
        return
        end
