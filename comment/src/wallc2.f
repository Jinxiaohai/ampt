        subroutine wallc2(i, i1, i2, i3, t, tmin)
c       this subroutine is used to get wall collision time
c       when particle is inside the cube, it sets the icsta at the same time
c       input i,i1,i2,i3,t
c       output tmin, icsta(i)
c       note the icsta is not finally set. we need further judgement in 
c       fixtim
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
cc      SAVE /prec4/
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
cc      SAVE /ilist3/
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        SAVE   
        x1p = gx(i)
        x2p = gy(i)
        x3p = gz(i)
        tf = ft(i)
        v1p = vx(i)
        v2p = vy(i)
        v3p = vz(i)
        if (v1p .gt. 0d0) then
           t1 = (5d0 * size1 - x1p) / v1p + tf
        else if (v1p .lt. 0d0) then
           t1 = (-5d0 * size1 - x1p) / v1p + tf
        else
           t1 = tlarge
        end if
        if (v2p .gt. 0d0) then
           t2 = (5d0 * size2 - x2p) / v2p + tf
        else if (v2p .lt. 0d0) then
           t2 = (- 5d0 * size2 - x2p) / v2p +tf
        else
           t2 = tlarge
        end if
        if (iconfg .eq. 5) then
           if (v3p .gt. 0d0) then
              t3 = (5d0 * size3 - x3p) / v3p + tf
           else if (v3p .lt. 0d0) then
              t3 = (- 5d0 * size3 - x3p) / v3p +tf
           else
              t3 = tlarge
           end if
        else
           t3 = tlarge
        end if
        tmin = min(t1, t2, t3)
c       set icsta,
c       after checking this is not an earlier collision comparing with 
c       a collision with another particle, we need to set icsta=0
c       after checking whether there is also a particle collision 
c       at the same time, we need to reset the second bit of icsta
        if (tmin .eq. t1) then
           if (v1p .gt. 0d0) then
              icsta(i) = 101
           else
              icsta(i) = 102
           end if
        end if
        if (tmin .eq. t2) then
           if (v2p .gt. 0d0) then
              icsta(i) = 103
           else
              icsta(i) = 104
           end if
        end if
        if (tmin .eq. t3) then
           if (v3p .gt. 0d0) then
              icsta(i) = 105
           else
              icsta(i) = 106
           end if
        end if
        return
        end
