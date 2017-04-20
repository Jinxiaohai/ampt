        subroutine wallc1(i, i1, i2, i3, t, tmin)
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
        if (t .lt. size .and. tf .lt. size) then
           if (v1p .gt. 0d0) then
              t1 = ((dble(i1) - 5d0) * size1 - x1p) / v1p + tf
           else if (v1p .lt. 0d0) then
              t1 = ((dble(i1) - 6d0) * size1 - x1p) / v1p + tf
           else
              t1 = tlarge
           end if
           if (v2p .gt. 0d0) then
              t2 = ((dble(i2) - 5d0) * size2 - x2p) / v2p + tf
           else if (v2p .lt. 0d0) then
              t2 = ((dble(i2) - 6d0) * size2 - x2p) / v2p + tf
           else
              t2 = tlarge
           end if
           if (v3p .gt. 0d0) then
              t3 = ((dble(i3) - 5d0) * size3 - x3p) / v3p + tf
           else if (v3p .lt. 0d0) then
              t3 = ((dble(i3) - 6d0) * size3 - x3p) / v3p + tf
           else
              t3 = tlarge
           end if
c       if a particle is on the wall, we don't collide it on the same wall
c        if (t1 .eq. 0d0) t1 = tlarge
c        if (t2 .eq. 0d0) t2 = tlarge
c        if (t3 .eq. 0d0) t3 = tlarge
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
        if (tmin .le. size) return
        end if
        if (v1p .gt. (i1 - 5) * v1) then
           t1 = ((i1 - 5) * (size1 - v1 * size) +
     &          v1p * tf - x1p) / (v1p - (i1 - 5) * v1)
        else if (v1p .lt. (i1 - 6) * v1) then
           t1 = ((i1 - 6) * (size1 - v1 * size) +
     &          v1p * tf - x1p) / (v1p - (i1 - 6) * v1)
        else
           t1 = tlarge
        end if
        if (v2p .gt. (i2 - 5) * v2) then
           t2 = ((i2 - 5) * (size2 - v2 * size) +
     &          v2p * tf - x2p) / (v2p - (i2 - 5) * v2)
        else if (v2p .lt. (i2 - 6) * v2) then
           t2 = ((i2 - 6) * (size2 - v2 * size) +
     &          v2p * tf - x2p) / (v2p - (i2 - 6) * v2)
        else
           t2 = tlarge
        end if
        if (v3p .gt. (i3 - 5) * v3) then
           t3 = ((i3 - 5) * (size3 - v3 * size) +
     &          v3p * tf - x3p) / (v3p - (i3 - 5) * v3)
        else if (v3p .lt. (i3 - 6) * v3) then
           t3 = ((i3 - 6) * (size3 - v3 * size) +
     &          v3p * tf - x3p) / (v3p - (i3 - 6) * v3)
        else
           t3 = tlarge
        end if
c       if a particle is on the wall, we don't collide it on the same wall
c        if (t1 .eq. 0d0) t1 = tlarge
c        if (t2 .eq. 0d0) t2 = tlarge
c        if (t3 .eq. 0d0) t3 = tlarge
        tmin = min(t1, t2, t3)
c       set icsta,
c       after checking this is not an earlier collision comparing with 
c       a collision with another particle, we need to set icsta=0
c       after checking whether there is also a particle collision 
c       at the same time, we need to reset the second bit of icsta
        if (tmin .eq. t1) then
           if (v1p .gt. (i1 - 5) * v1) then
              icsta(i) = 101
           else
              icsta(i) = 102
           end if
        end if
        if (tmin .eq. t2) then
           if (v2p .gt. (i2 - 5) * v2) then
              icsta(i) = 103
           else
              icsta(i) = 104
           end if
        end if
        if (tmin .eq. t3) then
           if (v3p .gt. (i3 - 5) * v3) then
              icsta(i) = 105
           else
              icsta(i) = 106
           end if
        end if
        return
        end
