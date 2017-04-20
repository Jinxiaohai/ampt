        subroutine wallcb(i, t, tmin)
c       this subroutine is used to calculate the wall collision time 
c       when the particle is outside the cube
c       input i,t
c       output tmin,icsta(i)
c       note the icsta is not finally set. we need further judgement in 
c       fixtim
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
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
c       check if there is a collision by looking at the closest approach point
c       and see if it's inside the cube
        if (size1 .eq. 0d0 .or. size2 .eq. 0d0 .or. 
     &     size3 .eq. 0d0) return
        x1p = gx(i)
        x2p = gy(i)
        x3p = gz(i)
        v1p = vx(i)
        v2p = vy(i)
        v3p = vz(i)
        tf = ft(i)
        if (t .lt. size .and. tf .lt. size) then
           if (x1p .lt. - 5d0 * size1 .and. v1p .gt. 0d0) then
              t1 = (- 5d0 * size1 - x1p) / v1p + tf
           else if(x1p .gt. 5d0 * size1 .and. v1p .lt. 0d0) then
              t1 = - (x1p - 5d0 * size1) / v1p + tf
           else
              t1 = tlarge 
           end if
           if (t1 .ne. tlarge) then
              x2pp = x2p + v2p * (t1 - tf)
              x3pp = x3p + v3p * (t1 - tf)
              if (x2pp .le. - 5d0 * size2 .or. x2pp .ge. 5d0 * size2
     &             .or. x3pp .le. - 5d0 * size3 
     &             .or. x3pp .ge. 5d0 * size3)
     &             t1 = tlarge
           end if
           if (x2p .lt. - 5d0 * size2 .and. v2p .gt. 0d0) then
              t2 = (- 5d0 * size2 - x2p) / v2p + tf
           else if(x2p .gt. 5d0 * size2 .and. v2p .lt. 0d0) then
              t2 = - (x2p - 5d0 * size2) / v2p + tf
           else
              t2 = tlarge 
           end if
           if (t2 .ne. tlarge) then
              x1pp = x1p + v1p * (t2 - tf)
              x3pp = x3p + v3p * (t2 - tf)
              if (x1pp .le. - 5d0 * size1 .or. x1pp .ge. 5d0 * size1
     &          .or. x3pp .le. - 5d0 * size3 .or. x3pp .ge. 5d0 * size3)
     &             t2 = tlarge
           end if
           if (x3p .lt. - 5d0 * size3 .and. v3p .gt. 0d0) then
              t3 = (- 5d0 * size3 - x3p) / v3p + tf
           else if(x3p .gt. 5d0 * size3 .and. v3p .lt. 0d0) then
              t3 = - (x3p - 5d0 * size3) / v3p + tf
           else
              t3 = tlarge 
           end if
           if (t3 .ne. tlarge) then
              x1pp = x1p + v1p * (t3 - tf)
              x2pp = x2p + v2p * (t3 - tf)
              if (x1pp .le. - 5d0 * size1 .or. x1pp .ge. 5d0 * size1
     &          .or. x2pp .le. - 5d0 * size2 .or. x2pp .ge. 5d0 * size2)
     &             t3 = tlarge
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
        if (tmin .le. size) return
        end if
c       notice now x1q, x2q, x3q are coordinates at time t
        x1q = x1p + v1p * (t - tf)
        x2q = x2p + v2p * (t - tf)
        x3q = x3p + v3p * (t - tf)
        if (x1q .lt. - 5d0 * (size1 + v1 * (t - size)) .and. 
     &      v1p .gt. - 5d0 * v1) then
           t1 = (- 5d0 * (size1 - v1 * size) + v1p * tf - x1p) /
     &          (v1p - (- 5d0) * v1)
           icsta1 = 101
        else if (x1q .gt. 5d0 * (size1 + v1 * (t-size)) .and. 
     &     v1p .lt. 5d0 * v1) then
           t1 = (5d0 * (size1 - v1 * size) + v1p * tf - x1p) /
     &          (v1p - 5d0 * v1)
           icsta1 = 102
        else
           t1 = tlarge 
        end if
        if (t1 .ne. tlarge) then
           x2pp = x2p + v2p * (t1 - tf)
           x3pp = x3p + v3p * (t1 - tf)
           if (x2pp .le. - 5d0 * (size2 + v2 * (t1 - size))
     &        .or. x2pp .ge. 5d0 * (size2 + v2 * (t1 - size))
     &        .or. x3pp .le. - 5d0 * (size3 + v3 * (t1 - size))
     &        .or. x3pp .ge. 5d0 * (size3 + v3 * (t1 - size)))
     &        t1 = tlarge
        end if
        if (x2q .lt. - 5d0 * (size2 + v2 * (t - size)) .and.
     &     v2p .gt. - 5d0 * v2) then
           t2 = (- 5d0 * (size2 - v2 * size) + v2p * tf - x2p) /
     &          (v2p - (- 5d0) * v2)
           icsta2 = 103
        else if (x2q .gt. 5d0 * (size2 + v2 * (t - size)) .and.
     &     v2p .lt. 5d0 * v2) then
           t2 = (5d0 * (size2 - v2 * size) + v2p * tf - x2p) / 
     &          (v2p - 5d0 * v2)
           icsta2 = 104
        else
           t2 = tlarge 
        end if
        if (t2 .ne. tlarge) then
           x1pp = x1p + v1p * (t2 - tf)
           x3pp = x3p + v3p * (t2 - tf)
           if (x1pp .le. - 5d0 * (size1 + v1 * (t2 - size))
     &        .or. x1pp .ge. 5d0 * (size1 + v1 * (t2 - size))
     &        .or. x3pp .le. - 5d0 * (size3 + v3 * (t2 - size))
     &        .or. x3pp .ge. 5d0 * (size3 + v3 * (t2 - size)))
     &        t2 = tlarge
        end if
        if (x3q .lt. - 5d0 * (size3 + v3 * (t - size)) .and. 
     &     v3p .gt. - 5d0 * v3) then
           t3 = (- 5d0 * (size3 - v3 * size) + v3p * tf - x3p) /
     &          (v3p - (- 5d0) * v3)
           icsta3 = 105
        else if (x3q .gt. 5d0 * (size3 + v3 * (t - size)) .and.
     &     v3p .lt. 5d0 * v3) then
           t3 = (5d0 * (size3 - v3 * size) + v3p * tf - x3p) /
     &          (v3p - 5d0 * v3)
           icsta3 = 106
        else
           t3 = tlarge 
        end if
        if (t3 .ne. tlarge) then
           x2pp = x2p + v2p * (t3 - tf)
           x1pp = x1p + v1p * (t3 - tf)
           if (x2pp .le. - 5d0 * (size2 + v2 * (t3 - size))
     &        .or. x2pp .ge. 5d0 * (size2 + v2 * (t3 - size))
     &        .or. x1pp .le. - 5d0 * (size1 + v1 * (t3 - size))
     &        .or. x1pp .ge. 5d0 * (size1 + v1 * (t3 - size)))
     &        t3 = tlarge
        end if
        tmin = min(t1, t2, t3)
c       set icsta,
c       after checking this is not an earlier collision comparing with 
c       a collision with another particle, we need to set icsta=0
c       after checking whether there is also a particle collision 
c       at the same time, we need to reset the second bit of icsta
        if (tmin .eq. t1) then
           icsta(i) = icsta1
        else if (tmin .eq. t2) then
           icsta(i) = icsta2
        else if (tmin .eq. t3) then
           icsta(i) = icsta3
        end if
        return
        end
