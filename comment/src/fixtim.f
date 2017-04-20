        subroutine fixtim(l, t, tmin1, tmin, nc)
c       this subroutine is used to compare the collision time with wall tmin1
c       and new collision time with particles for particle l
c       when used in ulist, input nc may be 0, which indicates no particle
c       collisions happen before wall collision, of course, then tmin=tmin1
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        SAVE   
        k = nc
        if (tmin .lt. tmin1) then
           ot(l) = tmin
           if (ct(l) .lt. tmin1) then
              icsta(l) = 0
           else
              icsta(l) = icsta(l) + 10
           end if
           next(l) = k
        else if (tmin .eq. tmin1) then
           ot(l) = tmin
           if (nc .eq. 0) then
              next(l) = 0
           else
              icsta(l) = icsta(l) + 10
              next(l) = k
           end if
        else
           ot(l) = tmin1
           next(l) = 0
        end if
        return
        end
