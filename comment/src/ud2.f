        subroutine ud2(i, j, t, tmin, nc)
c       this subroutine is used to update next(i), ct(i), ot(i),
c        and get tmin, nc for j
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        common /aurec1/ jxa, jya, jza
cc      SAVE /aurec1/
        common /aurec2/ dgxa(MAXPTN), dgya(MAXPTN), dgza(MAXPTN)
cc      SAVE /aurec2/
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
        logical allok
        call isco(i, j, allok, tm, t1, t2)
        if (allok) then
c       tm eq tmin, change nc to make sure fixtime get the collision with both 
c       wall and particle
             if (tm .lt. tmin) then
                tmin = tm
                ct(j) = t2
                nc = i
                if (iconfg .eq. 3 .or. iconfg .eq. 5) then
                   dgxa(j) = jxa * 10d0 * size1
                   dgya(j) = jya * 10d0 * size2
                   if (iconfg .eq. 5) then
                      dgza(j) = jza * 10d0 * size3
                   end if
                end if
             end if
             if (tm .le. ot(i)) then
                ct(i) = t1
                icels0 = icels(i)
                i1 = icels0 / 10000
                i2 = (icels0 - i1 * 10000) / 100
                i3 = icels0 - i1 * 10000 - i2 * 100
                call wallc(i, i1, i2, i3, t, tmin1)
                call fixtim(i, t, tmin1, tm, j)
                if (iconfg .eq. 3 .or. iconfg .eq. 5) then
                   dgxa(i) = - jxa * 10d0 * size1
                   dgya(i) = - jya * 10d0 * size2
                   if (iconfg .eq. 5) then
                      dgza(i) = - jza * 10d0 * size3
                   end if
                end if
             end if
             if (tm .gt. ot(i) .and. next(i) .eq. j) then
                ct(i) = t1
                call reor(t, tm, i, j)
             end if
           else if (next(i) .eq. j) then
             tm = tlarge
             call reor(t, tm, i, j)
          end if
        return
        end
