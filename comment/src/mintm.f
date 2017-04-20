        subroutine mintm(i, j, tmin, nc)
c       this subroutine is used to check whether particle j has smaller
c       collision time with particle i than other particles
c       or in other words, update next(i)
c       input i,j,tmin,nc
c       output tmin,nc
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
        if (allok .and. tm .lt. tmin) then
           tmin = tm
           ct(i) = t1
           nc = j
           if (iconfg .eq. 3 .or. iconfg .eq. 5) then
              dgxa(i) = - jxa * 10d0 * size1
              dgya(i) = - jya * 10d0 * size2
              if (iconfg .eq. 5) then
                 dgza(i) = - jza * 10d0 * size3
              end if
           end if
        end if
         return
        end
******************************************************************************
******************************************************************************
