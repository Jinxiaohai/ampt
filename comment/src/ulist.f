      subroutine ulist(t)
c     this subroutine is used to update a new collision time list
c       notice this t has been updated
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
cc      SAVE /ilist4/
        SAVE   

        if (ictype .eq. 1 .or. ictype .eq. 2 .or. ictype .eq. 5
     &     .or. ictype .eq. 6) then
           l = ifmpt
           call ulist1(l, t)
        end if

        if (ictype .ne. 1) then
           l = iscat
           call ulist1(l, t)
           if (jscat .ne. 0) then
              l = jscat
              call ulist1(l, t)
           end if
        end if

        return
        end
