        subroutine chcell(il, i1, i2, i3, last0, t, tmin, nc)
c       this program is used to check through all the particles, except last0
c       in the cell (i1,i2,i3) and see if we can get a particle collision 
c       with time less than the original input tmin ( the collision time of 
c       il with the wall
c       last0 cas be set to 0 if we don't want to exclude last0
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        common /ilist2/ icell, icel(10, 10, 10)
cc      SAVE /ilist2/
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
cc      SAVE /ilist4/
        SAVE   
        if (iconfg .eq. 3 .or. iconfg .eq. 5) then
           jj = ichkpt
           do 1001 j = 1, jj
c     10/24/02 get rid of argument usage mismatch in mintm():
              jmintm=j
              if (j .ne. il .and. j .ne. last0)
     &          call mintm(il, jmintm, tmin, nc)
c     &          call mintm(il, j, tmin, nc)
 1001         continue
           return
        end if
c       set l
        if (i1 .eq. 11 .and. i2 .eq. 11 .and. i3 .eq. 11) then
           l = icell
        else
           l = icel(i1 ,i2, i3)
        end if
        if (l .eq. 0) return
        j = nic(l)
c       if there is only one particle
        if (j .eq. 0) then
c       if it's not il or last0,when last is not wall
           if (l .eq. il .or. l .eq. last0) return
           call mintm(il, l, tmin, nc)
c       if there are many particles
        else
           if (l .ne. il .and. l .ne. last0)
     &        call mintm(il, l, tmin, nc)
           do 10 while(j .ne. l)
              if (j .ne. il .and. j .ne. last0)
     &             call mintm(il, j, tmin, nc)
              j = nic(j)
 10           continue
        end if
        return
        end
