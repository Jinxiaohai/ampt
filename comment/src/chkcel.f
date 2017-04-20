        subroutine chkcel(il, i1, i2, i3, t, tmin, nc)
c       this program is used to check through all the particles
c       in the cell (i1,i2,i3) and see if we can get a particle collision 
c       with time less than the original input tmin ( the collision time of 
c       il with the wall
c       and update the affected particles
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
              call ck(j, ick)
c     10/24/02 get rid of argument usage mismatch in ud2():
                            jud2=j
c              if (ick .eq. 1) call ud2(j, il, t, tmin, nc)
              if (ick .eq. 1) call ud2(jud2, il, t, tmin, nc)
 1001      continue
           return
        end if
        if (i1 .eq. 11 .and. i2 .eq. 11 .and. i3 .eq. 11) then
           l = icell
        else
           l = icel(i1, i2, i3)
        end if
c       if there is no particle
        if (l .eq. 0) then
           return
        end if
        j = nic(l)
c       if there is only one particle
        if (j .eq. 0) then
           call ck(l, ick)
           if (ick .eq. 1) call ud2(l, il, t, tmin, nc)
c       if there are many particles
        else
c       we don't worry about the other colliding particle because it's
c       set in last(), and will be checked in ud2
           call ck(l, ick)
           if (ick .eq. 1) call ud2(l, il, t, tmin, nc)
           do 10 while(j .ne. l)
              call ck(j, ick)
              if (ick .eq. 1) call ud2(j, il, t, tmin, nc)
              j = nic(j)
 10           continue
        end if
        return
        end
