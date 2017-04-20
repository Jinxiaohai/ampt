        subroutine newcre(i, k)
c       this subroutine is used to mk rearrange of the new cell a particle
c       enters,
c       input i
c       output nic(i)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        SAVE   
        if (k .eq. 0) then
           k = i
           nic(i) = 0
        else if (nic(k) .eq. 0) then
           nic(k) = i
           nic(i) = k
        else
           j = k
           do 10 while (nic(j) .ne. k)
              j = nic(j)
 10           continue
           nic(j) = i
           nic(i) = k
        end if
        return
        end
