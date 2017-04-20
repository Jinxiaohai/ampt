        subroutine chkin3(l, i1, i2, i3, t, tmin, nc)
c       this subroutine is used to check collisions for particle inside
c       the cube
c       and update the afftected particles through chkcel
        implicit double precision (a-h, o-z)
        SAVE   
c       itest is a flag to make sure the 111111 cell is checked only once
        itest = 0
        do 1003 i = i1 - 1, i1 + 1
           do 1002 j = i2 - 1, i2 + 1
              do 1001 k =  i3 - 1, i3 + 1
                 if (i .eq. 0) then
                    ia = 10
                 else if (i .eq. 11) then
                    ia = 1
                 else
                    ia = i
                 end if
                 if (j .eq. 0) then
                    ib = 10
                 else if (j .eq. 11) then
                    ib = 1
                 else
                    ib = j
                 end if
                 if (k .eq. 0) then
                    ic = 10
                 else if (k .eq. 11) then
                    ic = 1
                 else
                    ic = k
                 end if
                 call chkcel(l, ia, ib, ic, t, tmin, nc)
 1001         continue
 1002      continue
 1003   continue
        return
        end
