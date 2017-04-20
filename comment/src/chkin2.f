        subroutine chkin2(l, i1, i2, i3, t, tmin, nc)
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
                 ia = i
                 ib = j
                 ic = k
                 if (k .ge. 1 .and. k .le. 10) then
                    if (i .eq. 0) ia = 10
                    if (i .eq. 11) ia = 1
                    if (j .eq. 0) ib = 10
                    if (j .eq. 11) ib = 1
                    call chkcel(l, ia, ib, ic, t, tmin, nc)
                 end if
 1001         continue
 1002      continue
 1003   continue
        return
        end
