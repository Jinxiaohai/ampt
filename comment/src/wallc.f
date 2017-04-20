        subroutine wallc(i, i1, i2, i3, t, tmin)
c       this subroutine calculates the next time for collision with wall 
c       for particle i
c       input particle label i,t
c       output tmin collision time with wall, icsta(i) wall collision
c       information
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        SAVE   
        tmin = tlarge
        if (iconfg .le. 2 .or. iconfg .eq. 4) then
c       if particle is inside the cube
           if ((i1 .ge. 1 .and. i1 .le. 10)
     &          .or. (i2 .ge. 1 .and. i2 .le. 10)
     &          .or. (i3 .ge. 1 .and. i3 .le. 10)) then
              call wallc1(i, i1, i2, i3, t, tmin)
c       if particle is outside the cube
           else
              call wallcb(i, t, tmin)              
           end if
        else if (iconfg .eq. 3 .or. iconfg .eq. 5) then
           call wallc2(i, i1, i2, i3, t, tmin)
        end if
        return
        end
