        subroutine chkout(l, t, tmin, nc)
c       this subroutine is used to check the collisions with particles in 
c       surface cells to see if we can get a smaller collision time than tmin
c       with particle nc, when the colliding particle is outside the cube
c       input l,t,tmin,nc
c       output tmin, nc
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        SAVE   
        m1 = 11
        m2 = 11
        m3 = 11
        call chkcel(l, m1, m2, m3, t, tmin, nc)
        do 1003 i = 1, 10
           do 1002 j = 1, 10
              do 1001 k = 1, 10
                 if (i .eq. 1 .or. i .eq. 10 .or. j .eq. 1
     &              .or. j .eq. 10 .or. k .eq. 1 .or. k .eq. 10) 
     &                    call chkcel(l, i, j, k, t, tmin, nc)
 1001         continue
 1002      continue
 1003   continue
        return
        end
