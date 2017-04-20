        subroutine isco12(i, j, allok, tm, t1, t2)
c       this subroutine is used to decide whether there is a collision between
c       particle i and j, if there is one allok=1, and tm gives the 
c       collision time, t1 the collision time for i,
c       t2 the collision time for j
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
cc      SAVE /para2/
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
cc      SAVE /prec4/
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
        logical allok, allokp
c       preventing consecutive collisions
        allok = last(i) .ne. j .or. last(j) .ne. i
        if (ft(i) .ge. ft(j)) then
           i1 = j
           i2 = i
           isign = -1
        else 
           i1 = i
           i2 = j
           isign = 1
        end if
        if (allok) then
           tm = tlarge
           t1 = ft(i1)
           vx1 = vx(i1)
           vy1 = vy(i1)
           vz1 = vz(i1)
           t2 = ft(i2)
           dvx = vx(i2) - vx1
           dvy = vy(i2) - vy1
           dvz = vz(i2) - vz1
           dt = t2 - t1
           do 1000 ii = - 1, 1
              do 2000 jj = - 1, 1
                 do 3000 kk = -1, 1
                 allokp = .true.
                 dx = gx(i2) - gx(i1) - vx1 * dt
                 dy = gy(i2) - gy(i1) - vy1 * dt
                 dz = gz(i2) - gz(i1) - vz1 * dt
                 dx = dx + ii * 10d0 * size1
                 dy = dy + jj * 10d0 * size2
                 dz = dz + kk * 10d0 * size3
                 vp = dvx * dx + dvy * dy + dvz * dz
                 allokp = allokp .and. vp .lt. 0d0
                 if (allokp) then
                    v2 = dvx * dvx + dvy * dvy + dvz * dvz
                    if (v2 .eq. 0d0) then
                       tmp = tlarge
                    else
                       tmp = t2 - vp / v2
                    end if
c       note now tm is the absolute time
                    allokp = allokp .and. tmp .gt. t1 .and.
     &                         tmp .gt. t2
                 end if
                 if (allokp) then
                    dgx = dx - dvx * t2
                    dgy = dy - dvy * t2
                    dgz = dz - dvz * t2
                    dm2 = - v2 * tmp ** 2  + dgx * dgx +
     &                    dgy * dgy + dgz * dgz
                    allokp = allokp .and. dm2 .lt. cutof2
                 end if
                 if (allokp .and. tmp .lt. tm) then
                    tm = tmp
                    jxa = isign * ii
                    jya = isign * jj
                    jza = isign * kk
                 end if
 3000                 continue
 2000              continue
 1000           continue
           if (tm .eq. tlarge) then
              allok = .false.
           end if
        end if
        if (allok) then
           e1 = e(i1)
           px1 = px(i1)
           py1 = py(i1)
           pz1 = pz(i1)
           e2 = e(i2)
           px2 = px(i2)
           py2 = py(i2)
           pz2 = pz(i2)
           rts2 = (e1 + e2) ** 2 - (px1 + px2) ** 2
     &          - (py1 + py2) ** 2 - (pz1 + pz2) ** 2
           allok = allok .and. rts2 .gt. rscut2
        end if
        if (.not. allok) then
           tm = tlarge
           t1 = tlarge
           t2 = tlarge
        else
           t1 = tm
           t2 = tm
        end if
        return
        end
