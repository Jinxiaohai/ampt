      subroutine zprota(xn1, xn2, xn3, theta, v1, v2, v3)
cbz1/29/99end
c     this subroutine is used to rotate the vector (v1,v2,v3) by an angle theta
c     around the unit vector (xn1, xn2, xn3)
      implicit double precision (a-h, o-z)
      SAVE   
      vx = v1
      vy = v2
      vz = v3
      c = cos(theta)
      omc = 1d0 - c
      s = sin(theta)
      a11 = xn1 ** 2 * omc + c
      a12 = xn1 * xn2 * omc - s * xn3
      a13 = xn1 * xn3 * omc + s * xn2
      a21 = xn1 * xn2 * omc + s * xn3
      a22 = xn2 **2 * omc + c
      a23 = xn2 * xn3 * omc - s * xn1
      a31 = xn1 * xn3 * omc - s * xn2
      a32 = xn3 * xn2 * omc + s * xn1
      a33 = xn3 ** 2 * omc + c
      v1 = vx * a11 + vy * a12 + vz * a13
      v2 = vx * a21 + vy * a22 + vz * a23
      v3 = vx * a31 + vy * a32 + vz * a33
      return
      end
