      subroutine xnormv(vx, vy, vz)
c      this subroutine is used to get a normalized vector 
      implicit double precision (a-h, o-z)
      SAVE   
clin-7/20/01:
c      vv = sqrt(vx ** 2 + vy ** 2 + vz ** 2)
      vv = dsqrt(vx ** 2 + vy ** 2 + vz ** 2)
      vx = vx / vv
      vy = vy / vv
      vz = vz / vv
      return
      end
cbz1/29/99
c      subroutine rotate(xn1, xn2, xn3, theta, v1, v2, v3)
