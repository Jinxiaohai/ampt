      subroutine cropro(vx1, vy1, vz1, vx2, vy2, vz2)
c     this subroutine is used to calculate the cross product of 
c     (vx1,vy1,vz1) and (vx2,vy2,vz2) and get the result (vx3,vy3,vz3)
c     and put the vector into common /cprod/
      implicit double precision (a-h, o-z)
      common/cprod/ vx3, vy3, vz3
cc      SAVE /cprod/
      SAVE   
      vx3 = vy1 * vz2 - vz1 * vy2
      vy3 = vz1 * vx2 - vx1 * vz2
      vz3 = vx1 * vy2 - vy1 * vx2
      return
      end
