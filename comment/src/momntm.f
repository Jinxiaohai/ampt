        subroutine momntm(px, py, pz, e)
c       to generate the 3 components of the momentum px, py, pz,
c       from the magnitude of the momentum e
        implicit double precision (a-h,o-z)
        external ran1
        parameter (pi = 3.14159265358979d0)
        common /rndm3/ iseedp
cc      SAVE /rndm3/
        SAVE   
        iseed=iseedp
        cost = 2d0 * ran1(iseed) - 1d0
c     7/20/01:
c        sint = sqrt(1d0 - cost ** 2)
        sint = dsqrt(1d0 - cost ** 2)
        phi = 2d0 * pi * ran1(iseed)
        px = e * sint * cos(phi)
        py = e * sint * sin(phi)
        pz = e * cost
        return
        end
