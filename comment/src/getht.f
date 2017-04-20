        subroutine getht(iscat, jscat, pp2, that)
c       this subroutine is used to get \hat{t} for a particular processes
        implicit double precision (a-h, o-z)
        parameter (hbarc = 0.197327054d0)
        parameter (MAXPTN=400001)
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
cc      SAVE /para2/
        common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
        external ran1
        common /rndm3/ iseedp
cc      SAVE /rndm3/
        SAVE   
        iseed=iseedp
        xmu2 = (hbarc * xmu) ** 2
        xmp2 = xmp ** 2
        xm2 = xmu2 + xmp2
        rx=ran1(iseed)
        that = xm2*(1d0+1d0/((1d0-xm2/(4d0*pp2+xm2))*rx-1d0))
ctest off isotropic scattering:
c     &     + 1d0/((1d0 - xm2 / (4d0 * pp2 + xm2)) * ran1(2) - 1d0))
c        if(izpc.eq.100) that=-4d0*pp2*ran1(2)
        if(izpc.eq.100) that=-4d0*pp2*rx
        return
        end
