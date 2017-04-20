      real function ppk1(srt)
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in mb                                    *
*                                                                             *
******************************************
c      real*4   xarray(7), earray(7)
      real   xarray(7), earray(7)
      SAVE   
      data xarray /0.013,0.025,0.016,0.012,0.017,0.029,0.025/
      data earray /3.67,4.95,5.52,5.97,6.05,6.92,7.87/
           pmass=0.9383 
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
c      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
       ppk1=0.
       if(srt.le.2.63)return
       if(srt.gt.4.08)then
       ppk1=0.025
       return
       endif
        plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
        if (plab .lt. earray(1)) then
        ppk1 =xarray(1)
        return
      end if
*
* 2.Interpolate double logarithmically to find sigma(srt)
*
      do 1001 ie = 1,7
        if (earray(ie) .eq. plab) then
          ppk1 = xarray(ie)
          go to 10
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          ppk1 = exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &/(xmax-xmin) )
          go to 10
        end if
 1001 continue
10       continue
      return
        END
**********************************
*                                                                      *
*                                                                      *
