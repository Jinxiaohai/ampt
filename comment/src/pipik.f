      real function pipik(srt)
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in mb                                    *
*  NOTE: DEVIDE THE CROSS SECTION TO OBTAIN K+ PRODUCTION                     *
******************************************
c      real*4   xarray(5), earray(5)
      real   xarray(5), earray(5)
      SAVE   
      data xarray /0.001, 0.7,1.5,1.7,2.0/
      data earray /1.,1.2,1.6,2.0,2.4/
           pmass=0.9383 
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
c      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
       pipik=0.
       if(srt.le.1.)return
       if(srt.gt.2.4)then
           pipik=2.0/2.
           return
       endif
        if (srt .lt. earray(1)) then
           pipik =xarray(1)/2.
           return
        end if
*
* 2.Interpolate double logarithmically to find sigma(srt)
*
      do 1001 ie = 1,5
        if (earray(ie) .eq. srt) then
          pipik = xarray(ie)
          go to 10
        else if (earray(ie) .gt. srt) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          pipik = exp(ymin + (alog(srt)-xmin)*(ymax-ymin)
     &/(xmax-xmin) )
          go to 10
        end if
 1001 continue
10       PIPIK=PIPIK/2.
       continue
      return
        END
**********************************
* TOTAL PION-P INELASTIC CROSS SECTION 
*  from the CERN data book
*  date: Sept.2, 1994
*  for pion++p-->Delta+pion
c      real*4 function pionpp(srt)
