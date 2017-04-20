      real function X2pi(srt)
*  This function contains the experimental 
c     total pp-pp+pi(+)pi(-) Xsections    *
*  srt    = DSQRT(s) in GeV                                                  *
*  xsec   = production cross section in mb                                   *
*  earray = EXPerimental table with proton momentum in GeV/c                 *
*  xarray = EXPerimental table with cross sections in mb (curve to guide eye)*
*                                                                            *
******************************************
c      real*4   xarray(15), earray(15)
      real   xarray(15), earray(15)
      SAVE   
      data earray /2.23,2.81,3.67,4.0,4.95,5.52,5.97,6.04,
     &6.6,6.9,7.87,8.11,10.01,16.0,19./
      data xarray /1.22,2.51,2.67,2.95,2.96,2.84,2.8,3.2,
     &2.7,3.0,2.54,2.46,2.4,1.66,1.5/
           pmass=0.9383 
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
c      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
       x2pi=0.000001
       if(srt.le.2.2)return
      plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
      if (plab .lt. earray(1)) then
        x2pi = xarray(1)
        return
      end if
*
* 2.Interpolate double logarithmically to find sigma(srt)
*
      do 1001 ie = 1,15
        if (earray(ie) .eq. plab) then
          x2pi= xarray(ie)
          return
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          X2pi = exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &    /(xmax-xmin) )
          return
        end if
 1001 continue
      return
        END
******************************************
* for pp-->pn+pi(+)pi(+)pi(-)
c      real*4 function X3pi(srt)
