      real function X33pi(srt)
*  This function contains the experimental pp->pp+3pi cross sections          *
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in mb                                    *
*  earray = EXPerimental table with proton energies in MeV                    *
*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
*                                                                             *
******************************************
c      real*4   xarray(12), earray(12)
      real   xarray(12), earray(12)
      SAVE   
      data xarray /0.02,0.22,0.74,1.10,1.76,1.84,2.20,
     &2.40,2.15,2.60,2.30,1.70/
      data earray /2.23,2.81,3.67,4.00,4.95,5.52,5.97,
     &6.04,6.60,6.90,10.01,19./
           pmass=0.9383 
       x33pi=1.E-06
       if(srt.le.2.3)return
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
c      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
      plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
      if (plab .lt. earray(1)) then
        x33pi = xarray(1)
        return
      end if
*
* 2.Interpolate double logarithmically to find sigma(srt)
*
      do 1001 ie = 1,12
        if (earray(ie) .eq. plab) then
          x33pi= xarray(ie)
          return
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          x33pi= exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &    /(xmax-xmin))
          return
        end if
 1001   continue
        return
        END
******************************************
c       REAL*4 FUNCTION X4pi(SRT)
