      real function ppt(srt)
      SAVE   
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in mb                                    *
*  earray = EXPerimental table with proton energies in MeV                    *
*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
*                                                                             *
******************************************
           pmass=0.9383 
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
c      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
clin-9/2012: check argument in sqrt():
       scheck=((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2
       if(scheck.lt.0) then
          write(99,*) 'scheck34: ', scheck
          scheck=0.
       endif
       plab=sqrt(scheck)
c      plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
       pmin=3. 
       pmax=2100
      if ((plab .lt. pmin).or.(plab.gt.pmax)) then
        ppt = 55.
        return
      end if
c* fit parameters
       a=45.6
       b=219.0
       c=0.410
       d=-3.41
       an=-4.23
        ppt = a+b*(plab**an)+c*(alog(plab))**2+d*alog(plab)
       if(ppt.le.0)ppt=0.0
        return
        END
*************************
* cross section for N*(1535) production in PP collisions
* VARIABLES:
* LB1,LB2 ARE THE LABLES OF THE TWO COLLIDING PARTICLES
* SRT IS THE CMS ENERGY
* X1535 IS THE N*(1535) PRODUCTION CROSS SECTION
* NOTE THAT THE N*(1535) PRODUCTION CROSS SECTION IS 2 TIMES THE ETA 
* PRODUCTION CROSS SECTION
* DATE: Aug. 1 , 1994
* ********************************
