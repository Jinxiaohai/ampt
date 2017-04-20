      real function pionpp(srt)
      SAVE   
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in fm**2                                 *
*  earray = EXPerimental table with proton energies in MeV                    *
*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
*                                                                             *
******************************************
           pmass=0.14 
       pmass1=0.938
       PIONPP=0.00001
       IF(SRT.LE.1.22)RETURN
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
c      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
        plab=sqrt(((srt**2-pmass**2-pmass1**2)/(2.*pmass1))**2-pmass**2)
       pmin=0.3
       pmax=25.0
       if(plab.gt.pmax)then
       pionpp=20./10.
       return
       endif
        if(plab .lt. pmin)then
        pionpp = 0.
        return
        end if
c* fit parameters
       a=24.3
       b=-12.3
       c=0.324
       an=-1.91
       d=-2.44
        pionpp = a+b*(plab**an)+c*(alog(plab))**2+d*alog(plab)
       if(pionpp.le.0)pionpp=0
       pionpp=pionpp/10.
        return
        END
**********************************
* elementary cross sections
*  from the CERN data book
*  date: Sept.2, 1994
*  for pion-+p-->INELASTIC
c      real*4 function pipp1(srt)
