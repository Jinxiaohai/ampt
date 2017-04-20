      real function pipp1(srt)
      SAVE   
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in fm**2                                 *
*  earray = EXPerimental table with proton energies in MeV                    *
*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
*  UNITS: FM**2
******************************************
           pmass=0.14 
       pmass1=0.938
       PIPP1=0.0001
       IF(SRT.LE.1.22)RETURN
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
c      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
        plab=sqrt(((srt**2-pmass**2-pmass1**2)/(2.*pmass1))**2-pmass**2)
       pmin=0.3
       pmax=25.0
       if(plab.gt.pmax)then
       pipp1=20./10.
       return
       endif
        if(plab .lt. pmin)then
        pipp1 = 0.
        return
        end if
c* fit parameters
       a=26.6
       b=-7.18
       c=0.327
       an=-1.86
       d=-2.81
        pipp1 = a+b*(plab**an)+c*(alog(plab))**2+d*alog(plab)
       if(pipp1.le.0)pipp1=0
       PIPP1=PIPP1/10.
        return
        END
* *****************************
c       real*4 function xrho(srt)
