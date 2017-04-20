      real function pplpk(srt)
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
*   find the center of mass energy corresponding to the given pm as
*   if Lambda+N+K are produced
       pplpk=0.
clin-9/2012: check argument in sqrt():
       scheck=((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2
       if(scheck.lt.0) then
          write(99,*) 'scheck35: ', scheck
          scheck=0.
       endif
       plab=sqrt(scheck)
c        plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
       pmin=2.82
       pmax=25.0
       if(plab.gt.pmax)then
       pplpk=0.036
       return
       endif
        if(plab .lt. pmin)then
        pplpk = 0.
        return
        end if
c* fit parameters
       a=0.0654
       b=-3.16
       c=-0.0029
       an=-4.14
        pplpk = a+b*(plab**an)+c*(alog(plab))**2
       if(pplpk.le.0)pplpk=0
        return
        END
******************************************
* for pp-->pSigma+K0
c      real*4 function ppk0(srt)
