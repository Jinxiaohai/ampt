      real function pp2(srt)
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
          write(99,*) 'scheck33: ', scheck
          scheck=0.
       endif
       plab=sqrt(scheck)
c      plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
       pmin=2.
       pmax=2050
       if(plab.gt.pmax)then
       pp2=8.
       return
       endif
        if(plab .lt. pmin)then
        pp2 = 25.
        return
        end if
c* fit parameters
       a=11.2
       b=25.5
       c=0.151
       d=-1.62
       an=-1.12
        pp2 = a+b*(plab**an)+c*(alog(plab))**2+d*alog(plab)
       if(pp2.le.0)pp2=0
        return
        END
******************************************
* for pp-->total
c      real*4 function ppt(srt)
