      real function pp1(srt)
      SAVE   
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in mb                                    *
*  earray = EXPerimental table with proton energies in MeV                    *
*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
*                                                                             *
******************************************
           pmass=0.9383 
       PP1=0.
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
c      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
      plab2=((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2
       IF(PLAB2.LE.0)RETURN
      plab=sqrt(PLAB2)
       pmin=0.968
       pmax=2080
      if ((plab .lt. pmin).or.(plab.gt.pmax)) then
        pp1 = 0.
        return
      end if
c* fit parameters
       a=30.9
       b=-28.9
       c=0.192
       d=-0.835
       an=-2.46
        pp1 = a+b*(plab**an)+c*(alog(plab))**2
       if(pp1.le.0)pp1=0.0
        return
        END
******************************************
* for pp-->elastic
c      real*4 function pp2(srt)
