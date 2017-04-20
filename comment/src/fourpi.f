      real function FOURPI(srt)
*  This function contains the experimental pi+p-->DELTA+PION cross sections   *
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in mb                                    *
*  earray = EXPerimental table with proton energies in MeV                    *
*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
*                                                                             *
******************************************
c      real*4   xarray(10), earray(10)
      real   xarray(10), earray(10)
      SAVE   
      data xarray /0.0001,1.986597,6.411932,7.636956,    
     &9.598362,9.889740,10.24317,10.80138,11.86988,12.83925/    
      data earray /2.468,2.718,2.968,0.322E+01,
     &0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01,
     &0.472E+01/
           pmass=0.14 
       pmass1=0.938
       FOURPI=0.000001
       if(srt.le.1.52)return
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
        plab=SRT
      if (plab .lt. earray(1)) then
        FOURPI= 0.00001
        return
      end if
*
* 2.Interpolate double logarithmically to find sigma(srt)
*
      do 1001 ie = 1,10
        if (earray(ie) .eq. plab) then
          FOURPI= xarray(ie)
          return
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          FOURPI= exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &    /(xmax-xmin) )
          return
        end if
 1001   continue
      return
        END
******************************************
******************************************
* for pion (rho or omega)+baryon resonance collisions
c      real*4 function reab(i1,i2,srt,ictrl)
