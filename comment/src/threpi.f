      real function THREPI(srt)
*  This function contains the experimental pi+p-->DELTA + rho cross sections  *
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in mb                                    *
*  earray = EXPerimental table with proton energies in MeV                    *
*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
*                                                                             *
******************************************
c      real*4   xarray(15), earray(15)
      real   xarray(15), earray(15)
      SAVE   
      data xarray /8.0000000E-06,6.1999999E-05,1.881940,5.025690,    
     &11.80154,13.92114,15.07308,11.79571,11.53772,10.01197,9.792673,    
     &9.465264,8.970490,7.944254,6.886320/    
      data earray /0.122E+01, 0.147E+01, 0.172E+01, 0.197E+01,
     &0.222E+01, 0.247E+01, 0.272E+01, 0.297E+01, 0.322E+01,
     &0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01,
     &0.472E+01/
           pmass=0.14 
       pmass1=0.938
       THREPI=0.000001
       if(srt.le.1.36)return
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
        plab=SRT
      if (plab .lt. earray(1)) then
        THREPI = 0.00001
        return
      end if
*
* 2.Interpolate double logarithmically to find sigma(srt)
*
      do 1001 ie = 1,15
        if (earray(ie) .eq. plab) then
          THREPI= xarray(ie)
          return
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          THREPI = exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &    /(xmax-xmin) )
          return
        end if
 1001   continue
      return
        END
******************************************
******************************************
* for ppi(+)-->DELTA+omega
c      real*4 function FOURPI(srt)
