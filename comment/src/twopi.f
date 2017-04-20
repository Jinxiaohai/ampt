      real function TWOPI(srt)
*  This function contains the experimental pi+p-->DELTA+PION cross sections   *
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in mb                                    *
*  earray = EXPerimental table with proton energies in MeV                    *
*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
*                                                                             *
******************************************
c      real*4   xarray(19), earray(19)
      real   xarray(19), earray(19)
      SAVE   
      data xarray /0.300E-05,0.187E+01,0.110E+02,0.149E+02,0.935E+01,
     &0.765E+01,0.462E+01,0.345E+01,0.241E+01,0.185E+01,0.165E+01,
     &0.150E+01,0.132E+01,0.117E+01,0.116E+01,0.100E+01,0.856E+00,
     &0.745E+00,0.300E-05/
      data earray /0.122E+01, 0.147E+01, 0.172E+01, 0.197E+01,
     &0.222E+01, 0.247E+01, 0.272E+01, 0.297E+01, 0.322E+01,
     &0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01,
     &0.472E+01, 0.497E+01, 0.522E+01, 0.547E+01, 0.572E+01/
           pmass=0.14 
       pmass1=0.938
       TWOPI=0.000001
       if(srt.le.1.22)return
* 1.Calculate p(lab)  from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
        plab=SRT
      if (plab .lt. earray(1)) then
        TWOPI= 0.00001
        return
      end if
*
* 2.Interpolate double logarithmically to find sigma(srt)
*
      do 1001 ie = 1,19
        if (earray(ie) .eq. plab) then
          TWOPI= xarray(ie)
          return
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          TWOPI= exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &    /(xmax-xmin) )
          return
        end if
 1001   continue
      return
        END
******************************************
******************************************
* for ppi(+)-->DELTA+RHO
c      real*4 function THREPI(srt)
