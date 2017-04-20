      real function xrho(srt)
      SAVE   
*       xsection for pp-->pp+rho
* *****************************
       pmass=0.9383
       rmass=0.77
       trho=0.151
       xrho=0.000000001
       if(srt.le.2.67)return
       ESMIN=2.*0.9383+rmass-trho/2.
       ES=srt
* the cross section for tho0 production is
       xrho0=0.24*(es-esmin)/(1.4+(es-esmin)**2)
       xrho=3.*Xrho0
       return
       end
* *****************************
c       real*4 function omega(srt)
