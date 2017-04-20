       real function ptdis(x)
      SAVE   
* NUCLEON TRANSVERSE MOMENTUM DISTRIBUTION AT HIGH ENERGIES
* DATE: Aug. 11, 1994
*********************************
       b=3.78
       c=0.47
       d=3.60
c       b=b*3
c       d=d*3
       ptdis=1./(2.*b)*(1.-exp(-b*x**2))-c/d*x*exp(-d*x)
     1     -c/D**2*(exp(-d*x)-1.)
       return
       end
*****************************
