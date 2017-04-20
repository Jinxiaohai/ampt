       Real function fkaon(p,pmax)
      SAVE   
       fmax=0.148
       if(pmax.eq.0.)pmax=0.000001
       fkaon=(1.-p/pmax)*(p/pmax)**2
       if(fkaon.gt.fmax)fkaon=fmax
       fkaon=fkaon/fmax
       return
       end
*************************
* cross section for N*(1535) production in ND OR NN* collisions
* VARIABLES:
* LB1,LB2 ARE THE LABLES OF THE TWO COLLIDING PARTICLES
* SRT IS THE CMS ENERGY
* X1535 IS THE N*(1535) PRODUCTION CROSS SECTION
* NOTE THAT THE N*(1535) PRODUCTION CROSS SECTION IS 2 TIMES THE ETA 
* PRODUCTION CROSS SECTION
* DATE: MAY 18, 1994
* ***********************
