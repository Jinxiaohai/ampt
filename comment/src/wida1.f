        subroutine WIDA1(DMASS,rhomp,wa1,iseed)
      SAVE   
c
        PIMASS=0.137265
        coupa = 14.8
c
       RHOMAX = DMASS-PIMASS-0.02
       IF(RHOMAX.LE.0)then
         rhomp=0.
c   !! no decay
         wa1=-10.
        endif
        icount = 0
711       rhomp=RHOMAS(RHOMAX,ISEED)
      icount=icount+1
      if(dmass.le.(pimass+rhomp)) then
       if(icount.le.100) then
        goto 711
       else
         rhomp=0.
c   !! no decay
         wa1=-10.
        return
       endif
      endif
      qqp2=(dmass**2-(rhomp+pimass)**2)*(dmass**2-(rhomp-pimass)**2)
      qqp=sqrt(qqp2)/(2.0*dmass)
      epi=sqrt(pimass**2+qqp**2)
      erho=sqrt(rhomp**2+qqp**2)
      epirho=2.0*(epi*erho+qqp**2)**2+rhomp**2*epi**2
      wa1=coupa**2*qqp*epirho/(24.0*3.1416*dmass**2)
       return
       end
************************************       
* FUNCTION W1535(DMASS) GIVES THE N*(1535) DECAY WIDTH 
c     FOR A GIVEN N*(1535) MASS
* HERE THE FORMULA GIVEN BY KITAZOE IS USED
