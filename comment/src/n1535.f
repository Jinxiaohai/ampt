       Subroutine N1535(LB1,LB2,SRT,X1535)
      SAVE   
       S0=2.424
       x1535=0.
       IF(SRT.LE.S0)RETURN
       SIGMA=2.*0.102*(SRT-S0)/(0.058+(SRT-S0)**2)
* I N*(1535) PRODUCTION IN NUCLEON-NUCLEON COLLISIONS
*(1) pp->pN*(+)(1535), nn->nN*(0)(1535)
cbdbg11/25/98
c       IF((LB1*LB2.EQ.1).OR.(LB1*LB2.EQ.4))then
       IF((LB1*LB2.EQ.1).OR.
     &     (LB1.EQ.2.AND.LB2.EQ.2))then
cbz11/25/98end
       X1535=SIGMA
       return
       endif
*(2) pn->pN*(0)(1535),pn->nN*(+)(1535)
       IF(LB1*LB2.EQ.2)then
       X1535=3.*SIGMA
       return
       endif 
* III N*(1535) PRODUCTION IN DELTA+DELTA REACTIONS
* (5) D(++)+D(0), D(+)+D(+),D(+)+D(-),D(0)+D(0)
cbz11/25/98
c       IF((LB1*LB2.EQ.63).OR.(LB1*LB2.EQ.64).OR.(LB1*LB2.EQ.48).
c     1  OR.(LB1*LB2.EQ.49))then
       IF((LB1*LB2.EQ.63.AND.(LB1.EQ.7.OR.LB2.EQ.7)).OR.
     &     (LB1*LB2.EQ.64.AND.(LB1.EQ.8.OR.LB2.EQ.8)).OR.
     &     (LB1*LB2.EQ.48.AND.(LB1.EQ.6.OR.LB2.EQ.6)).OR.
     &     (LB1*LB2.EQ.49.AND.(LB1.EQ.7.OR.LB2.EQ.7)))then
cbz11/25/98end
       X1535=SIGMA
       return
       endif
* (6) D(++)+D(-),D(+)+D(0)
cbz11/25/98
c       IF((LB1*LB2.EQ.54).OR.(LB1*LB2.EQ.56))then
       IF((LB1*LB2.EQ.54.AND.(LB1.EQ.6.OR.LB2.EQ.6)).OR.
     &     (LB1*LB2.EQ.56.AND.(LB1.EQ.7.OR.LB2.EQ.7)))then
cbz11/25/98end
       X1535=3.*SIGMA
       return
       endif
* IV N*(1535) PRODUCTION IN N*(1440)+N*(1440) REACTIONS
cbz11/25/98
c       IF((LB1*LB2.EQ.100).OR.(LB1*LB2.EQ.11*11))X1535=SIGMA
       IF((LB1.EQ.10.AND.LB2.EQ.10).OR.
     &     (LB1.EQ.11.AND.LB2.EQ.11))X1535=SIGMA
c       IF(LB1*LB2.EQ.110)X1535=3.*SIGMA
       IF(LB1*LB2.EQ.110.AND.(LB1.EQ.10.OR.LB2.EQ.10))X1535=3.*SIGMA
cbdbg11/25/98end
       RETURN
       END
************************************       
* FUNCTION WA1(DMASS) GIVES THE A1 DECAY WIDTH
