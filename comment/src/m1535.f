       Subroutine M1535(LB1,LB2,SRT,X1535)
      SAVE   
       S0=2.424
       x1535=0.
       IF(SRT.LE.S0)RETURN
       SIGMA=2.*0.102*(SRT-S0)/(0.058+(SRT-S0)**2)
* I N*(1535) PRODUCTION IN NUCLEON-DELTA COLLISIONS
*(1) nD(++)->pN*(+)(1535), pD(-)->nN*(0)(1535),pD(+)-->N*(+)p
cbz11/25/98
c       IF((LB1*LB2.EQ.18).OR.(LB1*LB2.EQ.6).
c     1  or.(lb1*lb2).eq.8)then
       IF((LB1*LB2.EQ.18.AND.(LB1.EQ.2.OR.LB2.EQ.2)).OR.
     &     (LB1*LB2.EQ.6.AND.(LB1.EQ.1.OR.LB2.EQ.1)).or.
     &     (lb1*lb2.eq.8.AND.(LB1.EQ.1.OR.LB2.EQ.1)))then
cbz11/25/98end
       X1535=SIGMA
       return
       ENDIF
*(2) pD(0)->pN*(0)(1535),pD(0)->nN*(+)(1535)
       IF(LB1*LB2.EQ.7)THEN
       X1535=3.*SIGMA
       RETURN
       ENDIF 
* II N*(1535) PRODUCTION IN N*(1440)+NUCLEON REACTIONS
*(3) N*(+)(1440)p->N*(0+)(1535)p, N*(0)(1440)n->N*(0)(1535)
cbz11/25/98
c       IF((LB1*LB2.EQ.11).OR.(LB1*LB2.EQ.20))THEN
       IF((LB1*LB2.EQ.11).OR.
     &     (LB1*LB2.EQ.20.AND.(LB1.EQ.2.OR.LB2.EQ.2)))THEN
cbz11/25/98end
       X1535=SIGMA
       RETURN
       ENDIF
*(4) N*(0)(1440)p->N*(0+) or N*(+)(1440)n->N*(0+)(1535)
cbz11/25/98
c       IF((LB1*LB2.EQ.10).OR.(LB1*LB2.EQ.22))X1535=3.*SIGMA
       IF((LB1*LB2.EQ.10.AND.(LB1.EQ.1.OR.LB2.EQ.1)).OR.
     &     (LB1*LB2.EQ.22.AND.(LB1.EQ.2.OR.LB2.EQ.2)))
     &     X1535=3.*SIGMA
cbz11/25/98end
       RETURN
       END
*************************
* cross section for N*(1535) production in NN collisions
* VARIABLES:
* LB1,LB2 ARE THE LABLES OF THE TWO COLLIDING PARTICLES
* SRT IS THE CMS ENERGY
* X1535 IS THE N*(1535) PRODUCTION CROSS SECTION
* NOTE THAT THE N*(1535) PRODUCTION CROSS SECTION IS 2 TIMES THE ETA 
* PRODUCTION CROSS SECTION
* DATE: MAY 18, 1994
* ***********************
