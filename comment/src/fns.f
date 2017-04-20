        REAL FUNCTION FNS(DMASS,SRT,CON)
      SAVE   
        WIDTH=0.2
        AMN=0.938869
        AVPI=0.13803333
        AN0=1.43
        FN=4.*(AN0**2)*WIDTH/((DMASS**2-1.44**2)**2+AN0**2*WIDTH**2)
        IF(CON.EQ.1.)THEN
clin-9/2012: check argument in sqrt():
           scheck=(SRT**2+DMASS**2-AMN**2)**2/(4.*SRT**2)-DMASS**2
           if(scheck.lt.0) then
              write(99,*) 'scheck13: ', scheck
              scheck=0.
           endif
           P1=SQRT(scheck)
c        P1=SQRT((SRT**2+DMASS**2-AMN**2)**2
c     1  /(4.*SRT**2)-DMASS**2)
        ELSE
        DMASS=AMN+AVPI
clin-9/2012: check argument in sqrt():
        scheck=(SRT**2+DMASS**2-AMN**2)**2/(4.*SRT**2)-DMASS**2
        if(scheck.lt.0) then
           write(99,*) 'scheck14: ', scheck
           scheck=0.
        endif
        P1=SQRT(scheck)
c        P1=SQRT((SRT**2+DMASS**2-AMN**2)**2
c     1  /(4.*SRT**2)-DMASS**2)
        ENDIF
        FNS=FN*P1*DMASS
        RETURN
        END
*-----------------------------------------------------------------------------
*-----------------------------------------------------------------------------
* PURPOSE:1. SORT N*(1440) and N*(1535) 2-body DECAY PRODUCTS
*         2. DETERMINE THE MOMENTUM AND COORDINATES OF NUCLEON AND PION
*            AFTER THE DELTA OR N* DECAYING
* DATE   : JAN. 24,1990, MODIFIED ON MAY 17, 1994 TO INCLUDE ETA 
