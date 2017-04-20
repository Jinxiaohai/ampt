        REAL FUNCTION W1440(DMASS)
      SAVE   
        AVMASS=0.938868
        PIMASS=0.137265
           AUX = 0.25*(DMASS**2-AVMASS**2-PIMASS**2)**2
     &           -(AVMASS*PIMASS)**2
            IF (AUX .GT. 0.) THEN
              QAVAIL = SQRT(AUX)/DMASS
            ELSE
              QAVAIL = 1.E-06
            END IF
c              w1440=0.2 
           W1440 = 0.2* (QAVAIL/0.397)**3
        RETURN
        END
****************
* PURPOSE : CALCULATE THE PION(ETA)+NUCLEON CROSS SECTION 
*           ACCORDING TO THE BREIT-WIGNER FORMULA, 
*           NOTE THAT N*(1535) IS S_11
* VARIABLE : LA = 1 FOR PI+N
*            LA = 0 FOR ETA+N
* DATE    : MAY 16, 1994
****************
