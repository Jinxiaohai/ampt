        REAL FUNCTION FD5(DMASS,SRT,CON)
      SAVE   
        AMN=0.938869
        AVPI=0.13803333
        AM0=1.535
        FD=4.*(AM0**2)*W1535(DMASS)/((DMASS**2-1.535**2)**2
     1  +AM0**2*W1535(DMASS)**2)
        IF(CON.EQ.1.)THEN
clin-9/2012: check argument in sqrt():
           scheck=(SRT**2+DMASS**2-AMN**2)**2/(4.*SRT**2)-DMASS**2
           if(scheck.lt.0) then
              write(99,*) 'scheck11: ', scheck
              scheck=0.
           endif
           P1=SQRT(scheck)
c           P1=SQRT((SRT**2+DMASS**2-AMN**2)**2
c     1          /(4.*SRT**2)-DMASS**2)
        ELSE
        DMASS=AMN+AVPI
clin-9/2012: check argument in sqrt():
        scheck=(SRT**2+DMASS**2-AMN**2)**2/(4.*SRT**2)-DMASS**2
        if(scheck.lt.0) then
           write(99,*) 'scheck12: ', scheck
           scheck=0.
        endif
        P1=SQRT(scheck)
c        P1=SQRT((SRT**2+DMASS**2-AMN**2)**2
c     1  /(4.*SRT**2)-DMASS**2)
        ENDIF
        FD5=FD*P1*DMASS
        RETURN
        END
*--------------------------------------------------------------------------
*FUNCTION FNS(DMASS) GIVES N* MASS DISTRIBUTION 
c     BY USING OF BREIT-WIGNER FORMULA
