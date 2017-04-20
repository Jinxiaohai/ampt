      SUBROUTINE Crppba(PX,PY,PZ,SRT,I1,I2,IBLOCK)
*     PURPOSE:                                                         *
clin-8/29/00*             DEALING WITH anti-nucleon annihilation with 
*             DEALING WITH anti-baryon annihilation with 
*             nucleons or baryon resonances
*             Determine:                                               *
*             (1) no. of pions in the final state
*             (2) relable particles in the final state
*             (3) new momenta of final state particles                 *
*                  
*     QUANTITIES:                                                      *
*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
*           SRT      - SQRT OF S                                       *
*           IBLOCK   - INFORMATION about the reaction channel          *
*                
*           iblock   - 1902 annihilation-->pion(+)+pion(-)   (2 pion)
*           iblock   - 1903 annihilation-->pion(+)+rho(-)    (3 pion)
*           iblock   - 1904 annihilation-->rho(+)+rho(-)     (4 pion)
*           iblock   - 1905 annihilation-->rho(0)+omega      (5 pion)
*           iblock   - 1906 annihilation-->omega+omega       (6 pion)
*       charge conservation is enforced in relabling particles 
*       in the final state (note: at the momentum we don't check the
*       initial charges while dealing with annihilation, since some
*       annihilation channels between antinucleons and nucleons (baryon
*       resonances) might be forbiden by charge conservation, this effect
*       should be small, but keep it in mind.
**********************************
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AMRHO=0.769,AMOMGA=0.782,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974)
        parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
        COMMON /AA/ R(3,MAXSTR)
cc      SAVE /AA/
        COMMON /BB/ P(3,MAXSTR)
cc      SAVE /BB/
        COMMON /CC/ E(MAXSTR)
cc      SAVE /CC/
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
       PX0=PX
       PY0=PY
       PZ0=PZ
* determine the no. of pions in the final state using a 
* statistical model
       call pbarfs(srt,npion,iseed)
* find the masses of the final state particles before calculate 
* their momenta, and relable them. The masses of rho and omega 
* will be generated according to the Breit Wigner formula       (NOTE!!!
* NOT DONE YET, AT THE MOMENT LET US USE FIXED RHO AND OMEGA MAEES)
cbali2/22/99
* Here we generate two stes of integer random numbers (3,4,5)
* one or both of them are used directly as the lables of pions
* similarly, 22+nchrg1 and 22+nchrg2 are used directly 
* to label rhos  
       nchrg1=3+int(3*RANART(NSEED))
       nchrg2=3+int(3*RANART(NSEED))
* the corresponding masses of pions
      pmass1=ap1
       pmass2=ap1
       if(nchrg1.eq.3.or.nchrg1.eq.5)pmass1=ap2
       if(nchrg2.eq.3.or.nchrg2.eq.5)pmass2=ap2
* (1) for 2 pion production
       IF(NPION.EQ.2)THEN 
       IBLOCK=1902
* randomly generate the charges of final state particles,
       LB(I1)=nchrg1
       E(I1)=pmass1
       LB(I2)=nchrg2
       E(I2)=pmass2
* TO CALCULATE THE FINAL MOMENTA
       GO TO 50
       ENDIF
* (2) FOR 3 PION PRODUCTION
       IF(NPION.EQ.3)THEN 
       IBLOCK=1903
       LB(I1)=nchrg1
       E(I1)=pmass1
       LB(I2)=22+nchrg2
            E(I2)=AMRHO
       GO TO 50
       ENDIF
* (3) FOR 4 PION PRODUCTION
* we allow both rho+rho and pi+omega with 50-50% probability
        IF(NPION.EQ.4)THEN 
       IBLOCK=1904
* determine rho+rho or pi+omega
       if(RANART(NSEED).ge.0.5)then
* rho+rho  
       LB(I1)=22+nchrg1
       E(I1)=AMRHO
       LB(I2)=22+nchrg2
            E(I2)=AMRHO
       else
* pion+omega
       LB(I1)=nchrg1
       E(I1)=pmass1
       LB(I2)=28
            E(I2)=AMOMGA
       endif
       GO TO 50
       ENDIF
* (4) FOR 5 PION PRODUCTION
        IF(NPION.EQ.5)THEN 
       IBLOCK=1905
* RHO AND OMEGA
        LB(I1)=22+nchrg1
       E(I1)=AMRHO
       LB(I2)=28
       E(I2)=AMOMGA
       GO TO 50
       ENDIF
* (5) FOR 6 PION PRODUCTION
         IF(NPION.EQ.6)THEN 
       IBLOCK=1906
* OMEGA AND OMEGA
        LB(I1)=28
       E(I1)=AMOMGA
       LB(I2)=28
          E(I2)=AMOMGA
       ENDIF
cbali2/22/99
50    EM1=E(I1)
      EM2=E(I2)
*-----------------------------------------------------------------------
* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
* ENERGY CONSERVATION
          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.E-08
          PR=SQRT(PR2)/(2.*SRT)
* WE ASSUME AN ISOTROPIC ANGULAR DISTRIBUTION IN THE CMS 
          C1   = 1.0 - 2.0 * RANART(NSEED)
          T1   = 2.0 * PI * RANART(NSEED)
      S1   = SQRT( 1.0 - C1**2 )
      CT1  = COS(T1)
      ST1  = SIN(T1)
* THE MOMENTUM IN THE CMS IN THE FINAL STATE
      PZ   = PR * C1
      PX   = PR * S1*CT1 
      PY   = PR * S1*ST1
* ROTATE IT 
       CALL ROTATE(PX0,PY0,PZ0,PX,PY,PZ) 
      RETURN
      END
cbali2/7/99end
cbali3/5/99
**********************************
*     PURPOSE:                                                         *
*     assign final states for K+K- --> light mesons
*
