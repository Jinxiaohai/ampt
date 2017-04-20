      SUBROUTINE Crdir(PX,PY,PZ,SRT,I1,I2,IBLOCK)
*     PURPOSE:                                                         *
*             DEALING WITH pion+N-->pion+N PROCESS                   *
*     NOTE   :                                                         *
*          
*     QUANTITIES:                                                 *
*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
*           SRT      - SQRT OF S                                       *
*           IBLOCK   - THE INFORMATION BACK                            *
*                    
**********************************
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,
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
        IBLOCK=999
        NTAG=0
        EM1=E(I1)
        EM2=E(I2)
*-----------------------------------------------------------------------
* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
* ENERGY CONSERVATION
        PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.e-09
          PR=SQRT(PR2)/(2.*SRT)
clin-10/25/02 get rid of argument usage mismatch in PTR():
          xptr=0.33*pr
c         cc1=ptr(0.33*pr,iseed)
         cc1=ptr(xptr,iseed)
clin-10/25/02-end
clin-9/2012: check argument in sqrt():
         scheck=pr**2-cc1**2
         if(scheck.lt.0) then
            write(99,*) 'scheck37: ', scheck
            scheck=0.
         endif
         c1=sqrt(scheck)/pr
c         c1=sqrt(pr**2-cc1**2)/pr
           T1   = 2.0 * PI * RANART(NSEED)
      S1   = SQRT( 1.0 - C1**2 )
      CT1  = COS(T1)
      ST1  = SIN(T1)
* THE MOMENTUM IN THE CMS IN THE FINAL STATE
      PZ   = PR * C1
      PX   = PR * S1*CT1 
      PY   = PR * S1*ST1
* ROTATE the momentum
      call rotate(px0,py0,pz0,px,py,pz)
      RETURN
      END
**********************************
*                                                                      *
*                                                                      *
