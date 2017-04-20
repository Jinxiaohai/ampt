      SUBROUTINE GRADUN(IX,IY,IZ,GRADXN,GRADYN,GRADZN)
*                                                                      *
*       PURPOSE:     DETERMINE THE GRADIENT OF THE NEUTRON DENSITY     *
*       VARIABLES:                                                     *
*                                                                           *
*         IX, IY, IZ          - COORDINATES OF POINT   (INTEGER,INPUT) *
*         GRADXN, GRADYN, GRADZN - GRADIENT OF THE NEUTRON             *
*                                  DENSITY(REAL,OUTPUT)                *
*                                                                      *
**********************************
      PARAMETER         (MAXX =    20,  MAXZ =   24)
      PARAMETER         (RHO0 = 0.168)
*
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
cc      SAVE /DD/
      common  /ss/      inout(20)
cc      SAVE /ss/
      SAVE   
*
      RXPLUS   = RHON(IX+1,IY,  IZ  ) / RHO0
      RXMINS   = RHON(IX-1,IY,  IZ  ) / RHO0
      RYPLUS   = RHON(IX,  IY+1,IZ  ) / RHO0
      RYMINS   = RHON(IX,  IY-1,IZ  ) / RHO0
      RZPLUS   = RHON(IX,  IY,  IZ+1) / RHO0
      RZMINS   = RHON(IX,  IY,  IZ-1) / RHO0
*-----------------------------------------------------------------------
*
           GRADXN  = (RXPLUS - RXMINS)/2. 
           GRADYN  = (RYPLUS - RYMINS)/2.
           GRADZN  = (RZPLUS - RZMINS)/2.
           RETURN
      END
*-----------------------------------------------------------------------------
*FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF
*KITAZOE'S FORMULA
