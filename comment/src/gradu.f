      SUBROUTINE GRADU(IOPT,IX,IY,IZ,GRADX,GRADY,GRADZ)
*                                                                      *
*       PURPOSE:     DETERMINE GRAD(U(RHO(X,Y,Z)))                     *
*       VARIABLES:                                                     *
*         IOPT                - METHOD FOR EVALUATING THE GRADIENT     *
*                                                      (INTEGER,INPUT) *
*         IX, IY, IZ          - COORDINATES OF POINT   (INTEGER,INPUT) *
*         GRADX, GRADY, GRADZ - GRADIENT OF U            (REAL,OUTPUT) *
*                                                                      *
**********************************
      PARAMETER         (MAXX =    20,  MAXZ =   24)
      PARAMETER         (RHO0 = 0.167)
*
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                  RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                  RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
cc      SAVE /DD/
      common  /ss/      inout(20)
cc      SAVE /ss/
      common  /tt/  PEL(-maxx:maxx,-maxx:maxx,-maxz:maxz)
     &,rxy(-maxx:maxx,-maxx:maxx,-maxz:maxz)
cc      SAVE /tt/
      SAVE   
*
      RXPLUS   = RHO(IX+1,IY,  IZ  ) / RHO0
      RXMINS   = RHO(IX-1,IY,  IZ  ) / RHO0
      RYPLUS   = RHO(IX,  IY+1,IZ  ) / RHO0
      RYMINS   = RHO(IX,  IY-1,IZ  ) / RHO0
      RZPLUS   = RHO(IX,  IY,  IZ+1) / RHO0
      RZMINS   = RHO(IX,  IY,  IZ-1) / RHO0
      den0     = RHO(IX,  IY,  IZ) / RHO0
      ene0     = pel(IX,  IY,  IZ) 
*-----------------------------------------------------------------------
      GOTO (1,2,3,4,5) IOPT
       if(iopt.eq.6)go to 6
       if(iopt.eq.7)go to 7
*
    1 CONTINUE
*       POTENTIAL USED IN 1) (STIFF):
*       U = -.124 * RHO/RHO0 + .0705 (RHO/RHO0)**2 GEV
*
           GRADX  = -0.062 * (RXPLUS - RXMINS) + 0.03525 * (RXPLUS**2 -
     &                                                      RXMINS**2)
           GRADY  = -0.062 * (RYPLUS - RYMINS) + 0.03525 * (RYPLUS**2 -
     &                                                      RYMINS**2)
           GRADZ  = -0.062 * (RZPLUS - RZMINS) + 0.03525 * (RZPLUS**2 -
     &                                                      RZMINS**2)
           RETURN
*
    2 CONTINUE
*       POTENTIAL USED IN 2):
*       U = -.218 * RHO/RHO0 + .164 (RHO/RHO0)**(4/3) GEV
*
           EXPNT = 1.3333333
           GRADX = -0.109 * (RXPLUS - RXMINS) 
     &     + 0.082 * (RXPLUS**EXPNT-RXMINS**EXPNT)
           GRADY = -0.109 * (RYPLUS - RYMINS) 
     &     + 0.082 * (RYPLUS**EXPNT-RYMINS**EXPNT)
           GRADZ = -0.109 * (RZPLUS - RZMINS) 
     &     + 0.082 * (RZPLUS**EXPNT-RZMINS**EXPNT)
           RETURN
*
    3 CONTINUE
*       POTENTIAL USED IN 3) (SOFT):
*       U = -.356 * RHO/RHO0 + .303 * (RHO/RHO0)**(7/6)  GEV
*
           EXPNT = 1.1666667
          acoef = 0.178
           GRADX = -acoef * (RXPLUS - RXMINS) 
     &     + 0.1515 * (RXPLUS**EXPNT-RXMINS**EXPNT)
           GRADY = -acoef * (RYPLUS - RYMINS) 
     &     + 0.1515 * (RYPLUS**EXPNT-RYMINS**EXPNT)
           GRADZ = -acoef * (RZPLUS - RZMINS) 
     &     + 0.1515 * (RZPLUS**EXPNT-RZMINS**EXPNT)
                 RETURN
*
*
    4   CONTINUE
*       POTENTIAL USED IN 4) (super-soft in the mixed phase of 4 < rho/rho <7):
*       U1 = -.356 * RHO/RHO0 + .303 * (RHO/RHO0)**(7/6)  GEV
*       normal phase, soft eos of iopt=3
*       U2 = -.02 * (RHO/RHO0)**(2/3) -0.0253 * (RHO/RHO0)**(7/6)  GEV
*
       eh=4.
       eqgp=7.
           acoef=0.178
           EXPNT = 1.1666667
       denr=rho(ix,iy,iz)/rho0
       if(denr.le.eh.or.denr.ge.eqgp)then
           GRADX = -acoef * (RXPLUS - RXMINS) 
     &     + 0.1515 * (RXPLUS**EXPNT-RXMINS**EXPNT)
           GRADY = -acoef * (RYPLUS - RYMINS) 
     &     + 0.1515 * (RYPLUS**EXPNT-RYMINS**EXPNT)
           GRADZ = -acoef * (RZPLUS - RZMINS) 
     &     + 0.1515 * (RZPLUS**EXPNT-RZMINS**EXPNT)
       else
          acoef1=0.178
          acoef2=0.0
          expnt2=2./3.
           GRADX =-acoef1* (RXPLUS**EXPNT-RXMINS**EXPNT)
     &                 -acoef2* (RXPLUS**expnt2 - RXMINS**expnt2) 
           GRADy =-acoef1* (RyPLUS**EXPNT-RyMINS**EXPNT)
     &                 -acoef2* (RyPLUS**expnt2 - RyMINS**expnt2) 
           GRADz =-acoef1* (RzPLUS**EXPNT-RzMINS**EXPNT)
     &                 -acoef2* (RzPLUS**expnt2 - RzMINS**expnt2) 
       endif
       return
*     
    5   CONTINUE
*       POTENTIAL USED IN 5) (SUPER STIFF):
*       U = -.10322 * RHO/RHO0 + .04956 * (RHO/RHO0)**(2.77)  GEV
*
           EXPNT = 2.77
           GRADX = -0.0516 * (RXPLUS - RXMINS) 
     &     + 0.02498 * (RXPLUS**EXPNT-RXMINS**EXPNT)
           GRADY = -0.0516 * (RYPLUS - RYMINS) 
     &     + 0.02498 * (RYPLUS**EXPNT-RYMINS**EXPNT)
           GRADZ = -0.0516 * (RZPLUS - RZMINS) 
     &     + 0.02498 * (RZPLUS**EXPNT-RZMINS**EXPNT)
           RETURN
*
    6 CONTINUE
*       POTENTIAL USED IN 6) (STIFF-qgp):
*       U = -.124 * RHO/RHO0 + .0705 (RHO/RHO0)**2 GEV
*
       if(ene0.le.0.5)then
           GRADX  = -0.062 * (RXPLUS - RXMINS) + 0.03525 * (RXPLUS**2 -
     &                                                      RXMINS**2)
           GRADY  = -0.062 * (RYPLUS - RYMINS) + 0.03525 * (RYPLUS**2 -
     &                                                      RYMINS**2)
           GRADZ  = -0.062 * (RZPLUS - RZMINS) + 0.03525 * (RZPLUS**2 -
     &                                                      RZMINS**2)
           RETURN
       endif
       if(ene0.gt.0.5.and.ene0.le.1.5)then
*       U=c1-ef*rho/rho0**2/3
       ef=36./1000.
           GRADX  = -0.5*ef* (RXPLUS**0.67-RXMINS**0.67)
           GRADy  = -0.5*ef* (RyPLUS**0.67-RyMINS**0.67)
           GRADz  = -0.5*ef* (RzPLUS**0.67-RzMINS**0.67)
           RETURN
       endif
       if(ene0.gt.1.5)then
* U=800*(rho/rho0)**1/3.-Ef*(rho/rho0)**2/3.-c2
       ef=36./1000.
       cf0=0.8
        GRADX  =0.5*cf0*(rxplus**0.333-rxmins**0.333) 
     &         -0.5*ef* (RXPLUS**0.67-RXMINS**0.67)
        GRADy  =0.5*cf0*(ryplus**0.333-rymins**0.333) 
     &         -0.5*ef* (RyPLUS**0.67-RyMINS**0.67)
        GRADz  =0.5*cf0*(rzplus**0.333-rzmins**0.333) 
     &         -0.5*ef* (RzPLUS**0.67-RzMINS**0.67)
           RETURN
       endif
*
    7 CONTINUE
*       POTENTIAL USED IN 7) (Soft-qgp):
       if(den0.le.4.5)then
*       POTENTIAL USED is the same as IN 3) (SOFT):
*       U = -.356 * RHO/RHO0 + .303 * (RHO/RHO0)**(7/6)  GEV
*
           EXPNT = 1.1666667
          acoef = 0.178
           GRADX = -acoef * (RXPLUS - RXMINS) 
     &     + 0.1515 * (RXPLUS**EXPNT-RXMINS**EXPNT)
           GRADY = -acoef * (RYPLUS - RYMINS) 
     &     + 0.1515 * (RYPLUS**EXPNT-RYMINS**EXPNT)
           GRADZ = -acoef * (RZPLUS - RZMINS) 
     &     + 0.1515 * (RZPLUS**EXPNT-RZMINS**EXPNT)
       return
       endif
       if(den0.gt.4.5.and.den0.le.5.1)then
*       U=c1-ef*rho/rho0**2/3
       ef=36./1000.
           GRADX  = -0.5*ef* (RXPLUS**0.67-RXMINS**0.67)
           GRADy  = -0.5*ef* (RyPLUS**0.67-RyMINS**0.67)
           GRADz  = -0.5*ef* (RzPLUS**0.67-RzMINS**0.67)
           RETURN
       endif
       if(den0.gt.5.1)then
* U=800*(rho/rho0)**1/3.-Ef*(rho/rho0)**2/3.-c2
       ef=36./1000.
       cf0=0.8
        GRADX  =0.5*cf0*(rxplus**0.333-rxmins**0.333) 
     &         -0.5*ef* (RXPLUS**0.67-RXMINS**0.67)
        GRADy  =0.5*cf0*(ryplus**0.333-rymins**0.333) 
     &         -0.5*ef* (RyPLUS**0.67-RyMINS**0.67)
        GRADz  =0.5*cf0*(rzplus**0.333-rzmins**0.333) 
     &         -0.5*ef* (RzPLUS**0.67-RzMINS**0.67)
           RETURN
       endif
        END
**********************************
*                                                                      *
