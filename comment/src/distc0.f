            SUBROUTINE distc0(drmax,deltr0,DT,
     1     Ifirst,PX1CM,PY1CM,PZ1CM,
     2     x1,y1,z1,px1,py1,pz1,em1,x2,y2,z2,px2,py2,pz2,em2)
* PURPOSE : CHECK IF THE COLLISION BETWEEN TWO PARTICLES CAN HAPPEN
*           BY CHECKING
*                      (2) IF PARTICLE WILL PASS EACH OTHER WITHIN
*           TWO HARD CORE RADIUS.
*                      (3) IF PARTICLES WILL GET CLOSER.
* VARIABLES :
*           Ifirst=1 COLLISION may HAPPENED
*           Ifirst=-1 COLLISION CAN NOT HAPPEN
*****************************************
            COMMON   /BG/  BETAX,BETAY,BETAZ,GAMMA
cc      SAVE /BG/
      SAVE   
            Ifirst=-1
            E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
*NOW PARTICLES ARE CLOSE ENOUGH TO EACH OTHER !
            E2     = SQRT ( EM2**2 + PX2**2 + PY2**2 + PZ2**2 )
*NOW THERE IS ENOUGH ENERGY AVAILABLE !
*LORENTZ-TRANSFORMATION IN I1-I2-C.M. SYSTEM
* BETAX, BETAY, BETAZ AND GAMMA HAVE BEEN GIVEN IN THE SUBROUTINE CMS
*TRANSFORMATION OF MOMENTA (PX1CM = - PX2CM)
              P1BETA = PX1*BETAX + PY1*BETAY + PZ1 * BETAZ
              TRANSF = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) - E1 )
              PRCM   = SQRT (PX1CM**2 + PY1CM**2 + PZ1CM**2)
              IF (PRCM .LE. 0.00001) return
*TRANSFORMATION OF SPATIAL DISTANCE
              DRBETA = BETAX*(X1-X2) + BETAY*(Y1-Y2) + BETAZ*(Z1-Z2)
              TRANSF = GAMMA * GAMMA * DRBETA / (GAMMA + 1)
              DXCM   = BETAX * TRANSF + X1 - X2
              DYCM   = BETAY * TRANSF + Y1 - Y2
              DZCM   = BETAZ * TRANSF + Z1 - Z2
*DETERMINING IF THIS IS THE POINT OF CLOSEST APPROACH
              DRCM   = SQRT (DXCM**2  + DYCM**2  + DZCM**2 )
              DZZ    = (PX1CM*DXCM + PY1CM*DYCM + PZ1CM*DZCM) / PRCM
              if ((drcm**2 - dzz**2) .le. 0.) then
                BBB = 0.
              else
                BBB    = SQRT (DRCM**2 - DZZ**2)
              end if
*WILL PARTICLE PASS EACH OTHER WITHIN 2 * HARD CORE RADIUS ?
              IF (BBB .GT. drmax) return
              RELVEL = PRCM * (1.0/E1 + 1.0/E2)
              DDD    = RELVEL * DT * 0.5
*WILL PARTICLES GET CLOSER ?
              IF (ABS(DDD) .LT. ABS(DZZ)) return
              Ifirst=1
              RETURN
              END
*---------------------------------------------------------------------------
c
clin-8/2008 B+B->Deuteron+Meson cross section in mb:
