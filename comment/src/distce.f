            SUBROUTINE DISTCE(I1,I2,DELTAR,DS,DT,EC,SRT
     1      ,IC,PX1CM,PY1CM,PZ1CM)
* PURPOSE : CHECK IF THE COLLISION BETWEEN TWO PARTICLES CAN HAPPEN
*           BY CHECKING
*                      (1) IF THE DISTANCE BETWEEN THEM IS SMALLER
*           THAN THE MAXIMUM DISTANCE DETERMINED FROM THE CROSS SECTION.
*                      (2) IF PARTICLE WILL PASS EACH OTHER WITHIN
*           TWO HARD CORE RADIUS.
*                      (3) IF PARTICLES WILL GET CLOSER.
* VARIABLES :
*           IC=1 COLLISION HAPPENED
*           IC=-1 COLLISION CAN NOT HAPPEN
*****************************************
            PARAMETER (MAXSTR=150001)
            COMMON   /AA/  R(3,MAXSTR)
cc      SAVE /AA/
            COMMON   /BB/  P(3,MAXSTR)
cc      SAVE /BB/
            COMMON   /CC/  E(MAXSTR)
cc      SAVE /CC/
            COMMON   /BG/  BETAX,BETAY,BETAZ,GAMMA
            COMMON  /EE/      ID(MAXSTR),LB(MAXSTR)
cc      SAVE /BG/
            common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1           px1n,py1n,pz1n,dp1n
            common /dpi/em2,lb2
            SAVE   
            IC=0
            X1=R(1,I1)
            Y1=R(2,I1)
            Z1=R(3,I1)
            PX1=P(1,I1)
            PY1=P(2,I1)
            PZ1=P(3,I1)
            X2=R(1,I2)
            Y2=R(2,I2)
            Z2=R(3,I2)
            PX2=P(1,I2)
            PY2=P(2,I2)
            PZ2=P(3,I2)
            EM1=E(I1)
            EM2=E(I2)
            E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
c            IF (ABS(X1-X2) .GT. DELTAR) GO TO 400
c            IF (ABS(Y1-Y2) .GT. DELTAR) GO TO 400
c            IF (ABS(Z1-Z2) .GT. DELTAR) GO TO 400
            RSQARE = (X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2
            IF (RSQARE .GT. DELTAR**2) GO TO 400
*NOW PARTICLES ARE CLOSE ENOUGH TO EACH OTHER !
              E2     = SQRT ( EM2**2 + PX2**2 + PY2**2 + PZ2**2 )
              S      = SRT*SRT
            IF (S .LT. EC) GO TO 400
*NOW THERE IS ENOUGH ENERGY AVAILABLE !
*LORENTZ-TRANSFORMATION IN I1-I2-C.M. SYSTEM
* BETAX, BETAY, BETAZ AND GAMMA HAVE BEEN GIVEN IN THE SUBROUTINE CMS
*TRANSFORMATION OF MOMENTA (PX1CM = - PX2CM)
              P1BETA = PX1*BETAX + PY1*BETAY + PZ1 * BETAZ
              TRANSF = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) - E1 )
              PRCM   = SQRT (PX1CM**2 + PY1CM**2 + PZ1CM**2)
              IF (PRCM .LE. 0.00001) GO TO 400
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
              IF (BBB .GT. DS) GO TO 400
              RELVEL = PRCM * (1.0/E1 + 1.0/E2)
              DDD    = RELVEL * DT * 0.5
*WILL PARTICLES GET CLOSER ?
              IF (ABS(DDD) .LT. ABS(DZZ)) GO TO 400
              IC=1
              GO TO 500
400           IC=-1
500           CONTINUE
              RETURN
              END
****************************************
*                                                                      *
*                                                                      *
