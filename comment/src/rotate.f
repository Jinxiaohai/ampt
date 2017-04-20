      SUBROUTINE rotate(PX0,PY0,PZ0,px,py,pz)
      SAVE   
* purpose: rotate the momentum of a particle in the CMS of p1+p2 such that 
* the x' y' and z' in the cms of p1+p2 is the same as the fixed x y and z
* quantities:
*            px0,py0 and pz0 are the cms momentum of the incoming colliding
*            particles
*            px, py and pz are the cms momentum of any one of the particles 
*            after the collision to be rotated
***************************************
* the momentum, polar and azimuthal angles of the incoming momentm
      PR0  = SQRT( PX0**2 + PY0**2 + PZ0**2 )
      IF(PR0.EQ.0)PR0=0.00000001
      C2  = PZ0 / PR0
      IF(PX0 .EQ. 0.0 .AND. PY0 .EQ. 0.0) THEN
        T2 = 0.0
      ELSE
        T2=ATAN2(PY0,PX0)
      END IF
clin-9/2012: check argument in sqrt():
      scheck=1.0 - C2**2
      if(scheck.lt.0) then
         write(99,*) 'scheck45: ', scheck
         scheck=0.
      endif
      S2=sqrt(scheck)
c      S2  =  SQRT( 1.0 - C2**2 )
      CT2  = COS(T2)
      ST2  = SIN(T2)
* the momentum, polar and azimuthal angles of the momentum to be rotated
      PR=SQRT(PX**2+PY**2+PZ**2)
      IF(PR.EQ.0)PR=0.0000001
      C1=PZ/PR
      IF(PX.EQ.0.AND.PY.EQ.0)THEN
      T1=0.
      ELSE
      T1=ATAN2(PY,PX)
      ENDIF
clin-9/2012: check argument in sqrt():
      scheck=1.0 - C1**2
      if(scheck.lt.0) then
         write(99,*) 'scheck46: ', scheck
         scheck=0.
      endif
      S1=sqrt(scheck)
c      S1   = SQRT( 1.0 - C1**2 )
      CT1  = COS(T1)
      ST1  = SIN(T1)
      SS   = C2 * S1 * CT1  +  S2 * C1
* THE MOMENTUM AFTER ROTATION
      PX   = PR * ( SS*CT2 - S1*ST1*ST2 )
      PY   = PR * ( SS*ST2 + S1*ST1*CT2 )
      PZ   = PR * ( C1*C2 - S1*S2*CT1 )
      RETURN
      END
******************************************
c      real*4 function Xpp(srt)
