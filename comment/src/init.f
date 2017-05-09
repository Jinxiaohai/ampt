      SUBROUTINE INIT(MINNUM,MAXNUM,NUM,RADIUS,X0,Z0,P0,
     &                GAMMA,ISEED,MASS,IOPT)
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      目的:providing initial conditions for phase-space distribution
c$$$      of testparticles
c$$$      variables:   (all input)
c$$$      minnum   :   first testparticle treated in one run
c$$$      maxnum   :   last  testparticle treated in one run
c$$$      num      :   number of testparticles per nucleon
c$$$      radius   :   radius of nucleus "FM"
c$$$      x0, y0   :   displacement of center ofnucleus in x,z direction
c$$$                   "fm"
c$$$      p0       :   momentum-boost in cm frame
c$$$      gamma    :   relativistic gamma-factor
c$$$      iseed    :   seed for random-number generator
c$$$      mass     :   total mass of the system
c$$$      iopt     :   option for different occupation of momentum
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      PARAMETER     (MAXSTR=150001,  AMU   = 0.9383)
      PARAMETER     (MAXX   =   20,  MAXZ  =    24)
      PARAMETER     (PI=3.1415926)
*
      REAL              PTOT(3)
      COMMON  /AA/      R(3,MAXSTR)
cc      SAVE /AA/
      COMMON  /BB/      P(3,MAXSTR)
cc      SAVE /BB/
      COMMON  /CC/      E(MAXSTR)
cc      SAVE /CC/
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
cc      SAVE /DD/
      COMMON  /EE/      ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      common  /ss/      inout(20)
cc      SAVE /ss/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
*----------------------------------------------------------------------
*     PREPARATION FOR LORENTZ-TRANSFORMATIONS
*
      IF (P0 .NE. 0.) THEN
        SIGN = P0 / ABS(P0)
      ELSE
        SIGN = 0.
      END IF
clin-9/2012: check argument in sqrt():
      scheck=GAMMA**2-1.
      if(scheck.lt.0) then
         write(99,*) 'scheck10: ', scheck
         scheck=0.
      endif
      BETA=SIGN*SQRT(scheck)/GAMMA
c      BETA = SIGN * SQRT(GAMMA**2-1.)/GAMMA
*-----------------------------------------------------------------------
*     TARGET-ID = 1 AND PROJECTILE-ID = -1
*
      IF (MINNUM .EQ. 1) THEN
        IDNUM = 1
      ELSE
        IDNUM = -1
      END IF
*-----------------------------------------------------------------------
*     IDENTIFICATION OF TESTPARTICLES AND ASSIGMENT OF RESTMASS
*
*     LOOP OVER ALL PARALLEL RUNS:
      DO 400 IRUN = 1,NUM
        DO 100 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
          ID(I) = IDNUM
          E(I)  = AMU
  100   CONTINUE
*-----------------------------------------------------------------------
*       OCCUPATION OF COORDINATE-SPACE
*
        DO 300 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
  200     CONTINUE
            X = 1.0 - 2.0 * RANART(NSEED)
            Y = 1.0 - 2.0 * RANART(NSEED)
            Z = 1.0 - 2.0 * RANART(NSEED)
          IF ((X*X+Y*Y+Z*Z) .GT. 1.0) GOTO 200
          R(1,I) = X * RADIUS
          R(2,I) = Y * RADIUS
          R(3,I) = Z * RADIUS
  300   CONTINUE
  400 CONTINUE
*=======================================================================
      IF (IOPT .NE. 3) THEN
*-----
*     OPTION 1: USE WOODS-SAXON PARAMETRIZATION FOR DENSITY AND
*-----          CALCULATE LOCAL FERMI-MOMENTUM
*
        RHOW0  = 0.168
        DO 1000 IRUN = 1,NUM
          DO 600 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
  500       CONTINUE
              PX = 1.0 - 2.0 * RANART(NSEED)
              PY = 1.0 - 2.0 * RANART(NSEED)
              PZ = 1.0 - 2.0 * RANART(NSEED)
            IF (PX*PX+PY*PY+PZ*PZ .GT. 1.0) GOTO 500
            RDIST  = SQRT( R(1,I)**2 + R(2,I)**2 + R(3,I)**2 )
            RHOWS  = RHOW0 / (  1.0 + EXP( (RDIST-RADIUS) / 0.55 )  )
            PFERMI = 0.197 * (1.5 * PI*PI * RHOWS)**(1./3.)
*-----
*     OPTION 2: NUCLEAR MATTER CASE
            IF(IOPT.EQ.2) PFERMI=0.27
           if(iopt.eq.4) pfermi=0.
*-----
            P(1,I) = PFERMI * PX
            P(2,I) = PFERMI * PY
            P(3,I) = PFERMI * PZ
  600     CONTINUE
*
*         SET TOTAL MOMENTUM TO 0 IN REST FRAME AND BOOST
*
          DO 700 IDIR = 1,3
            PTOT(IDIR) = 0.0
  700     CONTINUE
          NPART = 0
          DO 900 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
            NPART = NPART + 1
            DO 800 IDIR = 1,3
              PTOT(IDIR) = PTOT(IDIR) + P(IDIR,I)
  800       CONTINUE
  900     CONTINUE
          DO 950 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
            DO 925 IDIR = 1,3
              P(IDIR,I) = P(IDIR,I) - PTOT(IDIR) / FLOAT(NPART)
  925       CONTINUE
*           BOOST
            IF ((IOPT .EQ. 1).or.(iopt.eq.2)) THEN
              EPART = SQRT(P(1,I)**2+P(2,I)**2+P(3,I)**2+AMU**2)
              P(3,I) = GAMMA*(P(3,I) + BETA*EPART)
            ELSE
              P(3,I) = P(3,I) + P0
            END IF
  950     CONTINUE
 1000   CONTINUE
*-----
      ELSE
*-----
*     OPTION 3: GIVE ALL NUCLEONS JUST A Z-MOMENTUM ACCORDING TO
*               THE BOOST OF THE NUCLEI
*
        DO 1200 IRUN = 1,NUM
          DO 1100 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
            P(1,I) = 0.0
            P(2,I) = 0.0
            P(3,I) = P0
 1100     CONTINUE
 1200   CONTINUE
*-----
      END IF
*=======================================================================
*     PUT PARTICLES IN THEIR POSITION IN COORDINATE-SPACE
*     (SHIFT AND RELATIVISTIC CONTRACTION)
*
      DO 1400 IRUN = 1,NUM
        DO 1300 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
          R(1,I) = R(1,I) + X0
* two nuclei in touch after contraction
          R(3,I) = (R(3,I)+Z0)/ GAMMA 
* two nuclei in touch before contraction
c          R(3,I) = R(3,I) / GAMMA + Z0
 1300   CONTINUE
 1400 CONTINUE
*
      RETURN
      END
**********************************
*                                                                      *
