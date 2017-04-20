      FUNCTION INVFLV(IART)
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
c.....anti-Delta-
      IF (IART .EQ. -6) THEN
         INVFLV = -1114
         RETURN
      END IF
c.....anti-Delta0
      IF (IART .EQ. -7) THEN
         INVFLV = -2114
         RETURN
      END IF
c.....anti-Delta+
      IF (IART .EQ. -8) THEN
         INVFLV = -2214
         RETURN
      END IF
c.....anti-Delta++
      IF (IART .EQ. -9) THEN
         INVFLV = -2224
         RETURN
      END IF
cbzdbg2/23/99
c.....anti-proton
      IF (IART .EQ. -1) THEN
         INVFLV = -2212
         RETURN
      END IF
c.....anti-neutron
      IF (IART .EQ. -2) THEN
         INVFLV = -2112
         RETURN
      END IF
cbzdbg2/23/99end
c.....eta
      IF (IART .EQ. 0) THEN
         INVFLV = 221
         RETURN
      END IF
c.....proton
      IF (IART .EQ. 1) THEN
         INVFLV = 2212
         RETURN
      END IF
c.....neutron
      IF (IART .EQ. 2) THEN
         INVFLV = 2112
         RETURN
      END IF
c.....pi-
      IF (IART .EQ. 3) THEN
         INVFLV = -211
         RETURN
      END IF
c.....pi0
      IF (IART .EQ. 4) THEN
         INVFLV = 111
         RETURN
      END IF
c.....pi+
      IF (IART .EQ. 5) THEN
         INVFLV = 211
         RETURN
      END IF
c.....Delta-
      IF (IART .EQ. 6) THEN
         INVFLV = 1114
         RETURN
      END IF
c.....Delta0
      IF (IART .EQ. 7) THEN
         INVFLV = 2114
         RETURN
      END IF
c.....Delta+
      IF (IART .EQ. 8) THEN
         INVFLV = 2214
         RETURN
      END IF
c.....Delta++
      IF (IART .EQ. 9) THEN
         INVFLV = 2224
         RETURN
      END IF
cc.....N*(1440), N*(1535) temporary entry
c      IF (IART .GE. 10 .AND. IART .LE.13) THEN
c         INVFLV = 0
c         RETURN
c      END IF
c.....Lambda
      IF (IART .EQ. 14) THEN
         INVFLV = 3122
         RETURN
      END IF
c.....Lambda-bar
      IF (IART .EQ. -14) THEN
         INVFLV = -3122
         RETURN
      END IF 
cbz3/12/99
c.....temporary entry for Sigma's
c      IF (IART .EQ. 15) THEN
c         R = RANART(NSEED)
c         IF (R .GT. 2. / 3.) THEN
c            INVFLV = 3112
c         ELSE IF (R .GT. 1./ 3. .AND. R .LE. 2. / 3.) THEN
c            INVFLV = 3212
c         ELSE
c            INVFLV = 3222
c         END IF
c         RETURN
c      END IF
c.....Sigma-
      IF (IART .EQ. 15) THEN
         INVFLV = 3112
         RETURN
      END IF
c.....Sigma- bar
      IF (IART .EQ. -15) THEN
         INVFLV = -3112
         RETURN
      END IF 
c.....Sigma0
      IF (IART .EQ. 16) THEN
         INVFLV = 3212
         RETURN
      END IF
c.....Sigma0 -bar
      IF (IART .EQ. -16) THEN
         INVFLV = -3212
         RETURN
      END IF
c.....Sigma+
      IF (IART .EQ. 17) THEN
         INVFLV = 3222
         RETURN
      END IF
c.....Sigma+ -bar
      IF (IART .EQ. -17) THEN
         INVFLV = -3222
         RETURN
      END IF 
clin-2/23/03 K0S and K0L are generated at the last timestep:
c.....temporary entry for K- and K0bar
      IF (IART .EQ. 21) THEN
c         R = RANART(NSEED)
c         IF (R .GT. 0.5) THEN
            INVFLV = -321
c         ELSE
c            INVFLV = -311
c            R = RANART(NSEED)
c            IF (R .GT. 0.5) THEN
c               INVFLV = 310
c            ELSE
c               INVFLV = 130
c            END IF
c         END IF
         RETURN
      END IF
c.....temporary entry for K+ and K0
      IF (IART .EQ. 23) THEN
c         R = RANART(NSEED)
c         IF (R .GT. 0.5) THEN
            INVFLV = 321
c         ELSE
c            INVFLV = 311
c            R = RANART(NSEED)
c            IF (R .GT. 0.5) THEN
c               INVFLV = 310
c            ELSE
c               INVFLV = 130
c            END IF
c         END IF
         RETURN
      END IF
c.....K0Long:
      IF (IART .EQ. 22) THEN
         INVFLV = 130
         RETURN
      ENDIF
c.....K0Short:
      IF (IART .EQ. 24) THEN
         INVFLV = 310
         RETURN
      ENDIF
c.....rho-
      IF (IART .EQ. 25) THEN
         INVFLV = -213
         RETURN
      END IF
c.....rho0
      IF (IART .EQ. 26) THEN
         INVFLV = 113
         RETURN
      END IF
c.....rho+
      IF (IART .EQ. 27) THEN
         INVFLV = 213
         RETURN
      END IF
c.....omega
      IF (IART .EQ. 28) THEN
         INVFLV = 223
         RETURN
      END IF
c.....phi
      IF (IART .EQ. 29) THEN
         INVFLV = 333
         RETURN
      END IF
c.....temporary entry for K*+ and K*0
      IF (IART .EQ. 30) THEN
         INVFLV = 323
         IF (RANART(NSEED).GT.0.5) INVFLV = 313
         RETURN
      END IF
c.....temporary entry for K*- and K*0bar
      IF (IART .EQ. -30) THEN
         INVFLV = -323
         IF (RANART(NSEED).GT.0.5) INVFLV = -313
         RETURN
      END IF
c... eta-prime (bar)
      IF (IART .EQ. 31) THEN
         INVFLV = 331
         RETURN
      END IF
c... a1
      IF (IART .EQ. 32) THEN
         INVFLV = 777
         RETURN
      END IF
c... cascade-
      IF (IART .EQ. 40) THEN
         INVFLV = 3312
         RETURN
      END IF                   
c... cascade+ (bar)
      IF (IART .EQ. -40) THEN
         INVFLV = -3312
         RETURN
      END IF
c... cascade0
      IF (IART .EQ. 41) THEN
         INVFLV = 3322
         RETURN
      END IF
c... cascade0 -bar
      IF (IART .EQ. -41) THEN
         INVFLV = -3322
         RETURN
      END IF
c... Omega-
      IF (IART .EQ. 45) THEN
         INVFLV = 3334
         RETURN
      END IF
c... Omega+ (bar)
      IF (IART .EQ. -45) THEN
         INVFLV = -3334
         RETURN
      END IF
c... Di-Omega
      IF (IART .EQ. 44) THEN
         INVFLV = 6666
         RETURN
      END IF
c sp 12/19/00 end           
clin-5/2008 deuteron ID numbers in ART and ampt.dat:
      IF (IART .EQ. 42) THEN
         INVFLV = 42
         RETURN
      ELSEIF (IART .EQ. -42) THEN         
         INVFLV = -42
         RETURN
      END IF
c
c.....other
      INVFLV = IART - 10000
      RETURN
      END
c=======================================================================
