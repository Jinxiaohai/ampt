       REAL FUNCTION RMASS(DMAX,ISEED)
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
* THE MINIMUM MASS FOR DELTA
          DMIN = 1.078
* Delta(1232) production
          IF(DMAX.LT.1.232) THEN
          FM=FDELTA(DMAX)
          ELSE
          FM=1.
          ENDIF
          IF(FM.EQ.0.)FM=1.E-06
          NTRY1=0
10        DM = RANART(NSEED) * (DMAX-DMIN) + DMIN
          NTRY1=NTRY1+1
          IF((RANART(NSEED) .GT. FDELTA(DM)/FM).AND.
     1    (NTRY1.LE.10)) GOTO 10
clin-2/26/03 sometimes Delta mass can reach very high values (e.g. 15.GeV),
c     thus violating the thresh of the collision which produces it 
c     and leads to large violation of energy conservation. 
c     To limit the above, limit the Delta mass below a certain value 
c     (here taken as its central value + 2* B-W fullwidth):
          if(dm.gt.1.47) goto 10
       RMASS=DM
       RETURN
       END
*------------------------------------------------------------------
* THE Breit Wigner FORMULA
