       REAL FUNCTION RHOMAS(DMAX,ISEED)
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
* THE MINIMUM MASS FOR DELTA
          DMIN = 0.28
* RHO(770) production
          IF(DMAX.LT.0.77) THEN
          FM=FRHO(DMAX)
          ELSE
          FM=1.
          ENDIF
          IF(FM.EQ.0.)FM=1.E-06
          NTRY1=0
10        DM = RANART(NSEED) * (DMAX-DMIN) + DMIN
          NTRY1=NTRY1+1
          IF((RANART(NSEED) .GT. FRHO(DM)/FM).AND.
     1    (NTRY1.LE.10)) GOTO 10
clin-2/26/03 limit the rho mass below a certain value
c     (here taken as its central value + 2* B-W fullwidth):
          if(dm.gt.1.07) goto 10
       RHOMAS=DM
       RETURN
       END
******************************************
* for pp-->pp+2pi
c      real*4 function X2pi(srt)
