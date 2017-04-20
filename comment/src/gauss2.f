        FUNCTION GAUSS2(F,A,B,EPS)
        EXTERNAL F
        DIMENSION W(12),X(12)
        SAVE   
        DATA CONST/1.0E-12/
        DATA W/0.1012285,.2223810,.3137067,.3623838,.0271525,
     &         .0622535,0.0951585,.1246290,.1495960,.1691565,
     &         .1826034,.1894506/
        DATA X/0.9602899,.7966665,.5255324,.1834346,.9894009,
     &         .9445750,0.8656312,.7554044,.6178762,.4580168,
     &         .2816036,.0950125/
        DELTA=CONST*ABS(A-B)
        GAUSS2=0.0
        AA=A
5        Y=B-AA
        IF(ABS(Y).LE.DELTA) RETURN
2        BB=AA+Y
        C1=0.5*(AA+BB)
        C2=C1-AA
        S8=0.0
        S16=0.0
        DO 1 I=1,4
        U=X(I)*C2
1        S8=S8+W(I)*(F(C1+U)+F(C1-U))
        DO 3 I=5,12
        U=X(I)*C2
3        S16=S16+W(I)*(F(C1+U)+F(C1-U))
        S8=S8*C2
        S16=S16*C2
        IF(ABS(S16-S8).GT.EPS*(1.+ABS(S16))) GOTO 4
        GAUSS2=GAUSS2+S16
        AA=BB
        GOTO 5
4        Y=0.5*Y
        IF(ABS(Y).GT.DELTA) GOTO 2
        WRITE(6,7)
        GAUSS2=0.0
        RETURN
7        FORMAT(1X,'GAUSS2....TOO HIGH ACURACY REQUIRED')
        END
C
C
C
C
C
