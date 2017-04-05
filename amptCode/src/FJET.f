        FUNCTION FJET(X,WGT)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        REAL HIPR1(100),HINT1(100)
        COMMON/HPARNT/HIPR1,IHPR2(50),HINT1,IHNT2(50)
        DIMENSION X(10)
        SAVE   
        PT2=dble(HINT1(1)**2/4.0-HIPR1(8)**2)*X(1)+dble(HIPR1(8))**2
        XT=2.0d0*DSQRT(PT2)/dble(HINT1(1))
        YMX1=DLOG(1.0d0/XT+DSQRT(1.0d0/XT**2-1.0d0))
        Y1=2.0d0*YMX1*X(2)-YMX1
        YMX2=DLOG(2.0d0/XT-DEXP(Y1))
        YMN2=DLOG(2.0d0/XT-DEXP(-Y1))
        Y2=(YMX2+YMN2)*X(3)-YMN2
        FJET=2.0d0*YMX1*(YMX2+YMN2)*dble(HINT1(1)**2/4.0-HIPR1(8)**2)
     &                *G(Y1,Y2,PT2)/2.0d0
        RETURN
        END
