        FUNCTION FNKICK(X)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
        SAVE   
        FNKICK=1.0/(X+HIPR1(19)**2)/(X+HIPR1(20)**2)
     &                /(1+EXP((SQRT(X)-HIPR1(20))/0.4))
        RETURN
        END
C
C
