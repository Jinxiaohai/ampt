        FUNCTION FTOT(X)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
        SAVE   
        OMG=OMG0(X)*(HIPR1(30)+HINT1(11))/HIPR1(31)/2.0
        FTOT=2.0*(1.0-EXP(-OMG))
        RETURN
        END
