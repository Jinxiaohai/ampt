c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    该函数得到FTOTJT
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        FUNCTION FTOTJT(X)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
        SAVE   
        OMG=OMG0(X)*HINT1(11)/HIPR1(31)/2.0
        FTOTJT=1.0-EXP(-2.0*OMG)
        RETURN
        END
C
C
C
