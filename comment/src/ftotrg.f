c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     该函数得到FTOTRG的值。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        FUNCTION FTOTRG(X)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
        SAVE   
        OMG=OMG0(X)*HINT1(60)/HIPR1(31)/2.0
        FTOTRG=1.0-EXP(-2.0*OMG)
        RETURN
        END
C
C
C
C
