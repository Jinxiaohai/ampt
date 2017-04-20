c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    更新QRAN(I)的数值。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      SUBROUTINE ARAN9(QRAN,NDIM)
      DIMENSION QRAN(10)
      COMMON/SEDVAX/NUM1
      SAVE   
      DO 1 I=1,NDIM
    1 QRAN(I)=RANART(NUM1)
      RETURN
      END
C
C
C*********GAUSSIAN ONE-DIMENSIONAL INTEGRATION PROGRAM*************
C
