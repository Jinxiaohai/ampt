        SUBROUTINE CRSJET
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        REAL HIPR1(100),HINT1(100)
        COMMON/HPARNT/HIPR1,IHPR2(50),HINT1,IHNT2(50)
cc      SAVE /HPARNT/
        COMMON/NJET/N,ipcrs
cc      SAVE /NJET/
        COMMON/BVEG1/XL(10),XU(10),ACC,NDIM,NCALL,ITMX,NPRN
cc      SAVE /BVEG1/
        COMMON/BVEG2/XI(50,10),SI,SI2,SWGT,SCHI,NDO,IT
cc      SAVE /BVEG2/
        COMMON/BVEG3/F,TI,TSI
cc      SAVE /BVEG3/
        COMMON/SEDVAX/NUM1
cc      SAVE /SEDVAX/
        EXTERNAL FJET,FJETRG
        SAVE   
C
c************************
c        NCALL give the number of inner-iteration, ITMX 
C       gives the limit of out-iteration. Nprn is an option
C       ( 1: print the integration process. 0: do not print)
C
        NDIM=3
        ipcrs=0
        CALL VEGAS(FJET,AVGI,SD,CHI2A)
        HINT1(14)=sngl(AVGI)/2.5682
        IF(IHPR2(6).EQ.1 .AND. IHNT2(1).GT.1) THEN
                ipcrs=1
                CALL VEGAS(FJET,AVGI,SD,CHI2A)
                HINT1(15)=sngl(AVGI)/2.5682
        ENDIF
        IF(IHPR2(6).EQ.1 .AND. IHNT2(3).GT.1) THEN
                ipcrs=2
                CALL VEGAS(FJET,AVGI,SD,CHI2A)
                HINT1(16)=sngl(AVGI)/2.5682
        ENDIF
        IF(IHPR2(6).EQ.1.AND.IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
                ipcrs=3
                CALL VEGAS(FJET,AVGI,SD,CHI2A)
                HINT1(17)=sngl(AVGI)/2.5682
        ENDIF
C                ********Total inclusive jet cross section(Pt>P0) 
C
        IF(IHPR2(3).NE.0) THEN
           ipcrs=0
           CALL VEGAS(FJETRG,AVGI,SD,CHI2A)
           HINT1(61)=sngl(AVGI)/2.5682
           IF(IHPR2(6).EQ.1 .AND. IHNT2(1).GT.1) THEN
              ipcrs=1
              CALL VEGAS(FJETRG,AVGI,SD,CHI2A)
              HINT1(62)=sngl(AVGI)/2.5682
           ENDIF
           IF(IHPR2(6).EQ.1 .AND. IHNT2(3).GT.1) THEN
              ipcrs=2
              CALL VEGAS(FJETRG,AVGI,SD,CHI2A)
              HINT1(63)=sngl(AVGI)/2.5682
           ENDIF
           IF(IHPR2(6).EQ.1.AND.IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
              ipcrs=3
              CALL VEGAS(FJETRG,AVGI,SD,CHI2A)
              HINT1(64)=sngl(AVGI)/2.5682
           ENDIF
        ENDIF
C                        ********cross section of trigger jet
C
        RETURN
        END
C
C
C
