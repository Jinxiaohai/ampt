      REAL FUNCTION FDR(DMASS,aj,al,width,widb0,srt,em1,em2)
      SAVE   
        AMd=em1
        AmP=em2
           Ak02= 0.25*(DMASS**2-amd**2-amp**2)**2
     &           -(Amp*amd)**2
            IF (ak02 .GT. 0.) THEN
              Q0 = SQRT(ak02/DMASS)
            ELSE
              Q0= 0.0
             fdR=0
           return
            END IF
           Ak2= 0.25*(srt**2-amd**2-amp**2)**2
     &           -(Amp*amd)**2
            IF (ak2 .GT. 0.) THEN
              Q = SQRT(ak2/DMASS)
            ELSE
              Q= 0.00
             fdR=0
             return
            END IF
       b=widb0*1.2*dmass/srt*(q/q0)**(2.*al+1)
     &  /(1.+0.2*(q/q0)**(2*al))
        FDR=(2.*aj+1)*WIDTH**2*b/((srt-dmass)**2
     1  +0.25*WIDTH**2)/(6.*q**2)
        RETURN
        END
******************************
* this program calculates the elastic cross section for pion+delta
* through higher resonances
c       REAL*4 FUNCTION DIRCT3(SRT)
