       SUBROUTINE XphiB(LB1, LB2, EM1, EM2, SRT,
     &                  XSK1, XSK2, XSK3, XSK4, XSK5, SIGP)
c
* ***************************
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AP2=0.13957,AM0=1.232,PI=3.1415926)
          PARAMETER (AKA=0.498, ALA = 1.1157, PIMASS=0.140, APHI=1.02)
        parameter (arho=0.77)
      SAVE   
       SIGP = 1.E-08
        XSK1 = 0.0
        XSK2 = 0.0
        XSK3 = 0.0
        XSK4 = 0.0
        XSK5 = 0.0
        XSK6 = 0.0
          srrt = srt - (em1+em2)
c* phi + N(D) -> elastic scattering
c            XSK1 = 0.56  !! mb
c  !! mb  (photo-production xsecn used)
            XSK1 = 8.00
c
c* phi + N(D) -> pi + N
        IF (srt  .GT. (ap1+amn)) THEN
             XSK2 = 0.0235*srrt**(-0.519) 
        END IF
c
c* phi + N(D) -> pi + D
        IF (srt  .GT. (ap1+am0)) THEN
            if(srrt .lt. 0.7)then
             XSK3 = 0.0119*srrt**(-0.534)
            else
             XSK3 = 0.0130*srrt**(-0.304)
            endif      
        END IF
c
c* phi + N(D) -> rho + N
        IF (srt  .GT. (arho+amn)) THEN
           if(srrt .lt. 0.7)then
             XSK4 = 0.0166*srrt**(-0.786)
            else
             XSK4 = 0.0189*srrt**(-0.277)
            endif
        END IF
c
c* phi + N(D) -> rho + D   (same as pi + D)
        IF (srt  .GT. (arho+am0)) THEN
            if(srrt .lt. 0.7)then
             XSK5 = 0.0119*srrt**(-0.534)
            else
             XSK5 = 0.0130*srrt**(-0.304)
            endif      
        END IF
c
c* phi + N -> K+ + La
       IF( (lb1.ge.1.and.lb1.le.2) .or. (lb2.ge.1.and.lb2.le.2) )THEN
        IF (srt  .GT. (aka+ala)) THEN
           XSK6 = 1.715/((srrt+3.508)**2-12.138)  
        END IF
       END IF
        SIGP = XSK1 + XSK2 + XSK3 + XSK4 + XSK5 + XSK6
       RETURN
        END
c
**********************************
*
