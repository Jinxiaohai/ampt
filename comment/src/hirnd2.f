        FUNCTION HIRND2(I,XMIN,XMAX)
        COMMON/HIJHB/RR(10,201),XX(10,201)
cc      SAVE /HIJHB/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
        SAVE   
        IF(XMIN.LT.XX(I,1)) XMIN=XX(I,1)
        IF(XMAX.GT.XX(I,201)) XMAX=XX(I,201)
        JMIN=1+int(200*(XMIN-XX(I,1))/(XX(I,201)-XX(I,1)))
        JMAX=1+int(200*(XMAX-XX(I,1))/(XX(I,201)-XX(I,1)))
        RX=RR(I,JMIN)+(RR(I,JMAX)-RR(I,JMIN))*RANART(NSEED)
        JL=0
        JU=202
10        IF(JU-JL.GT.1) THEN
           JM=(JU+JL)/2
           IF((RR(I,201).GT.RR(I,1)).EQV.(RX.GT.RR(I,JM))) THEN
              JL=JM
           ELSE
              JU=JM
           ENDIF
        GO TO 10
        ENDIF
        J=JL
        IF(J.LT.1) J=1
        IF(J.GE.201) J=200
        HIRND2=(XX(I,J)+XX(I,J+1))/2.0
        RETURN
        END        
C
C
C
C
