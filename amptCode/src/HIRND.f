        FUNCTION HIRND(I)
        COMMON/HIJHB/RR(10,201),XX(10,201)
      COMMON/RNDF77/NSEED
        SAVE   
        RX=RANART(NSEED)
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
        HIRND=(XX(I,J)+XX(I,J+1))/2.0
        RETURN
        END        
