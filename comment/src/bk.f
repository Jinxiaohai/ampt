        FUNCTION BK(X)
        COMMON /BESEL/X4
cc      SAVE /BESEL/
        SAVE   
        BK=EXP(-X)*(X**2-X4**2)**2.50/15.0
        RETURN
        END
C
C
C        THIS PROGRAM IS TO CALCULATE THE JET CROSS SECTION
C        THE INTEGRATION IS DONE BY USING VEGAS
C
