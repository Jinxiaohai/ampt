        FUNCTION RWDSAX(X)
        SAVE   
              RWDSAX=X*X*WDSAX(X)
              RETURN
              END
C
C
C
C
C The next three subroutines are for Monte Carlo generation 
C according to a given function FHB. One calls first HIFUN 
C with assigned channel number I, low and up limits. Then to 
C generate the distribution one can call HIRND(I) which gives 
C you a random number generated according to the given function.
C 
