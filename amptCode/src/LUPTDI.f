      SUBROUTINE LUPTDI(KFL,PX,PY)  
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      KFLA=IABS(KFL)    
      PT=PARJ(21)*SQRT(-LOG(MAX(1E-10,RLU(0)))) 
      IF(MSTJ(91).EQ.1) PT=PARJ(22)*PT  
      IF(KFLA.EQ.0.AND.MSTJ(13).LE.0) PT=0. 
      PHI=PARU(2)*RLU(0)    
      PX=PT*COS(PHI)    
      PY=PT*SIN(PHI)    
      RETURN    
      END   
