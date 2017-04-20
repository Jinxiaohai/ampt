      real function PNLKA(srt)
      SAVE   
* units: fm**2
***********************************C
      ala=1.116
      aka=0.498
      ana=0.939
      t1=ala+aka      
      if(srt.le.t1) THEN
      Pnlka=0
      Else
      IF(SRT.LT.1.7)sbbk=(0.9/0.091)*(SRT-T1)
      IF(SRT.GE.1.7)sbbk=0.09/(SRT-1.6)
      Pnlka=0.25*sbbk
* give the cross section in units of fm**2
       pnlka=pnlka/10.
      endif     
      return
      end
*-------------------------------------------------------------------------
*****subprogram * kaon production from pi+B collisions *******************
