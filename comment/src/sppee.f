        SUBROUTINE sppee(lb1,lb2,srt)
        parameter (ETAM=0.5475)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
cc      SAVE /ppb1/
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
cc      SAVE /ppmm/
      SAVE   
        ppee=0.
        if((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.3.and.lb2.le.5)) then
           if(srt.gt.(2*ETAM)) ppee=ptoe(srt)
        elseif(lb1.eq.0.and.lb2.eq.0) then
           ppee=etop(srt)
        endif
        return
        END
*****************************************
* for pi pi -> eta eta, determined from detailed balance, spin-isospin averaged
