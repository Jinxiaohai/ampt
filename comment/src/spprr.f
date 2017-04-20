        SUBROUTINE spprr(lb1,lb2,srt)
        parameter (arho=0.77)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
cc      SAVE /ppb1/
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
cc      SAVE /ppmm/
      SAVE   
        pprr=0.
        if((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.3.and.lb2.le.5)) then
c     for now, rho mass taken to be the central value in these two processes
           if(srt.gt.(2*arho)) pprr=ptor(srt)
        elseif((lb1.ge.25.and.lb1.le.27).and.(lb2.ge.25.and.lb2.le.27)) 
     1          then
           pprr=rtop(srt)
        endif
c
        return
        END
*****************************************
* for pi pi -> rho rho, determined from detailed balance
