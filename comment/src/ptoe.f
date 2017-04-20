      real function ptoe(srt)
*****************************************
      parameter (pimass=0.140,ETAM=0.5475)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
cc      SAVE /ppb1/
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
cc      SAVE /ppmm/
      SAVE   
      s2=srt**2
      ptoe=1./9.*(s2-4*etam**2)/(s2-4*pimass**2)*etop(srt)
      return
      END
*****************************************
* for eta eta -> pi pi, assumed a constant cross section (in mb)
