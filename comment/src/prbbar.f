      real function prbbar(srt)
*****************************************
      parameter (pimass=0.140,arho=0.77,aomega=0.782)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
cc      SAVE /ppb1/
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
cc      SAVE /ppmm/
      SAVE   
      sppb3p=xppbar(srt)*factr2(3)*ene/fsum
      pi2=(s-(pimass+arho)**2)*(s-(pimass-arho)**2)/4/s
      prbbar=4./27.*sppb3p/pi2*wtot
      return
      END
*****************************************
* for rho+rho-->Bbar B                                                      *
c      real*4 function rrbbar(srt)
