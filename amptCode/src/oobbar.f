      real function oobbar(srt)
      parameter (pimass=0.140,arho=0.77,aomega=0.782)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
      sppb6p=xppbar(srt)*factr2(6)*ene**4/fsum
      pi2=(s-4*aomega**2)/4
      oobbar=4./9.*sppb6p/pi2*wtot
      return
      END
