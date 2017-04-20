      BLOCK DATA PPBDAT 
      parameter (AMP=0.93828,AMN=0.939457,
     1     AM0=1.232,AM1440 = 1.44, AM1535 = 1.535)
c     to give default values to parameters for BbarB production from mesons
      COMMON/ppbmas/niso(15),nstate,ppbm(15,2),thresh(15),weight(15)
cc      SAVE /ppbmas/
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
cc      SAVE /ppb1/
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
cc      SAVE /ppmm/
      SAVE   
c     thresh(i) gives the mass thresh for final channel i:
      DATA thresh/1.87656,1.877737,1.878914,2.17028,
     1     2.171457,2.37828,2.379457,2.464,2.47328,2.474457,
     2     2.672,2.767,2.88,2.975,3.07/
c     ppbm(i,j=1,2) gives masses for the two final baryons of channel i,
c     with j=1 for the lighter baryon:
      DATA (ppbm(i,1),i=1,15)/amp,amp,amn,amp,amn,amp,amn,
     1     am0,amp,amn,am0,am0,am1440,am1440,am1535/
      DATA (ppbm(i,2),i=1,15)/amp,amn,amn,am0,am0,am1440,am1440,
     1     am0,am1535,am1535,am1440,am1535,am1440,am1535,am1535/
c     factr2(i) gives weights for producing i pions from ppbar annihilation:
      DATA factr2/0,1,1.17e-01,3.27e-03,3.58e-05,1.93e-07/
c     niso(i) gives the degeneracy factor for final channel i:
      DATA niso/1,2,1,16,16,4,4,64,4,4,32,32,4,8,4/
      END   
*****************************************
* get the number of BbarB states available for mm collisions of energy srt 
