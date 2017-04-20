      subroutine dbelangle(pxn,pyn,pzn,pfinal)
      PARAMETER (PI=3.1415926)
      COMMON/RNDF77/NSEED
      SAVE   
c     take isotropic distribution for now:
      C1=1.0-2.0*RANART(NSEED)
      T1=2.0*PI*RANART(NSEED)
      S1=SQRT(1.0-C1**2)
      CT1=COS(T1)
      ST1=SIN(T1)
* THE MOMENTUM IN THE CMS IN THE FINAL STATE
      Pzn=pfinal*C1
      Pxn=pfinal*S1*CT1 
      Pyn=pfinal*S1*ST1
      return
      end
c
c     Cross section of Deuteron+Pi elastic (in mb):
