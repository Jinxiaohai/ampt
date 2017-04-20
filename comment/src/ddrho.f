        SUBROUTINE ddrho(SRT,ISEED,PX,PY,PZ,DM1,PNX,
     &  PNY,PNZ,DM2,PPX,PPY,PPZ,amp,icou1)
* PURPOSE : CALCULATE MOMENTUM OF PARTICLES IN THE FINAL SATAT FROM
* THE PROCESS N+N--->D1+D2+rho
*       DATE : Nov.5, 1994
* Generate the masses and momentum for particles in the NN-->DDrho process
* for a given center of mass energy srt, the momenta are given in the center
* of mass of the NN
*****************************************
        COMMON/TABLE/ xarray(0:1000),earray(0:1000)
cc      SAVE /TABLE/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
       icou1=0
       pi=3.1415926
        AMN=938.925/1000.
        AMP=770./1000.
* (1) GENGRATE THE MASS OF DELTA1 AND DELTA2 USING
       srt1=srt-amp-0.02
       ntrym=0
8       call Rmasdd(srt1,1.232,1.232,1.08,
     &  1.08,ISEED,1,dm1,dm2)
       ntrym=ntrym+1
* GENERATE THE MASS FOR THE RHO
       RHOMAX = SRT-DM1-DM2-0.02
       IF(RHOMAX.LE.0.and.ntrym.le.20)go to 8
       AMP=RHOMAS(RHOMAX,ISEED)
* CONSTANTS FOR GENERATING THE LONGITUDINAL MOMENTUM 
* FOR ONE OF THE RESONANCES
       V=0.43
       W=-0.84
* (2) Generate the transverse momentum
*     OF DELTA1
* (2.1) estimate the maximum transverse momentum
       PTMAX2=(SRT**2-(DM1+DM2+AMP)**2)*
     1  (SRT**2-(DM1-AMP-DM2)**2)/4./SRT**2
clin-9/2012: check argument in sqrt():
       scheck=PTMAX2
       if(scheck.lt.0) then
          write(99,*) 'scheck24: ', scheck
          scheck=0.
       endif
       PTMAX=SQRT(scheck)*1./3.
c       PTMAX=SQRT(PTMAX2)*1./3.
7       PT=PTR(PTMAX,ISEED)
* (3) GENGRATE THE LONGITUDINAL MOMENTUM FOR DM1
*     USING THE GIVEN DISTRIBUTION
* (3.1) THE MAXIMUM LONGITUDINAL MOMENTUM IS
       PZMAX2=(SRT**2-(DM1+DM2+AMP)**2)*
     1  (SRT**2-(DM1-AMP-DM2)**2)/4./SRT**2-PT**2
       IF((PZMAX2.LT.0.).and.ntrym.le.100)then 
       go to 7
       else
       pzmax2=1.E-06
       endif
       PZMAX=SQRT(PZMAX2)
       XMAX=2.*PZMAX/SRT
* (3.2) THE GENERATED X IS
* THE DSTRIBUTION HAS A MAXIMUM AT X0=-V/(2*w), f(X0)=1.056
       ntryx=0
       fmax00=1.056
       x00=0.26
       if(abs(xmax).gt.0.26)then
       f00=fmax00
       else
       f00=1.+v*abs(xmax)+w*xmax**2
       endif
9       X=XMAX*(1.-2.*RANART(NSEED))
       ntryx=ntryx+1
       xratio=(1.+V*ABS(X)+W*X**2)/f00       
clin-8/17/00       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9       
       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9       
* (3.5) THE PZ IS
       PZ=0.5*SRT*X
* The x and y components of the delta1
       fai=2.*pi*RANART(NSEED)
       Px=pt*cos(fai)
       Py=pt*sin(fai)
* find the momentum of delta2 and rho
* the energy of the delta1
       ek=sqrt(dm1**2+PT**2+Pz**2)
* (1) Generate the momentum of the delta2 in the cms of delta2 and rho
*     the energy of the cms of Drho
        eln=srt-ek
       IF(ELN.lE.0)then
       icou1=-1
       return
       endif
* beta and gamma of the cms of delta2 and rho
       bx=-Px/eln
       by=-Py/eln
       bz=-Pz/eln
clin-9/2012: check argument in sqrt():
       scheck=1.-bx**2-by**2-bz**2
       if(scheck.le.0) then
          write(99,*) 'scheck25: ', scheck
          stop
       endif
       ga=1./sqrt(scheck)
c       ga=1./sqrt(1.-bx**2-by**2-bz**2)
       elnc=eln/ga
       pn2=((elnc**2+dm2**2-amp**2)/(2.*elnc))**2-dm2**2
       if(pn2.le.0)then
       icou1=-1
       return
       endif
       pn=sqrt(pn2)
clin-10/25/02 get rid of argument usage mismatch in PTR():
        xptr=0.33*PN
c       PNT=PTR(0.33*PN,ISEED)
       PNT=PTR(xptr,ISEED)
clin-10/25/02-end
       fain=2.*pi*RANART(NSEED)
       pnx=pnT*cos(fain)
       pny=pnT*sin(fain)
       SIG=1
       IF(X.GT.0)SIG=-1
clin-9/2012: check argument in sqrt():
       scheck=pn**2-PNT**2
       if(scheck.lt.0) then
          write(99,*) 'scheck26: ', scheck
          scheck=0.
       endif
       pnz=SIG*SQRT(scheck)
c       pnz=SIG*SQRT(pn**2-PNT**2)
       en=sqrt(dm2**2+pnx**2+pny**2+pnz**2)
* (2) the momentum for the rho
       ppx=-pnx
       ppy=-pny
       ppz=-pnz
       ep=sqrt(amp**2+ppx**2+ppy**2+ppz**2)
* (3) for the delta2, LORENTZ-TRANSFORMATION INTO nn cms FRAME
        PBETA  = PnX*BX + PnY*By+ PnZ*Bz
              TRANS0  = GA * ( GA * PBETA / (GA + 1.) + En )
              Pnx = BX * TRANS0 + PnX
              Pny = BY * TRANS0 + PnY
              Pnz = BZ * TRANS0 + PnZ
* (4) for the rho, LORENTZ-TRANSFORMATION INTO nn cms FRAME
             if(ep.eq.0.)ep=1.e-09
              PBETA  = PPX*BX + PPY*By+ PPZ*Bz
              TRANS0  = GA * ( GA * PBETA / (GA + 1.) + EP )
              PPx = BX * TRANS0 + PPX
              PPy = BY * TRANS0 + PPY
              PPz = BZ * TRANS0 + PPZ
       return
       end
****************************************
