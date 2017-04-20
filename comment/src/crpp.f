      SUBROUTINE CRPP(PX,PY,PZ,SRT,I1,I2,IBLOCK,
     &ppel,ppin,spprho,ipp)
*     PURPOSE:                                                         *
*             DEALING WITH PION-PION COLLISIONS                        *
*     NOTE   :                                                         *
*           VALID ONLY FOR PION-PION-DISTANCES LESS THAN 2.5 FM        *
*     QUANTITIES:                                                 *
*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
*           SRT      - SQRT OF S                                       *
*           IBLOCK   - THE INFORMATION BACK                            *
*                     6-> Meson+Meson elastic
*                     66-> Meson+meson-->K+K-
**********************************
      PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1     AMP=0.93828,AP1=0.13496,
     2 AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
      PARAMETER      (AKA=0.498,aks=0.895)
      parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
      COMMON /AA/ R(3,MAXSTR)
cc      SAVE /AA/
      COMMON /BB/ P(3,MAXSTR)
cc      SAVE /BB/
      COMMON /CC/ E(MAXSTR)
cc      SAVE /CC/
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
cc      SAVE /ppb1/
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
cc      SAVE /ppmm/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
      lb1i=lb(i1)
      lb2i=lb(i2)
       PX0=PX
       PY0=PY
       PZ0=PZ
        iblock=1
*-----------------------------------------------------------------------
* check Meson+Meson inelastic collisions
clin-9/28/00
c        if((srt.gt.1.).and.(ppin/(ppin+ppel).gt.RANART(NSEED)))then
c        iblock=66
c        e(i1)=0.498
c        e(i2)=0.498
c        lb(i1)=21
c        lb(i2)=23
c        go to 10
clin-11/07/00
c        if(srt.gt.1.and.(ppin/(ppin+ppel)).gt.RANART(NSEED)) then
clin-4/03/02
        if(srt.gt.(2*aka).and.(ppin/(ppin+ppel)).gt.RANART(NSEED)) then
c        if(ppin/(ppin+ppel).gt.RANART(NSEED)) then
clin-10/08/00
           ranpi=RANART(NSEED)
           if((pprr/ppin).ge.ranpi) then
c     1) pi pi <-> rho rho:
              call pi2ro2(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
clin-4/03/02 eta equilibration:
           elseif((pprr+ppee)/ppin.ge.ranpi) then
c     4) pi pi <-> eta eta:
              call pi2et2(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
           elseif(((pprr+ppee+pppe)/ppin).ge.ranpi) then
c     5) pi pi <-> pi eta:
              call pi3eta(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
           elseif(((pprr+ppee+pppe+rpre)/ppin).ge.ranpi) then
c     6) rho pi <-> pi eta:
              call rpiret(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
           elseif(((pprr+ppee+pppe+rpre+xopoe)/ppin).ge.ranpi) then
c     7) omega pi <-> omega eta:
              call opioet(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
           elseif(((pprr+ppee+pppe+rpre+xopoe+rree)
     1             /ppin).ge.ranpi) then
c     8) rho rho <-> eta eta:
              call ro2et2(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
clin-4/03/02-end
c     2) BBbar production:
           elseif(((pprr+ppee+pppe+rpre+xopoe+rree+ppinnb)/ppin)
     1             .ge.ranpi) then
              call bbarfs(lbb1,lbb2,ei1,ei2,iblock,iseed)
c     3) KKbar production:
           else
              iblock=66
              ei1=aka
              ei2=aka
              lbb1=21
              lbb2=23
clin-11/07/00 pi rho -> K* Kbar and K*bar K productions:
              lb1=lb(i1)
              lb2=lb(i2)
clin-2/13/03 include omega the same as rho, eta the same as pi:
c        if(((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.25.and.lb2.le.27))
c     1  .or.((lb2.ge.3.and.lb2.le.5).and.(lb1.ge.25.and.lb1.le.27)))
        if( ( (lb1.eq.0.or.(lb1.ge.3.and.lb1.le.5))
     1       .and.(lb2.ge.25.and.lb2.le.28))
     2       .or. ( (lb2.eq.0.or.(lb2.ge.3.and.lb2.le.5))
     3       .and.(lb1.ge.25.and.lb1.le.28))) then
           ei1=aks
           ei2=aka
           if(RANART(NSEED).ge.0.5) then
              iblock=366
              lbb1=30
              lbb2=21
           else
              iblock=367
              lbb1=-30
              lbb2=23
           endif
        endif
clin-11/07/00-end
           endif
clin-ppbar-8/25/00
           e(i1)=ei1
           e(i2)=ei2
           lb(i1)=lbb1
           lb(i2)=lbb2
clin-10/08/00-end
       else
cbzdbg10/15/99
c.....for meson+meson elastic srt.le.2Mk, if not pi+pi collision return
         if ((lb(i1).lt.3.or.lb(i1).gt.5).and.
     &        (lb(i2).lt.3.or.lb(i2).gt.5)) return
cbzdbg10/15/99 end
* check Meson+Meson elastic collisions
        IBLOCK=6
* direct process
       if(ipp.eq.1.or.ipp.eq.4.or.ipp.eq.6)go to 10
       if(spprho/ppel.gt.RANART(NSEED))go to 20
       endif
10      NTAG=0
        EM1=E(I1)
        EM2=E(I2)
*-----------------------------------------------------------------------
* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
* ENERGY CONSERVATION
          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.e-09
          PR=SQRT(PR2)/(2.*SRT)
          C1   = 1.0 - 2.0 * RANART(NSEED)
          T1   = 2.0 * PI * RANART(NSEED)
          S1   = SQRT( 1.0 - C1**2 )
      CT1  = COS(T1)
      ST1  = SIN(T1)
      PZ   = PR * C1
      PX   = PR * S1*CT1 
      PY   = PR * S1*ST1
* for isotropic distribution no need to ROTATE THE MOMENTUM
* ROTATE IT 
      CALL ROTATE(PX0,PY0,PZ0,PX,PY,PZ) 
      RETURN
20       continue
       iblock=666
* treat rho formation in pion+pion collisions
* calculate the mass and momentum of rho in the nucleus-nucleus frame
       call rhores(i1,i2)
       if(ipp.eq.2)lb(i1)=27
       if(ipp.eq.3)lb(i1)=26
       if(ipp.eq.5)lb(i1)=25
       return       
      END
**********************************
**********************************
*                                                                      *
*                                                                      *
