        subroutine resdec(i1,nt,nnn,wid,idecay,ipion)
        PARAMETER (hbarc=0.19733)
        PARAMETER (AK0=0.498,APICH=0.140,API0=0.135,AN=0.940,ADDM=0.02)
        PARAMETER (MAXSTR=150001, MAXR=1)
        COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &       IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
cc      SAVE /INPUT2/
        COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
cc      SAVE /LUJETS/
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
cc      SAVE /LUDAT1/
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
cc      SAVE /LUDAT2/
        COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
cc      SAVE /LUDAT3/
        COMMON /CC/ E(MAXSTR)
cc      SAVE /CC/
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
        COMMON   /PA/RPION(3,MAXSTR,MAXR)
cc      SAVE /PA/
        COMMON   /PB/PPION(3,MAXSTR,MAXR)
cc      SAVE /PB/
        COMMON   /PC/EPION(MAXSTR,MAXR)
cc      SAVE /PC/
        COMMON   /PD/LPION(MAXSTR,MAXR)
cc      SAVE /PD/
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
        common/resdcy/NSAV,iksdcy
cc      SAVE /resdcy/
        common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1       px1n,py1n,pz1n,dp1n
cc      SAVE /leadng/
        EXTERNAL IARFLV, INVFLV
        COMMON/tdecay/tfdcy(MAXSTR),tfdpi(MAXSTR,MAXR),tft(MAXSTR)
cc      SAVE /tdecay/
        COMMON/RNDF77/NSEED
        COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1       dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2       dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
cc      SAVE /RNDF77/
        common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
        SAVE   
        irun=idecay
clin-4/2012 for option of pi0 decay:
        if(nt.eq.ntmax.and.ipi0dcy.eq.1
     &       .and.((lb1.eq.4.and.ipion.eq.0).or.ipion.ge.1)) then
           kf=111
c        if(lb1.eq.0.or.lb1.eq.25.or.lb1.eq.26.or.lb1.eq.27
        elseif(lb1.eq.0.or.lb1.eq.25.or.lb1.eq.26.or.lb1.eq.27
     &       .or.lb1.eq.28.or.lb1.eq.29.or.iabs(lb1).eq.30
     &       .or.lb1.eq.24.or.(iabs(lb1).ge.6.and.iabs(lb1).le.9) 
     &       .or.iabs(lb1).eq.16) then
           kf=INVFLV(lb1)
        else
           return
        endif
c
        IP=1
c     label as undecayed and the only particle in the record:
        N=1
        K(IP,1)=1
        K(IP,3)=0
        K(IP,4)=0
        K(IP,5)=0
c
        K(IP,2)=kf
clin-4/2012 for option of pi0 decay:
        if(ipion.eq.0) then
c
        P(IP,1)=px1
        P(IP,2)=py1
        P(IP,3)=pz1
c        em1a=em1
c     eta or omega in ART may be below or too close to (pi+pi-pi0) mass, 
c     causing LUDECY error,thus increase their mass ADDM above this thresh,
c     noting that rho (m=0.281) too close to 2pi thrshold fails to decay:
        if((lb1.eq.0.or.lb1.eq.28).and.em1.lt.(2*APICH+API0+ADDM)) then
           em1=2*APICH+API0+ADDM
c     rho
        elseif(lb1.ge.25.and.lb1.le.27.and.em1.lt.(2*APICH+ADDM)) then
           em1=2*APICH+ADDM
c     K*
        elseif(iabs(lb1).eq.30.and.em1.lt.(APICH+AK0+ADDM)) then
           em1=APICH+AK0+ADDM
c     Delta created in ART may be below (n+pich) mass, causing LUDECY error:
        elseif(iabs(lb1).ge.6.and.iabs(lb1).le.9
     1          .and.em1.lt.(APICH+AN+ADDM)) then
           em1=APICH+AN+ADDM
        endif
c        if(em1.ge.(em1a+0.01)) write (6,*) 
c     1       'Mass increase in resdec():',nt,em1-em1a,lb1
        e1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
        P(IP,4)=e1
        P(IP,5)=em1
clin-5/2008:
        dpdecp=dpertp(i1)
clin-4/2012 for option of pi0 decay:
        elseif(nt.eq.ntmax.and.ipi0dcy.eq.1.and.ipion.ge.1) then        
           P(IP,1)=PPION(1,ipion,IRUN)
           P(IP,2)=PPION(2,ipion,IRUN)
           P(IP,3)=PPION(3,ipion,IRUN)
           P(IP,5)=EPION(ipion,IRUN)
           P(IP,4)=SQRT(P(IP,5)**2+P(IP,1)**2+P(IP,2)**2+P(IP,3)**2)
           dpdecp=dppion(ipion,IRUN)
ctest off
c           write(99,*) P(IP,4), P(IP,5), dpdecp, ipion, wid
        else
           print *, 'stopped in resdec() a'
           stop
        endif
c
        call ludecy(IP)
c     add decay time to daughter's formation time at the last timestep:
        if(nt.eq.ntmax) then
           tau0=hbarc/wid
           taudcy=tau0*(-1.)*alog(1.-RANART(NSEED))
           ndaut=n-nsav
           if(ndaut.le.1) then
              write(10,*) 'note: ndaut(<1)=',ndaut
              call lulist(2)
              stop
            endif
c     lorentz boost:
clin-4/2012 for option of pi0 decay:
            if(ipion.eq.0) then
               taudcy=taudcy*e1/em1
               tfnl=tfnl+taudcy
               xfnl=xfnl+px1/e1*taudcy
               yfnl=yfnl+py1/e1*taudcy
               zfnl=zfnl+pz1/e1*taudcy
            elseif(ipion.ge.1) then
               taudcy=taudcy*P(IP,4)/P(IP,5)
               tfnl=tfdpi(ipion,IRUN)+taudcy
               xfnl=RPION(1,ipion,IRUN)+P(IP,1)/P(IP,4)*taudcy
               yfnl=RPION(2,ipion,IRUN)+P(IP,2)/P(IP,4)*taudcy
               zfnl=RPION(3,ipion,IRUN)+P(IP,3)/P(IP,4)*taudcy
            else
               print *, 'stopped in resdec() b',ipion,wid,P(ip,4)
               stop
            endif
c     at the last timestep, assign rho, K0S or eta (decay daughter)
c     to lb(i1) only (not to lpion) in order to decay them again:
clin-4/2012 for option of pi0 decay:
c           if(n.ge.(nsav+2)) then
           if(n.ge.(nsav+2).and.ipion.eq.0) then
              do 1001 idau=nsav+2,n
                 kdaut=K(idau,2)
                 if(kdaut.eq.221.or.kdaut.eq.113
     1                .or.kdaut.eq.213.or.kdaut.eq.-213
     2                .or.kdaut.eq.310) then
c     switch idau and i1(nsav+1):
                    ksave=kdaut
                    pxsave=p(idau,1)
                    pysave=p(idau,2)
                    pzsave=p(idau,3)
                    esave=p(idau,4)
                    xmsave=p(idau,5)
                    K(idau,2)=K(nsav+1,2)
                    p(idau,1)=p(nsav+1,1)
                    p(idau,2)=p(nsav+1,2)
                    p(idau,3)=p(nsav+1,3)
                    p(idau,4)=p(nsav+1,4)
                    p(idau,5)=p(nsav+1,5)
                    K(nsav+1,2)=ksave
                    p(nsav+1,1)=pxsave
                    p(nsav+1,2)=pysave
                    p(nsav+1,3)=pzsave
                    p(nsav+1,4)=esave
                    p(nsav+1,5)=xmsave
c     note: phi decay may produce rho, K0s or eta, N*(1535) decay may produce 
c     eta, but only one daughter may be rho, K0s or eta:
                    goto 111
                 endif
 1001         continue
           endif
 111       continue
c     
           enet=0.
           do 1002 idau=nsav+1,n
              enet=enet+p(idau,4)
 1002      continue
c           if(abs(enet-e1).gt.0.02) 
c     1          write(93,*) 'resdec(): nt=',nt,enet-e1,lb1
        endif
        do 1003 idau=nsav+1,n
           kdaut=K(idau,2)
           lbdaut=IARFLV(kdaut)
c     K0S and K0L are named K+/K- during hadron cascade, and only 
c     at the last timestep they keep their real LB # before output;
c     K0/K0bar (from K* decay) converted to K0S and K0L at the last timestep:
           if(nt.eq.ntmax.and.(kdaut.eq.130.or.kdaut.eq.310
     1          .or.iabs(kdaut).eq.311)) then
              if(kdaut.eq.130) then
                 lbdaut=22
              elseif(kdaut.eq.310) then
                 lbdaut=24
              elseif(iabs(kdaut).eq.311) then
                 if(RANART(NSEED).lt.0.5) then
                    lbdaut=22
                 else
                    lbdaut=24
                 endif
              endif
           endif
c
           if(idau.eq.(nsav+1)) then
clin-4/2012 for option of pi0 decay:
              if(ipion.eq.0) then
                 LB(i1)=lbdaut
                 E(i1)=p(idau,5)
                 px1n=p(idau,1)
                 py1n=p(idau,2)
                 pz1n=p(idau,3)
clin-5/2008:
                 dp1n=dpdecp
              elseif(ipion.ge.1) then
                 LPION(ipion,IRUN)=lbdaut
                 EPION(ipion,IRUN)=p(idau,5)
                 PPION(1,ipion,IRUN)=p(idau,1)
                 PPION(2,ipion,IRUN)=p(idau,2)
                 PPION(3,ipion,IRUN)=p(idau,3)
                 RPION(1,ipion,IRUN)=xfnl
                 RPION(2,ipion,IRUN)=yfnl
                 RPION(3,ipion,IRUN)=zfnl
                 tfdpi(ipion,IRUN)=tfnl
                 dppion(ipion,IRUN)=dpdecp
              endif
c
           else
              nnn=nnn+1
              LPION(NNN,IRUN)=lbdaut
              EPION(NNN,IRUN)=p(idau,5)
              PPION(1,NNN,IRUN)=p(idau,1)
              PPION(2,NNN,IRUN)=p(idau,2)
              PPION(3,NNN,IRUN)=p(idau,3)
              RPION(1,NNN,IRUN)=xfnl
              RPION(2,NNN,IRUN)=yfnl
              RPION(3,NNN,IRUN)=zfnl
              tfdpi(NNN,IRUN)=tfnl
clin-5/2008:
              dppion(NNN,IRUN)=dpdecp
           endif
 1003   continue
        return
        end
c=======================================================================
