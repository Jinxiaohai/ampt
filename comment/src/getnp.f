        subroutine getnp
        PARAMETER (MAXSTR=150001)
        COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
cc      SAVE /HMAIN1/
        COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
cc      SAVE /HMAIN2/
        COMMON /HPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)
cc      SAVE /HPARNT/
        common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
cc      SAVE /snn/
        SAVE   
        if(NATT.eq.0) then
           npart1=0
           npart2=0
           return
        endif
c
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9971,*)"beforePZPROJ = ",PZPROJ,"beforePZTARG = ",PZTARG
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        PZPROJ=SQRT(HINT1(6)**2-HINT1(8)**2)
        PZTARG=SQRT(HINT1(7)**2-HINT1(9)**2)
        epsiPz=0.01
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9971,*)"afterPZPROJ = ",PZPROJ,"afterPZTARG = ",PZTARG
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
clin-9/2011-add Pt tolerance in determining spectator nucleons
c     (affect string melting runs when LAB frame is used):
        epsiPt=1e-6
c
        nspec1=0
        nspec2=0
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9971,*)"The value of NATT is =====> ", NATT
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 1000 I = 1, NATT
clin-9/2011 determine spectator nucleons consistently
c           if((KATT(I,1).eq.2112.or.KATT(I,1).eq.2212)
c     1          .and.PATT(I, 1).eq.0.and.PATT(I, 2).eq.0) then
           if((KATT(I,1).eq.2112.or.KATT(I,1).eq.2212)
     1          .and.abs(PATT(I, 1)).le.epsiPt
     2          .and.abs(PATT(I, 2)).le.epsiPt) then
              if(PATT(I, 3).gt.amax1(0.,PZPROJ-epsiPz)) then
                 nspec1=nspec1+1
              elseif(PATT(I, 3).lt.(-PZTARG+epsiPz)) then
                 nspec2=nspec2+1
              endif
           endif
 1000   CONTINUE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9971, *)"beforenpart1 is ", npart1, "beforenpart2 is ",
     &       npart2
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        npart1=IHNT2(1)-nspec1
        npart2=IHNT2(3)-nspec2
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9971, *)"afternpart1 is ", npart1, "afternpart2 is ",
     &       npart2
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        return
        end
c=======================================================================
c     2/18/03 use PYTHIA to decay eta,rho,omega,k*,phi and Delta
c     4/2012 added pi0 decay flag: 
c       ipion=0: resonance or pi0 in lb(i1); >0: pi0 in lpion(ipion).
