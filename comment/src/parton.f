        SUBROUTINE PARTON(F,X1,X2,QQ)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        REAL HIPR1(100),HINT1(100)
        COMMON/HPARNT/HIPR1,IHPR2(50),HINT1,IHNT2(50)
cc      SAVE /HPARNT/
        COMMON/NJET/N,ipcrs
cc      SAVE /NJET/
clin-7/2009:
        common/cmsflag/dshadow,ishadow
        DIMENSION F(2,7) 
        SAVE   
        DLAM=dble(HIPR1(15))
        Q0=dble(HIPR1(16))
        S=DLOG(DLOG(QQ/DLAM**2)/DLOG(Q0**2/DLAM**2))
        IF(IHPR2(7).EQ.2) GO TO 200
C*******************************************************
        AT1=0.419d0+0.004d0*S-0.007d0*S**2
        AT2=3.460d0+0.724d0*S-0.066d0*S**2
        GMUD=4.40d0-4.86d0*S+1.33d0*S**2
        AT3=0.763d0-0.237d0*S+0.026d0*S**2
        AT4=4.00d0+0.627d0*S-0.019d0*S**2
        GMD=-0.421d0*S+0.033d0*S**2
C*******************************************************
        CAS=1.265d0-1.132d0*S+0.293d0*S**2
        AS=-0.372d0*S-0.029d0*S**2
        BS=8.05d0+1.59d0*S-0.153d0*S**2
        APHS=6.31d0*S-0.273d0*S**2
        BTAS=-10.5d0*S-3.17d0*S**2
        GMS=14.7d0*S+9.80d0*S**2
C********************************************************
C        CAC=0.135*S-0.075*S**2
C        AC=-0.036-0.222*S-0.058*S**2
C        BC=6.35+3.26*S-0.909*S**2
C        APHC=-3.03*S+1.50*S**2
C        BTAC=17.4*S-11.3*S**2
C        GMC=-17.9*S+15.6*S**2
C***********************************************************
        CAG=1.56d0-1.71d0*S+0.638d0*S**2
        AG=-0.949d0*S+0.325d0*S**2
        BG=6.0d0+1.44d0*S-1.05d0*S**2
        APHG=9.0d0-7.19d0*S+0.255d0*S**2
        BTAG=-16.5d0*S+10.9d0*S**2
        GMG=15.3d0*S-10.1d0*S**2
        GO TO 300
C********************************************************
200        AT1=0.374d0+0.014d0*S
        AT2=3.33d0+0.753d0*S-0.076d0*S**2
        GMUD=6.03d0-6.22d0*S+1.56d0*S**2
        AT3=0.761d0-0.232d0*S+0.023d0*S**2
        AT4=3.83d0+0.627d0*S-0.019d0*S**2
        GMD=-0.418d0*S+0.036d0*S**2
C************************************
        CAS=1.67d0-1.92d0*S+0.582d0*S**2
        AS=-0.273d0*S-0.164d0*S**2
        BS=9.15d0+0.530d0*S-0.763d0*S**2
        APHS=15.7d0*S-2.83d0*S**2
        BTAS=-101.0d0*S+44.7d0*S**2
        GMS=223.0d0*S-117.0d0*S**2
C*********************************
C        CAC=0.067*S-0.031*S**2
C        AC=-0.120-0.233*S-0.023*S**2
C        BC=3.51+3.66*S-0.453*S**2
C        APHC=-0.474*S+0.358*S**2
C        BTAC=9.50*S-5.43*S**2
C        GMC=-16.6*S+15.5*S**2
C**********************************
        CAG=0.879d0-0.971d0*S+0.434d0*S**2
        AG=-1.16d0*S+0.476d0*S**2
        BG=4.0d0+1.23d0*S-0.254d0*S**2
        APHG=9.0d0-5.64d0*S-0.817d0*S**2
        BTAG=-7.54d0*S+5.50d0*S**2
        GMG=-0.596d0*S+1.26d0*S**2
C*********************************
300        B12=DEXP(GMRE(AT1)+GMRE(AT2+1.D0)-GMRE(AT1+AT2+1.D0))
        B34=DEXP(GMRE(AT3)+GMRE(AT4+1.D0)-GMRE(AT3+AT4+1.D0))
        CNUD=3.D0/B12/(1.D0+GMUD*AT1/(AT1+AT2+1.D0))
        CND=1.D0/B34/(1.D0+GMD*AT3/(AT3+AT4+1.D0))
C********************************************************
C        FUD=X*(U+D)
C        FS=X*2(UBAR+DBAR+SBAR)  AND UBAR=DBAR=SBAR
C*******************************************************
        FUD1=CNUD*X1**AT1*(1.D0-X1)**AT2*(1.D0+GMUD*X1)
        FS1=CAS*X1**AS*(1.D0-X1)**BS*(1.D0+APHS*X1
     &      +BTAS*X1**2+GMS*X1**3)
        F(1,3)=CND*X1**AT3*(1.D0-X1)**AT4*(1.D0+GMD*X1)+FS1/6.D0
        F(1,1)=FUD1-F(1,3)+FS1/3.D0
        F(1,2)=FS1/6.D0
        F(1,4)=FS1/6.D0
        F(1,5)=FS1/6.D0
        F(1,6)=FS1/6.D0
        F(1,7)=CAG*X1**AG*(1.D0-X1)**BG*(1.D0+APHG*X1
     &         +BTAG*X1**2+GMG*X1**3)
C
        FUD2=CNUD*X2**AT1*(1.D0-X2)**AT2*(1.D0+GMUD*X2)
        FS2=CAS*X2**AS*(1.D0-X2)**BS*(1.D0+APHS*X2
     &      +BTAS*X2**2+GMS*X2**3)
        F(2,3)=CND*X2**AT3*(1.D0-X2)**AT4*(1.D0+GMD*X2)+FS2/6.D0
        F(2,1)=FUD2-F(2,3)+FS2/3.D0
        F(2,2)=FS2/6.D0
        F(2,4)=FS2/6.D0
        F(2,5)=FS2/6.D0
        F(2,6)=FS2/6.D0
        F(2,7)=CAG*X2**AG*(1.D0-X2)**BG*(1.D0+APHG*X2
     &         +BTAG*X2**2+GMG*X2**3)
C***********Nuclear effect on the structure function****************
C
        IF(IHPR2(6).EQ.1 .AND. IHNT2(1).GT.1) THEN
           AAX=1.193d0*dble(ALOG(FLOAT(IHNT2(1)))**0.16666666)
           RRX=AAX*(X1**3-1.2d0*X1**2+0.21d0*X1)+1.d0
     &               +dble(1.079*(FLOAT(IHNT2(1))**0.33333333-1.0))
     &          /dble(ALOG(float(IHNT2(1))+1.0))*DSQRT(X1)
     &          *DEXP(-X1**2/0.01d0)
clin-8/2015 DEXP() above may cause IEEE_UNDERFLOW, 
c     does not seem to affect results.
c     &          /DLOG(IHNT2(1)+1.0D0)*(DSQRT(X1)*DEXP(-X1**2/0.01)
clin-7/2009 enable users to modify nuclear shadowing:
           if(ishadow.eq.1) RRX=1.d0+dshadow*(RRX-1.d0)
           IF(ipcrs.EQ.1 .OR.ipcrs.EQ.3) RRX=DEXP(-X1**2/0.01d0)
clin-7/2009:
           if((ipcrs.EQ.1.OR.ipcrs.EQ.3).and.ishadow.eq.1) 
     1          RRX=DEXP(-X1**2/0.01d0)*dshadow
           DO 400 I=1,7
              F(1,I)=RRX*F(1,I)
 400           CONTINUE
        ENDIF
        IF(IHPR2(6).EQ.1 .AND. IHNT2(3).GT.1) THEN
           AAX=1.193d0*dble(ALOG(FLOAT(IHNT2(3)))**0.16666666)
           RRX=AAX*(X2**3-1.2d0*X2**2+0.21d0*X2)+1.d0
     &               +dble(1.079*(FLOAT(IHNT2(3))**0.33333-1.0))
     &          /dble(ALOG(float(IHNT2(3))+1.0))*DSQRT(X2)
     &          *DEXP(-X2**2/0.01d0)
c     &         /DLOG(IHNT2(3)+1.0D0)*DSQRT(X2)*DEXP(-X2**2/0.01)
clin-7/2009:
           if(ishadow.eq.1) RRX=1.d0+dshadow*(RRX-1.d0)
           IF(ipcrs.EQ.2 .OR. ipcrs.EQ.3) RRX=DEXP(-X2**2/0.01d0)
clin-7/2009:
           if((ipcrs.EQ.2.OR.ipcrs.EQ.3).and.ishadow.eq.1) 
     1          RRX=DEXP(-X2**2/0.01d0)*dshadow
           DO 500 I=1,7
              F(2,I)=RRX*F(2,I)
 500           CONTINUE
        ENDIF
c
        RETURN
        END
C
C
C
