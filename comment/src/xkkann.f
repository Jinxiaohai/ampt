      SUBROUTINE XKKANN(SRT, XSK1, XSK2, XSK3, XSK4, XSK5,
     &     XSK6, XSK7, XSK8, XSK9, XSK10, XSK11, SIGK, rrkk)
*  srt    = DSQRT(s) in GeV                                       *
*  xsk1   = annihilation into pi pi                               *
*  xsk2   = annihilation into pi rho (shifted to XKKSAN)         *
*  xsk3   = annihilation into pi omega (shifted to XKKSAN)       *
*  xsk4   = annihilation into pi eta                              *
*  xsk5   = annihilation into rho rho                             *
*  xsk6   = annihilation into rho omega                           *
*  xsk7   = annihilation into rho eta (shifted to XKKSAN)        *
*  xsk8   = annihilation into omega omega                         *
*  xsk9   = annihilation into omega eta (shifted to XKKSAN)      *
*  xsk10  = annihilation into eta eta                             *
*  sigk   = xsection in mb obtained from                          *
*           the detailed balance                                  *
* ***************************
      PARAMETER  (MAXSTR=150001, MAXX=20,  MAXZ=24)
          PARAMETER (AKA=0.498, PIMASS=0.140, RHOM = 0.770, 
     &     OMEGAM = 0.7819, ETAM = 0.5473, APHI=1.02)
      COMMON  /AA/ R(3,MAXSTR)
cc      SAVE /AA/
      COMMON /BB/  P(3,MAXSTR)
cc      SAVE /BB/
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
cc      SAVE /DD/
      SAVE   
        S = SRT ** 2
       SIGK = 1.E-08
        XSK1 = 0.0
        XSK2 = 0.0
        XSK3 = 0.0
        XSK4 = 0.0
        XSK5 = 0.0
        XSK6 = 0.0
        XSK7 = 0.0
        XSK8 = 0.0
        XSK9 = 0.0
        XSK10 = 0.0
        XSK11 = 0.0
        XPION0 = PIPIK(SRT)
c.....take into account both K+ and K0
        XPION0 = 2.0 * XPION0
        PI2 = S * (S - 4.0 * AKA ** 2)
         if(PI2 .le. 0.0)return
        XM1 = PIMASS
        XM2 = PIMASS
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK1 = 9.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
clin-8/28/00 (pi eta) eta -> K+K- is assumed the same as pi pi -> K+K-:
        XM1 = PIMASS
        XM2 = ETAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK4 = 3.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XM1 = ETAM
        XM2 = ETAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK10 = 1.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XPION0 = rrkk
clin-11/07/00: (pi eta) (rho omega) -> K* Kbar (or K*bar K) instead to K Kbar:
c        XM1 = PIMASS
c        XM2 = RHOM
c        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
c        IF (PF2 .GT. 0.0) THEN
c           XSK2 = 27.0 / 4.0 * PF2 / PI2 * XPION0
c        END IF
c        XM1 = PIMASS
c        XM2 = OMEGAM
c        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
c        IF (PF2 .GT. 0.0) THEN
c           XSK3 = 9.0 / 4.0 * PF2 / PI2 * XPION0
c        END IF
        XM1 = RHOM
        XM2 = RHOM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK5 = 81.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XM1 = RHOM
        XM2 = OMEGAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK6 = 27.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
c        XM1 = RHOM
c        XM2 = ETAM
c        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
c        IF (PF2 .GT. 0.0) THEN
c           XSK7 = 9.0 / 4.0 * PF2 / PI2 * XPION0
c        END IF
        XM1 = OMEGAM
        XM2 = OMEGAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK8 = 9.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
c        XM1 = OMEGAM
c        XM2 = ETAM
c        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
c        IF (PF2 .GT. 0.0) THEN
c           XSK9 = 3.0 / 4.0 * PF2 / PI2 * XPION0
c        END IF
c* K+ + K- --> phi
          fwdp = 1.68*(aphi**2-4.*aka**2)**1.5/6./aphi/aphi     
clin-9/2012: check argument in sqrt():
          scheck=srt**2-4.0*aka**2
          if(scheck.le.0) then
             write(99,*) 'scheck47: ', scheck
             stop
          endif
          pkaon=0.5*sqrt(scheck)
c          pkaon=0.5*sqrt(srt**2-4.0*aka**2)
          XSK11 = 30.*3.14159*0.1973**2*(aphi*fwdp)**2/
     &             ((srt**2-aphi**2)**2+(aphi*fwdp)**2)/pkaon**2
c
        SIGK = XSK1 + XSK2 + XSK3 + XSK4 + XSK5 + 
     &     XSK6 + XSK7 + XSK8 + XSK9 + XSK10 + XSK11
       RETURN
        END
cbz3/9/99 kkbar end
*****************************
* purpose: Xsection for Phi + B 
