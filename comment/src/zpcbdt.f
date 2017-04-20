        block data zpcbdt
c       set initial values in block data
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        PARAMETER (MAXSTR=150001)
        common /para1/ mul
cc      SAVE /para1/
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
cc      SAVE /para2/
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
cc      SAVE /para3/
        common /para4/ iftflg, ireflg, igeflg, ibstfg
cc      SAVE /para4/
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        common /para6/ centy
cc      SAVE /para6/
clin-6/2009 nsmbbbar and nsmmeson respectively give the total number of 
c     baryons/anti-baryons and mesons for each event:
c        common /para7/ ioscar
        common /para7/ ioscar,nsmbbbar,nsmmeson
cc      SAVE /para7/
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
cc      SAVE /prec1/
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        common /prec3/gxs(MAXPTN),gys(MAXPTN),gzs(MAXPTN),fts(MAXPTN),
     &     pxs(MAXPTN), pys(MAXPTN), pzs(MAXPTN), es(MAXPTN),
     &     xmasss(MAXPTN), ityps(MAXPTN)
cc      SAVE /prec3/
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
cc      SAVE /prec4/
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
cc      SAVE /prec5/
        common /prec6/ etas(MAXPTN), raps(MAXPTN), taus(MAXPTN)
cc      SAVE /prec6/
        common /aurec1/ jxa, jya, jza
cc      SAVE /aurec1/
        common /aurec2/ dgxa(MAXPTN), dgya(MAXPTN), dgza(MAXPTN)
cc      SAVE /aurec2/
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        common /ilist2/ icell, icel(10,10,10)
cc      SAVE /ilist2/
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
cc      SAVE /ilist3/
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
cc      SAVE /ilist4/
c     6/07/02 initialize in ftime to expedite compiling:
c        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        common /ilist6/ t, iopern, icolln
cc      SAVE /ilist6/
        COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
cc      SAVE /ilist7/
        COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)
cc      SAVE /ilist8/
        common /rndm1/ number
cc      SAVE /rndm1/
        common /rndm2/ iff
cc      SAVE /rndm2/
        common /rndm3/ iseedp
cc      SAVE /rndm3/
        common /ana1/ ts(12)
cc      SAVE /ana1/
        common /ana2/
     &     det(12), dn(12), detdy(12), detdn(12), dndy(12),
     &     det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12),
     &     det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
cc      SAVE /ana2/
        common /ana3/ em(4, 4, 12)
cc      SAVE /ana3/
        common /ana4/ fdetdy(24), fdndy(24), fdndpt(12)
cc      SAVE /ana4/
        SAVE   
        data centy/0d0/
c     6/07/02 initialize in ftime to expedite compiling:
c        data (ct(i), i = 1, MAXPTN)/MAXPTN*0d0/
c        data (ot(i), i = 1, MAXPTN)/MAXPTN*0d0/
c        data tlarge/1000000.d0/
        data number/0/
        data ts/0.11d0, 0.12d0, 0.15d0, 0.2d0, 0.3d0, 0.4d0, 0.6d0,
     &     0.8d0, 1d0, 2d0, 4d0, 6d0/
c
        end
******************************************************************************
******************************************************************************
