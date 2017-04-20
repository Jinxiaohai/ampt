        subroutine readpa
        implicit double precision (a-h, o-z)
        external ran1
        character*50 str
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
cc      SAVE /para2/
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
cc      SAVE /para3/
        common /para4/ iftflg, ireflg, igeflg, ibstfg
cc      SAVE /para4/
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        common /para7/ ioscar,nsmbbbar,nsmmeson
cc      SAVE /para7/
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
cc      SAVE /ilist3/
        common /rndm1/ number
cc      SAVE /rndm1/
        common /rndm2/ iff
cc      SAVE /rndm2/
        common /rndm3/ iseedp
        common /xiaohai/xiaohaiflag
cc      SAVE /rndm3/
        SAVE   
        iseed=iseedp
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(9963,*) "call readpa end"
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c       this is the initialization file containing the initial values of 
c          the parameters
cbz1/31/99
c        open (5, file = 'zpc.ini', status = 'unknown')
cbz1/31/99end
c       this is the final data file containing general info about the cascade
cbz1/31/99
c        open (6, file = 'zpc.res', status = 'unknown')
        open (25, file = 'ana/zpc.res', status = 'unknown')
cbz1/31/99end
c       this is the input file containing initial particle records
cbz1/25/99
c        open (7, file = 'zpc.inp', status = 'unknown')
cbz1/25/99end
c       this gives the optional OSCAR standard output
cbz1/31/99
c        open (8, file = 'zpc.oscar', status = 'unknown')
        if(ioscar.eq.1) then
           open (26, file = 'ana/parton.oscar', status = 'unknown')
           open (19, file = 'ana/hadron.oscar', status = 'unknown')
        endif
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$         WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$          WW     WW WW     WW  RR  RR   II     TT     EE
c$$$           WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$            WW WW     WW WW    RR RR    II     TT     EE
c$$$             WW        WW      RR  RR   II     TT     EEEEEE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
           open(9997, file='testOut/parton.oscar', status = 'unknown')
           open(9996, file='testOut/hadron.oscar', status = 'unknown')
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$         WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$          WW     WW WW     WW  RR  RR   II     TT     EE
c$$$           WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$            WW WW     WW WW    RR RR    II     TT     EE
c$$$             WW        WW      RR  RR   II     TT     EEEEEE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
cbz1/31/99end
c     2/11/03 combine zpc initialization into ampt.ini:
c        open (29, file = 'zpc.ini', status = 'unknown')
c        read (29, *) str, xmp
        xmp=0d0
c        read (29, *) str, xmu
c        read (29, *) str, alpha
        cutof2 = 4.5d0 * (alpha / xmu) ** 2
c        read (29, *) str, rscut2
        rscut2=0.01d0
c        read (29, *) str, nsevt
        nsevt=1
c        read (29, *) str, nevnt
        nevnt=1
c        read (29, *) str, nsbrun
        nsbrun=1
c        read (29, *) str, iftflg
        iftflg=0
c        read (29, *) str, ireflg
        ireflg=1
cbz1/31/99
        IF (ireflg .EQ. 0) THEN
           OPEN (27, FILE = 'zpc.inp', STATUS = 'UNKNOWN')
        END IF
cbz1/31/99end
c        read (29, *) str, igeflg
        igeflg=0
c        read (29, *) str, ibstfg
        ibstfg=0
c        read (29, *) str, iconfg
        iconfg=1
c        read (29, *) str, iordsc
        iordsc=11
c        read (29, *) str, ioscar
c        read (29, *) str, v1, v2, v3
        v1=0.2d0
        v2=0.2d0
        v3=0.2d0
c        read (29, *) str, size1, size2, size3
        size1=1.5d0
        size2=1.5d0
        size3=0.7d0
        if (size1 .eq. 0d0 .or. size2 .eq. 0d0 .or. 
     &     size3 .eq. 0d0) then
           if (size1 .ne. 0d0 .or. size2 .ne. 0d0 .or. size3 .ne. 0d0
     &        .or. v1 .ne. 0d0 .or. v2 .ne. 0d0 .or. v3 .ne. 0d0) then
              print *, 'to get rid of space division:'
              print *, 'set all sizes and vs to 0'
              stop 'chker'
           end if
        end if
        size = min(size1, size2, size3)
c        read (29, *) str, iff
        iff=-1
c        read (29, *) str, iseed
c     10/24/02 get rid of argument usage mismatch in ran1():
        isedng=-iseed
c        a = ran1(-iseed)
        a = ran1(isedng)
c        read (29, *) str, irused
        irused=2
        do 1001 i = 1, irused - 1
c           a = ran1(2)
           iseed2=2
           a = ran1(iseed2)
 1001   continue
c     10/24/02-end
        if (iconfg .eq. 2 .or. iconfg .eq. 3) then
           v1 = 0d0
           v2 = 0d0
        end if
        if (iconfg .eq. 4 .or. iconfg .eq. 5) then
           v1 = 0d0
           v2 = 0d0
           v3 = 0d0
        end if
        close(5)
        return
        end
