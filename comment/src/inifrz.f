c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    此程序无调用
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        subroutine inifrz
c$$$c
c$$$        implicit double precision  (a-h, o-z)
c$$$        PARAMETER (MAXPTN=400001)
c$$$        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
c$$$cc      SAVE /ilist5/
c$$$        common /frzprc/ 
c$$$     &       gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
c$$$     &       pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
c$$$     &       xmfrz(MAXPTN), 
c$$$     &       tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
c$$$cc      SAVE /frzprc/
c$$$        SAVE   
c$$$c
c$$$c     for freezeout time 0-10fm, use interval of 0.1fm; 
c$$$c     for 10-100fm, use interval of 1fm; 
c$$$c     for 100-1000fm, use interval of 10fm; 
c$$$c     for 1000-3000fm, use interval of 100fm: 
c$$$        step1=0.1d0
c$$$        step2=1d0
c$$$        step3=10d0
c$$$        step4=100d0
c$$$c     
c$$$        do 1001 it=1,101
c$$$           tfrz(it)=0d0+dble(it-1)*step1
c$$$ 1001 continue
c$$$        do 1002 it=102,191
c$$$           tfrz(it)=10d0+dble(it-101)*step2
c$$$ 1002   continue
c$$$        do 1003 it=192,281
c$$$           tfrz(it)=100d0+dble(it-191)*step3
c$$$ 1003   continue
c$$$        do 1004 it=282,301
c$$$           tfrz(it)=1000d0+dble(it-281)*step4
c$$$ 1004   continue
c$$$        tfrz(302)=tlarge
c$$$c
c$$$        return
c$$$        end
c$$$clin-5/2009 v2 analysis
c$$$c=======================================================================
c$$$c     idd=0,1,2,3 specifies different subroutines for partonic flow analysis.
