      subroutine hbtout(nnew,nt,ntmax)
c
      PARAMETER  (MAXSTR=150001,MAXR=1)
clin-5/2008 give tolerance to regular particles (perturbative probability 1):
      PARAMETER  (oneminus=0.99999,oneplus=1.00001)
      dimension lastkp(MAXSTR), newkp(MAXSTR),xnew(3)
      common /para7/ ioscar,nsmbbbar,nsmmeson
cc      SAVE /para7/
      COMMON/hbt/lblast(MAXSTR),xlast(4,MAXSTR),plast(4,MAXSTR),nlast
cc      SAVE /hbt/
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      COMMON   /AA/  R(3,MAXSTR)
cc      SAVE /AA/
      COMMON   /BB/  P(3,MAXSTR)
cc      SAVE /BB/
      COMMON   /CC/  E(MAXSTR)
cc      SAVE /CC/
      COMMON   /EE/  ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      common /lastt/itimeh,bimp
cc      SAVE /lastt/
      COMMON/tdecay/tfdcy(MAXSTR),tfdpi(MAXSTR,MAXR),tft(MAXSTR)
cc      SAVE /tdecay/
      COMMON /AREVT/ IAEVT, IARUN, MISS
cc      SAVE /AREVT/
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
cc      SAVE /snn/
      COMMON/HJGLBR/NELT,NINTHJ,NELP,NINP
cc      SAVE /HJGLBR/
      COMMON/FTMAX/ftsv(MAXSTR),ftsvt(MAXSTR, MAXR)
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
clin-12/14/03:
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      EXTERNAL IARFLV, INVFLV
      common /para8/ idpert,npertd,idxsec
clin-2/2012:
      common /phiHJ/iphirp,phiRP
      SAVE   
c
      do 1001 i=1,max0(nlast,nnew)
         lastkp(i)=0
 1001 continue
      do 1002 i=1,nnew
         newkp(i)=0
 1002 continue
c     for each of the particles, search the freezeout record (common /hbt/) 
c     to find & keep those which do not have interactions during this timestep:
      do 100 ip=1,nnew
         do 1004 iplast=1,nlast
            if(p(1,ip).eq.plast(1,iplast).and.
     1           p(2,ip).eq.plast(2,iplast).and.
     2           p(3,ip).eq.plast(3,iplast).and.
     3           e(ip).eq.plast(4,iplast).and.
     4           lb(ip).eq.lblast(iplast).and.
     5      dpertp(ip).eq.dplast(iplast).and.lastkp(iplast).eq.0) then
clin-5/2008 modified below to the above in case we have perturbative particles:
c     5           lastkp(iplast).eq.0) then
               deltat=nt*dt-xlast(4,iplast)
               ene=sqrt(plast(1,iplast)**2+plast(2,iplast)**2
     1              +plast(3,iplast)**2+plast(4,iplast)**2)
c     xnew gives the coordinate if a particle free-streams to current time:
               do 1003 ii=1,3
                  xnew(ii)=xlast(ii,iplast)+plast(ii,iplast)/ene*deltat
 1003          continue
                  dr=sqrt((r(1,ip)-xnew(1))**2+(r(2,ip)-xnew(2))**2
     1              +(r(3,ip)-xnew(3))**2)
c     find particles with dp=0 and dr<0.01, considered to be those 
c     without any interactions during this timestep, 
c     thus keep their last positions and time:
               if(dr.le.0.01) then
                  lastkp(iplast)=1
                  newkp(ip)=1
c                  if(lb(ip).eq.41) then
c                write(95,*) 'nt,ip,px,x=',nt,ip,p(1,ip),r(1,ip),ftsv(ip)
c                write(95,*) 'xnew=',xnew(1),xnew(2),xnew(3),xlast(4,ip)
c                  endif
clin-5/2009 Take care of formation time of particles read in at nt=ntmax-1:
                  if(nt.eq.ntmax.and.ftsv(ip).gt.((ntmax-1)*dt)) 
     1                 xlast(4,iplast)=ftsv(ip)
                  goto 100
               endif
            endif
 1004    continue
 100  continue
c     for current particles with interactions, fill their current info in 
c     the freezeout record (if that record entry needs not to be kept):
      do 150 ip=1,nnew
         if(newkp(ip).eq.0) then
            do 1005 iplast=1,nnew
               if(lastkp(iplast).eq.0) then
ctest off: write collision info
c                  if(lb(ip).eq.41) then
c                     write(95,*) 'nt,lb(ip)=',nt,lb(ip)
c                  write(95,*) '  last p=',plast(1,iplast),
c     1 plast(2,iplast),plast(3,iplast),plast(4,iplast)
c                  write(95,*) '  after p=',p(1,ip),p(2,ip),p(3,ip),e(ip)
c                  write(95,*) 'after x=',r(1,ip),r(2,ip),r(3,ip),ftsv(ip)
c                  endif
c
                  xlast(1,iplast)=r(1,ip)
                  xlast(2,iplast)=r(2,ip)
                  xlast(3,iplast)=r(3,ip)
                  xlast(4,iplast)=nt*dt
c
                  if(nt.eq.ntmax) then
c     freezeout time for decay daughters at the last timestep 
c     needs to include the decay time of the parent:
                     if(tfdcy(ip).gt.(ntmax*dt+0.001)) then
                        xlast(4,iplast)=tfdcy(ip)
c     freezeout time for particles unformed at the next-to-last timestep 
c     needs to be their formation time instead of (ntmax*dt):
                     elseif(ftsv(ip).gt.((ntmax-1)*dt)) then
                        xlast(4,iplast)=ftsv(ip)
                     endif
                  endif
                  plast(1,iplast)=p(1,ip)
                  plast(2,iplast)=p(2,ip)
                  plast(3,iplast)=p(3,ip)
                  plast(4,iplast)=e(ip)
                  lblast(iplast)=lb(ip)
                  lastkp(iplast)=1
clin-5/2008:
                  dplast(iplast)=dpertp(ip)
                  goto 150
               endif
 1005       continue
         endif
 150  continue
c     if the current particle list is shorter than the freezeout record,
c     condense the last-collision record by filling new record from 1 to nnew, 
c     and label these entries as keep:
      if(nnew.lt.nlast) then
         do 170 iplast=1,nlast
            if(lastkp(iplast).eq.0) then
               do 1006 ip2=iplast+1,nlast
                  if(lastkp(ip2).eq.1) then
                     xlast(1,iplast)=xlast(1,ip2)
                     xlast(2,iplast)=xlast(2,ip2)
                     xlast(3,iplast)=xlast(3,ip2)
                     xlast(4,iplast)=xlast(4,ip2)
                     plast(1,iplast)=plast(1,ip2)
                     plast(2,iplast)=plast(2,ip2)
                     plast(3,iplast)=plast(3,ip2)
                     plast(4,iplast)=plast(4,ip2)
                     lblast(iplast)=lblast(ip2)
                     lastkp(iplast)=1
clin-5/2008:
                     dplast(iplast)=dplast(ip2)
                     goto 170
                  endif
 1006          continue
            endif
 170     continue
      endif
      nlast=nnew
ctest off look inside each NT timestep (for debugging purpose):
c      do ip=1,nlast
c         write(99,*) ' p ',nt,ip,lblast(ip),plast(1,ip),
c     1        plast(2,ip),plast(3,ip),plast(4,ip),dplast(ip)
c         write(99,*) '  x ',nt,ip,lblast(ip),xlast(1,ip),
c     1        xlast(2,ip),xlast(3,ip),xlast(4,ip),dplast(ip)
c      enddo
c
      if(nt.eq.ntmax) then
clin-5/2008 find final number of perturbative particles (deuterons only):
         ndpert=0
         do ip=1,nlast
            if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
            else
               ndpert=ndpert+1
            endif
         enddo
c
c         write(16,190) IAEVT,IARUN,nlast,bimp,npart1,npart2,
c     1 NELP,NINP,NELT,NINTHJ
clin-2/2012:
c         write(16,190) IAEVT,IARUN,nlast-ndpert,bimp,npart1,npart2,
c     1 NELP,NINP,NELT,NINTHJ
         write(16,191) IAEVT,IARUN,nlast-ndpert,bimp,npart1,npart2,
     1 NELP,NINP,NELT,NINTHJ,phiRP
clin-5/2008 write out perturbatively-produced particles (deuterons only):
         if(idpert.eq.1.or.idpert.eq.2)
     1        write(90,190) IAEVT,IARUN,ndpert,bimp,npart1,npart2,
     2        NELP,NINP,NELT,NINTHJ
         do 1007 ip=1,nlast
clin-12/14/03   No formation time for spectator projectile or target nucleons,
c     see ARINI1 in 'amptsub.f':
clin-3/2009 To be consistent with new particles produced in hadron cascade
c     that are limited by the time-resolution (DT) of the hadron cascade, 
c     freezeout time of spectator projectile or target nucleons is written as 
c     DT as they are read at the 1st timestep and then propagated to time DT: 
c
clin-9/2011 determine spectator nucleons consistently
c            if(plast(1,ip).eq.0.and.plast(2,ip).eq.0
c     1           .and.(sqrt(plast(3,ip)**2+plast(4,ip)**2)*2/HINT1(1))
c     2           .gt.0.99.and.(lblast(ip).eq.1.or.lblast(ip).eq.2)) then
            if(abs(plast(1,ip)).le.epsiPt.and.abs(plast(2,ip)).le.epsiPt
     1           .and.(plast(3,ip).gt.amax1(0.,PZPROJ-epsiPz)
     2                .or.plast(3,ip).lt.(-PZTARG+epsiPz))
     3           .and.(lblast(ip).eq.1.or.lblast(ip).eq.2)) then
clin-5/2008 perturbatively-produced particles (currently only deuterons) 
c     are written to ana/ampt_pert.dat (without the column for the mass); 
c     ana/ampt.dat has regularly-produced particles (including deuterons);
c     these two sets of deuteron data are close to each other(but not the same 
c     because of the bias from triggering the perturbative production); 
c     ONLY use one data set for analysis to avoid double-counting:
               if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
                  write(16,200) INVFLV(lblast(ip)), plast(1,ip),
     1                 plast(2,ip),plast(3,ip),plast(4,ip),
     2                 xlast(1,ip),xlast(2,ip),xlast(3,ip),
     3                 xlast(4,ip)
clin-12/14/03-end
               else
                  if(idpert.eq.1.or.idpert.eq.2) then
                     write(90,250) INVFLV(lblast(ip)), plast(1,ip),
     1                 plast(2,ip),plast(3,ip),
     2                 xlast(1,ip),xlast(2,ip),xlast(3,ip),
     3                 xlast(4,ip)
                  else
                     write(99,*) 'Unexpected perturbative particles'
                  endif
               endif
            elseif(amax1(abs(xlast(1,ip)),abs(xlast(2,ip)),
     1              abs(xlast(3,ip)),abs(xlast(4,ip))).lt.9999) then
               if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
            write(16,200) INVFLV(lblast(ip)), plast(1,ip),
     1           plast(2,ip),plast(3,ip),plast(4,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip)
               else
                  if(idpert.eq.1.or.idpert.eq.2) then
            write(90,250) INVFLV(lblast(ip)),plast(1,ip),
     1           plast(2,ip),plast(3,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip),
     3           dplast(ip)
                  else
                     write(99,*) 'Unexpected perturbative particles'
                  endif
               endif
            else
c     change format for large numbers:
               if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
            write(16,201) INVFLV(lblast(ip)), plast(1,ip),
     1           plast(2,ip),plast(3,ip),plast(4,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip)
               else
                  if(idpert.eq.1.or.idpert.eq.2) then
                     write(90,251) INVFLV(lblast(ip)), plast(1,ip),
     1           plast(2,ip),plast(3,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip),
     3           dplast(ip)
                  else
                     write(99,*) 'Unexpected perturbative particles'
                  endif
               endif
            endif
 1007    continue
         if(ioscar.eq.1) call hoscar
      endif
 190  format(3(i7),f10.4,5x,6(i4))
 191  format(3(i7),f10.4,5x,6(i4),5x,f7.4)
clin-3/2009 improve the output accuracy of Pz
 200  format(I6,2(1x,f8.3),1x,f11.4,1x,f6.3,4(1x,f8.2))
 201  format(I6,2(1x,f8.3),1x,f11.4,1x,f6.3,4(1x,e8.2))
 250  format(I5,2(1x,f8.3),1x,f10.3,2(1x,f7.1),1x,f8.2,1x,f7.2,1x,e10.4)
 251  format(I5,2(1x,f8.3),1x,f10.3,4(1x,e8.2),1x,e10.4)
c     
        return
        end
c=======================================================================
