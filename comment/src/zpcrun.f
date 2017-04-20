        subroutine zpcrun(*)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        parameter (tend1 = 250d0)
        parameter (tend2 = 6.1d0)
        common /para1/ mul
cc      SAVE /para1/
        common /para5/ iconfg, iordsc
        common /para7/ ioscar,nsmbbbar,nsmmeson
cc      SAVE /para5/
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
cc      SAVE /prec4/
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
cc      SAVE /prec5/
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
cc      SAVE /ilist4/
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        common /ilist6/ t, iopern, icolln
cc      SAVE /ilist6/
        common /ana1/ ts(12)
cc      SAVE /ana1/
        common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
        COMMON /AREVT/ IAEVT, IARUN, MISS
        SAVE   
c       save last collision info
        if (mod(ictype, 2) .eq. 0) then
           call savrec(iscat)
           call savrec(jscat)
        end if
c1      get operation type
        call getict(t1)
c2      check freezeout condition
        if (iconfg .eq. 1 .and. t1 .gt. tlarge / 2d0) return 1
        if (iconfg .eq. 2 .or. iconfg .eq. 3) then
           if (t1 .gt. tend1) return 1
c           if (ichkpt .eq. mul) then
c              ii = 0
c              do i = 1, mul
c                 gztemp = gz(i) + vz(i) * (t1 - ft(i))
c                 if (sqrt(t1 ** 2 - gztemp ** 2) .lt. tend) then
c                    ii = 1
c                    goto 1000
c                 end if
c              end do
c 1000              continue
c              if (ii .eq. 0) return 1
c           end if
        end if
        if (iconfg .eq. 4 .or. iconfg .eq. 5) then
           if (t1 .gt. tend2) return 1
        end if
clin-6/06/02 local freezeout for string melting,
c     decide what partons have frozen out at time t1:
      if(isoft.eq.5) then
         call local(t1)
      endif
c3      update iopern, t
        iopern = iopern + 1
        t = t1
        if (mod(ictype, 2) .eq. 0) then
           icolln = icolln + 1
c     4/18/01-ctest off
c           write (2006, 1233) 'iscat=', iscat, 'jscat=', jscat,
c           write (2006, *) 'iscat=', iscat, ' jscat=', jscat,
c     1 ityp(iscat), ityp(jscat)
c           write (2006, 1233) 'iscat=', max(indx(iscat), indx(jscat)),
c     &        'jscat=', min(indx(iscat), indx(jscat))
c           write (2006, 1234) ' icolln=', icolln, 't=', t
c 1233           format (a10, i10, a10, i10)
c 1234           format (a15, i10, a5, f23.17, a5, f23.17)
        end if
c4.1    deal with formation
        if (iconfg .eq. 1
     &     .or. iconfg .eq. 2
     &     .or. iconfg .eq. 4) then
           if (ictype .eq. 1 .or. ictype .eq. 2 .or. 
     &        ictype .eq. 5 .or. ictype .eq. 6) then
              call celasn
           end if
        end if
c4.2    deal with collisions
        if (ictype .ne. 1) then
           iscat0 = iscat
           jscat0 = jscat
c        iscat is the larger one so that if it's a wall collision,
c       it's still ok
           iscat = max0(iscat0, jscat0)
           jscat = min0(iscat0, jscat0)
ctest off check icsta(i): 0 with f77 compiler
c        write(9,*) 'BB:ictype,t1,iscat,jscat,icsta(i)=',
c     1 ictype,t1,iscat,jscat,icsta(iscat)
c       check collision time table error 'tterr'
clin-4/2008 to avoid out-of-bound error in next():
c           if (jscat .ne. 0 .and. next(jscat) .ne. iscat)
c     &        then
c              print *, 'iscat=', iscat, 'jscat=', jscat,
c     &             'next(', jscat, ')=', next(jscat)
c
c              if (ct(iscat) .lt. tlarge / 2d0) stop 'tterr'
c              if (ct(jscat) .lt. tlarge / 2d0) stop 'tterr'
c           end if 
           if (jscat .ne. 0) then
              if(next(jscat) .ne. iscat) then
                 print *, 'iscat=', iscat, 'jscat=', jscat,
     &                'next(', jscat, ')=', next(jscat)
                 if (ct(iscat) .lt. tlarge / 2d0) stop 'tterr'
                 if (ct(jscat) .lt. tlarge / 2d0) stop 'tterr'
              endif
           end if 
clin-4/2008-end
c4.2.1     collisions with wall
c     8/19/02 avoid actual argument in common blocks of cellre:
         niscat=iscat
         njscat=jscat
c           if (icsta(iscat) .ne. 0) call cellre(iscat, t)
c           if (jscat .ne. 0) then
c              if (icsta(jscat) .ne. 0) call cellre(jscat, t)
c           end if
           if (icsta(iscat) .ne. 0) call cellre(niscat, t)
           if (jscat .ne. 0) then
              if (icsta(jscat) .ne. 0) call cellre(njscat, t)
           end if
c4.2.2     collision between particles     
clin-6/2009 write out info for each collision:
c           if (mod(ictype, 2) .eq. 0) call scat(t, iscat, jscat)
           if (mod(ictype, 2) .eq. 0) then
              if(ioscar.eq.3) then
            write(95,*) 'event,miss,iscat,jscat=',iaevt,miss,iscat,jscat
                 if(dmax1(abs(gx(iscat)),abs(gy(iscat)),
     1                abs(gz(iscat)),abs(ft(iscat)),abs(gx(jscat)),
     2                abs(gy(jscat)),abs(gz(jscat)),abs(ft(jscat)))
     3                .lt.9999) then
                    write(95,200) ityp(iscat),px(iscat),py(iscat),
     1                   pz(iscat),xmass(iscat),gx(iscat),gy(iscat),
     2                   gz(iscat),ft(iscat)
                    write(95,200) ityp(jscat),px(jscat),py(jscat),
     1                   pz(jscat),xmass(jscat),gx(jscat),gy(jscat),
     2                   gz(jscat),ft(jscat)
                 else
                    write(95,201) ityp(iscat),px(iscat),py(iscat),
     1                   pz(iscat),xmass(iscat),gx(iscat),gy(iscat),
     2                   gz(iscat),ft(iscat)
                    write(95,201) ityp(jscat),px(jscat),py(jscat),
     1                   pz(jscat),xmass(jscat),gx(jscat),gy(jscat),
     2                   gz(jscat),ft(jscat)
                 endif
              endif
c     
              call scat(t, iscat, jscat)
c     
              if(ioscar.eq.3) then
                 if(dmax1(abs(gx(iscat)),abs(gy(iscat)),
     1                abs(gz(iscat)),abs(ft(iscat)),abs(gx(jscat)),
     2                abs(gy(jscat)),abs(gz(jscat)),abs(ft(jscat)))
     3                .lt.9999) then
                    write(95,200) ityp(iscat),px(iscat),py(iscat),
     1                   pz(iscat),xmass(iscat),gx(iscat),gy(iscat),
     2                   gz(iscat),ft(iscat)
                    write(95,200) ityp(jscat),px(jscat),py(jscat),
     1                   pz(jscat),xmass(jscat),gx(jscat),gy(jscat),
     2                   gz(jscat),ft(jscat)
                 else
                    write(95,201) ityp(iscat),px(iscat),py(iscat),
     1                   pz(iscat),xmass(iscat),gx(iscat),gy(iscat),
     2                   gz(iscat),ft(iscat)
                    write(95,201) ityp(jscat),px(jscat),py(jscat),
     1                   pz(jscat),xmass(jscat),gx(jscat),gy(jscat),
     2                   gz(jscat),ft(jscat)
                 endif
              endif
           endif
 200       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
 201       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
        end if
c5      update the interaction list
        call ulist(t)
c6      update ifmpt. ichkpt
c       old ichkpt and ifmpt are more conveniently used in ulist
        if (ifmpt .le. mul) then
           if (ictype .ne. 0 .and. ictype .ne. 3 
     &        .and. ictype .ne. 4) then
              ichkpt = ichkpt + 1
              ifmpt = ifmpt + 1
           end if
        end if
        return
        end
