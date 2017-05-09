        subroutine zpca2a
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para1/ mul
cc      SAVE /para1/
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
cc      SAVE /para2/
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
cc      SAVE /para3/
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        common /para6/ centy
cc      SAVE /para6/
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
cc      SAVE /prec5/
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
cc      SAVE /ilist4/
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        common /ilist6/ t, iopern, icolln
cc      SAVE /ilist6/
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
        common /ana4/ fdetdy(24), fdndy(24), fdndpt(12)
cc      SAVE /ana4/
        SAVE   
        do 1004 i = 1, ichkpt
           rapi = rap(i)
c     7/20/01:
c           et = sqrt(px(i) ** 2 + py(i) ** 2 + xmp ** 2)
           et = dsqrt(px(i) ** 2 + py(i) ** 2 + xmp ** 2)

           do 1001 j = 1, 24
              if (rapi .gt. j + centy - 13d0 
     &           .and. rapi .lt. j  + centy - 12d0) then
                 fdetdy(j) = fdetdy(j) + et
                 fdndy(j) = fdndy(j) + 1d0
              end if
 1001      continue
           do 1002 j = 1, 12
              if (et .gt. 0.5d0 * (j - 1) .and.
     &           et .lt. 0.5d0 * j ) then
                 fdndpt(j) = fdndpt(j) + 1d0
              end if
 1002      continue
           if (iconfg .eq. 1) then
              t1 = ft(i)
              t2 = tlarge
              ipic = 11
           else
              t1 = tau(i)
              t2 = tlarge
              ipic = 12
           end if

           do 1003 ian = 1, ipic
              if (t1 .le. ts(ian) .and.
     &           t2 .gt. ts(ian)) then
                 call zpca1b(rapi, et, ian)
              end if
 1003      continue

           if (iconfg .eq. 1) then
              call zpca1b(rapi, et, 12)
           end if
 1004   continue

        do 1005 ian = 1, 12
           if (dn(ian) .eq. 0d0 .or. dn1(ian) .eq. 0d0 .or.
     &        dn2(ian) .eq. 0d0) then
clin-9/2012 suppress output:
c              print *, 'event=', ievt
c              print *, 'dn(', ian, ')=', dn(ian), 'dn1(', ian,
c     &           ')=', dn1(ian), 'dn2(', ian, ')=', dn2(ian)
           end if
           detdy(ian) = detdy(ian) + det(ian)
           if (dn(ian) .ne. 0) then
              detdn(ian) = detdn(ian) + det(ian) / dn(ian)
           end if
           dndy(ian) = dndy(ian) + dn(ian)
           detdy1(ian) = detdy1(ian) + det1(ian)
           if (dn1(ian) .ne. 0) then
              detdn1(ian) = detdn1(ian) + det1(ian) / dn1(ian)
           end if
           dndy1(ian) = dndy1(ian) + dn1(ian)
           detdy2(ian) = detdy2(ian) + det2(ian)
           if (dn2(ian) .ne. 0) then
              detdn2(ian) = detdn2(ian) + det2(ian) / dn2(ian)
           end if
           dndy2(ian) = dndy2(ian) + dn2(ian)
 1005   continue
        return
        end
