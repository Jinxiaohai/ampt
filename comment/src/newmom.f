        subroutine newmom(t)
c       this subroutine is used to calculate the 2 particle scattering
        implicit double precision (a-h, o-z)
        parameter (hbarc = 0.197327054d0)
        parameter (MAXPTN=400001)
        parameter (pi = 3.14159265358979d0)
        COMMON /para1/ mul
cc      SAVE /para1/
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
cc      SAVE /para2/
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
ctrans
        common /para6/ centy
cc      SAVE /para6/
ctransend
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
cc      SAVE /prec2/
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
cc      SAVE /prec4/
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
cc      SAVE /prec5/
        common /aurec1/ jxa, jya, jza
cc      SAVE /aurec1/
        common /aurec2/ dgxa(MAXPTN), dgya(MAXPTN), dgza(MAXPTN)
cc      SAVE /aurec2/
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
cc      SAVE /ilist3/
        common /lor/ enenew, pxnew, pynew, pznew
cc      SAVE /lor/
        common /cprod/ xn1, xn2, xn3
cc      SAVE /cprod/
        common /rndm2/ iff
cc      SAVE /rndm2/
        common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
        common /frzprc/ 
     &       gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &       pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &       xmfrz(MAXPTN), 
     &       tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
cc      SAVE /frzprc/
        SAVE   
clin-6/06/02 no momentum change for partons already frozen out,
c     however, spatial upgrade is needed to ensure overall system freezeout:
      if(isoft.eq.5) then
         if(ifrz(iscat).eq.1.or.ifrz(jscat).eq.1) then
            last(iscat) = jscat
            last(jscat) = iscat
            return
         endif
      endif
clin-6/06/02-end
c       iff is used to randomize the interaction to have both attractive and
c        repulsive
        iff = - iff
        if (iconfg .eq. 2 .or. iconfg .eq. 4) then
           icels1 = icels(iscat)
           i1 = icels1 / 10000
           j1 = (icels1 - i1 * 10000) / 100
           icels2 = icels(jscat)
           i2 = icels2 / 10000
           j2 = (icels2 - i2 * 10000) / 100
           if (iconfg .eq. 4) then
              k1 = icels1 - i1 * 10000 - j1 * 100
              k2 = icels2 - i2 * 10000 - j2 * 100
           end if
        end if
        px1 = px(iscat)
        py1 = py(iscat)
        pz1 = pz(iscat)
        e1 = e(iscat)
        x1 = gx(iscat)
        y1 = gy(iscat)
        z1 = gz(iscat)
        t1 = ft(iscat)
        px2 = px(jscat)
        py2 = py(jscat)
        pz2 = pz(jscat)
        e2 = e(jscat)
        if (iconfg .eq. 1) then
           x2 = gx(jscat)
           y2 = gy(jscat)
           z2 = gz(jscat)
        else if (iconfg .eq. 2 .or. iconfg .eq. 4) then
           if (i1 - i2 .gt. 5) then
              x2 = gx(jscat) + 10d0 * size1
           else if (i1 - i2 .lt. -5) then
              x2 = gx(jscat) - 10d0 * size1
           else
              x2 = gx(jscat)
           end if
           if (j1 - j2 .gt. 5) then
              y2 = gy(jscat) + 10d0 * size2
           else if (j1 - j2 .lt. -5) then
              y2 = gy(jscat) - 10d0 * size2
           else
              y2 = gy(jscat)
           end if
           if (iconfg .eq. 4) then
              if (k1 - k2 .gt. 5) then
                 z2 = gz(jscat) + 10d0 * size3
              else if (k1 - k2 .lt. -5) then
                 z2 = gz(jscat) - 10d0 * size3
              else
                 z2 = gz(jscat)
              end if
           else
              z2 = gz(jscat)
           end if
        else if (iconfg .eq. 3 .or. iconfg .eq. 5) then
           x2 = gx(jscat) + dgxa(jscat)
           y2 = gy(jscat) + dgya(jscat)
           if (iconfg .eq. 5) then
              z2 = gz(jscat) + dgza(jscat)
           else
              z2 = gz(jscat)
           end if
        end if
        t2 = ft(jscat)
ctrans
        rts2 = (e1 + e2) ** 2 - (px1 + px2) ** 2 -
     &     (py1 + py2) ** 2 - (pz1 + pz2) ** 2
ctransend
        bex = (px1 + px2) / (e1 + e2)
        bey = (py1 + py2) / (e1 + e2)
        bez = (pz1 + pz2) / (e1 + e2)
clin-11/2015-ctest off
c        write(99,*) 'iscat,jscat,etotalA=',iscat,jscat,e1+e2
        call lorenz(e1, px1, py1, pz1, bex, bey, bez)
cc      SAVE pxnew, ..., values for later use.
        px1 = pxnew
        py1 = pynew
        pz1 = pznew
        e1 = enenew
        pp2 = pxnew ** 2 + pynew ** 2 + pznew ** 2
        call getht(iscat, jscat, pp2, that)
        theta = dacos(that / (2d0 * pp2) + 1d0)
        theta = dble(iff) * theta
c       we boost to the cm frame, get rotation axis, and rotate 1 particle 
c       momentum
        call lorenz(t1, x1, y1, z1, bex, bey, bez)
        x1 = pxnew
        y1 = pynew
        z1 = pznew
        call lorenz(t2, x2, y2, z2, bex, bey, bez)
        x2 = pxnew
        y2 = pynew
        z2 = pznew
c       notice now pxnew, ..., are new positions
        call cropro(x1-x2, y1-y2, z1-z2, px1, py1, pz1)
        call xnormv(xn1, xn2, xn3)
cbz1/29/99
c        call rotate(xn1, xn2, xn3, theta, px1, py1, pz1)
        call zprota(xn1, xn2, xn3, theta, px1, py1, pz1)
cbz1/29/99end
c       we invert the momentum to get the other particle's momentum
        px2 = -px1
        py2 = -py1
        pz2 = -pz1
clin-4/13/01: modify in case m1, m2 are different:
c        e2 = e1
        e2 = dsqrt(px2**2+py2**2+pz2**2+xmass(jscat)**2)
clin-11/2015-ctest off
c        write(99,*) 'iscat,jscat,masses= ',iscat,jscat,
c     1       xmass(iscat),xmass(jscat)
c       boost the 2 particle 4 momentum back to lab frame
        call lorenz(e1, px1, py1, pz1, -bex, -bey, -bez)
        px(iscat) = pxnew
        py(iscat) = pynew
        pz(iscat) = pznew
        e(iscat) = enenew
        call lorenz(e2, px2, py2, pz2, -bex, -bey, -bez)        
        px(jscat) = pxnew
        py(jscat) = pynew
        pz(jscat) = pznew
        e(jscat) = enenew
clin-11/2015-ctest off
c        write(99,*) 'iscat,jscat,etotalB= ',iscat,jscat,
c     1       e(iscat)+e(jscat)
        vx(iscat) = px(iscat) / e(iscat)
        vy(iscat) = py(iscat) / e(iscat)
        vz(iscat) = pz(iscat) / e(iscat)
        vx(jscat) = px(jscat) / e(jscat)
        vy(jscat) = py(jscat) / e(jscat)
        vz(jscat) = pz(jscat) / e(jscat)
        last(iscat) = jscat
        last(jscat) = iscat
        if (iconfg .le. 3) then
           if (e(iscat) .le. abs(pz(iscat))) then
              rap(iscat) = 1000000.d0
           else
              rap(iscat) = 0.5d0 * log((e(iscat) + pz(iscat)) /
     &           (e(iscat) - pz(iscat)))
           end if
           if (e(jscat) .le. abs(pz(jscat))) then
              rap(jscat) = 1000000.d0
           else
              rap(jscat) = 0.5d0 * log((e(jscat) + pz(jscat)) /
     &           (e(jscat) - pz(jscat)))
           end if
ctrans
           rap1 = rap(iscat)
           rap2 = rap(jscat)
           if ((rap1 .lt. centy + 0.5d0 .and.
     &        rap1 .gt. centy - 0.5d0)) then
c              write (9, *) sqrt(ft(iscat) ** 2 - gz(iscat) ** 2), rts2
           end if
           if ((rap2 .lt. centy + 0.5d0 .and.
     &        rap2 .gt. centy - 0.5d0)) then
c              write (9, *) sqrt(ft(jscat) ** 2 - gz(jscat) ** 2), rts2
           end if
ctransend
        end if
clin-11/2015-ctest off
c        write(99,*) 'iscat,jscat,xmp,xmu,that=',iscat,jscat,xmp,xmu,that
        return
        end
