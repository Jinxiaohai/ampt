c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$     对粒子进行排序。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine arindx(n, m, arrin, indx)
cbz1/29/99end
c     indexes the first m elements of ARRIN of length n, i.e., outputs INDX
c     such that ARRIN(INDEX(J)) is in ascending order for J=1,...,m
c      implicit real*4 (a-h, o-z)
      dimension arrin(n), indx(n)
      SAVE   
      do 1001 j = 1, m
         indx(j) = j
 1001 continue
      l = m / 2 + 1
      ir = m
 10   continue
      if (l .gt. 1) then
         l = l - 1
         indxt = indx(l)
         q = arrin(indxt)
      else
         indxt = indx(ir)
         q = arrin(indxt)
         indx(ir) = indx(1)
         ir = ir - 1
         if (ir .eq. 1) then
            indx(1) = indxt
            return
         end if
      end if
      i = l
      j = l + l
 20   if (j .le. ir) then
         if (j .lt. ir) then
            if (arrin(indx(j)) .lt. arrin(indx(j + 1))) j = j + 1
         end if
         if (q .lt. arrin(indx(j))) then
            indx(i) = indx(j)
            i = j
            j = j + j
         else
            j = ir + 1
         end if
      goto 20
      end if
      indx(i) = indxt
      goto 10
      end
c-----------------------------------------------------------------------
c.....extracted from G. Song's ART expasion including K- interactions
c.....file `NEWKAON.FOR'
c     5/01/03 send iblock value into art1f.f, necessary for resonance studies:
c        subroutine newka(icase,irun,iseed,dt,nt,ictrl,i1,i2,
c     &                                   srt,pcx,pcy,pcz)
