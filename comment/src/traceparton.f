      subroutine traceparton
      implicit double precision (a-h, o-z)
      parameter (MAXPTN=400001)
      common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &     px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &     xmass(MAXPTN), ityp(MAXPTN)
      common /ilist6/ t, iopern, icolln
      common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
      common /para1/ mul
      common /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
      common /tracexiaohai/ tracetime(30), indextime(30)
      SAVE

      do 457 i=1, 30
         if(indextime(i) .eq. 0) then
            if(t .lt. tracetime(i+1) .and. t .gt. tracetime(i)) then
               write(9887,*)MUL
               do 455 ihai=1, mul
                  write(9887,200)ityp(ihai),px(ihai),py(ihai),
     &                 pz(ihai),xmass(ihai),gx(ihai),gy(ihai),
     &                 gz(ihai),ft(ihai), indx(ihai),LSTRG0(indx(ihai))
 455           continue
            indextime(i) = 1
            endif
         endif
 457  continue

 200  format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2), 2(2x,i8))
      end
