       real function ang(srt,iseed)
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
c        if(srt.le.2.14)then
c       b1s=0.5
c       b2s=0.
c      endif
      if((srt.gt.2.14).and.(srt.le.2.4))then
       b1s=29.03-23.75*srt+4.865*srt**2
         b2s=-30.33+25.53*srt-5.301*srt**2
      endif
      if(srt.gt.2.4)then
       b1s=0.06
         b2s=0.4
      endif
        x=RANART(NSEED)
       p=b1s/b2s
       q=(2.*x-1.)*(b1s+b2s)/b2s
       IF((-q/2.+sqrt((q/2.)**2+(p/3.)**3)).GE.0.)THEN
       ang1=(-q/2.+sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
       ELSE
       ang1=-(q/2.-sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
       ENDIF
       IF((-q/2.-sqrt((q/2.)**2+(p/3.)**3).GE.0.))THEN
       ang2=(-q/2.-sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
       ELSE
       ang2=-(q/2.+sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
       ENDIF
       ANG=ANG1+ANG2
       return
       end
*--------------------------------------------------------------------------
*****subprogram * kaon production from pi+B collisions *******************
