        subroutine zpca1
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        SAVE   
        if (mod(ictype,2) .eq. 0) then
           call zpca1a(iscat)
           call zpca1a(jscat)
clin-5/2009 ctest off v2 for parton:
c           call flowp(1)
        end if
        return
        end
