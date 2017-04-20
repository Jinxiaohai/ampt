        subroutine zpca2
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
cc      SAVE /para3/
        common /para5/ iconfg, iordsc
cc      SAVE /para5/
        common /para7/ ioscar,nsmbbbar,nsmmeson
cc      SAVE /para7/
        common /ilist6/ t, iopern, icolln
cc      SAVE /ilist6/
        common /rndm1/ number
cc      SAVE /rndm1/
        common /rndm2/ iff
cc      SAVE /rndm2/
        common /rndm3/ iseedp
cc      SAVE /rndm3/
        COMMON /AREVT/ IAEVT, IARUN, MISS
cc      SAVE /AREVT/
        SAVE   
        if (iconfg .le. 3) then
           call zpca2a
        else
           call zpca2b
        end if
        if (ioscar .eq. 1) then
           call zpca2c
        end if
cbzdbg2/17/99
c        write (25, *) 'Event', nsevt - 1 + ievt, 
c    &         ', run', isbrun,
c        WRITE (25, *) ' Event ', IAEVT, ', run ', IARUN,
c     &     ',\n\t number of operations = ', iopern,
c     &     ',\n\t number of collisions between particles = ', 
c     &         icolln,
c     &     ',\n\t freezeout time=', t,
c     &     ',\n\t ending at the ', number, 'th random number',
c     &     ',\n\t ending collision iff=', iff
        WRITE (25, *) ' Event ', IAEVT, ', run ', IARUN
        WRITE (25, *) '    number of operations = ', iopern
        WRITE (25, *) '    number of collisions between particles = ', 
     &       icolln
        WRITE (25, *) '    freezeout time=', t
        WRITE (25, *) '    ending at the ', number, 'th random number'
        WRITE (25, *) '    ending collision iff=', iff
        return
        end
