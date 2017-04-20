        SUBROUTINE TITLE
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
        SAVE   
        WRITE(6,200)
clin-8/15/02 f77:
c200        FORMAT(//10X,
c     &        '**************************************************'/10X,
c     &  '*     |      \       _______      /  ------/     *'/10X,
c     &        '*   ----- ------     |_____|     /_/     /       *'/10X,
c     &        '*    ||\    /        |_____|      /    / \       *'/10X,
c     &        '*    /| \  /_/       /_______    /_  /    \_     *'/10X,
c     &        '*   / |     / /     /  /  / |        -------     *'/10X,
c     &        '*     |    / /\       /  /  |     /     |        *'/10X,
c     &        '*     |   / /  \     /  / \_|    /   -------     *'/10X,
200        FORMAT(//10X,
     &        '**************************************************'/10X,
     &  '*     |      |       _______      /  ------/     *'/10X,
     &        '*   ----- ------     |_____|     /_/     /       *'/10X,
     &        '*    |||    /        |_____|      /    / |       *'/10X,
     &        '*    /| |  /_/       /_______    /_  /    |      *'/10X,
     &        '*   / |     / /     /  /  / |        -------     *'/10X,
     &        '*     |    / /|       /  /  |     /     |        *'/10X,
     &        '*     |   / /  |     /  /  _|    /   -------     *'/10X,
     &        '*                                                *'/10X,
     &        '**************************************************'/10X,
     &        '                      HIJING                      '/10X,
     &        '       Heavy Ion Jet INteraction Generator        '/10X,
     &        '                        by                        '/10X,
     &  '            X. N. Wang  and  M. Gyulassy           '/10X,
     &  '             Lawrence Berkeley Laboratory           '//)        
        RETURN
        END
c.................... hipyset1.35.f
C
C
C
C     Modified for HIJING program
c
c    modification July 22, 1997  In pyremnn put an upper limit
c     on the total pt kick the parton can accumulate via multiple
C     scattering. Set the upper limit to be the sqrt(s)/2,
c     this is fix cronin bug for Pb+Pb events at SPS energy.
c
C
C Last modification Oct. 1993 to comply with non-vax
C machines' compiler 
C
C*********************************************************************  
