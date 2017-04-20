        subroutine zpcou1
        implicit double precision (a-h, o-z)
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
cc      SAVE /para3/
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
c
        dpt = 0.5d0
        dy2 = 1d0
        dy1 = 0.5d0
        dy = 0.2d0
        ntotal = nevnt * nsbrun
c
        return
        end
