c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      完成两粒子的散射。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        subroutine scat(t, iscat, jscat)
c       this subroutine is used to calculate the 2 particle scattering
        implicit double precision (a-h, o-z)
        SAVE   
        call newpos(t, iscat)
        call newpos(t, jscat)
        call newmom(t)
        return
        end
