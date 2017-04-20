        subroutine inirun
        SAVE   
c       sort prec2 according to increasing formation time
        call ftime
        call inirec
        call iilist
        call inian2
        return
        end
