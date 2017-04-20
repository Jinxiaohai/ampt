        subroutine inievt
        implicit double precision (a-h, o-z)
        COMMON /para1/ mul
cc      SAVE /para1/
        common /para4/ iftflg, ireflg, igeflg, ibstfg
cc      SAVE /para4/
        SAVE   
cbz1/25/99
c        mul = 0
cbz1/25/99
        if (ireflg .eq. 0) call readi
        if (igeflg .ne. 0) call genei
        if (ibstfg .ne. 0) call boosti
        return
        end
