c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  此函数生成一个时间信息进行返回。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        double precision function ftime1(iseed)
c       this program is used to generate formation time
c       the calling program needs a common /par1/
c       and declare external ftime1
clin-8/19/02
        implicit double precision (a-h, o-z)
        external ran1
        parameter (hbarc = 0.197327054d0)
        common /par1/ formt
cc      SAVE /par1/
        SAVE   
        aa = hbarc / formt
clin7/20/01:
c        ftime1 = aa * sqrt(1d0 / ran1(iseed) - 1d0)
        ftime1 = aa * dsqrt(1d0 / ran1(iseed) - 1d0)
        return
        end
