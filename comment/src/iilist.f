        subroutine iilist
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para1/ mul
cc      SAVE /para1/
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        common /ilist2/ icell, icel(10,10,10)
cc      SAVE /ilist2/
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
cc      SAVE /ilist4/
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        common /ilist6/ t, iopern, icolln
cc      SAVE /ilist6/
        SAVE   
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /I_list1/ ISCAT,JSCAT,NEXT(NMAXGL), LAST(NMAXGL),
c$$$        ICTYPES, ICSTA(NMAXGL), NIC(NMAXGL), ICELSTA(NMAXGL)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        ISCAT:particle index. if the operation involves only one
c$$$        particle, then iscat is the particle index, if it involves 2
c$$$        particles, iscat is the larger particle index.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        JSCAT:particle index. if the operation involves only one
c$$$        particle, then iscat is 0, if it involves 2
c$$$        particles, iscat is the smaller particle index.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NEXT(I):the next operation partner of particle I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        LAST(I):the last operation partner of particle I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        ICTYPE:the operation type.
c$$$        = 0: a collision between particles.
c$$$        = 1: the formation of a particle.
c$$$        = 2: both a collision between particles and the formation of a
c$$$             particle.
c$$$        = 3: a collision with wall.
c$$$        = 4: both a collision between particles and a wall collision.
c$$$        = 5: both a wall collision and a formation.
c$$$        = 6: a formation, collision between particles and a wall collision
c$$$             at the same time.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        ICSTA(I):the operation type for particle I.
c$$$        = 0: an ordinary collision.
c$$$        = 101:a collision with the wall with larger x.
c$$$        = 102:a collision with the wall with smaller x.
c$$$        = 103:a collision with the wall with larger y.
c$$$        = 104:a collision with the wall with smaller y.
c$$$        = 105:a collision with the wall with larger z.
c$$$        = 106:a collision with the wall with smaller z.
c$$$        = 111:a collision with another particle and the wall with lartger x
c$$$        = 112:a collision with another particle and the wall with smaller x
c$$$        = 113:a collision with another particle and the wall with lartger y
c$$$        = 114:a collision with another particle and the wall with lartger y
c$$$        = 115:a collision with another particle and the wall with lartger z
c$$$        = 116:a collision with another particle and the wall with lartger z
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NIC(I):the next particle index in the same cell as I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        ICELSTA(I):the encoded information of the cell number
c$$$        particle I is in. If particle is in cell (i1, i2, i3), then it
c$$$        equals i1 * 10000 + i2 * 100 + i3 for particles inside the cube.
c$$$        when a particle is outside the 10 by 10 by 10 box. its value is 111111.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        iscat = MAXPTN
        jscat = MAXPTN
        do 1001 i = 1, mul
           next(i) = 0
           last(i) = 0
           icsta(i) = 0
           nic(i) = 0
           icels(i) = 0
 1001   continue
        icell = 0
        do 1004 i1 = 1, 10
           do 1003 i2 = 1, 10
              do 1002 i3 = 1, 10
                 icel(i1, i2, i3) = 0
 1002         continue
 1003      continue
 1004   continue
        ichkpt = 0
        ifmpt = 1
        do 1005 i = 1, mul
           ct(i) = tlarge
           ot(i) = tlarge
 1005   continue
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  ct(i),ot(i)全为0.
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO 478 ihai = 1, MUL
           write(9925,*)ct(i), ot(i)
 478       continue
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$  WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$  WW     WW WW     WW  RR  RR   II     TT     EE
c$$$  WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$  WW WW     WW WW    RR RR    II     TT     EE
c$$$  WW        WW      RR  RR   II     TT     EEEEEE
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
        iopern = 0
        icolln = 0
        t = 0.d0
        return
        end
