c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$    此程序不参与调用，故全注释掉了。
c$$$  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        subroutine frztm(NEVNT,idd)
c$$$c
c$$$        implicit double precision  (a-h, o-z)
c$$$        PARAMETER (MAXPTN=400001)
c$$$        dimension tsf(31)
c$$$        COMMON /PARA1/ MUL
c$$$cc      SAVE /PARA1/
c$$$        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
c$$$     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
c$$$     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
c$$$cc      SAVE /prec1/
c$$$        COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
c$$$     &       PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
c$$$     &       XMASS5(MAXPTN), ITYP5(MAXPTN)
c$$$cc      SAVE /prec2/
c$$$        COMMON /frzout/ xnprod(30),etprod(30),xnfrz(30),etfrz(30),
c$$$     & dnprod(30),detpro(30),dnfrz(30),detfrz(30)
c$$$cc      SAVE /frzout/ 
c$$$        SAVE   
c$$$        data tsf/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
c$$$     &       1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
c$$$     &       2  , 3,   4,   5,   6,   7,   8,   9,   10,  20,  30/
c$$$c
c$$$        if(idd.eq.0) then
c$$$           do 1001 ii=1,30
c$$$              xnprod(ii)=0d0
c$$$              etprod(ii)=0d0
c$$$              xnfrz(ii)=0d0
c$$$              etfrz(ii)=0d0
c$$$              dnprod(ii)=0d0
c$$$              detpro(ii)=0d0
c$$$              dnfrz(ii)=0d0
c$$$              detfrz(ii)=0d0
c$$$ 1001      continue
c$$$           OPEN (86, FILE = 'ana1/production.dat', STATUS = 'UNKNOWN')
c$$$           OPEN (87, FILE = 'ana1/freezeout.dat', STATUS = 'UNKNOWN')
c$$$        else if(idd.eq.1) then
c$$$           DO 100 ip = 1, MUL
c$$$              do 1002 ii=1,30
c$$$                 eth0=dSQRT(PX0(ip)**2+PY0(ip)**2+XMASS0(ip)**2)
c$$$                 eth2=dSQRT(PX5(ip)**2+PY5(ip)**2+XMASS5(ip)**2)
c$$$c     total number and Et produced by time tsf(ii):
c$$$                 if (ft0(ip).lt.tsf(ii+1)) then
c$$$                    xnprod(ii)=xnprod(ii)+1d0
c$$$                    etprod(ii)=etprod(ii)+eth0
c$$$c     number and Et produced from time tsf(ii) to tsf(ii+1):
c$$$                    if (ft0(ip).ge.tsf(ii)) then
c$$$                       dnprod(ii)=dnprod(ii)+1d0
c$$$                       detpro(ii)=detpro(ii)+eth0
c$$$                    endif
c$$$                 endif
c$$$c     total number and Et freezed out by time tsf(ii):
c$$$                 if (FT5(ip).lt.tsf(ii+1)) then
c$$$                    xnfrz(ii)=xnfrz(ii)+1d0
c$$$                    etfrz(ii)=etfrz(ii)+eth2
c$$$c     number and Et freezed out from time tsf(ii) to tsf(ii+1):
c$$$                    if (FT5(ip).ge.tsf(ii)) then
c$$$                       dnfrz(ii)=dnfrz(ii)+1d0
c$$$                       detfrz(ii)=detfrz(ii)+eth2
c$$$                    endif
c$$$                 endif
c$$$ 1002         continue
c$$$ 100           continue
c$$$        else if(idd.eq.2) then
c$$$           write (86,*) '       t,       np,       dnp/dt,      etp '//
c$$$     1 ' detp/dt'
c$$$           write (87,*) '       t,       nf,       dnf/dt,      etf '//
c$$$     1 ' detf/dt'
c$$$           do 1003 ii=1,30
c$$$              xnp=xnprod(ii)/dble(NEVNT)
c$$$              xnf=xnfrz(ii)/dble(NEVNT)
c$$$              etp=etprod(ii)/dble(NEVNT)
c$$$              etf=etfrz(ii)/dble(NEVNT)
c$$$              dxnp=dnprod(ii)/dble(NEVNT)/(tsf(ii+1)-tsf(ii))
c$$$              dxnf=dnfrz(ii)/dble(NEVNT)/(tsf(ii+1)-tsf(ii))
c$$$              detp=detpro(ii)/dble(NEVNT)/(tsf(ii+1)-tsf(ii))
c$$$              detf=detfrz(ii)/dble(NEVNT)/(tsf(ii+1)-tsf(ii))
c$$$              write (86, 200) 
c$$$     1        tsf(ii+1),xnp,dxnp,etp,detp
c$$$              write (87, 200) 
c$$$     1        tsf(ii+1),xnf,dxnf,etf,detf
c$$$ 1003      continue
c$$$        endif
c$$$ 200    format(2x,f9.2,4(2x,f10.2))
c$$$c
c$$$        return
c$$$        end
c=======================================================================
clin-6/2009 write out initial minijet information 
c     before propagating to its formation time:
clin-2/2012:
c        subroutine minijet_out(BB)
