       subroutine ppxS(lb1,lb2,srt,ppsig,spprho,ipp)
* purpose: this subroutine gives the cross section for pion+pion 
*          elastic collision
* variables: 
*       input: lb1,lb2 and srt are the labels and srt for I1 and I2
*       output: ppsig: pp xsection
*               ipp: label for the pion+pion channel
*               Ipp=0 NOTHING HAPPEND 
*                  1 for Pi(+)+PI(+) DIRECT
*                   2     PI(+)+PI(0) FORMING RHO(+)
*                  3     PI(+)+PI(-) FORMING RHO(0)
*                   4     PI(0)+PI(O) DIRECT
*                  5     PI(0)+PI(-) FORMING RHO(-)
*                  6     PI(-)+PI(-) DIRECT
* reference: G.F. Bertsch, Phys. Rev. D37 (1988) 1202.
* date : Aug 29, 1994
*****************************
       parameter (amp=0.14,pi=3.1415926)
      SAVE   
       PPSIG=0.0
cbzdbg10/15/99
        spprho=0.0
cbzdbg10/15/99 end
       IPP=0
       IF(SRT.LE.0.3)RETURN
       q=sqrt((srt/2)**2-amp**2)
       esigma=5.8*amp
       tsigma=2.06*q
       erho=0.77
       trho=0.095*q*(q/amp/(1.+(q/erho)**2))**2
       esi=esigma-srt
       if(esi.eq.0)then
       d00=pi/2.
       go to 10
       endif
       d00=atan(tsigma/2./esi)
10       erh=erho-srt
       if(erh.eq.0.)then
       d11=pi/2.
       go to 20
       endif
       d11=atan(trho/2./erh)
20       d20=-0.12*q/amp
       s0=8.*pi*sin(d00)**2/q**2
       s1=8*pi*3*sin(d11)**2/q**2
       s2=8*pi*5*sin(d20)**2/q**2
c    !! GeV^-2 to mb
        s0=s0*0.197**2*10.
        s1=s1*0.197**2*10.
        s2=s2*0.197**2*10.
C       ppXS=s0/9.+s1/3.+s2*0.56
C       if(ppxs.le.0)ppxs=0.00001
       spprho=s1/2.
* (1) PI(+)+PI(+)
       IF(LB1.EQ.5.AND.LB2.EQ.5)THEN
       IPP=1
       PPSIG=S2
       RETURN
       ENDIF
* (2) PI(+)+PI(0)
       IF((LB1.EQ.5.AND.LB2.EQ.4).OR.(LB1.EQ.4.AND.LB2.EQ.5))THEN
       IPP=2
       PPSIG=S2/2.+S1/2.
       RETURN
       ENDIF
* (3) PI(+)+PI(-)
       IF((LB1.EQ.5.AND.LB2.EQ.3).OR.(LB1.EQ.3.AND.LB2.EQ.5))THEN
       IPP=3
       PPSIG=S2/6.+S1/2.+S0/3.
       RETURN
       ENDIF
* (4) PI(0)+PI(0)
       IF(LB1.EQ.4.AND.LB2.EQ.4)THEN
       IPP=4
       PPSIG=2*S2/3.+S0/3.
       RETURN
       ENDIF
* (5) PI(0)+PI(-)
       IF((LB1.EQ.4.AND.LB2.EQ.3).OR.(LB1.EQ.3.AND.LB2.EQ.4))THEN
       IPP=5
       PPSIG=S2/2.+S1/2.
       RETURN
       ENDIF
* (6) PI(-)+PI(-)
       IF(LB1.EQ.3.AND.LB2.EQ.3)THEN
       IPP=6
       PPSIG=S2
       ENDIF
       return
       end
**********************************
* elementary kaon production cross sections
*  from the CERN data book
*  date: Sept.2, 1994
*  for pp-->pLK+
c      real*4 function pplpk(srt)
