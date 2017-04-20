      REAL FUNCTION X4pi(SRT)
      SAVE   
*       CROSS SECTION FOR NN-->DD+rho PROCESS
* *****************************
       akp=0.498
       ak0=0.498
       ana=0.94
       ada=1.232
       al=1.1157
       as=1.1197
       pmass=0.9383
       ES=SRT
       IF(ES.LE.4)THEN
       X4pi=0.
       ELSE
* cross section for two resonance pp-->DD+DN*+N*N*
       xpp2pi=4.*x2pi(es)
* cross section for pp-->pp+spi
       xpp3pi=3.*(x3pi(es)+x33pi(es))
* cross section for pp-->pD+ and nD++
       pps1=sigma(es,1,1,0)+0.5*sigma(es,1,1,1)
       pps2=1.5*sigma(es,1,1,1)
       ppsngl=pps1+pps2+s1535(es)
* CROSS SECTION FOR KAON PRODUCTION from the four channels
* for NLK channel
       xk1=0
       xk2=0
       xk3=0
       xk4=0
       t1nlk=ana+al+akp
       t2nlk=ana+al-akp
       if(es.le.t1nlk)go to 333
       pmnlk2=(es**2-t1nlk**2)*(es**2-t2nlk**2)/(4.*es**2)
       pmnlk=sqrt(pmnlk2)
       xk1=pplpk(es)
* for DLK channel
       t1dlk=ada+al+akp
       t2dlk=ada+al-akp
       if(es.le.t1dlk)go to 333
       pmdlk2=(es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
       pmdlk=sqrt(pmdlk2)
       xk3=pplpk(es)
* for NSK channel
       t1nsk=ana+as+akp
       t2nsk=ana+as-akp
       if(es.le.t1nsk)go to 333
       pmnsk2=(es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
       pmnsk=sqrt(pmnsk2)
       xk2=ppk1(es)+ppk0(es)
* for DSK channel
       t1DSk=aDa+aS+akp
       t2DSk=aDa+aS-akp
       if(es.le.t1dsk)go to 333
       pmDSk2=(es**2-t1DSk**2)*(es**2-t2DSk**2)/(4.*es**2)
       pmDSk=sqrt(pmDSk2)
       xk4=ppk1(es)+ppk0(es)
* THE TOTAL KAON+ AND KAON0 PRODUCTION CROSS SECTION IS THEN
333       XKAON=3.*(xk1+xk2+xk3+xk4)
* cross section for pp-->DD+rho
       x4pi=pp1(es)-ppsngl-xpp2pi-xpp3pi-XKAON
       if(x4pi.le.0)x4pi=1.E-06
       ENDIF
       RETURN
       END
******************************************
* for pp-->inelastic
c      real*4 function pp1(srt)
