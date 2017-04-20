      subroutine pbarfs(srt,npion,iseed)
cbz2/25/99end
* Quantities: 
*  srt: DSQRT(s) in GeV                                                    *
*  npion: No. of pions produced in the annihilation of ppbar at srt        *
*  nmax=6, cutoff of the maximum no. of n the code can handle     
*                                             
*  Reference: C.M. Ko and R. Yuan, Phys. Lett. B192 (1987) 31      *
*
******************************************
       parameter (pimass=0.140,pi=3.1415926) 
       Dimension factor(6),pnpi(6) 
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      SAVE   
C the factorial coefficients in the pion no. distribution 
* from n=2 to 6 calculated use the formula in the reference
       factor(2)=1.
       factor(3)=1.17e-01
       factor(4)=3.27e-03
       factor(5)=3.58e-05
       factor(6)=1.93e-07
       ene=(srt/pimass)**3/(6.*pi**2)
c the relative probability from n=2 to 6
       do 1001 n=2,6 
           pnpi(n)=ene**n*factor(n)
 1001   continue
c find the maximum of the probabilities, I checked a 
c Fortan manual: max() returns the maximum value of 
c the same type as in the argument list
       pmax=max(pnpi(2),pnpi(3),pnpi(4),pnpi(5),pnpi(6))
c randomly generate n between 2 and 6
       ntry=0
 10    npion=2+int(5*RANART(NSEED))
clin-4/2008 check bounds:
       if(npion.gt.6) goto 10
       thisp=pnpi(npion)/pmax  
       ntry=ntry+1 
c decide whether to take this npion according to the distribution
c using rejection method.
       if((thisp.lt.RANART(NSEED)).and.(ntry.le.20)) go to 10
c now take the last generated npion and return
       return
       END
**********************************
cbali2/6/99 end
cbz3/9/99 kkbar
cbali3/5/99
******************************************
* purpose: Xsection for K+ K- to pi+ pi-
c      real*4 function xkkpi(srt)
*  srt    = DSQRT(s) in GeV                                  *
*  xkkpi   = xsection in mb obtained from
*           the detailed balance                             *
* ******************************************
c          parameter (pimass=0.140,aka=0.498)
c       xkkpi=1.e-08 
c       ppi2=(srt/2)**2-pimass**2
c       pk2=(srt/2)**2-aka**2
c       if(ppi2.le.0.or.pk2.le.0)return
cbz3/9/99 kkbar
c       xkkpi=ppi2/pk2*pipik(srt)
c       xkkpi=9.0 / 4.0 * ppi2/pk2*pipik(srt)
c        xkkpi = 2.0 * xkkpi
cbz3/9/99 kkbar end
cbz3/9/99 kkbar
c       end
c       return
c        END
cbz3/9/99 kkbar end
cbali3/5/99 end
cbz3/9/99 kkbar end
cbz3/9/99 kkbar
*****************************
* purpose: Xsection for K+ K- to pi+ pi-
