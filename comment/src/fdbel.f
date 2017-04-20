      function fdbel(s)
      parameter(srt0=2.012)
      if(s.le.srt0**2) then
         fdbel=0.
      else
         fdbel=2500.*exp(-(s-7.93)**2/0.003)
     1        +300.*exp(-(s-7.93)**2/0.1)+10.
      endif
      return
      end
c.................... hijing1.383_ampt.f
c     Version 1.383
c     The variables isng in HIJSFT and JL in ATTRAD were not initialized.
c     The version initialize them. (as found by Fernando Marroquim)
c
c
c
c     Version 1.382
c     Nuclear distribution for deuteron is taken as the Hulthen wave
c     function as provided by Brian Cole (Columbia)
clin     used my own implementation of impact parameter 
clin     & proton-neutron distance within a deuteron.
c
c
c     Version 1.381
c
c     The parameters for Wood-Saxon distribution for deuteron are
c     constrained to give the right rms ratius 2.116 fm
c     (R=0.0, D=0.5882)
c
c
c     Version 1.38
c
c     The following common block is added to record the number of elastic
c     (NELT, NELP) and inelastic (NINT, NINP) participants
c
c        COMMON/HJGLBR/NELT,NINT,NELP,NINP
c        SAVE /HJGLBR/
c
c     Version 1.37
c
c     A bug in the quenching subroutine is corrected. When calculating the
c     distance between two wounded nucleons, the displacement of the
c     impact parameter was not inculded. This bug was discovered by
c     Dr. V.Uzhinskii JINR, Dubna, Russia
c
c
C     Version 1.36
c
c     Modification Oct. 8, 1998. In hijing, log(ran(nseed)) occasionally
c     causes overfloat. It is modified to log(max(ran(nseed),1.0e-20)).
c
c
C     Nothing important has been changed here. A few 'garbage' has been
C     cleaned up here, like common block HJJET3 for the sea quark strings
C     which were originally created to implement the DPM scheme which
C     later was abadoned in the final version. The lines which operate
C     on these data are also deleted in the program.
C
C
C     Version 1.35
C     There are some changes in the program: subroutine HARDJET is now
C     consolidated with HIJHRD. HARDJET is used to re-initiate PYTHIA
C     for the triggered hard processes. Now that is done  altogether
C     with other normal hard processes in modified JETINI. In the new
C     version one calls JETINI every time one calls HIJHRD. In the new
C     version the effect of the isospin of the nucleon on hard processes,
C     especially direct photons is correctly considered.
C     For A+A collisions, one has to initilize pythia
C     separately for each type of collisions, pp, pn,np and nn,
C     or hp and hn for hA collisions. In JETINI we use the following
C     catalogue for different types of collisions:
C     h+h: h+h (itype=1)
C     h+A: h+p (itype=1), h+n (itype=2)
C     A+h: p+h (itype=1), n+h (itype=2)
C     A+A: p+p (itype=1), p+n (itype=2), n+p (itype=3), n+n (itype=4)
C*****************************************************************
c
C
C     Version 1.34
C     Last modification on January 5, 1998. Two mistakes are corrected in
C     function G. A Mistake in the subroutine Parton is also corrected.
C     (These are pointed out by Ysushi Nara).
C
C
C       Last modifcation on April 10, 1996. To conduct final
C       state radiation, PYTHIA reorganize the two scattered
C       partons and their final momenta will be a little
C       different. The summed total momenta of the partons
C       from the final state radiation are stored in HINT1(26-29)
C       and HINT1(36-39) which are little different from 
C       HINT1(21-24) and HINT1(41-44).
C
C       Version 1.33
C
C       Last modfication  on September 11, 1995. When HIJING and
C       PYTHIA are initialized, the shadowing is evaluated at
C       b=0 which is the maximum. This will cause overestimate
C       of shadowing for peripheral interactions. To correct this
C       problem, shadowing is set to zero when initializing. Then
C       use these maximum  cross section without shadowing as a
C       normalization of the Monte Carlo. This however increase
C       the computing time. IHNT2(16) is used to indicate whether
C       the sturcture function is called for (IHNT2(16)=1) initialization
C       or for (IHNT2(16)=0)normal collisions simulation
C
C       Last modification on Aagust 28, 1994. Two bugs associate
C       with the impact parameter dependence of the shadowing is
C       corrected.
C
C
c       Last modification on October 14, 1994. One bug is corrected
c       in the direct photon production option in subroutine
C       HIJHRD.( this problem was reported by Jim Carroll and Mike Beddo).
C       Another bug associated with keeping the decay history
C       in the particle information is also corrected.(this problem
C       was reported by Matt Bloomer)
C
C
C       Last modification on July 15, 1994. The option to trig on
C       heavy quark production (charm IHPR2(18)=0 or beauty IHPR2(18)=1) 
C       is added. To do this, set IHPR2(3)=3. For inclusive production,
C       one should reset HIPR1(10)=0.0. One can also trig larger pt
C       QQbar production by giving HIPR1(10) a nonvanishing value.
C       The mass of the heavy quark in the calculation of the cross
C       section (HINT1(59)--HINT1(65)) is given by HIPR1(7) (the
C       default is the charm mass D=1.5). We also include a separate
C       K-factor for heavy quark and direct photon production by
C       HIPR1(23)(D=2.0).
C
C       Last modification on May 24, 1994.  The option to
C       retain the information of all particles including those
C       who have decayed is IHPR(21)=1 (default=0). KATT(I,3) is 
C       added to contain the line number of the parent particle 
C       of the current line which is produced via a decay. 
C       KATT(I,4) is the status number of the particle: 11=particle
C       which has decayed; 1=finally produced particle.
C
C
C       Last modification on May 24, 1994( in HIJSFT when valence quark
C       is quenched, the following error is corrected. 1.2*IHNT2(1) --> 
C       1.2*IHNT2(1)**0.333333, 1.2*IHNT2(3) -->1.2*IHNT(3)**0.333333)
C
C
C       Last modification on March 16, 1994 (heavy flavor production
C       processes MSUB(81)=1 MSUB(82)=1 have been switched on,
C       charm production is the default, B-quark option is
C       IHPR2(18), when it is switched on, charm quark is 
C       automatically off)
C
C
C       Last modification on March 23, 1994 (an error is corrected
C       in the impact parameter dependence of the jet cross section)
C
C       Last modification Oct. 1993 to comply with non-vax
C       machines' compiler 
C
C*********************************************
C	LAST MODIFICATION April 5, 1991
CQUARK DISTRIBUTIOIN (1-X)**A/(X**2+C**2/S)**B 
C(A=HIPR1(44),B=HIPR1(46),C=HIPR1(45))
C STRING FLIP, VENUS OPTION IHPR2(15)=1,IN WHICH ONE CAN HAVE ONE AND
C TWO COLOR CHANGES, (1-W)**2,W*(1-W),W*(1-W),AND W*2, W=HIPR1(18), 
C AMONG PT DISTRIBUTION OF SEA QUARKS IS CONTROLLED BY HIPR1(42)
C
C	gluon jets can form a single string system
C
C	initial state radiation is included
C	
C	all QCD subprocesses are included
c
c	direct particles production is included(currently only direct
C		photon)
c
C	Effect of high P_T trigger bias on multiple jets distribution
c
C******************************************************************
C	                        HIJING.10                         *
C	          Heavy Ion Jet INteraction Generator        	  *
C	                           by                       	  *
C		   X. N. Wang      and   M. Gyulassy           	  *
C	 	      Lawrence Berkeley Laboratory		  *
C								  *
C******************************************************************
C
C******************************************************************
C NFP(K,1),NFP(K,2)=flavor of q and di-q, NFP(K,3)=present ID of  *
C proj, NFP(K,4) original ID of proj.  NFP(K,5)=colli status(0=no,*
C 1=elastic,2=the diffrac one in single-diffrac,3= excited string.*
C |NFP(K,6)| is the total # of jet production, if NFP(K,6)<0 it   *
C can not produce jet anymore. NFP(K,10)=valence quarks scattering*
C (0=has not been,1=is going to be, -1=has already been scattered *
C NFP(k,11) total number of interactions this proj has suffered   *
C PP(K,1)=PX,PP(K,2)=PY,PP(K,3)=PZ,PP(K,4)=E,PP(K,5)=M(invariant  *
C mass), PP(K,6,7),PP(K,8,9)=transverse momentum of quark and     *
C diquark,PP(K,10)=PT of the hard scattering between the valence  *
C quarks; PP(K,14,15)=the mass of quark,diquark.       		  * 
C******************************************************************
C
C****************************************************************
C
C	SUBROUTINE HIJING
C
C****************************************************************
