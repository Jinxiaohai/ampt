      FUNCTION asinh(x)
      SAVE
      if(x.gt.0) then
         ASINH=alog(x+sqrt(x**2+1.))
      else
c     a la suggestion de YP Liu:
         ASINH=-alog(-x+sqrt(x**2+1.))
      endif
      return
      end
c....................art1f.f
**************************************
*
*                           PROGRAM ART1.0 
*
*        A relativistic transport (ART) model for heavy-ion collisions
*
*   sp/01/04/2002
*   calculates K+K- from phi decay, dimuons from phi decay
*   has finite baryon density & possibilites of varying Kaon 
*   in-medium mass in phiproduction-annhilation channel only.
*
*
* RELEASING DATE: JAN., 1997 
***************************************
* 
* Bao-An Li & Che Ming Ko
* Cyclotron Institute, Texas A&M University.
* Phone: (409) 845-1411
* e-mail: Bali@comp.tamu.edu & Ko@comp.tamu.edu 
* http://wwwcyc.tamu.edu/~bali
***************************************
* Speical notice on the limitation of the code:
* 
* (1) ART is a hadronic transport model
* 
* (2) E_beam/A <= 15 GeV
* 
* (3) The mass of the colliding system is limited by the dimensions of arrays
*    which can be extended purposely. Presently the dimensions are large enough
*     for running Au+Au at 15 GeV/A.
*
* (4) The production and absorption of antiparticles (e.g., ki-, anti-nucleons,
*     etc) are not fully included in this version of the model. They, however, 
*     have essentially no effect on the reaction dynamics and observables 
*     related to nucleons, pions and kaons (K+) at and below AGS energies.
* 
* (5) Bose enhancement for mesons and Pauli blocking for fermions are 
*     turned off.
* 
*********************************
*
* USEFUL REFERENCES ON PHYSICS AND NUMERICS OF NUCLEAR TRANSPORT MODELS:
*     G.F. BERTSCH AND DAS GUPTA, PHYS. REP. 160 (1988) 189.
*     B.A. LI AND W. BAUER, PHYS. REV. C44 (1991) 450.
*     B.A. LI, W. BAUER AND G.F. BERTSCH, PHYS. REV. C44 (1991) 2095.
*     P. DANIELEWICZ AND G.F. BERTSCH, NUCL. PHYS. A533 (1991) 712.
* 
* MAIN REFERENCES ON THIS VERSION OF ART MODEL:
*     B.A. LI AND C.M. KO, PHYS. REV. C52 (1995) 2037; 
*                          NUCL. PHYS. A601 (1996) 457. 
*
**********************************
**********************************
*  VARIABLES IN INPUT-SECTION:                                               * 
*                                                                      *
*  1) TARGET-RELATED QUANTITIES                                        *
*       MASSTA, ZTA -  TARGET MASS IN AMU, TARGET CHARGE  (INTEGER)    *
*                                                                      *
*  2) PROJECTILE-RELATED QUANTITIES                                    *
*       MASSPR, ZPR -  PROJECTILE MASS IN AMU, PROJ. CHARGE(INTEGER)   *
*       ELAB     -  BEAM ENERGY IN [MEV/NUCLEON]               (REAL)  *
*       ZEROPT   -  DISPLACEMENT OF THE SYSTEM IN Z-DIREC. [FM](REAL)  *
*       B        -  IMPACT PARAMETER [FM]                      (REAL)  *
*                                                                      *
*  3) PROGRAM-CONTROL PARAMETERS                                       *
*       ISEED    -  SEED FOR RANDOM NUMBER GENERATOR        (INTEGER)  *
*       DT       -  TIME-STEP-SIZE [FM/C]                      (REAL)  *
*       NTMAX    -  TOTAL NUMBER OF TIMESTEPS               (INTEGER)  *
*       ICOLL    -  (= 1 -> MEAN FIELD ONLY,                           *
*                -   =-1 -> CACADE ONLY, ELSE FULL ART)     (INTEGER)  *
*       NUM      -  NUMBER OF TESTPARTICLES PER NUCLEON     (INTEGER)  *
*       INSYS    -  (=0 -> LAB-SYSTEM, ELSE C.M. SYSTEM)    (INTEGER)  *
*       IPOT     -  1 -> SIGMA=2; 2 -> SIGMA=4/3; 3 -> SIGMA=7/6       *
*                   IN MEAN FIELD POTENTIAL                 (INTEGER)  *
*       MODE     -  (=1 -> interpolation for pauli-blocking,           *
*                    =2 -> local lookup, other -> unblocked)(integer)  *
*       DX,DY,DZ -  widths of cell for paulat in coor. sp. [fm](real)  *
*       DPX,DPY,DPZ-widths of cell for paulat in mom. sp.[GeV/c](real) *
*       IAVOID   -  (=1 -> AVOID FIRST COLL. WITHIN SAME NUCL.         *
*                    =0 -> ALLOW THEM)                      (INTEGER)  *
*       IMOMEN   -  FLAG FOR CHOICE OF INITIAL MOMENTUM DISTRIBUTION   *
*                   (=1 -> WOODS-SAXON DENSITY AND LOCAL THOMAS-FERMI  *
*                    =2 -> NUCLEAR MATTER DEN. AND LOCAL THOMAS-FERMI  *
*                    =3 -> COHERENT BOOST IN Z-DIRECTION)   (INTEGER)  *
*  4) CONTROL-PRINTOUT OPTIONS                                         *
*       NFREQ    -  NUMBER OF TIMSTEPS AFTER WHICH PRINTOUT            *
*                   IS REQUIRED OR ON-LINE ANALYSIS IS PERFORMED       *
*       ICFLOW      =1 PERFORM ON-LINE FLOW ANALYSIS EVERY NFREQ STEPS *
*       ICRHO       =1 PRINT OUT THE BARYON,PION AND ENERGY MATRIX IN  *
*                      THE REACTION PLANE EVERY NFREQ TIME-STEPS       *
*  5)
*       CYCBOX   -  ne.0 => cyclic boundary conditions;boxsize CYCBOX  *
*
**********************************
*               Lables of particles used in this code                     *
**********************************
*         
*         LB(I) IS USED TO LABEL PARTICLE'S CHARGE STATE
*    
*         LB(I)   =
clin-11/07/00:
*                -30 K*-
clin-8/29/00
*                -13 anti-N*(+1)(1535),s_11
*                -12 anti-N*0(1535),s_11
*                 -11 anti-N*(+1)(1440),p_11
*                 -10 anti-N*0(1440), p_11
*                  -9 anti-DELTA+2
*                  -8 anti-DELTA+1
*                  -7 anti-DELTA0
*                  -6 anti-DELTA-1
clin-8/29/00-end
cbali2/7/99 
*                  -2 antineutron 
*                             -1       antiproton
cbali2/7/99 end 
*                   0 eta
*                        1 PROTON
*                   2 NUETRON
*                   3 PION-
*                   4 PION0
*                   5 PION+
*                   6 DELTA-1
*                   7 DELTA0
*                   8 DELTA+1
*                   9 DELTA+2
*                   10 N*0(1440), p_11
*                   11 N*(+1)(1440),p_11
*                  12 N*0(1535),s_11
*                  13 N*(+1)(1535),s_11
*                  14 LAMBDA
*                   15 sigma-, since we used isospin averaged xsection for
*                   16 sigma0  sigma associated K+ production, sigma0 and 
*                   17 sigma+  sigma+ are counted as sigma-
*                   21 kaon-
*                   23 KAON+
*                   24 kaon0
*                   25 rho-
*                         26 rho0
*                   27 rho+
*                   28 omega meson
*                   29 phi
clin-11/07/00:
*                  30 K*+
* sp01/03/01
*                 -14 LAMBDA(bar)
*                  -15 sigma-(bar)
*                  -16 sigma0(bar)
*                  -17 sigma+(bar)
*                   31 eta-prime
*                   40 cascade-
*                  -40 cascade-(bar)
*                   41 cascade0
*                  -41 cascade0(bar)
*                   45 Omega baryon
*                  -45 Omega baryon(bar)
* sp01/03/01 end
clin-5/2008:
*                   42 Deuteron (same in ampt.dat)
*                  -42 anti-Deuteron (same in ampt.dat)
c
*                   ++  ------- SEE BAO-AN LI'S NOTE BOOK
**********************************
cbz11/16/98
c      PROGRAM ART
