      function RANART(NSEED)
      SAVE   
clin-4/2008 ran(nseed) is renamed to avoid conflict with system functions:
c      ran=rand()
      ranart=rand(0)
c     one may also use the following random number generator in PYTHIA/JETSET:
c      ranart=rlu(0)
      return
      end
clin-3/2009
c     Initialize hadron weights; 
c     Can add initial hadrons before the hadron cascade starts (but after ZPC).
