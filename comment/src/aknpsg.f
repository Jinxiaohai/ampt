        real function akNPsg(pkaon)
      SAVE   
*cross section in mb for K- + N reactions.
c       sigma1: x section for K- + p/n -> sigma0 + PI0
         if(pkaon.le.0.345)then
           sigma1=0.624*pkaon**(-1.83)
         else
           sigma1=0.7*pkaon**(-2.09)
         endif
        akNPsg=sigma1
        return
        end   
c-----------------------------------------------------------------------
c.....extracted from G. Song's ART expasion including K- interactions
c.....file `NEWNNK.FOR'
