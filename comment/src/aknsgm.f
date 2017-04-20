        real function akNsgm(pkaon)
      SAVE   
*cross section in mb for K- + N reactions.
c        sigma2: x section for K- + n -> sigma0 + PI-
        if(pkaon.lt.0.5.or. pkaon.ge.6.0) sigma2=0.
        if(pkaon.ge.0.5.and.pkaon.lt.1.0) sigma2=1.2*pkaon**(-1.3)
        if(pkaon.ge.1.0.and.pkaon.lt.6.0) sigma2=1.2*pkaon**(-2.3)
        akNsgm=sigma2
        return
        end
