        real function akNel(pkaon)
      SAVE   
*cross section in mb for K- + N reactions.
c        the following data come from PRC 41 (1701)
c        sigma1: K(-) + neutron elastic
        if(pkaon.lt.0.5.or. pkaon.ge.4.0) sigma1=0.
        if(pkaon.ge.0.5.and.pkaon.lt.1.0) sigma1=20.*pkaon**2.74
        if(pkaon.ge.1.0.and.pkaon.lt.4.0) sigma1=20.*pkaon**(-1.8)
        akNel=sigma1
        return
        end
