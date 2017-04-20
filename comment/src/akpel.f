        real function akPel(pkaon)
      SAVE   
*cross section in mb for K- + N reactions.
c        the following data come from PRC 41 (1701)
c        sigma2: K(-) + proton elastic
        if(pkaon.lt.0.25.or. pkaon.ge.4.0) sigma2=0.
        if(pkaon.ge.0.25.and.pkaon.lt.4.0) sigma2=13.*pkaon**(-0.9)
        akPel=sigma2
        return
        end
