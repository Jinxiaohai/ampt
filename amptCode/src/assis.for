c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        声明四个数组变量HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50),
c$$$        来储存几乎所有的物理量。
c$$$        Purpose: contains input parameters (HIPR1,IHPR2) for event
c$$$        options and some extra information(HINT1,IHNT2) of current event
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(1):(D=1.5GeV/c2)minimum value for the invariant mass
c$$$        of the excited string system in a hadron-hadron interaction.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(2):(D=0.35GeV)width of the Gaussian Pt disstribution
c$$$        of produced hadron in Lund string fragmentation.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(3), HIPR1(4): (D=0.5, 0.9 GeV −2 ) give the a and b
c$$$        parameters of the symmetric Lund fragmentation
c$$$        function (PARJ(41), PARJ(42) in JETSET 7.2).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(5): (D=2.0 GeV/c 2 ) invariant mass cut-off for
c$$$        the dipole radiation of a string system
c$$$        below which soft gluon radiations are terminated.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(6): (D=0.1) the depth of shadowing of structure
c$$$        functions at x = 0 as defined in
c$$$        α A = HIPR1(6) × (A 1/3 − 1).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(8): (D=2.0 GeV/c) minimum P T transfer in hard or
c$$$        semihard scatterings.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(9): (D=−1.0 GeV/c) maximum P T transfer in hard
c$$$        or semihard scatterings. If
c$$$        negative, the limit is set by the colliding energy.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(10): (D=−2.25 GeV/c) specifies the value of P T for
c$$$        each triggered hard scattering generated per event
c$$$        (see Section ??). If HIPR1(10) is negative,
c$$$        its absolute value gives the low limit of the P T of the
c$$$        triggered jets.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(11): (D=2.0 GeV/c) minimum P T of a jet which will
c$$$        interact with excited nuclear matter. When the Pt
c$$$        of a jet is smaller than HIPR1(11) it will stop interacting
c$$$        further.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(12): (D=1.0 fm) transverse distance between a
c$$$        traversing jet and an excited nucleon (string system)
c$$$        below which they will interact and the jet
c$$$        will lose energy and momentum to that string system.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(13): (D=1.0 fm) the mean free path of a jet when
c$$$        it goes through the excited nuclear matter.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(14): (D=2.0 GeV/fm) the energy loss dE/dz of a
c$$$        gluon jet inside the excited nuclear
c$$$        matter. The energy loss for a quark jet is half
c$$$        of the energy loss of a gluon.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(15): (D=0.2 GeV/c) the scale Λ in the calculation of
c$$$        α s .
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(19),HIPR1(20):parameters in the distribution for the
c$$$        PT kick from soft interactions.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(30): the inclusive cross section sigma for soft
c$$$        interactions. the default value is used to ensure the
c$$$        geometrical scaling of PP interaction cross sections at low
c$$$        eneries
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(31): the cross section sigma which characterizes the
c$$$        geometrical size of a nucleon(此处看文献上的等式).The default
c$$$        value is only for high-energy limit > 200GeV. At lower energies,
c$$$        a slight decrease which depends on energy is parametrized in the
c$$$        program. The default values of the two parameters HIPR1(30),
c$$$        HIPR1(31) are only for NN type interactions.for other kinds of
c$$$        projectile or target hadrons, users should change these values
c$$$        so that correct inelastic and total cross sections are obtained
c$$$        by the program.(散射截面和核子大小的关系).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(34): maximum radial coordinate for projectile nucleons
c$$$        to be given by the initialization program HIJSET.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(35): maximum radial coordiante for target nucleons
c$$$        to be given by the initialization program HIJSET.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(40): value of PION.(pion 的具体数值).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(43):(D=0.01)fractional energy error relative to the
c$$$        colliding energy permitted per nucleon-necleon collision.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(44),HIPR1(45),HIPR1(46):(D=1.5, 0.1GeV, 0.25)
c$$$        parameters alpha , c and beta in the valence quark distributions
c$$$        for soft string excitation.具体参看第20页的内容。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HIPR1(47),HIPR1(48):(D=0.0, 0.5)parameters alpha and beta
c$$$        in valence quark distribution,() for the disassociated
c$$$        excitation in a single diffractive collision.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<










        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(1):(D=1)switch for dipole-approximated QCD radiation
c$$$        of the string system in soft interactions.(偶极近似的开关)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(2):(D=3)option for initial and final state radiation
c$$$        in the hard scattering.
c$$$        = 0 : both initial and final radiation are off.
c$$$        = 1 : initial radiation on and final radiation off.
c$$$        = 2 : initial radiation off and final radiation on.
c$$$        = 3 : both initial and final radiation are on.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(3):(D=3)switch for triggered hard scattering with
c$$$        specified Pt >= HIPR1(10).
c$$$        = 0 : no triggered jet production.
c$$$        = 1 : ordinary hard processes
c$$$        = 2 : only direct photon production
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(4):(D=1)switch for jet quenching in the excited
c$$$        nuclear matter.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(5):(D=1)switch for the PT kick due to soft interactions.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(6):(D=1)switch for the nuclear effect ona the parton
c$$$        distribution function such as shadowing.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(7):(D=1)selection of Duke-Owens set (1 or 2)of
c$$$        parameterization of nucleon structure functions.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(8):(D=10)maximum number of hard scattering per nucleon
c$$$        nucleon interaction. when IHPR2(8)=0, jet production will be
c$$$        turned off. when IHPR2(8)<0, the number of jet production will
c$$$        be fixed at its absolute value for each NN collision.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(9):(D=0)switch to guarantee at least one pair of
c$$$        minijets production per event (pp, pA or AB)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(10):(D=0)option to print warning messages about errors
c$$$        that might happen. when a fatal error happens the current event
c$$$        will abandoned and a new one is generated.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(10):(D=0)option to print waring message about errors
c$$$        that might happen. when a fatal error happens the current
c$$$        event will be abandoned and a new one is generated.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(11):(D=1)choice of baryon production model.
c$$$        = 0 : no baryon-antibaryon pair production, initial diquark
c$$$        treated as a unit.
c$$$        = 1 : diquark-antidiquark pair production allowed, initial
c$$$        diquark treated as a unit.
c$$$        = 2 : diquark-antidiquark pair production allowed, with the
c$$$        possibility for diquark to split according to the "popcorn"
c$$$        scheme(选择重子)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(12):(D=1) option to turn off the automatic decay of
c$$$        the following particles : pion0, K, D, sigma.(关掉衰变选项)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(13):(D=1)option to turn on single diffractive
c$$$        reactions.(单衍射反应????????????)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(14):(D=1)option to turn on elastic scattering.(打开
c$$$        弹性散射)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(18):(D=0)option to switch on B-quark production. charm
c$$$        production is the default.when B-quark production is on, charm
c$$$        quark production is automatically off.(B-quark产生的开关)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(19):(D=1)option to turn on initial state soft
c$$$        interaction.(打开软相互作用)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(20):(D=1)switch for the final fragmentation.(末态碎片
c$$$        的开关)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHPR2(21):(D=0)option to keep the information of all
c$$$        particles including those which have decayed and the decay
c$$$        history in the common block HIMAIN2. the line number of the
c$$$        parent particle is KATT(I,3).the status of a particle, whether
c$$$        it is a finally produced particle (KATT(I,4)=1) or a decayed
c$$$        particle (KATT(I,4)=11) is also kept.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



        


c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(1):colliding energy in the c.m frame of
c$$$        nucleon-nucleon collisions.(质心系的能量)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(2):Lorentz trandformation varialble Beta from
c$$$        laboratory to c.m.(洛伦兹变化的Beta)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(3):rapidity y of the c.m frame(质心系的快度).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(4):弹核的快度.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(5):靶核的快度.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(6):弹核的能量.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(7):靶核的能量.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(8):弹核的静止质量.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(9):靶核的静止质量.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(10):核核碰撞JET产生的平均截面.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(11):the averaged inclusive cross section for Jet 的产
c$$$        生再核核碰撞中.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(12):核核碰撞的平均非弹性散射截面.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(13):核核碰撞的平均总的截面.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(14):JET的产生截面不考虑核的shadowing效应.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(15):JET的产生截面考虑弹核的shadowing效应修正后.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(16):JET的产生截面考虑靶核的shadowing效应修正后.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(17):
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(18):JET产生的有效截面.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(19):最近一次事件的碰撞参数.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(20):the azimuthal angle phi of the impact parameter
c$$$        vector in the transverse plane of the latest event.最近一次
c$$$        事件中,碰撞参数在横平面的方位角。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(21)-HINT1(25):the four momentum and mass(px, py, pz,
c$$$        E, M)of the first scattered parton in the triggered hard
c$$$        scattering. This is before the final state radiation but
c$$$        after the initial state radiation.(第一次部分子散射后的
c$$$        四动量和质量)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(31)-HINT1(35):the four momentum and mass(px,py,pz,E,M)
c$$$        of the second scattered parton in the triggered hard scattering.
c$$$        this is before the final state radiation but after the initial
c$$$        state radiation.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(41)-HINT1(45):the four momentum and mass(px,py,pz,E,M)
c$$$        of the scattered parton in the latest hard scattering of the
c$$$        latest event.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(46):PT of the first scattered parton in the latest
c$$$        hard scattering of the latest event.(第一次散射后的横动量)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(51)-HINT1(55):the four momentum and mass(px,py,pz,E,M)
c$$$        of the second scattered parton in the latest hard scattering
c$$$        of the latest event.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(56):PT of the second scattered parton in the latest
c$$$        hard scattering of the latest event.(第二次散射后的横动量)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(59):the averaged cross section of the triggered jet
c$$$        production per nucleon-nucleon collision.(JET产生的平均散射
c$$$        截面)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(60):the averaged inclusive cross section of the
c$$$        triggered jet production per nucleon-nucleon.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(61):the triggered jet production cross section without
c$$$        nuclear shadowing effect.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(62):the cross section to account for the projectile
c$$$        shadowing correction term int the triggered jet cross section.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(63):the cross section to account for target shadowing
c$$$        correction term in the triggered jet cross section.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(64):the cross section to account for the term of
c$$$        shadowing correction in the triggered jet cross section.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(65):the inclusive cross section for latest triggered
c$$$        jet production which depends on the transverse coordinates of
c$$$        the colliding nucleons.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(72)-HINT1(75):three parameters for the Wood-Saxon
c$$$        projectile nuclear distribution and the normalization read from
c$$$        a table inside the program.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(76)-HINT1(79):three parameters for the Wood-Saxon
c$$$        target nuclear distribution and the normalization read from
c$$$        a table inside the program.(详细的请看25页的具体参数的介绍)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        HINT1(80)-HINT1(100):the probability of j = 0-20 number of
c$$$        hard scatterings per nucleon-nucleon collisions.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        




c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(1):The mass number of the projectile nucleus.(弹核
c$$$        核子的质量数)(1 for a hadron).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(2):the charge number of the projectile nucleus.(弹
c$$$        核的电荷数)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(3):The mass number of the target nucleus.(靶核
c$$$        核子的质量数)(1 for a hadron).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(4):the charge number of the taget nucleus.(靶核的电荷
c$$$        数)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(5):the flavor code of the projectile hadron(0 for
c$$$        nucleus).弹核的味道代码。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(6):the flavor code of the target hadron(0 for nucleus)
c$$$        靶核的夸克味道代码。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(9):the flavor code of the first scattered parton in
c$$$        the triggered hard scattering.(第一次散射时部分子的味道)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(10):the flavor code of the second scattered parton in
c$$$        the triggered hard scattering.(第二次散射时部分子的味道)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(11):the sequence(次序) number of the projectile nucleon in
c$$$        the latest nucleon-nucleon interaction of the latest event.
c$$$        最近一次核核(弹核)相互作用的次序数。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(12):the sequence(次序) number of the target nucleon in
c$$$        the latest nucleon-nucleon interaction of the latest event.
c$$$        最近一次核核(靶核)相互作用的次序数。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(13):status of the latest soft string excitation.
c$$$        = 1 : double diffractive.
c$$$        = 2 : single diffractive.
c$$$        = 3 : non-single diffractive.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(14):the flavor code of the first scattered parton in
c$$$        the latest hard scattering of the latest event.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IHNT2(15):the flavor code of the second scattered parton in
c$$$        the latest hard scattering of the latest event.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
C
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        分别为弹核和靶核的三坐标YP(3,300), YT(3,300)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/hjcrdn/YP(3,300),YT(3,300)
cc      SAVE /hjcrdn/
clin-7/16/03 NINT is a intrinsic fortran function, rename it to NINTHJ
c        COMMON/HJGLBR/NELT,NINT,NELP,NINP
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        弹核和靶核核子碰撞的状态(弹性和非弹性)。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HJGLBR/NELT,NINTHJ,NELP,NINP
cc      SAVE /HJGLBR/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        EATT:总的能量,用来检验能量守恒。
c$$$        JATT:the total number of jet-pairs in the current event.
c$$$        NATT:当前事件产生的稳定的粒子的总数。
c$$$        NT:靶核受伤的核子数。
c$$$        NP:弹核受伤的核子数。
c$$$        N0:       number of N-N, N-N_wounded, N_wounded-N, and 
c$$$        N01:      N_wounded-N_wounded collisions in the current
c$$$        N10:      events.(这几个变量似乎是受伤的核子)。
c$$$        N11:
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
cc      SAVE /HMAIN1/
clin-4/26/01
c        COMMON/HMAIN2/KATT(130000,4),PATT(130000,4)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KATT两维数组的介绍！！！！！
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KATT(I,1): particle ID of the I_th produced particle. Users
c$$$        have to refer to JETSET7.2 for identifying particles with
c$$$        their ID's.(产生的第I个粒子的ID).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KATT(I,2): This is a code to identify the sources from which
c$$$        the particle comes.
c$$$        = 0 : projectile which has not interacted at all.
c$$$        = 1 : projectile nucleon (or hadron) which only suffers an
c$$$              elastic collision.
c$$$        = 2 : from a diffractive projectile nucleon (or hadron) in a
c$$$              single diffractive interaction.
c$$$        = 3 : from the fragmentation of a projectile string
c$$$              system(including gluon jets).
c$$$        = 10: target nucleon (or hadron) which has not interacted at
c$$$              all.
c$$$        = 11: target nucleon (or hadron) which only suffers an elastic
c$$$              collision.
c$$$        = 12: from a diffractive target nucleon (or hadron) in a single
c$$$              diffractive interaction.
c$$$        = 13: from the fragmentation of a target string system
c$$$              (including gluon jets).
c$$$        = 20: from scattered partons which form string systems
c$$$              themeselves.
c$$$        = 40: from direct production in the hard process (currently,
c$$$              only direct photons are included).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KATT(I,3): (I=1,...,NATT)line number of the parent particle.
c$$$        For finally produced or directly produced (not from the decay of
c$$$        another particle)particles, it is set to 0(The option to keep
c$$$        the information of all particles including the decayed ones is
c$$$        IHPR2(21)=1).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KATT(I,4): (I=1,...,NATT)status number of the particle.
c$$$        = 1 : finally of directly produced particles.
c$$$        = 11: particles which has already decayed.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PATT(I,1-4): (I=1,...,NATT) four-momentum(px, py, pz, E) of
c$$$        the produced particles.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
cc      SAVE /HMAIN2/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,1):(I=1,...,IAP)flavor code of the valence quark in
c$$$        projectile  nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,2):flavor code of diquark in projectile nucleon
c$$$       (anti_quark in projectile meson) I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,3):present flavor code of the projectile nucleon
c$$$        (hadron) I(a necleon of meson can be excired to its vector
c$$$        resonance).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,4):original flavor code of projectile
c$$$        nucleon(hadron)I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,5):collision status of projectile nucleon(hadron) I.
c$$$        = 0 : suffered no collision.
c$$$        = 1 : suffered an elastic collision.
c$$$        = 2 : being the diffractive one in a single-diffractive
c$$$        collision.
c$$$        = 3 : became an excited string after an inelastic collision. 
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,6): the total number of hard scattering associated
c$$$        with projectile nucleon(hadron) I. if NEP(I,6)<0, it can not
c$$$        produce jets any more due to energy conservation.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,10): to indicate whether the valence quarks of
c$$$        diquarks (anti_quarks) in projectile nucleon (hadron) I suffered
c$$$        a hard scattering.
c$$$        = 0 : has not suffered a hard scattering.
c$$$        = 1 : suffered one or more hard scatterings in current binary
c$$$        nucleon-nucleon collision.
c$$$        = -1: suffered one or more hard scatterings in previous binary
c$$$        nucleon-nucleon collisions.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NFP(I,11): total number of interactions projectile nucleon
c$$$        (hadron) I has suffered so far.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,1),PP(I,2),PP(I,3),PP(I,4),PP(I,5):four momentum and
c$$$        the invariant mass (px,py,pz,E,M) of projectile nucleon(hadron)I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,6), PP(I,7): transverse momentum (px,py)of the valence
c$$$        quark in projectile nucleon(hadron)I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,8), PP(I,9):transverse momentum (px,py)of the diquark
c$$$        (anti_quark) in projectile nucleon(hadron)I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,10), PP(I,11), PP(I,12):three momentum (px,py,pz)
c$$$        transferred to the quark or diquark (anti_quark) in projectile
c$$$        nucleon (hadron) I from the last hard scattering.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,14):mass of the quark in projectile nucleon(hadron) I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PP(I,15):mass of the diquark (anti_quark) in projectile
c$$$        nucleon (hadron) I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
cc      SAVE /HSTRNG/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NPJ(I):(I=1,...,IAP)number of partons associated with
c$$$        projectile nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KFPJ(I,J):(I=1,...,IAP,J=1,...,NPJ(I))parton flavor code of
c$$$        the parton J associated with projectile nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PJPX(I,J),PJPY(I,J),PJPZ(I,J),PJPM(I,J):the four momentum
c$$$        and mass (px,py,pz,E,M) of parton J associated with the
c$$$        projectile nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NTJ(I):(I=1,...,IAT)number of partons associated with
c$$$        target nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KFTJ(I,J):(I=1,...,IAT,J=1,...,NTJ(I)):parton flavor code
c$$$        of the parton J associated with target nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PJTX(I,J),PJTY(I,J),PJTZ(I,J),PJTE(I,J),PJTM(I,J):
c$$$        the four momentum and mass (px,py,pz,E,M) of parton J associated
c$$$        with the target nucleon I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &                PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &                PJPM(300,500),NTJ(300),KFTJ(300,500),
     &                PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &                PJTE(300,500),PJTM(300,500)
cc      SAVE /HJJET1/
clin-4/2008
c        COMMON/HJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100),
c     &       K2SG(900,100),PXSG(900,100),PYSG(900,100),
c     &       PZSG(900,100),PESG(900,100),PMSG(900,100)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NSG:the total number of such string systems.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NJSG(I):(I=1,...,NSG)number of partons in the string system I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IASG(I,1),IASG(I,2): to specify which projectile and target
c$$$        nucleons produce string system I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        K1SG(I,J):(J=1,...,NJSG(I))color flow information of parton
c$$$        J in string system I (see JETSET 7.2 for detailed explanation).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        K2SG(I,J):(J=1,...,NJSG(I))flavor code of parton J in
c$$$        string system I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PXSG(I,J),PYSG(I,J),PZSG(I,J),PESG(I,J),PMSG(I,J):four
c$$$        momentum and mass (px,py,pz,E,M) of parton J in string system I
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &       K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &       PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
cc      SAVE /HJJET2/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NDR:total number of directly produced particles.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        IADR(I,1),IADR(I,2):the sequence numbers of projectile and
c$$$        target nucleons which produce particle I during their
c$$$        interaction.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        KFDR(I):the flavor code of directly produced particle I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        PDR(I,1,...,5):four momentum and mass (px,py,pz,E,M)of
c$$$        particle I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        COMMON/HJJET4/NDR,IADR(MAXSTR,2),KFDR(MAXSTR),PDR(MAXSTR,5)
clin-4/2008:
c        common/xydr/rtdr(900,2)
        common/xydr/rtdr(MAXSTR,2)
cc      SAVE /HJJET4/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
C
        COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)   
cc      SAVE /LUJETS/
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
cc      SAVE /LUDAT1/
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        changed name in order to distinguish from /prec2/(进行区分)
c$$$        参看下面下面下面的数据块进行区分。
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
clin-9/29/03 changed name in order to distinguish from /prec2/
        COMMON /ARPRC/ ITYPAR(MAXSTR),
     &       GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &       PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &       XMAR(MAXSTR)
ccbz11/11/98




c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /PARA1/ MUL
c$$$        mul : the multiplicity of one particular event.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /PARA2/ xmp, xmu, alpha, rts_cut2, cutoff2
c$$$        xmp : the mass of the particle.
c$$$        xmu : the screening mass.
c$$$        alpha : the strong interaction coupling constant.
c$$$        rts_cut2 : lower square root of s cut for parton collisions.
c$$$        cutoff2 : the total parton scattering cross section divided by pion.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /PARA3/ n_sta_evt, n_event, n_run, i_event, i_run
c$$$        n_sta_evt : the starting event number in the input file.
c$$$        n_event : the total number of events.
c$$$        n_run : the total number of runs per event.
c$$$        i_event : the current event number that goes from 1 to n_event.
c$$$        i_run : the current run number.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /PARA4/ ift_flag, ire_flag, igen_flag, ibst_flag
c$$$        ift_flag : switch for the generation of particle formation time.
c$$$          = 0 : generation formation time.
c$$$          = 1 : read in formation time from zpc.inp.
c$$$        ire_flag : switch to read initial conditions from zpc.inp.
c$$$          = 0 : read the initial conditions.
c$$$          = 1 : do not read the initial conditions.
c$$$        igen_flag : switch to generate initial conditions.
c$$$          = 0 : do not generate initial conditions.
c$$$          = 1 : generate initial conditions.
c$$$        ibst_flag : choice of global frame.
c$$$          = 0: the global frame (the ordering frame) is the colliser center
c$$$               of mass frame for Au on Au collisons.
c$$$          = 1 : means the global frame is taken to be the target frame.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /PARA5/ i_config, iord_sch
c$$$        i_config : choice of geometric configuration and space cell
c$$$                   dividion optimization.
c$$$          = 1 : the system is undergoing 3-d expansion.
c$$$          = 2 : the system is undergoing 1-d expansion with space cell
c$$$                division.
c$$$          = 3 : the system is undergoing 1-d expansion without space cell
c$$$                division.
c$$$          = 4 : the system is confined in a box with space cell division.
c$$$          = 5 : the system is confined in a box without space cell division.
c$$$        iord_sch : choice of collision scheme: iord_sch = 10*i1 + i2
c$$$        i1 : choice of collision prescription.
c$$$          = 0 : the collison frame is the two parton center frame. the
c$$$                scattering space point is the center point of the 2 parton
c$$$                positions at closest approach in the collision frame.
c$$$          = 1 : the collision frame is the two parton center of mass frame.
c$$$                the scattering point is the position of the parton at the closest
c$$$                approach in the collision frame. this gives two generally
c$$$                different collision times for the two colliding partons.
c$$$        = 2 : the collisons frame is the same as the ordering global
c$$$              frame.
c$$$        i2 : choice of ordering time.
c$$$          = 0 : the earlier of the two collision times in the global frame
c$$$                is the ordering time for parton collisions.
c$$$          = 1 : the average of the two collision times is taken as the
c$$$                ordering time.
c$$$          = 2 : the later collision time is the ordering time.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /para6/ cent_rap
c$$$        cent_rap : the central rapidity of the rapidity plateau.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /para7/ i_oscar
c$$$        i_oscar : option to write oscat standard output.
c$$$          = 0 : write oscar standard output.
c$$$          = 1 : do not write oscar standard output.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /prec1/
c$$$        purpose:stores the initial conditions for the event.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /prec2/
c$$$        purpose:main particle record that has particle types, the
c$$$        positions at formation times, the corresponding formation times
c$$$        and the four momenta of the particle, plus the masses of the
c$$$        particles.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /prec3/
c$$$        purpose:the same as /prec2/ expect the two colliding particles
c$$$        do not update their records until the next collision.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /prec4/
c$$$        purpose:stores vecolities of particles.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /prec5/
c$$$        purpose:stores formation time pseudo-rapidities, rapidities and
c$$$        proper times of particles.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /prec6/
c$$$        purpose:the same as /prec5/ except the two colliding particles do
c$$$        not update their records until the next collision
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<






      
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        six interaction list common block i_list* give us particle
c$$$        run time information and cell information
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /I_list1/ ISCAT,JSCAT,NEXT(NMAXGL), LAST(NMAXGL),
c$$$        ICTYPES, ICSTA(NMAXGL), NIC(NMAXGL), ICELSTA(NMAXGL)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        ISCAT:particle index. if the operation involves only one
c$$$        particle, then iscat is the particle index, if it involves 2
c$$$        particles, iscat is the larger particle index.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        JSCAT:particle index. if the operation involves only one
c$$$        particle, then iscat is 0, if it involves 2
c$$$        particles, iscat is the smaller particle index.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NEXT(I):the next operation partner of particle I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        LAST(I):the last operation partner of particle I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        ICTYPE:the operation type.
c$$$        = 0: a collision between particles.
c$$$        = 1: the formation of a particle.
c$$$        = 2: both a collision between particles and the formation of a
c$$$             particle.
c$$$        = 3: a collision with wall.
c$$$        = 4: both a collision between particles and a wall collision.
c$$$        = 5: both a wall collision and a formation.
c$$$        = 6: a formation, collision between particles and a wall collision
c$$$             at the same time.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        ICSTA(I):the operation type for particle I.
c$$$        = 0: an ordinary collision.
c$$$        = 101:a collision with the wall with larger x.
c$$$        = 102:a collision with the wall with smaller x.
c$$$        = 103:a collision with the wall with larger y.
c$$$        = 104:a collision with the wall with smaller y.
c$$$        = 105:a collision with the wall with larger z.
c$$$        = 106:a collision with the wall with smaller z.
c$$$        = 111:a collision with another particle and the wall with lartger x
c$$$        = 112:a collision with another particle and the wall with smaller x
c$$$        = 113:a collision with another particle and the wall with lartger y
c$$$        = 114:a collision with another particle and the wall with lartger y
c$$$        = 115:a collision with another particle and the wall with lartger z
c$$$        = 116:a collision with another particle and the wall with lartger z
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        NIC(I):the next particle index in the same cell as I.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        ICELSTA(I):the encoded information of the cell number
c$$$        particle I is in. If particle is in cell (i1, i2, i3), then it
c$$$        equals i1 * 10000 + i2 * 100 + i3 for particles inside the cube.
c$$$        when a particle is outside the 10 by 10 by 10 box. its value is 111111.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /I_LIST2/ ICELL, ICEL(10,10,10)
c$$$        ICELL:pointer to one particle outside the 10 by 10 by 10
c$$$        cell volume.
c$$$        ICEL(I,J,K):one particle in the cell(I,J,K).
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /I_LIST3/ SIZE1, SIZE2, SIZE3, V1, V2, V3, SIZE
c$$$        SIZE1: cell size in the x direction. it must be larger than the
c$$$        inverse square root of cutoff2.
c$$$        SIZE2: cell size in the y direction.
c$$$        SIZE3: cell size in the z direction.
c$$$        V1: cell expanding velocity in x direction.
c$$$        V2: cell expanding velocity in y direction.
c$$$        V3: cell expanding velocity in z direction.
c$$$        SIZE: the time that cells begin to expand
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /I_LIST4/ IFMPT, ICHKPT, INDX(NMAXGL)
c$$$        IFMPT: the partilce index of the parton that is to be formed next.
c$$$        ICHKPT: the last formed particle index.
c$$$        INDX(I): the index of particle I of /INI_REC/ in an incresing
c$$$        formation time order.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /I_LIST5/ CT(NMAXGL), OT(NMAXGL), TLARGE
c$$$        CT(I): the collision time for particle I.
c$$$        OT(I): The ordering time for particle I.
c$$$        TLARGE: A large number for time cutoff.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /I_LIST6/ I_OPERATION, I_COLLISION, T
c$$$        N_OPERATION: the current nember of operations.
c$$$        I_COLLISION: the current number of collision between particles
c$$$        that have been performed.
c$$$        T: the current operation time.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        COMMON /RANDOM/ NUMBER
c$$$        NUMBER: the number of random numbers that have been used.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        COMMON /RANDOM2/ IFF
c$$$        IFF: choice  of attractive or repulsive force. It alternates
c$$$        between 1 and -1.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        Currentlt, there are four common blocks ana* for analysis.
c$$$        COMMON /ANA1/ TS(12)
c$$$        TS: stores the time points that we want to sample data for analysis.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        COMMON /ANA2/
c$$$        purpose: to record the dE/dy, dE/dN, and dN/dy time evolution.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        COMMON /ANA3/ EM(4,4,12)
c$$$        EM: energy momentum tensor for 12 different time points.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        COMMON /ANA4/ FDETDY(24), FDNDY(24), FDNDPT(12)
c$$$        FDETDY: final dE/dy distribution.
c$$$        FDNDY: final dN/dy distribution.
c$$$        FDNDPT: final dN/dp distribution.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
c.................... zpc.f
c	PROGRAM ZPC
      SUBROUTINE ZPCMN
c       Version: 1.0.1
c       Author: Bin Zhang 
c       (suggestions, problems -> bzhang@nt1.phys.columbia.edu)
        implicit double precision (a-h, o-z)
clin-4/20/01        PARAMETER (NMAXGL = 16000)
        parameter (MAXPTN=400001)
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        common /PARA3/ n_sta_evt, n_event, n_run, i_event, i_run
c$$$        n_sta_evt : the starting event number in the input file.
c$$$        n_event : the total number of events.
c$$$        n_run : the total number of runs per event.
c$$$        i_event : the current event number that goes from 1 to n_event.
c$$$        i_run : the current run number.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
cc      SAVE /para3/
        SAVE   
c
c       loop over events
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        遍历所有的事件.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do 1000 i = 1, nevnt
           ievt = i
c       generation of the initial condition for one event
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        调用了此函数,但是里面什么屁功能都没有.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$           call inievt              !!!!!!初始事件。
c      loop over many runs of the same event
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
           write(*,*) "nsbrun = ", nsbrun
c$$$       这里run的总数好像永远为1，也就是说下面的循环只进行一次。           
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
           do 2000 j = 1, nsbrun
              isbrun = j
c       initialization for one run of an event
              call inirun               
clin-4/2008 not used:
c             CALL HJAN1A
 3000         continue
c       do one collision
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        进行一次碰撞.
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              call zpcrun(*4000)
              call zpca1
              goto 3000
 4000         continue
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        子程序zpca2的功能和zpca1的功能很相近,这里对比下.


c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$         WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$          WW     WW WW     WW  RR  RR   II     TT     EE
c$$$           WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$            WW WW     WW WW    RR RR    II     TT     EE
c$$$             WW        WW      RR  RR   II     TT     EEEEEE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<
            
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$         WW       WW        WW RRRRRRR  II TTTTTTTTTT EEEEEE
c$$$          WW     WW WW     WW  RR  RR   II     TT     EE
c$$$           WW   WW   WW   WW   RRRRR    II     TT     EEEEEE
c$$$            WW WW     WW WW    RR RR    II     TT     EE
c$$$             WW        WW      RR  RR   II     TT     EEEEEE
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<


c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



      
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c$$$        
c$$$      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
