*********************
rmgpy.kinetics.Wigner
*********************

.. autoclass:: rmgpy.kinetics.Wigner

    An early formulation for incorporating the effect of tunneling is that of 
    Wigner [1932Wigner]_:

    .. math:: \kappa(T) = 1 + \frac{1}{24} \left( \frac{h \left| \nu_\mathrm{TS} \right|}{ k_\mathrm{B} T} \right)^2

    where :math:`h` is the Planck constant, :math:`\nu_\mathrm{TS}` is the
    negative frequency, :math:`k_\mathrm{B}` is the Boltzmann constant, and
    :math:`T` is the absolute temperature. 

    The Wigner formula represents the first correction term in a perturbative
    expansion for a parabolic barrier [1959Bell]_, and is therefore only accurate
    in the limit of a small tunneling correction. There are many cases for which
    the tunneling correction is very large; for these cases the Wigner model is
    inappropriate.
    
.. [1932Wigner] E.\ Wigner. *Phys. Rev.* **40**, p. 749-759 (1932). `doi:10.1103/PhysRev.40.749 <http://dx.doi.org/10.1103/PhysRev.40.749>`_

.. [1959Bell] R.\  P.\  Bell. *Trans. Faraday Soc.* **55**, p. 1-4 (1959). `doi:10.1039/TF9595500001 <http://dx.doi.org/10.1039/TF9595500001>`_
