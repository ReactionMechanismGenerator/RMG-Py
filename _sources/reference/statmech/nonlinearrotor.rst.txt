*****************************
rmgpy.statmech.NonlinearRotor
*****************************

.. autoclass:: rmgpy.statmech.NonlinearRotor

    A nonlinear rigid rotor is the generalization of the linear rotor to a 
    nonlinear polyatomic system. Such a system is characterized by three moments
    of inertia :math:`I_\mathrm{A}`, :math:`I_\mathrm{B}`, and :math:`I_\mathrm{C}` instead of just one.
    The solution to the Schrodinger equation for the quantum nonlinear rotor is not 
    well defined, so we will simply show the classical result instead:

    .. math:: Q_\mathrm{rot}^\mathrm{cl}(T) = \frac{\pi^{1/2}}{\sigma} \left( \frac{8 k_\mathrm{B} T}{h^2} \right)^{3/2} \sqrt{I_\mathrm{A} I_\mathrm{B} I_\mathrm{C}}
