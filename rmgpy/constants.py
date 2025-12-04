#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

r"""
The :mod:`rmgpy.constants` module contains module-level variables defining
relevant physical constants relevant in chemistry applications. The recommended
method of importing this module is ::

    import rmgpy.constants as constants

so as to not place the constants in the importing module's global namespace.

The constants defined in this module are listed in the table below:

.. table:: Physical constants defined in the :mod:`rmgpy.constants` module

Values updated November 2025 using exact values (when available) sourced from https://arxiv.org/pdf/2409.03787
    
    ======================= =================== =========================================================== ============================================================
    Symbol                  Constant            Value                                                       Description
    ======================= =================== =========================================================== ============================================================
    :math:`E_\mathrm{h}`    :data:`E_h`         :math:`4.3597447222060 \times 10^{-18} \ \mathrm{J}`        Hartree energy
    :math:`F`               :data:`F`           :math:`96 48533212 \ldots \ \mathrm{C/mol}`                 Faraday constant (exact)
    :math:`G`               :data:`G`           :math:`6.67430 \times 10^{-11} \ \mathrm{m^3/kg \cdot s^2}` Newtonian gravitational constant
    :math:`N_\mathrm{A}`    :data:`Na`          :math:`6.02214076 \times 10^{23} \ \mathrm{mol^{-1}}`       Avogadro constant (exact)
    :math:`R`               :data:`R`           :math:`8.314462618 \ldots \ \mathrm{J/mol \cdot K}`         gas law constant (exact) 
    :math:`a_0`             :data:`a0`          :math:`5.29177210544 \ldots \times 10^{-11} \ \mathrm{m}`   Bohr radius
    :math:`c`               :data:`c`           :math:`299792458 \ \mathrm{m/s}`                            speed of light in a vacuum (exact)
    :math:`e`               :data:`e`           :math:`1.602176634 \times 10^{-19} \ \mathrm{C}`            elementary charge (exact)
    :math:`g`               :data:`g`           :math:`9.80665 \ \mathrm{m/s^2}`                            standard acceleration due to gravity (exact)
    :math:`h`               :data:`h`           :math:`6.62607015 \times 10^{-34} \ \mathrm{J \cdot s}`     Planck constant (exact)
    :math:`\hbar`           :data:`hbar`        :math:`1.054571817 \times 10^{-34} \ \mathrm{J \cdot s}`    reduced Planck constant (exact)
    :math:`k_\mathrm{B}`    :data:`kB`          :math:`1.380649 \times 10^{-23} \ \mathrm{J/K}`             Boltzmann constant (exact)
    :math:`m_\mathrm{e}`    :data:`m_e`         :math:`9.1093837139 \times 10^{-31} \ \mathrm{kg}`          electron rest mass
    :math:`m_\mathrm{n}`    :data:`m_n`         :math:`1.67492750056 \times 10^{-27} \ \mathrm{kg}`         neutron rest mass
    :math:`m_\mathrm{p}`    :data:`m_p`         :math:`1.67262192595 \times 10^{-27} \ \mathrm{kg}`         proton rest mass
    :math:`M_\mathrm{u}`    :data:`M_u`         :math:`1.0000000010531 \times 10^{-3} \ \mathrm{kg/mol}`    molar mass constant
    :math:`m_\mathrm{u}`    :data:`amu`         :math:`1.66053906892 \times 10^{-27} \ \mathrm{kg}`         atomic mass unit
    :math:`\pi`             :data:`pi`          :math:`3.14159 \ldots`      
    ======================= =================== =========================================================== ============================================================

"""

import math

################################################################################

#: :math:`\pi = 3.14159 \ldots`      
pi = float(math.pi)

#: The Avogadro constant :math:`N_\mathrm{A}` in :math:`\mathrm{mol^{-1}}`       
Na = 6.02214076e23

#: The Avogadro constant :math:`N_\mathrm{A}` in :math:`\mathrm{mol^{-1}}`       
Mu = 1.0000000010531e-3

#: The Boltzmann constant :math:`k_\mathrm{B}` in :math:`\mathrm{J/K}`
kB = 1.380649e-23

#: The speed of light in a vacuum :math:`c` in :math:`\mathrm{m/s}`
c = 299792458

#: The elementary charge :math:`e` in :math:`\mathrm{C}`
e = 1.602176634e-19

#: The Planck constant :math:`h` in :math:`\mathrm{J \cdot s}`
h = 6.62607015e-34

#: The reduced Planck constant :math:`\hbar` in :math:`\mathrm{J \cdot s}`
hbar = h / (2 * pi)

#: The gas law constant :math:`R` in :math:`\mathrm{J/mol \cdot K}`
R = kB * Na 

#: Vacuum permittivity 
epsilon_0 = 8.8541878188e-12

#: The Bohr radius :math:`a_0` in :math:`\mathrm{m}`
a0 = 5.29177210544e-11

#: The Hartree energy :math:`E_\mathrm{h}` in :math:`\mathrm{J}`
E_h = e**2 / (4 * pi * epsilon_0 * a0)

#: Faradays Constant F in C/mol
F = Na * e 

#: The atomic mass unit in :math:`\mathrm{kg}`
amu = Mu / Na

#: The mass of an electron :math:`m_\mathrm{e}` in :math:`\mathrm{kg}`
m_e = 9.1093837139e-31

#: The mass of a neutron :math:`m_\mathrm{n}` in :math:`\mathrm{kg}`
m_n = 1.67492750056e-27

#: The mass of a proton :math:`m_\mathrm{p}` in :math:`\mathrm{kg}`
m_p = 1.67262192595e-27


################################################################################

# Cython does not automatically place module-level variables into the module
# symbol table when in compiled mode, so we must do this manually so that we
# can use the constants from both Cython and regular Python code
globals().update({
    'E_h': E_h,
    'Na': Na,
    'R': R,
    'a0': a0,
    'amu': amu,
    'c': c,
    'e': e,
    'h': h,
    'hbar': hbar,
    'kB': kB,
    'm_e': m_e,
    'm_n': m_n,
    'm_p': m_p,
    'pi': pi,
    'F': F,
    'epsilon_0': epsilon_0,
})
