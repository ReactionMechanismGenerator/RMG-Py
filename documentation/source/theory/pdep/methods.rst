**********************************************************
Methods for Determining Phenomenological Rate Coefficients
**********************************************************

Solving the energy-grained master equation is often prohibitively expensive for 
chemical reaction mechanisms of even modest size. Instead, we seek to reduce
the master equation matrix down to a set of phenomenological rate coefficients
:math:`k(T, P)`. In particular, we need to replace the isomer population 
distributions :math:`p_i(E, t)` with the corresponding time-dependent total 
isomer populations :math:`x_i(t)`.

Three methods of varying rigor, computational cost, and robustness will be
discussed in the upcoming sections. The modified strong collision (MSC) method
is the fastest and most robust, but utilizes the least realistic approximations.
The reservoir state (RS) method uses better approximations, which leads to
increased accuracy, but requires more computational effort. Finally, the
chemically-significant eigenvalues (CSE) method is the most theoretically sound,
but is very computationally expensive and not very robust. Your choice of
method will depend on the particular balance between expense, robustness, and
rigor that is required for your intended application.

A Common Formalism
==================

All of the methods discussed here can be expressed in terms of a common 
formalism. Each method seeks to express the population distribution vector 
:math:`p_i(E, t)` for each unimolecular isomer :math:`i` as a linear 
combination of the total populations :math:`x_j(t)` and 
:math:`y_{m\mathrm{A}}(t) y_{m\mathrm{B}}` of unimolecular isomers 
:math:`\ce{A}_j` and reactant channels :math:`\ce{A}_m + \ce{B}_m`:

.. math:: p_i(E, t) = \sum_{j=1}^{N_\mathrm{isom}} x_j(t) u_{ij}(E) + \sum_{m=1}^{N_\mathrm{reac}} y_{m\mathrm{A}}(t) y_{m\mathrm{B}} v_{im}(E)

The function :math:`u_{ij}(E)` represents the portion of the population 
distribution of unimolecular isomer :math:`i` at energy :math:`E` that tracks 
the population of isomer :math:`j`. In the modified strong collision and 
reservoir state methods, this is because the energy levels of isomer :math:`i`
are in pseudo-steady-state relationships with isomer :math:`j`. The 
interpretation is a bit different for the chemically-significant eigenvalues 
method, but the form of the equations is the same. Similarly, the function
:math:`v_{im}(E)` represents the population distribution of unimolecular isomer
:math:`i` at energy :math:`E` that tracks the population of reactant channel 
:math:`m`. Both functions :math:`u_{ij}(E)` and :math:`v_{im}(E)` are functions
of energy only, and not of time.

After discretizing the energy domain, the above becomes

.. math:: 

    \renewcommand{\vector}[1]{\boldsymbol{\mathbf{#1}}}
    \renewcommand{\matrix}[1]{\boldsymbol{\mathbf{#1}}}
    \vector{p}_i(t) = \sum_{j=1}^{N_\mathrm{isom}} x_j(t) \vector{u}_{ij} + \sum_{m=1}^{N_\mathrm{reac}} y_{m\mathrm{A}}(t) y_{m\mathrm{B}} \vector{v}_{im}

The phenomenological rate coefficients can be constructed from the 
energy-grained master equation matrix and the vectors :math:`\vector{u}_{ij}` 
and :math:`\vector{v}_{im}`:

.. math:: 

    k_{ij}(T,P) &= \sum_{s=1}^{N_\mathrm{grains}} \left( \matrix{M}_i \vector{u}_{ij} \right)_s + 
                   \sum_{\ell \ne i}^{N_\mathrm{isom}} \sum_{s=1}^{N_\mathrm{grains}} \left( \matrix{K}_{i \ell} \vector{u}_{\ell j} \right)_s \\
    k_{im}(T,P) &= \sum_{s=1}^{N_\mathrm{grains}} \left( \matrix{M}_i \vector{v}_{im} \right)_s + 
                   \sum_{\ell \ne i}^{N_\mathrm{isom}} \sum_{s=1}^{N_\mathrm{grains}} \left( \matrix{K}_{i \ell} \vector{v}_{\ell m} \right)_s + \sum_{s=1}^{N_\mathrm{grains}} \left( \matrix{F}_{im} \vector{b}_m \right)_s \\
    k_{nj}(T,P) &= \sum_{\ell=1}^{N_\mathrm{isom}} \vector{g}_{n \ell} \cdot \vector{u}_{\ell j} \\
    k_{nm}(T,P) &= \sum_{\ell=1}^{N_\mathrm{isom}} \vector{g}_{n \ell} \cdot \vector{v}_{\ell m} 

Above, the indices :math:`i` and :math:`j` represent unimolecular isomers of 
the initial adduct, :math:`m` represents bimolecular reactants, :math:`n` 
represents bimolecular reactants and products, and :math:`s` represents an 
energy grain. Thus, the rate coefficients above are for isomerization, 
association, dissociation, and bimolecular reactions, respectively.

The output from each of the three methods is a set of phenomenological rate
coefficients :math:`k(T, P)` and the vectors :math:`\vector{u}_{ij}` and
:math:`\vector{v}_{im}` which can be used to construct the approximate
population distribution predicted by that method.

The Modified Strong Collision Method
====================================

The Reservoir State Method
==========================

The Chemically-Signficant Eigenvalues Method
============================================
