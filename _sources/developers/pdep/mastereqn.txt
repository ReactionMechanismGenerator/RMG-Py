*******************
The master equation
*******************

Isomers, reactants, and products
================================

Throughout this document we will utilize the following terminology:

* An **isomer** is a unimolecular configuration on the potential energy surface.

* A **reactant channel** is a bimolecular configuration that associates to form
  an isomer. Dissociation from the isomer back to reactants is allowed.

* A **product channel** is a bimolecular configuration that is formed by
  dissociation of an isomer. Reassociation of products to the isomer is *not*
  allowed.

The isomers are the configurations for which we must model the energy states.
We designate :math:`p_i(E, J, t)` as the population of isomer :math:`i` having
total energy :math:`E` and total angular momentum quantum number :math:`J` at
time :math:`t`. At long times, statistical mechanics requires that the
population of each isomer approach a Boltzmann distribution :math:`b_i(E, J)`:

.. math:: \lim_{t \rightarrow \infty} p_i(E, J, t) \propto b_i(E, J)

We can simplify by eliminating the angular momentum quantum number to get

.. math:: p_i(E, t) = \sum_J p_i(E, J, t)

Let us also denote the (time-dependent) total population of isomer :math:`i` 
by :math:`x_i(t)`:

.. math:: x_i(t) \equiv \sum_J \int_0^\infty p_i(E, J, t) \ dE

The two molecules of a reactant or product channel are free to move apart from
one another and interact independently with other molecules in the system. 
Accordingly, we treat these channels as fully thermalized, leaving as the only
variable the total concentrations :math:`y_{n\mathrm{A}}(t)` and 
:math:`y_{n\mathrm{B}}(t)` of the molecules :math:`\mathrm{A}_n` and :math:`\mathrm{B}_n`
of reactant channel :math:`n`. (Since the product channels act as infinite 
sinks, their populations do not need to be considered explicitly.)

Finally, we will use :math:`N_\mathrm{isom}`, :math:`N_\mathrm{reac}`, and
:math:`N_\mathrm{prod}` as the numbers of isomers, reactant channels, and
product channels, respectively, in the system.

The master equation
===================

The governing equation for the population distributions :math:`p_i(E, J, t)`
of each isomer :math:`i` and the reactant concentrations 
:math:`y_{n\mathrm{A}}(t)` and :math:`y_{n\mathrm{B}}(t)` combines the
collision and reaction models to give a linear integro-differential equation:

.. math::

    \frac{d}{dt} p_i(E, J, t) &= \omega_i(T, P) \sum_{J'} \int_0^\infty P_i(E, J, E', J') p_i(E', J', t) \ dE' - \omega_i(T, P) p_i(E, J, t) \\
                           &  \mbox{} + \sum_{j \ne i}^{N_\mathrm{isom}} k_{ij}(E, J) p_j(E, J, t) - \sum_{j \ne i}^{N_\mathrm{isom}} k_{ji}(E, J) p_i(E, J, t) \\
                           &  \mbox{} + \sum_{n=1}^{N_\mathrm{reac}} y_{n\mathrm{A}}(t) y_{n\mathrm{B}}(t) f_{in}(E, J) b_n(E, J, t) - \sum_{n=1}^{N_\mathrm{reac} + N_\mathrm{prod}} g_{ni}(E, J) p_i(E, J, t) \\
    
    \frac{d}{dt} y_{n\mathrm{A}}(t) = \frac{d}{dt} y_{n\mathrm{B}}(t) &= \sum_{i=1}^{N_\mathrm{isom}} \int_0^\infty g_{ni}(E, J) p_i(E, J, t) \ dE \\
                                                                      &  \mbox{} - \sum_{i=1}^{N_\mathrm{isom}} y_{n\mathrm{A}}(t) y_{n\mathrm{B}}(t) \int_0^\infty f_{in}(E, J) b_n(E, J, t) \  dE

A summary of the variables is given below:

.. table::

    ======================================= ========================================
    Variable                                Meaning
    ======================================= ========================================
    :math:`p_i(E, J, t)`                    Population distribution of isomer :math:`i`
    :math:`y_{n\mathrm{A}}(t)`              Total population of species :math:`\mathrm{A}_n` in reactant channel :math:`n`
    :math:`\omega_i(T, P)`                  Collision frequency of isomer :math:`i`
    :math:`P_i(E, J, E', J')`               Collisional transfer probability from :math:`(E', J')` to :math:`(E, J)` for isomer :math:`i`
    :math:`k_{ij}(E, J)`                    Microcanonical rate coefficient for isomerization from isomer :math:`j` to isomer :math:`i`
    :math:`f_{im}(E, J)`                    Microcanonical rate coefficient for association from reactant channel :math:`m` to isomer :math:`i`
    :math:`g_{nj}(E, J)`                    Microcanonical rate coefficient for dissociation from isomer :math:`j` to reactant or product channel :math:`n`
    :math:`b_n(E, J, t)`                    Boltzmann distribution for reactant channel :math:`n`
    :math:`N_\mathrm{isom}`                 Total number of isomers
    :math:`N_\mathrm{reac}`                 Total number of reactant channels
    :math:`N_\mathrm{prod}`                 Total number of product channels
    ======================================= ========================================

The above is called the two-dimensional master equation because it contains two
dimensions: total energy :math:`E` and total angular momentum quantum number 
:math:`J`. In the first equation (for isomers), the first pair of terms 
correspond to collision, the second pair to isomerization, and the final pair
to association/dissociation. Similarly, in the second equation above (for 
reactant channels), the pair of terms refer to dissociation/association.

We can also simplify the above to the one-dimensional form, which 
only has :math:`E` as a dimension:

.. math::

    \frac{d}{dt} p_i(E, t) &= \omega_i(T, P) \int_0^\infty P_i(E, E') p_i(E', t) \ dE' - \omega_i(T, P) p_i(E, t) \\
                           &  \mbox{} + \sum_{j \ne i}^{N_\mathrm{isom}} k_{ij}(E) p_j(E, t) - \sum_{j \ne i}^{N_\mathrm{isom}} k_{ji}(E) p_i(E, t) \\
                           &  \mbox{} + \sum_{n=1}^{N_\mathrm{reac}} y_{n\mathrm{A}}(t) y_{n\mathrm{B}}(t) f_{in}(E) b_n(E, t) - \sum_{n=1}^{N_\mathrm{reac} + N_\mathrm{prod}} g_{ni}(E) p_i(E, t) \\
    
    \frac{d}{dt} y_{n\mathrm{A}}(t) = \frac{d}{dt} y_{n\mathrm{B}}(t) &= \sum_{i=1}^{N_\mathrm{isom}} \int_0^\infty g_{ni}(E) p_i(E, t) \ dE \\
                                                                      &  \mbox{} - \sum_{i=1}^{N_\mathrm{isom}} y_{n\mathrm{A}}(t) y_{n\mathrm{B}}(t) \int_0^\infty f_{in}(E) b_n(E, t) \ dE


The equations as given are nonlinear, both due to the presence of the 
bimolecular reactants and because both :math:`\omega_i` and 
:math:`P_i(E, E')` depend on the composition, which is changing with time.
The rate coefficients can be derived from considering the pseudo-first-order 
situation where :math:`y_{n\mathrm{A}}(t) \ll y_{n\mathrm{B}}(t)`, and all 
:math:`y(t)` are negligible compared to the bath gas :math:`\mathrm{M}`. From these 
assumptions the changes in :math:`\omega_i`, :math:`P_i(E, E')`, and all 
:math:`y_{n\mathrm{B}}` can be neglected, which yields a linear equation system.

The energy-grained master equation
==================================

Except for the simplest of unimolecular reaction networks, both the
one-dimensional and two-dimensional master equation must be solved numerically.
To do this we must discretize and truncate the energy domain into a finite
number of discrete bins called *grains*. This converts the linear 
integro-differential equation into a system of first-order ordinary
differential equations:

.. math::

    \frac{d}{dt} \begin{bmatrix}
    \mathbf{p}_1 \\
    \mathbf{p}_2 \\
    \vdots \\
    y_{1\mathrm{A}} \\
    y_{2\mathrm{A}} \\
    \vdots
    \end{bmatrix} = \begin{bmatrix}
    \mathbf{M}_1 & \mathbf{K}_{12} & \ldots & \mathbf{F}_{11} \mathbf{b}_1 y_{1\mathrm{B}} & \mathbf{F}_{12} \mathbf{b}_2 y_{2\mathrm{B}} & \ldots \\
    \mathbf{K}_{21} & \mathbf{M}_2 & \ldots & \mathbf{F}_{21} \mathbf{b}_1 y_{1\mathrm{B}} & \mathbf{F}_{22} \mathbf{b}_2 y_{2\mathrm{B}} & \ldots \\
    \vdots & \vdots & \ddots & \vdots & \vdots & \ddots \\
    (\mathbf{g}_{11})^T & (\mathbf{g}_{12})^T & \ldots & h_1 & 0 & \ldots \\
    (\mathbf{g}_{21})^T & (\mathbf{g}_{22})^T & \ldots & 0 & h_2 & \ldots \\
    \vdots & \vdots & \ddots & \vdots & \vdots & \ddots 
    \end{bmatrix} \begin{bmatrix}
    \mathbf{p}_1 \\
    \mathbf{p}_2 \\
    \vdots \\
    y_{1\mathrm{A}} \\
    y_{2\mathrm{A}} \\
    \vdots
    \end{bmatrix}

The diagonal matrices :math:`\mathrm{K}_{ij}` and :math:`\mathrm{F}_{in}` and 
the vector :math:`\mathrm{g}_{ni}` contain the microcanonical rate coefficients 
for isomerization, association, and dissociation, respectively:

.. math::
    
    (\mathbf{K}_{ij})_{rs} &= \begin{cases}
    \frac{1}{\Delta E_r} \int_{E_r - \Delta E_r/2}^{E_r + \Delta E_r/2} k_{ij}(E) \, dE & r = s \\
    0 & r \ne s
    \end{cases} \\
    (\mathbf{F}_{in})_{rs} &= \begin{cases}
    \frac{1}{\Delta E_r} \int_{E_r - \Delta E_r/2}^{E_r + \Delta E_r/2} f_{in}(E) \, dE & r = s \\
    0 & r \ne s
    \end{cases} \\
    (\mathbf{g}_{ni})_r &= \frac{1}{\Delta E_r} \int_{E_r - \Delta E_r/2}^{E_r + \Delta E_r/2} g_{ni}(E) \, dE
    
The matrices :math:`\mathbf{M}_i` represent the collisional transfer 
probabilities minus the rates of reactive loss to other isomers and to 
reactants and products:

.. math::
    
    (\mathbf{M}_i)_{rs} = \begin{cases}
    \omega_i \left[ P_i(E_r, E_r) - 1 \right] - \sum_{j \ne i}^{N_\mathrm{isom}} k_{ij}(E_r) - \sum_{n=1}^{N_\mathrm{reac} + N_\mathrm{prod}} g_{ni}(E_r) & r = s \\
    \omega_i P_i(E_r, E_s) & r \ne s
    \end{cases}

The scalars :math:`h_n` are simply the total rate coefficient for loss of 
reactant channel :math:`n` due to chemical reactions:

.. math:: h_n = - \sum_{i=1}^{N_\mathrm{isom}} \sum_{r=1}^{N_\mathrm{grains}} y_{n\mathrm{B}} f_{in}(E_r) b_n(E_r)

The interested reader is referred to any of a variety of other sources for
alternative presentations, of which an illustrative sampling is given here
[Gilbert1990]_ [Baer1996]_ [Holbrook1996]_ [Forst2003]_ [Pilling2003]_.

.. [Gilbert1990] R. G. Gilbert and S. C. Smith. *Theory of Unimolecular and 
   Recombination Reactions*. Blackwell Sci. (1990).

.. [Baer1996] T. Baer and W. L. Hase. *Unimolecular Reaction Dynamics*.
   Oxford University Press (1996).

.. [Holbrook1996] K. A. Holbrook, M. J. Pilling, and S. H. Robertson.
   *Unimolecular Reactions*. Second Edition. John Wiley and Sons (1996).

.. [Forst2003] W. Forst. *Unimolecular Reactions: A Concise Introduction*.
   Cambridge University Press (2003).

.. [Pilling2003] M. J. Pilling and S. H. Robertson. *Annu. Rev. Phys. Chem.* 
   **54**, p. 245-275 (2003).
   `doi:10.1146/annurev.physchem.54.011002.103822 <http://dx.doi.org/10.1146/annurev.physchem.54.011002.103822>`_

Generating the master equation matrix
=====================================

.. autofunction:: rmgpy.pdep.me.generateFullMEMatrix
