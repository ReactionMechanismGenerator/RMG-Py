*******************
The Master Equation
*******************

A full treatment of the energy states of each molecule is unfeasible for
molecules larger than diatomics, as there are simply too many states. To 
simplify things we apply the RRKM approximation, which leaves the state of a 
molecule as a function of two quantities: the total energy :math:`E` and total  
angular momentum quantum number :math:`J`. Frequently we will find that even 
this is too difficult, and will only keep the total energy :math:`E` as an 
independent variable.

Isomers, Reactants, and Products
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
variable the total concentrations :math:`y_{n\ce{A}}(t)` and 
:math:`y_{n\ce{B}}(t)` of the molecules :math:`\ce{A}_n` and :math:`\ce{B}_n`
of reactant channel :math:`n`. (Since the product channels act as infinite 
sinks, their populations do not need to be considered explicitly.)

Finally, we will use :math:`N_\mathrm{isom}`, :math:`N_\mathrm{reac}`, and
:math:`N_\mathrm{prod}` as the numbers of isomers, reactant channels, and
product channels, respectively, in the system.

Collision Models
================

Bimolecular collisions with an inert species :math:`\ce{M}` are the primary
means by which an isomer molecule changes its energy. A reasonable estimate -- 
although generally a bit of an underestimate -- of the total rate of collisions
:math:`k_{\mathrm{coll},i}(T)` for each isomer :math:`i` comes from 
Lennard-Jones collision theory:

.. math:: k_{\mathrm{coll},i}(T) = \sqrt{\frac{8 k_\mathrm{B} T}{\pi \mu_i}} \pi d_i^2 \Omega_i^{(2,2)\ast}

Above, :math:`\mu_i` is the reduced mass, :math:`d_i` is the collision diameter, 
and :math:`k_\mathrm{B}` is the Boltzmann constant. The collision diameter is
generally taken as :math:`d \approx \frac{1}{2} (\sigma_i + \sigma_\ce{M})`,
the arithmetic average of the Lennard-Jones :math:`\sigma` parameter for the
isomer and the bath gas. The parameter :math:`\Omega_i^{(2,2)\ast}` represents 
a configurational integral, which is well-approximated by the expression

.. math:: \Omega_i^{(2,2)\ast} = 1.16145 \tilde{T}^{-0.14874} + 0.52487 e^{-0.7732 \tilde{T}}  + 2.16178 e^{-2.437887 \tilde{T}}

where :math:`\tilde{T} \equiv k_\mathrm{B} T / \sqrt{\epsilon_i \epsilon_\ce{M}}`
is a reduced temperature and :math:`\epsilon_i` is the Lennard-Jones 
:math:`\epsilon` parameter. Note that we have used a geometric average for the
:math:`\epsilon` parameters of the isomer and the bath gas in this expression.
Assuming the total gas concentration to be constant and that the gas is ideal, 
we obtain an expression for the collision frequency :math:`\omega_i(T,P)`, 
which makes explicit the pressure dependence:

.. math:: \omega_i(T,P) = k_{\mathrm{coll},i}(T) \frac{P}{k_\mathrm{B} T}

Now that we have an estimate for the total rate of collisions, we need to
develop a model of the effect that these collisions have on the state of the
isomer distribution. To this end, we define 
:math:`P(E, J, E^\prime, J^\prime)` as the probability of a collision
resulting in a transfer of a molecule from state :math:`(E^\prime, J^\prime)`
to state :math:`(E, J)`. There are two mathematical constraints on
:math:`P(E, J, E^\prime, J^\prime)`. The first of these is normalization:

.. math:: \sum_{J^\prime} (2 J^\prime + 1) \int_0^\infty P(E, J, E^\prime, J^\prime) \ dE^\prime = 1

The second of these is detailed balance, required in order to obtain the
Boltzmann distribution at long times:

.. math:: P(E, J, E^\prime, J^\prime) b(E^\prime, J^\prime) = P(E^\prime, J^\prime, E, J) b(E, J)

.. math:: P_i(E^\prime, J^\prime, E, J) = \frac{\rho_i(E^\prime, J^\prime)}{\rho_i(E, J)} \exp \left(- \frac{E^\prime - E}{k_\mathrm{B} T} \right) P_i(E, J, E^\prime, J^\prime) \hspace{40pt} E < E^\prime

Rather than define models directly for :math:`P(E, J, E^\prime, J^\prime)`, 
we usually eliminate the angular momentum contribution and instead define
:math:`P(E, E^\prime)`. This can be related to 
:math:`P(E, J, E^\prime, J^\prime)` via

.. math:: P(E, J, E^\prime, J^\prime) = P(E, E^\prime) \phi(E, J) = P(E, E^\prime) (2J + 1) \frac{\rho(E, J)}{\rho(E)}

where :math:`\rho(E) \equiv \sum_J (2J+1) \rho(E, J)`.

There are a variety of models used for :math:`P(E, E^\prime)`. 
By far the most common is the single exponential down model

.. math:: P(E, E^\prime) = C(E^\prime) \exp \left( -\frac{E^\prime - E}{\alpha} \right) \hspace{40pt} E < E^\prime

where :math:`C(E^\prime)` is determined from the normalization constraint. Note
that this function has been defined for the deactivating direction 
(:math:`E < E^\prime`) only, as the activating direction (:math:`E > E^\prime`)
is then set from detailed balance. The parameter :math:`\alpha` corresponds to
the average energy transferred in a deactivating collision 
:math:`\left< \Delta E_\mathrm{d} \right>`, which itself is a weak function of
temperature.

Other models for :math:`P(E, J, E^\prime, J^\prime)` include the Gaussian down

.. math:: P(E, E^\prime) = C(E^\prime) \exp \left[ - \frac{(E^\prime - E)^2}{\alpha^2} \right] \hspace{40pt} E < E^\prime

and the double exponential down

.. math:: P(E, E^\prime) = C(E^\prime) \left[ (1 - f) \exp \left( -\frac{E^\prime - E}{\alpha_1} \right) + f \exp \left( -\frac{E^\prime - E}{\alpha_2} \right) \right] \hspace{40pt} E < E^\prime

The parameters for these simple models generally contain so much uncertainty
that more complex functional forms are generally not used.

Reaction Models
===============

Chemical reaction events cause a change in molecular configuration at constant
energy. The rate coefficient for this process must be determined as a function
of energy rather than the usual temperature. Such a quantity is called a
*microcanonical rate coefficient* and written as :math:`k(E, J)`. In the master
equation we will differentiate between microcanonical rate coefficients for
isomerization, dissociation, and association by using different letters:
:math:`k_{ij}(E, J)` for isomerization, :math:`g_{nj}(E, J)` for dissociation, 
and :math:`f_{im}(E, J)` for association. (By convention, we use indices
:math:`i` and :math:`j` to refer to unimolecular isomers, :math:`m` and
:math:`n` to refer to bimolecular reactant and product channels, and, later,
:math:`r` and :math:`s` to refer to energy grains.)

As with collision models, the values of the microcanonical rate coefficients
are constrained by detailed balance so that the proper equilibrium is 
obtained. The detailed balance expressions have the form

.. math:: k_{ij}(E, J) \rho_j(E, J) = k_{ji}(E, J) \rho_i(E, J)

for isomerization and

.. math:: f_{in}(E, J) \rho_n(E, J) = g_{ni}(E, J) \rho_i(E, J)

for association/dissociation, where :math:`\rho_i(E, J)` is the density of 
states of the appropriate unimolecular or bimolecular configuration. 

An alternative formulation incorporates the macroscopic equilibrium coefficient
:math:`K_\mathrm{eq}(T)` and equilibrium distributions :math:`b_i(E, J, T)` at 
each temperature:

.. math:: k_{ij}(E, J) b_j(E, J, T) = K_\mathrm{eq}(T) k_{ji}(E, J) b_i(E, J, T)

for isomerization and

.. math:: f_{in}(E, J) b_n(E, J, T) = K_\mathrm{eq}(T) g_{ni}(E, J) b_i(E, J, T)

for association/dissociation. Note that these two formulations are equivalent if
the molecular degrees of freedom are consistent with the macroscopic 
thermodynamic parameters. There are multiple reasons to use the latter 
formulation:

* Only the density of states of the unimolecular isomers need be computed. This 
  is a result of the assumption of thermalized bimolecular channels, which 
  means that we only need to compute the product :math:`f_{in} b_n`, and not 
  the individual values of :math:`f_{in}` and :math:`b_n`.

* Only the reactive rovibrational modes need be included in the density of 
  states. Missing modes will not affect the observed equilibrium because we are 
  imposing the macroscopic equilibrium via :math:`K_\mathrm{eq}(T)`.

* Constants of proportionality in the density of states become unimportant, as 
  they cancel when taking the ratio :math:`\rho(E,J)/Q(\beta)`. For example, if
  the external rotational constants are unknown then we will include an active
  K-rotor in the density of states; this property means that the rotational
  constant of this active K-rotor cancels and is therefore arbitrary.

There are two common ways of determining values for :math:`k(E, J)`: the 
inverse Laplace transform method and RRKM theory. The latter requires detailed
information about the transition state, while the former only requires the 
high-pressure limit rate coefficient :math:`k_\infty(T)`.

Inverse Laplace Transform 
-------------------------

The microcanonical rate coefficient :math:`k(E)` is related to the canonical
high-pressure limit rate coefficient :math:`k_\infty(T)` via a Boltzmann
averaging

.. math:: k_\infty(T) = \frac{\sum_J \int_0^\infty k(E) \rho(E, J) e^{-\beta E} \ dE}{\sum_J \int_0^\infty \rho(E, J) e^{-\beta E} \ dE}

where :math:`\rho(E, J)` is the rovibrational density of states for the 
reactants and :math:`\beta \equiv (k_\mathrm{B} T)^{-1}`. Neglecting the
angular momentum dependence, the above can be written in terms of Laplace 
transforms as

.. math:: k_\infty(T) = \frac{\mathcal{L} \left[ k(E) \rho(E) \right]} {\mathcal{L} \left[ \rho(E) \right]} = \frac{\mathcal{L} \left[ k(E) \rho(E) \right]} {Q(\beta)}

where :math:`Q(\beta)` is the rovibrational partition function for the 
reactants. The above implies that :math:`E` and :math:`\beta` are the transform
variables. We can take an inverse Laplace transform in order to solve for
:math:`k(E)`:

.. math:: k(E) = \frac{\mathcal{L}^{-1} \left[ k_\infty(\beta) Q(\beta) \right] }{\rho(E)}

Hidden in the above manipulation is the assumption that :math:`k_\infty(\beta)`
is valid over a temperature range from zero to positive infinity.

The most common form of :math:`k_\infty(T)` is the modified Arrhenius expression

.. math:: k(T) = A T^n \exp \left( - \frac{E_\mathrm{a}}{k_\mathrm{B} T} \right)

where :math:`A`, :math:`n`, and :math:`E_\mathrm{a}` are the Arrhenius 
preexpoential, temperature exponent, and activation energy, respectively. For
:math:`n = 0` and :math:`E_\mathrm{a} > 0` the inverse Laplace transform can be
easily evaluated to give

.. math:: k(E) = A \frac{\rho(E - E_\mathrm{a})}{\rho(E)} \hspace{40pt} E > E_\mathrm{a}

We can also determine an expression when :math:`n > 0` and 
:math:`E_\mathrm{a} > 0` using a convolution integral:

.. math:: k(E) = A \frac{\phi(E - E_\mathrm{a})}{\rho(E)} \hspace{40pt} E > E_\mathrm{a}

.. math:: \phi(E) = \mathcal{L}^{-1} \left[ T^n Q(\beta) \right] = \frac{1}{k_\mathrm{B}^n \Gamma(n)} \int_0^E (E - x)^{n-1} \rho(x) \ dx

Finally, for cases where :math:`n < 0` and/or :math:`E_\mathrm{a} < 0` we 
obtain a rough estimate by lumping these contributions into the preexponential
at the temperature we are working at. By redoing this at each temperature being
considered we minimize the error introduced, at the expense of not being able
to identify a single :math:`k(E)`.

RRKM Theory 
-----------

RRKM theory -- named for Rice, Ramsperger, Kassel, and Marcus -- is a 
microcanonical transition state theory. Like canonical transition state theory,
detailed information about the transition state and reactants are required,
e.g. from a quantum chemistry calculation. If such information is available,
then the microcanonical rate coefficient can be evaluated via the equation

.. math:: k(E, J) = \frac{N^\ddagger(E, J)}{h \rho(E, J)}

where :math:`N^\ddagger(E, J)` is the sum of states of the transition state,
:math:`\rho(E, J)` is the density of states of the reactant, and :math:`h` is
the Planck constant. Both the transition state and the reactants have been
referenced to the same zero of energy. The sum of states is related to the 
density of states via

.. math:: N(E, J) = \int_0^E \rho(x, J) \ dx

The angular momentum quantum number dependence can be removed via

.. math:: k(E) = \sum_J (2J+1) k(E, J)

The Full Master Equation
========================

The governing equation for the population distributions :math:`p_i(E, J, t)`
of each isomer :math:`i` and the reactant concentrations 
:math:`y_{n\mathrm{A}}(t)` and :math:`y_{n\mathrm{B}}(t)` combines the
collision and reaction models to give a linear integro-differential equation:

.. math::

    \frac{d}{dt} p_i(E, J, t) &= \omega_i(T, P) \sum_{J^\prime} \int_0^\infty P_i(E, J, E^\prime, J^\prime) p_i(E^\prime, J^\prime, t) \ dE^\prime - \omega_i(T, P) p_i(E, J, t) \\
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
    :math:`y_{n\mathrm{A}}(t)`              Total population of species :math:`\ce{A}_n` in reactant channel :math:`n`
    :math:`\omega_i(T, P)`                  Collision frequency of isomer :math:`i`
    :math:`P_i(E, J, E^\prime, J^\prime)`   Collisional transfer probability from :math:`(E^\prime, J^\prime)` to :math:`(E, J)` for isomer :math:`i`
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

    \frac{d}{dt} p_i(E, t) &= \omega_i(T, P) \int_0^\infty P_i(E, E^\prime) p_i(E^\prime, t) \ dE^\prime - \omega_i(T, P) p_i(E, t) \\
                           &  \mbox{} + \sum_{j \ne i}^{N_\mathrm{isom}} k_{ij}(E) p_j(E, t) - \sum_{j \ne i}^{N_\mathrm{isom}} k_{ji}(E) p_i(E, t) \\
                           &  \mbox{} + \sum_{n=1}^{N_\mathrm{reac}} y_{n\mathrm{A}}(t) y_{n\mathrm{B}}(t) f_{in}(E) b_n(E, t) - \sum_{n=1}^{N_\mathrm{reac} + N_\mathrm{prod}} g_{ni}(E) p_i(E, t) \\
    
    \frac{d}{dt} y_{n\mathrm{A}}(t) = \frac{d}{dt} y_{n\mathrm{B}}(t) &= \sum_{i=1}^{N_\mathrm{isom}} \int_0^\infty g_{ni}(E) p_i(E, t) \ dE \\
                                                                      &  \mbox{} - \sum_{i=1}^{N_\mathrm{isom}} y_{n\mathrm{A}}(t) y_{n\mathrm{B}}(t) \int_0^\infty f_{in}(E) b_n(E, t) \ dE


The equations as given are nonlinear, both due to the presence of the 
bimolecular reactants and because both :math:`\omega_i` and 
:math:`P_i(E, E^\prime)` depend on the composition, which is changing with time.
The rate coefficients can be derived from considering the pseudo-first-order 
situation where :math:`y_{n\mathrm{A}}(t) \ll y_{n\mathrm{B}}(t)`, and all 
:math:`y(t)` are negligible compared to the bath gas :math:`\ce{M}`. From these 
assumptions the changes in :math:`\omega_i`, :math:`P_i(E, E^\prime)`, and all 
:math:`y_{n\mathrm{B}}` can be neglected, which yields a linear equation system.

The Energy-Grained Master Equation
==================================

Except for the simplest of unimolecular reaction networks, both the
one-dimensional and two-dimensional master equation must be solved numerically.
To do this we must discretize and truncate the energy domain into a finite
number of discrete bins called *grains*. This converts the linear 
integro-differential equation into a system of first-order ordinary
differential equations:

.. math::

    \frac{d}{dt} \begin{bmatrix}
    \vector{p}_1 \\
    \vector{p}_2 \\
    \vdots \\
    y_{1\mathrm{A}} \\
    y_{2\mathrm{A}} \\
    \vdots
    \end{bmatrix} = \begin{bmatrix}
    \matrix{M}_1 & \matrix{K}_{12} & \ldots & \matrix{F}_{11} \vector{b}_1 y_{1\mathrm{B}} & \matrix{F}_{12} \vector{b}_2 y_{2\mathrm{B}} & \ldots \\
    \matrix{K}_{21} & \matrix{M}_2 & \ldots & \matrix{F}_{21} \vector{b}_1 y_{1\mathrm{B}} & \matrix{F}_{22} \vector{b}_2 y_{2\mathrm{B}} & \ldots \\
    \vdots & \vdots & \ddots & \vdots & \vdots & \ddots \\
    (\vector{g}_{11})^T & (\vector{g}_{12})^T & \ldots & h_1 & 0 & \ldots \\
    (\vector{g}_{21})^T & (\vector{g}_{22})^T & \ldots & 0 & h_2 & \ldots \\
    \vdots & \vdots & \ddots & \vdots & \vdots & \ddots 
    \end{bmatrix} \begin{bmatrix}
    \vector{p}_1 \\
    \vector{p}_2 \\
    \vdots \\
    y_{1\mathrm{A}} \\
    y_{2\mathrm{A}} \\
    \vdots
    \end{bmatrix}

The diagonal matrices :math:`\matrix{K}_{ij}` and :math:`\matrix{F}_{in}` and 
the vector :math:`\vector{g}_{ni}` contain the microcanonical rate coefficients 
for isomerization, association, and dissociation, respectively:

.. math::
    
    (\matrix{K}_{ij})_{rs} &= \begin{cases}
    \frac{1}{\Delta E_r} \int_{E_r - \Delta E_r/2}^{E_r + \Delta E_r/2} k_{ij}(E) \, dE & r = s \\
    0 & r \ne s
    \end{cases} \\
    (\matrix{F}_{in})_{rs} &= \begin{cases}
    \frac{1}{\Delta E_r} \int_{E_r - \Delta E_r/2}^{E_r + \Delta E_r/2} f_{in}(E) \, dE & r = s \\
    0 & r \ne s
    \end{cases} \\
    (\matrix{g}_{ni})_r &= \frac{1}{\Delta E_r} \int_{E_r - \Delta E_r/2}^{E_r + \Delta E_r/2} g_{ni}(E) \, dE
    
The matrices :math:`\matrix{M}_i` represent the collisional transfer 
probabilities minus the rates of reactive loss to other isomers and to 
reactants and products:

.. math::
    
    (\matrix{M}_i)_{rs} = \begin{cases}
    \omega_i \left[ P_i(E_r, E_r) - 1 \right] - \sum_{j \ne i}^{N_\mathrm{isom}} k_{ij}(E_r) - \sum_{n=1}^{N_\mathrm{reac} + N_\mathrm{prod}} g_{ni}(E_r) & r = s \\
    \omega_i P_i(E_r, E_s) & r \ne s
    \end{cases}

The scalars :math:`h_n` are simply the total rate coefficient for loss of 
reactant channel :math:`n` due to chemical reactions:

.. math:: h_n = - \sum_{i=1}^{N_\mathrm{isom}} \sum_{r=1}^{N_\mathrm{grains}} y_{n\mathrm{B}} f_{in}(E_r) b_n(E_r)

Further Reading
===============

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

