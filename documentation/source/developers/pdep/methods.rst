************************************
Methods for estimating k(T,P) values
************************************

The objective of each of the methods described in this section is to reduce 
the master equation into a small number of phenomenological rate coefficients 
:math:`k(T,P)`. All of the methods share a common formalism in that they seek
to express the population distribution vector :math:`\mathbf{p}_i` for each
unimolecular isomer :math:`i` as a linear combination of the total populations
of all unimolecular isomers and bimolecular reactant channels.

The modified strong collision method
====================================

The modified strong collision method utilizes a greatly simplified collision
model that allows for a decoupling of the energy grains. In the simplified
collision model, collisional stabilization of a reactive isomer is treated as a
single-step process, ignoring the effects of collisional energy redistribution
within the reactive energy space. An attempt to correct for the effect of
collisional energy redistribution is made by modifying the collision frequency
:math:`\omega_i(T,P)` with a collision efficiency :math:`\beta_i(T)` estimated
from the low-pressure limit fall-off of a single isomer. 

By approximating the reactive populations as existing in pseudo-steady state, 
the master equation is converted to a matrix equation is at each energy. Solving
these small matrix equations gives the pseudo-steady state populations of each
isomer as a function of the total population of each isomer and reactant 
channel, which are then applied to determine the :math:`k(T,P)` values.

In practice, the modified strong collision method is the fastest and most
robust of the methods, and is reasonably accurate over a wide range of
temperatures and pressures.

.. autofunction:: rmgpy.pdep.msc.applyModifiedStrongCollisionMethod

The reservoir state method
==========================

In the reservoir state method, the population distribution of each isomer is
partitioned into the low-energy grains (called the *reservoir*) and the
high-energy grains (called the *active space*). The partition generally occurs
at or near the lowest transition state energy for each isomer. The reservoir
population is assumed to be thermalized, while the active-space population is
assumed to be in pseudo-steady state. Applying these approximations converts
the master equation into a single large matrix equation. Solving this matrix 
equation gives the pseudo-steady state populations of each isomer as a function 
of the total population of each isomer and reactant channel, which are then
applied to determine the :math:`k(T,P)` values.

The reservoir state method is only slightly more expensive than the modified
strong collision method. At low temperatures the approximations used are very
good, and the resulting :math:`k(T,P)` values are more accurate than the
modified strong collision values. However, at high temperatures the thermalized
reservoir approximation breaks down, resulting in very inaccurate 
:math:`k(T,P)` values. Thus, the reservoir state method is not robustly
applicable over a wide range of temperatures and pressures.

.. autofunction:: rmgpy.pdep.rs.applyReservoirStateMethod

The chemically-significant eigenvalues method
=============================================

In the chemically-significant eigenvalues method, the master equation matrix
is diagonized to determine its eigenmodes. Only the slowest of these modes are
relevant to the chemistry; the rest involve internal energy relaxation due to
collisions. Keeping only these "chemically-significant" eigenmodes allows for
reduction to :math:`k(T,P)` values.

The chemically-significant eigenvalues method is the most accurate method,
and is considered to be exact as long as the chemically-significant eigenmodes
are separable and distinct from the internal energy relaxation eigenmodes.
However, this is often only the case near the high-pressure limit, even for
networks of only modest size. The chemically-significant eigenvalues method
is also substantially more expensive to apply than the other methods.

.. autofunction:: rmgpy.pdep.cse.applyChemicallySignificantEigenvaluesMethod

