.. _ratebasedmodelenlarger:

Rate-Based Model Enlarger
=========================

To construct a mechanism, the user must specify an initial set of species and
the initial conditions (temperature, pressure, species concentrations, etc.).
RMG reacts the initial species in all possible ways according to its known
reaction families, and it integrates the model in time. RMG tracks the rate
(flux) at which each new "edge" species is produced, and species (and the
reactions producing them) that are produced with significant fluxes are
incorporated into the model (the "core"). These new core species are reacted
with all other core species in the model, to generate a new set of edge species
and reactions. The time-integration restarts, and the expanded list of edge
species is monitored for significant species to be included in the core. The
process continues until all significant species and reactions have been
included in the model. The user is free to vary the definition of a
"significant" rate to refine the mechanism as desired. For a more detailed
description on rate-based model enlargement, please see [Susnow1997]_.

.. [Susnow1997] \ R. G. Susnow, A. M. Dean, W. H. Green, P. K. Peczak, and L. J. Broadbelt. "Rate-Based Construction of Kinetic Models for Complex Systems." *J. Phys. Chem. A* **101**, p. 3731-3740 (1997).