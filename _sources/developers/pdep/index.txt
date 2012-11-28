**************************
Pressure dependence in RMG
**************************

.. module:: rmgpy.pdep

The :mod:`rmgpy.pdep` subpackage provides functionality for calcuating the
pressure-dependent rate coefficients :math:`k(T,P)` for unimolecular reaction
networks.

A unimolecular reaction network is defined by a set of chemically reactive 
molecular configurations - local minima on a potential energy surface - divided
into unimolecular isomers and bimolecular reactants or products. In our 
vernacular, reactants can associate to form an isomer, while such association
is neglected for products. These configurations are connected by chemical 
reactions to form a network; these are referred to as *path* reactions. The
system also consists of an excess of inert gas M, representing a thermal bath;
this allows for neglecting all collisions other than those between an isomer
and the bath gas. 

An isomer molecule at sufficiently high internal energy can be transformed by a
number of possible events:

* The isomer molecule can collide with any other molecule, resulting in an
  increase or decrease in energy

* The isomer molecule can isomerize to an adjacent isomer at the same energy

* The isomer molecule can dissociate into any directly connected bimolecular
  reactant or product channel

It is this competition between collision and reaction events that gives rise
to pressure-dependent kinetics.
  
.. toctree::
    :maxdepth: 2
    
    collision
    reaction
    configuration
    network
    mastereqn
    methods
