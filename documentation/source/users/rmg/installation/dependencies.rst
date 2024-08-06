.. _dependencies:

************
Dependencies
************


List of Dependencies
====================

A list of RMG's dependencies can be found in the ``environment.yml`` file on the `RMG GitHub page: <https://github.com/ReactionMechanismGenerator/RMG-Py/blob/main/environment.yml>`_. 


.. _dependenciesRestrictions:

License Restrictions on Dependencies
====================================

All of RMG's dependencies except the ones listed below are freely available and compatible with RMG's open source MIT license (though the specific nature of their licenses vary). 

* **pydas**: The DAE solvers used in the simulations come from `Linda Petzold's research group <https://cse.cs.ucsb.edu/software/>`_ at UCSB.  For running sensitivity analysis in RMG, the DASPK 3.1 solver is required, which "is subject to copyright restrictions” for non-academic use. Please visit their website for more details. To run RMG without this restriction, one may switch to compiling with the DASSL solver instead in RMG, which is "available in the public domain.”
