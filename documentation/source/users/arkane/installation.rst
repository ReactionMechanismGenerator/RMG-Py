************
Installation
************

Installing Arkane
=================

Arkane can be obtained by installing the `RMG-Py <https://rmg.mit.edu/>`_ software, which
includes all neccesary dependencies.

Instructions to install RMG-Py can be found at the :ref:`RMG-Py Installation page <installation>`.

Note that you'll need to choose between the Basic User or Developer installation instructions
that are specific to your operating system. Modifying Arkane source code will
require Developer installation. If you are only looking to run the code, the
Basic User installation will work.

Installing Q2DTor
=================

Q2DTor is a software for calculating the partition functions and themodynamic properties of molecular systems with two or more
torsional modes developed by David Ferro Costas (david.ferro@usc.es) and Antonio Fernandez Ramos (qf.ramos@usc.es) at
the Universidade de Santiago de Compostela. Arkane can integrate Q2DTor to compute the quantum mechanical partition function 
of 2D rotors.  

For use of Q2DTor and HinderedRotor2D within Arkane please cite:  
D. Ferro-Costas, M. N. D. S. Cordeiro, D. G. Truhlar, A. Fernández-Ramos, Comput. Phys. Commun. 232, 190-205, 2018.

To install a Q2DTor version compatible with RMG run in the RMG-Py directory::

    make q2dtor


