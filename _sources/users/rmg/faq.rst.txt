.. _faq:

**************************
Frequently Asked Questions
**************************


For any other questions related to RMG and its usage and installation, please
post an issue at https://github.com/ReactionMechanismGenerator/RMG-Py/issues and the RMG
developers will get back to you as soon as we can.  You can also search for your problem on the issues
page to see if there are already solutions in development.  Alternatively, you can email us at
rmg_dev@mit.edu.

Why can't my adjacency lists be read any more?
==============================================

The adjacency list syntax changed in July 2014.
The minimal requirement for most translations is to prefix the number
of unpaired electrons with the letter `u`.

Example old syntax::

    HXD13
    1    C 0       {2,D}
    2    C 0 {1,D} {3,S}
    3    C 0 {2,S} {4,D}
    4    C 0 {3,D} {5,S}
    5 *1 C 0 {4,S} {6,S}
    6 *2 C 0 {5,S}

Example new syntax::

    HXD13
    1    C u0       {2,D}
    2    C u0 {1,D} {3,S}
    3    C u0 {2,S} {4,D}
    4    C u0 {3,D} {5,S}
    5 *1 C u0 {4,S} {6,S}
    6 *2 C u0 {5,S}
    
The new syntax, however, allows much
greater flexibility, including definition of lone pairs, partial charges, 
wildcards, and molecule multiplicities, and was necessary to allow us to 
add Nitrogen chemistry.
See :ref:`rmgpy.molecule.adjlist` for details of the new syntax.