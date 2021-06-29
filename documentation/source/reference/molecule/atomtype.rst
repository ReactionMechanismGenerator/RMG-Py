***********************
rmgpy.molecule.AtomType
***********************

.. _atom-types:

.. autoclass:: rmgpy.molecule.AtomType
    :members:

.. autofunction:: rmgpy.molecule.get_atomtype

The atom type of an atom describes the atom itself and (often) something about
the local bond structure around that atom. This is a useful semantic tool for
accelerating graph isomorphism queries, and a useful shorthand when specifying
molecular substructure patterns via an RMG-style adjacency list.

We define the following basic atom types:

=============== ==============================================================================================================================================================
Atom type       Description
=============== ==============================================================================================================================================================
*General atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``R``           any atom with any local bond structure
``R!H``         any non-hydrogen atom with any local bond structure
*Hydrogen atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``H``           hydrogen atom with up to one single bond
*Carbon atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``C``           carbon atom with any local bond structure
``Ca``          carbon atom with two lone pairs and no bonds
``Cs``          carbon atom with up to four single bonds
``Csc``         charged carbon atom with up to three single bonds
``Cd``          carbon atom with one double bond (not to O or S) and up to two single bonds
``Cdc``         charged carbon atom with one double bond and up to one single bond
``CO``          carbon atom with one double bond to oxygen and up to two single bonds
``CS``          carbon atom with one double bond to sulfur and up to two single bonds
``Cdd``         carbon atom with two double bonds
``Ct``          carbon atom with one triple bond and up to one single bond
``Cb``          carbon atom with up to two benzene bonds and up to one single bond
``Cbf``         carbon atom with three benzene bonds
``C2s``         carbon atom with one lone pair (valance 2) and up to two single bonds
``C2sc``        charged carbon atom with one lone pair (valance 2) and up to three single bonds
``C2d``         carbon atom with one lone pair (valance 2) and one double bond
``C2dc``        charged carbon atom with one lone pair (valance 2), one double bond and up to one single bond
``C2tc``        charged carbon atom with one lone pair (valance 2), one triple bond
*Nitrogen atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``N``           nitrogen atom with any local bond structure
``N0sc``        charged nitrogen atom with three lone pairs (valance 0) with up to one single bond
``N1s``         nitrogen atom with two lone pairs (valance 1) and up to one single bond
``N1sc``        charged nitrogen atom with two lone pairs (valance 1) up to two single bonds
``N1dc``        charged nitrogen atom with two lone pairs (valance 1), one double bond
``N3s``         nitrogen atom with one lone pair (valance 3) with up to three single bonds
``N3d``         nitrogen atom with one lone pair (valance 3), one double bond and up to one single bond
``N3t``         nitrogen atom with one lone pair (valance 3) and one triple bond
``N3b``         nitrogen atom with one lone pair (valance 3) and two benzene bonds
``N5sc``        charged nitrogen atom with no lone pairs (valance 5) with up to four single bonds
``N5dc``        charged nitrogen atom with no lone pairs (valance 5), one double bond and up to two single bonds
``N5ddc``       charged nitrogen atom with with no lone pairs (valance 5) and two double bonds
``N5dddc``      charged nitrogen atom with with no lone pairs (valance 5) and three double bonds
``N5tc``        charged nitrogen atom with with no lone pairs (valance 5), one triple bond and up to one single bond
``N5b``         nitrogen atom with with no lone pairs (valance 5) and two benzene bonds (one of the lone pairs also participates in the aromatic bond) and up to one single bond
``N5bd``        nitrogen atom with with no lone pairs (valance 5), two benzene bonds, and one double bond
*Oxygen atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``O``           oxygen atom with any local bond structure
``Oa``          oxygen atom with three lone pairs and no bonds
``O0sc``        charged oxygen with three lone pairs (valance 0) and up to one single bond
``O0dc``        charged oxygen atom with three lone pairs (valance 0) and one double bond
``O2s``         oxygen atom with two lone pairs (valance 2) and up to two single bonds
``O2sc``        charged oxygen atom with two lone pairs (valance 2) and up to one single bond
``O2d``         oxygen atom with two lone pairs (valance 2) and one doubel bond
``O4sc``        charged oxygen atom with one one pair (valance 4) and up to three single bonds
``O4dc``        charged oxygen atom with one one pair (valance 4), one double bond and up to one single bond
``O4tc``        charged oxygen atom with one one pair (valance 4) and one triple bond
``O4b``         oxygen atom with one one pair (valance 4) and and two benzene bonds
*Silicon atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``Si``          silicon atom with any local bond structure
``Sis``         silicon atom with four single bonds
``Sid``         silicon atom with one double bond (to carbon) and two single bonds
``SiO``         silicon atom with one double bond (to oxygen) and two single bonds
``Sidd``        silicon atom with two double bonds
``Sit``         silicon atom with one triple bond and one single bond
``Sib``         silicon atom with two benzene bonds and one single bond
``Sibf``        silicon atom with three benzene bonds
*Phosphorus atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``P``           phosphorus atom with any local bond structure
``P0sc``        charged phosphorus atom with three lone pairs (valence 0) and up to 1 single bond
``P1s``         phosphorus atom with two lone pairs (valence 1) and up to 1 single bond
``P1sc``        charged phosphorus atom with two lone pairs (valence 1) and up to 2 single bonds
``P1dc``        charged phosphorus atom with two lone pairs (valence 1) and 1 double bond
``P3s``         phosphorus atom with one lone pair (valence 3) and up to 3 single bonds
``P3d``         phosphorus atom with one lone pair (valence 3), 1 double bond and up to 1 single bond
``P3t``         phosphorus atom with one lone pair (valence 3) and 1 triple bond
``P3b``         phosphorus atom with one lone pair (valence 3) and 2 benzene bonds
``P5s``         phosphorus atom with no lone pairs (valence 5) and up to 5 single bonds
``P5sc``        charged phosphorus atom with no lone pairs (valence 5) and up to 6 single bonds
``P5d``         phosphorus atom with no lone pairs (valence 5), 1 double bond and up to 3 single bonds
``P5dd``        phosphorus atom with no lone pairs (valence 5), 2 double bonds and up to 1 single bond
``P5dc``        charged phosphorus atom with no lone pairs (valence 5), 1 double bond and up to 2 single bonds
``P5ddc``       charged phosphorus atom with no lone pairs (valence 5) and 2 double bonds
``P5t``         phosphorus atom with no lone pairs (valence 5), 1 triple bond and up to 2 single bonds
``P5td``        phosphorus atom with no lone pairs (valence 5), 1 triple bond and 1 double bond
``P5tc``        charged phosphorus atom with no lone pairs (valence 5), 1 triple bond and up to 1 single bond
``P5b``         phosphorus atom with no lone pairs (valence 5), 2 benzene bonds and up to 1 single bond
``P5bd``        phosphorus atom with no lone pairs (valence 5), 2 benzene bonds and 1 double bond
*Sulfur atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``S``           sulfur atom with any local bond structure
``Sa``          sulfur atom with three lone pairs and no bonds
``S0sc``        charged sulfur atom with three lone pairs (valance 0) and up to one single bonds
``S2s``         sulfur atom with two lone pairs (valance 2) and up to two single bonds
``S2sc``        charged sulfur atom with two lone pairs (valance 2) and up to three single bonds
``S2d``         sulfur atom with two lone pairs (valance 2) and one double bond
``S2dc``        charged sulfur atom with two lone pairs (valance 2), one double bond and up to one single bond
``S2tc``        charged sulfur atom with two lone pairs (valance 2) and one triple bond
``S4s``         sulfur atom with one lone pair (valance 4) and up to four single bonds
``S4sc``        charged sulfur atom with one lone pair (valance 4) and up to five single bonds
``S4d``         sulfur atom with one lone pair (valance 4), one double bond and up to two single bonds
``S4dd``        sulfur atom with one lone pair (valance 4) and two double bonds
``S4dc``        charged sulfur atom with one lone pair (valance 4), one to three double bonds and up to three single bonds
``S4b``         sulfur atom with one lone pair (valance 4) and two benzene bonds (one of the lone pairs also participates in the aromatic bond)
``S4t``         sulfur atom with one lone pair (valance 4), one triple bond and up to one single bond
``S4tdc``       charged sulfur atom with one lone pair (valance 4) one to two triple bonds, up to two double bonds, and up to three single bonds
``S6s``         sulfur atom with no lone pairs (valance 6) and up to six single bonds
``S6sc``        charged sulfur atom with no lone pairs (valance 6) and up to seven single bonds
``S6d``         sulfur atom with no lone pairs (valance 6), one double bond and up to four single bonds
``S6dd``        sulfur atom with no lone pairs (valance 6), two double bonds and up to two single bonds
``S6ddd``       sulfur atom with no lone pairs (valance 6) and three double bonds
``S6dc``        charged sulfur atom with no lone pairs (valance 6), one to three double bonds and up to five single bonds
``S6t``         sulfur atom with no lone pairs (valance 6), one triple bond and up to three single bonds
``S6td``        sulfur atom with no lone pairs (valance 6), one triple bond, one double bond and up to one single bond
``S6tt``        sulfur atom with no lone pairs (valance 6) and two triple bonds
``S6tdc``       charged sulfur atom with no lone pairs (valance 6), one to two triple bonds, up to two double bonds, and up to four single bonds
*Chlorine atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``Cl``          chlorine atom with any local bond structure
``Cl1s``        chlorine atom with three lone pairs and zero to one single bonds
*Bromine atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``Br``          bromine atom with any local bond structure
``Br1s``        bromine atom with three lone pairs and zero to one single bonds
*Iodine atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``I``           iodine atom with any local bond structure
``I1s``         iodine atom with three lone pairs and zero to one single bonds
*Fluorine atom types*
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``F``           fluorine atom with any local bond structure
``F1s``         fluorine atom with three lone pairs and zero to one single bonds
=============== ==============================================================================================================================================================
