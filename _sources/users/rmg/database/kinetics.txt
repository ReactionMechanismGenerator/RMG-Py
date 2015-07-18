.. _kineticsDatabase:

*****************
Kinetics Database
*****************
This section describes the general usage of RMG's kinetic database. See :ref:`kinetic-database-modification` for 
instructions on modifying the database.

Pressure independent reaction rates in RMG are calculated using a modified 
Arrhenius equation, designating the reaction coefficient as :math:`k(T)` at 
temperature :math:`T`.

.. math:: k(T) = A\left(\frac{T}{T_0}\right)^ne^{-(E_a + \alpha \Delta H_{rxn})/(RT)}

:math:`R` is the universal gas constant. The **kinetic parameters** determining 
the rate coefficient are:

* :math:`A`:	the pre-exponential A-factor 

* :math:`T_0`:	the reference temperature

* :math:`n`:	the temperature exponent

* :math:`E_a`:	the activation energy 

* :math:`\alpha`:	the Evans-Polanyi coefficient

* :math:`\Delta H_{rxn}`: the enthalpy of reaction

When Evans-Polanyi corrections are used, :math:`\Delta H_{rxn}` is calculated
using RMG's thermo database, instead of being specified in the kinetic database.  

Libraries
=========
Kinetic libraries delineate kinetic parameters for specific reactions. 
RMG always chooses to use kinetics from libraries over families. If multiple libraries
contain the same reaction, then precedence is given to whichever library is
listed first in the input.py file.

For combustion mechanisms, you should always use *one* small-molecule 
combustion library, such as the pre-packaged ERC-Foundation Fuel. 
The reactions contained in these libraries are poorly estimated by kinetic 
families and are universally important to combustion systems.

Kinetic libraries should also be used in the cases where:

* A set of reaction rates were optimized together
* You know the reaction rate is not generalizable to similar species (perhaps due to catalysis or aromatic structures)
* No family exists for the class of reaction
* You are not confident about the accuracy of kinetic parameters

.. _kineticsFamilies:

Families
========
Allowable reactions in RMG are divided up into classes called **reaction families**.
All reactions not listed in a kinetic library have their kinetic parameters 
estimated from the reaction families. 

Each reaction family contains the files:

* groups.py containing the recipe, group definitions, and hierarchical trees
* training.py containing a training set for the family
* rules.py containing kinetic parameters for rules

There are currently 45 reaction families in RMG:

**1+2_Cycloaddition**     

.. image:: images/kinetics_families/1+2_Cycloaddition.png 
	:scale: 40% 

**1,2-Birad_to_alkene**     

.. image:: images/kinetics_families/1,2-Birad_to_alkene.png 
	:scale: 40% 

**1,2_Insertion_carbene**     

.. image:: images/kinetics_families/1,2_Insertion_carbene.png 
	:scale: 40%  

**1,2_Insertion_CO**     

.. image:: images/kinetics_families/1,2_Insertion_CO.png 
	:scale: 40% 

**1,2_shiftS**     

.. image:: images/kinetics_families/1,2_shiftS.png 
	:scale: 40% 

**1,3_Insertion_CO2**     

.. image:: images/kinetics_families/1,3_Insertion_CO2.png 
	:scale: 40% 

**1,3_Insertion_ROR**     

.. image:: images/kinetics_families/1,3_Insertion_ROR.png 
	:scale: 40% 

**1,3_Insertion_RSR**     

.. image:: images/kinetics_families/1,3_Insertion_RSR.png 
	:scale: 40% 

**1,4_Cyclic_birad_scission**     

.. image:: images/kinetics_families/1,4_Cyclic_birad_scission.png 
	:scale: 40% 

**1,4_Linear_birad_scission**     

.. image:: images/kinetics_families/1,4_Linear_birad_scission.png 
	:scale: 40% 

**2+2_cycloaddition_CCO**     

.. image:: images/kinetics_families/2+2_cycloaddition_CCO.png 
	:scale: 40% 

**2+2_cycloaddition_Cd**     

.. image:: images/kinetics_families/2+2_cycloaddition_Cd.png 
	:scale: 40% 

**2+2_cycloaddition_CO**     

.. image:: images/kinetics_families/2+2_cycloaddition_CO.png 
	:scale: 40% 

**Birad_recombination**     

.. image:: images/kinetics_families/Birad_recombination.png 
	:scale: 40% 

**Cyclic_Ether_Formation**     

.. image:: images/kinetics_families/Cyclic_Ether_Formation.png 
	:scale: 40% 

**Diels_alder_addition**     

.. image:: images/kinetics_families/Diels_alder_addition.png 
	:scale: 40% 

**Disproportionation**     

.. image:: images/kinetics_families/Disproportionation.png 
	:scale: 40% 

**H_Abstraction**     

.. image:: images/kinetics_families/H_Abstraction.png 
	:scale: 40% 

**H_shift_cyclopentadiene**     

.. image:: images/kinetics_families/H_shift_cyclopentadiene.png 
	:scale: 40% 

**HO2_Elimination_from_PeroxyRadical**     

.. image:: images/kinetics_families/HO2_Elimination_from_PeroxyRadical.png 
	:scale: 40% 

**Intra_Diels_alder**     

.. image:: images/kinetics_families/Intra_Diels_alder.png 
	:scale: 40% 

**Intra_Disproportionation**     

.. image:: images/kinetics_families/Intra_Disproportionation.png 
	:scale: 40% 

**intra_H_migration**     

.. image:: images/kinetics_families/intra_H_migration.png 
	:scale: 40% 

**intra_NO2_ONO_conversion**     

.. image:: images/kinetics_families/intra_NO2_ONO_conversion.png 
	:scale: 40% 

**intra_OH_migration**     

.. image:: images/kinetics_families/intra_OH_migration.png 
	:scale: 40% 

**Intra_R_Add_Endocyclic**     

.. image:: images/kinetics_families/Intra_R_Add_Endocyclic.png 
	:scale: 40% 

**Intra_R_Add_Exocyclic**     

.. image:: images/kinetics_families/Intra_R_Add_Exocyclic.png 
	:scale: 40% 

**Intra_R_Add_ExoTetCyclic**     

.. image:: images/kinetics_families/Intra_R_Add_ExoTetCyclic.png 
	:scale: 40% 

**Intra_RH_Add_Endocyclic**     

.. image:: images/kinetics_families/Intra_RH_Add_Endocyclic.png 
	:scale: 40% 

**Intra_RH_Add_Exocyclic**     

.. image:: images/kinetics_families/Intra_RH_Add_Exocyclic.png 
	:scale: 40% 

**intra_substitutionCS_cyclization**     

.. image:: images/kinetics_families/intra_substitutionCS_cyclization.png 
	:scale: 40% 

**intra_substitutionCS_isomerization**     

.. image:: images/kinetics_families/intra_substitutionCS_isomerization.png 
	:scale: 40% 

**intra_substitutionS_cyclization**     

.. image:: images/kinetics_families/intra_substitutionS_cyclization.png 
	:scale: 40% 

**intra_substitutionS_isomerization**     

.. image:: images/kinetics_families/intra_substitutionS_isomerization.png 
	:scale: 40% 

**ketoenol**     

.. image:: images/kinetics_families/ketoenol.png 
	:scale: 40% 

**Korcek_step1**     

.. image:: images/kinetics_families/Korcek_step1.png 
	:scale: 40% 

**Korcek_step2**     

.. image:: images/kinetics_families/Korcek_step2.png 
	:scale: 40% 

**lone_electron_pair_bond**     

.. image:: images/kinetics_families/lone_electron_pair_bond.png 
	:scale: 40% 

**Oa_R_Recombination**     

.. image:: images/kinetics_families/Oa_R_Recombination.png 
	:scale: 40% 

**R_Addition_COm**     

.. image:: images/kinetics_families/R_Addition_COm.png 
	:scale: 40% 

**R_Addition_CSm**     

.. image:: images/kinetics_families/R_Addition_CSm.png 
	:scale: 40% 

**R_Addition_MultipleBond**     

.. image:: images/kinetics_families/R_Addition_MultipleBond.png 
	:scale: 40% 

**R_Recombination**     

.. image:: images/kinetics_families/R_Recombination.png 
	:scale: 40% 

**Substitution_O**     

.. image:: images/kinetics_families/Substitution_O.png 
	:scale: 40% 

**SubstitutionS**     

.. image:: images/kinetics_families/SubstitutionS.png 
	:scale: 40% 




Recipe
------
The recipe can be found near the top of groups.py and describes the changes in
bond order and radicals that occur during the reaction. Reacting atoms are
labelled with a starred number. Shown below is the recipe for the H-abstraction 
family.

.. image:: images/Recipe.png
	:scale: 65%
	:align: center

The table below shows the possible actions for recipes. The arguments are given 
in the curly braces as shown above. For the order of bond change in the 
Change_Bond action, a -1 could represent a triple bond changing to a double 
bond while a +1 could represent a single bond changing to a double bond. 

+------------+-----------------+---------------------+------------------+
|Action      |Argument1        |Argument2            |Argument3         |
+============+=================+=====================+==================+
|Break_Bond  |First bonded atom|Type of bond         |Second bonded atom|
+------------+-----------------+---------------------+------------------+
|Form_Bond   |First bonded atom|Type of bond         |Second bonded atom|
+------------+-----------------+---------------------+------------------+
|Change_Bond |First bonded atom|Order of bond change |Second bonded atom|
+------------+-----------------+---------------------+------------------+
|Gain_Radical|Specified atom   |Number of radicals   |                  |
+------------+-----------------+---------------------+------------------+
|Lose_Radical|Specified atom   |Number of radicals   |                  |
+------------+-----------------+---------------------+------------------+

Change_Bond order cannot be directly used on benzene bonds. During generation,
aromatic species are kekulized to alternating double and single bonds such that
reaction families can be applied. However, RMG cannot properly handle benzene bonds 
written in the kinetic group definitions.

Training Set vs Rules
---------------------
The training set and rules both contain trusted kinetics that are used to fill in
templates in a family. The **training set** contains kinetics for specific reactions,
which are then matched to a template. The kinetic **rules** contain kinetic 
parameters that do not necessarily correspond to a specific reaction, but have 
been generalized for a template.

When determining the kinetics for a reaction, a match for the template
is searched for in the kinetic database. The three cases in order
of decreasing reliability are:

#. Reaction match from training set
#. Node template exact match using either training set or rules
#. Node template estimate averaged from children nodes

The reaction match from training set is accurate within the documented uncertainty for that
reaction. A template exact match is usually accurate within about one order
of magnitude. When there is no kinetics available for for the template in
either the training set or rules, the kinetics are averaged from the children
nodes as an estimate. In these cases, the kinetic parameters are much less reliable.
For more information on the estimation algorithm see :ref:`kinetics`. 

The training set can be modified in training.py and the rules can be modified in
rules.py. For more information on modification see :ref:`kinetic-training-set` and :ref:`kinetic-rules`.
