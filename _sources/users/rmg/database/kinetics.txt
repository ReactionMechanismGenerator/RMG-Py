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

Families
========
Allowable reactions in RMG are divided up into classes called **reaction families**.
All reactions not listed in a kinetic library have their kinetic parameters 
estimated from the reaction families. 

Each reaction family contains the files:

* groups.py containing the recipe, group definitions, and hierarchical trees
* training.py containing a training set for the family
* rules.py containing kinetic parameters for rules

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

Group Definitions
-----------------
The main section in groups.py are the group definitions. Groups are 
:ref:`adjacency lists <rmgpy.molecule.adjlist>`
that describe structures around the reacting atoms. Between the adjacency
list's index number and atom type, a starred number is inserted if the
atom is a reacting atom.

Because groups typically do not describe entire molecules, atoms may appear to 
be lacking full valency. When this occurs, the omitted bonds are allowed to be 
anything. An example of a primary carbon group from H-Abstraction is shown below.
The adjacency list defined on the left matches any of the three drawn structures
on the right (the numbers correspond to the index from the adjacency list).

.. image:: images/Group.png
	:scale: 70%
	:align: center

New **atom types** are also introduced to describe atoms in group definitions. The 
table below shows all atoms types in RMG.

+----------+-------------------+----------------------------------------------+
|Atom Type |Chemical Element   |Bonding                                       |
+==========+===================+==============================================+
|R         |Any                |No requirements                               |
+----------+-------------------+----------------------------------------------+
|R!H       |Any except hydrogen|No requirements                               |
+----------+-------------------+----------------------------------------------+
|H         |Hydrogen           |No requirements                               |
+----------+-------------------+----------------------------------------------+
|C         |Carbon             |No requirements                               |
+----------+-------------------+----------------------------------------------+
|Cs        |Carbon             |No double or triple bonds                     |
+----------+-------------------+----------------------------------------------+
|Cd        |Carbon             |Exactly one double bond to any non-oxygen atom|
+----------+-------------------+----------------------------------------------+
|Cdd       |Carbon             |Two double bonds to any atoms                 |
+----------+-------------------+----------------------------------------------+
|Ct        |Carbon             |One triple bond                               |
+----------+-------------------+----------------------------------------------+
|CO        |Carbon             |Exactly one double bond to an oxygen atom     |
+----------+-------------------+----------------------------------------------+
|Cb        |Carbon             |Exactly two benzene bonds                     |
+----------+-------------------+----------------------------------------------+
|Cbf       |Carbon             |Three benzene bonds (Fused aromatics)         |
+----------+-------------------+----------------------------------------------+
|O         |Oxygen             |No requirements                               |
+----------+-------------------+----------------------------------------------+
|Os        |Oxygen             |No double bonds                               |
+----------+-------------------+----------------------------------------------+
|Od        |Oxygen             |One double bond                               |
+----------+-------------------+----------------------------------------------+
|Oa        |Oxygen             |No bonds (Oxygen atom)                        |
+----------+-------------------+----------------------------------------------+
|S         |Sulfur             |No requirements                               |
+----------+-------------------+----------------------------------------------+
|Ss        |Sulfur             |No double bond                                |
+----------+-------------------+----------------------------------------------+
|Sd        |Sulfur             |One double bond                               |
+----------+-------------------+----------------------------------------------+
|Sa        |Sulfur             |No bonds (Sulfur atom)                        |
+----------+-------------------+----------------------------------------------+

Additionally, groups can also be defined as unions of other groups. For example,::

	label="X_H_or_Xrad_H",
	group=OR{X_H, Xrad_H}, 
    

Forbidden Groups
----------------
Forbidden groups can be defined to ban structures globally in RMG or to
ban pathways in a specific kinetic family.

Globally forbidden structures will ban all reactions containing either reactants
or products that are forbidden.  These groups are stored in in the file located at
``RMG-database/input/forbiddenStructures.py``. 


To ban certain specific pathways in the kinetics 
families, a `forbidden` group must be created, like the following group
in the ``intra_H_migration`` family ::

    forbidden(
        label = "bridged56_1254",
    group =
    """""""
    1 *1 C 1 {2,S} {6,S}
    2 *4 C 0 {1,S} {3,S} {7,S}
    3    C 0 {2,S} {4,S}
    4 *2 C 0 {3,S} {5,S} {8,S}
    5 *5 C 0 {4,S} {6,S} {7,S}
    6    C 0 {1,S} {5,S}
    7    C 0 {2,S} {5,S}
    8 *3 H 0 {4,S}
    """,
        shortDesc = u"""""",
        longDesc = 
    u"""
    
    """,
    )

Forbidden groups should be placed inside the groups.py file located inside the
specific kinetics family's folder ``RMG-database/input/kinetics/family_name/`` 
alongside normal group entries. The starred atoms in the forbidden group
ban the specified reaction recipe from occurring in matched products and reactants.

Hierarchical Trees
------------------
Kinetic groups are ordered into the nodes of a hierarchical trees which is written 
at the end of groups.py. The root node of each tree is the most general group with 
the reacting atoms required for the family. Descending from the root node are 
more specific groups. Each child node is a subset of the parent node above it.

A simplified example of the trees for H-abstraction is shown below. The indented
text shows the syntax in groups.py and a schematic is given underneath.

.. image:: images/Trees.png
	:align: center

Individual groups only describe part of the reaction. To describe an entire reaction
we need one group from each tree, which we call **node templates** or simply templates. 
(C_pri, O_pri_rad), (H2, O_sec_rad), and (X_H, Y_rad) are all valid examples of templates. 
Templates can be filled in with kinetic parameters from the training set or rules.

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