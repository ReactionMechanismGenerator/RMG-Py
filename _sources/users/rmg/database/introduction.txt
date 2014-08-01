.. _introDatabase:

************
Introduction
************
This section describes some of the general characteristics of RMG's databases.

Group Definitions
-----------------
The main section in many of RMG's databases are the 'group' definitions. Groups are 
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
Groups are ordered into the nodes of a hierarchical trees which is written 
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