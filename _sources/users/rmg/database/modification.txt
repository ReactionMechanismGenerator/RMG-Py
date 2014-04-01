.. _databaseModification:

*********************
Database Modification
*********************
Note that the RMG-Py database is written in Python code where line indentions
determine the scope. When modifying the database, be sure to preserve all 
line indentions shown in the examples.

Modifying the Thermo Database
=============================

Creating Thermo Libraries
-------------------------


Adding Thermo Groups
--------------------


Adding Thermo to the Depository
-------------------------------

.. _kinetic-database-modification:

Modifying the Kinetics Database
===============================

For the casual user, it is recommended to use either a kinetic library or 
add to the training set instead of modifying the kinetic groups. 

Put kinetic parameters into a kinetic library when:

* A set of reaction rates were optimized together
* You know the reaction rate is not generalizable to similar species (perhaps due to catalysis or aromatic structures)
* No family exists for the class of reaction
* You are not confident about the accuracy of kinetic parameters

Put kinetic parameters into the training set when:

* You are confident on the accurcy of the kinetic parameter
* You wish for the reaction to be generalized to similar reactions in your mechanism

Creating Kinetics Libraries
---------------------------

Adding New Kinetic Groups and Rate Rules
----------------------------------------

Decide on a Template
--------------------
First you need to know the template for your reaction to decide whether or not
to create new groups: 

#. Type your reaction into the kinetics search at http://rmg.mit.edu/database/kinetics/search/
#. Select the correct reaction
#. In the results search for "(RMG-Py rate rules)" and select that link. The kinetic family listed is the family of interest.
#. Scroll to the bottom and look at the end of the long description. There may be very long description of the averaging scheme, but the template for the reaction is the very last one listed:

.. image:: images/GroupSearch.png
	:align: center

Now you must determine whether the chosen template is appropriate.  A good rule
of thumb is to see if the all neighbours of the reacting atoms are as specified
as possible. For example, assume your species is ethanol

.. image:: images/ethanol.png
	:scale: 150%
	:align: center

and RMG suggests the group::

	label = "C_sec",
	group = 
	"""
	1 *1 Cs  0 {2,S} {3,S} {4,S}
	2 *2 H   0 {1,S}
	3    R!H 0 {1,S}
	4    R!H 0 {1,S}
	""",

If you use the suggested groups you will not capture the effect of the alcohol 
group. Therefore it is better to make a new group. ::

	label = "C/H2/CsO",
	group = 
	"""
	1 *1 Cs  0 {2,S} {3,S} {4,S} {5,S}
	2 *2 H  0 {1,S}
	3    H  0 {1,S}
	4    O  0 {1,S}
	5    Cs 0 {1,S}
	""",

If you have determined the suggested groups is appropriate, skip to 
:ref:`kinetic-training-set` or :ref:`kinetic-rules`. Otherwise proceed to the 
next section for instructions on creating the new group.

Creating a New Group
--------------------

In the family's groups.py, you will need to add an entry of the format::

	entry(
		index = 61,
		label = "C_sec",
		group = 
	"""
	1 *1 Cs   0 {2,S} {3,S} {4,S} {5,S}
	2 *2 H   0 {1,S}
	3    C   0 {1,S}
	4    H   0 {1,S}
	5    R!H 0 {1,S}
	""",
		kinetics = None,
		reference = None,
		referenceType = "",
		shortDesc = u"""""",
		longDesc = u"""""",
	)

* The index can be any number not already present in the set
* The label is the name of the group.
* The group is the group adjacency list with the starred reacting atoms.
* The other attributes do not need to be filled for a group

Next, you must enter your new group into the tree. At the bottom of groups.py
you will find the trees. Place your group in the appropriate position. In the 
example given in the previous section, the new group would be added under the C_sec. ::

	L1: X_H
		L2: H2
		L2: Cs_H
			L3: C_pri
			L3: C_sec
				L4: C/H2/CsO
			L3: C_ter

.. _kinetic-rules:
			
Adding Kinetic Rules
--------------------
Rules give generalized kinetic parameters for a specific node template. In most
cases, your kinetic parameters describe a specific reaction in which case you
will want to add your reaction to the training set.
 
The rule must be added into rules.py in the form::

	entry(
		index = 150,
		label = "C/H/Cs3;O_rad/NonDeO",
		group1 = 
	"""
	1 *1 Cs  0 {2,S} {3,S} {4,S} {5,S}
	2 *2 H  0 {1,S}
	3    Cs 0 {1,S}
	4    Cs 0 {1,S}
	5    Cs 0 {1,S}
	""",
		group2 = 
	"""
	1 *3 O 1 {2,S}
	2    O 0 {1,S}
	""",
		kinetics = ArrheniusEP(
			A = (2800000000000.0, 'cm^3/(mol*s)', '*|/', 5),
			n = 0,
			alpha = 0,
			E0 = (16.013, 'kcal/mol', '+|-', 1),
			Tmin = (300, 'K'),
			Tmax = (1500, 'K'),
		),
		reference = None,
		referenceType = "",
		rank = 5,
		shortDesc = u"""Curran et al. [8] Rate expressions for H atom abstraction from fuels.""",
		longDesc = 
	u"""
	[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
	Rate expressions for H atom abstraction from fuels.

	pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:HO2, Site: tertiary (c)
	
	Verified by Karma James
	""",
	) 

* The index can be any number not already used in rules.py.
* The label is the name of the rule.
* The groups must have the adjacency list of the respective groups. Between them they should have all starred atoms from the recipe.
* The value and units of kinetic parameters must be given. 
	* Multiplicative uncertainty is given as "'*\|/,' 5" meaning within a factor of 5 
	* Additive uncertainty is given as "'+\|/-', 2" meaning plus or minus 2.
* Rank determines the priority of the rule when compared with other rules.
* The short description will appear in the annotated chemkin file.
* The long description only appears in the database.

.. _kinetic-training-set:

Adding Training Reactions
-------------------------

If you know the kinetics of a specific reaction, rather than a rate rule for a template, you can
add the kinetics to the database training set.  By default, RMG creates new rate rules from this 
training set, which in turn benefits the kinetics of similar reactions.  The new rate rules
are formed by matching the reaction to the most most specific template nodes within
the reaction's respective family. If you do not want the
training depository reactions to create new rate rules in the database, set the option for 
``kineticsDepositories`` within the ``database`` field in your input file to ::

    kineticsDepositories = ['!training'],


Currently, RMG's rate rule estimates overrides all kinetics depository kinetics, including training
reactions.  Unless the training reaction's rate rule ranks higher than the existing node, it 
will not be used.  If you want the training reaction to override the rate rule estimates, you should put the reaction into
a reaction library or seed mechanism.  

The easiest way to add training reactions to the database is via the RMG website.  First, search for 
the reaction using http://rmg.mit.edu/database/kinetics/search/ . This will automatically search 
the existing RMG database for the reaction, as well as identify the reaction family template
that this reaction matches.  If the reaction does not match any family, then it cannot be added to the 
training reactions.  Click the 'Create training rate from average' button underneath the kinetics plot 
for the reaction and edit the kinetics and reference descriptions for the reaction.  The atom labels
marking the reaction recipe actions (lose bond, add radical, etc.) will already be automatically 
labeled for you.  After editing the reaction data, write a short message for the reaction added under 
the 'Summary of changes' field, then click 'Save.'  You will need an account for the RMG website to 
make an entry.

.. note::

	If you are entering the reaction in the reverse direction of the family, you must still label the
	reactants and products with the atomLabels of the original reaction template.  Otherwise, RMG
	will not be able to locate the nodes in the group tree to match the reaction.
	
	Entries added in the reverse direction of the original template will use the current RMG job's 
	thermo database	to estimate the kinetics in the forward direction.  Therefore this value can differ
	depending on the order of thermo libraries used when running a job.
 
If adding the training reaction manually, first identify the reaction family of the reaction, then
go to the family's folder in ``RMG-database/input/kinetics/families/``.  Create a new kinetics entry
in the ``training.py`` file.  Make sure to apply the reaction recipe labels properly for the
reactants and products.

Pitfalls
--------
Be careful with the specificity when naming neighbouring atoms. On upper nodes,
you should try to be general so that you do not exclude reactions. 

Sibling nodes must be exclusive from one another so that there is no question
which group a molecule qualifies as. However, you do not need to be exhaustive and
list out every possibility.

Make sure your nodes are actually children of their parents. Currently RMG does
no atom-by-atom checking and assumes whatever is put into the tree is correct.

Be sure to give errors whenever adding rules. If you don't know the uncertainty,
why do you trust the kinetics?

After you are done always check via populate reactions or the website, that your
modifications are behaving the way you expect.

Caveat regarding how rate rules are used by RMG and the rate parameters you input: because tunneling is
important for many chemical reactions, the rate of a reaction may not be easily represented by
a bi-Arrhenius fit. 3-parameter fits are more common. However, the resulting fit may report an
'activation energy' that is much different (possibly by 10+ kcals) than the the true barrier height. 
When RMG is assembling pressure-dependent networks, it will use barrier heights from rate rules. This can 
lead to very inaccurate rate calculations. To avoid this issue, try to ensure that your fitted arrhenius 
activation energy truly does reflect the reaction barrier height. 