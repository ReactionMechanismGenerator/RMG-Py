.. _databaseModification:

*********************
Database Modification
*********************


Modifying the Thermo Database
=============================

Creating Thermo Libraries
----------------------------------


Adding Thermo Groups
----------------------------------


Adding Thermo to the Depository
----------------------------------

Modifying the Kinetics Database
===============================


Creating Kinetics Libraries
----------------------------------

Adding New Kinetic Groups and Rate Rules
----------------------------------------

Adding Training Reactions 
----------------------------------------

If you know the kinetics of a specific reaction, rather than a rate rule for a template, you can
add the kinetics to the database training set.  By default, RMG creates new rate rules from this 
training set, which in turn benefits the kinetics of similar reactions.  The new rate rules
are formed by matching the reaction to the most most specific template nodes within
the reaction's respective family.  The new rate rule will have a rank of 3, and s kinetics
will be used if there are no other kinetics for the node, or if the other node has lower rank, (i.e. the same node
selection method applies for rate rules added through training reactions.)  If you do not want the
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
