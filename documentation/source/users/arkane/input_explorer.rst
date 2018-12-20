*************************************************************************
Creating Input Files for Automated Pressure Dependent Network Exploration
*************************************************************************

Syntax
======

Network exploration starts from a pressure dependent calculation job so naturally it requires an input file
containing everything required for a pressure-dependent calculation input file.  In addition all species 
blocks must have a structure input and there must be a database block in the job.  At the end of the 
pressure-dependent calculation input file you append an explorer block.  For example ::

    explorer(
        source=['methoxy'],
        explore_tol=(1e-2,'s^-1'),
        energy_tol=1e4,
        flux_tol=1e-6,
    )

The ``source`` is a list containing either the label for a single isomer or the labels corresponding to a bimolecular
source channel.  The network is expanded starting from this starting isomer/channel.   

Network Exploration
===================

The ``explore_tol`` is largest acceptable total rate at which isomers inside the network isomerize to become 
isomers that are not part of the network.  Network expansion is done starting from just the network source using
values from the rest of the Arkane job when available, otherwise from RMG.  It cycles through all of the
temperature and pressure points specified for fitting in the pressure dependence job and checks the total network
leak rate at each one.  Whenever this rate is greater than ``explore_tol`` the outside isomer with the most leak is
added to the network and reacted and the loop is flagged to cycle through all of the temperatures and pressures
again.  Once this loop is finished a network_full.py file is generated in the pdep directory that has the full
explored network.  

Network Reduction
=================

``energy_tol`` and ``flux_tol`` control reduction of the full network.  This is highly recommended as the full model
will often contain many unimportant reactions.  How strict you decide to set these parameters should be related
to how much you trust the thermo and rate information being used in the job.  

``energy_tol`` is related to the maximum difference in E0 allowed between a reaction TS and the source channel/isomer.  
A given reaction is deemed filterable at a given condition if ``energy_tol`` * RT < E0_TS - E0_source.  For example, 
if T = 1000 K and ``energy_tol`` is set to 100 and the source energy is ``100 (kJ mol^-1)``, all reactions with a
transition state energy more than ``100 * 8.314 (J K^-1 mol^-1) * 1000 K + 100 (kJ mol^-1) = 931.4 (kJ mol^-1)`` will
not be considered in the final network.

``flux_tol`` is related to the solution of a steady state problem.  In this problem a constant flux of 1.0 (choice of 
units irrelevant) is applied to the source channel/isomer while all A => B + C reactions are assumed to be irreversible.  
Under these conditions steady state concentrations of each isomer and the rates of every unimolecular reaction can be 
calculated.  Since the input flux to the system is 1.0 to some degree these rates represent the fraction of the input 
flux passing through that reaction.  When the flux through a reaction is less than ``flux_tol`` it is deemed filterable 
at that condition.  

Overall a reaction is removed from the network if it is deemed filterable by both methods (if both specified) at every 
temperature and pressure combination.  After reduction a network_reduced.py file is generated in the pdep directory 
containing the reduced network.  
