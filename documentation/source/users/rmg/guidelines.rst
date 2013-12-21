.. _guidelines:

*******************************
Guidelines for Building a Model
*******************************

RMG has been designed to build kinetic models for gas phase pyrolysis and combustion of 
organic molecules made of C, H, O and S. By kinetic model, we mean a set 
of reactions and associated kinetics that represent the chemical transformations 
occurring in the system of interest. These systems could be the combustion of fuels, 
pyrolysis of hydrocarbon feedstocks, etc. The total number of reactions and species 
typically required to describe some of these processes can run into the thousands making 
these models difficult and error-prone to build manually. This is the main motivation 
behind using software like RMG that build such models automatically in a systematic 
reproducible manner.

In RMG, the user is expected to provide an input file specifying the conditions
(temperature, pressure, etc.) under which one desires to develop kinetic models. 
The following are some tips for setting up your input/condition file.


Start with a good seed mechanism
--------------------------------
RMG is a useful tool in elucidating important pathways in a given process but may not 
capture certain special reaction types which may be specific to the system you are 
interested in. However, if you already have a good idea of these reactions that are 
important and are not available in the standard RMG library, you can create a 
‘seed mechanism’  and include it in the input file to RMG. This will directly 
include these in the model core and add other reactions from the RMG library on 
top of it using our rate based algorithm. (Similarly, you can specify your own 
thermodynamic parameters for species using thermochemistry libraries which are similar 
in concept to seed mechanisms. In order to build these libraries, you will need to 
specify all species in the RMG adjacency list format.) In a combustion system, RMG 
tends to do a decent job filling in the termination and propagation steps of a mechanism 
if it is guided with the initiation and chain branching steps using a seed mechanism. 
Ideally, RMG should be able to find all the right chemistry through our kinetics 
database but holes in current kinetic databases can make this task difficult. A 
good seed mechanism can address this issue for the system of interest and also 
reduce the size, cost and time taken to arrive at a converged model.

Setting up the right termination criterion
------------------------------------------
Start with a relatively large tolerance (such as 0.1) when building your first model 
to make sure that RMG can converge the model to completion without any hiccups, then 
begin tightening the tolerance if you are able to converge the initial model.  For 
large molecules such as tetradecane (C14), even a tolerance of 0.1 may be too tight 
for RMG to work with and lead to convergence problems. Note that a good seed mechanism 
allows for faster convergence.

Restricting the number of carbon atoms, oxygen atoms, and radical sites per species
-----------------------------------------------------------------------------------
Options to tune the maximum number of carbon or oxygen atoms, or number of radical 
sites per species can be specified at the beginning of the condition file. In most 
systems, we do not expect large contributions from species with more than 1 radical 
center (i.e. biradicals, etc.) to affect the overall chemistry, thus it may be useful 
to limit the maximum number of radicals to 2 (to allow for O2). The same applies for 
the maximum number of oxygens you want to allow per species. Restricting the number 
of carbon atoms in each species may also be worthwhile to prevent very large molecules 
from being generated if many such species appear in your model.  Using any of these 
options requires some prior knowledge of the chemistry in your system. It is 
recommended that an initial model be generated without turning these options on. 
If many unlikely species show up in your model (or if your model has trouble 
converging and is generating many unlikely species on the edge), you can begin 
tuning these options to produce a better model. 

Adding key species into the initial condition file
--------------------------------------------------
Sometimes, chain branching reactions like dissociation of ROOH species do not make 
it to the core directly because if their fluxes are very small and the tolerance is 
not tight enough. In these cases, seeding the condition file with these species 
(with zero concentration) is helpful. By adding these species to the initial set 
of species in the condition file, the reactions involving those species will be 
automatically added to the core. (Putting these reactions in the seed mechanism 
has the same effect.)  Thus, if a species is known to be a part of your system 
and RMG is having trouble incorporating it within your model, it should be added 
to the condition file with 0.0 set as the concentration.

Starting with a single molecule when generating a model for a mixture
---------------------------------------------------------------------
For modeling the combustion of fuel mixtures, you may want to start with determining 
their composition and starting with a kinetic study of the dominant compound. It is 
possible to model the combustion of fuel mixtures but they are more challenging as 
well as harder to converge in RMG because RMG will automatically generate all cross 
reactions between the reacting species and intermedites. Starting with single species 
is always a good idea and is also useful when thinking about fuel mixtures. In order 
to build a better background in chemical kinetic model development and validation, 
please look at a recent paper from our group on butanol combustion available
`here <http://www.sciencedirect.com/science/article/pii/S0010218010001586>`_. 
This should give you some idea about how RMG can be put to use for the species 
of interest to you.
