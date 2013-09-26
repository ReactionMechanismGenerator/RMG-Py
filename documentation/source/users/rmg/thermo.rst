.. _thermo:

**************************
Thermochemistry Estimation
**************************

Benson's Group Additivity approach, used by RMG, divides a molecule into functional 
groups, and the contribution of each functional group to the overall thermodynamics is 
included. RMG's kinetic rules use a similar strategy: the reaction rate coefficient
depends on the two (or more) functional groups involved in the reaction. The rate
coefficient rule database uses an innovative tree structure that allows users to 
include ever-more specific defintions of the functional groups and the kinetics 
involving them. In both cases, users can easily modify and expand these databases to 
suit their needs. In principle, the trees could be expanded until all possible molecules
and data for the reactions between all possible pairs of molecules have been included 
(of course, in practice, the database is much more limited). The kinetics database 
includes 19 reaction types ("families"), including hydrogen abstractions, cycloadditions,
disproportionation reactions, radical recombination reactions, intramolecular 
rearrangements, etc. 

Each reaction family includes data on many possible functional group pairings 
involved in the reaction, for a total of roughly 1000 rate rules. When a perfect 
match for a given reaction is not found in the database, RMG estimates the rate 
coefficient parameters based on "nearby" rate coefficient rules that are similar 
to the desired one. In addition to modifying the RMG databases, the user may also 
enter thermodynamic or kinetic information for specific molecules/reactions into 
Primary Thermo/Kinetic Libraries which supersede the estimated values. In this manner, 
precise data for well-studied reactions/species can be included. 

