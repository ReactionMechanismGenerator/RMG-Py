.. _features:

********************
Overview of Features
********************

Thermodynamics estimation using group additivity.
	Group additivity based on Benson's groups provide fast and reliable thermochemistry estimates. A standalone utility for estimating heat of formation, entropy, and heat capacity is also included.

Rate-based model enlargement 
 	Reactions are added to the model based on their rate, fastest first.

Rate-based termination. 
	The model enlargement stops when all excluded reactions are slower than a given threshold.
	This provides a controllable error bound on the kinetic model that is generated.

Extensible libraries
	Ability to include reaction models on top of the provided reaction families.

Pressure-dependent reaction networks. 
	Dissociation, combination, and isomerization reactions have the potential to have rate coefficients that are dependent on both temperature and pressure, and RMG is able to estimate both for networks of arbitrary complexity with a bounded error.
	
Simultaneous mechanism generation for several conditions.
	Concurrent generation of a reaction mechanism over multiple temperature and pressure conditions. 
	Mechanisms generated this way are valid over a range of reaction conditions.

Dynamic simulation to a target conversion or time.
	Often the desired simulation time is not known *a priori*, so a target conversion is preferred.

Transport properties estimation using group additivity (in development)
	The Lennard-Jones sigma and epsilon parameters are estimated using empirical correlations (based on a species' critical properties and acentric factor).
	The critical properties are estimated using a group-additivity approach; the acentric factor is also estimated using empirical correlations.
	A standalone application for estimating these parameters is provided, and the output is stored in CHEMKIN-readable format.
