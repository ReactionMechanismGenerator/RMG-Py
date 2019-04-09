.. _ratebasedmodelenlarger:

Rate-based Model Enlarging Algorithm
====================================


.. figure:: fluxAlgorithm.png

To construct a mechanism, the user must specify an initial set of species and
the initial conditions (temperature, pressure, species concentrations, etc.).
RMG reacts the initial species in all possible ways according to its known
reaction families, and it integrates the model in time. RMG tracks the rate
(flux) at which each new "edge" species is produced, and species (and the
reactions producing them) that are produced with significant fluxes are
incorporated into the model (the "core"). These new core species are reacted
with all other core species in the model to generate a new set of edge species
and reactions. The time-integration restarts, and the expanded list of edge
species is monitored for significant species to be included in the core. The
process continues until all significant species and reactions have been
included in the model. The definition of a "significant" rate can be specified by the user
by taking the following definition for a single species rate:

:math:`R_i = \frac{dC_i}{dt}`

and the following definition for the reaction system's characteristic rate, which is the sum of
all **core** species rates:

:math:`R_{char} = \sqrt{\sum\limits_{j} R_{j}^2}\quad    \quad  \textrm{species $j$ $\in$  core}`

When a :math:`\textrm{species $i$ $\in$ edge}`  exceeds a "significant" rate equal to :math:`\epsilon R_{char}`,
it is added to the core. The parameter :math:`\epsilon` is the user-specified
``toleranceMoveToCore`` that can be adjusted under the :ref:`model tolerances <modeltolerances>`
in the :ref:`RMG Input File <input>`.



For more information on rate-based model enlargement, please refer to the papers [Gao2016]_ or [Susnow1997]_. 

.. [Gao2016] \ C. W. Gao, J. W. Allen, W. H. Green, R. H. West, "Reaction Mechanism Generator: automatic construction of chemical kinetic mechanisms." *Computer Physics Communications* (2016).
.. [Susnow1997] \ R. G. Susnow, A. M. Dean, W. H. Green, P. K. Peczak, and L. J. Broadbelt. "Rate-Based Construction of Kinetic Models for Complex Systems." *J. Phys. Chem. A* **101**, p. 3731-3740 (1997).


.. _filterReactionsTheory:

Filtering Reactions within the Rate-based Algorithm
---------------------------------------------------

Filtering reactions in the react step in the flux-based algorithm attempts to speed up model generation by attacking the pain point.  RMG has trouble 
converging when generating models for large molecules because it searches for reactions on the order of :math:`(n_{reaction\: sites})^{{n_{species}}}`.  

The original algorithm performs in the following manner:

1. Reacts species together (slow) 
2. Determines which reactions are negligible (fast)

By filtering reactions we add a pre-filtering step before step 1 which prevents species from reacting together when the reactions are expected to be negligible
throughout the simulation.


.. figure:: fluxAlgorithmWithReactionFiltering.png

Here, ``unimolecularThreshold``, ``bimolecularThreshold``, and ``trimolecularThreshold`` are binary arrays storing flags for whether a species or a pair of species are above a reaction threshold.
For a unimolecular rate, this threshold is set to ``True`` if the unimolecular rate of :math:`\textrm{reaction $k$}` for a species A 

:math:`R_{unimolecular} = k_{threshold}C_A > \epsilon R_{char}` 

at any given time :math:`t` in the reaction system, where :math:`k_{threshold} = \frac{k_B T}{h}`

For a bimolecular reaction occuring between species A and B, this threshold is set to ``True`` if the bimolecular rate 

:math:`R_{bimolecular} = k_{threshold}C_A C_B > \epsilon R_{char}` 

where :math:`k_{threshold} = filterThreshold`. ``filterThreshold`` is adjusted for each reaction family in RMG, ensuring that very fast reactions are not accidentally filtered out. The fits shown in the table below are obtained by extracting the highest reaction rates from all training reactions of one family, for a discrete set of temperatures, and fitting these highest reaction rates using an Arrhenius fit. The used ipython notebook is RMG-Py/ipython/kinetics_training_analysis.ipynb.

Similarly, for a trimolecular reaction, the following expression is used:

:math:`R_{trimolecular} = k_{threshold}C_A C_B C_C > \epsilon R_{char}`

where :math:`k_{threshold} = 10^{-3} \cdot filterThreshold \frac{m^6}{mol^2\cdot s}`. Based on extending Smoluchowski theory to multiple molecules, the diffusion limit rate constant for trimolecular reactions (in :math:`\frac{m^6}{mol^2\cdot s}`) is approximately three orders of magnitude smaller than the rate constant for bimolecular reactions (in :math:`\frac{m^3}{mol\cdot s}`). It is assumed here that Smoluchowski theory gives a sufficient approximation to collision theory in the gas phase.

When the liquid-phase reactor is used, the diffusion limits are calculated using the Stokes-Einstein equation instead.
For bimolecular reactions, this results in

:math:`k_{threshold}[m^3/mol/s] = 22.2\frac{T[K]}{\mu[Pa\cdot s]}`

and for trimolecular

:math:`k_{threshold}[m^6/mol^2/s] = 0.11\frac{T[K]}{\mu[Pa\cdot s]}`

where :math:`\mu` is the solvent viscosity. The coefficients in the above equations were obtained by using a
representative value of the molecular radius of 2 Angstrom. More details on the calculation of diffusion limits in the
liquid phase can be found in :ref:`the description of liquid-phase systems <liquids>` under
:ref:`diffusion-limited kinetics <diffusionLimited>`.

Three additional binary arrays ``unimolecularReact``, ``bimolecularReact``, and ``trimolecularReact`` store flags for
when the ``unimolecularThreshold``, ``bimolecularThreshold``, or ``trimolecularThreshold`` flag
shifts from ``False`` to ``True``.  RMG reacts species when the flag is set to ``True``.

Table: filterThreshold for each reaction family.

==================================================== ================================================= ====================================================
                     Reaction family                     		Unimolecular                                   	Bimolecular 
==================================================== ================================================= ====================================================
Intra_R_Add_Exocyclic                                    Arrhenius(A=(2.1548e-15,'s^-1'),n=8.65061, 			:math:`10^8`
							  Ea=(-122.419,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA = 20.5283, dn = +|- 0.383077, 
							  dEa = +|- 2.61431 kJ/mol""")
Cyclopentadiene_scission                                 Arrhenius(A=(2.51056e+17,'s^-1'), 				:math:`10^8`
							  n=-1.48728, Ea=(95.6898,'kJ/mol'),  
							  T0=(1,'K'), Tmin=(300,'K'),  
							  Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA = / 1, dn = +|- 2.11795e-15, 
							  dEa = +|- 1.4454e-14 kJ/mol""")   
2+2_cycloaddition_CO                                     Arrhenius(A=(4.19097e+11,'s^-1'),n=0.542031,	Arrhenius(A=(2.319e-07,'m^3/(mol*s)'), n=3.416,
							  Ea=(48.0397,'kJ/mol'), T0=(1,'K'), 		 Ea=(322.616,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 		 comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 2.67845e-15,		 dA =  / 1, dn = +|- 5.81451e-15, 
							  dEa = +|- 1.82791e-14 kJ/mol""")		 dEa = +|- 3.96812e-14 kJ/mol""")					 
R_Addition_CSm                                           Arrhenius(A=(3.06643e+07,'s^-1'),n=4.05506, 	Arrhenius(A=(1.2e+07,'m^3/(mol*s)'), n=2.11, 
							  Ea=(364.729,'kJ/mol'), T0=(1,'K'),		 Ea=(10.2926,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 		 comment="""Fitted to 30 data points; 
							  dA =  / 3.71715, dn = +|- 0.166445,		 dA =  / 1, dn = +|- 1.66533e-15, 
							  dEa = +|- 1.1359 kJ/mol""")			 dEa = +|- 1.13651e-14 kJ/mol""")
Disproportionation                                          		:math:`10^8`			Arrhenius(A=(3.71358e-08,'m^3/(mol*s)'), n=4.90833, 
													 Ea=(-21.5849,'kJ/mol'), T0=(1,'K'), 
													 Tmin=(300,'K'), Tmax=(2500,'K'), 
													 comment="""Fitted to 30 data points; 
													 dA =  / 8.35058, dn = +|- 0.26905, 
													 dEa = +|- 1.83613 kJ/mol""")
1,2-Birad_to_alkene    					 Arrhenius(A=(1e+10,'s^-1'),n=-6.64137e-15, 			:math:`10^8`
							  Ea=(6.11426e-14,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 1.69925e-15, 
							  dEa = +|- 1.15965e-14 kJ/mol""")                                     
Intra_R_Add_Exo_scission 				 Arrhenius(A=(7.809e+07,'s^-1'), n=1.057, 			:math:`10^8`
							  Ea=(63.0152,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 1.62696e-15, 
							  dEa = +|- 1.11032e-14 kJ/mol""")                                   
H2_Loss 						 Arrhenius(A=(4.82588e+09,'s^-1'),n=0.803687,	Arrhenius(A=(477137,'m^3/(mol*s)'), n=2.9449, 
							  Ea=(72.0667,'kJ/mol'), T0=(1,'K'),		 Ea=(-96.2703,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 1.15661, dn = +|- 0.0184441,		 dA =  / 1, dn = +|- 3.48379e-15, 
							  dEa = +|- 0.125872 kJ/mol""")			 dEa = +|- 2.37752e-14 kJ/mol""") 
1,3_Insertion_ROR      					 Arrhenius(A=(4.49138e+06,'s^-1'),n=3.19054, 	Arrhenius(A=(4.86e-07,'m^3/(mol*s)'), n=3.55, 								  Ea=(123.447,'kJ/mol'), T0=(1,'K'),		 Ea=(101.797,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 		 comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 3.08198e-15,		 dA =  / 1, dn = +|- 2.24939e-15, 
							  dEa = +|- 2.1033e-14 kJ/mol""")		 dEa = +|- 1.5351e-14 kJ/mol""")
Baeyer-Villiger_step1_cat   						:math:`10^8`			Arrhenius(A=(3.46333e+06,'m^3/(mol*s)'), 
													 n=-0.929161, Ea=(42.315,'kJ/mol'), T0=(1,'K'), 													 Tmin=(300,'K'), Tmax=(2500,'K'), 
													 comment="""Fitted to 30 data points; 
													 dA =  / 1, dn = +|- 9.48759e-16, 
													 dEa = +|- 6.47482e-15 kJ/mol""")
Intra_RH_Add_Endocyclic    						:math:`10^8`					:math:`10^8`                                 
Baeyer-Villiger_step2_cat                                   		:math:`10^8`			Arrhenius(A=(4.74858e-09,'m^3/(mol*s)'), 
													 n=4.24247, Ea=(83.3223,'kJ/mol'), T0=(1,'K'), 														 Tmin=(300,'K'), Tmax=(2500,'K'), 
													 comment="""Fitted to 30 data points; 
													 dA =  / 1, dn = +|- 2.32974e-15, 
													 dEa = +|- 1.58993e-14 kJ/mol""")
Korcek_step2       							:math:`10^8`					:math:`10^8`                                         
Singlet_Val6_to_triplet 				 Arrhenius(A=(2.69922e+09,'s^-1'),n=2.23108, 			:math:`10^8`
							  Ea=(-23.2466,'kJ/mol'), T0=(1,'K'), 
 							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 5.18462, dn = +|- 0.208626, 
							  dEa = +|- 1.42377 kJ/mol""")          
Intra_Retro_Diels_alder_bicyclic 					:math:`10^8`					:math:`10^8`			                           
R_Addition_MultipleBond                			 Arrhenius(A=(8.36742e-37,'s^-1'),n=15.6692, 	Arrhenius(A=(3.92165e-11,'m^3/(mol*s)'), 
							  Ea=(-79.9819,'kJ/mol'), T0=(1,'K'),		 n=5.41917, Ea=(-92.6718,'kJ/mol'), T0=(1,'K'), 							  Tmin=(300,'K'), Tmax=(2500,'K'), 		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 37.755, dn = +|- 0.46032,	 	 dA =  / 17.0778, dn = +|- 0.359748, 
							  dEa = +|- 3.14146 kJ/mol""")			 dEa = +|- 2.4551 kJ/mol""")
Concerted_Intra_Diels_alder_monocyclic_1,2_shiftH        Arrhenius(A=(3.60572e+16,'s^-1'),				:math:`10^8`
							  n=-0.490221, 			
							  Ea=(118.723,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 3.56237e-15, 
							  dEa = +|- 2.43114e-14 kJ/mol""")   
Cyclic_Thioether_Formation 						:math:`10^8`					:math:`10^8`                                 
Intra_R_Add_Endocyclic 					 Arrhenius(A=(31.1508,'s^-1'), n=4.20415, 			:math:`10^8`
							  Ea=(-67.1521,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 4.40734, dn = +|- 0.188036, 
							  dEa = +|- 1.28325 kJ/mol""")                                     
1,3_Insertion_CO2					 Arrhenius(A=(22135.6,'s^-1'), n=3.10576, 	Arrhenius(A=(1.33741,'m^3/(mol*s)'), n=2.15746, 							  Ea=(287.824,'kJ/mol'), T0=(1,'K'),		 Ea=(304.387,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 4.16477, dn = +|- 0.180859, 		 dA =  / 4.07671, dn = +|- 0.17815, 
							  dEa = +|- 1.23428 kJ/mol""")			 dEa = +|- 1.21579 kJ/mol""")
1+2_Cycloaddition          				 Arrhenius(A=(1.25061e+25,'s^-1'),n=-2.6694,   	Arrhenius(A=(8.19836e+06,'m^3/(mol*s)'), 
							  Ea=(362.294,'kJ/mol'), T0=(1,'K'),		 n=0.145041, Ea=(-3.96815,'kJ/mol'), T0=(1,'K'), 							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 1.10756, dn = +|- 0.0129508,		 dA =  / 1.50944, dn = +|- 0.0521967, 
							  dEa = +|- 0.088383 kJ/mol""")			 dEa = +|- 0.356217 kJ/mol""")			                                   
Bimolec_Hydroperoxide_Decomposition 					:math:`10^8`			Arrhenius(A=(127965,'m^3/(mol*s)'), n=0.933342, 													 Ea=(111.519,'kJ/mol'), T0=(1,'K'), 
													 Tmin=(300,'K'), Tmax=(2500,'K'), 
													 comment="""Fitted to 30 data points; 
													 dA =  / 6.41962, dn = +|- 0.235713, 
													 dEa = +|- 1.60862 kJ/mol""")                        
Intra_R_Add_ExoTetCyclic 						:math:`10^8`					:math:`10^8`                                   
Peroxyl_Termination   							:math:`10^8`			Arrhenius(A=(120000,'m^3/(mol*s)'), 
													 n=-1.31735e-14, Ea=(-4.184,'kJ/mol'), 
													 T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), 														 comment="""Fitted to 30 data points; 
													 dA =  / 1, dn = +|- 1.94625e-15, 
													 dEa = +|- 1.32822e-14 kJ/mol""")                                      
CO_Disproportionation							:math:`10^8`			Arrhenius(A=(1.2e+08,'m^3/(mol*s)'), 
													 n=-1.20033e-14, Ea=(1.11491e-13,'kJ/mol'), 														 T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), 														 comment="""Fitted to 30 data points; 
													 dA =  / 1, dn = +|- 3.09852e-15, 
													 dEa = +|- 2.11459e-14 kJ/mol""")                                       
Intra_Disproportionation				 Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, 			:math:`10^8`
							  Ea=(22.8614,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 3.78153e-15, 
							  dEa = +|- 2.58071e-14 kJ/mol""")                                    
SubstitutionS  								:math:`10^8`			Arrhenius(A=(0.000160472,'m^3/(mol*s)'), 
													 n=3.8024, Ea=(-10.9764,'kJ/mol'), T0=(1,'K'), 														 Tmin=(300,'K'), Tmax=(2500,'K'), 
													 comment="""Fitted to 30 data points; 
													 dA =  / 15.0279, dn = +|- 0.343538, 
													 dEa = +|- 2.34448 kJ/mol""")			                                      
Korcek_step1  								:math:`10^8`					:math:`10^8`                                              
intra_substitutionS_cyclization 			 Arrhenius(A=(7.15311,'s^-1'), n=3.63158, 	Arrhenius(A=(17580.9,'m^3/(mol*s)'), n=1.98548, 							  Ea=(32.7853,'kJ/mol'), T0=(1,'K'),		 Ea=(-106.532,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 51.476, dn = +|- 0.499619,		 dA =  / 1, dn = +|- 4.80061e-15, 
							  dEa = +|- 3.40965 kJ/mol""")			 dEa = +|- 3.27618e-14 kJ/mol""")		                              
Korcek_step1_cat 							:math:`10^8`					:math:`10^8`                                          
1,4_Linear_birad_scission                                   		:math:`10^8`					:math:`10^8`
1,2_Insertion_carbene      				Arrhenius(A=(9.09546e+17,'s^-1'), 		Arrhenius(A=(51795.3,'m^3/(mol*s)'), n=0.681821, 							  n=-0.747035, Ea=(447.041,'kJ/mol'),		 Ea=(-9.92113,'kJ/mol'), T0=(1,'K'), 
							  T0=(1,'K'), Tmin=(300,'K'),			 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  Tmax=(2500,'K'), 				 comment="""Fitted to 30 data points; 
							  comment="""Fitted to 30 data points;		 dA =  / 2.01108, dn = +|- 0.0885713, 
							  dA =  / 1.42527, dn = +|- 0.0449229,		 dEa = +|- 0.604456 kJ/mol""")
							  dEa = +|- 0.306577 kJ/mol""")				                                   
H_Abstraction 								:math:`10^8`			Arrhenius(A=(1.46107e-12,'m^3/(mol*s)'), 
													 n=4.73222, Ea=(-278.229,'kJ/mol'), T0=(1,'K'), 													 Tmin=(300,'K'), Tmax=(2500,'K'), 
													 comment="""Fitted to 30 data points; 
													 dA =  / 16576.3, dn = +|- 1.23167, 
													 dEa = +|- 8.40556 kJ/mol""")                                              
Intra_5_membered_conjugated_C=C_C=C_addition   		Arrhenius(A=(9.86304e+10,'s^-1'),n=0.836047, 			:math:`10^8`
							  Ea=(79.2187,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 1.74493e-15, 
							  dEa = +|- 1.19083e-14 kJ/mol""")             
Intra_ene_reaction   					Arrhenius(A=(27.788,'s^-1'), n=3.56981, 			:math:`10^8`
							  Ea=(17.5282,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 2.95262, dn = +|- 0.137254, 
							  dEa = +|- 0.936691 kJ/mol""")                                       
intra_H_migration  					Arrhenius(A=(32804.1,'s^-1'), n=2.27586, 			:math:`10^8`
							  Ea=(-3.34049,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 5.22479, dn = +|- 0.209605, 
							  dEa = +|- 1.43045 kJ/mol""")                                         
Baeyer-Villiger_step2  					Arrhenius(A=(2.8594e+09,'s^-1'), n=1.09585, 	Arrhenius(A=(2.03275e-17,'m^3/(mol*s)'), n=6.6413, 							  Ea=(90.8931,'kJ/mol'), T0=(1,'K'),		 Ea=(205.734,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 2.42346, dn = +|- 0.112217,		 dA =  / 1, dn = +|- 6.88308e-15, 
							  dEa = +|- 0.765828 kJ/mol""")			 dEa = +|- 4.69737e-14 kJ/mol""")	1,2_Insertion_CO 					Arrhenius(A=(3.33623e+08,'s^-1'), n=1.58566,    Arrhenius(A=(1.27e-07,'m^3/(mol*s)'), n=3.7, 								  Ea=(274.448,'kJ/mol'), T0=(1,'K'),		 Ea=(223.258,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 4.34887e-15,		 dA =  / 1, dn = +|- 2.66187e-15, 
							  dEa = +|- 2.96789e-14 kJ/mol""")		 dEa = +|- 1.81659e-14 kJ/mol""")
Substitution_O   							:math:`10^8`			Arrhenius(A=(1.34785e-07,'m^3/(mol*s)'), n=4.7731, 													 Ea=(-85.8887,'kJ/mol'), T0=(1,'K'), 
													 Tmin=(300,'K'), Tmax=(2500,'K'), 
													 comment="""Fitted to 30 data points; 
													 dA =  / 1, dn = +|- 2.54359e-15, 
													 dEa = +|- 1.73587e-14 kJ/mol""")                                           
Intra_RH_Add_Exocyclic  						:math:`10^8`					:math:`10^8`                                    
Cyclic_Ether_Formation 					Arrhenius(A=(0.0211122,'s^-1'), n=4.09341, 	Arrhenius(A=(8.18963e-31,'m^3/(mol*s)'), n=11.428, 							  Ea=(4.27861,'kJ/mol'), T0=(1,'K'),		 Ea=(-13.2083,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 33.1388, dn = +|- 0.443787,		 dA =  / 140.706, dn = +|- 0.627094, 
							  dEa = +|- 3.02863 kJ/mol""")			 dEa = +|- 4.27961 kJ/mol""")			                                    
1,2_shiftC    						Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, 			:math:`10^8`
							  Ea=(94.4747,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 1.63603e-15, 
							  dEa = +|- 1.11651e-14 kJ/mol""")                                              
lone_electron_pair_bond     						:math:`10^8`					:math:`10^8`                               
HO2_Elimination_from_PeroxyRadical  			Arrhenius(A=(21583.8,'s^-1'), n=2.90923, 	Arrhenius(A=(7.37936e-08,'m^3/(mol*s)'), 
							  Ea=(106.401,'kJ/mol'), T0=(1,'K'),		 n=3.90078, Ea=(20.1194,'kJ/mol'), T0=(1,'K'), 								  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 		
							  dA =  / 6.91976, dn = +|- 0.245223, 		 dA =  / 3.75901, dn = +|- 0.167864, 
							  dEa = +|- 1.67353 kJ/mol""") 			 dEa = +|- 1.14559 kJ/mol""")			 
Birad_recombination    					Arrhenius(A=(2.18e+16,'s^-1'), n=-1.0059e-14, 			:math:`10^8`
							  Ea=(2.9288,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 5.97248e-15, 
							  dEa = +|- 4.07593e-14 kJ/mol""")                                    
Diels_alder_addition    				Arrhenius(A=(1.17242e-21,'s^-1'), n=11.0661, 	Arrhenius(A=(8.35312e-24,'m^3/(mol*s)'), 
							  Ea=(17.1225,'kJ/mol'), T0=(1,'K'),		 n=9.17707, Ea=(-1.16914,'kJ/mol'), T0=(1,'K'), 							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 4.46687e-15,		 dA =  / 29.3773, dn = +|- 0.428514, 
							  dEa = +|- 3.04842e-14 kJ/mol""") 		 dEa = +|- 2.9244 kJ/mol""")			                                     
R_Addition_COm  					Arrhenius(A=(1.5378e+14,'s^-1'), n=0.264564, 	Arrhenius(A=(0.0202143,'m^3/(mol*s)'), 
							  Ea=(20.2142,'kJ/mol'), T0=(1,'K'),		 n=2.45362, Ea=(2.95866,'kJ/mol'), T0=(1,'K'), 								  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 3.09632e-15,		 dA =  / 3.00278, dn = +|- 0.13939, 
							  dEa = +|- 2.11309e-14 kJ/mol""")		 dEa = +|- 0.951266 kJ/mol""")			                                              
intra_substitutionCS_cyclization   					:math:`10^8`					:math:`10^8`                         
2+2_cycloaddition_CS   							:math:`10^8`					:math:`10^8`                                     
1,2_shiftS     								:math:`10^8`					:math:`10^8`                                             
intra_OH_migration  					Arrhenius(A=(181.219,'s^-1'), n=2.3668, 			:math:`10^8`
							  Ea=(50.9862,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 5.67183, dn = +|- 0.220012, 
							  dEa = +|- 1.50148 kJ/mol""")                                        
Birad_R_Recombination      				Arrhenius(A=(9.40883e+22,'s^-1'), n=-0.830214, 	Arrhenius(A=(1.29008e+06,'m^3/(mol*s)'), 
							  Ea=(159.721,'kJ/mol'), T0=(1,'K'),		 n=0.806257, Ea=(-5.31319,'kJ/mol'), T0=(1,'K'), 							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 5.2588e-15,		 dA =  / 2.22228, dn = +|- 0.101231, 
							  dEa = +|- 3.58887e-14 kJ/mol""")		 dEa = +|- 0.690853 kJ/mol""")		                                   
Singlet_Carbene_Intra_Disproportionation   		Arrhenius(A=(1.454e+12,'s^-1'), n=0.178, 			:math:`10^8`
							  Ea=(0.85772,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 3.373e-15, 
							  dEa = +|- 2.30191e-14 kJ/mol""")                 
6_membered_central_C-C_shift                       	Arrhenius(A=(3.53521e+20,'s^-1'), n=-2.14941, 			:math:`10^8`
							  Ea=(84.4898,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 2.0056e-15, 
							  dEa = +|- 1.36872e-14 kJ/mol""")         
intra_substitutionCS_isomerization  			Arrhenius(A=(607.614,'s^-1'), n=2.9594, 			:math:`10^8`
							  Ea=(180.721,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 6.44969, dn = +|- 0.236305, 
							  dEa = +|- 1.61267 kJ/mol""")                        
2+2_cycloaddition_CCO    						:math:`10^8`					:math:`10^8`                                   
intra_substitutionS_isomerization  			Arrhenius(A=(2228.54,'s^-1'), n=2.59523, 			:math:`10^8`
							  Ea=(79.3728,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 9.24147, dn = +|- 0.281901, 
							  dEa = +|- 1.92383 kJ/mol""")                         
2+2_cycloaddition_Cd      				Arrhenius(A=(3.47729e+15,'s^-1'),n=-0.0133205,	Arrhenius(A=(4.66,'m^3/(mol*s)'), n=1.65, 								  Ea=(253.025,'kJ/mol'), T0=(1,'K'),		 Ea=(226.564,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 3.42822e-15,		 dA =  / 1, dn = +|- 4.19081e-15, 
							  dEa = +|- 2.33959e-14 kJ/mol""")		 dEa = +|- 2.86002e-14 kJ/mol""") 		                                    
Intra_2+2_cycloaddition_Cd      			Arrhenius(A=(0.0435533,'s^-1'), n=4.08745, 			:math:`10^8`
							  Ea=(103.277,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 15.3826, dn = +|- 0.346495, 
							  dEa = +|- 2.36466 kJ/mol""")                            
Intra_Diels_alder_monocyclic 				Arrhenius(A=(1.4544e+12,'s^-1'), n=0.301801, 			:math:`10^8`
							  Ea=(-1.2548,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 1.47428e-15, 
							  dEa = +|- 1.00612e-14 kJ/mol""")                               
R_Recombination    					Arrhenius(A=(6.85921e-07,'s^-1'), n=5.65681, 	Arrhenius(A=(5.6874e+10,'m^3/(mol*s)'), 								  Ea=(-64.6036,'kJ/mol'), T0=(1,'K'),		 n=-0.0206709, Ea=(-191.106,'kJ/mol'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'),		 T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), 								  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 21.0499, dn = +|- 0.386257, 		 dA =  / 1, dn = +|- 3.38094e-15, 
							  dEa = +|- 2.63602 kJ/mol""") 			 dEa = +|- 2.30733e-14 kJ/mol""")		                                          
1,3_Insertion_RSR     					Arrhenius(A=(6.13062e+10,'s^-1'), n=2.21279, 	Arrhenius(A=(3.80823e-13,'m^3/(mol*s)'),
							  Ea=(128.085,'kJ/mol'), T0=(1,'K'),		 n=5.9792, Ea=(141.605,'kJ/mol'), T0=(1,'K'), 								  Tmin=(300,'K'), Tmax=(2500,'K'),		 Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points;		 comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 2.304e-15,		 dA =  / 5.14419, dn = +|- 0.207634, 
							  dEa = +|- 1.57237e-14 kJ/mol""")		 dEa = +|- 1.417 kJ/mol""")			                                        
1,4_Cyclic_birad_scission     				Arrhenius(A=(1.67986e+13,'s^-1'), n=0.420292, 			:math:`10^8`
							  Ea=(21.9887,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 2.88394e-15, 
							  dEa = +|- 1.96815e-14 kJ/mol""")                              
intra_NO2_ONO_conversion      				Arrhenius(A=(722095,'s^-1'), n=2.22801, 			:math:`10^8`
							  Ea=(238.399,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
 							  comment="""Fitted to 30 data points; 
							  dA =  / 1.94083, dn = +|- 0.0840635, 
							  dEa = +|- 0.573692 kJ/mol""")                              
ketoenol                                     		Arrhenius(A=(104,'s^-1'), n=3.21, 				:math:`10^8`
							  Ea=(82.0482,'kJ/mol'), T0=(1,'K'), 
							  Tmin=(300,'K'), Tmax=(2500,'K'), 
							  comment="""Fitted to 30 data points; 
							  dA =  / 1, dn = +|- 7.65184e-15, 
							  dEa = +|- 5.22201e-14 kJ/mol""")               
Peroxyl_Disproportionation                                		:math:`10^8`			Arrhenius(A=(1.1e+06,'m^3/(mol*s)'), 
													 n=-1.05506e-14, Ea=(-4.184,'kJ/mol'), 
													 T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), 														 comment="""Fitted to 30 data points; 
													 dA =  / 1, dn = +|- 1.19083e-15, 
													 dEa = +|- 8.12686e-15 kJ/mol""")  
==================================================== ================================================= ====================================================

