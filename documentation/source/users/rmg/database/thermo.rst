.. _thermoDatabase:

***************
Thermo Database
***************

This section describes the general usage of RMG’s thermochemistry databases.
Thermochemical data in RMG is reported using three different quantities:

#. Standard heat capacity data :math:`C_p^o(T)` as a function of temperature :math:`T`
#. Standard enthalpy of formation at 298K :math:`\Delta_fH^{o}(298K)`
#. Standard entropy at 298K :math:`S^{o}(298K)`

A heat capacity model based on the Wilhoit equation is used for inter- and 
extrapolation of the heat capacity data as a function of temperature.

Libraries
=========

Species thermochemistry libraries
---------------------------------

The folder ``RMG-database/input/thermo/libraries/`` in RMG-database is the location to store
species thermochemistry libraries. Each particularly library is stored in a file
with the extension .py, e.g. 'DFT_QCI_thermo.py'.

An example of a species thermochemistry entry is shown here below::

	entry(
	    index = 1,
	    label = "H2",
	    molecule = 
	"""
	1 H 0 0 {2,S}
	2 H 0 0 {1,S}
	""",
	    thermo = ThermoData(
	        Tdata = ([300,400,500,600,800,1000,1500],'K'),
	        Cpdata = ([6.948,6.948,6.949,6.954,6.995,7.095,7.493],'cal/(mol*K)'),
	        H298 = (0,'kcal/mol'),
	        S298 = (31.095,'cal/(mol*K)'),
	    ),
	    shortDesc = u"""""",
	    longDesc = 
	u"""
	
	""",
	)
	
The text above describes the first entry in the library (index = 1), 
labeled 'H2', through the adjacency list representation. Heat capacity data ('Cpdata') is described
at 7 different temperatures, along with the standard enthalpy of formation at 298K ('H298'), and 
the standard entropy at 298K ('S298').

Groups
======

The folder ``RMG-database/input/thermo/groups/`` in RMG-database is the location to store
group contribution databases. Each particularly type of group contribution is stored in a file
with the extension .py, e.g. 'groups.py':

.. table::

    ======================================= ======================================================
    file		                            Type of group contribution
    ======================================= ======================================================
    gauche.py				                1,4-gauche non-nearest neighbor interactions (NNIs)
    group.py		              			group additive values (GAVs)
    int15.py			                	1,5-repulsion non-nearest neighbor interactions (NNIs)
    other.py				                other non-nearest neighbor interactions (NNIs)
    polycyclic.py		              		polycyclic ring corrections (RSCs)
    radical.py			                	hydrogen bond increments (HBIs)
    ring.py					                monocyclic ring corrections (RSCs)
    ======================================= ======================================================


Like many other entities in RMG, the database of each type of group contribution 
is organized in a hierarchical tree, and is defined at the bottom of the database file. E.g.::
	
	tree(
	"""
	L1: R
	    L2: C
	        L3: Cbf
	            L4: Cbf-CbCbCbf
	            L4: Cbf-CbCbfCbf
	            L4: Cbf-CbfCbfCbf
	        L3: Cb
	            L4: Cb-H
	            L4: Cb-Os
	            L4: Cb-Ss
	            L4: Cb-C
	                L5: Cb-Cs
	                L5: Cb-Cds
	                    L6: Cb-(Cds-Od)
	                    ...
  
More information on hierarchical tree structures in RMG can be found here:
:ref:`introDatabase`.

Group additive values (GAV)
---------------------------

An example of a GAV entry in group.py is shown here below::

	entry(
	    index = 3,
	    label = "Cbf-CbCbCbf",
	    group = 
	"""
	1 * Cbf 0 {2,B} {3,B} {4,B}
	2   Cb  0 {1,B}
	3   Cb  0 {1,B}
	4   Cbf 0 {1,B}
	""",
	    thermo = ThermoData(
	        Tdata = ([300,400,500,600,800,1000,1500],'K'),
	        Cpdata = ([3.01,3.68,4.2,4.61,5.2,5.7,6.2],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
	        H298 = (4.8,'kcal/mol','+|-',0.17),
	        S298 = (-5,'cal/(mol*K)','+|-',0.1),
	    ),
	    shortDesc = u"""Cbf-CbfCbCb STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
	    longDesc = 
	u"""
	
	""",
	)

The text above describes a GAV "Cbf-CbCbCbf", with the central atom denoted by the asterisk in 
the adjacency list representation. Uncertainty margins are added in the data, after the unit
specification. A short description 'shortDesc' specifies the origin of the data.


Ring Strain Corrections (RSC)
-----------------------------
RMG distinguishes between monocyclic and polycyclic ring correction databases. 

Monocyclic RSCs are used for molecules that contain one single ring.
An example of a  monocyclic RSC entry in ring.py is shown here below::

	entry(
	    index = 1,
	    label = "Cyclopropane",
	    group = 
	"""
	1 * Cs 0 {2,S} {3,S}
	2   Cs 0 {1,S} {3,S}
	3   Cs 0 {1,S} {2,S}
	""",
	    thermo = ThermoData(
	        Tdata = ([300,400,500,600,800,1000,1500],'K'),
	        Cpdata = ([-3.227,-2.849,-2.536,-2.35,-2.191,-2.111,-1.76],'cal/(mol*K)'),
	        H298 = (27.53,'kcal/mol'),
	        S298 = (32.0088,'cal/(mol*K)'),
	    ),
	    shortDesc = u"""Cyclopropane ring BENSON""",
	    longDesc = 
	u"""
	
	""",
	)

A molecule may have two or more fused rings that mutually interact. In that case, a polycyclic ring
strain correction may be more adequate. RMG identifies molecules with fused ring systems and subsequently
searches through polycyclic.py to identify an adequate RSC.
 
An example of a  polycyclic RSC entry in polycyclic.py is shown here below::

	entry(
	    index = 2,
	    label = "norbornane",
	    group = 
	"""
	1 * Cs 0 {3,S} {4,S} {7,S}
	2   Cs 0 {3,S} {5,S} {6,S}
	3   Cs 0 {1,S} {2,S}
	4   Cs 0 {1,S} {5,S}
	5   Cs 0 {2,S} {4,S}
	6   Cs 0 {2,S} {7,S}
	7   Cs 0 {1,S} {6,S}
	""",
	    thermo = ThermoData(
	        Tdata = ([300,400,500,600,800,1000,1500],'K'),
	        Cpdata = ([-4.5,-3.942,-3.291,-2.759,-2.08,-1.628,-0.898],'cal/(mol*K)'),
	        H298 = (16.14,'kcal/mol'),
	        S298 = (53.47,'cal/(mol*K)'),
	    ),
	    shortDesc = u"""""",
	    longDesc = 
	u"""
	
	""",
	)

Hydrogen Bond Increments (HBI)
------------------------------

An example of a HBI entry in radical.py is shown here below::

	entry(
	    index = 4,
	    label = "CH3",
	    group = 
	"""
	1 * C 1 {2,S} {3,S} {4,S}
	2   H 0 {1,S}
	3   H 0 {1,S}
	4   H 0 {1,S}
	""",
	    thermo = ThermoData(
	        Tdata = ([300,400,500,600,800,1000,1500],'K'),
	        Cpdata = ([0.71,0.34,-0.33,-1.07,-2.43,-3.54,-5.43],'cal/(mol*K)'),
	        H298 = (104.81,'kcal/mol','+|-',0.1),
	        S298 = (0.52,'cal/(mol*K)'),
	    ),
	    shortDesc = u"""Calculated in relation to methane from NIST values""",
	    longDesc = 
	u"""
	
	""",
	)

Non-nearest neighbor interactions
--------------------------------- 

The majority of the NNIs groups pertain to small enthalpy of formation corrections. Only a very limited
number include entropy or heat capacity corrections. The database other.py contains
cis-, ortho- and ketene-corrections.

An example of a NNI entry in gauche.py is shown here below::

	entry(
	    index = 11,
	    label = "Cs(Cs(CsCsR)Cs(CsCsR)RR)",
	    group = 
	"""
	1  * Cs                         0 {2,S} {3,S} {4,S} {5,S}
	2    Cs                         0 {1,S} {6,S} {7,S} {8,S}
	3    Cs                         0 {1,S} {9,S} {10,S} {11,S}
	4    {Cd,Cdd,Ct,Cb,Cbf,Os,CO,H} 0 {1,S}
	5    {Cd,Cdd,Ct,Cb,Cbf,Os,CO,H} 0 {1,S}
	6    Cs                         0 {2,S}
	7    Cs                         0 {2,S}
	8    {Cd,Cdd,Ct,Cb,Cbf,Os,CO,H} 0 {2,S}
	9    Cs                         0 {3,S}
	10   Cs                         0 {3,S}
	11   {Cd,Cdd,Ct,Cb,Cbf,Os,CO,H} 0 {3,S}
	""",
	    thermo = ThermoData(
	        Tdata = ([300,400,500,600,800,1000,1500],'K'),
	        Cpdata = ([0,0,0,0,0,0,0],'cal/(mol*K)'),
	        H298 = (0.8,'kcal/mol'),
	        S298 = (0,'cal/(mol*K)'),
	    ),
	    shortDesc = u"""""",
	    longDesc = 
	u"""
	
	""",
	)
