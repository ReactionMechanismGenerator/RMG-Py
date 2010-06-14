# encoding: utf-8
header = """
H_Abstraction

X_H + Y_rad_birad -> X_rad + Y_H


(1) BREAK_BOND		{*1,S,*2}
(2) FORM_BOND		{*2,S,*3}
(3) GAIN_RADICAL	{*1,1}
(4) LOSE_RADICAL 	{*3,1}


Generated on 7th April 2010 at 17:08
Generated on 22nd June 2010 at 11:22
Generated on 22nd June 2010 at 12:28
Generated on 22nd June 2010 at 12:58
"""

reaction_family_name = "H_Abstraction"

# determines permitted units for rate expression:
reaction_order = 2

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f01: HAbstraction reaction
// original from rate library.txt, CDW 03/08/01
// SR and JS rename some nodes according to the tree and dictionary correction, Nov., 20, 2002
// JS, comment out 10th, change XH in 23 from C/Cd/H3 to C/H3/Cd
// JS, remove CO_birad to form a new family later: CO + RH -> HCO + R.  Aug, 26, 2003

// JS, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// Catherina Wijaya thesis, pg 151 -154.

//f01_intermolecular_HA
//No.	XH			Y_rad				Temp.		A			n		a		E0		DA		Dn		Da		DE0		Rank	Comments
//142.	C/H3/Cs		O_pri_rad			300-1500	2.92E+06	1.80	0		0.278	0		0		0		0		5		Curran et al. [8] Rate expressions for H atom abstraction from fuels. (changed to per H)	
//208.	Cd_pri		CO_birad			300-2500	3.78E+13	0		0		90.62	*5.0	0		0		0		4		Tsang [89] literature review.
//213.	Cd/H/NonDeC	H_rad				300-2500	2.37E+00	0		0		7.31	*5.0	0		0		0		4		Tsang [93] literature review.
//214.	Cd/H/NonDeC	C_methyl			300-2500	3.8E-01		0		0		5.98	*100.0	0		0		0		4		Tsang [93] literature review.
//215.	Cd/H/NonDeC	Cd_pri_rad			300-2500	3.8E-01		0		0		4.99	*10.0	0		0		0		4		Tsang [93] literature review.
//217.	Cd/H/NonDeC	O_pri_rad			300-2500	1.11E+06	2.00	0		1.45	0		0		0		0		4		Tsang [93] literature review.
//219.	Ct_H		CO_birad			300-2500	2.41E+14	0		0		107		*10.0	0		0		0		4		Tsang [89] literature review.
//245.	CO/H/NonDe	O_pri_rad			295-600		2.0E+06		1.80	0		1.30	0		0		0		0		3		Taylor et al. [127] Transition state theory.
//257.	O_pri		O_pri_rad			200-700		6.45E-06	5.01	0		0.61	0		0		0		0		3		Masgrau et al. [141] Transition state theory w/tunneling correction.
//263.	O/H/NonDeC	C_rad/H/NonDeC		300-2500	1.45E+01	3.10	0		10.33	*5.0	0		0		0		4		Tsang [90] literature review.
//264.	O/H/NonDeC	C_rad/Cs3			300-2500	1.51E+03	1.80	0		9.36	*10.0	0		0		0		4		Tsang [90] literature review.
// (Reverse) Rates for nButanol+HO2=H2O2+radicals
// (Reverse) Rates for sButanol+HO2=H2O2+radicals
// (Reverse) Rates for tButanol+HO2=H2O2+radicals

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
X_H_or_Xrad_H
Union {X_H, Xrad_H}
""",
  group2 = 
"""
Y_rad_birad
Union {Y_2centeradjbirad, Y_1centerbirad, Y_rad}
""",
  kf = Arrhenius(A=(1E+05,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "0",
  short_comment = "Default",
  long_comment = 
"""
If a biradical CH2JJ can abstract from RCH4 to make RCH3J and CH3J 
then a Y_rad CH3J should be able to abstract from RCH3J which means X_H needs 
to include Xrad_H. I.e. you can abstract from a radical. To make this possible
a head node has been created X_H_or_Xrad_H which is a union of X_H and Xrad_H.
The kinetics for it have just been copied from X_H and are only defined for 
abstraction by Y_rad_birad. I.e. the top level very approximate guess.

Do better kinetics for this exist? Do we in fact use the reverse kinetics anyway?
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 2
rate(
  group1 = 
"""
X_H
1 *1 R 0 {2,S}
2 *2 H 0 {1,S}
""",
  group2 = 
"""
Y_rad_birad
Union {Y_2centeradjbirad, Y_1centerbirad, Y_rad}
""",
  kf = Arrhenius(A=(1E+05,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "1",
  short_comment = "Default",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
X_H
1 *1 R 0 {2,S}
2 *2 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.4E+08,A_UNITS,"+-",0.0),
                 n=(1.5,None,"+-",0.0),
                 alpha=(0.65,None,"+-",0.0),
                 E0=(9.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "2",
  short_comment = "Dean, A. M. [118]",
  long_comment = 
"""
[118] Dean, A.M. Development and application of Detailed Kinetic Mechanisms for Free Radical Systems.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
X_H
1 *1 R 0 {2,S}
2 *2 H 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(1.7E+08,A_UNITS,"+-",0.0),
                 n=(1.5,None,"+-",0.0),
                 alpha=(0.75,None,"+-",0.0),
                 E0=(6.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "3",
  short_comment = "Dean, A. M. [118]",
  long_comment = 
"""
[118] Dean, A.M. Development and application of Detailed Kinetic Mechanisms for Free Radical Systems.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
X_H
1 *1 R 0 {2,S}
2 *2 H 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.2E+06,A_UNITS,"+-",0.0),
                 n=(2.0,None,"+-",0.0),
                 alpha=(0.50,None,"+-",0.0),
                 E0=(10.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "4",
  short_comment = "Dean, A. M. [118]",
  long_comment = 
"""
[118] Dean, A.M. Development and application of Detailed Kinetic Mechanisms for Free Radical Systems.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
X_H
1 *1 R 0 {2,S}
2 *2 H 0 {1,S}
""",
  group2 = 
"""
O_sec_rad
1 *3 O 1 {2,S}
2 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.4E+04,A_UNITS,"+-",0.0),
                 n=(2.69,None,"+-",0.0),
                 alpha=(0.60,None,"+-",0.0),
                 E0=(11.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "5",
  short_comment = "Dean, A. M. [118]",
  long_comment = 
"""
[118] Dean, A.M. Development and application of Detailed Kinetic Mechanisms for Free Radical Systems.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
X_H
1 *1 R 0 {2,S}
2 *2 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.1E+05,A_UNITS,"+-",0.0),
                 n=(1.87,None,"+-",0.0),
                 alpha=(0.65,None,"+-",0.0),
                 E0=(13.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "6",
  short_comment = "Dean, A. M. [118]",
  long_comment = 
"""
[118] Dean, A.M. Development and application of Detailed Kinetic Mechanisms for Free Radical Systems.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.30E-01,A_UNITS,"+-",0.0),
                 n=(3.71,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.85,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "7",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 9
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.26E+00,A_UNITS,"+-",0.0),
                 n=(3.55,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.31,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "8",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 10
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.85E+00,A_UNITS,"+-",0.0),
                 n=(3.62,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "9",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 11
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.58E+01,A_UNITS,"+-",0.0),
                 n=(3.01,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.34,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "10",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 12
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.52E+01,A_UNITS,"+-",0.0),
                 n=(3.19,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.31,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "11",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 13
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.71E+01,A_UNITS,"+-",0.0),
                 n=(3.23,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.27,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "12",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 14
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.22E+03,A_UNITS,"+-",0.0),
                 n=(2.51,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.06,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "13",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 15
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.54E+03,A_UNITS,"+-",0.0),
                 n=(2.66,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "14",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 16
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(6.59E+02,A_UNITS,"+-",0.0),
                 n=(2.71,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.92,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "15",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 17
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.74E+05,A_UNITS,"+-",0.0),
                 n=(1.83,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.94,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "16",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 18
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.45E+06,A_UNITS,"+-",0.0),
                 n=(1.77,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.53,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "17",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 19
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.78E+05,A_UNITS,"+-",0.0),
                 n=(1.9,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.05,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "18",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 20
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.01E+04,A_UNITS,"+-",0.0),
                 n=(2.47,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.96,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "19",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 21
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(4.83E+08,A_UNITS,"+-",0.0),
                 n=(1.54,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.98,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "20",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 22
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.30E+08,A_UNITS,"+-",0.0),
                 n=(1.69,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.78,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "21",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 23
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(6.28E+07,A_UNITS,"+-",0.0),
                 n=(1.75,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.51,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "22",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 24
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(3.06E+07,A_UNITS,"+-",0.0),
                 n=(1.87,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.59,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "23",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 25
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(8.70E+08,A_UNITS,"+-",0.0),
                 n=(1.39,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.07,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "24",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

ROH + H --> RO + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 26
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(5.48E+07,A_UNITS,"+-",0.0),
                 n=(1.82,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.44,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "25",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

HCHO + H --> HCO + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 27
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(8.07E+07,A_UNITS,"+-",0.0),
                 n=(1.76,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.67,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "26",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

RCHO + H --> RCO + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 28
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.53E+07,A_UNITS,"+-",0.0),
                 n=(1.98,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.78,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "27",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Primary vinylic {Cd/H2}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

R2C=CH2 + H --> R2C=CH + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 29
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(4.53E+03,A_UNITS,"+-",0.0),
                 n=(2.43,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.85,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "28",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Ketene hydrogen {CCO/H2}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 30
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.46E+07,A_UNITS,"+-",0.0),
                 n=(2.09,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.49,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "29",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Allene hydrogen {Cd/H2/Ca}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 31
rate(
  group1 = 
"""
Cd/H/NonDeC
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.98E+07,A_UNITS,"+-",0.0),
                 n=(1.95,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.65,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "30",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

RCH=CR2 + H --> RC=CR2 + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 32
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(4.33E+05,A_UNITS,"+-",0.0),
                 n=(2.38,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "31",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

R2C=CRCH3 + H --> R2C=CRCH2 + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 33
rate(
  group1 = 
"""
C/H2/OneDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 {Cd,Ct,CO,Cb} 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(6.99E+05,A_UNITS,"+-",0.0),
                 n=(2.36,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.11,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "32",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Primary allylic hydrogen {C/Cd/C/H2}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

RR2C=CRCH2R + H --> R2C=CRCHR + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 34
rate(
  group1 = 
"""
C/H2/OneDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 {Cd,Ct,CO,Cb} 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(7.79E+07,A_UNITS,"+-",0.0),
                 n=(1.78,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.11,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "33",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Primary propergylic {C/Ct/C/H2}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

RCCCH2R + H --> RCCCHR + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 35
rate(
  group1 = 
"""
C/H/Cs2
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(3.02E+06,A_UNITS,"+-",0.0),
                 n=(2.16,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.45,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "34",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Secondary allylic hydrogen {C/Cd/C2/H}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

R2C=CRCHR2 + H --> R2C=CRCR2 + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 36
rate(
  group1 = 
"""
C/H/Cs2
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.21E+08,A_UNITS,"+-",0.0),
                 n=(1.72,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.73,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "35",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Secondary propergylic {C/Ct/C2/H}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

RCCCHR2 + H --> RCCCR2 + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 37
rate(
  group1 = 
"""
C/H2/TwoDe
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 {Cd,Ct,CO,Cb} 0 {1,S}
5 {Cd,Ct,CO,Cb} 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(7.09E+03,A_UNITS,"+-",0.0),
                 n=(2.85,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-1.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "36",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

R2C=CH-CH2-CH=CR2 + H --> R2C=CH-CH-CH=CR2 + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 38
rate(
  group1 = 
"""
Cd/H/OneDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.93E+08,A_UNITS,"+-",0.0),
                 n=(1.74,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.28,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "37",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Dienylic {Cd/Cd/H}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

R2C=CRCH=CR2 + H --> R2C=CRC=CR2 + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 39
rate(
  group1 = 
"""
Cd/H/OneDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.18E+06,A_UNITS,"+-",0.0),
                 n=(2.4,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.11,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "38",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Eneynic {Cd/Ct/H}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

RCC-CH=CR2 + H --> RCC-C=CR2 + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 40
rate(
  group1 = 
"""
Ct_H
1 *1 C 0 {2,T}, {3,S}
2 C 0 {1,T}
3 *2 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.65E+08,A_UNITS,"+-",0.0),
                 n=(1.85,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.52,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "39",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

RCCH + H --> RCC + H2

NOTE: There was a discrepancy in the rate values. The published values were: A = 1.30E+08, n = 1.88, 

E0 = 1.34E+04

RMG values: A=1.65E+08, n=1.85, E0=	26.52.

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 41
rate(
  group1 = 
"""
C/H3/Ct
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Ct 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.70E+07,A_UNITS,"+-",0.0),
                 n=(1.91,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.99,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "40",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

RCCCH3 + H --> RCCCH2 + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 42
rate(
  group1 = 
"""
Cd/H/NonDeO
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 O 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(9.67E+09,A_UNITS,"+-",0.0),
                 n=(1.23,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.69,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "41",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 43
rate(
  group1 = 
"""
O/H/OneDe
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.30E+11,A_UNITS,"+-",0.0),
                 n=(0.82,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.75,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "42",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 44
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.24E+03,A_UNITS,"+-",0.0),
                 n=(2.58,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.04,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "43",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 45
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.04E+01,A_UNITS,"+-",0.0),
                 n=(2.92,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.16,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "44",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 46
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.54E+05,A_UNITS,"+-",0.0),
                 n=(1.89,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.97,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "45",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 47
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.92E+04,A_UNITS,"+-",0.0),
                 n=(2.29,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.44,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "46",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 48
rate(
  group1 = 
"""
O/H/OneDe
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.07E+04,A_UNITS,"+-",0.0),
                 n=(2.04,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.85,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "47",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Olefinic alcohol {O/Cd/H}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 49
rate(
  group1 = 
"""
O/H/OneDe
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(3.30E+08,A_UNITS,"+-",0.0),
                 n=(1.56,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.94,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "48",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom. Acid O-H {O/CO/H}",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
Sumathi, R.; Carstensen, H.-H.; Green, W.H. Jr.; J. Phys. Chem. A. 2001, 105, 8978

Table 8: Modified ArrHenius Fitted Parameters for GA Predicted Rates.

RCOOOH + H --> RCOOO + H2

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 50
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.53E+03,A_UNITS,"+-",0.0),
                 n=(2.43,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.85,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "49",
  short_comment = "Sumathy CBS-Q calculations. Rate expression per H atom.",
  long_comment = 
"""
Sumathy CBS-Q calculations. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 51
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.03E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "50",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom. 
Same reaction as #19. 

Saeys, M.; Reyniers, M.-F.; Marin, G.B.; Van Speybroeck, V.; Waroquier, M. J. Phys. Chem. A 2003, 107, 9147 - 9159.

CH3 + CH4 --> CH4 + CH3

pg 9156 Table 6: Calculated and Experimental Activation Energies(kJ/mol) at 0 K, deltaE (0 k), 

for Three Families of Radical Reactions from Various Levels of Theory.

From reference: E0 = 71.0/4.185 = 16.97, @ 0 K, from database: E0 = 19.0 @ 300 - 1500 K

Experimental values from reference @ 0 K = 55.4 kJ/mol, 60.7 kJ/mol, 61.9 kJ/mol
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 52
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.80E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "51",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 53
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.15E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "52",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 54
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.71E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "53",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 55
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.02E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "54",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 56
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.02E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "55",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 57
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.42E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "56",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 58
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.62E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "57",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 59
rate(
  group1 = 
"""
C/H2/OneDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 {Cd,Ct,CO,Cb} 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.73E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "58",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 60
rate(
  group1 = 
"""
C/H/Cs2
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.73E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "59",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 61
rate(
  group1 = 
"""
C/H2/TwoDe
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 {Cd,Ct,CO,Cb} 0 {1,S}
5 {Cd,Ct,CO,Cb} 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.54E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "60",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 62
rate(
  group1 = 
"""
C/H/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.32E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "61",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 63
rate(
  group1 = 
"""
C/H3/Cb
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cb 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.90E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "62",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 64
rate(
  group1 = 
"""
C/H2/OneDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 {Cd,Ct,CO,Cb} 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.70E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "63",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 65
rate(
  group1 = 
"""
C/H/Cs2
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.17E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "64",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 66
rate(
  group1 = 
"""
C/H3/Ct
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Ct 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.60E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "65",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 67
rate(
  group1 = 
"""
C/H2/OneDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 {Cd,Ct,CO,Cb} 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.81E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "66",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 68
rate(
  group1 = 
"""
C/H2/TwoDe
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 {Cd,Ct,CO,Cb} 0 {1,S}
5 {Cd,Ct,CO,Cb} 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.38E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "67",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 69
rate(
  group1 = 
"""
C/H/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.43E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "68",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 70
rate(
  group1 = 
"""
C/H/Cs2
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.03E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "69",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 71
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.88E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "70",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 72
rate(
  group1 = 
"""
Cd/H/NonDeC
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.43E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "71",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 73
rate(
  group1 = 
"""
Cd/H/OneDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.78E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "72",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 74
rate(
  group1 = 
"""
Cd/H/OneDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.26E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "73",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 75
rate(
  group1 = 
"""
Cd/H/OneDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.37E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "74",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 76
rate(
  group1 = 
"""
Ct_H
1 *1 C 0 {2,T}, {3,S}
2 C 0 {1,T}
3 *2 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.33E+18,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(28.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "75",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 77
rate(
  group1 = 
"""
Cb_H
1 *1 Cb 0 {2,B}, {3,B}, {4,S}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
4 *2 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.17E+15,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "76",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 78
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.87E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "77",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 79
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.51E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "78",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 80
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.46E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "79",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 81
rate(
  group1 = 
"""
C/H2/OneDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 {Cd,Ct,CO,Cb} 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.98E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "80",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 82
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.87E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "81",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 83
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.03E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "82",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 84
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cd
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  kf = Arrhenius(A=(1.53E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "83",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 85
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.59E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "84",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 86
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.86E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "85",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 87
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(5.31E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "86",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 88
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.02E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "87",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 89
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cd
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  kf = Arrhenius(A=(5.02E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "88",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 90
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.95E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "89",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 91
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.43E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "90",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 92
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.03E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "91",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 93
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.44E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "92",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 94
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H/OneDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.46E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(21.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "93",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 95
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs2
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.10E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(21.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "94",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 96
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
Cd_rad/NonDeC
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.02E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "95",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 97
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
Cd_rad/OneDe
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  kf = Arrhenius(A=(4.16E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "96",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 98
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Ct
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Ct 0 {1,S}
""",
  kf = Arrhenius(A=(3.14E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(17.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "97",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 99
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H/OneDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.53E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "98",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 100
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs2
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(9.78E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "99",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 101
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(3.39E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "100",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 102
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.88E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "101",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 103
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cd
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  kf = Arrhenius(A=(5.36E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "102",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 104
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.47E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "103",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 105
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.93E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "104",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 106
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(7.82E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "105",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 107
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.25E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "106",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 108
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H/OneDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.51E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(31.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "107",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 109
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs2
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.39E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "108",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 110
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
Cd_rad/NonDeC
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(6.04E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "109",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 111
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
Cd_rad/OneDe
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  kf = Arrhenius(A=(1.30E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "110",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 112
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Ct
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Ct 0 {1,S}
""",
  kf = Arrhenius(A=(2.55E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "111",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 113
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H/OneDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.41E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(28.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "112",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 114
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs2
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.59E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "113",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 115
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(3.13E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "114",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 116
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.62E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "115",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 117
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cd
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  kf = Arrhenius(A=(5.78E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(21.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "116",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 118
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.73E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "117",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 119
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.18E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "118",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 120
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.60E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "119",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 121
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.87E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "120",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 122
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_rad/H/OneDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.52E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(21.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "121",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 123
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs2
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.13E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(21.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "122",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 124
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
Cd_rad/NonDeC
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.52E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "123",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 125
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
Cd_rad/OneDe
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  kf = Arrhenius(A=(1.95E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "124",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 126
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Ct
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Ct 0 {1,S}
""",
  kf = Arrhenius(A=(7.58E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "125",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 127
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_rad/H/OneDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.09E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "126",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 128
rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs2
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.32E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "127",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 129
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.37E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "128",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 130
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.95E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "129",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 131
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cd
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  kf = Arrhenius(A=(2.87E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(25.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "130",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 132
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.49E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "131",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 133
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.08E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "132",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 134
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.57E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "133",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 135
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.04E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "134",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 136
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H/OneDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.82E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(27.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "135",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 137
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs2
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.30E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(27.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "136",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 138
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
Cd_rad/NonDeC
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.01E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "137",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 139
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
Cd_rad/OneDe
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  kf = Arrhenius(A=(2.37E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "138",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 140
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Ct
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Ct 0 {1,S}
""",
  kf = Arrhenius(A=(2.91E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(23.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "139",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 141
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H/OneDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.34E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(24.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "140",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 142
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs2
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.20E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(24.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "141",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 143
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.93E+06,A_UNITS,"+-",0.0),
                 n=(1.80,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.431,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "142",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels. Fixed by RWest (changed to per H)",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253. http://dx.doi.org/10.1016/S0010-2180(01)00373-X

Rate expressions for H atom abstraction from fuels. 
pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:OH, Site: primary (a)
Verified by Karma James

**HOWEVER** This entry should probably use the numbers for primary(d) not primary(a).
Primary(a) is for a primary on neopentane; primary(d) is for a primary on propane.
Richard West. (Updated accordingly).

These numbers reported by Curran et al. were apparently taken from
N. Cohen, *Intl. J. Chem. Kinet.* 14 (1982), p. 1339 http://dx.doi.org/10.1002/kin.550141206

Rate expression is changed to per H.(divided by 3)
Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 144
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.50E+05,A_UNITS,"+-",0.0),
                 n=(2.00,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-1.133,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "143",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels. (changed to per H)",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
http://dx.doi.org/10.1016/S0010-2180(01)00373-X

Rate expressions for H atom abstraction from fuels. 
pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:OH, Site: secondary (b)

Verified by Karma James

These numbers reported by Curran et al. were apparently taken from
N. Cohen, *Intl. J. Chem. Kinet.* 14 (1982), p. 1339 http://dx.doi.org/10.1002/kin.550141206


Rate expression is changed to per H.(divided by 2)
Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 145
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.70E+06,A_UNITS,"+-",0.0),
                 n=(1.90,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-1.451,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "144",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
http://dx.doi.org/10.1016/S0010-2180(01)00373-X

Rate expressions for H atom abstraction from fuels.
pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:OH, Site: tertiary (c)

Verified by Karma James

These numbers reported by Curran et al. were apparently taken from
N. Cohen, *Intl. J. Chem. Kinet.* 14 (1982), p. 1339 http://dx.doi.org/10.1002/kin.550141206
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 146
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(9.50E+02,A_UNITS,"+-",0.0),
                 n=(3.05,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.123,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "145",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels. (changed to per H)",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:O, Site: primary (a)

Verified by Karma James

Rate expression is changed to per H.(divided by 9)
Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 147
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(2.39E+04,A_UNITS,"+-",0.0),
                 n=(2.71,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.106,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "146",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels. (changed to per H)",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:O, Site: secondary (b)

Verified by Karma James


Rate expression is changed to per H.(divided by 2)
Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 148
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(3.83E+05,A_UNITS,"+-",0.0),
                 n=(2.41,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.140,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "147",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:O, Site: tertiary (c)

Verified by Karma James


This rate parameter actually comes from following new mechanism for PRF.

https://www-pls.llnl.gov/data/docs/science_and_technology/chemistry/combustion/prf_2d_mech.txt

Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 149
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeO
1 *3 O 1 {2,S}
2    O 0 {1,S}
""",
  kf = Arrhenius(A=(2.80E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.435,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "148",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels. (changed to per H)",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:HO2, Site: primary (a)
Verified by Karma James

Rate expression is changed to per H.(divided by 9)
Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 150
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeO
1 *3 O 1 {2,S}
2    O 0 {1,S}
""",
  kf = Arrhenius(A=(2.80E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(17.686,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "149",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels. (changed to per H)",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:HO2, Site: secondary (b)

Verified by Karma James

Rate expression is changed to per H.(divided by 2)
Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 151
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeO
1 *3 O 1 {2,S}
2    O 0 {1,S}
""",
  kf = Arrhenius(A=(2.80E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.013,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "150",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:HO2, Site: tertiary (c)

Verified by Karma James
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 152
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeC
1 *3 O 1 {2,S}
2 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.27E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.000,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "151",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels. (changed to per H)",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:CH3O, Site: primary (a)

Verified by Karma James

Rate expression is changed to per H.(divided by 9)
Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 153
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeC
1 *3 O 1 {2,S}
2 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.50E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.000,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "152",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels. (changed to per H)",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:CH3O, Site: secondary (b)

Verified by Karma James

Rate expression is changed to per H.(divided by 2)
Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 154
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeC
1 *3 O 1 {2,S}
2 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.90E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.800,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "153",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:CH3O, Site: tertiary (c)

Verified by Karma James	 
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 155
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(7.00E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(50.76,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "154",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels. (changed to per H)",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:O2, Site: primary (a)

Verified by Karma James

Rate expression is changed to per H.(divided by 9)
Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 156
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(7.00E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(48.21,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "155",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels. (changed to per H)",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:O2, Site: secondary (b)

Verified by Karma James

Rate expression is changed to per H.(divided by 2)
Yushi Suzuki
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 157
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(7.00E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(46.06,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "156",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:O2, Site: tertiary (c)

Verified by Karma James	 
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 158
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(7.25E+13,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(56.64,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,800),
  rank = 4,
  old_id = "157",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
H2 + O2 --> H + HO2 C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 1091, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 3,2.

Verified by Karma James

pg. 1109: Discussion of evaluated data

Recommended value computed using reverse rate and thermodynamics

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 159
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.73E+03,A_UNITS,"+-",0.0),
                 n=(2.56,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.03,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (200,3000),
  rank = 3,
  old_id = "158",
  short_comment = "Knyazev et al. [119] Transition state theory.",
  long_comment = 
"""
[119] Knyazev, V.D; Bencsura, A.; Stoliarov, S.I.; Slagle, I.R. J. Phys. Chem. 1996, 100, 11346.
H2 + C2H3 --> H + C2H4 C.D.W divided original rate expression by 2 ( from A = 9.45E+03), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 160
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.11E+04,A_UNITS,"+-",0.0),
                 n=(2.48,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.13,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,3500),
  rank = 3,
  old_id = "159",
  short_comment = "Mebel et al. [120] Transition state theory.",
  long_comment = 
"""
[120] Mebel, A.M.; Morokuma, K.; Lin, M.C. J Chem. Phys. 1995, 103, 3440.
H2 + C2H3 --> H + C2H4 C.D.W divided original rate expression by 2, to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 161
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.58E+09,A_UNITS,"+-",0.0),
                 n=(0.70,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.11,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "160",
  short_comment = "Weissman et al. [121] Transition state theory.",
  long_comment = 
"""
[121] Weissman, M.A.; Benson, S.W. J. Phys. Chem. 1988, 92, 4080.
H2 + C2H3 --> H + C2H4 C.D.W divided original rate expression by 2 ( from A = 3.15E+09), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 162
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
Ct_rad
1 *3 C 1 {2,T}
2 C 0 {1,T}
""",
  kf = Arrhenius(A=(5.4E+12,A_UNITS,"*/",3.16),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.17,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "161",
  short_comment = "Baulch et al. [94] literature review.",
  long_comment = 
"""
[94] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Frank, P.; Hayman, G,; Just, T.; Kerr, J.A.; Murrells, T.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1994, 23, 847.

H2 + C2H --> H + C2H2 C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 863 Evaluated Kinetic Data for Combustion Modelling Supplement 1, Table 1. Bimolecular reactions - C2H Radical Reactions.

Verified by Karma James

pg.1013-1014: Discussion on evaluated data

C2H+H2-->C2H2+H: Recommended rate coefficient is that reported by Koshi et al.  Rate

coefficient was computed for low temperatures, but extrapolation to higher temperatures
fits other reported data reasonably well.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 163
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
Cb_rad
1 *3 Cb 1 {2,B}, {3,B}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
""",
  kf = Arrhenius(A=(2.86E+04,A_UNITS,"+-",0.0),
                 n=(2.43,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.28,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,5000),
  rank = 3,
  old_id = "162",
  short_comment = "Mebel et al. [122] Transition state theory.",
  long_comment = 
"""
[122] Mebel, A.M.; Lin, M.C.; Yu, T.; Morokuma, K. J. Phys. Chem. A. 1997, 101, 3189.
H2 + phenyl --> H + benzene C.D.W divided original rate expression by 2 ( from A = 5.71E+04), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 164
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
CO_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.0E+05,A_UNITS,"*/",5.0),
                 n=(2.00,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(17.83,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "163",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
H2 + HCO --> H + CH2O C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 1094, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 15,2.

Verified by Karma James

pg. 1147: Discussion of evaluated data

Recommended value computed using reverse rate and thermodynamics

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 165
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
CO_rad/NonDe
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(2.06E+06,A_UNITS,"*/",3.0),
                 n=(1.82,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(17.61,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "164",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
H2 + CH3CO --> H + CH3CHO C.D.W divided original rate expression by 2, to get rate expression per H atom.

//WAS UNABLE TO VERIFY DATA!!! DATA NOT FOUND IN REFERENCE.

pg. 1229: Discussion on evaluated data

No experimental data for forward rxn, at the time

Reviewers noticed that k(H+HCHO=H2+HCO) / k(H+CH3CHO=H2+CH3CO) ~ 2, due to double the number of H atoms available

Used 0.5*k(H+HCHO=H2+HCO) and equilibrium constant to compute recommended rate expression

Verified by MRH on 10Aug2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 166
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.10E+08,A_UNITS,"+-",0.0),
                 n=(1.21,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.71,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (200,2400),
  rank = 3,
  old_id = "165",
  short_comment = "Isaacson [123] Transition state theory.",
  long_comment = 
"""
[123] Isaacson, A.D. J. Chem. Phys. 1997, 107, 3832.
H2 + O2 --> H + H2O C.D.W divided original rate expression by 2, to get rate expression per H atom.

166. [100] Jodkowski, J.T.; Rauez, M.-T.; Rayez, J.-C. J. Phys. Chem. A. 1999, 103, 3750.

H2 + CH3O --> H + CH3OH The calculated reverse rate constants are in good agreement with experiment. (This is -R1 in the paper)

C.D.W divided original rate expression by 2, to get rate expression per H atom.

Verified by Greg Magoon; maximum error of fitted expression from tabular data for forward rate constant, kr1 is 15% (cf. p. 3758)
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 167
rate(
  group1 = 
"""
H2
1  *1 H 0 {2,S}
2  *2 H 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeC
1 *3 O 1 {2,S}
2 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(6.32E-02,A_UNITS,"+-",0.0),
                 n=(4.00,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.91,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 2,
  old_id = "166",
  short_comment = "Jodkowski et al. [100] ab initio calculations.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 168
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(9.925E+12,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(56.83,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (500,2000),
  rank = 4,
  old_id = "167",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

CH4 + O2 --> CH3 + HO2 C.D.W divided original rate expression by 4, to get rate expression per H atom.

pg 417 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - O2 Reactions.

Verified by Karma James

pg.483: Discussion on evaluated data

O2+CH4 --> HO2+CH3: Recommended data based on experimental value for CH2O + O2 -->

HO2 + HCO.  Assumes equal A factor per C-H bond and Ea = deltaH.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 169
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.16E-02,A_UNITS,"*/",2.0),
                 n=(4.14,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.56,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "168",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
CH4 + C2H5 --> CH3 + C2H6 C.D.W divided original rate expression by 4, to get rate expression per H atom.

//WAS UNABLE TO VERIFY DATA!!! DATA NOT FOUND IN REFERENCE.

pg. 1177: Discussion on evaluated data

No experimental data for forward rxn, at the time

Recommended data from reverse rate and equilibrium constant

Verified by MRH on 10Aug2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 170
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.81E-04,A_UNITS,"*/",2.0),
                 n=(4.40,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.79,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "169",
  short_comment = "Tsang et al. [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
CH4 + iso-C3H7 --> CH3 + C3H8 C.D.W divided original rate expression by 4, to get rate expression per H atom.

pg 894, Chemical Kinetic Database For Combustion Chemistry, Part 3. Index of Reactions and Summary of Recommended Rate Expressions. No. 42,10.

Verified by Karma James

pg. 935: Discussion on evaluated data

Entry 42,10: No data available at the time.  Author recommends rate coefficient

expression based on reverse rate and equilibrium constant.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 171
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
Ct_rad
1 *3 C 1 {2,T}
2 C 0 {1,T}
""",
  kf = Arrhenius(A=(4.53E+11,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "170",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
CH4 + C2H --> CH3 + C2H2 C.D.W divided original rate expression by 4, to get rate expression per H atom.

pg 1101, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 21,10.

Verified by Karma James

pg. 1220: Discussion of evaluated data

Recommended data is expression given by Brown and Laufer (1981).

They computed the pre-exponential factor by the bond energy-bond order (BEBO) method

and combined that with experimental k at room temperature to yield Arrhenius expression
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 172
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
Cb_rad
1 *3 Cb 1 {2,B}, {3,B}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
""",
  kf = Arrhenius(A=(5.0E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (560,1410),
  rank = 2,
  old_id = "171",
  short_comment = "Heckmann et al. [124]",
  long_comment = 
"""
[124] Heckmann, E.; Hippler, H. Troe, J. Sypm. Int. Combust. Proc. 1996, 26, 543.
Absolute value measured directly (excitation technique: thermal, analytical technique: vis-UV absorption) CH4 + phenyl --> benzene

C.D.W divided original rate expression by 4, to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 173
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
CO_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.82E+03,A_UNITS,"*/",5.0),
                 n=(2.85,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(22.46,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "172",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
CH4 + HCO --> CH3 + CH2O C.D.W divided original rate expression by 4, to get rate expression per H atom.

pg 1094, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 15,10.

Verified by Karma James

pg. 1150: Discussion on evaluated data

Recommended data computed using reverse rate and equilibrium constant

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 174
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
CO_rad/NonDe
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(5.43E+02,A_UNITS,"*/",5.0),
                 n=(2.88,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(21.46,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "173",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
CH4 + CH3CO --> CH3 + CH3CHO C.D.W divided original rate expression by 4, to get rate expression per H atom.

pg 1102, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 22,10.

Verified by Karma James

pg. 1231: Discussion on evaluated data

Recommended number computed from reverse rate and equilibrium constant

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 175
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.85E-01,A_UNITS,"+-",0.0),
                 n=(3.95,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.55,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (223,2400),
  rank = 3,
  old_id = "174",
  short_comment = "Melissas and Truhlar [125] Transition state theory.",
  long_comment = 
"""
[125] Melissas, V.S.; Truhlar, D.G. J. Chem. Phys. 1993,99,1010.
CH4 + OH --> CH3 + H2O C.D.W divided original rate expression by 4, to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 176
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.93E+06,A_UNITS,"*/",1.41),
                 n=(1.83,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.78,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (240,2500),
  rank = 4,
  old_id = "175",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

CH4 + OH --> CH3 + H2O C.D.W divided original rate expression by 4, to get rate expression per H atom.

pg 419 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - OH Radical Reactions.

Verified by Karma James

pg.571-572: Discussion on evaluated data

OH+CH4 --> H2O+CH3: \"The preferred value of k is that obtained experimentally by

Madronich and Felder which predicts very precisely the data obtained between
240 and 2000K.\"
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 177
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.55E+07,A_UNITS,"+-",0.0),
                 n=(1.60,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.12,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,1510),
  rank = 3,
  old_id = "176",
  short_comment = "Cohen et al. [101] Transition state theory.",
  long_comment = 
"""
[101] Cohen, N. Int. J. Chem. Kinet. 1991, 23, 397.
CH4 + OH --> CH3 + H2O C.D.W divided original rate expression by 4, to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 178
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeC
1 *3 O 1 {2,S}
2 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.55E-04,A_UNITS,"+-",0.0),
                 n=(5.00,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.58,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 2,
  old_id = "177",
  short_comment = "Jodkowski et al. [100] ab initio calculations.",
  long_comment = 
"""
[100] Jodkowski, J.T.; Rauez, M.-T.; Rayez, J.-C. J. Phys. Chem. A. 1999, 103, 3750.
CH4 + CH3O --> CH3 + CH3OH The calculated reverse rate constants are in good agreement with experiment. (Rxn. -R3 in paper)

C.D.W divided original rate expression by 4 ( from A= 1.51E+09), to get rate expression per H atom.

Verified by Greg Magoon; cf. reverse reaction, #261, below
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 179
rate(
  group1 = 
"""
C_methane
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeO
1 *3 O 1 {2,S}
2    O 0 {1,S}
""",
  kf = Arrhenius(A=(4.53E+10,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.58,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "178",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
CH4 + HO2 --> CH3 + H2O2 C.D.W divided original rate expression by 4, to get rate expression per H atom.

pg 1093, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 10,7.

Verified by Karma James

pg. 1131: Discussion on evaluated data

Recommended data is based on expression for HO2 attach on alkanes (Walker)

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 180
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(1.005E+13,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(51.87,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (500,2000),
  rank = 4,
  old_id = "179",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

C2H6 + O2 --> C2H5 + HO2 C.D.W divided original rate expression by 6, to get rate expression per H atom.

pg 417 Evaluated Kinetic Data for Combustion Modelling  Table 1. Bimolecular reactions - O2 Reactions. (The value for E0 does not 

match the value in the reference, E0 RMG = 1.87; E0 Reference = 51.86)

Verified by Karma James

pg.484: Discussion on evaluated data

O2+C2H6 --> HO2+C2H5: \"The value given in the Walker review has been modified slightly

to allow for the higher heat of formation of the C2H5 radical now recommended
and for an assumed equal A factor per C-H bond in CH2O+O2 and C2H6+O2.\"
*** NOTE: MRH agrees with KJ on discrepancy in RMG-stored E0.  MRH is changing the value

of E0 in RMG from 1.87 kcal/mol to 51.87 kcal/mol. ***
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 181
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
Ct_rad
1 *3 C 1 {2,T}
2 C 0 {1,T}
""",
  kf = Arrhenius(A=(6.02E+11,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 4,
  old_id = "180",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
C2H6 + C2H --> C2H5 + C2H2 C.D.W divided original rate expression by 6, to get rate expression per H atom.

pg 1101, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 21,11.

Verified by Karma James

pg. 1221: Discussion on evaluated data

Recommended data is based on expression given by Brown and Laufer (1981).

Brown and Laufer calculated pre-exponential factor by BEBO method and
combined calculation with experimental measurement of k at room temperature.
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 182
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
Cb_rad
1 *3 Cb 1 {2,B}, {3,B}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
""",
  kf = Arrhenius(A=(3.48E+10,A_UNITS,"*/",2.35),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.44,E_UNITS,"+-",0.18)
                 ),
  temperature_range = (565,1000),
  rank = 2,
  old_id = "181",
  short_comment = "Park et al. [126]",
  long_comment = 
"""
[126] Park, J.; Gheyas, S.; Lin, M.C. Int. J. Chem. Kinet. 2001, 33, 64.
Absolute value measured directly. Static or low flow, flash photolysis excitation, Vis-UV absoprtion analysis. 

Phenyl radicals are produced from 193 nm photolysis of C6H5COCH3. The cavity ringdown spectroscopy and/or mass spectroscopy

have been used to monitor reactant and/or products. C2H6 + phenyl --> C2H5 + benzene.

C.D.W divided original rate expression by 6 ( from A= 2.09E+11), to get rate expression per H atom. Original delta A = 2.0E+10.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 183
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
CO_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.82E+03,A_UNITS,"*/",5.0),
                 n=(2.72,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.24,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "182",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
C2H6 + HCO --> C2H5 + CH2O C.D.W divided original rate expression by 6(from A = 4.69E+04), to get rate expression per H atom.

pg 1094, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 15,11.

Verified by Karma James

pg. 1150: Discussion on evaluated data

Recommended data computed from reverse rate and equilibrium constant

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 184
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
CO_rad/NonDe
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(3.02E+03,A_UNITS,"*/",5.0),
                 n=(2.75,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(17.53,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "183",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
C2H6 + CH3CO --> C2H5 + CH3CHO C.D.W divided original rate expression by 6(from A = 1.81E+04), to get rate expression per H atom.

pg 1102, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 22,11.

Verified by Karma James

pg. 1231: Discussion on evaluated data

Recommended data computed using rate of C2H5+CH2O divided by 2 (since only one O=C-H

hydrogen is present in CH3CHO) and equilibrium constant
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 185
rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.20E+06,A_UNITS,"*/",1.41),
                 n=(2.00,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.86,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (250,2000),
  rank = 4,
  old_id = "184",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

C2H6 + OH --> C2H5 + H2O C.D.W divided original rate expression by 6, to get rate expression per H atom.

pg 420 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - OH Radical Reactions.

Verified by Karma James

pg.589-590: Discussion on evaluated data

OH+C2H6 --> H2O+C2H5: \"The preferred value of k is almost indistinguishable from the

value obtained by Cohen from transition state calculations carried out for
temperatures between 300 and 2000K.\"
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 186
rate(
  group1 = 
"""
C/H3/CO
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 CO 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.17E+05,A_UNITS,"+-",0.0),
                 n=(2.20,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (295,600),
  rank = 3,
  old_id = "185",
  short_comment = "Taylor et al. [127] Transition state theory.",
  long_comment = 
"""
[127] Taylor, P.H.; Rahman, M.S.; Arif, M.; Dellinger, B.; Marshall, P. Sypm. Int. Combust. Proc. 1996, 26, 497.
CH3CHO + OH --> CH2CHO + H2O Rate constant is high pressure limit (pressure 0.13-0.97atm?) 

C.D.W divided original rate expression by 3(from A = 1.55E+06), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 187
rate(
  group1 = 
"""
C/H3/O
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 O 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.05E-04,A_UNITS,"+-",0.0),
                 n=(4.90,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.72,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 2,
  old_id = "186",
  short_comment = "Jodkowski et al. [100] ab initio calculations.",
  long_comment = 
"""
[100] Jodkowski, J.T.; Rauez, M.-T.; Rayez, J.-C. J. Phys. Chem. A. 1999, 103, 3750.
CH3OH + CH3 --> CH2OH + CH4 The calculated rate constants are in good agreement with experiment. (Rxn. R4 in paper)

C.D.W divided original rate expression by 3 ( from A= 8.43E+08), to get rate expression per H atom.

Verified by Greg Magoon
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 188
rate(
  group1 = 
"""
C/H3/O
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 O 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.14E+03,A_UNITS,"+-",0.0),
                 n=(2.80,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.42,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 2,
  old_id = "187",
  short_comment = "Jodkowski et al. [100] ab initio calculations.",
  long_comment = 
"""
[100] Jodkowski, J.T.; Rauez, M.-T.; Rayez, J.-C. J. Phys. Chem. A. 1999, 103, 3750.
CH3OH + OH --> CH2OH + H2O The calculated rate constants are in good agreement with experiment. (Rxn. R6 in paper)

C.D.W divided original rate expression by 3 ( from A= 2.11E+11), to get rate expression per H atom.

Verified by Greg Magoon
**Note that R2 from this paper appears to be missing from the RMG library, so I have added it as 100_R2**

100_R2: [100] Jodkowski, J.T.; Rauez, M.-T.; Rayez, J.-C. J. Phys. Chem. A. 1999, 103, 3750.

CH3OH + H --> CH2OH + H2 (Rxn. R2 in paper)

divided original rate expression by 3 to get rate expression per H atom.

Created by Greg Magoon; maximum error of fitted expression from tabular data for kr2 is 20% (cf. p. 3758); rank of 2 assigned based on rank for other values reported in the paper in the rateLibrary (also 2)
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 189
rate(
  group1 = 
"""
C/H3/O
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 O 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(4.51E+02,A_UNITS,"+-",0.0),
                 n=(3.20,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.49,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 2,
  old_id = "100_R2",
  short_comment = "Jodkowski et al. [100] ab initio calculations. added by Greg Magoon 08/25/09",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 190
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(1.985E+13,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(47.69,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "188",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
C3H8 + O2 --> iso-C3H7 + HO2  C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 891, Chemical Kinetic Database For Combustion Chemistry, Part 3. Index of Reactions and Summary of Recommended Rate Expressions. No. 40,3.

//NOTE: For A value, Database value = 1.985E+13 and Reference value = 1.65E+13

Verified by Karma James

NOTE: MRH computed Reference A value of 1.99E+13 (11Aug2009)

pg. 899: Discussion on evaluated data

Entry 40,3 (b): No data available at the time.  The author \"estimates\" the rate

coefficient expressions (no indication of how).
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 191
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
CH2_triplet
1 *3 C 2T {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.55E-01,A_UNITS,"*/",10.0),
                 n=(3.46,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.47,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "189",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
C3H8 + CH2 --> iso-C3H7 + CH3  C.D.W divided original rate expression by 2(from A = 1.51), to get rate expression per H atom.

pg 892, Chemical Kinetic Database For Combustion Chemistry, Part 3. Index of Reactions and Summary of Recommended Rate Expressions. No. 40,26.
Verified by Karma James

pg. 910: Discussion on evaluated data

Entry 40,26 (b): No data available at the time.  Author estimates the rate coefficient

expression as that of CH3+C3H8=i-C3H7+CH4.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 192
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(2.39E+04,A_UNITS,"*/",2.0),
                 n=(2.71,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.11,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "190",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
C3H8 + O --> iso-C3H7 + OH  C.D.W divided original rate expression by 2(from A = 4.77E+04), to get rate expression per H atom.

pg 891, Chemical Kinetic Database For Combustion Chemistry, Part 3. Index of Reactions and Summary of Recommended Rate Expressions. No. 40,5.

Verified by Karma James

pg. 901: Discussion on evaluated data

Entry 40,5 (b): The author notes \"considerable scatter\" among the existing data.  The

author computed Arrhenius A and n parameters using a BEBO calculation and performed
a \"fit\" on the data reported by Herron and Huie to obtain the Arrhenius E.  This
rate coefficient expression is stated to fit 3 (of the 5) raw data reported.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 193
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/O
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 O 0 {1,S}
""",
  kf = Arrhenius(A=(3.02E+01,A_UNITS,"*/",5.0),
                 n=(2.95,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.98,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "191",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
C3H8 + CH2OH --> iso-C3H7 + CH3OH  C.D.W divided original rate expression by 2(from A = 6.03E+01), to get rate expression per H atom.

//WAS UNABLE TO VERIFY DATA!!! DATA NOT FOUND IN REFERENCE.

pg. 910: Discussion on evaluated data

Entry 40,39 (b)

No experimental data, at the time

Recommended value for C3H8+CH2OH-->n-C3H7+CH3OH comes from rate for C2H6+CH2OH-->C2H5+CH3OH

No discussion on where rate for C3H8+CH2OH-->i-C3H7+CH3OH comes from:

A is ~ factor of 3 smaller (6 hydrogens vs 2 ... seems reasonable to MRH)
E is 1 kcal/mol smaller (more stable to form secondary radical than primary)
Verified by MRH on 10Aug2009

MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 194
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.10E+02,A_UNITS,"*/",10.0),
                 n=(3.10,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.82,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "192",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
C3H8 + C2H3 --> iso-C3H7 + C2H4  C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 891, Chemical Kinetic Database For Combustion Chemistry, Part 3. Index of Reactions and Summary of Recommended Rate Expressions. No. 40,19.

Verified by Karma James

pg. 906: Discussion on evaluated data

Entry 40,19 (b): No data available at the time.  The author recommends the rate coefficient

expression of C2H3+C2H6=C2H5+C2H4 for the rxn C2H3+C3H8=n-C3H7+C2H4.  The author
assumes the ratio of secondary-to-primary H-atom abstraction for the rxn CH3+C3H8
to obtain the recommended rate coefficient expression.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 195
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
Ct_rad
1 *3 C 1 {2,T}
2 C 0 {1,T}
""",
  kf = Arrhenius(A=(6.05E+11,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "193",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
C3H8 + C2H --> iso-C3H7 + C2H2  C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 891, Chemical Kinetic Database For Combustion Chemistry, Part 3. Index of Reactions and Summary of Recommended Rate Expressions. No. 40,21.

Verified by Karma James

pg. 906-907: Discussion on evaluated data

Entry 40,21 (b): No data available at the time.  The author recommends the rate coefficient

of C2H6+C2H=C2H2+C2H5 for the rxn C3H8+C2H=C2H2+n-C3H7.  Due to the high exothermicity
of the rxn, the author assumes the H-atom abstraction rxn is limited to the number
of H-atoms available, thus recommedning a rate coefficient equal to one-third that
recommended for C3H8+C2H=C2H2+n-C3H7.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 196
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
CO_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.4E+06,A_UNITS,"*/",3.0),
                 n=(1.90,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(17.01,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "194",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
C3H8 + HCO --> iso-C3H7 + CH2O  C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 891, Chemical Kinetic Database For Combustion Chemistry, Part 3. Index of Reactions and Summary of Recommended Rate Expressions. No. 40,15.

Verified by Karma James

pg. 904: Discussion on evaluated data

Entry 40,15 (b): No data available at the time.  The author recommends a rate coefficient

expression based on the reverse rxn (note: the author uses the rate of the rxn
n-C3H7+CH2O=HCO+C3H8 instead of i-C3H7+CH2O=HCO+C3H8) and equilibrium constant.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 197
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
CO_rad/NonDe
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(2.65E+06,A_UNITS,"*/",3.0),
                 n=(2.00,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.24,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "195",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
C3H8 + CH3CO --> iso-C3H7 + CH3CHO  C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 891, Chemical Kinetic Database For Combustion Chemistry, Part 3. Index of Reactions and Summary of Recommended Rate Expressions. No. 40,22.

Verified by Karma James

pg. 907: Discussion on evaluated data

Entry 40,22 (b): No data available at the time.  The author recommends a rate coefficient

expression based on the equilibrium constant and the reverse rate (note: the author
estimates this reverse rate using the suggestions of Kerr, J.A. and Trotman-Dickenson, A.F.).
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 198
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.95E+06,A_UNITS,"+-",0.0),
                 n=(1.90,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.16,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (295,1220),
  rank = 3,
  old_id = "196",
  short_comment = "Cohen et al. [101] Transition state theory.",
  long_comment = 
"""
[101] Cohen, N. Int. J. Chem. Kinet. 1991, 23, 397.
C3H8 + OH --> iso-C3H7 + H20  C.D.W divided original rate expression by 2, to get rate expression per H atom.

Not yet checked
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 199
rate(
  group1 = 
"""
C/H2/NonDeC
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeC
1 *3 O 1 {2,S}
2 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(7.25E+10,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.57,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "197",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
C3H8 + CH3O --> iso-C3H7 + CH3OH  C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 891, Chemical Kinetic Database For Combustion Chemistry, Part 3. Index of Reactions and Summary of Recommended Rate Expressions. No. 40,24.

Verified by Karma James

pg. 908: Discussion on evaluated data

Entry 40,24 (b): The author assumes the Arrhenius A parameter should follow:

A(C3H8+CH3O=CH3OH+n-C3H7) / A(C3H8+CH3O=CH3OH+i-C3H7) = 3
The author also demands that the recommended data fit both sets of experiments
reported.  These assumptions (plus one unwritten one, as we still have 3
unknown parameters [A1, E1, E2 ... A2=f(A1)]) produce the reported rate
coefficient expression.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 200
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(3.97E+13,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(43.92,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "198",
  short_comment = "Tsang [92] literature review.",
  long_comment = 
"""
[92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
Iso-C4H10 + O2 --> tert-C4H9 + HO2

pg 5, Chemical Kinetic Database For Combustion Chemistry, Part 4 - Isobutane. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 43,3.

Verified by Karma James

pg.14: Discussion on evaluated data

Entry 43,3(b): No direct measurements at the time.  A review article reported a rate

coefficient expression for iC4H10+O2-->tC4H9+HO2.  An experiment performed by
Baldwin, Drewery, and Walker reported a rate coefficient expression for O2 abstracting
a tertiary H-atom from 2,3-dimethylbutane.  The experiment\'s value matched well
with the review\'s value, so Tsang recommended the review\'s value.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 201
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(1.57E+05,A_UNITS,"*/",2.0),
                 n=(2.50,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.11,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "199",
  short_comment = "Tsang [92] literature review.",
  long_comment = 
"""
[92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
Iso-C4H10 + O --> tert-C4H9 + OH

pg 5, Chemical Kinetic Database For Combustion Chemistry, Part 4 - Isobutane. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 43,5.

Verified by Karma James

pg.15: Discussion on evaluated data

Entry 43,5(b): Michael et al. reported the rate coefficient expression for iC4H10+O=OH+C4H9 isomer.

Tsang subtracted from this expression the contributions from iC4H10+O=OH+iC4H9 (What
rate expression was used?  The one recommended here?) to obtain an expression for
iC4H10+O=OH+tC4H9.  Tsang then adjusted the rate expression such that the A-factor\'s
temperature dependence was 2.5 (was this 2.5 based on the review by Cohen and Westberg?).
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 202
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
CH2_triplet
1 *3 C 2T {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.09E+12,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.91,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "200",
  short_comment = "Tsang [92] literature review.",
  long_comment = 
"""
[92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
Iso-C4H10 + CH2 --> tert-C4H9 + CH3

pg 6, Chemical Kinetic Database For Combustion Chemistry, Part 4 - Isobutane. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 43,25.

Verified by Karma James

pg.23-24: Discussion on evaluated data

Entry 43,25(b): Tsang recommends the rate coefficient expression reported by Bohland et al.

Tsang notes that the rate for CH2_triplet abstracting a H-atom is faster than
the recommended value for CH3 abstracting a H-atom.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 203
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.04E-01,A_UNITS,"*/",5.0),
                 n=(3.46,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "201",
  short_comment = "Tsang [92] literature review.",
  long_comment = 
"""
[92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
Iso-C4H10 + C2H3 --> tert-C4H9 + C2H4

pg 5, Chemical Kinetic Database For Combustion Chemistry, Part 4 - Isobutane. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 43,19.

Verified by Karma James

pg.20: Discussion on evaluated data

Entry 43,19(b): No data available at the time.  Author recommends rate coefficient expression

based on the rxn CH3+iC4H10=CH4+tC4H9: same Arrhenius A and n parameters, Ea decreased
by 8.5 kJ/mol.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 204
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
Ct_rad
1 *3 C 1 {2,T}
2 C 0 {1,T}
""",
  kf = Arrhenius(A=(6.62E+11,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "202",
  short_comment = "Tsang [92] literature review.",
  long_comment = 
"""
[92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
Iso-C4H10 + C2H --> tert-C4H9 + C2H2

pg 5, Chemical Kinetic Database For Combustion Chemistry, Part 4 - Isobutane. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 43,21.

Verified by Karma James

pg.21: Discussion on evaluated data

Entry 43,21(b): No data available at the time.  For the rxn iC4H10+C2H=C2H2+iC4H9, author

recommends 1.5x the rate of the rxn C2H6+C2H=C2H2+C2H5 (9 vs. 6 primary H-atoms).
The author then recommends a rate coefficient for iC4H10+C2H=C2H2+tC4H9 that appears
to be 1/9 the rate of iC4H10+C2H=C2H2+iC4H9 (9 vs. 1 H-atom).
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 205
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
CO_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.43E+04,A_UNITS,"*/",5.0),
                 n=(2.50,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.51,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "203",
  short_comment = "Tsang [92] literature review.",
  long_comment = 
"""
[92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
Iso-C4H10 + HCO --> tert-C4H9 + CH2O

pg 5, Chemical Kinetic Database For Combustion Chemistry, Part 4 - Isobutane. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 43,15.

Verified by Karma James

pg.18: Discussion on evaluated data

Entry 43,15(b): No data available at the time.  For the rxn iC4H10+HCO=CH2O+iC4H9, author

recommends 1.5x the rate of the rxn C3H8+HCO+CH2O+nC3H7 (9 vs. 6 primary H-atoms).
The author then recommends the rate coefficient for iC4H10+HCO=CH2O+tC4H9 to be the 
rate coefficient of iC4H10+HCO=CH2O+iC4H9, with the A factor divided by 9 (9 vs. 1
H-atoms) and the Ea decreased by 20 kJ/mol.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 206
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
CO_rad/NonDe
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(3.43E+04,A_UNITS,"*/",10.0),
                 n=(2.50,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.51,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "204",
  short_comment = "Tsang [92] literature review.",
  long_comment = 
"""
[92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
Iso-C4H10 + CH3CO --> tert-C4H9 + CH3CHO

pg 5, Chemical Kinetic Database For Combustion Chemistry, Part 4 - Isobutane. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 43,22.

Verified by Karma James

pg.21: Discussion on evaluated data

Entry 43,22(b): No data available at the time.  Author recommends rate coefficient expression

based on the rxn iC4H10+HCO=CH2O+tC4H9.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 207
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.57E+06,A_UNITS,"+-",0.0),
                 n=(1.90,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.45,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,1150),
  rank = 3,
  old_id = "205",
  short_comment = "Cohen et al. [101] Transition state theory.",
  long_comment = 
"""
[101] Cohen, N. Int. J. Chem. Kinet. 1991, 23, 397.
Iso-C4H10 + OH --> tert-C4H9 + H2O

Not yet checked
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 208
rate(
  group1 = 
"""
C/H/Cs3
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeC
1 *3 O 1 {2,S}
2 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.29E+10,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.88,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "206",
  short_comment = "Tsang [92] literature review.",
  long_comment = 
"""
[92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
Iso-C4H10 + CH3O --> tert-C4H9 + CH3OH

pg 6, Chemical Kinetic Database For Combustion Chemistry, Part 4 - Isobutane. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 43,24.

Verified by Karma James

pg.22: Discussion on evaluated data

Entry 43,24(b): A study by Berces and Trotman-Dickenson reported a rate coefficient

for the rxn iC4H10+CH3O=CH3OH+C4H9 isomer.  Tsang used the rate coefficient for
the rxn CH3O+C(CH3)4=CH3OH+H2C*-C(CH3)3 determined by Shaw and Trotman-Dickenson
as a characteristic for CH3O+primary H-atom abstraction to recommend a rate coefficient
for iC4H10+CH3O=CH3OH+iC4H9.  This rate expression was subtracted from the rate
coefficient reported by Berces and Trotman-Dickenson to yield the rate coefficient
for iC4H10+CH3O=CH3OH+tC4H9.  Lastly, the pre-exponential factor was cut in half,
due to Tsang\'s correcting an arithmetic error by Kerr and Parsonage (perhaps this
work was referenced in the Berces and Trotman-Dickenson study?)
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 209
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(1.055E+13,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(57.63,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "207",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
C2H4 + O2 --> C2H3 + HO2 C.D.W divided original rate expression by 4, to get rate expression per H atom.

pg 1097, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 18,3.

Verified by Karma James

pg. 1184: Discussion on evaluated data

Recommended data follows Walker\'s estimates for O2+alkane

Note: The authors note that a lower lying channel, involving addition and
rearrangement prior to decomposition, may exist.
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 210
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(3.78E+06,A_UNITS,"+-",0.0),
                 n=(1.91,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.74,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (290,1510),
  rank = 3,
  old_id = "209",
  short_comment = "Mahmud et al. [128] Transition state theory",
  long_comment = 
"""
[128] Mahmud, K.; Marshall, P.; Fontijn, A. J Phys. Chem. 1987, 91, `568.
C2H4 + O --> C2H3 + OH C.D.W divided original rate expression by 4(from A= 1.51E+07), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 211
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.58E+02,A_UNITS,"*/",10.0),
                 n=(3.13,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "210",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
C2H4 + C2H5 --> C2H3 + C2H6 C.D.W divided original rate expression by 4, to get rate expression per H atom.

pg 1098, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 18,17.

Verified by Karma James

pgs. 1191-1192: Discussion on evaluated data

Recommended data based on study performed by Halstead and Quinn

Tsang fit the data against BEBO calculations (to attain the Arrhenius A, n)
and manually adjusted the E.
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 212
rate(
  group1 = 
"""
Cd_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.13E+12,A_UNITS,"*/",3.16),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.94,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (650,1500),
  rank = 4,
  old_id = "211",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

C2H4 + OH --> C2H3 + H2O C.D.W divided original rate expression by 4(from A= 2.05E+13), to get rate expression per H atom.

pg 420 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - OH Radical Reactions.

Verified by Karma James

pg.586-587: Discussion on evaluated data

OH+C2H4 --> H2O+C2H3: Recommended rate taken from expression reported by Tully (1988).

MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 213
rate(
  group1 = 
"""
Cd/H/NonDeC
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(6.02E+10,A_UNITS,"*/",3.0),
                 n=(0.70,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.63,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "212",
  short_comment = "Tsang [93] literature review.",
  long_comment = 
"""
[93] Tsang, W. J. Phys. Chem. Ref. Data 1991, 20, 221.
CH3CH=CH2 + O --> CH3C=CH2 + OH

pg 233-234: Discussion on evaluated data

Verified by MRH on 6Aug2009

Entry 46,5(f): No measurements on H-atom abstraction rxns. Recommended rate coefficient

is computed as follows:

The rate of O + C3H6 --> OH + H2C=CH-*CH2 is computed using the expression:
[k(O+C2H6-->C2H5+HO)/k(OH+C2H6-->C2H5+H2O)] * k(OH+C3H6-->H2C=CH-*CH2+H2O).
The author uses this expression because he notes that OH and O H-atom abstraction
rxns generally follow the same trend.  The O+C2H6, OH+C2H6, and OH+C3H6
are from other Tsang review articles.
The rate of O+C3H6 --> OH+CH3C=CH2 is computed by adjusting the O+C3H6 --> OH+H2C=CH-*CH2
rate coefficient: increasing the Ea/R by 880 Kelvin and multiplying the A
by ~0.345; these values come from the relative difference between the rxns
OH+C3H6 --> H2O+H2C=CH-*CH2 and OH+C3H6 --> H2O+CH3C=CH2
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 214
rate(
  group1 = 
"""
Cd/H/NonDeC
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(4.09e+05,A_UNITS,"*/",4.0),
                 n=(2.5,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.79,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "213",
  short_comment = "Tsang [93] literature review.",
  long_comment = 
"""
[93] Tsang, W. J. Phys. Chem. Ref. Data 1991, 20, 221.
CH3CH=CH2 + H --> CH3C=CH2 + H2

pg 231: Discussion on evaluated data

Previous modified Arrhenius parameters were for RELATIVE rate (kc/ka)

Multipled kc/ka by ka to get kc (only one H to abstract, so no division necessary)

Certified by MRH on 6Aug2009

Entry 46,4(c): No data available for H-atom abstraction rxns.  The recommended rate

coefficient is based on the author\'s assumption that methyl substitution has the
same influence in olefins as in alkanes.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 215
rate(
  group1 = 
"""
Cd/H/NonDeC
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(0.842,A_UNITS,"*/",6.0),
                 n=(3.5,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(11.66,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "214",
  short_comment = "Tsang [93] literature review.",
  long_comment = 
"""
[93] Tsang, W. J. Phys. Chem. Ref. Data 1991, 20, 221.
CH3CH=CH2 + CH3 --> CH3C=CH2 + CH4

pg 237-239

Previous modified Arrhenius parameters were for RELATIVE rate (ke/kc)

Multiplied ke/kc by kc to get ke (only one H to abstract, so no division necessary)

Certified by MRH on 6Aug2009

Entry 46,16(e): Recommended rate coefficient is based on the author\'s assumption

that methyl substitution has the same influence in olefins as in alkanes.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 216
rate(
  group1 = 
"""
Cd/H/NonDeC
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(0.842,A_UNITS,"*/",6.0),
                 n=(3.5,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(9.67,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "215",
  short_comment = "Tsang [93] literature review.",
  long_comment = 
"""
[93] Tsang, W. J. Phys. Chem. Ref. Data 1991, 20, 221.
CH3CH=CH2 + C2H3 --> CH3C=CH2 + C2H4

pg 240-241

Previous modified Arrhenius parameters were for RELATIVE rate (kc/ka)

Multiplied kc/ka by ka to get kc (only one H to abstract, so no division necessary)

Certified by MRH on 6Aug2009

Entry 46,19(c): No data available at the time.  The recommended rate coefficient

is based on the rate expressions for CH3 abstracting a H-atom from C3H6; all of
the Ea\'s have been decreased by 4kJ/mol.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 217
rate(
  group1 = 
"""
Cd/H/NonDeC
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
Ct_rad
1 *3 C 1 {2,T}
2 C 0 {1,T}
""",
  kf = Arrhenius(A=(1.21E+12,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "216",
  short_comment = "Tsang [93] literature review.",
  long_comment = 
"""
[93] Tsang, W. J. Phys. Chem. Ref. Data 1991, 20, 221.
CH3CH=CH2 + C2H --> CH3C=CH2 + C2H2 

pg 241

Verified by MRH on 6Aug2009

Entry 46,21(e): No data available at the time.  Recommended rate expression is \"somewhat

smaller\" than the rate of rxn C3H6+C2H --> C2H2+H2C=CH-*CH2.  The rate of this rxn
is assumed to be the rate of the rxn C2H+C2H6 --> C2H2+C2H5.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 218
rate(
  group1 = 
"""
Cd/H/NonDeC
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 C 0 {1,D}
3 *2 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.11E+06,A_UNITS,"*/",2.0),
                 n=(2.00,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.45,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "217",
  short_comment = "Tsang [93] literature review.",
  long_comment = 
"""
[93] Tsang, W. J. Phys. Chem. Ref. Data 1991, 20, 221.
CH3CH=CH2 + OH --> CH3C=CH2 + H2O

pg 235-236

Valid T range in reference suggested 700-2500K

Uncertainty stated in reference was *2.0

Verified by MRH on 6Aug2009

Entry 46,6(d): No direct measurements of H-atom abstraction rxns.  The recommended

H-atom abstraction rxns are based on \"the results of Tully (1988) for the rxn
of OH + C2H4 and the rate constant ratio of OH + primary Hydrogens in ethane by
Tully et al. (1986) to OH + secondary Hydrogens by Droege and Tully (1986)\".  The
author has also introduced a T^2 dependence in the A-factor.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 219
rate(
  group1 = 
"""
Ct_H
1 *1 C 0 {2,T}, {3,S}
2 C 0 {1,T}
3 *2 H 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(6.05E+12,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(74.52,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "218",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
C2H2 + O2 --> C2H + HO2 C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 1100, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 20,3.

Verified by Karma James

pg. 1209: Discussion on evaluated data

Recommended data based on report by Walker

NOTE: Authors note that a lower-lying channel of O2 addition, rearrangement,
and decomposition may exist.
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 220
rate(
  group1 = 
"""
Ct_H
1 *1 C 0 {2,T}, {3,S}
2 C 0 {1,T}
3 *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.36E+11,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(23.45,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "220",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
C2H2 + C2H5 --> C2H + C2H6 C.D.W divided original rate expression by 2 (from A= 2.71E+11), to get rate expression per H atom.

pg 1100, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 20,17.

Verified by Karma James

pg. 1215: Discussion on evaluated data

Recommended data based on reverse rate and equilibrium constant

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 221
rate(
  group1 = 
"""
Ct_H
1 *1 C 0 {2,T}, {3,S}
2 C 0 {1,T}
3 *2 H 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.25E+03,A_UNITS,"*/",10.0),
                 n=(2.68,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.04,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "221",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
C2H2 + OH --> C2H + H2O C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 1100, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 20,6.

Verified by Karma James

pg. 1213: Discussion on evaluated data

Recommended data is derived from BEBO method calculation

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 222
rate(
  group1 = 
"""
Cb_H
1 *1 Cb 0 {2,B}, {3,B}, {4,S}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
4 *2 H 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(1.052E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(60.01,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (1200,1700),
  rank = 4,
  old_id = "222",
  short_comment = "Asaba et al. [129]. Data are estimated.",
  long_comment = 
"""
[129] Asaba, T.; Fujii, N.; Proc. Int. Symp. Shock Tubes Waves 1971, 8, 1.
Benzene + O2 --> phenyl + HO2 C.D.W divided original rate expression by 6(from A = 6.31E+13), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 223
rate(
  group1 = 
"""
Cb_H
1 *1 Cb 0 {2,B}, {3,B}, {4,S}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
4 *2 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.00E+08,A_UNITS,"+-",0.0),
                 n=(1.80,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.35,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (500,800),
  rank = 4,
  old_id = "223",
  short_comment = "Mebel et al. [122] RRK(M) extrapolation.",
  long_comment = 
"""
[122] Mebel, A.M.; Lin, M.C.; Yu, T.; Morokuma, K. J. Phys. Chem. A. 1997, 101, 3189.
Rate constant is high pressure limit. Benzene + H --> phenyl + H2

C.D.W divided original rate expression by 6(from A = 6.02E+08), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 224
rate(
  group1 = 
"""
Cb_H
1 *1 Cb 0 {2,B}, {3,B}, {4,S}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
4 *2 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(5.02E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.11,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,1000),
  rank = 5,
  old_id = "224",
  short_comment = "Nicovich et al. [130]",
  long_comment = 
"""
[130] Nicovich, J.M.; Ravishankara, A.R. J. Phys. Chem. 1984, 88, 2534.
Pressure 0.01-0.26 atm Excitation: flash photolysis, analysis: resonance fluorescence. Benzene + H --> phenyl + H2

C.D.W divided original rate expression by 6(from A = 3.01E+12), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 225
rate(
  group1 = 
"""
Cb_H
1 *1 Cb 0 {2,B}, {3,B}, {4,S}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
4 *2 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.05E+11,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.86,E_UNITS,"+-",1.19)
                 ),
  temperature_range = (650,770),
  rank = 5,
  old_id = "225",
  short_comment = "Zhang et al. [131]",
  long_comment = 
"""
[131] Zhang, H.X.; Ahonkhai, S.I. Back, H.M. Can. J. Chem. 1989, 67, 1541.
Pressure 0.30-0.50 atm Excitation: thermal, analysis: GC. Benzene + C2H5 --> phenyl + C2H6

C.D.W divided original rate expression by 6(from A = 6.31E+11), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 226
rate(
  group1 = 
"""
Cb_H
1 *1 Cb 0 {2,B}, {3,B}, {4,S}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
4 *2 H 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.72E+07,A_UNITS,"*/",2.0),
                 n=(1.42,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.45,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (400,1500),
  rank = 4,
  old_id = "226",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

Benzene + OH --> phenyl + H2O  C.D.W divided original rate expression by 6(from A = 1.63E+08), to get rate expression per H atom.

pg 420 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - OH Radical Reactions.

Verified by Karma James

pg.595-597: Discussion on evaluated data

OH+C6H6 --> H2O+C6H5: Authors note that this rxn should be dominant at temperatures

above 500K.  No other comment on where the recommended rate expression comes from
(although MRH believes it is a best-fit to the available data, based on graph).
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 227
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(2.34E+07,A_UNITS,"+-",0.0),
                 n=(2.05,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.93,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2200),
  rank = 3,
  old_id = "227",
  short_comment = "Michael et al. [132] Transition state theory.",
  long_comment = 
"""
[132] Michael, J.V.; Kumaran, S.S.; Su, M.-C. J. Phys. Chem. A. 1999, 103, 5942.
CH2O + O2 --> HCO + HO2 C.D.W divided original rate expression by 2, to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 228
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(2.08E+11,A_UNITS,"*/",2.0),
                 n=(0.57,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.76,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (250,2200),
  rank = 4,
  old_id = "228",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

CH2O + O --> HCO + OH C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 416 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - O Atom Reactions.

Verified by Karma James

pg.449-450: Discussion on evaluated data

O+CH2O --> OH+HCO: \"The preferred values are based on the low temperature data which are

all in good agreement, and on the higher temperature value of Bowman.\"
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 229
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
CH2_triplet
1 *3 C 2T {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.02E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "229",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
Rate constant is an upper limit. CH2O + CH2 --> HCO + CH3

C.D.W divided original rate expression by 2 (from A= 6.03E+09), to get rate expression per H atom.

pg 1106, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 26,12.

Verified by Karma James

pg. 1267: Discussion on evaluated data

Recommended data based on triplet methylene\'s general lack of reactivity in H-atom abstractions

NOTE: Rate coefficient is an upper limit
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 230
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.89E-08,A_UNITS,"*/",1.58),
                 n=(6.10,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.97,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 4,
  old_id = "230",
  short_comment = "Baulch et al. [94] literature review.",
  long_comment = 
"""
[94] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Frank, P.; Hayman, G,; Just, T.; Kerr, J.A.; Murrells, T.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1994, 23, 847.

CH2O + CH3 --> HCO + CH4 C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 862 Evaluated Kinetic Data for Combustion Modelling Supplement 1, Table 1. Bimolecular reactions - CH3 Radical Reactions.

Verified by Karma James

pg.989-990: Discussion on evaluated data

CH3+HCHO --> CH4+HCO: The recommended value is a \"best fit to the data of Choudhury et al.,

the reworked data from Anastasi, together with those at lower temperatures from
Refs. 4, 5, and 7.\"
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 231
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.75E+03,A_UNITS,"*/",5.0),
                 n=(2.81,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.86,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "231",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
CH2O + C2H5 --> HCO + C2H6 C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 1096, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 17,12.

Verified by Karma James

pg. 1178: Discussion on evaluated data

Recommended data is the rate of CH2O+CH3-->HCO+CH4.

Authors note that rate coefficients for alkyl radicals w/aldehydic H-atoms are
similar (as noted by Kerr, J.A. and Trotman-Dickenson, A.F.
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 232
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.40E+10,A_UNITS,"*/",2.5),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.96,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "232",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
CH2O + iso-C3H7 --> HCO + C3H8 C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 894, Chemical Kinetic Database For Combustion Chemistry, Part 3. Index of Reactions and Summary of Recommended Rate Expressions. No. 42,12.

Verified by Karma James

pg. 936: Discussion on evaluated data

Entry 42,12: No data available at the time.  The author recommends a rate coefficient

expression that is twice that of the rxn i-C3H7+(CH3)2CHCHO, taken from a study
by Kerr, J.A. and Trotman-Dickenson, A.F. (1959).  The author states that a correction
was made to the 1959 report, taking the recommended rate of i-C3H7 recombination
(reported by Tsang) into consideration.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 233
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.63E+09,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.56,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "233",
  short_comment = "Tsang [92] literature review.",
  long_comment = 
"""
[92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
CH2O + tert-C4H9 --> HCO + iso-C4H10 C.D.W divided original rate expression by 2 (from A= 3.25E+09), to get rate expression per H atom.

pg 7, Chemical Kinetic Database For Combustion Chemistry, Part 4 - Isobutane. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 44,12.

Verified by Karma James

pg.35: Discussion on evaluated data

Entry 44,12: No data available at the time.  The author recommends 2x the rate coefficient

of the rxn tC4H9+tC4H9-CHO=iC4H10+tC4H9-CO (2 vs. 1 aldehydic H-atoms); this value
was reported by Birrell and Trotman-Dickenson.  The author also notes that he has
taken into account tC4H9 combination (perhaps meaning he used a geometric mean rule
to derive the final form of the expression?)
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 234
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.71E+03,A_UNITS,"*/",5.0),
                 n=(2.81,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.86,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "234",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
CH2O + C2H3 --> HCO + C2H4 C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 1099, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 19,12.

Verified by Karma James

pg. 1197: Discussion on evaluated data

Recommended data is the rate of CH2O+CH3-->HCO+CH4.

Authors note that rate coefficients for alkyl radicals w/aldehydic H-atoms are
similar (as noted by Kerr, J.A. and Trotman-Dickenson, A.F.
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 235
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
CO_rad/NonDe
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(9.05E+10,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.92,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "235",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
CH2O + CH3CO --> HCO + CH3CHO C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 1102, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 22,12.

Verified by Karma James

pg. 1231: Discussion on evaluated data

Recommended data based on \"analogous systems\"

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 236
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.72E+09,A_UNITS,"*/",5.0),
                 n=(1.18,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.45,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,3000),
  rank = 4,
  old_id = "236",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

CH2O + OH --> HCO + H2O C.D.W divided original rate expression by 2 (from A= 3.43E+09), to get rate expression per H atom.

pg 419 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - OH Radical Reactions.

Verified by Karma James

pg.575-576: Discussion on evaluated data

OH+CH2O --> H2O+HCO: The recommended rate coefficient is the value reported by Tsang

and Hampson.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 237
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeC
1 *3 O 1 {2,S}
2 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.10E+10,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.98,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "237",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
CH2O + CH3O --> HCO + CH3OH C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 1104, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 24,12.

Verified by Karma James

pg. 1245: Discussion on evaluated data

Recommended data based on review by Gray, based on experiments performed by Hoare and Wellington.

Authors note that experimental conditions were such that rxn of interest was
in competition with the disproportionation of two CH3O radicals (CH3O+CH3O-->CH3OH+CH2O)
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 238
rate(
  group1 = 
"""
CO_pri
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeO
1 *3 O 1 {2,S}
2    O 0 {1,S}
""",
  kf = Arrhenius(A=(2.06E+04,A_UNITS,"+-",0.0),
                 n=(2.50,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.21,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (641,1600),
  rank = 4,
  old_id = "238",
  short_comment = "Eiteneer et al. [133] literature review.",
  long_comment = 
"""
[133] Eiteneer, B.; Yu, C.-L.; Goldenberg, M.; Frenklach, M. J. Phys. Chem. A. 1998, 102, 5196.
CH2O + HO2 --> HCO + H2O2 C.D.W divided original rate expression by 2 (from A= 4.11E+04), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 239
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(3.01E+13,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.15,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,1100),
  rank = 4,
  old_id = "239",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

CH3CHO + O2 --> CH3CO + HO2

pg 417 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - O2 Reactions.

Verified by Karma James

pg.485: Discussion on evaluated data

O2+CH3CHO --> HO2+CH3CO: \"For this evaluation we prefer the approach of Walker and

the recommended value is based on the best current deltaH298 value (=163.8 kJ/mol
using deltaHf(CH3CO)=11.0 kJ/mol and deltaHf(HO2)=14.6 kJ/mol) and A=5.0x10^-11
cm3/molecule/s.\"
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 240
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(5.0E+12,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.79,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 4,
  old_id = "240",
  short_comment = "Warnatz [134] literature review",
  long_comment = 
"""
[134] Warnatz, J. Rate coefficeints in the C/H/O system. In Combustion Chemistry, 1984; pp 197.
CH3CHO + O --> CH3CO + OH
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 241
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(4.0E+13,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.21,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 4,
  old_id = "241",
  short_comment = "Warnatz [134] literature review",
  long_comment = 
"""
[134] Warnatz, J. Rate coefficeints in the C/H/O system. In Combustion Chemistry, 1984; pp 197.
CH3CHO + H --> CH3CO + H2
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 242
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.99E-06,A_UNITS,"*/",2.0),
                 n=(5.64,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.46,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1250),
  rank = 4,
  old_id = "242",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

CH3CHO + CH3 --> CH3CO + CH4

pg 423 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - CH3 Radical Reactions.

Verified by Karma James

pg.671: Discussion on evaluated data

CH3+CH3CHO --> CH4+CH3CO: \"There are no direct studies of the kinetics of this reaction

and all of the k values are relative to methyl recombination ... The preferred values
are based on a line constructed through the mean of the low temperature data and the
data of Liu and Laidler and Colket et al.\"
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 243
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cd
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  kf = Arrhenius(A=(3.8E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.21,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (790,810),
  rank = 5,
  old_id = "243",
  short_comment = "Loser et al. [135] bond strength-bond length method.",
  long_comment = 
"""
[135] Loser, U.; Scherzer, K.; Weber, K. Z. Phys. Chem. (Leipzig) 1989, 270, 237.
CH3CHO + CH2CH=CH2 --> CH3CO + CH3CH=CH2
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 244
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.13E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.68,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (480,520),
  rank = 5,
  old_id = "244",
  short_comment = "Scherzer et al. [136] bond energy-bond order method.",
  long_comment = 
"""
[136] Scherzer, K.; Loser, U.; Stiller, W. Z. Phys. Chem. 1987, 27, 300.
CH3CHO + C2H3 --> CH3CO + C2H4
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 245
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.0E+06,A_UNITS,"+-",0.0),
                 n=(1.80,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-1.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (295,600),
  rank = 3,
  old_id = "245",
  short_comment = "Taylor et al. [127] Transition state theory.",
  long_comment = 
"""
[127] Taylor, P.H.; Rahman, M.S.; Arif, M.; Dellinger, B.; Marshall, P. Sypm. Int. Combust. Proc. 1996, 26, 497.
CH3CHO + OH --> CH3CO + H2O Pressure 0.13-0.97 atm. Rate constant is high pressure limit.

pg 501, Table 1, k2 = 2.00x10^6 T^1.8 exp(1300/RT)

Previous modified Arrhenius parameters had E=1.3 kcal/mol; it should be E=-1.3 kcal/mol

Certified by MRH on 6Aug2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 246
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.0E+13,A_UNITS,"*/",2.51),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 4,
  old_id = "246",
  short_comment = "Warnatz [134] literature review",
  long_comment = 
"""
[134] Warnatz, J. Rate coefficeints in the C/H/O system. In Combustion Chemistry, 1984; pp 197.
CH3CHO + OH --> CH3CO + H2O
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 247
rate(
  group1 = 
"""
CO/H/NonDe
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 O 0 {1,D}
3 *2 H 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeO
1 *3 O 1 {2,S}
2    O 0 {1,S}
""",
  kf = Arrhenius(A=(3.01E+12,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.92,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (900,1200),
  rank = 4,
  old_id = "247",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

CH3CHO + HO2 --> CH3CO + H2O2

pg 421 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - HO2 Radical Reactions.

Verified by Karma James

pg.614-615: Discussion on evaluated data

HO2+CH3CHO --> CH3CO+H2O2: \"The preferred expression is based on a value of 1.7x10^-14

cm3/molecule/s at 1050K from a study performed by Colket et al. and an assumed A
factor of 5.0x10^-12 cm3/molecule/s.\"
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 248
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
O2b
1 *3 O 1 {2,S}
2    O 1 {1,S}
""",
  kf = Arrhenius(A=(2.325E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(74.12,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1000),
  rank = 5,
  old_id = "248",
  short_comment = "Mayer et al. [137] Bond energy-bond order.",
  long_comment = 
"""
[137] Mayer, S.W.; Schieler, L. J. Phys. Chem. 1968, 72, 2628.
http://dx.doi.org/10.1021/j100853a066

H2O + O2 --> OH + HO2. 
C.D.W divided original rate expression by 2, to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 249
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(2.63E+09,A_UNITS,"+-",0.0),
                 n=(1.20,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(17.83,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,1000),
  rank = 3,
  old_id = "249",
  short_comment = "Karach et al. [138] Transition state theory.",
  long_comment = 
"""
[138] Karach, S.P.; Oscherov, V.I. J. Phys. Chem. 1999, 110, 11918.
H2O + O --> OH + OH. C.D.W divided original rate expression by 2 (from A= 2.95E+39), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 250
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(7.40E+04,A_UNITS,"+-",0.0),
                 n=(2.60,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.18,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 3,
  old_id = "250",
  short_comment = "Harding et al. [139] Transition state theory.",
  long_comment = 
"""
[139] Harding, L.B.; Wagner, A.F. Symp. Int. Combust. proc. 1989, 22, 983.
H2O + O --> OH + OH. C.D.W divided original rate expression by 2 (from A= 1.48E+05), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 251
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.26E+08,A_UNITS,"*/",1.6),
                 n=(1.60,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.32,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "251",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.

H2O + H --> OH + H2. C.D.W divided original rate expression by 2, to get rate expression per H atom.

pg 418 Evaluated Kinetic Data for Combustion Modelling, Table 1. Bimolecular reactions - H Atom Reactions.

NOTE: E0 Rference = 18.4, E0 RMG database = 19.32

Verified by Karma James

pg.504: Discussion on evaluated data

H+H2O --> OH+H2: \"The recommended rate coefficient is based on the spare high temperature

measurements and rate data of the reverse rxn combined with the equilibrium constant.\"
MRH agrees with Karma.  However, the discrepancy is small and NIST\'s online database Webbook

has an E = 19.32 kcal/mol.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 252
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.20E+00,A_UNITS,"+-",0.0),
                 n=(3.31,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.56,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "252",
  short_comment = "Ma et al. [140] Transition state theory.",
  long_comment = 
"""
[140] Ma, S.; Liu, R.; Sci. China Ser. S: 1996, 39, 37.
H2O + CH3 --> OH + CH4. C.D.W divided original rate expression by 2 (from A= 6.39), to get rate expression per H atom.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 253
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.42E+02,A_UNITS,"*/",1.6),
                 n=(2.90,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.86,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "253",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
H2O + CH3 --> OH + CH4. C.D.W divided original rate expression by 2 (from A= 4.83E+02), to get rate expression per H atom.

pg 1095, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 16,9.

Verified by Karma James

pg. 1163: Discussion on evaluated data

Recommended data based on reverse rate and equilibrium constant

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 254
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.70E+06,A_UNITS,"*/",2.0),
                 n=(1.44,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.27,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "254",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
H2O + C2H5 --> OH + C2H6. C.D.W divided original rate expression by 2 (from A= 3.39E+06), to get rate expression per H atom.

pg 1096, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 17,9.

Verified by Karma James

pg. 1177: Discussion on evaluated data

Recommended data based on reverse rate and equilibrium constant

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 255
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.42E+02,A_UNITS,"*/",5.0),
                 n=(2.90,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.86,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "255",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
H2O + C2H3 --> OH + C2H4. C.D.W divided original rate expression by 2 (from A= 4.83E+02), to get rate expression per H atom.

pg 1098, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 19,9.

Verified by Karma James

pg. 1196: Discussion on evaluated data

Recommended data based on expression for CH3+H2O=CH4+OH

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 256
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
CO_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.18E+08,A_UNITS,"*/",5.0),
                 n=(1.35,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.03,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "256",
  short_comment = "Tsang [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087.
H2O + HCO --> OH + CH2O. C.D.W divided original rate expression by 2 (from A= 2.35E+08), to get rate expression per H atom.

pg 1094, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 15,9.

Verified by Karma James

pg. 1150: Discussion on evaluated data

Recommended data based on reverse rate and equilibrium constant

MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 257
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.417E-07,A_UNITS,"+-",0.0),
                 n=(5.48,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.274,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (200,700),
  rank = 3,
  old_id = "257",
  short_comment = "Masgrau et al. [141] Transition state theory.",
  long_comment = 
"""
[141] Masgrau, L.; Gonzalez-Lafont, A.; Lluch, J.M. J. Phys. Chem. A. 1999, 103, 1044.
H2O + OH --> OH + H2O . C.D.W refitted their k(T) to get A, n, and Ea, and divided original rate expression by 2, to get rate expression per H atom.

pg 1050, Table 4, Section: HO + HOH = HOH + OH(1), Column k_ab_CVT/SCT

MRH computed modified Arrhenius parameters using least-squares regression: ln(k) = ln(A) + n*ln(T) - (E/R)*(1/T)

E: Multiplied (E/R) expression by 1.987e-3

A: exp(ln(A)), multiplied by 6.02e23 (to convert /molecule to /mol) and divided by 2 (to get rate per H atom)

Certified by MRH on 7Aug2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 258
rate(
  group1 = 
"""
O_pri
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeC
1 *3 O 1 {2,S}
2 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.74E-01,A_UNITS,"+-",0.0),
                 n=(3.80,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.49,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 2,
  old_id = "258",
  short_comment = "Jodkowski et al. [100] ab initio calculations.",
  long_comment = 
"""
[100] Jodkowski, J.T.; Rauez, M.-T.; Rayez, J.-C. J. Phys. Chem. A. 1999, 103, 3750.
H2O + CH3O --> OH + CH3OH C.D.W divided original rate expression by 2 (from A= 9.03E+08), to get rate expression per H atom.; This is Rxn. -R5 from mpaper

Verified by Greg Magoon: note that this reaction is endothermic; the reverse (R5), appears as #267, below
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 259
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
O_atom_triplet
1 *3 O 2T
""",
  kf = Arrhenius(A=(1.0E+13,A_UNITS,"*/",2.51),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.69,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1000),
  rank = 4,
  old_id = "259",
  short_comment = "Warnatz [134] literature review",
  long_comment = 
"""
[134] Warnatz, J. Rate coefficeints in the C/H/O system. In Combustion Chemistry, 1984; pp 197.
CH3OH + O --> CH3O + OH
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 260
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
CH2_triplet
1 *3 C 2T {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.44E+01,A_UNITS,"*/",3.0),
                 n=(3.10,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.94,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "260",
  short_comment = "Tsang [90] literature review.",
  long_comment = 
"""
[90] Tsang, W. J. Phys. Chem. Ref. Data 1987, 16, 471.
CH3OH + CH2 --> CH3O + CH3

pg 475, Chemical Kinetic Database For Combustion Chemistry, Part 2 - Methanol. 

//Index of Reactions and Summary of Recommended Rate Expressions. No. 38,25.

Verified by Karma James

Data for Rate Expression 38,26 (pg. 493)

Stated uncertainty factor is 3

Verified by MRH on 11Aug2009

Entry 38,26 (b): No data available at the time.  Author suggests the rate coefficient

expression for CH3+CH3OH=CH4+CH3O
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 261
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.70E-04,A_UNITS,"+-",0.0),
                 n=(4.70,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.78,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 2,
  old_id = "261",
  short_comment = "Jodkowski et al. [100] ab initio calculations.",
  long_comment = 
"""
[100] Jodkowski, J.T.; Rauez, M.-T.; Rayez, J.-C. J. Phys. Chem. A. 1999, 103, 3750.
The calculated rate constants are in good agreement with experiment. CH3OH + CH3 --> CH3O + CH4 (Rxn. R3 from paper)

Verified by Greg Magoon: I changed upper temperature to 2000 K (was 2500) in line with other reactions from same paper; note that according to the paper, this reaction is very slightly endothermic; the exothermic reverse (-R3) is included above as #177
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 262
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.44E+01,A_UNITS,"*/",3.0),
                 n=(3.10,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.94,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "262",
  short_comment = "Tsang [90] literature review.",
  long_comment = 
"""
[90] Tsang, W. J. Phys. Chem. Ref. Data 1987, 16, 471.
CH3OH + C2H5 --> CH3O + C2H6

pg 475, Chemical Kinetic Database For Combustion Chemistry, Part 2 - Methanol. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 38,17.

Verified by Karma James

pg. 489: Discussion on evaluated data

Entry 38,17 (b): No data available at the time.  Author notes ethyl radicals are known

to be considerably less reactive than methyl.  Author recommends the rate coefficient
expression of CH3+CH3OH=CH4+CH3O, with a slight adjustment (based on the observed
trends in methyl vs. ethyl radical reactivity with hydrocarbons).
MRH 30-Aug-2009

//263: [90] Tsang, W. J. Phys. Chem. Ref. Data 1987, 16, 471.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 263
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.45E+01,A_UNITS,"*/",5.0),
                 n=(3.10,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.33,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "263",
  short_comment = "Tsang [91] literature review.",
  long_comment = 
"""
[91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
CH3OH + iso-C3H7 --> CH3O + C3H8

//WAS UNABLE TO VERIFY DATA!!! DATA NOT FOUND IN REFERENCE.

Ref[90] was incorrect; rate matches that reported in Ref[91].

pg. 944: Discussion on evaluated data

Entry 42,38 (b)

No experimental data, at the time

Recommended rate is based on C2H5+CH3OH=C2H6+CH3O

Verified by MRH on 10Aug2009

MRH 30-Aug-2009

//264: [90] Tsang, W. J. Phys. Chem. Ref. Data 1987, 16, 471.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 264
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.51E+03,A_UNITS,"*/",10.0),
                 n=(1.80,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.36,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "264",
  short_comment = "Tsang [92] literature review.",
  long_comment = 
"""
[92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
CH3OH + tert-C4H9 --> CH3O + iso-C4H10

//WAS UNABLE TO VERIFY DATA!!! DATA NOT FOUND IN REFERENCE.

Ref[90] was incorrect; rate matches that reported in Ref[92].

pg.44: Discussion on evaluated data

Entry 44,38(b)

Reference reports reaction as: t-C4H9+CH3OH=t-C4H10+CH3O

This is a typo: no t-C4H10 molecule exists (should be i-C4H10)
No experimental data, at the time

Recommended rate is based on reverse rxn and equilibrium constant

Verified by MRH on 10Aug2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 265
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.44E+01,A_UNITS,"*/",10.0),
                 n=(3.10,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.94,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "265",
  short_comment = "Tsang [90] literature review.",
  long_comment = 
"""
[90] Tsang, W. J. Phys. Chem. Ref. Data 1987, 16, 471.
CH3OH + C2H3 --> CH3O + C2H4

pg 475, Chemical Kinetic Database For Combustion Chemistry, Part 2 - Methanol. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 38,19.

Verified by Karma James

pg. 489: Discussion on evaluated data

Entry 38,19 (b): No data available at the time.  Author recommends the rate coefficient

expression for CH3+CH3OH=CH4+CH3O.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 266
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
Ct_rad
1 *3 C 1 {2,T}
2 C 0 {1,T}
""",
  kf = Arrhenius(A=(1.21E+12,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "266",
  short_comment = "Tsang [90] literature review.",
  long_comment = 
"""
[90] Tsang, W. J. Phys. Chem. Ref. Data 1987, 16, 471.
CH3OH + C2H --> CH3O + C2H2

pg 475, Chemical Kinetic Database For Combustion Chemistry, Part 2 - Methanol. 

Index of Reactions and Summary of Recommended Rate Expressions. No. 38,21.

Verified by Karma James

pg. 490: Discussion on evaluated data

Entry 38,21 (b): No data available at the time.  Author recommends a rate coefficient

expression based on measurements of C2H+CH4 and C2H+C2H6 rxns
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 267
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.73E+01,A_UNITS,"+-",0.0),
                 n=(3.40,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-1.14,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 2,
  old_id = "267",
  short_comment = "Jodkowski et al. [100] ab initio calculations.",
  long_comment = 
"""
[100] Jodkowski, J.T.; Rauez, M.-T.; Rayez, J.-C. J. Phys. Chem. A. 1999, 103, 3750.
The calculated rate constants are in good agreement with experiment. CH3OH + OH --> CH3O + H2O (Rxn. R5 from paper)

Verified by Greg Magoon (cf. reverse, #258, above)
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 268
rate(
  group1 = 
"""
O/H/NonDeC
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.0E+13,A_UNITS,"*/",3.16),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 4,
  old_id = "268",
  short_comment = "Warnatz [134] literature review",
  long_comment = 
"""
[134] Warnatz, J. Rate coefficeints in the C/H/O system. In Combustion Chemistry, 1984; pp 197.
CH3OH + OH --> CH3O + H2O
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 269
rate(
  group1 = 
"""
H2O2
1 *1 O 0 {2,S} {3,S}
2 O 0 {1,S} {4,S}
3 *2 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
InChI=1/C4H9O/c1-2-3-4-5/h5H,1-4H2
1 *3 C 1 {2,S} {6,S} {7,S}
2 C 0 {1,S} {3,S} {8,S} {9,S}
3 C 0 {2,S} {4,S} {10,S} {11,S}
4 C 0 {3,S} {5,S} {12,S} {13,S}
5 O 0 {4,S} {14,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {2,S}
9 H 0 {2,S}
10 H 0 {3,S}
11 H 0 {3,S}
12 H 0 {4,S}
13 H 0 {4,S}
14 H 0 {5,S}
""",
  kf = Arrhenius(A=(2.88E+00,A_UNITS,"+-",0.0),
                 n=(3.16,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(0.75,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "301",
  short_comment = "MRH CBS-QB3 calculations w/o HR corrections",
  long_comment = 
"""
MRH CBS-QB3 calculations w/o HR corrections
H2O2 + *CH2CH2CH2CH2OH = nButanol + HO2

CBS-QB3 method was used to calculate electronic energy of reactants, products, and TS; frequencies were
calculated using B3LYP/CBSB7 method.  Arrhenius expression was computed using CanTherm: an asymmetric Eckart
tunneling correction was employed and the frequencies were scaled by 0.99 (as suggested by Montgomery et al.
J.Chem.Phys. 110 (1999) 2822-2827).  The external symmetry number for H2O2 was 2; the external symmetry number
for the remaining species and TS were set to 1.  The rate coefficient was computed at 600-2000K (in 200 K increments).
The computed pre-exponential factor was divided by 2 and this is the reported value.

For nButanol+HO2=H2O2+*CH2CH2CH2CH2OH:
Moc et al. (AIP Conference Proceedings (2009) 1148 161-164 \"The Unimolecular Decomposition
and H Abstraction Reactions by HO and HO2 from n-Butanol\") report reaction barriers and
enthalpies(0 K); our CBS-QB3 calculations are shown in comparison (all units are kcal/mol).
				G3		CCSD(T)/cc-pVTZ		CBS-QB3
Barrier:		18.8		19.62			17.57
Enthalpy:		14.25		14.66			13.70
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 270
rate(
  group1 = 
"""
H2O2
1 *1 O 0 {2,S} {3,S}
2 O 0 {1,S} {4,S}
3 *2 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
InChI=1/C4H9O/c1-2-3-4-5/h2,5H,3-4H2,1H3
1 C 0 {2,S} {6,S} {7,S} {8,S}
2 *3 C 1 {1,S} {3,S} {9,S}
3 C 0 {2,S} {4,S} {10,S} {11,S}
4 C 0 {3,S} {5,S} {12,S} {13,S}
5 O 0 {4,S} {14,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {1,S}
9 H 0 {2,S}
10 H 0 {3,S}
11 H 0 {3,S}
12 H 0 {4,S}
13 H 0 {4,S}
14 H 0 {5,S}
""",
  kf = Arrhenius(A=(6.75E-01,A_UNITS,"+-",0.0),
                 n=(3.42,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(1.43,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "302",
  short_comment = "MRH CBS-QB3 calculations w/o HR corrections",
  long_comment = 
"""
MRH CBS-QB3 calculations w/o HR corrections
H2O2 + CH3*CHCH2CH2OH = nButanol + HO2

CBS-QB3 method was used to calculate electronic energy of reactants, products, and TS; frequencies were
calculated using B3LYP/CBSB7 method.  Arrhenius expression was computed using CanTherm: an asymmetric Eckart
tunneling correction was employed and the frequencies were scaled by 0.99 (as suggested by Montgomery et al.
J.Chem.Phys. 110 (1999) 2822-2827).  The external symmetry number for H2O2 was 2; the external symmetry number
for the remaining species and TS were set to 1.  The rate coefficient was computed at 600-2000K (in 200 K increments).
The computed pre-exponential factor was divided by 2 and this is the reported value.

For nButanol+HO2=H2O2+CH3*CHCH2CH2OH:
Moc et al. (AIP Conference Proceedings (2009) 1148 161-164 \"The Unimolecular Decomposition
and H Abstraction Reactions by HO and HO2 from n-Butanol\") report reaction barriers and
enthalpies(0 K); our CBS-QB3 calculations are shown in comparison (all units are kcal/mol).
				G3		CCSD(T)/cc-pVTZ		CBS-QB3
Barrier:		14.64		15.47			14.72
Enthalpy:		11.05		12.41			10.11
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 271
rate(
  group1 = 
"""
H2O2
1 *1 O 0 {2,S} {3,S}
2 O 0 {1,S} {4,S}
3 *2 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
InChI=1/C4H9O/c1-2-3-4-5/h3,5H,2,4H2,1H3
1 C 0 {2,S} {6,S} {7,S} {8,S}
2 C 0 {1,S} {3,S} {9,S} {10,S}
3 *3 C 1 {2,S} {4,S} {11,S}
4 C 0 {3,S} {5,S} {12,S} {13,S}
5 O 0 {4,S} {14,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {1,S}
9 H 0 {2,S}
10 H 0 {2,S}
11 H 0 {3,S}
12 H 0 {4,S}
13 H 0 {4,S}
14 H 0 {5,S}
""",
  kf = Arrhenius(A=(3.145E-01,A_UNITS,"+-",0.0),
                 n=(3.52,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(1.61,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "303",
  short_comment = "MRH CBS-QB3 calculations w/o HR corrections",
  long_comment = 
"""
MRH CBS-QB3 calculations w/o HR corrections
H2O2 + CH3CH2*CHCH2OH = nButanol + HO2

CBS-QB3 method was used to calculate electronic energy of reactants, products, and TS; frequencies were
calculated using B3LYP/CBSB7 method.  Arrhenius expression was computed using CanTherm: an asymmetric Eckart
tunneling correction was employed and the frequencies were scaled by 0.99 (as suggested by Montgomery et al.
J.Chem.Phys. 110 (1999) 2822-2827).  The external symmetry number for H2O2 was 2; the external symmetry number
for the remaining species and TS were set to 1.  The rate coefficient was computed at 600-2000K (in 200 K increments).
The computed pre-exponential factor was divided by 2 and this is the reported value.

For nButanol+HO2=H2O2+CH3CH2*CHCH2OH:
Moc et al. (AIP Conference Proceedings (2009) 1148 161-164 \"The Unimolecular Decomposition
and H Abstraction Reactions by HO and HO2 from n-Butanol\") report reaction barriers and
enthalpies(0 K); our CBS-QB3 calculations are shown in comparison (all units are kcal/mol).
				G3		CCSD(T)/cc-pVTZ		CBS-QB3
Barrier:		15.43		16.37			16.33
Enthalpy:		13.53		14.02			11.48
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 272
rate(
  group1 = 
"""
H2O2
1 *1 O 0 {2,S} {3,S}
2 O 0 {1,S} {4,S}
3 *2 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
InChI=1/C4H9O/c1-2-3-4-5/h4-5H,2-3H2,1H3
1 C 0 {2,S} {6,S} {7,S} {8,S}
2 C 0 {1,S} {3,S} {9,S} {10,S}
3 C 0 {2,S} {4,S} {11,S} {12,S}
4 *3 C 1 {3,S} {5,S} {13,S}
5 O 0 {4,S} {14,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {1,S}
9 H 0 {2,S}
10 H 0 {2,S}
11 H 0 {3,S}
12 H 0 {3,S}
13 H 0 {4,S}
14 H 0 {5,S}
""",
  kf = Arrhenius(A=(1.485E+00,A_UNITS,"+-",0.0),
                 n=(3.39,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(1.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "304",
  short_comment = "MRH CBS-QB3 calculations w/o HR corrections",
  long_comment = 
"""
MRH CBS-QB3 calculations w/o HR corrections
H2O2 + CH3CH2CH2*CHOH = nButanol + HO2

CBS-QB3 method was used to calculate electronic energy of reactants, products, and TS; frequencies were
calculated using B3LYP/CBSB7 method.  Arrhenius expression was computed using CanTherm: an asymmetric Eckart
tunneling correction was employed and the frequencies were scaled by 0.99 (as suggested by Montgomery et al.
J.Chem.Phys. 110 (1999) 2822-2827).  The external symmetry number for H2O2 was 2; the external symmetry number
for the remaining species and TS were set to 1.  The rate coefficient was computed at 600-2000K (in 200 K increments).
The computed pre-exponential factor was divided by 2 and this is the reported value.

For nButanol+HO2=H2O2+CH3CH2CH2*CHOH:
Moc et al. (AIP Conference Proceedings (2009) 1148 161-164 \"The Unimolecular Decomposition
and H Abstraction Reactions by HO and HO2 from n-Butanol\") report reaction barriers and
enthalpies(0 K); our CBS-QB3 calculations are shown in comparison (all units are kcal/mol).
				G3		CCSD(T)/cc-pVTZ		CBS-QB3
Barrier:		12.62		13.23			11.74
Enthalpy:		 8.35		 8.63			 7.17
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 273
rate(
  group1 = 
"""
H2O2
1 *1 O 0 {2,S} {3,S}
2 O 0 {1,S} {4,S}
3 *2 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
InChI=1/C4H9O/c1-3-4(2)5/h4-5H,1,3H2,2H3
1 *3 C 1 {3,S} {6,S} {7,S}
2 C 0 {4,S} {8,S} {9,S} {10,S}
3 C 0 {1,S} {4,S} {11,S} {12,S}
4 C 0 {2,S} {3,S} {5,S} {13,S}
5 O 0 {4,S} {14,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {2,S}
9 H 0 {2,S}
10 H 0 {2,S}
11 H 0 {3,S}
12 H 0 {3,S}
13 H 0 {4,S}
14 H 0 {5,S}
""",
  kf = Arrhenius(A=(5.75E+00,A_UNITS,"+-",0.0),
                 n=(2.94,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(0.46,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "305",
  short_comment = "MRH CBS-QB3 calculations w/o HR corrections",
  long_comment = 
"""
MRH CBS-QB3 calculations w/o HR corrections
H2O2 + *CH2CH2CH[OH]CH3 = 2-Butanol + HO2

CBS-QB3 method was used to calculate electronic energy of reactants, products, and TS; frequencies were
calculated using B3LYP/CBSB7 method.  Arrhenius expression was computed using CanTherm: an asymmetric Eckart
tunneling correction was employed and the frequencies were scaled by 0.99 (as suggested by Montgomery et al.
J.Chem.Phys. 110 (1999) 2822-2827).  The external symmetry number for H2O2 was 2; the external symmetry number
for the remaining species and TS were set to 1.  The rate coefficient was computed at 600-2000K (in 200 K increments).
The computed pre-exponential factor was divided by 2 and this is the reported value.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 274
rate(
  group1 = 
"""
H2O2
1 *1 O 0 {2,S} {3,S}
2 O 0 {1,S} {4,S}
3 *2 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
InChI=1/C4H9O/c1-3-4(2)5/h3-5H,1-2H3
1 C 0 {3,S} {6,S} {7,S} {8,S}
2 C 0 {4,S} {9,S} {10,S} {11,S}
3 *3 C 1 {1,S} {4,S} {12,S}
4 C 0 {2,S} {3,S} {5,S} {13,S}
5 O 0 {4,S} {14,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {1,S}
9 H 0 {2,S}
10 H 0 {2,S}
11 H 0 {2,S}
12 H 0 {3,S}
13 H 0 {4,S}
14 H 0 {5,S}
""",
  kf = Arrhenius(A=(8.75E-01,A_UNITS,"+-",0.0),
                 n=(2.91,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(-0.41,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "306",
  short_comment = "MRH CBS-QB3 calculations w/o HR corrections",
  long_comment = 
"""
MRH CBS-QB3 calculations w/o HR corrections
H2O2 + CH3*CHCH[OH]CH3 = 2-Butanol + HO2

CBS-QB3 method was used to calculate electronic energy of reactants, products, and TS; frequencies were
calculated using B3LYP/CBSB7 method.  Arrhenius expression was computed using CanTherm: an asymmetric Eckart
tunneling correction was employed and the frequencies were scaled by 0.99 (as suggested by Montgomery et al.
J.Chem.Phys. 110 (1999) 2822-2827).  The external symmetry number for H2O2 was 2; the external symmetry number
for the remaining species and TS were set to 1.  The rate coefficient was computed at 600-2000K (in 200 K increments).
The computed pre-exponential factor was divided by 2 and this is the reported value.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 275
rate(
  group1 = 
"""
H2O2
1 *1 O 0 {2,S} {3,S}
2 O 0 {1,S} {4,S}
3 *2 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
InChI=1/C4H9O/c1-3-4(2)5/h5H,3H2,1-2H3
1 C 0 {3,S} {6,S} {7,S} {8,S}
2 C 0 {4,S} {9,S} {10,S} {11,S}
3 C 0 {1,S} {4,S} {12,S} {13,S}
4 *3 C 1 {2,S} {3,S} {5,S}
5 O 0 {4,S} {14,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {1,S}
9 H 0 {2,S}
10 H 0 {2,S}
11 H 0 {2,S}
12 H 0 {3,S}
13 H 0 {3,S}
14 H 0 {5,S}
""",
  kf = Arrhenius(A=(1.73E+01,A_UNITS,"+-",0.0),
                 n=(3.05,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(1.02,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "307",
  short_comment = "MRH CBS-QB3 calculations w/o HR corrections",
  long_comment = 
"""
MRH CBS-QB3 calculations w/o HR corrections
H2O2 + CH3CH2*C[OH]CH3 = 2-Butanol + HO2

CBS-QB3 method was used to calculate electronic energy of reactants, products, and TS; frequencies were
calculated using B3LYP/CBSB7 method.  Arrhenius expression was computed using CanTherm: an asymmetric Eckart
tunneling correction was employed and the frequencies were scaled by 0.99 (as suggested by Montgomery et al.
J.Chem.Phys. 110 (1999) 2822-2827).  The external symmetry number for H2O2 was 2; the external symmetry number
for the remaining species and TS were set to 1.  The rate coefficient was computed at 600-2000K (in 200 K increments).
The computed pre-exponential factor was divided by 2 and this is the reported value.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 276
rate(
  group1 = 
"""
H2O2
1 *1 O 0 {2,S} {3,S}
2 O 0 {1,S} {4,S}
3 *2 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
InChI=1/C4H9O/c1-3-4(2)5/h4-5H,2-3H2,1H3
1 C 0 {3,S} {6,S} {7,S} {8,S}
2 *3 C 1 {4,S} {9,S} {10,S}
3 C 0 {1,S} {4,S} {11,S} {12,S}
4 C 0 {2,S} {3,S} {5,S} {13,S}
5 O 0 {4,S} {14,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {1,S}
9 H 0 {2,S}
10 H 0 {2,S}
11 H 0 {3,S}
12 H 0 {3,S}
13 H 0 {4,S}
14 H 0 {5,S}
""",
  kf = Arrhenius(A=(3.055E-01,A_UNITS,"+-",0.0),
                 n=(3.53,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(1.52,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "308",
  short_comment = "MRH CBS-QB3 calculations w/o HR corrections",
  long_comment = 
"""
MRH CBS-QB3 calculations w/o HR corrections
H2O2 + CH3CH2CH[OH]*CH2 = 2-Butanol + HO2

CBS-QB3 method was used to calculate electronic energy of reactants, products, and TS; frequencies were
calculated using B3LYP/CBSB7 method.  Arrhenius expression was computed using CanTherm: an asymmetric Eckart
tunneling correction was employed and the frequencies were scaled by 0.99 (as suggested by Montgomery et al.
J.Chem.Phys. 110 (1999) 2822-2827).  The external symmetry number for H2O2 was 2; the external symmetry number
for the remaining species and TS were set to 1.  The rate coefficient was computed at 600-2000K (in 200 K increments).
The computed pre-exponential factor was divided by 2 and this is the reported value.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 277
rate(
  group1 = 
"""
H2O2
1 *1 O 0 {2,S} {3,S}
2 O 0 {1,S} {4,S}
3 *2 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
InChI=1/C4H9O/c1-4(2,3)5/h5H,1H2,2-3H3
1 *3 C 1 {4,S} {6,S} {7,S}
2 C 0 {4,S} {8,S} {9,S} {10,S}
3 C 0 {4,S} {11,S} {12,S} {13,S}
4 C 0 {1,S} {2,S} {3,S} {5,S}
5 O 0 {4,S} {14,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {2,S}
9 H 0 {2,S}
10 H 0 {2,S}
11 H 0 {3,S}
12 H 0 {3,S}
13 H 0 {3,S}
14 H 0 {5,S}
""",
  kf = Arrhenius(A=(2.100E-01,A_UNITS,"+-",0.0),
                 n=(3.53,None,"+-",0.0),
                 alpha=(0.0,None,"+-",0.0),
                 E0=(1.56,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "309",
  short_comment = "MRH CBS-QB3 calculations w/o HR corrections",
  long_comment = 
"""
MRH CBS-QB3 calculations w/o HR corrections
H2O2 + HOC[*CH2][CH3][CH3] = tert-Butanol + HO2

CBS-QB3 method was used to calculate electronic energy of reactants, products, and TS; frequencies were
calculated using B3LYP/CBSB7 method.  Arrhenius expression was computed using CanTherm: an asymmetric Eckart
tunneling correction was employed and the frequencies were scaled by 0.99 (as suggested by Montgomery et al.
J.Chem.Phys. 110 (1999) 2822-2827).  The external symmetry number for H2O2 was 2; the external symmetry number
for the remaining species and TS were set to 1.  The rate coefficient was computed at 600-2000K (in 200 K increments).
The computed pre-exponential factor was divided by 2 and this is the reported value.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 278
rate(
  group1 = 
"""
Orad_O_H
1 *1 O 0 {2,S}, {3,S}
2 *2 H 0 {1,S}
3    O 1 {1,S}
""",
  group2 = 
"""
O_rad/NonDeO
1 *3 O 1 {2,S}
2    O 0 {1,S}
""",
  kf = Arrhenius(A=(1.75E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-3.275,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "486",
  short_comment = "[8] Curran\'s estimation in reaction type 13, RO2 + HO2 --> RO2H + O2",
  long_comment = 
"""
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)


