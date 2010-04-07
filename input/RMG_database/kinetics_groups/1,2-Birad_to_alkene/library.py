# encoding: utf-8
header = """
1,2-Birad_to_alkene

Y_12birad -> Y_alkene


Reverse name: Alkene_to_1,2-birad

(1) CHANGE_BOND		{*1,1,*2}
(2) LOSE_RADICAL 	{*1,1}
(3) LOSE_RADICAL	{*2,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "1,2-Birad_to_alkene"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
//f24 : 1,2-biradical to alkene 

Arrhenius_EP
//reaction family added by gmagoon 3/1/10
//correlation used for #1-16: log10(tau(s)) = -8.0+0.2*m+0.3*n; m = number of alkyl substituents, n = number of aryl/vinyl substituents
//references for this correlation and for treatment of alkene triplets as 1,2-biradicals:
//1.	R. A. Caldwell, Intersystem crossing in organic photochemical intermediates, Pure Appl. Chem., 56 (1984) 1167-1177.
//2.	R. A. Caldwell, Laser flash photolysis studies of intersystem crossing in biradicals and alkene triplets, p. 110, in Kinetics and spectroscopy of carbenes and biradicals, ed. M. S. Platz, 1990.
//3.	D. J. Unett and R. A. Caldwell, The triplet state of alkenes: structure, dynamics, energetics and chemistry, Res. Chem. Intermed., 21 (1995) 665-709
//4.	T. Ni, R. A. Caldwell, and L. A. Melton, The relaxed and spectroscopic energies of olefin triplets, J. Am. Chem. Soc., 111 (1989) 457-464.
//to extrapolate to other groups, I am basing my (rough) assignments on the author\'s argument in Ref. 1, that effect is mainly related to number of hydrogens and mass of substituents, rather than electronic stabilization/polar effects; note however, that this will not correctly account for large mass of extended groups like large alkyl chains (in any case, the effect is relatively small, and the resulting estimate should still be within an order or magnitude or so of what we would obtain if we assigned groups differently)
//the assignments I use are:  no slow-down: H
//					m slow-down: Cs, Os
//					n slow-down: Cd, Ct, Cb, CO
//it is assumed that k = 1/tau(s)
//nomenclature used is Y_12_mn where m and n are defined as above

//No.	Y_12birad		Temp.		A		n		a		E0		DA		Dn		Da		DE0		Rank	Comments


"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
Y_12birad
Union {Y_12_00,Y_12_10,Y_12_20,Y_12_30,Y_12_40,Y_12_01,Y_12_02,Y_12_03,Y_12_04,Y_12_11,Y_12_12,Y_12_21,Y_12_22,Y_12_13,Y_12_31}
""",
  kf = Arrhenius(A=(1.0E+8,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
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

# Number 2
rate(
  group1 = 
"""
Y_12_00
1 *1 Cs 1 {2,S}, {3,S}, {4,S}
2 *2 Cs 1 {1,S}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.0E+8,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "2",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
Y_12_10
1 *1 Cs 1 {2,S}, {3,S}, {4,S}
2 *2 Cs 1 {1,S}, {5,S}, {6,S}
3 {Cs,Os} 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(6.31E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "3",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
Y_12_20
Union {Y_12_20a, Y_12_20b}
""",
  kf = Arrhenius(A=(3.98E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "4",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
Y_12_30
1 *1 Cs 1 {2,S}, {3,S}, {4,S}
2 *2 Cs 1 {1,S}, {5,S}, {6,S}
3 {Cs,Os} 0 {1,S}
4 {Cs,Os} 0 {1,S}
5 {Cs,Os} 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(2.51E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "5",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
Y_12_40
1 *1 Cs 1 {2,S}, {3,S}, {4,S}
2 *2 Cs 1 {1,S}, {5,S}, {6,S}
3 {Cs,Os} 0 {1,S}
4 {Cs,Os} 0 {1,S}
5 {Cs,Os} 0 {2,S}
6 {Cs,Os} 0 {2,S}
""",
  kf = Arrhenius(A=(1.58E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "6",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
Y_12_01
1 *1 Cs 1 {2,S}, {3,S}, {4,S}
2 *2 Cs 1 {1,S}, {5,S}, {6,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(5.01E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "7",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
Y_12_02
Union {Y_12_02a, Y_12_02b}
""",
  kf = Arrhenius(A=(2.51E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "8",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 9
rate(
  group1 = 
"""
Y_12_03
1 *1 Cs 1 {2,S}, {3,S}, {4,S}
2 *2 Cs 1 {1,S}, {5,S}, {6,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 {Cd,Ct,Cb,CO} 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.26E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "9",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 10
rate(
  group1 = 
"""
Y_12_04
1 *1 Cs 1 {2,S}, {3,S}, {4,S}
2 *2 Cs 1 {1,S}, {5,S}, {6,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 {Cd,Ct,Cb,CO} 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  kf = Arrhenius(A=(6.31E+6,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "10",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 11
rate(
  group1 = 
"""
Y_12_11
Union {Y_12_11a, Y_12_11b}
""",
  kf = Arrhenius(A=(3.16E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "11",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 12
rate(
  group1 = 
"""
Y_12_12
Union {Y_12_12a, Y_12_12b}
""",
  kf = Arrhenius(A=(1.58E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "12",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 13
rate(
  group1 = 
"""
Y_12_21
Union {Y_12_21a, Y_12_21b}
""",
  kf = Arrhenius(A=(2.00E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "13",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 14
rate(
  group1 = 
"""
Y_12_22
Union {Y_12_22a, Y_12_22b}
""",
  kf = Arrhenius(A=(1.0E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "14",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 15
rate(
  group1 = 
"""
Y_12_13
1 *1 Cs 1 {2,S}, {3,S}, {4,S}
2 *2 Cs 1 {1,S}, {5,S}, {6,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 {Cd,Ct,Cb,CO} 0 {2,S}
6 {Cs,Os} 0 {2,S}
""",
  kf = Arrhenius(A=(7.94E+6,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "15",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 16
rate(
  group1 = 
"""
Y_12_31
1 *1 Cs 1 {2,S}, {3,S}, {4,S}
2 *2 Cs 1 {1,S}, {5,S}, {6,S}
3 {Cs,Os} 0 {1,S}
4 {Cs,Os} 0 {1,S}
5 {Cs,Os} 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  kf = Arrhenius(A=(1.26E+7,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "16",
  short_comment = "see references header of 1,2-Birad_to_alkene/rateLibrary.txt",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


