# encoding: utf-8
header = """
Diels_alder_addition

ene + diene_out -> Six_Ring


Reverse name: Retro_Diels_Alder_Addition

(1) CHANGE_BOND		{*1,-1,*2}
(2) CHANGE_BOND		{*3,-1,*4}
(3) CHANGE_BOND		{*4,1,*5}
(4) CHANGE_BOND		{*5,-1,*6}
(5) FORM_BOND		{*1,S,*3}
(6) FORM_BOND		{*2,S,*6}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "Diels_alder_addition"

# determines permitted units for rate expression:
reaction_order = 2

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f19: diels-Alder reaction

// JS, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// f19_diels_Alde
// rate constants from rate_library_4.txt, Cath, 03/07/28

// Catherina Wijaya thesis, pg 132

// [108] Huybrechts, G.; Poppelsdorf, H.; Maesschalck, L.; Van Mele, B. Int. J. Chem. Kinet. 1984, 16, 93.
// [109] Huybrechts, G.; Rigaux, D.; Vankeerberghen, J.; Van Mele, B. Int. J. Chem. Kinet. 1980, 12, 253.
// [110] Van Mele, B.; Tybaert, C.; Huybrechts, G.  Int. J. Chem. Kinet. 1987, 19, 1063.
// [111] Huybrechts, G.;Hubin, Y.; Narmon, M.; Van Mele, B. Int. J. Chem. Kinet. 1982, 14, 259.
// [112] Kistiakowsky, G. B.; Lacher, J. R. J. Am. Chem. Soc. 1936, 58, 123.
// [198] Huybrechts, G.; Luyckx, L.; Vandenboom, T.; Van Mele, B. Int. J. Chem. Kinet. 1977, 9, 283.
// [199] Simmie, J. M. Int. J. Chem. Kinet. 1978, 10, 227.
// [200] Benford, G. A.; Wassermann, A. J. Chem. Soc. 1939, 362. 

//No.	diene_out						diene_in			ene				Temp.		A			N		a		E0		DA		Dn		Da		DE0		Rank	Comments

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
diene_out
Union {diene_unsub_unsub_out,diene_unsub_monosub_out,diene_unsub_disub_out,diene_monosub_monosub_out,diene_monosub_disub_out,diene_disub_disub_out}
""",
  group2 = 
"""
diene_in
1 *4 Cd 0 {2,S}
2 *5 Cd 0 {1,S}
""",
  group3 = 
"""
ene
1 *1 Cd 0 {2,D}
2 *2 Cd 0 {1,D}
""",
  kf = Arrhenius(A=(5E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "589",
  short_comment = "\"default\"",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 2
rate(
  group1 = 
"""
diene_unsub_unsub_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HDe
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  kf = Arrhenius(A=(8.91E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(24.44,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (464,557),
  rank = 3,
  old_id = "590",
  short_comment = "Huybrechts et al. [198]",
  long_comment = 
"""
[198] Huybrechts, G.; Luyckx, L.; Vandenboom, T.; Van Mele, B. Int. J. Chem. Kinet. 1977, 9, 283.
(E)-CH2=CHCH=CH2 + (E)-CH2=CHCH=CH2 --> 4-vinylcyclohexene

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.06-0.59 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
diene_unsub_unsub_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_HDe_2H
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(8.91E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(24.44,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (464,557),
  rank = 3,
  old_id = "591",
  short_comment = "Huybrechts et al. [198]",
  long_comment = 
"""
[198] Huybrechts, G.; Luyckx, L.; Vandenboom, T.; Van Mele, B. Int. J. Chem. Kinet. 1977, 9, 283.
(E)-CH2=CHCH=CH2 + (E)-CH2=CHCH=CH2 --> 4-vinylcyclohexene

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.06-0.59 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
diene_unsub_unsub_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_HNd_HDe
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cs,O} 0 {1,S}
5 H 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  kf = Arrhenius(A=(8.99E+08,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(22.06,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (515,572),
  rank = 3,
  old_id = "592",
  short_comment = "Kistiakowsky et al [112]",
  long_comment = 
"""
[112] Kistiakowsky, G. B.; Lacher, J. R. J. Am. Chem. Soc. 1936, 58, 123.
(Z)-CH3CH=CHCHO + (E)-CH2=CHCH=CH2 --> 3-cyclohexene-1-carboxaldehyde,6-methyl

Absolute value measured directly. Excitation: thermal. Pressure 0.64-0.78 atm	
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
diene_unsub_unsub_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_HDe_HNd
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(8.99E+08,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(22.06,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (515,572),
  rank = 3,
  old_id = "593",
  short_comment = "Kistiakowsky et al [112]",
  long_comment = 
"""
[112] Kistiakowsky, G. B.; Lacher, J. R. J. Am. Chem. Soc. 1936, 58, 123.
(Z)-CH3CH=CHCHO + (E)-CH2=CHCH=CH2 --> 3-cyclohexene-1-carboxaldehyde,6-methyl

Absolute value measured directly. Excitation: thermal. Pressure 0.64-0.78 atm	
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
diene_unsub_unsub_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_HNd
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 {Cs,O} 0 {2,S}
""",
  group3 = 
"""
ene_unsub_unsub
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.32E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.61,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (1000,1180),
  rank = 4,
  old_id = "594",
  short_comment = "Simmie [199]",
  long_comment = 
"""
[199] Simmie, J. M. Int. J. Chem. Kinet. 1978, 10, 227.
CH2=C(CH3)CH=CH2 + C2H4 --> 1-methyl-cyclohexane

Data derived from fitting to a complex mechanism. Excitation: thermal. Analysis: GC. Presure 2.50 atm
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
diene_unsub_unsub_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_NdH
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 {Cs,O} 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_unsub_unsub
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.32E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.61,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (1000,1180),
  rank = 4,
  old_id = "595",
  short_comment = "Simmie [199]",
  long_comment = 
"""
[199] Simmie, J. M. Int. J. Chem. Kinet. 1978, 10, 227.
CH2=C(CH3)CH=CH2 + C2H4 --> 1-methyl-cyclohexane

Data derived from fitting to a complex mechanism. Excitation: thermal. Analysis: GC. Presure 2.50 atm
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
diene_unsub_unsub_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_HNd
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 {Cs,O} 0 {2,S}
""",
  group3 = 
"""
ene_HDe_2H
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.02E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (492,606),
  rank = 3,
  old_id = "596",
  short_comment = "Kistiakowsky et al [112]",
  long_comment = 
"""
[112] Kistiakowsky, G. B.; Lacher, J. R. J. Am. Chem. Soc. 1936, 58, 123.
CH2=CHCHO + CH2=C(CH3)CH=CH2 --> 3-cyclohexene-1-carboxaldehyde,4-methyl

Absolute value measured directly. Excitation: thermal. Pressure 0.64-0.78 atm
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 9
rate(
  group1 = 
"""
diene_unsub_unsub_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_NdH
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 {Cs,O} 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HDe
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  kf = Arrhenius(A=(1.02E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (492,606),
  rank = 3,
  old_id = "597",
  short_comment = "Kistiakowsky et al [112]",
  long_comment = 
"""
[112] Kistiakowsky, G. B.; Lacher, J. R. J. Am. Chem. Soc. 1936, 58, 123.
CH2=CHCHO + CH2=C(CH3)CH=CH2 --> 3-cyclohexene-1-carboxaldehyde,4-methyl

Absolute value measured directly. Excitation: thermal. Pressure 0.64-0.78 atm
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 10
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_unsub_unsub
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(4.57E+09,A_UNITS,"*/",1.05),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.03,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (450,592),
  rank = 3,
  old_id = "598",
  short_comment = "Huybrechts et al. [109]",
  long_comment = 
"""
[109] Huybrechts, G.; Rigaux, D.; Vankeerberghen, J.; Van Mele, B. Int. J. Chem. Kinet. 1980, 12, 253.
1,3-cyclohexadiene + C2H4 --> bicyclo[2.2.2]oct-2-ene

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.05-0.25 atm.	
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 11
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HNd
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(1.12E+09,A_UNITS,"*/",1.12),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.63,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (488,606),
  rank = 3,
  old_id = "599",
  short_comment = "Huybrechts et al. [108]",
  long_comment = 
"""
[108] Huybrechts, G.; Poppelsdorf, H.; Maesschalck, L.; Van Mele, B. Int. J. Chem. Kinet. 1984, 16, 93.
1,3-cyclohexadiene + CH3CH=CH2 --> bicyclo[2.2.2]oct-2-ene,5-METHYL-(1alpha, 4alpha, 5beta)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.07-0.82 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 12
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HNd
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(2.09E+09,A_UNITS,"*/",1.12),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(28.81,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (488,606),
  rank = 3,
  old_id = "600",
  short_comment = "Huybrechts et al. [108]",
  long_comment = 
"""
[108] Huybrechts, G.; Poppelsdorf, H.; Maesschalck, L.; Van Mele, B. Int. J. Chem. Kinet. 1984, 16, 93.
1,3-cyclohexadiene + CH3CH=CH2 --> bicyclo[2.2.2]oct-2-ene,5-METHYL-(1alpha, 4alpha, 5alpha)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.07-0.82 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 13
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HNd
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(7.08E+08,A_UNITS,"*/",1.12),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.23,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (488,606),
  rank = 3,
  old_id = "601",
  short_comment = "Huybrechts et al. [108]",
  long_comment = 
"""
[108] Huybrechts, G.; Poppelsdorf, H.; Maesschalck, L.; Van Mele, B. Int. J. Chem. Kinet. 1984, 16, 93.
1,3-cyclohexadiene + 1-C4H8 --> bicyclo[2.2.2]oct-2-ene,5-ETHYL-(1alpha, 4alpha, 5beta)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.07-0.82 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 14
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HNd
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(1.17E+09,A_UNITS,"*/",1.12),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(28.62,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (488,606),
  rank = 3,
  old_id = "602",
  short_comment = "Huybrechts et al. [108]",
  long_comment = 
"""
[108] Huybrechts, G.; Poppelsdorf, H.; Maesschalck, L.; Van Mele, B. Int. J. Chem. Kinet. 1984, 16, 93.
1,3-cyclohexadiene + 1-C4H8 --> bicyclo[2.2.2]oct-2-ene,5-ETHYL-(1alpha, 4alpha, 5alpha)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.07-0.82 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 15
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HNd
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(3.72E+08,A_UNITS,"*/",1.07),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.63,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (488,600),
  rank = 3,
  old_id = "603",
  short_comment = "Huybrechts et al. [108]",
  long_comment = 
"""
[108] Huybrechts, G.; Poppelsdorf, H.; Maesschalck, L.; Van Mele, B. Int. J. Chem. Kinet. 1984, 16, 93.
1,3-cyclohexadiene + (CH3)2CHCH=CH2 --> bicyclo[2.2.2]oct-2-ene,5-(1-METHYLETHYL)-(1alpha, 4alpha, 5beta)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.07-0.82 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 16
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HNd
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(2.95E+08,A_UNITS,"*/",1.1),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(28.42,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (486,600),
  rank = 3,
  old_id = "604",
  short_comment = "Huybrechts et al. [108]",
  long_comment = 
"""
[108] Huybrechts, G.; Poppelsdorf, H.; Maesschalck, L.; Van Mele, B. Int. J. Chem. Kinet. 1984, 16, 93.
1,3-cyclohexadiene + (CH3)2CHCH=CH2 --> bicyclo[2.2.2]oct-2-ene,5-(1-METHYLETHYL)-(1alpha, 4alpha, 5alpha)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.07-0.82 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 17
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_HNd_2H
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cs,O} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.12E+09,A_UNITS,"*/",1.12),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.63,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (488,606),
  rank = 3,
  old_id = "605",
  short_comment = "Huybrechts et al. [108]",
  long_comment = 
"""
[108] Huybrechts, G.; Poppelsdorf, H.; Maesschalck, L.; Van Mele, B. Int. J. Chem. Kinet. 1984, 16, 93.
1,3-cyclohexadiene + CH3CH=CH2 --> bicyclo[2.2.2]oct-2-ene,5-METHYL-(1alpha, 4alpha, 5beta)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.07-0.82 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 18
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HDe
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  kf = Arrhenius(A=(1.02E+09,A_UNITS,"*/",1.07),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.07,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (379,581),
  rank = 3,
  old_id = "606",
  short_comment = "Van Mele et al [110]",
  long_comment = 
"""
[110] Van Mele, B.; Tybaert, C.; Huybrechts, G.  Int. J. Chem. Kinet. 1987, 19, 1063.
1,3-cyclohexadiene + CH2=CHCHO --> bicyclo[2.2.2]oct-2-ene,2-carboxaldehyde(1alpha, 2alpha, 4alpha)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.06-0.27 atm.	
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 19
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HDe
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  kf = Arrhenius(A=(6.03E+08,A_UNITS,"*/",1.07),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.87,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (379,581),
  rank = 3,
  old_id = "607",
  short_comment = "Van Mele et al [110]",
  long_comment = 
"""
[110] Van Mele, B.; Tybaert, C.; Huybrechts, G.  Int. J. Chem. Kinet. 1987, 19, 1063.
1,3-cyclohexadiene + CH2=CHCHO --> bicyclo[2.2.2]oct-2-ene,2-carboxaldehyde(1alpha, 2beta, 4alpha)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.06-0.27 atm.	
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 20
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HDe
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  kf = Arrhenius(A=(1.15E+10,A_UNITS,"*/",1.05),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.83,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (437,526),
  rank = 3,
  old_id = "608",
  short_comment = "Huybrechts et al. [111]",
  long_comment = 
"""
[111] Huybrechts, G.;Hubin, Y.; Narmon, M.; Van Mele, B. Int. J. Chem. Kinet. 1982, 14, 259.
1,3-cyclohexadiene + (E)CH2=CHCH=CH2 --> bicyclo[2.2.2]oct-2-ene,5-ethenyl(1alpha, 4alpha, 5alpha)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.15-0.64 atm.	
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 21
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_2H_HDe
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  kf = Arrhenius(A=(3.8E+09,A_UNITS,"*/",1.05),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(24.84,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (437,526),
  rank = 3,
  old_id = "609",
  short_comment = "Huybrechts et al. [111]",
  long_comment = 
"""
[111] Huybrechts, G.;Hubin, Y.; Narmon, M.; Van Mele, B. Int. J. Chem. Kinet. 1982, 14, 259.
1,3-cyclohexadiene + (E)CH2=CHCH=CH2 --> bicyclo[2.2.2]oct-2-ene,5-ethenyl(1alpha, 4alpha, 5beta)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.15-0.64 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 22
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_HDe_2H
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.02E+09,A_UNITS,"*/",1.07),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.07,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (379,581),
  rank = 3,
  old_id = "610",
  short_comment = "Van Mele et al [110]",
  long_comment = 
"""
[110] Van Mele, B.; Tybaert, C.; Huybrechts, G.  Int. J. Chem. Kinet. 1987, 19, 1063.
1,3-cyclohexadiene + CH2=CHCHO --> bicyclo[2.2.2]oct-2-ene,2-carboxaldehyde(1alpha, 2alpha, 4alpha)

Absolute value measured directly using thermal excitation technique and GC. Pressure 0.06-0.27 atm.	
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 23
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_HNd_HDe
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cs,O} 0 {1,S}
5 H 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  kf = Arrhenius(A=(1.26E+9,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.69,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (352,423),
  rank = 3,
  old_id = "611",
  short_comment = "Benford et al [200]",
  long_comment = 
"""
[200] Benford, G. A.; Wassermann, A. J. Chem. Soc. 1939, 362. 
Cyclopentadiene + cyclopentadiene --> Tricyclo[5.2.1.02,6]deca-c,8-diene.

Absolute value measured directly using thermal excitation technique and mass spectrometry. Pressure 0.20-0.97 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 24
rate(
  group1 = 
"""
diene_monosubNd_monosubNd_out
1 *3 Cd 0 {2,D}, {5,S}, {6,S}
2 *4 Cd 0 {1,D}, {3,S}
3 *5 Cd 0 {2,S}, {4,D}
4 *6 Cd 0 {3,D}, {7,S}, {8,S}
5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {4,S}
8 H 0 {4,S}
""",
  group2 = 
"""
diene_in_2H
1 *4 Cd 0 {2,S} {3,S}
2 *5 Cd 0 {1,S} {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group3 = 
"""
ene_HDe_HNd
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(1.26E+9,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.69,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (352,423),
  rank = 3,
  old_id = "612",
  short_comment = "Benford et al [200]",
  long_comment = 
"""
[200] Benford, G. A.; Wassermann, A. J. Chem. Soc. 1939, 362. 
Cyclopentadiene + cyclopentadiene --> Tricyclo[5.2.1.02,6]deca-c,8-diene.

Absolute value measured directly using thermal excitation technique and mass spectrometry. Pressure 0.20-0.97 atm.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


