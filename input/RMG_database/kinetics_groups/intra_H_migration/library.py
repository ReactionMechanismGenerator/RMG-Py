# encoding: utf-8
header = """
intra_H_migration

RnH -> HR.


(1) BREAK_BOND		{*2,S,*3}
(2) FORM_BOND		{*1,S,*3}
(3) GAIN_RADICAL	{*2,1}
(4) LOSE_RADICAL 	{*1,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "intra_H_migration"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f33: intra hydroxyl migration
// Originally from rate library.doc by CDW

// \"jing,\" define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// #690 has a conflict of the name and the catogery in sumathy\'s original excel \"file,\" should ask her about it. Cath and \"Jing,\" 3/7/2028
// make change of \"#673,\" \"#675,\" #676 according to sumathy\'s \"correction,\" \"JS,\" 3/8/2026

// Catherine Wijaya Thesis, pg 133, 159.

//f25_intra_H_migration
//No	RnH				Y_rad_out				XH_out					Temp		A			n		Alpha	E		DA	DN	DAlpha	DE	Rank	Comments
//690.	R2H_S_cy4		C_rad_out_				Cs_H_out_				300-1500	2.38E+11	0.5		0		40.60	0	0	0		0	2		Sumathy B3LYP/CCPVDZ calculations
//From 850 to 869 added by Sandeep. The rate rules are from DFT/CBSB7 level of calculations.
//// rwest 2008/10/22: these commented out because functional group R6H_SSSSS_(Cs/Cs)OO is not recognised
//867 R6H_SSSSS_(Cs/Cs)OO		O_rad_out		Cs_H_out_H/NonDeO		300-1500	1.6e4		1.89	0		14.38	0	0	0		0	3		Sandeep\'s DFT/CBSB7 level of calculations.
//867	R6H_SSSSS_(Cs/Cs)OO		O_rad_out		Cs_H_out_NDMustO		300-1500	1.03e5		1.26	0		12.02	0	0	0		0	3		Sandeep\'s DFT/CBSB7 level of calculations.
//From 870 to 875 added by sandeep based on the work done on 1,3-hexadiene chemistry cbs-qb3

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
RnH
Union {R2H, R3H, R4H, R5H, R6H, R7H}
""",
  group2 = 
"""
Y_rad_out
1 *1 {R!H} 1
""",
  group3 = 
"""
XH_out
1 *2 {R!H} 0 {2,S}
2 *3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.00E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(25.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "614",
  short_comment = "default",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 2
rate(
  group1 = 
"""
R2H_S
1 *1 {R!H} 1 {2,S}
2 *2 {R!H} 0 {1,S}, {3,S}
3 *3 H 0 {2,S}
""",
  group2 = 
"""
C_rad_out_single
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.48E+08,A_UNITS,"+-",0.0),
                 n=(1.62,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.76,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "615",
  short_comment = "Currans\'s estimation [5] in his reaction type 5.",
  long_comment = 
"""
[5] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 1998, 114, 149. 
Currans\'s estimation in his reaction type 5. C7H15

Checked by Paul Green. 
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
R2H_S
1 *1 {R!H} 1 {2,S}
2 *2 {R!H} 0 {1,S}, {3,S}
3 *3 H 0 {2,S}
""",
  group2 = 
"""
C_rad_out_single
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_1H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.59E+08,A_UNITS,"+-",0.0),
                 n=(1.39,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "616",
  short_comment = "Currans\'s estimation [5] in his reaction type 5.",
  long_comment = 
"""
[5] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 1998, 114, 149. 
Currans\'s estimation in his reaction type 5. C7H15

Checked by Paul Green.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
R3H_SS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_single
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.39E+09,A_UNITS,"+-",0.0),
                 n=(0.98,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.76,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "617",
  short_comment = "Currans\'s estimation [5] in his reaction type 5.",
  long_comment = 
"""
[5] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 1998, 114, 149. 
Currans\'s estimation in his reaction type 5. C7H15

Checked By Paul Green
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
R3H_SS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_single
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_1H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.76E+09,A_UNITS,"+-",0.0),
                 n=(0.76,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "618",
  short_comment = "Currans\'s estimation [5] in his reaction type 5.",
  long_comment = 
"""
[5] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 1998, 114, 149. 
Currans\'s estimation in his reaction type 5. C7H15

Checked By Paul Green.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
C_rad_out_single
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.54E+09,A_UNITS,"+-",0.0),
                 n=(0.35,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.76,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "619",
  short_comment = "Currans\'s estimation [5] in his reaction type 5.",
  long_comment = 
"""
[5] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 1998, 114, 149. 
Currans\'s estimation in his reaction type 5. C7H15

Checked By Paul Green.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
C_rad_out_single
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_1H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.22E+09,A_UNITS,"+-",0.0),
                 n=(0.13,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "620",
  short_comment = "Currans\'s estimation [5] in his reaction type 5.",
  long_comment = 
"""
[5] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 1998, 114, 149. 
Currans\'s estimation in his reaction type 5. C7H15

Checked By Paul Green.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
C_rad_out_single
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_noH
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 {R!H} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.86E+10,A_UNITS,"+-",0.0),
                 n=(0.58,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.19,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "621",
  short_comment = "Currans\'s estimation [5] in his reaction type 5.",
  long_comment = 
"""
[5] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 1998, 114, 149. 
Currans\'s estimation in his reaction type 5. 

NEEDS TO BE CHECKED	
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 9
rate(
  group1 = 
"""
R5H_SSSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,S}
4 *5 {R!H} 0 {3,S}, {5,S}
5 *2 {R!H} 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
C_rad_out_single
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.28E+11,A_UNITS,"+-",0.0),
                 n=(-1.05,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.76,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "622",
  short_comment = "Currans\'s estimation [5] in his reaction type 5.",
  long_comment = 
"""
[5] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 1998, 114, 149. 
Currans\'s estimation in his reaction type 5. C7H15

Checked By Paul Green
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 10
rate(
  group1 = 
"""
R5H_SSSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,S}
4 *5 {R!H} 0 {3,S}, {5,S}
5 *2 {R!H} 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
C_rad_out_single
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_1H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.36E+10,A_UNITS,"+-",0.0),
                 n=(-0.66,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.28,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "623",
  short_comment = "Currans\'s estimation [5] in his reaction type 5.",
  long_comment = 
"""
[5] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 1998, 114, 149. 
Currans\'s estimation in his reaction type 5. C7H15

Checked by Paul Green
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 11
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.00E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "624",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 12
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_1H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.00E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.85,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "625",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 13
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_noH
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 {R!H} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.00E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(24.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "626",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 14
rate(
  group1 = 
"""
R5H_SSSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,S}
4 *5 {R!H} 0 {3,S}, {5,S}
5 *2 {R!H} 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.25E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(24.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "627",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 15
rate(
  group1 = 
"""
R5H_SSSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,S}
4 *5 {R!H} 0 {3,S}, {5,S}
5 *2 {R!H} 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_1H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.25E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.85,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "628",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 16
rate(
  group1 = 
"""
R5H_SSSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,S}
4 *5 {R!H} 0 {3,S}, {5,S}
5 *2 {R!H} 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_noH
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 {R!H} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.25E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "629",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 17
rate(
  group1 = 
"""
R6H_SSSSS
1 *1 {R!H} 1        {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3    {R!H} 0 {2,S}, {4,S}
4    {R!H} 0 {3,S}, {5,S}
5 *5 {R!H} 0 {4,S}, {6,S}
6 *2 {R!H} 0 {5,S}, {7,S}
7 *3  H    0 {6,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.56E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(22.35,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "630",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 18
rate(
  group1 = 
"""
R6H_SSSSS
1 *1 {R!H} 1        {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3    {R!H} 0 {2,S}, {4,S}
4    {R!H} 0 {3,S}, {5,S}
5 *5 {R!H} 0 {4,S}, {6,S}
6 *2 {R!H} 0 {5,S}, {7,S}
7 *3  H    0 {6,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_1H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.56E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.05,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "631",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 19
rate(
  group1 = 
"""
R6H_SSSSS
1 *1 {R!H} 1        {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3    {R!H} 0 {2,S}, {4,S}
4    {R!H} 0 {3,S}, {5,S}
5 *5 {R!H} 0 {4,S}, {6,S}
6 *2 {R!H} 0 {5,S}, {7,S}
7 *3  H    0 {6,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_noH
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 {R!H} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.56E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(17.05,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "632",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 20
rate(
  group1 = 
"""
R7H
1 *1 {R!H} 1 {2,{S,D,T,B}}
2 *4 {R!H} 0 {1,{S,D,T,B}}, {3,{S,D,T,B}}
3 {R!H} 0 {2,{S,D,T,B}}, {4,{S,D,T,B}}
4 {R!H} 0 {3,{S,D,T,B}}, {5,{S,D,T,B}}
5 {R!H} 0 {4,{S,D,T,B}}, {6,{S,D,T,B}}
6 *5 {R!H} 0 {5,{S,D,T,B}}, {7,{S,D,T,B}}
7 *2 {R!H} 0 {6,{S,D,T,B}}, {8,S}
8 *3 H 0 {7,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.95E+08,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(25.55,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "633",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 21
rate(
  group1 = 
"""
R7H
1 *1 {R!H} 1 {2,{S,D,T,B}}
2 *4 {R!H} 0 {1,{S,D,T,B}}, {3,{S,D,T,B}}
3 {R!H} 0 {2,{S,D,T,B}}, {4,{S,D,T,B}}
4 {R!H} 0 {3,{S,D,T,B}}, {5,{S,D,T,B}}
5 {R!H} 0 {4,{S,D,T,B}}, {6,{S,D,T,B}}
6 *5 {R!H} 0 {5,{S,D,T,B}}, {7,{S,D,T,B}}
7 *2 {R!H} 0 {6,{S,D,T,B}}, {8,S}
8 *3 H 0 {7,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_1H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.95E+08,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(25.55,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "634",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 22
rate(
  group1 = 
"""
R7H
1 *1 {R!H} 1 {2,{S,D,T,B}}
2 *4 {R!H} 0 {1,{S,D,T,B}}, {3,{S,D,T,B}}
3 {R!H} 0 {2,{S,D,T,B}}, {4,{S,D,T,B}}
4 {R!H} 0 {3,{S,D,T,B}}, {5,{S,D,T,B}}
5 {R!H} 0 {4,{S,D,T,B}}, {6,{S,D,T,B}}
6 *5 {R!H} 0 {5,{S,D,T,B}}, {7,{S,D,T,B}}
7 *2 {R!H} 0 {6,{S,D,T,B}}, {8,S}
8 *3 H 0 {7,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_noH
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 {R!H} 0 {1,S}
4 {R!H} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.95E+08,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(25.55,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "635",
  short_comment = "Curran\'s estimstion [8] in his reaction type 12 RO2 isomerization.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimstion in his reaction type 12 RO2 isomerization.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 23
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.45E+09,A_UNITS,"+-",0.0),
                 n=(1.12,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "636",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 24
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.10E+08,A_UNITS,"+-",0.0),
                 n=(1.32,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(40.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "637",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 25
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.69E+09,A_UNITS,"+-",0.0),
                 n=(0.89,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "638",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 26
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.12E+07,A_UNITS,"+-",0.0),
                 n=(1.66,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(40.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "639",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 27
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(4.04E+10,A_UNITS,"+-",0.0),
                 n=(0.64,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "640",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 28
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.28E+10,A_UNITS,"+-",0.0),
                 n=(0.97,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "641",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 29
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.38E+09,A_UNITS,"+-",0.0),
                 n=(0.88,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "642",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 30
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(7.25E+10,A_UNITS,"+-",0.0),
                 n=(0.6,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "643",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 31
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.12E+09,A_UNITS,"+-",0.0),
                 n=(1.19,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "644",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 32
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Cd_rad_out_double
1 *1 Cd 1 {2,D}
2 Cd 0 {1,D}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.44E+09,A_UNITS,"+-",0.0),
                 n=(1.12,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "645",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 33
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_doubleC
1 *2 Cd 0 {2,S} {3,D}
2 *3 H 0 {1,S}
3 Cd 0 {1,D}
""",
  kf = Arrhenius(A=(2.68E+11,A_UNITS,"+-",0.0),
                 n=(0.63,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(62.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "646",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 34
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Cd_rad_out_double
1 *1 Cd 1 {2,D}
2 Cd 0 {1,D}
""",
  group3 = 
"""
Cs_H_out_H/OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.24E+09,A_UNITS,"+-",0.0),
                 n=(0.82,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "647",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 35
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/OneDe
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_doubleC
1 *2 Cd 0 {2,S} {3,D}
2 *3 H 0 {1,S}
3 Cd 0 {1,D}
""",
  kf = Arrhenius(A=(9.38E+10,A_UNITS,"+-",0.0),
                 n=(0.71,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(62.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "648",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 36
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Cd_rad_out_double
1 *1 Cd 1 {2,D}
2 Cd 0 {1,D}
""",
  group3 = 
"""
Cs_H_out_OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 {Cs,O} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.67E+10,A_UNITS,"+-",0.0),
                 n=(0.79,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "649",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 37
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_OneDe/Cs
1 *1 C 1 {2,S}, {3,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_doubleC
1 *2 Cd 0 {2,S} {3,D}
2 *3 H 0 {1,S}
3 Cd 0 {1,D}
""",
  kf = Arrhenius(A=(1.03E+09,A_UNITS,"+-",0.0),
                 n=(1.31,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(62.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "650",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 38
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/OneDe
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.06E+09,A_UNITS,"+-",0.0),
                 n=(1.22,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(47.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "651",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 39
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.41E+08,A_UNITS,"+-",0.0),
                 n=(1.28,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(27.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "652",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 40
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/OneDe
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.45E+10,A_UNITS,"+-",0.0),
                 n=(0.75,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(45.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "653",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 41
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.41E+09,A_UNITS,"+-",0.0),
                 n=(0.35,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "654",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 42
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/OneDe
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.01E+12,A_UNITS,"+-",0.0),
                 n=(0.33,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(42.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "655",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 43
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.47E+08,A_UNITS,"+-",0.0),
                 n=(1.27,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "656",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 44
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_OneDe/Cs
1 *1 C 1 {2,S}, {3,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.69E+08,A_UNITS,"+-",0.0),
                 n=(1.31,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(48.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "657",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 45
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 {Cs,O} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(4.89E+09,A_UNITS,"+-",0.0),
                 n=(0.81,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(25.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "658",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 46
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_OneDe/Cs
1 *1 C 1 {2,S}, {3,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.13E+10,A_UNITS,"+-",0.0),
                 n=(0.77,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(46.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "659",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 47
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 {Cs,O} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(8.83E+10,A_UNITS,"+-",0.0),
                 n=(0.3,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "660",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 48
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_OneDe/Cs
1 *1 C 1 {2,S}, {3,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(3.62E+13,A_UNITS,"+-",0.0),
                 n=(-0.14,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(44.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "661",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 49
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 {Cs,O} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(8.20E+09,A_UNITS,"+-",0.0),
                 n=(0.65,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(31.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "662",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 50
rate(
  group1 = 
"""
R2H_D
1 *1 Cd 1 {2,D}
2 *2 Cd 0 {1,D}, {3,S}
3 *3 H 0 {2,S}
""",
  group2 = 
"""
Cd_rad_out_singleH
1 *1 Cd 1 {2,S}
2 H 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_singleH
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.28E+10,A_UNITS,"+-",0.0),
                 n=(0.86,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(45.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "663",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 51
rate(
  group1 = 
"""
R2H_D
1 *1 Cd 1 {2,D}
2 *2 Cd 0 {1,D}, {3,S}
3 *3 H 0 {2,S}
""",
  group2 = 
"""
Cd_rad_out_singleH
1 *1 Cd 1 {2,S}
2 H 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_singleNd
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(3.24E+11,A_UNITS,"+-",0.0),
                 n=(0.73,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(42.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "664",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 52
rate(
  group1 = 
"""
R2H_D
1 *1 Cd 1 {2,D}
2 *2 Cd 0 {1,D}, {3,S}
3 *3 H 0 {2,S}
""",
  group2 = 
"""
Cd_rad_out_singleNd
1 *1 Cd 1 {2,S}
2 {Cs,O} 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_singleH
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.62E+11,A_UNITS,"+-",0.0),
                 n=(0.8,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(47.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "665",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 53
rate(
  group1 = 
"""
R2H_D
1 *1 Cd 1 {2,D}
2 *2 Cd 0 {1,D}, {3,S}
3 *3 H 0 {2,S}
""",
  group2 = 
"""
Cd_rad_out_singleNd
1 *1 Cd 1 {2,S}
2 {Cs,O} 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_singleNd
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(3.94E+11,A_UNITS,"+-",0.0),
                 n=(0.69,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(44.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "666",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 54
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_Cs2_cy3
1 *1 C 1 {2,S}, {3,S}
2   Cs 0 {1,S}, {3,S}
3   Cs 0 {1,S}, {2,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.58E+09,A_UNITS,"+-",0.0),
                 n=(1.08,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(40.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "667",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 55
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy3
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {1,S}, {3,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.14E+10,A_UNITS,"+-",0.0),
                 n=(0.81,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(46.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "668",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 56
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_Cs2_cy3
1 *1 C 1 {2,S}, {3,S}
2   Cs 0 {1,S}, {3,S}
3   Cs 0 {1,S}, {2,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.33E+10,A_UNITS,"+-",0.0),
                 n=(0.65,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "669",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 57
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy3
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {1,S}, {3,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(2.74E+09,A_UNITS,"+-",0.0),
                 n=(0.98,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(46.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "670",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 58
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_Cs2_cy3
1 *1 C 1 {2,S}, {3,S}
2   Cs 0 {1,S}, {3,S}
3   Cs 0 {1,S}, {2,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(5.90E+11,A_UNITS,"+-",0.0),
                 n=(0.36,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "671",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 59
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy3
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {1,S}, {3,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.44E+08,A_UNITS,"+-",0.0),
                 n=(1.39,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(47.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "672",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 60
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy4
1 *2 Cs 0 {2,S} {3,S} {4,S} {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {5,S}
5 Cs 0 {3,S}, {4,S}
6 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(9.75E+09,A_UNITS,"+-",0.0),
                 n=(0.98,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "673",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 61
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_Cs2_cy4
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {4,S}
4 {Cs,Cd} 0 {2,S}, {3,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.44E+08,A_UNITS,"+-",0.0),
                 n=(1.2,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "674",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 62
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy4
1 *2 Cs 0 {2,S} {3,S} {4,S} {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {5,S}
5 Cs 0 {3,S}, {4,S}
6 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(5.64E+09,A_UNITS,"+-",0.0),
                 n=(1,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "675",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 63
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_Cs2_cy4
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {4,S}
4 {Cs,Cd} 0 {2,S}, {3,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.56E+09,A_UNITS,"+-",0.0),
                 n=(0.81,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "676",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 64
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy4
1 *2 Cs 0 {2,S} {3,S} {4,S} {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {5,S}
5 Cs 0 {3,S}, {4,S}
6 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(9.31E+08,A_UNITS,"+-",0.0),
                 n=(1.21,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "677",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 65
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_Cs2_cy4
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {4,S}
4 {Cs,Cd} 0 {2,S}, {3,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(4.86E+10,A_UNITS,"+-",0.0),
                 n=(0.58,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "678",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 66
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_Cs2_cy5
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {5,S}
4 {Cs,Cd,Cb,Ct} 0 {2,S}, {5,{S,D,T,B}}
5 {Cs,Cd,Cb,Ct} 0 {3,S}, {4,{S,D,T,B}}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.07E+09,A_UNITS,"+-",0.0),
                 n=(1.19,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(42.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "679",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 67
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy5
1 *2 Cs 0 {2,S} {3,S} {4,S} {7,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {6,S}
5 Cs 0 {3,S}, {6,S}
6 Cs 0 {4,S}, {5,S}
7 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(3.35E+09,A_UNITS,"+-",0.0),
                 n=(0.99,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "680",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 68
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy5
1 *2 Cs 0 {2,S} {3,S} {4,S} {7,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {6,S}
5 Cs 0 {3,S}, {6,S}
6 Cs 0 {4,S}, {5,S}
7 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(3.29E+09,A_UNITS,"+-",0.0),
                 n=(0.89,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "681",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 69
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_Cs2_cy5
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {5,S}
4 {Cs,Cd,Cb,Ct} 0 {2,S}, {5,{S,D,T,B}}
5 {Cs,Cd,Cb,Ct} 0 {3,S}, {4,{S,D,T,B}}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.08E+10,A_UNITS,"+-",0.0),
                 n=(0.81,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "682",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 70
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy5
1 *2 Cs 0 {2,S} {3,S} {4,S} {7,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {6,S}
5 Cs 0 {3,S}, {6,S}
6 Cs 0 {4,S}, {5,S}
7 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(7.48E+07,A_UNITS,"+-",0.0),
                 n=(1.45,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "683",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 71
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_Cs2_cy5
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {5,S}
4 {Cs,Cd,Cb,Ct} 0 {2,S}, {5,{S,D,T,B}}
5 {Cs,Cd,Cb,Ct} 0 {3,S}, {4,{S,D,T,B}}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.24E+11,A_UNITS,"+-",0.0),
                 n=(1.47,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "684",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 72
rate(
  group1 = 
"""
R2H_S_cy3
1 *1 {R!H} 1 {2,S}, {4,{S,D,B}}
2 *2 {R!H} 0 {1,S}, {3,S}, {4,{S,D,B}}
3 *3 H 0 {2,S}
4 {R!H} 0 {1,{S,D,B}}, {2,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.25E+11,A_UNITS,"+-",0.0),
                 n=(0.6,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(44.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "685",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 73
rate(
  group1 = 
"""
R2H_S_cy3
1 *1 {R!H} 1 {2,S}, {4,{S,D,B}}
2 *2 {R!H} 0 {1,S}, {3,S}, {4,{S,D,B}}
3 *3 H 0 {2,S}
4 {R!H} 0 {1,{S,D,B}}, {2,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.72E+12,A_UNITS,"+-",0.0),
                 n=(0.37,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "686",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 74
rate(
  group1 = 
"""
R2H_S_cy3
1 *1 {R!H} 1 {2,S}, {4,{S,D,B}}
2 *2 {R!H} 0 {1,S}, {3,S}, {4,{S,D,B}}
3 *3 H 0 {2,S}
4 {R!H} 0 {1,{S,D,B}}, {2,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(5.69E+11,A_UNITS,"+-",0.0),
                 n=(0.51,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(43.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "687",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 75
rate(
  group1 = 
"""
R2H_S_cy4
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *2 {R!H} 0 {1,S}, {3,S}, {4,{S,D,B}}
3 *3 H 0 {2,S}
4 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
5 {R!H} 0 {1,{S,D,B}}, {4,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.56E+12,A_UNITS,"+-",0.0),
                 n=(0.24,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "688",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 76
rate(
  group1 = 
"""
R2H_S_cy4
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *2 {R!H} 0 {1,S}, {3,S}, {4,{S,D,B}}
3 *3 H 0 {2,S}
4 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
5 {R!H} 0 {1,{S,D,B}}, {4,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.49E+10,A_UNITS,"+-",0.0),
                 n=(0.79,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(42.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "689",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 77
rate(
  group1 = 
"""
R2H_S_cy5
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *2 {R!H} 0 {1,S}, {3,S}, {4,{S,D,B}}
3 *3 H 0 {2,S}
4 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
5 {R!H} 0 {4,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {1,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.71E+11,A_UNITS,"+-",0.0),
                 n=(0.61,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "691",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 78
rate(
  group1 = 
"""
R2H_S_cy5
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *2 {R!H} 0 {1,S}, {3,S}, {4,{S,D,B}}
3 *3 H 0 {2,S}
4 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
5 {R!H} 0 {4,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {1,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(3.72E+12,A_UNITS,"+-",0.0),
                 n=(0.26,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "692",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 79
rate(
  group1 = 
"""
R2H_S_cy5
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *2 {R!H} 0 {1,S}, {3,S}, {4,{S,D,B}}
3 *3 H 0 {2,S}
4 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
5 {R!H} 0 {4,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {1,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.88E+11,A_UNITS,"+-",0.0),
                 n=(0.51,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "693",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 80
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.76E+08,A_UNITS,"+-",0.0),
                 n=(1.17,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "694",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 81
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.90E+09,A_UNITS,"+-",0.0),
                 n=(0.82,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "695",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 82
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.19E+08,A_UNITS,"+-",0.0),
                 n=(1.32,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "696",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 83
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.50E+08,A_UNITS,"+-",0.0),
                 n=(1.01,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "697",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 84
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(2.25E+10,A_UNITS,"+-",0.0),
                 n=(0.66,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(32.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "698",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 85
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.89E+06,A_UNITS,"+-",0.0),
                 n=(1.77,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "699",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 86
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(7.27E+09,A_UNITS,"+-",0.0),
                 n=(0.66,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "700",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 87
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.71E+07,A_UNITS,"+-",0.0),
                 n=(1.41,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "701",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 88
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(6.78E+08,A_UNITS,"+-",0.0),
                 n=(1,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "702",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 89
rate(
  group1 = 
"""
R3H_DS
1 *1 Cd 1 {2,D}
2 *4 Cd 0 {1,D}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Cd_rad_out_singleH
1 *1 Cd 1 {2,S}
2 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.10E+09,A_UNITS,"+-",0.0),
                 n=(0.97,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "703",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 90
rate(
  group1 = 
"""
R3H_SD
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,D}
3 *2 Cd 0 {2,D}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_singleH
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.11E+11,A_UNITS,"+-",0.0),
                 n=(0.58,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "704",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 91
rate(
  group1 = 
"""
R3H_DS
1 *1 Cd 1 {2,D}
2 *4 Cd 0 {1,D}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Cd_rad_out_singleH
1 *1 Cd 1 {2,S}
2 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.23E+09,A_UNITS,"+-",0.0),
                 n=(0.74,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "705",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 92
rate(
  group1 = 
"""
R3H_SD
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,D}
3 *2 Cd 0 {2,D}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_singleH
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.16E+10,A_UNITS,"+-",0.0),
                 n=(0.77,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(64.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "706",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 93
rate(
  group1 = 
"""
R3H_DS
1 *1 Cd 1 {2,D}
2 *4 Cd 0 {1,D}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Cd_rad_out_singleH
1 *1 Cd 1 {2,S}
2 H 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(6.04E+10,A_UNITS,"+-",0.0),
                 n=(0.59,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(32.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "707",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 94
rate(
  group1 = 
"""
R3H_SD
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,D}
3 *2 Cd 0 {2,D}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cd_H_out_singleH
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.53E+08,A_UNITS,"+-",0.0),
                 n=(1.27,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(63.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "708",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 95
rate(
  group1 = 
"""
R3H_DS
1 *1 Cd 1 {2,D}
2 *4 Cd 0 {1,D}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Cd_rad_out_singleNd
1 *1 Cd 1 {2,S}
2 {Cs,O} 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.58E+09,A_UNITS,"+-",0.0),
                 n=(1.08,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "709",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 96
rate(
  group1 = 
"""
R3H_SD
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,D}
3 *2 Cd 0 {2,D}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_singleNd
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.91E+11,A_UNITS,"+-",0.0),
                 n=(0.63,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(61.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "710",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 97
rate(
  group1 = 
"""
R3H_DS
1 *1 Cd 1 {2,D}
2 *4 Cd 0 {1,D}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Cd_rad_out_singleNd
1 *1 Cd 1 {2,S}
2 {Cs,O} 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.91E+09,A_UNITS,"+-",0.0),
                 n=(0.86,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "711",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 98
rate(
  group1 = 
"""
R3H_SD
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,D}
3 *2 Cd 0 {2,D}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_singleNd
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(3.96E+10,A_UNITS,"+-",0.0),
                 n=(0.83,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(61.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "712",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 99
rate(
  group1 = 
"""
R3H_DS
1 *1 Cd 1 {2,D}
2 *4 Cd 0 {1,D}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Cd_rad_out_singleNd
1 *1 Cd 1 {2,S}
2 {Cs,O} 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(8.05E+09,A_UNITS,"+-",0.0),
                 n=(0.86,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "713",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 100
rate(
  group1 = 
"""
R3H_SD
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,D}
3 *2 Cd 0 {2,D}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cd_H_out_singleNd
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(6.05E+10,A_UNITS,"+-",0.0),
                 n=(0.79,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(61.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "714",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 101
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Cd_rad_out_double
1 *1 Cd 1 {2,D}
2 Cd 0 {1,D}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.68E+08,A_UNITS,"+-",0.0),
                 n=(1.24,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "715",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 102
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_doubleC
1 *2 Cd 0 {2,S} {3,D}
2 *3 H 0 {1,S}
3 Cd 0 {1,D}
""",
  kf = Arrhenius(A=(3.24E+08,A_UNITS,"+-",0.0),
                 n=(1.14,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "716",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 103
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Cd_rad_out_double
1 *1 Cd 1 {2,D}
2 Cd 0 {1,D}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.66E+09,A_UNITS,"+-",0.0),
                 n=(0.99,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "717",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 104
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_doubleC
1 *2 Cd 0 {2,S} {3,D}
2 *3 H 0 {1,S}
3 Cd 0 {1,D}
""",
  kf = Arrhenius(A=(3.37E+07,A_UNITS,"+-",0.0),
                 n=(1.41,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(42.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "718",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 105
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Cd_rad_out_double
1 *1 Cd 1 {2,D}
2 Cd 0 {1,D}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.10E+10,A_UNITS,"+-",0.0),
                 n=(0.78,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(31.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "719",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 106
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cd_H_out_doubleC
1 *2 Cd 0 {2,S} {3,D}
2 *3 H 0 {1,S}
3 Cd 0 {1,D}
""",
  kf = Arrhenius(A=(3.50E+06,A_UNITS,"+-",0.0),
                 n=(1.68,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(42.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "720",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 107
rate(
  group1 = 
"""
R3H_SS_2Cd
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.93E+09,A_UNITS,"+-",0.0),
                 n=(1.26,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(52.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "721",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 108
rate(
  group1 = 
"""
R3H_SS_2Cd
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.20E+10,A_UNITS,"+-",0.0),
                 n=(0.82,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(49.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "722",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 109
rate(
  group1 = 
"""
R3H_SS_2Cd
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.64E+08,A_UNITS,"+-",0.0),
                 n=(1.47,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(53.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "723",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 110
rate(
  group1 = 
"""
R3H_SS_2Cd
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.43E+11,A_UNITS,"+-",0.0),
                 n=(0.65,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(47.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "724",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 111
rate(
  group1 = 
"""
R3H_SS_2Cd
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.12E+07,A_UNITS,"+-",0.0),
                 n=(1.72,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(51.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "725",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\')) 
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 112
rate(
  group1 = 
"""
R3H_SS_2Cd
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.21E+10,A_UNITS,"+-",0.0),
                 n=(0.91,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(49.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "726",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 113
rate(
  group1 = 
"""
R3H_SS_2Cd
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(3.46E+10,A_UNITS,"+-",0.0),
                 n=(0.76,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(47.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "727",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 114
rate(
  group1 = 
"""
R3H_SS_2Cd
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.14E+10,A_UNITS,"+-",0.0),
                 n=(0.8,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(48.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "728",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 115
rate(
  group1 = 
"""
R3H_SS_2Cd
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.76E+09,A_UNITS,"+-",0.0),
                 n=(1.18,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(43.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "729",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 116
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/OneDe
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.80E+09,A_UNITS,"+-",0.0),
                 n=(0.99,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(48.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "730",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 117
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.66E+08,A_UNITS,"+-",0.0),
                 n=(1.1,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "731",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 118
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/OneDe
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.77E+09,A_UNITS,"+-",0.0),
                 n=(0.74,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(46.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "732",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 119
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.41E+09,A_UNITS,"+-",0.0),
                 n=(0.73,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "733",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 120
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/OneDe
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(9.06E+10,A_UNITS,"+-",0.0),
                 n=(0.44,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(43.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "734",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 121
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/OneDe
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.40E+06,A_UNITS,"+-",0.0),
                 n=(1.56,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "735",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 122
rate(
  group1 = 
"""
R3H_SS_12cy3
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {1,{S,D,B}}, {2,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.62E+10,A_UNITS,"+-",0.0),
                 n=(0.69,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "736",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 123
rate(
  group1 = 
"""
R3H_SS_23cy3
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {3,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.04E+10,A_UNITS,"+-",0.0),
                 n=(0.71,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "737",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 124
rate(
  group1 = 
"""
R3H_SS_12cy3
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {1,{S,D,B}}, {2,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.16E+11,A_UNITS,"+-",0.0),
                 n=(0.26,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "738",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 125
rate(
  group1 = 
"""
R3H_SS_23cy3
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {3,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.63E+08,A_UNITS,"+-",0.0),
                 n=(1.01,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(45.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "739",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 126
rate(
  group1 = 
"""
R3H_SS_12cy3
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {1,{S,D,B}}, {2,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(2.26E+13,A_UNITS,"+-",0.0),
                 n=(0.26,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(32.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "740",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 127
rate(
  group1 = 
"""
R3H_SS_23cy3
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {3,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.68E+07,A_UNITS,"+-",0.0),
                 n=(1.42,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(46.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "741",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 128
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_Cs2_cy3
1 *1 C 1 {2,S}, {3,S}
2   Cs 0 {1,S}, {3,S}
3   Cs 0 {1,S}, {2,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.78E+09,A_UNITS,"+-",0.0),
                 n=(1.04,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "742",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 129
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy3
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {1,S}, {3,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(9.72E+09,A_UNITS,"+-",0.0),
                 n=(0.78,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "743",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 130
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_Cs2_cy3
1 *1 C 1 {2,S}, {3,S}
2   Cs 0 {1,S}, {3,S}
3   Cs 0 {1,S}, {2,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.39E+09,A_UNITS,"+-",0.0),
                 n=(0.77,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "744",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 131
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy3
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {1,S}, {3,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.73E+08,A_UNITS,"+-",0.0),
                 n=(1.14,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(40.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "745",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 132
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_Cs2_cy3
1 *1 C 1 {2,S}, {3,S}
2   Cs 0 {1,S}, {3,S}
3   Cs 0 {1,S}, {2,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(9.08E+10,A_UNITS,"+-",0.0),
                 n=(0.36,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(31.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "746",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 133
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy3
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {1,S}, {3,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(3.86E+06,A_UNITS,"+-",0.0),
                 n=(1.65,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "747",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 134
rate(
  group1 = 
"""
R3H_SS_12cy4
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {1,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.90E+10,A_UNITS,"+-",0.0),
                 n=(0.57,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "748",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 135
rate(
  group1 = 
"""
R3H_SS_23cy4
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {6,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.43E+09,A_UNITS,"+-",0.0),
                 n=(0.93,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "749",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 136
rate(
  group1 = 
"""
R3H_SS_12cy4
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {1,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.90E+11,A_UNITS,"+-",0.0),
                 n=(0.27,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "750",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 137
rate(
  group1 = 
"""
R3H_SS_23cy4
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {6,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.59E+08,A_UNITS,"+-",0.0),
                 n=(1.2,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "751",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 138
rate(
  group1 = 
"""
R3H_SS_12cy4
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {1,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.63E+12,A_UNITS,"+-",0.0),
                 n=(-0.04,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "752",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 139
rate(
  group1 = 
"""
R3H_SS_23cy4
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {6,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.19E+07,A_UNITS,"+-",0.0),
                 n=(1.55,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(40.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "753",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 140
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_Cs2_cy4
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {4,S}
4 {Cs,Cd} 0 {2,S}, {3,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.08E+08,A_UNITS,"+-",0.0),
                 n=(1.25,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "754",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 141
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy4
1 *2 Cs 0 {2,S} {3,S} {4,S} {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {5,S}
5 Cs 0 {3,S}, {4,S}
6 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(5.09E+09,A_UNITS,"+-",0.0),
                 n=(0.84,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "755",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 142
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_Cs2_cy4
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {4,S}
4 {Cs,Cd} 0 {2,S}, {3,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.05E+08,A_UNITS,"+-",0.0),
                 n=(0.99,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "756",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 143
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy4
1 *2 Cs 0 {2,S} {3,S} {4,S} {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {5,S}
5 Cs 0 {3,S}, {4,S}
6 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(5.69E+08,A_UNITS,"+-",0.0),
                 n=(0.97,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "757",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 144
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_Cs2_cy4
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {4,S}
4 {Cs,Cd} 0 {2,S}, {3,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(8.20E+09,A_UNITS,"+-",0.0),
                 n=(0.54,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "758",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 145
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy4
1 *2 Cs 0 {2,S} {3,S} {4,S} {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {5,S}
5 Cs 0 {3,S}, {4,S}
6 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(3.49E+07,A_UNITS,"+-",0.0),
                 n=(1.38,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "759",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 146
rate(
  group1 = 
"""
R3H_SS_12cy5
1 *1 {R!H} 1 {2,S}, {7,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {1,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.85E+10,A_UNITS,"+-",0.0),
                 n=(0.6,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(40.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "760",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 147
rate(
  group1 = 
"""
R3H_SS_23cy5
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {7,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.25E+09,A_UNITS,"+-",0.0),
                 n=(0.99,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "761",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 148
rate(
  group1 = 
"""
R3H_SS_12cy5
1 *1 {R!H} 1 {2,S}, {7,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {1,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.67E+11,A_UNITS,"+-",0.0),
                 n=(0.29,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "762",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 149
rate(
  group1 = 
"""
R3H_SS_23cy5
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {7,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.42E+08,A_UNITS,"+-",0.0),
                 n=(1.14,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "763",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 150
rate(
  group1 = 
"""
R3H_SS_12cy5
1 *1 {R!H} 1 {2,S}, {7,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {1,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(9.42E+11,A_UNITS,"+-",0.0),
                 n=(0.12,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "764",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 151
rate(
  group1 = 
"""
R3H_SS_23cy5
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {7,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.58E+06,A_UNITS,"+-",0.0),
                 n=(1.78,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "765",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 152
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_Cs2_cy5
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {5,S}
4 {Cs,Cd,Cb,Ct} 0 {2,S}, {5,{S,D,T,B}}
5 {Cs,Cd,Cb,Ct} 0 {3,S}, {4,{S,D,T,B}}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.14E+08,A_UNITS,"+-",0.0),
                 n=(1.26,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "766",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 153
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy5
1 *2 Cs 0 {2,S} {3,S} {4,S} {7,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {6,S}
5 Cs 0 {3,S}, {6,S}
6 Cs 0 {4,S}, {5,S}
7 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(6.90E+09,A_UNITS,"+-",0.0),
                 n=(0.82,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(32.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "767",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 154
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_Cs2_cy5
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {5,S}
4 {Cs,Cd,Cb,Ct} 0 {2,S}, {5,{S,D,T,B}}
5 {Cs,Cd,Cb,Ct} 0 {3,S}, {4,{S,D,T,B}}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.25E+08,A_UNITS,"+-",0.0),
                 n=(1.01,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "768",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 155
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy5
1 *2 Cs 0 {2,S} {3,S} {4,S} {7,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {6,S}
5 Cs 0 {3,S}, {6,S}
6 Cs 0 {4,S}, {5,S}
7 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(7.50E+08,A_UNITS,"+-",0.0),
                 n=(0.9,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "769",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 156
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_Cs2_cy5
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}, {4,S}
3 Cs 0 {1,S}, {5,S}
4 {Cs,Cd,Cb,Ct} 0 {2,S}, {5,{S,D,T,B}}
5 {Cs,Cd,Cb,Ct} 0 {3,S}, {4,{S,D,T,B}}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.97E+10,A_UNITS,"+-",0.0),
                 n=(0.46,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "770",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 157
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_Cs2_cy5
1 *2 Cs 0 {2,S} {3,S} {4,S} {7,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {5,S}
4 Cs 0 {1,S}, {6,S}
5 Cs 0 {3,S}, {6,S}
6 Cs 0 {4,S}, {5,S}
7 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(2.21E+08,A_UNITS,"+-",0.0),
                 n=(1.04,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "771",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 158
rate(
  group1 = 
"""
R3H_SS_12cy3
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {1,{S,D,B}}, {2,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.64E+09,A_UNITS,"+-",0.0),
                 n=(0.84,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "772",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 159
rate(
  group1 = 
"""
R3H_SS_23cy3
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {3,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(5.02E+10,A_UNITS,"+-",0.0),
                 n=(0.56,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(42.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "773",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 160
rate(
  group1 = 
"""
R3H_SS_12cy3
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {1,{S,D,B}}, {2,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.22E+11,A_UNITS,"+-",0.0),
                 n=(0.4,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "774",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 161
rate(
  group1 = 
"""
R3H_SS_23cy3
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {3,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(4.34E+09,A_UNITS,"+-",0.0),
                 n=(0.81,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(43.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "775",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 162
rate(
  group1 = 
"""
R3H_SS_12cy3
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {1,{S,D,B}}, {2,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(2.72E+12,A_UNITS,"+-",0.0),
                 n=(-0.04,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "776",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 163
rate(
  group1 = 
"""
R3H_SS_23cy3
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {3,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.61E+08,A_UNITS,"+-",0.0),
                 n=(1.26,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(42.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "777",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 164
rate(
  group1 = 
"""
R3H_SS_13cy4
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {1,{S,D,B}}, {3,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.78E+11,A_UNITS,"+-",0.0),
                 n=(0.29,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(54.30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "778",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 165
rate(
  group1 = 
"""
R3H_SS_13cy4
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {1,{S,D,B}}, {3,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.48E+10,A_UNITS,"+-",0.0),
                 n=(0.6,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(54.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "779",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 166
rate(
  group1 = 
"""
R3H_SS_13cy4
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {1,{S,D,B}}, {3,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(2.66E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(51.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "780",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 167
rate(
  group1 = 
"""
R3H_SS_13cy4
1 *1 {R!H} 1 {2,S}, {5,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {1,{S,D,B}}, {3,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(3.55E+11,A_UNITS,"+-",0.0),
                 n=(0.37,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(51.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "781",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 168
rate(
  group1 = 
"""
R3H_SS_13cy5
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {1,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.36E+11,A_UNITS,"+-",0.0),
                 n=(0.46,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(47.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "782",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 169
rate(
  group1 = 
"""
R3H_SS_13cy5
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {1,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.72E+09,A_UNITS,"+-",0.0),
                 n=(0.86,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(47.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "783",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 170
rate(
  group1 = 
"""
R3H_SS_13cy5
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {1,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.10E+12,A_UNITS,"+-",0.0),
                 n=(0.23,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(44.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "784",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 171
rate(
  group1 = 
"""
R3H_SS_13cy5
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {1,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(6.07E+10,A_UNITS,"+-",0.0),
                 n=(0.62,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(43.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "785",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 172
rate(
  group1 = 
"""
R3H_SS_12cy5
1 *1 {R!H} 1 {2,S}, {7,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {1,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.84E+09,A_UNITS,"+-",0.0),
                 n=(1.05,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "786",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 173
rate(
  group1 = 
"""
R3H_SS_23cy5
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {7,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(4.51E+09,A_UNITS,"+-",0.0),
                 n=(0.86,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "787",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 174
rate(
  group1 = 
"""
R3H_SS_12cy5
1 *1 {R!H} 1 {2,S}, {7,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {1,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.04E+09,A_UNITS,"+-",0.0),
                 n=(0.74,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "788",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 175
rate(
  group1 = 
"""
R3H_SS_23cy5
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {7,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.95E+09,A_UNITS,"+-",0.0),
                 n=(0.88,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "789",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 176
rate(
  group1 = 
"""
R3H_SS_12cy5
1 *1 {R!H} 1 {2,S}, {7,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {1,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.44E+10,A_UNITS,"+-",0.0),
                 n=(0.74,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "790",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 177
rate(
  group1 = 
"""
R3H_SS_23cy5
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {7,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {5,{S,D,B}}, {7,{S,D,B}}
7 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(2.85E+07,A_UNITS,"+-",0.0),
                 n=(1.46,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "791",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 178
rate(
  group1 = 
"""
R3H_SS_12cy4
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {1,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.28E+08,A_UNITS,"+-",0.0),
                 n=(1.07,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(40.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "792",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 179
rate(
  group1 = 
"""
R3H_SS_23cy4
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {6,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.93E+10,A_UNITS,"+-",0.0),
                 n=(0.75,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "793",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 180
rate(
  group1 = 
"""
R3H_SS_12cy4
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {1,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.41E+09,A_UNITS,"+-",0.0),
                 n=(0.77,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "794",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 181
rate(
  group1 = 
"""
R3H_SS_23cy4
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {6,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.96E+09,A_UNITS,"+-",0.0),
                 n=(0.96,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "795",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 182
rate(
  group1 = 
"""
R3H_SS_12cy4
1 *1 {R!H} 1 {2,S}, {6,{S,D,B}}
2 *4 {R!H} 0 {1,S}, {3,S}, {5,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
5 {R!H} 0 {2,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {1,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(2.37E+10,A_UNITS,"+-",0.0),
                 n=(0.62,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "796",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 183
rate(
  group1 = 
"""
R3H_SS_23cy4
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}, {6,{S,D,B}}
3 *2 {R!H} 0 {2,S}, {4,S}, {5,{S,D,B}}
4 *3 H 0 {3,S}
5 {R!H} 0 {3,{S,D,B}}, {6,{S,D,B}}
6 {R!H} 0 {2,{S,D,B}}, {5,{S,D,B}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(4.96E+07,A_UNITS,"+-",0.0),
                 n=(1.46,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "797",
  short_comment = "Sumathy B3LYP/CCPVDZ calculations",
  long_comment = 
"""
Sumathy B3LYP/CCPVDZ calculations (hindered rotor potential barrier calculations at B3LYP/6-31G(d\'))
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 184
rate(
  group1 = 
"""
R3H_SS_OC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *2 Cs 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.71E+08,A_UNITS,"+-",0.0),
                 n=(1.45,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(42.27,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "798",
  short_comment = "Sumathy CBS-Q calculations",
  long_comment = 
"""
Sumathy CBS-Q calculations
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 185
rate(
  group1 = 
"""
R3H_SS_OC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *2 Cs 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.66E+08,A_UNITS,"+-",0.0),
                 n=(1.28,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.74,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "799",
  short_comment = "Sumathy CBS-Q calculations",
  long_comment = 
"""
Sumathy CBS-Q calculations
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 186
rate(
  group1 = 
"""
R3H_SS_OC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *2 Cs 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/(NonDeC/Cs)
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S} {5,S}
4 H 0 {1,S}
5 Cs 0 {3,S}
""",
  kf = Arrhenius(A=(2.43E+09,A_UNITS,"+-",0.0),
                 n=(1.17,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.66,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "800",
  short_comment = "Sumathy CBS-Q calculations",
  long_comment = 
"""
Sumathy CBS-Q calculations
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 187
rate(
  group1 = 
"""
R3H_SS_OC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *2 Cs 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/(NonDeC/Cs/Cs)
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S} {5,S} {6,S}
4 H 0 {1,S}
5 Cs 0 {3,S}
6 Cs 0 {3,S}
""",
  kf = Arrhenius(A=(1.07E+10,A_UNITS,"+-",0.0),
                 n=(0.98,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.58,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "801",
  short_comment = "Sumathy CBS-Q calculations",
  long_comment = 
"""
Sumathy CBS-Q calculations
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 188
rate(
  group1 = 
"""
R3H_SS_OC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *2 Cs 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/(NonDeC/Cs/Cs/Cs)
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S} {5,S} {6,S} {7,S}
4 H 0 {1,S}
5 Cs 0 {3,S}
6 Cs 0 {3,S}
7 Cs 0 {3,S}
""",
  kf = Arrhenius(A=(6.62E+09,A_UNITS,"+-",0.0),
                 n=(1.04,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.34,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "802",
  short_comment = "Sumathy CBS-Q calculations",
  long_comment = 
"""
Sumathy CBS-Q calculations
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 189
rate(
  group1 = 
"""
R3H_SS_OC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *2 Cs 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(4.97E+09,A_UNITS,"+-",0.0),
                 n=(1.01,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.47,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "803",
  short_comment = "Sumathy CBS-Q calculations",
  long_comment = 
"""
Sumathy CBS-Q calculations
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 190
rate(
  group1 = 
"""
R4H_SSS_OOCsCs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.99E+11,A_UNITS,"+-",0.0),
                 n=(0.15,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.2107177,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "804",
  short_comment = "Sumathy CBS-Q calculations",
  long_comment = 
"""
Sumathy CBS-Q calculations
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 191
rate(
  group1 = 
"""
R4H_SSS_OO(Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}, {6,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
6    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.38E+08,A_UNITS,"+-",0.0),
                 n=(1.06,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.51,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "805",
  short_comment = "Sumathy CBS-Q calculations",
  long_comment = 
"""
Sumathy CBS-Q calculations
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 192
rate(
  group1 = 
"""
R4H_SSS_OO(Cs/Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}, {6,S}, {7,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
6    Cs 0 {3,S}
7    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.06E+08,A_UNITS,"+-",0.0),
                 n=(1.20,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(33.53,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "806",
  short_comment = "Sumathy CBS-Q calculations",
  long_comment = 
"""
Sumathy CBS-Q calculations
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 193
rate(
  group1 = 
"""
R4H_SSS_OOCsCs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.00E+08,A_UNITS,"+-",0.0),
                 n=(1.10,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30.09,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "807",
  short_comment = "",
  long_comment = 
"""

""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 194
rate(
  group1 = 
"""
R4H_SSS_OO(Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}, {6,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
6    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.81E+08,A_UNITS,"+-",0.0),
                 n=(0.88,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.48,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "808",
  short_comment = "",
  long_comment = 
"""

""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 195
rate(
  group1 = 
"""
R4H_SSS_OO(Cs/Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}, {6,S}, {7,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
6    Cs 0 {3,S}
7    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.53E+09,A_UNITS,"+-",0.0),
                 n=(0.69,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30.11,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "809",
  short_comment = "",
  long_comment = 
"""

""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 196
rate(
  group1 = 
"""
R4H_SSS_OOCsCs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(2.64E+09,A_UNITS,"+-",0.0),
                 n=(0.78,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(27.11,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "810",
  short_comment = "",
  long_comment = 
"""

""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 197
rate(
  group1 = 
"""
R4H_SSS_OO(Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}, {6,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
6    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(9.25E+09,A_UNITS,"+-",0.0),
                 n=(0.57,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(27.31,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "811",
  short_comment = "",
  long_comment = 
"""

""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 198
rate(
  group1 = 
"""
R4H_SSS_OO(Cs/Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}, {6,S}, {7,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
6    Cs 0 {3,S}
7    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(4.87E+10,A_UNITS,"+-",0.0),
                 n=(0.35,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.39,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "812",
  short_comment = "",
  long_comment = 
"""

""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 199
rate(
  group1 = 
"""
R5H_SSSS_OOCCC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 Cs 0 {2,S}, {4,S}
4 *5 Cs 0 {3,S}, {5,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1690000,A_UNITS,"+-",0.0),
                 n=(1.55,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(21.02,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "813",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman).",
  long_comment = 
"""
CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 200
rate(
  group1 = 
"""
R5H_SSSS_OO(Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}, {7,S}
4 *5 Cs 0 {3,S}, {5,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H  0 {5,S}
7    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6780000,A_UNITS,"+-",0.0),
                 n=(1.35,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.84,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "814",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman).",
  long_comment = 
"""
CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 201
rate(
  group1 = 
"""
R5H_SSSS_OO(Cs/Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}, {7,S}, {8,S}
4 *5 Cs 0 {3,S}, {5,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H  0 {5,S}
7    Cs 0 {3,S}
8    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.35E+07,A_UNITS,"+-",0.0),
                 n=(1.12,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(21.88,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "815",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman).",
  long_comment = 
"""
CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 202
rate(
  group1 = 
"""
R5H_SSSS_OOCs(Cs/Cs)
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}
4 *5 Cs 0 {3,S}, {5,S}, {7,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H  0 {5,S}
7    Cs 0 {4,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.41E+07,A_UNITS,"+-",0.0),
                 n=(1.32,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(21.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "816",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman).",
  long_comment = 
"""
CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 203
rate(
  group1 = 
"""
R5H_SSSS_OOCs(Cs/Cs/Cs)
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}
4 *5 Cs 0 {3,S}, {5,S}, {7,S}, {8,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H  0 {5,S}
7    Cs 0 {4,S}
8    Cs 0 {4,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.09E+08,A_UNITS,"+-",0.0),
                 n=(1.23,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(21.62,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "817",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman).",
  long_comment = 
"""
CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 204
rate(
  group1 = 
"""
R5H_SSSS_OOCCC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 Cs 0 {2,S}, {4,S}
4 *5 Cs 0 {3,S}, {5,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.94E+06,A_UNITS,"+-",0.0),
                 n=(1.26,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.17,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "818",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman).",
  long_comment = 
"""
CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 205
rate(
  group1 = 
"""
R5H_SSSS_OO(Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}, {7,S}
4 *5 Cs 0 {3,S}, {5,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H  0 {5,S}
7    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.38E+10,A_UNITS,"+-",0.0),
                 n=(0.21,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "819",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman).",
  long_comment = 
"""
CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 206
rate(
  group1 = 
"""
R5H_SSSS_OOCCC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 Cs 0 {2,S}, {4,S}
4 *5 Cs 0 {3,S}, {5,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(31.74E+07,A_UNITS,"+-",0.0),
                 n=(1.15,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.37,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "820",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman).",
  long_comment = 
"""
CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 207
rate(
  group1 = 
"""
R6H_SSSSS_OO
1 *1 Os 1        {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}
4    Cs 0 {3,S}, {5,S}
5 *5 Cs 0 {4,S}, {6,S}
6 *2 Cs 0 {5,S}, {7,S}
7 *3 H 0 {6,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.69E+05,A_UNITS,"+-",0.0),
                 n=(1.52,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.05,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "821",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman).",
  long_comment = 
"""
CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 208
rate(
  group1 = 
"""
R6H_SSSSS_OO
1 *1 Os 1        {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}
4    Cs 0 {3,S}, {5,S}
5 *5 Cs 0 {4,S}, {6,S}
6 *2 Cs 0 {5,S}, {7,S}
7 *3 H 0 {6,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.62E+06,A_UNITS,"+-",0.0),
                 n=(1.22,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "822",
  short_comment = "Curran\'s [8] estimation in reaction type 19, QOOH = cyclic ether + OH",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253. 
Curran\'s estimation in reaction type 19, QOOH = cyclic ether + OH
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 209
rate(
  group1 = 
"""
R6H_SSSSS_OO
1 *1 Os 1        {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}
4    Cs 0 {3,S}, {5,S}
5 *5 Cs 0 {4,S}, {6,S}
6 *2 Cs 0 {5,S}, {7,S}
7 *3 H 0 {6,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(1.48E+06,A_UNITS,"+-",0.0),
                 n=(1.22,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.84,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "823",
  short_comment = "Curran\'s [8] estimation in reaction type 19, QOOH = cyclic ether + OH",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimation in reaction type 19, QOOH = cyclic ether + OH
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 210
rate(
  group1 = 
"""
R7H_OOCs4
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,{S,D,T,B}}
4 {R!H} 0 {3,{S,D,T,B}}, {5,{S,D,T,B}}
5 {R!H} 0 {4,{S,D,T,B}}, {6,{S,D,T,B}}
6 *5 {R!H} 0 {5,{S,D,T,B}}, {7,{S,D,T,B}}
7 *2 {R!H} 0 {6,{S,D,T,B}}, {8,S}
8 *3 H 0 {7,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.06E+04,A_UNITS,"+-",0.0),
                 n=(1.51,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.95,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "824",
  short_comment = "Curran\'s [8] estimation in reaction type 19, QOOH = cyclic ether + OH",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Curran\'s estimation in reaction type 19, QOOH = cyclic ether + OH
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 211
rate(
  group1 = 
"""
R7H_OOCs4
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,{S,D,T,B}}
4 {R!H} 0 {3,{S,D,T,B}}, {5,{S,D,T,B}}
5 {R!H} 0 {4,{S,D,T,B}}, {6,{S,D,T,B}}
6 *5 {R!H} 0 {5,{S,D,T,B}}, {7,{S,D,T,B}}
7 *2 {R!H} 0 {6,{S,D,T,B}}, {8,S}
8 *3 H 0 {7,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeC
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.37E+06,A_UNITS,"+-",0.0),
                 n=(0.99,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.17,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "825",
  short_comment = "Curran\'s [8] estimation in reaction type 19, QOOH = cyclic ether + OH",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253. 
Curran\'s estimation in reaction type 19, QOOH = cyclic ether + OH
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 212
rate(
  group1 = 
"""
R7H_OOCs4
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,{S,D,T,B}}
4 {R!H} 0 {3,{S,D,T,B}}, {5,{S,D,T,B}}
5 {R!H} 0 {4,{S,D,T,B}}, {6,{S,D,T,B}}
6 *5 {R!H} 0 {5,{S,D,T,B}}, {7,{S,D,T,B}}
7 *2 {R!H} 0 {6,{S,D,T,B}}, {8,S}
8 *3 H 0 {7,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Others-Cs_H_out_Cs2
AND{ Cs_H_out_Cs2, NOT OR{Cs_H_out_Cs2_cy3, Cs_H_out_Cs2_cy4, Cs_H_out_Cs2_cy5}}
""",
  kf = Arrhenius(A=(5.62E+05,A_UNITS,"+-",0.0),
                 n=(1.09,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.28,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "826",
  short_comment = "",
  long_comment = 
"""

""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 213
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/(CCCOOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {3,S}, {5,S}
5 Os 0 {4,S}, {6,S}
6 Os 0 {5,S}
""",
  kf = Arrhenius(A=(4.07E+09,A_UNITS,"+-",0.0),
                 n=(0.99,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.33,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "827",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya).",
  long_comment = 
"""
CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 214
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/(CCCOOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {3,S}, {5,S}
5 Os 0 {4,S}, {6,S}
6 Os 0 {5,S}
""",
  kf = Arrhenius(A=(6.70E+08,A_UNITS,"+-",0.0),
                 n=(1.15,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.04,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "828",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya).",
  long_comment = 
"""
CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 215
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/(CCCOOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {3,S}, {5,S}
5 Os 0 {4,S}, {6,S}
6 Os 0 {5,S}
""",
  kf = Arrhenius(A=(3.56E+07,A_UNITS,"+-",0.0),
                 n=(1.53,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(40.58,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "829",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya).",
  long_comment = 
"""
CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 216
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)CCOOH)
1 *2 Cs 0 {2,S}, {3,S}, {7,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {3,S}, {5,S}
5 Os 0 {4,S}, {6,S}
6 Os 0 {5,S}
7 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.83E+11,A_UNITS,"+-",0.0),
                 n=(0.45,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.92,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "830",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya).",
  long_comment = 
"""
CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 217
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)CCOOH)
1 *2 Cs 0 {2,S}, {3,S}, {7,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Cs 0 {3,S}, {5,S}
5 Os 0 {4,S}, {6,S}
6 Os 0 {5,S}
7 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(6.89E+10,A_UNITS,"+-",0.0),
                 n=(0.43,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.88,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "831",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya).",
  long_comment = 
"""
CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 218
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/(CCOOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
""",
  kf = Arrhenius(A=(5.87E+08,A_UNITS,"+-",0.0),
                 n=(1.28,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "832",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya).",
  long_comment = 
"""
CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 219
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/(CCOOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
""",
  kf = Arrhenius(A=(1.75E+08,A_UNITS,"+-",0.0),
                 n=(1.29,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.93,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "833",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya).",
  long_comment = 
"""
CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 220
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/(CCOOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
""",
  kf = Arrhenius(A=(2.14E+08,A_UNITS,"+-",0.0),
                 n=(1.42,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.71,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "834",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya).",
  long_comment = 
"""
CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 221
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)COOH)
1 *2 Cs 0 {2,S}, {3,S}, {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
6 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.28E+09,A_UNITS,"+-",0.0),
                 n=(1.12,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.69,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "835",
  short_comment = "Sumathy\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 222
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)COOH)
1 *2 Cs 0 {2,S}, {3,S}, {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
6 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.53E+10,A_UNITS,"+-",0.0),
                 n=(0.68,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.43,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "836",
  short_comment = "Sumathy\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 223
rate(
  group1 = 
"""
Others-R2H_S
AND{R2H_S, NOT OR{R2H_S,R2H_S_cy3,R2H_S_cy4,R2H_S_cy5}}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)COOH)
1 *2 Cs 0 {2,S}, {3,S}, {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
6 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.24E+09,A_UNITS,"+-",0.0),
                 n=(1.11,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.38,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "837",
  short_comment = "Sumathy\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 224
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/(CCOOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
""",
  kf = Arrhenius(A=(6.02E+08,A_UNITS,"+-",0.0),
                 n=(1.11,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.56,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "838",
  short_comment = "Sumathy\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 225
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/(CCOOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
""",
  kf = Arrhenius(A=(1.63E+07,A_UNITS,"+-",0.0),
                 n=(1.54,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.27,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "839",
  short_comment = "Sumathy\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 226
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/(CCOOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
""",
  kf = Arrhenius(A=(3.13E+05,A_UNITS,"+-",0.0),
                 n=(2.04,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.64,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "840",
  short_comment = "Sumathy\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 227
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)COOH)
1 *2 Cs 0 {2,S}, {3,S}, {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
6 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(6.95E+09,A_UNITS,"+-",0.0),
                 n=(0.79,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.71,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "841",
  short_comment = "Sumathy\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 228
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)COOH)
1 *2 Cs 0 {2,S}, {3,S}, {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
6 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.53E+10,A_UNITS,"+-",0.0),
                 n=(0.68,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.43,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "842",
  short_comment = "Sumathy\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 229
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)COOH)
1 *2 Cs 0 {2,S}, {3,S}, {6,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}, {4,S}
4 Os 0 {3,S}, {5,S}
5 Os 0 {4,S}
6 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.44E+09,A_UNITS,"+-",0.0),
                 n=(0.8,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.84,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "843",
  short_comment = "Sumathy\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 230
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/(COOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
""",
  kf = Arrhenius(A=(1.51E+08,A_UNITS,"+-",0.0),
                 n=(1.16,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.24,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "844",
  short_comment = "Sumathy\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 231
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/(COOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
""",
  kf = Arrhenius(A=(1.37E+07,A_UNITS,"+-",0.0),
                 n=(1.36,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(37.15,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "845",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 232
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/(COOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
""",
  kf = Arrhenius(A=(4.08E+06,A_UNITS,"+-",0.0),
                 n=(1.55,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.68,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "846",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 233
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)OOH)
1 *2 Cs 0 {2,S}, {3,S}, {5,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
5 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.69E+09,A_UNITS,"+-",0.0),
                 n=(0.68,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.81,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "847",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 234
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)OOH)
1 *2 Cs 0 {2,S}, {3,S}, {5,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
5 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(5.18E+08,A_UNITS,"+-",0.0),
                 n=(0.87,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.12,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "848",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 235
rate(
  group1 = 
"""
Others-R3H_SS
AND{ R3H_SS, NOT OR{R3H_SS_12cy3, R3H_SS_23cy3, R3H_SS_12cy4, R3H_SS_23cy4, R3H_SS_13cy4, R3H_SS_12cy5, R3H_SS_23cy5, R3H_SS_13cy5, R3H_SS_2Cd, R3H_SS_OC }}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)OOH)
1 *2 Cs 0 {2,S}, {3,S}, {5,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
5 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.07E+07,A_UNITS,"+-",0.0),
                 n=(1.37,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(35.66,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "849",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 236
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/(COOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
""",
  kf = Arrhenius(A=(8.73E+06,A_UNITS,"+-",0.0),
                 n=(1.13,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.93,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "8441",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 237
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/(COOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
""",
  kf = Arrhenius(A=(6.95E+05,A_UNITS,"+-",0.0),
                 n=(1.38,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.67,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "8451",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 238
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/(COOH)
1 *2 Cs 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
""",
  kf = Arrhenius(A=(1.74E+04,A_UNITS,"+-",0.0),
                 n=(1.89,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(18.51,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "8461",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 239
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)OOH)
1 *2 Cs 0 {2,S}, {3,S}, {5,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
5 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.09E+09,A_UNITS,"+-",0.0),
                 n=(0.55,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(17.43,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "8471",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 240
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S},  {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)OOH)
1 *2 Cs 0 {2,S}, {3,S}, {5,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
5 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(7.86E+07,A_UNITS,"+-",0.0),
                 n=(0.84,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.38,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "8481",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 241
rate(
  group1 = 
"""
R4H_SSS
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
Others-C_rad_out_Cs2
AND{C_rad_out_Cs2, NOT OR{C_rad_out_Cs2_cy3, C_rad_out_Cs2_cy4, C_rad_out_Cs2_cy5 }}
""",
  group3 = 
"""
Cs_H_out_H/((C/C)OOH)
1 *2 Cs 0 {2,S}, {3,S}, {5,S}
2 *3 H 0 {1,S}
3 Os 0 {1,S}, {4,S}
4 Os 0 {3,S}
5 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.05E+06,A_UNITS,"+-",0.0),
                 n=(1.39,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.19,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "8491",
  short_comment = "",
  long_comment = 
"""

The rate was added by Sandeep Sharma on Feb 01 2006 (713a0f98f91) with message \"I have added the rate rules given to me by Sumathy which are present in Pitz.xls file in the Acads directory of my laptop.\"
The node was commented out of the tree, disabling the rate, by Sandeep Sharma on Feb 13 2006 (2e7b38d367c9) with the message \"Removed nodes Cs_H_out_H/(CCCOOH) and the others from under node Cs_H_out_H/NonDeC as it is not a subnode anyway.\"

On 6 April 2010, Josh Allen, Mike Harper and Richard West spent quite a while trying to put these nodes in the right place in the tree and to make the definitions valid and consistent.
Unfortunately it was not clear what they were intended to mean because many of the definitions overlap. We gave up, and they remain commented out.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 242
rate(
  group1 = 
"""
R3H_SS_OC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *2 Cs 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeO
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.00e8,A_UNITS,"+-",0.0),
                 n=(1.23,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.85,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "850",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 243
rate(
  group1 = 
"""
R3H_SS_OC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *2 Cs 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_NDMustO
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 {Cs,O} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(3.00e8,A_UNITS,"+-",0.0),
                 n=(1.23,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(36.85,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "851",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 244
rate(
  group1 = 
"""
R4H_SSS_OOCsCs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeO
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.61e8,A_UNITS,"+-",0.0),
                 n=(1.09,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.14,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "852",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 245
rate(
  group1 = 
"""
R4H_SSS_OOCsCs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_NDMustO
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 {Cs,O} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(5.29e9,A_UNITS,"+-",0.0),
                 n=(0.75,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(24.82,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "855",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 246
rate(
  group1 = 
"""
R4H_SSS_OO(Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}, {6,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
6    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeO
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.20e8,A_UNITS,"+-",0.0),
                 n=(0.82,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.28,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "855",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 247
rate(
  group1 = 
"""
R4H_SSS_OO(Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 *5 Cs 0 {2,S}, {4,S}, {6,S}
4 *2 Cs 0 {3,S}, {5,S}
5 *3 H  0 {4,S}
6    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_NDMustO
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 {Cs,O} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.18e11,A_UNITS,"+-",0.0),
                 n=(0.51,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.20,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "856",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 248
rate(
  group1 = 
"""
R5H_SSSS_OOCCC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 Cs 0 {2,S}, {4,S}
4 *5 Cs 0 {3,S}, {5,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeO
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.15e4,A_UNITS,"+-",0.0),
                 n=(2.11,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.47,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "858",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 249
rate(
  group1 = 
"""
R5H_SSSS_OOCCC
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 Cs 0 {2,S}, {4,S}
4 *5 Cs 0 {3,S}, {5,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_NDMustO
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 {Cs,O} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.9e7,A_UNITS,"+-",0.0),
                 n=(1.1,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "863",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 250
rate(
  group1 = 
"""
R5H_SSSS_OO(Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}, {7,S}
4 *5 Cs 0 {3,S}, {5,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H  0 {5,S}
7    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeO
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.29e8,A_UNITS,"+-",0.0),
                 n=(1.12,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.38,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "864",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 251
rate(
  group1 = 
"""
R5H_SSSS_OO(Cs/Cs)Cs
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}, {7,S}
4 *5 Cs 0 {3,S}, {5,S}
5 *2 Cs 0 {4,S}, {6,S}
6 *3 H  0 {5,S}
7    Cs 0 {3,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_NDMustO
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 {Cs,O} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.17e11,A_UNITS,"+-",0.0),
                 n=(0.43,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.40,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "865",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 252
rate(
  group1 = 
"""
R6H_SSSSS_OO
1 *1 Os 1        {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}
4    Cs 0 {3,S}, {5,S}
5 *5 Cs 0 {4,S}, {6,S}
6 *2 Cs 0 {5,S}, {7,S}
7 *3 H 0 {6,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeO
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.49e2,A_UNITS,"+-",0.0),
                 n=(2.21,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.38,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "866",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 253
rate(
  group1 = 
"""
R6H_SSSSS_OO
1 *1 Os 1        {2,S}
2 *4 Os 0 {1,S}, {3,S}
3    Cs 0 {2,S}, {4,S}
4    Cs 0 {3,S}, {5,S}
5 *5 Cs 0 {4,S}, {6,S}
6 *2 Cs 0 {5,S}, {7,S}
7 *3 H 0 {6,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_NDMustO
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 {Cs,O} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(7.23e4,A_UNITS,"+-",0.0),
                 n=(1.65,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.02,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "867",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 254
rate(
  group1 = 
"""
R7H_OOCs4
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,{S,D,T,B}}
4 {R!H} 0 {3,{S,D,T,B}}, {5,{S,D,T,B}}
5 {R!H} 0 {4,{S,D,T,B}}, {6,{S,D,T,B}}
6 *5 {R!H} 0 {5,{S,D,T,B}}, {7,{S,D,T,B}}
7 *2 {R!H} 0 {6,{S,D,T,B}}, {8,S}
8 *3 H 0 {7,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_H/NonDeO
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.33e5,A_UNITS,"+-",0.0),
                 n=(0.75,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.82,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "868",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 255
rate(
  group1 = 
"""
R7H_OOCs4
1 *1 Os 1 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,{S,D,T,B}}
4 {R!H} 0 {3,{S,D,T,B}}, {5,{S,D,T,B}}
5 {R!H} 0 {4,{S,D,T,B}}, {6,{S,D,T,B}}
6 *5 {R!H} 0 {5,{S,D,T,B}}, {7,{S,D,T,B}}
7 *2 {R!H} 0 {6,{S,D,T,B}}, {8,S}
8 *3 H 0 {7,S}
""",
  group2 = 
"""
O_rad_out
1 *1 O 1
""",
  group3 = 
"""
Cs_H_out_NDMustO
1 *2 Cs 0 {2,S} {3,S} {4,S} {5,S}
2 *3 H 0 {1,S}
3 O 0 {1,S}
4 {Cs,O} 0 {1,S}
5 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(3.41e6,A_UNITS,"+-",0.0),
                 n=(1.09,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "869",
  short_comment = "Sandeep\'s DFT/CBSB7 level of calculations.",
  long_comment = 
"""
Sandeep and Sumathy paper (submitted to JPCA 2009), intra_H_migration of ROO & HOOQOO.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 256
rate(
  group1 = 
"""
R4H_SDS
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,D}
3 *5 Cd 0 {2,D}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_doubleC
1 *2 Cd 0 {2,S} {3,D}
2 *3 H 0 {1,S}
3 Cd 0 {1,D}
""",
  kf = Arrhenius(A=(1.32E+06,A_UNITS,"+-",0.0),
                 n=(1.6229,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(44.071,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1600),
  rank = 3,
  old_id = "870",
  short_comment = "Sandeep\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sandeep\'s CBS-QB3 calculations.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 257
rate(
  group1 = 
"""
R4H_SDS
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,D}
3 *5 Cd 0 {2,D}, {4,S}
4 *2 {R!H} 0 {3,S}, {5,S}
5 *3 H 0 {4,S}
""",
  group2 = 
"""
Cd_rad_out_double
1 *1 Cd 1 {2,D}
2 Cd 0 {1,D}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.11E+08,A_UNITS,"+-",0.0),
                 n=(1.1915,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(24.7623,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1600),
  rank = 3,
  old_id = "871",
  short_comment = "Sandeep\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sandeep\'s CBS-QB3 calculations.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 258
rate(
  group1 = 
"""
R5H_SMSD
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,{D,T,B}}
3 {R!H} 0 {2,{D,T,B}}, {4,S}
4 *5 Cd 0 {3,S}, {5,D}
5 *2 Cd 0 {4,D}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_singleH
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.19E+05,A_UNITS,"+-",0.0),
                 n=(1.7613,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(38.275,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1600),
  rank = 3,
  old_id = "872",
  short_comment = "Sandeep\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sandeep\'s CBS-QB3 calculations.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 259
rate(
  group1 = 
"""
R5H_DSMS
1 *1 Cd 1 {2,D}
2 *4 Cd 0 {1,D}, {3,S}
3 {R!H} 0 {2,S}, {4,{D,T,B}}
4 *5 {R!H} 0 {3,{D,T,B}}, {5,S}
5 *2 {R!H} 0 {4,S}, {6,S}
6 *3 H 0 {5,S}
""",
  group2 = 
"""
Cd_rad_out_singleH
1 *1 Cd 1 {2,S}
2 H 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.36E+05,A_UNITS,"+-",0.0),
                 n=(1.9199,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.8968,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1600),
  rank = 3,
  old_id = "873",
  short_comment = "Sandeep\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sandeep\'s CBS-QB3 calculations.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 260
rate(
  group1 = 
"""
R3H_SD
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,D}
3 *2 Cd 0 {2,D}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group3 = 
"""
Cd_H_out_singleDe
1 *2 Cd 0 {2,S} {3,S}
2 *3 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  kf = Arrhenius(A=(1.59E+07,A_UNITS,"+-",0.0),
                 n=(1.4638,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(66.3163,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1600),
  rank = 3,
  old_id = "874",
  short_comment = "Sandeep\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sandeep\'s CBS-QB3 calculations.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 261
rate(
  group1 = 
"""
R3H_DS
1 *1 Cd 1 {2,D}
2 *4 Cd 0 {1,D}, {3,S}
3 *2 {R!H} 0 {2,S}, {4,S}
4 *3 H 0 {3,S}
""",
  group2 = 
"""
Cd_rad_out_singleDe
1 *1 Cd 1 {2,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  group3 = 
"""
Cs_H_out_2H
1 *2 Cs 0 {2,S} {3,S} {4,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1283712039,A_UNITS,"+-",0.0),
                 n=(1.0541,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(46.1467,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1600),
  rank = 3,
  old_id = "875",
  short_comment = "Sandeep\'s CBS-QB3 calculations.",
  long_comment = 
"""
Sandeep\'s CBS-QB3 calculations.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


