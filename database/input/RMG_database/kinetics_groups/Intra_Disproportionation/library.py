# encoding: utf-8
header = """
Intra_Disproportionation

Y_birad -> Y


Reverse name: BiradFromMultipleBond

(1) FORM_BOND		{*1,S,*4}
(2) BREAK_BOND		{*2,S,*4}
(3) CHANGE_BOND		{*2,1,*3}
(4) LOSE_RADICAL 	{*1,1}
(5) LOSE_RADICAL	{*3,1}


Generated on 28th July 2010 at 12:55
"""

reaction_family_name = "Intra_Disproportionation"

# determines permitted units for rate expression:
reaction_order = 1

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
//Intra_Disproportionation	
//gmagoon 08/06/09: estimates 1-5 below are based on:
//(1) reference Ea values from Herbinet, Sirjean, Bounaceur, Fournet, Battin-Leclerc, Scacchi, and Marquaire, \"Primary Mechanism of the Thermal Decomposition of Tricyclodecane\", J. Phys. Chem. A, 2006, 110 (39), 11298-11314 DOI: 10.1021/jp0623802
//(2) A (and n) factors from Eq. 1 (with deltan_int = -1) of Warth, Stef, Glaude, Battin-Leclerc, Scacchi, and Come, \"Computer-Aided Derivation of Gas-Phase Oxidation Mechanisms: Application to the Modeling of the Oxidation of n-Butane\", Comb. and Flame 114:81-102 (1998) doi:10.1016/S0010-2180(97)00273-3  
//these are more likely to be overestimates than underestimates																							

Arrhenius_EP

//No	Y_birad		Y_rad		XH_Rrad		Temp		A	N	Alpha	E	DA	DN	DAlpha	DE	Rank	Comments											

"""

# Set some units for all the rates in this file

A_UNITS = "1/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
Y_biCyc3
1 *1 {R!H} 1 {2,{S,D,B}}
2 *2 {R!H} 0 {1,{S,D,B}} {3,{S,D}} {4,S}
3 *3 {R!H} 1 {2,{S,D}}
4 *4 H 0 {2,S}
""",
  group2 = 
"""
Y_rad
1 *1 {R!H} 1
""",
  group3 = 
"""
XH_Rrad
1 *2 {R!H} 0 {2,{S,D}} {3,S}
2 *3 {R!H} 1 {1,{S,D}}
3 *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.66E10,A_UNITS,"+-",0.0),
                 n=(1.0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "1",
  short_comment = "Herbinet et al.(2006) reference Ea and Warth et al.(1998) prefactor with deltan_int=-1",
  long_comment = 
"""
""",
   history = [("2010-07-28","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 2
rate(
  group1 = 
"""
Y_biCyc4
1 *1 {R!H} 1 {5,{S,D,B,T}}
2 *2 {R!H} 0 {5,{S,D,B}} {3,{S,D}} {4,S}
3 *3 {R!H} 1 {2,{S,D}}
4 *4 H 0 {2,S}
5 {R!H} 0 {1,{S,D,B,T}} {2,{S,D,B}}
""",
  group2 = 
"""
Y_rad
1 *1 {R!H} 1
""",
  group3 = 
"""
XH_Rrad
1 *2 {R!H} 0 {2,{S,D}} {3,S}
2 *3 {R!H} 1 {1,{S,D}}
3 *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.66E10,A_UNITS,"+-",0.0),
                 n=(1.0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "2",
  short_comment = "Herbinet et al.(2006) reference Ea and Warth et al.(1998) prefactor with deltan_int=-1",
  long_comment = 
"""
""",
   history = [("2010-07-28","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
Y_biCyc5
Union {Y_biCyc5radEndo,Y_biCyc5radExo}
""",
  group2 = 
"""
Y_rad
1 *1 {R!H} 1
""",
  group3 = 
"""
XH_Rrad
1 *2 {R!H} 0 {2,{S,D}} {3,S}
2 *3 {R!H} 1 {1,{S,D}}
3 *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.66E10,A_UNITS,"+-",0.0),
                 n=(1.0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.75,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "3",
  short_comment = "Herbinet et al.(2006) reference Ea and Warth et al.(1998) prefactor with deltan_int=-1",
  long_comment = 
"""
""",
   history = [("2010-07-28","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
Y_biCyc6
Union {Y_biCyc6radEndo,Y_biCyc6radExo}
""",
  group2 = 
"""
Y_rad
1 *1 {R!H} 1
""",
  group3 = 
"""
XH_Rrad
1 *2 {R!H} 0 {2,{S,D}} {3,S}
2 *3 {R!H} 1 {1,{S,D}}
3 *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.66E10,A_UNITS,"+-",0.0),
                 n=(1.0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.85,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "4",
  short_comment = "Herbinet et al.(2006) reference Ea and Warth et al.(1998) prefactor with deltan_int=-1",
  long_comment = 
"""
""",
   history = [("2010-07-28","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
Y_biCyc7
Union {Y_biCyc7radEndo,Y_biCyc7radExo}
""",
  group2 = 
"""
Y_rad
1 *1 {R!H} 1
""",
  group3 = 
"""
XH_Rrad
1 *2 {R!H} 0 {2,{S,D}} {3,S}
2 *3 {R!H} 1 {1,{S,D}}
3 *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.66E10,A_UNITS,"+-",0.0),
                 n=(1.0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.75,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "5",
  short_comment = "Herbinet et al.(2006) reference Ea and Warth et al.(1998) prefactor with deltan_int=-1",
  long_comment = 
"""
""",
   history = [("2010-07-28","Generated from current RMG library.","rwest@mit.edu")]
)


