# encoding: utf-8
header = """
Intra_R_Add_Endocyclic

Rn -> RnCyclic


Reverse name: Ring_Open_Endo_Cycli_Radical

(1) CHANGE_BOND		{*2,-1,*3}
(2) FORM_BOND		{*1,S,*3}
(3) GAIN_RADICAL	{*2,1}
(4) LOSE_RADICAL 	{*1,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "Intra_R_Add_Endocyclic"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f29: intra radical addition to form endocyclic radical
// original from rate library.txt, CDW 10/20/2002

// JS, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// f29 intra_radical_addition_form_endocyclic_radical
// JS, add 144., a dummy number for debugging, need to be deleleted!!!!!!!!!!!!!!!!

//		Rn			multiplebond_intra			radadd_intra		Temp.		A			n		a	E0		DA	Dn	Da	DE0	Rank
//added by sandeep after performing simulations on 13HXD system
//very crude rough guess
//815.	R6_SSS_D    doublebond_intra_pri_NdNd   radadd_intra_cs		300-1500	1.00E+10	0		0	50.9	0	0	0	0	2


"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
Rn
Union {R3, R4, R5, R6}
""",
  group2 = 
"""
multiplebond_intra
1 *2 {Cd,Ct,CO} 0 {2,{D,T}}
2 *3 {Cd,Ct,Od} 0 {1,{D,T}}
""",
  group3 = 
"""
radadd_intra
1 *1 {R!H} 1
""",
  kf = Arrhenius(A=(1E+08,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "809",
  short_comment = "",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 2
rate(
  group1 = 
"""
R5_SS_D
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 *5 {R!H} 0 {2,S}, {4,S}
4 *2 Cd 0 {3,S}, {5,D}
5 *3 Cd 0 {4,D}
""",
  group2 = 
"""
doublebond_intra_pri_2H
1 *2 Cd 0 {2,D}, {3,S}
2 *3 Cd 0 {1,D}, {4,S}, {5,S}
3 H 0 {1,S}
4 H 0 {2,S}
5 H 0 {2,S}
""",
  group3 = 
"""
radadd_intra_cs2H
1 *1 Cs 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.22E+08,A_UNITS,"+-",0.0),
                 n=(1.05,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.82,E_UNITS,"+-",0.0)
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

# Number 3
rate(
  group1 = 
"""
R6_SSS_D
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S}, {3,S}
3 {R!H} 0 {2,S}, {4,S}
4 *5 {R!H} 0 {3,S}, {5,S}
5 *2 Cd 0 {4,S}, {6,D}
6 *3 Cd 0 {5,D}
""",
  group2 = 
"""
doublebond_intra_pri_2H
1 *2 Cd 0 {2,D}, {3,S}
2 *3 Cd 0 {1,D}, {4,S}, {5,S}
3 H 0 {1,S}
4 H 0 {2,S}
5 H 0 {2,S}
""",
  group3 = 
"""
radadd_intra_cs2H
1 *1 Cs 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.00E+08,A_UNITS,"+-",0.0),
                 n=(0.855,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.9,E_UNITS,"+-",0.0)
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

# Number 4
rate(
  group1 = 
"""
R5_SD_D
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S}, {3,D}
3 *5 Cd 0 {2,D}, {4,S}
4 *2 Cd 0 {3,S}, {5,D}
5 *3 Cd 0 {4,D}
""",
  group2 = 
"""
doublebond_intra_pri_2H
1 *2 Cd 0 {2,D}, {3,S}
2 *3 Cd 0 {1,D}, {4,S}, {5,S}
3 H 0 {1,S}
4 H 0 {2,S}
5 H 0 {2,S}
""",
  group3 = 
"""
radadd_intra_cs2H
1 *1 Cs 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.05E+07,A_UNITS,"+-",0.0),
                 n=(1.192,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(34.9116,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1600),
  rank = 2,
  old_id = "812",
  short_comment = "",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
R3_D
1 *1 {R!H} 1 {2,S}
2 *2 Cd 0 {1,S}, {3,D}
3 *3 Cd 0 {2,D}
""",
  group2 = 
"""
doublebond_intra_pri_HDe
1 *2 Cd 0 {2,D}, {3,S}
2 *3 Cd 0 {1,D}, {4,S}, {5,S}
3 H 0 {1,S}
4 H 0 {2,S}
5 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group3 = 
"""
radadd_intra_cs2H
1 *1 Cs 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.05E+08,A_UNITS,"+-",0.0),
                 n=(1.192,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(54,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1600),
  rank = 5,
  old_id = "813",
  short_comment = "",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
R3_T
1 *1 {R!H} 1 {2,S}
2 *2 Ct 0 {1,S}, {3,T}
3 *3 Ct 0 {2,T}
""",
  group2 = 
"""
triplebond_intra_H
1 *2 Ct 0 {2,T}
2 *3 Ct 0 {1,T}, {3,S}
3 H 0 {2,S}
""",
  group3 = 
"""
radadd_intra_cs2H
1 *1 Cs 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.05E+08,A_UNITS,"+-",0.0),
                 n=(1.192,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(54,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1600),
  rank = 5,
  old_id = "814",
  short_comment = "",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


