# encoding: utf-8
header = """
Intra_R_Add_Exocyclic

Rn -> RnCycle


Reverse name: Ring_Open_Exo_Cycli_Radical

(1) CHANGE_BOND		{*2,-1,*3}
(2) FORM_BOND		{*1,S,*2}
(2) LOSE_RADICAL 	{*1,1}
(3) GAIN_RADICAL 	{*3,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "Intra_R_Add_Exocyclic"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f27: intra radical addition to form exocyclic radical
// original from rate library.txt, CDW 10/20/2002

// JS, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// f27 intra_radical_addition_form_exocyclic_radical

//		Rn			multiplebond_intra			radadd_intra		Temp.		A			n	a		E0		DA		Dn		Da		DE0		Rank	Comments


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
Union {R4, R5, R6, R7}
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
  kf = Arrhenius(A=(1E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "807",
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
R6_SSS_D
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S} {3,S}
3 {R!H} 0 {2,S} {4,S}
4 *5 {R!H} 0 {3,S} {5,S}
5 *2 Cd 0 {4,S} {6,D}
6 *3 Cd 0 {5,D}
""",
  group2 = 
"""
doublebond_intra_2H_pri
1 *2 Cd 0 {2,D} {3,S}
2 *3 Cd 0 {1,D} {4,S} {5,S}
3 H 0 {1,S}
4 H 0 {2,S}
5 H 0 {2,S}
""",
  group3 = 
"""
radadd_intra_cs2H
1 *1 Cs 1 {2,S} {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.51E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.85,E_UNITS,"+-",0.0)
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

# Number 3
rate(
  group1 = 
"""
R4_S_D
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S} {3,S}
3 *2 Cd 0 {2,S} {4,D}
4 *3 Cd 0 {3,D}
""",
  group2 = 
"""
doublebond_intra_2H_pri
1 *2 Cd 0 {2,D} {3,S}
2 *3 Cd 0 {1,D} {4,S} {5,S}
3 H 0 {1,S}
4 H 0 {2,S}
5 H 0 {2,S}
""",
  group3 = 
"""
radadd_intra_csHDe
1 *1 Cs 1 {2,S} {3,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  kf = Arrhenius(A=(1.00E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(25.85,E_UNITS,"+-",0.0)
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

# Number 4
rate(
  group1 = 
"""
R6_SMS_D
1 *1 {R!H} 1 {2,S}
2 *4 {Cd,Ct,Cb} 0 {1,S} {3,{D,T,B}}
3 {Cd,Ct,Cb} 0 {2,{D,T,B}} {4,S}
4 *5 {R!H} 0 {3,S} {5,S}
5 *2 Cd 0 {4,S} {6,D}
6 *3 Cd 0 {5,D}
""",
  group2 = 
"""
doublebond_intra
1 *2 Cd 0 {2,D}
2 *3 Cd 0 {1,D}
""",
  group3 = 
"""
radadd_intra_cs
1 *1 Cs 1
""",
  kf = Arrhenius(A=(1E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(46.85,E_UNITS,"+-",0.0)
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

# Number 5
rate(
  group1 = 
"""
R6_SSM_D
1 *1 {R!H} 1 {2,S}
2 *4 {R!H} 0 {1,S} {3,S}
3 {Cd,Ct,Cb} 0 {2,S} {4,{D,T,B}}
4 *5 {Cd,Ct,Cb} 0 {3,{D,T,B}} {5,S}
5 *2 Cd 0 {4,S} {6,D}
6 *3 Cd 0 {5,D}
""",
  group2 = 
"""
doublebond_intra
1 *2 Cd 0 {2,D}
2 *3 Cd 0 {1,D}
""",
  group3 = 
"""
radadd_intra_cs
1 *1 Cs 1
""",
  kf = Arrhenius(A=(1E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(46.85,E_UNITS,"+-",0.0)
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

# Number 6
rate(
  group1 = 
"""
R5_SD_D
1 *1 {R!H} 1 {2,S}
2 *4 Cd 0 {1,S} {3,D}
3 *5 Cd 0 {2,D} {4,S}
4 *2 Cd 0 {3,S} {5,D}
5 *3 Cd 0 {4,D}
""",
  group2 = 
"""
doublebond_intra_HNd_pri
1 *2 Cd 0 {2,D} {3,S}
2 *3 Cd 0 {1,D} {4,S} {5,S}
3 H 0 {1,S}
4 H 0 {2,S}
5 {Cs,O} 0 {2,S}
""",
  group3 = 
"""
radadd_intra_csHNd
1 *1 Cs 1 {2,S} {3,S}
2 H 0 {1,S}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(46.85,E_UNITS,"+-",0.0)
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


