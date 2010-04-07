# encoding: utf-8
header = """
2+2_cycloaddition_Cd

db + doublebond -> four_ring


Reverse name: Four_Ring_Cleavage_Cd

(1) CHANGE_BOND		{*1,-1,*2}
(2) CHANGE_BOND		{*3,-1,*4}
(3) FORM_BOND		{*1,S,*3}
(4) FORM_BOND		{*2,S,*4}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "2+2_cycloaddition_Cd"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f17a: 2+2-cycloaddition_Cd

// jing, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

//f15a_2+2-cycloaddition_Cd
// from rate_library_4.txt, Cath, 03/07/28

// Catherina Wijaya thesis, pg 131

// [107] Quick, L. M. Int. J. Chem. Kinet. 1972, 4, 61. 

//No.	db		doublebond		Temp.		A			n		a		E0		DA		Dn		Da		DE0		Rank	Comments

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
db
Union {db_2H,db_HNd,db_HDe,db_Nd2,db_NdDe,db_De2}
""",
  group2 = 
"""
doublebond
Union {mb_db,mb_CO,mb_OC,mb_CCO,mb_COC}
""",
  kf = Arrhenius(A=(6.92E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(43.72,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (723,786),
  rank = 3,
  old_id = "586",
  short_comment = "Quick et al. [107]",
  long_comment = 
"""
[107] Quick, L. M. *Int. J. Chem. Kinet.* 1972, 4, 61. 

C2H4 + C2H4 --> cyclobutane, absolute value measured directly using thermal excitation technique 
and mass spectrometry. Pressure  0.40 - 1.73 bar.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 2
rate(
  group1 = 
"""
db_2H
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
mb_OC
1 *3 Od 0 {2,D}
2 *4 CO 0 {1,D}
""",
  kf = Arrhenius(A=(2.33E+06,A_UNITS,"*/",5.0),
                 n=(1.65,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(54.15,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "6000",
  short_comment = "",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


