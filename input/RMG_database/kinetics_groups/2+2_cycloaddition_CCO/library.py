# encoding: utf-8
header = """
2+2_cycloaddition_CCO

CCO + doublebond -> four_ring


Reverse name: Four_Ring_Cleavage_CCO

(1) CHANGE_BOND		{*1,-1,*2}
(2) CHANGE_BOND		{*3,-1,*4}
(3) FORM_BOND		{*1,S,*3}
(4) FORM_BOND		{*2,S,*4}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "2+2_cycloaddition_CCO"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f17c: 2+2-cycloaddition_CCO

// jing, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

//f15c_2+2-cycloaddition_CCO
// from rate_library_4.txt, Cath, 03/07/28

// Catherina Wijaya thesis, pg 131

// [107] Quick, L. M. Int. J. Chem. Kinet. 1972, 4, 61.

//No.	CCO	doublebond		Temp.		A			n		a		E0		DA		Dn		Da		DE0		Rank	Comments

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
CCO
1 *1 Cd 0 {2,D}
2 *2 Cdd 0 {1,D} {3,D}
3 Od 0 {2,D}
""",
  group2 = 
"""
doublebond
Union {mb_CCO, mb_COC}
""",
  kf = Arrhenius(A=(6.92E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(43.72,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "588",
  short_comment = "Quick et al. [107]",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


