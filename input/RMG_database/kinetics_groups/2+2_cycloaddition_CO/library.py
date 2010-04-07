# encoding: utf-8
header = """
2+2_cycloaddition_CO

CO + doublebond -> four_ring


Reverse name: Four_Ring_Cleavage_CO

(1) CHANGE_BOND		{*1,-1,*2}
(2) CHANGE_BOND		{*3,-1,*4}
(3) FORM_BOND		{*1,S,*3}
(4) FORM_BOND		{*2,S,*4}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "2+2_cycloaddition_CO"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f17b: 2+2-cycloaddition_CO

// jing, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

//f15b_2+2-cycloaddition_CO
// from rate_library_4.txt, Cath, 03/07/28

//No.	CO	doublebond		Temp.		A			N		a		E0		DA		Dn		Da		DE0		Rank	Comments

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
CO
1 *1 CO 0 {2,D}
2 *2 Od 0 {1,D}
""",
  group2 = 
"""
doublebond
Union {mb_CO, mb_OC, mb_CCO, mb_COC}
""",
  kf = Arrhenius(A=(6.92E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(43.72,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "587",
  short_comment = "",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


