# encoding: utf-8
header = """
Oa_R_Recombination

Y_rad + Oa -> YO.


Reverse name: RO_Bond_Dissociation

(1) FORM_BOND		{*1,S,*2}
(2) LOSE_RADICAL 	{*1,1}
(3) LOSE_RADICAL 	{*2,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "Oa_R_Recombination"

# determines permitted units for rate expression:
reaction_order = 2

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// Oa recomb with R.
// JS, July,22, 2003 

// jing, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

//No.	Y_rad	Oa	Temp.		A		N	a	E0	DA	Dn	Da	DE0	Rank

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
Y_rad
1 *1 R 1
""",
  group2 = 
"""
Oa
1 *2 O 2T
""",
  kf = Arrhenius(A=(1E13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "1000",
  short_comment = "",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


