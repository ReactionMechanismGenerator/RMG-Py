# encoding: utf-8
header = """
Birad_recombination

Rn -> Rncycle


Reverse name: Ring_Open

(1) FORM_BOND		{*1,S,*2}
(2) LOSE_RADICAL 	{*1,1}
(3) LOSE_RADICAL 	{*2,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "Birad_recombination"

# determines permitted units for rate expression:
reaction_order = 1

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f07: biradical recombination to form cycle

// jing, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// f07_biradical_recombination_to_form_cycle
// from rate_library_4.txt, Cath, 03/07/28

// Catherina Wijaya Thesis, pg 157.

// [185] Roth, W.R.; Scholz, B.P.; Breuckmann, R.; Jelich, K.; Lennartz, H-W. Chem. Ber. 1982, 115, 1934. 
// [186] Benson, S.W. J. Chem. Phys. 1967, 46, 4920.
// [187] Grimme, W.; Schumachers, L.; Roth, W.R.; Breuckmann, R. Chem. Ber. 1981, 114, 3197. 
// [188] Lewis, K.E.; Steiner, H. J. Chem. Soc. 1964, 3080.
// [x] Sirjean, B.; Glaude, P. A.; Ruiz-Lopez, M. F.; Fournet, R.; J. Phys. Chem. A. 2006, 110, 12693-12704. 

//No.	Rn			Y_rad_out			Ypri_rad_out		Temp.		A			n	a		E0		DA		Dn		Da		DE0		Rank	Comments
//481.	R4_SDS			C_rad_out_2H			Cpri_rad_out_2H		495-549		2.5E+11		0		0		35.7	1.3E+11 0		0		0.5	5	[185] Roth et al.
//483.	R6_SDSDS		C_rad_out_2H			Cpri_rad_out_2H		555-606		2.77E+13	0		0		45.51	*1.2	0		0		0		3	[187] Grimme et al.
//484.	R6_SDSDS		C_rad_out_2H			Cpri_rad_out_2H		390-463		7.14E+11	0		0		29.01	*3.7	0		0		0.58	3	[188] Lewis et al.

"""

# Set some units for all the rates in this file

A_UNITS = "1/s"
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
Y_rad_out
1 *1 {R!H} 1
""",
  group3 = 
"""
Ypri_rad_out
1 *2 {R!H} 1
""",
  kf = Arrhenius(A=(5E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "480",
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
R6_SSSDS
1 *1 {Cs,Cd,CO,Os} 1 {2,S}
2 *3 {Cs,Cd,CO,Os} 0 {1,S}, {3,S}
3 {Cs,Cd,CO,Os} 0 {2,S}, {4,S}
4 Cd 0 {3,S}, {5,D}
5 *4 Cd 0 {4,D}, {6,S}
6 *2 {Cs,Cd,CO,Os} 1 {5,S}
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
Cpri_rad_out_2H
1 *2 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.0E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.80,E_UNITS,"+-",1.0)
                 ),
  temperature_range = (550,650),
  rank = 4,
  old_id = "482",
  short_comment = "[186] Benson et al.",
  long_comment = 
"""

[186] Benson, S.W. J. Chem. Phys. 1967, 46, 4920.

CH2=CHCH(.)CH2CH2CH(.)CH=CH2 --> 4-vinylcyclohexene. (Rxn. -c); arises from birad recombination of resonance isomer: .CH2CH=CHCH2CH2CH(.)CH=CH2

Data are estimated.

***this only considers cis-cis isomer reaction*** cis-trans A prefactor is 50% of the value used here; also, on p. 4923, it is stated that cis trans rate is 5/6 of the overall rate, so maybe the k that should be used is 0.6 of the value currently in place?
Verified by Greg Magoon: Rxn. -d. also looks to be of interest here; whether the rate is high-pressure limit was not investigated, but p. 4922 says that pressures involved were low; DE0 uncertainty added; regarding temperature range, I considered dropping lower temperature limit from 550 K to 400 K (based on p. 4923), but it seems that experiments were performed at or around 600 K, so I will leave it at 550-650 K

Note: after some preliminary confusion on my part, it looks like the existing groups are correct (the correct resonance form for the CH2 radical is taken into account with the Ypri_rad_out (i.e. Cpri_rad_out_H2)), but arguably, another, a more-specific group (C_rad_out_H2/OneDe and Cpri_rad_out_H2/OneDe) should be specified to account for delocalizing group at this site
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
R4_SSS
1 *1 {Cs,Cd,CO,Os} 1 {2,S}
2 *3 {Cs,Cd,CO,Os} 0 {1,S}, {3,S}
3 *4 {Cs,Cd,CO,Os} 0 {2,S}, {4,S}
4 *2 {Cs,Cd,CO,Os} 1 {3,S}
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
Cpri_rad_out_2H
1 *2 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.62E+12,A_UNITS,"+-",0.0),
                 n=(-0.305,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.98,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "485",
  short_comment = "[x] Sirjean et al.",
  long_comment = 
"""

[x] Sirjean, B.; Glaude, P. A.; Ruiz-Lopez, M. F.; Fournet, R.; J. Phys. Chem. A. 2006, 110, 12693-12704. 
http://dx.doi.org/10.1021/jp0651081
.CH2CH2CH2CH2CH2. -> cyclopentane (k4-1 in Scheme 5/Table 7)

TST calculation

Added by Greg Magoon: Stated pressure is 1 atm, but I believe they are actually calculating the high-pressure limit rate constant
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
R5_SSSS
1 *1 {Cs,Cd,CO,Os} 1 {2,S}
2 *3 {Cs,Cd,CO,Os} 0 {1,S}, {3,S}
3 {Cs,Cd,CO,Os} 0 {2,S}, {4,S}
4 *4 {Cs,Cd,CO,Os} 0 {3,S}, {5,S}
5 *2 {Cs,Cd,CO,Os} 1 {4,S}
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
Cpri_rad_out_2H
1 *2 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.76E+9,A_UNITS,"+-",0.0),
                 n=(0.311,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "486",
  short_comment = "[x] Sirjean et al.",
  long_comment = 
"""

[x] Sirjean, B.; Glaude, P. A.; Ruiz-Lopez, M. F.; Fournet, R.; J. Phys. Chem. A. 2006, 110, 12693-12704. 
http://dx.doi.org/10.1021/jp0651081
.CH2CH2CH2CH2CH2CH2. -> cyclohexane (k5-1+k5-2 in Scheme 7/Table 10) (includes formation of both boat and chair conformations)

TST calculation

Added by Greg Magoon: Stated pressure is 1 atm, but I believe they are actually calculating the high-pressure limit rate constant; the rate constant added here was found my performing least squares fit for log(ktot) from 600-2000 K (where ktot is the sum of k5-1 and k5-2)

Note: Recent experimental/RRKM study by Kiefer, Gupte, Harding, and Klippenstein (http://www.combustion.org.uk/ECM_2009/P810069.pdf) (stated uncertainty is +/- 30%) appears to agree with Sirjean et al. results, but they only report forward rate constant
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
R6_SSSSS
1 *1 {Cs,Cd,CO,Os} 1 {2,S}
2 *3 {Cs,Cd,CO,Os} 0 {1,S}, {3,S}
3 {Cs,Cd,CO,Os} 0 {2,S}, {4,S}
4 {Cs,Cd,CO,Os} 0 {3,S}, {5,S}
5 *4 {Cs,Cd,CO,Os} 0 {4,S}, {6,S}
6 *2 {Cs,Cd,CO,Os} 1 {5,S}
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
Cpri_rad_out_2H
1 *2 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.21E+10,A_UNITS,"+-",0.0),
                 n=(0.137,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.12,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,2000),
  rank = 3,
  old_id = "486",
  short_comment = "[x] Sirjean et al.",
  long_comment = 
"""

[x] Sirjean, B.; Glaude, P. A.; Ruiz-Lopez, M. F.; Fournet, R.; J. Phys. Chem. A. 2006, 110, 12693-12704. 
http://dx.doi.org/10.1021/jp0651081
.CH2CH2CH2CH2CH2CH2. -> cyclohexane (k5-1+k5-2 in Scheme 7/Table 10) (includes formation of both boat and chair conformations)

TST calculation

Added by Greg Magoon: Stated pressure is 1 atm, but I believe they are actually calculating the high-pressure limit rate constant; the rate constant added here was found my performing least squares fit for log(ktot) from 600-2000 K (where ktot is the sum of k5-1 and k5-2)

Note: Recent experimental/RRKM study by Kiefer, Gupte, Harding, and Klippenstein (http://www.combustion.org.uk/ECM_2009/P810069.pdf) (stated uncertainty is +/- 30%) appears to agree with Sirjean et al. results, but they only report forward rate constant
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


