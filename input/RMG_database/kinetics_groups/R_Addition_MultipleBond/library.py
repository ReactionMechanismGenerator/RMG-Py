# encoding: utf-8
header = """
R_Addition_MultipleBond

XZ + Y_rad_birad -> YXZ.


Reverse name: Beta_Scission

(1) CHANGE_BOND		{*1,-1,*2}
(2) FORM_BOND		{*1,S,*3}
(2) GAIN_RADICAL 	{*2,1}
(3) LOSE_RADICAL 	{*3,1}


Generated on 7th April 2010 at 17:08
Generated on 22nd June 2010 at 12:58
"""

reaction_family_name = "R_Addition_MultipleBond"

# determines permitted units for rate expression:
reaction_order = 2

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f02: Radical addition to multiple bond
// original from rate_library_4.txt, Cath, 03/07/28

// JS, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

//f02_Radical_Addition
//No.		XZ				Y_rad			Temp.		A			N		a	E0		DA		Dn		Da	DE0		Rank	Comments
//416 Added by sandeep taken from the stufy of 1,3-Hexadiene dont at CBS-QB3 level

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
XZ
Union {CZ,OCO,OCddO,OSi,OSiddO}
""",
  group2 = 
"""
Y_rad_birad
Union {Y_rad, Y_birad}
""",
  kf = Arrhenius(A=(1E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "269",
  short_comment = "Default",
  long_comment = 
"""
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 2
rate(
  group1 = 
"""
Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}
3    H  0 {1,S}
4    H  0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "281",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
in his reaction type 3. Based on the recommendations of
[146] Allara, D.L.; Shaw, R. J Phys. Chem. Ref. Data 1980,9,523.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
Cd/H/Nd
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}
3       H   0 {1,S}
4  {Cs,O}   0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "282",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
in his reaction type 3. Based on the recommendations of
[146] Allara, D.L.; Shaw, R. J Phys. Chem. Ref. Data 1980,9,523.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
Cd/Nd2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "283",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
in his reaction type 3. Based on the recommendations of
[146] Allara, D.L.; Shaw, R. J Phys. Chem. Ref. Data 1980,9,523.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}
3    H  0 {1,S}
4    H  0 {1,S}
""",
  group2 = 
"""
Cs_rad
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 R 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(8.5E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "284",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
in his reaction type 3. Based on the recommendations of
[146] Allara, D.L.; Shaw, R. J Phys. Chem. Ref. Data 1980,9,523.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
Cd/H/Nd
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}
3       H   0 {1,S}
4  {Cs,O}   0 {1,S}
""",
  group2 = 
"""
Cs_rad
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 R 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(8.5E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "285",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
in his reaction type 3. Based on the recommendations of
[146] Allara, D.L.; Shaw, R. J Phys. Chem. Ref. Data 1980,9,523.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
Cd/Nd2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
Cs_rad
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 R 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(8.5E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "286",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
in his reaction type 3. Based on the recommendations of
[146] Allara, D.L.; Shaw, R. J Phys. Chem. Ref. Data 1980,9,523.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}
3    H  0 {1,S}
4    H  0 {1,S}
""",
  group2 = 
"""
O_rad/NonDe
1 *3 O 1 {2,S}
2 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.0E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "287",
  short_comment = "Curran et al. [8] in his reaction type 20. Based on recommendations of Chen and Bozzelli [57]",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
in his reaction type 20. Based on the recommendations of
[146] Allara, D.L.; Shaw, R. J Phys. Chem. Ref. Data 1980,9,523.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 9
rate(
  group1 = 
"""
Cd/H/Nd
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}
3       H   0 {1,S}
4  {Cs,O}   0 {1,S}
""",
  group2 = 
"""
O_rad/NonDe
1 *3 O 1 {2,S}
2 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.0E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "288",
  short_comment = "Curran et al. [8] in his reaction type 20. Based on recommendations of Chen and Bozzelli [57]",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
in his reaction type 20. Based on the recommendations of
[146] Allara, D.L.; Shaw, R. J Phys. Chem. Ref. Data 1980,9,523.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 10
rate(
  group1 = 
"""
Cd/Nd2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDe
1 *3 O 1 {2,S}
2 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.0E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "289",
  short_comment = "Curran et al. [8] in his reaction type 20. Based on recommendations of Chen and Bozzelli [57]",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
in his reaction type 20. Based on the recommendations of
[146] Allara, D.L.; Shaw, R. J Phys. Chem. Ref. Data 1980,9,523.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 11
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.985E+09,A_UNITS,"*/",2.0),
                 n=(1.28,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.29,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (200,1100),
  rank = 4,
  old_id = "290",
  short_comment = "Baulch et al. [94] literature review.",
  long_comment = 
"""
[94] Baulch,D.L.; Cobos,C.J.;Cox,R.A;Frank,P.;Hayman,G.;Just,T.;Kerr,J.A.;Murells,T.;Philling,M.J.;Troe,J.;Walker,R.W.; Warnatz, J. J Phys Chem. Ref. Data 1994,23,847.
literature review. C2H4 + H --> C2H5. C.D.W. divided rate expression by 2, to get rate of addition per site 
pg.916-920: Discussion on evaluated data

H+C2H4(+m) --> C2H5(+m): \"The analysis of the rxn is based on theoretical fall-off

curves and strong collision low pressure rate coefficients which were calculated
using a rxn threshold of 154.78 kJ/mol.\"  The rate coefficient stored in RMG
is the high-pressure limit, k_inf.
MRH 31-Aug-2009
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 12
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.655E+11,A_UNITS,"*/",1.3),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.71,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "291",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J. Phys. Chem. Ref. Data 1986, 15, 1087. 
literature review. C2H4 + CH3 --> n-C3H7. C.D.W. divided rate expression by 2, to get rate of addition per site
pg. 1191: Discussion on evaluated data

Entry 18,16 (b)

Recommended data is from other Review paper by Kerr and Parsonage (1972)

MRH 28-Aug-2009
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 13
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.99E+03,A_UNITS,"+-",0.0),
                 n=(2.44,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.37,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,1500),
  rank = 2,
  old_id = "292",
  short_comment = "Knyazev et al. [147]",
  long_comment = 
"""
[147] Knyazev,V.D.;Slagle,I.R. J Phys. Chem. 1996 100, 5318.
Pressure up to 10 atm. Excitation; thermal, analysis: mass spectrometry. C2H4 + C2H5--> n-C4H9. C.D.W. divided rate expression by 2, to get rate of addtion per site
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 14
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/O
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 O 0 {1,S}
""",
  kf = Arrhenius(A=(2.41E+10,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.96,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "293",
  short_comment = "Tsang et al. [90] literature review.",
  long_comment = 
"""
[90] Tsang,W.J. Phys. Chem. Ref. Data 1987,16,471.
literature review. C2H4+ CH2OH --> CH2CH2CH2OH C.D.W. divided rate expression by 2, to get rate of addition per site
pg. 502: Discussion on evaluated data

Entry 39,18 (a): No data available at the time.  Author suggests rate coefficient expression

of 8.0x10^-14 * exp(-3500/T) cm3/molecule/s noting rates of alkyl radical addition
to ethylene are similar (Kerr, J.A., Trotman-Dickenson, A.F.)
MRH 30-Aug-2009
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 15
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 Cd 1 {2,D}, {3,S}
2 Cd 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.0E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.01,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (1260,1310),
  rank = 5,
  old_id = "294",
  short_comment = "Weissman and Benson [148] Estimated values.",
  long_comment = 
"""
[148] Weissman and Benson. Estimated values. Activation energy is a lower limit. Pressure 1.00 atm. 
C2H4 + C2H3 --> CH2=CHCH2CH2 C.D.W. divided rate expression by 2, to get rate of addition per site
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 16
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.71E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 4,
  old_id = "295",
  short_comment = "Tsang et al. [89] literature review.",
  long_comment = 
"""
[89] Tsang et al. Literature Review.  
C2H4 + OH --> CH2CH2OH  C.D.W. divided rate expression by 2, to get rate of addition per site

pg. 1189: Discussion on evaluated data (in theory)

Online reference does not have pages 1188-1189; pages 1198-1199 come between
pages 1187&1190 and between 1197&1200
Following discussion is only based on table (pg. 1097) that summarizes all evaluated

data in the reference
Entry 18,6 (b)

Table states rxn is pressure-dependent: C2H4+OH(+M)=C2H4OH(+M)

Only data available in table is k=9.0x10^-12
MRH 28-Aug-2009
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 17
rate(
  group1 = 
"""
Cd/H2_Cd/H/Nd
1 *1 Cd     0 {2,D}, {3,S}, {4,S}
2 *2 Cd     0 {1,D}, {5,S}, {6,S}
3    H      0 {1,S}
4    H      0 {1,S}
5    H      0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.3E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.56,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (500,2500),
  rank = 4,
  old_id = "296",
  short_comment = "Tsang [149] experiments and limited review.",
  long_comment = 
"""
[149] Tsang experiments and limited review. CH3CH=CH2 + H --> iso-C3H7
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 18
rate(
  group1 = 
"""
Cd/H2_Cd/H/Nd
1 *1 Cd     0 {2,D}, {3,S}, {4,S}
2 *2 Cd     0 {1,D}, {5,S}, {6,S}
3    H      0 {1,S}
4    H      0 {1,S}
5    H      0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.28E+05,A_UNITS,"+-",0.0),
                 n=(2.28,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.60,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,1500),
  rank = 4,
  old_id = "297",
  short_comment = "Knyazev et al. [150]",
  long_comment = 
"""
[150] Knayzev et al. Data derived from fitting to a complex mechanism. Pressure up to 10 atm. Excitation : flash photolysis, analysis : mass spectrometry
CH3CH=CH2 + CH3 --> sec-C4H9 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 19
rate(
  group1 = 
"""
Cd/H2_Cd/H/Nd
1 *1 Cd     0 {2,D}, {3,S}, {4,S}
2 *2 Cd     0 {1,D}, {5,S}, {6,S}
3    H      0 {1,S}
4    H      0 {1,S}
5    H      0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.69E+11,A_UNITS,"*/",1.4),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.41,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "298",
  short_comment = "Tsang [93] literature review.",
  long_comment = 
"""
[93] Tsang literature review. CH3CH=CH2 + CH3 --> sec-C4H9 
pg.237-239: Discussion on evaluated data

Entry 46,16(a): Recommended rate coefficient is that reported by Kerr and Parsonage (1972).

Author notes that rxn is pressure dependent and lists fall-off ratios and
collision efficiencies; these are not stored in RMG.
MRH 31-Aug-2009
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 20
rate(
  group1 = 
"""
Cd/H2_Cd/H/Nd
1 *1 Cd     0 {2,D}, {3,S}, {4,S}
2 *2 Cd     0 {1,D}, {5,S}, {6,S}
3    H      0 {1,S}
4    H      0 {1,S}
5    H      0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cd
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  kf = Arrhenius(A=(3.55E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.89,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (762,811),
  rank = 4,
  old_id = "299",
  short_comment = "Barbe et al. [151] Data are estimated.",
  long_comment = 
"""
[151] Barbe et al. Data is estimated. Pressure 0.04-0.26 atm. CH3CH=CH2 + .CH2CH=CH2 --> CH3CH(.)CH2CH2CH=CH2
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 21
rate(
  group1 = 
"""
Cd/H2_Cd/H/Nd
1 *1 Cd     0 {2,D}, {3,S}, {4,S}
2 *2 Cd     0 {1,D}, {5,S}, {6,S}
3    H      0 {1,S}
4    H      0 {1,S}
5    H      0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.07E+09,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.88,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "300",
  short_comment = "Tsang [93] literature review.",
  long_comment = 
"""
[93] Tsang literature review. CH3CH=CH2 + tert-C4H9 --> (CH3)3CCH2CH(.)CH3
pg.247: Discussion on evaluated data

Entry 46,44(terminal): Recommended rate coefficient is based on summary of data on alkyl

radical addition to olefins (Kerr and Parsonage, 1972).
MRH 31-Aug-2009
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 22
rate(
  group1 = 
"""
Cd/H2_Cd/H/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       H       0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.155E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.49,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (743,772),
  rank = 5,
  old_id = "301",
  short_comment = "Perrin et al. [152] Data are estimated.",
  long_comment = 
"""
[152] Perrin et al. Data is estimated. Pressure 0.01-0.13 atm. 
CH2=CHCH=CH2 + .CH3 --> CH2CH=CHCH2CH3 C.D.W. divied rate expression by 2, to get rate of addition per site.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 23
rate(
  group1 = 
"""
Cd/H2_Cd/Nd2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4       H   0 {1,S}
5    {Cs,O} 0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(6.21E+12,A_UNITS,"+-",0.0),
                 n=(0.25,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.46,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (712,779),
  rank = 4,
  old_id = "302",
  short_comment = "Knyazev et al. [153]",
  long_comment = 
"""
[153] Knayzev et al. Pressure ~ 0.01 atm. Excitation : thermal, analysis : GC Iso-C4H8 + CH3 --> (CH3)2CCH2CH3
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 24
rate(
  group1 = 
"""
Cd/H2_Cd/Nd2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4       H   0 {1,S}
5    {Cs,O} 0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.51E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.70,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (391,449),
  rank = 4,
  old_id = "303",
  short_comment = "Seres et al. [154] Data derived from fitting a complex mechanism.",
  long_comment = 
"""
[303] Seres et al. Data derived from fitting to a complex mechanism. Excitation : thermal, analysis : GC Iso-C4H8 + CH3 --> (CH3)2CCH2CH3
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 25
rate(
  group1 = 
"""
Cd/H/Nd_Cd/H2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4   {Cs,O}  0 {1,S}
5       H   0 {2,S}
6       H   0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.3E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.26,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (500,2500),
  rank = 4,
  old_id = "304",
  short_comment = "Tsang [149] experiments and limited review.",
  long_comment = 
"""
[149] Tsang experiments and limited review. CH3CH=CH2 + H --> n-C3H7
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 26
rate(
  group1 = 
"""
Cd/H/Nd_Cd/H2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4   {Cs,O}  0 {1,S}
5       H   0 {2,S}
6       H   0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.0E+04,A_UNITS,"+-",0.0),
                 n=(2.57,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.71,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,1500),
  rank = 4,
  old_id = "305",
  short_comment = "Knyazev et al. [147]",
  long_comment = 
"""
[147] Knyazev et al. Pressure up to 10 atm. Excitation : thermal, analysis : mass spectrometry. 
CH3CH=CH2 + CH3 --> iso-C4H9
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 27
rate(
  group1 = 
"""
Cd/H/Nd_Cd/H2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4   {Cs,O}  0 {1,S}
5       H   0 {2,S}
6       H   0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.64E+10,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.01,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "306",
  short_comment = "Tsang [93] literature review.",
  long_comment = 
"""
[93] literature review. CH3CH=CH2 + CH3 --> iso-C4H9
pg.237-239: Discussion on evaluated data

Entry 46,16(b): Recommended rate coefficient is from reverse rate and equilibrium constant.

Author notes that rxn is pressure dependent and lists fall-off ratios and
collision efficiencies; these are not stored in RMG.
MRH 31-Aug-2009
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 28
rate(
  group1 = 
"""
Cd/Nd2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.23E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.59,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (560,650),
  rank = 5,
  old_id = "307",
  short_comment = "Slagle et al. [155] Data derived from detailed balance/reverse rate.",
  long_comment = 
"""
[155] Slagle et al. Data deriver from detailed balance/reverse rate. Pressure ~ 0.01 atm. 
Iso-C4H8 + .CH3 --> (CH3)3CCH2
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 29
rate(
  group1 = 
"""
Cd/H2_Ca
1 *1 Cd  0 {2,D}, {3,S}, {4,S}
2 *2 Cdd  0 {1,D}, {5,D}
3 H  0 {1,S}
4 H  0 {1,S}
5 C 0 {2,D}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "308",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran et al. in his reaction type 3. Based on recommendations of Allara and Shaw. [146] 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 30
rate(
  group1 = 
"""
Cd/H/Nd_Ca
1 *1 Cd  0 {2,D}, {3,S}, {4,S}
2 *2 Cdd  0 {1,D}, {5,D}
3 H  0 {1,S}
4 {Cs,O} 0 {1,S}
5 C 0 {2,D}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "309",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran et al. in his reaction type 3. Based on recommendations of Allara and Shaw. [146] 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 31
rate(
  group1 = 
"""
Cd/Nd2_Ca
1 *1 Cd  0 {2,D}, {3,S}, {4,S}
2 *2 Cdd  0 {1,D}, {5,D}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
5 C 0 {2,D}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "310",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran et al. in his reaction type 3. Based on recommendations of Allara and Shaw. [146] 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 32
rate(
  group1 = 
"""
Cd/H2_Ca
1 *1 Cd  0 {2,D}, {3,S}, {4,S}
2 *2 Cdd  0 {1,D}, {5,D}
3 H  0 {1,S}
4 H  0 {1,S}
5 C 0 {2,D}
""",
  group2 = 
"""
Cs_rad
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 R 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(8.5E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "311",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran et al. in his reaction type 3. Based on recommendations of Allara and Shaw. [146] 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 33
rate(
  group1 = 
"""
Cd/H/Nd_Ca
1 *1 Cd  0 {2,D}, {3,S}, {4,S}
2 *2 Cdd  0 {1,D}, {5,D}
3 H  0 {1,S}
4 {Cs,O} 0 {1,S}
5 C 0 {2,D}
""",
  group2 = 
"""
Cs_rad
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 R 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(8.5E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "312",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
[8] Curran et al. in his reaction type 3. Based on recommendations of Allara and Shaw. [146] 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 34
rate(
  group1 = 
"""
Cd/Nd2_Ca
1 *1 Cd  0 {2,D}, {3,S}, {4,S}
2 *2 Cdd  0 {1,D}, {5,D}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
5 C 0 {2,D}
""",
  group2 = 
"""
Cs_rad
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 R 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(8.5E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "313",
  short_comment = "Curran et al. [8] in his reaction type 3. Based on recommendations of Allara and Shaw [146]",
  long_comment = 
"""
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 35
rate(
  group1 = 
"""
Cd/H2_Ca
1 *1 Cd  0 {2,D}, {3,S}, {4,S}
2 *2 Cdd  0 {1,D}, {5,D}
3 H  0 {1,S}
4 H  0 {1,S}
5 C 0 {2,D}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.75E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.84,E_UNITS,"+-",0.2)
                 ),
  temperature_range = (573,595),
  rank = 5,
  old_id = "314",
  short_comment = "Scherzer et al. [156] Data derived from fitting a complex mechanism.",
  long_comment = 
"""
[156] Scherzer et al. Data derived from fitting to a complex mechanism. Pressure 0.04 atm. Excitation: thermal, analysis: GC.
CH2=C=CH2 + .CH3 --> CH3CH2C=CH2
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 36
rate(
  group1 = 
"""
Ca_Cd/H2
1 *1 Cdd 0 {2,D}, {3,D}
2 *2 Cd 0 {1,D} {4,S} {5,S}
3 C 0 {1,D}
4 H 0 {2,S}
5 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.2E+11,A_UNITS,"*/",2.5),
                 n=(0.69,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (350,1200),
  rank = 4,
  old_id = "315",
  short_comment = "Tsang et al. [157]",
  long_comment = 
"""
[157] Tsang et al. Absolute Value Measured directly. Pressure 2 - 7 atm. Excitation: thermal, analysis : GC. 
CH2=C=CH2 + H --> .CH2CH=CH2
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 37
rate(
  group1 = 
"""
Ca_Cd/H2
1 *1 Cdd 0 {2,D}, {3,D}
2 *2 Cd 0 {1,D} {4,S} {5,S}
3 C 0 {1,D}
4 H 0 {2,S}
5 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.58E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.97,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (996,1180),
  rank = 5,
  old_id = "316",
  short_comment = "Tsang [158] Data is estimated.",
  long_comment = 
"""
[158] Tsang. Data is estimated. Pressure 1.50-5.00 atm. CH2=C=CH2 + CH3 --> CH2C(CH3)=CH2
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 38
rate(
  group1 = 
"""
CO_O
1 *1 CO 0 {2,D}
2 *2 Od 0 {1,D}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.0E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "317",
  short_comment = "Curran et al. [8] in his reaction type 18.",
  long_comment = 
"""
[8] Curran et al. In his reaction type 18. 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 39
rate(
  group1 = 
"""
CO_O
1 *1 CO 0 {2,D}
2 *2 Od 0 {1,D}
""",
  group2 = 
"""
Cs_rad
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 R 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(1.0E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "318",
  short_comment = "Curran et al. [8] in his reaction type 18.",
  long_comment = 
"""
[8] Curran et al. In his reaction type 18. 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 40
rate(
  group1 = 
"""
CO_O
1 *1 CO 0 {2,D}
2 *2 Od 0 {1,D}
""",
  group2 = 
"""
CO_pri_rad
1 *3 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.2E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.560,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 4,
  old_id = "319",
  short_comment = "Bozzelli et al. [144] based on CH3 addition to CO (Anastasi and Maw)",
  long_comment = 
"""
[144] Bozzelli et al. Based upon CH3 addition to CO (Anastasi and Maw)
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 41
rate(
  group1 = 
"""
CO_O
1 *1 CO 0 {2,D}
2 *2 Od 0 {1,D}
""",
  group2 = 
"""
O_rad/OneDe
1 *3 O 1 {2,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  kf = Arrhenius(A=(1.3E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "320",
  short_comment = "Curran esitmation [159] in DME oxidation modeling for ketohydroperoxide decomposition.",
  long_comment = 
"""
[159] Curran et al. His estimation in DME oxidation modeling for ketohydroperoxide decomposition. 
H2CO + HCO2. (formic acid radical) --> +  .OCH2OCHO (ester) (Rxn. 338, p. 234)

Verified by Greg Magoon; it is not immediately clear whether this rate constant is for high pressure limit, but based on other references to high pressure limit in the paper, I suspect that it is a high pressure limit value; also, note that CO_O group is used for H2CO...MRH and I have interpreted CO_O as referring to any carbonyl group
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 42
rate(
  group1 = 
"""
CO/H2_O
1 *1 CO  0 {2,D}, {3,S}, {4,S}
2 *2 Od  0 {1,D}
3 H  0 {1,S}
4 H  0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(7.94E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.70,E_UNITS,"+-",0.47)
                 ),
  temperature_range = (333,363),
  rank = 5,
  old_id = "321",
  short_comment = "Knoll et al. [160] Data derived from fitting a complex mechanism.",
  long_comment = 
"""
[160] Knoll et al. Data derived from fitting to a complex mechanism. Pressure 0.08 atm. Excitation : direct photolysis, analysis : mass spectrometry.
N-C3H7 + C2HO --> N-C4H9O 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 43
rate(
  group1 = 
"""
CO/Nd2_O
1 *1 CO  0 {2,D}, {3,S}, {4,S}
2 *2 Od  0 {1,D}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.16E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.51,E_UNITS,"+-",1.15)
                 ),
  temperature_range = (413,563),
  rank = 3,
  old_id = "322",
  short_comment = "Knoll et al. [161]",
  long_comment = 
"""
[161] Knoll et al. Absolute value measured directly. Pressure 0.28 - 1.17 atm. Excitation : thermal, analysis : mass spectrometry. 
(CH3)2CO + .CH3 --> (CH3)3CO
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 44
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.75E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.42,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2000),
  rank = 4,
  old_id = "323",
  short_comment = "Warnatz [134] literature review.",
  long_comment = 
"""
[134] Warnatz literature review. C.D.W divided rate expression by 2, to get rate of addition per site.
C2H2 + H --> C2H3
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 45
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.875E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.77,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (370,478),
  rank = 4,
  old_id = "324",
  short_comment = "E.W. Diau and M.C. Lin [162] RRK(M) extrapolation.",
  long_comment = 
"""
[162] E.W.Diau and M.C.Lin. RRK(M) extrapolation. C.D.W divided rate expression by 2, to get rate of addition per site. 
C2H2 + CH3 --> CH3CH=CH
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 46
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.505E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.99,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (373,473),
  rank = 4,
  old_id = "325",
  short_comment = "Kerr et al. [163] literature review.",
  long_comment = 
"""
[163] Kerr et al. literature review. Pressure 0.03-0.20 atm. C.D.W divided rate expression by 2, to get rate of addition per site.
C2H2 + .C2H5 --> CH3CH2CH=CH 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 47
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cd
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  kf = Arrhenius(A=(1.595E+10,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.96,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "326",
  short_comment = "Tsang [93] literature review.",
  long_comment = 
"""
[93] Tsang et al. literature review. Pressure 0.03-0.20 atm. C.D.W divided rate expression by 2, to get rate of addition per site.
C2H2 + .CH2CH=CH2 --> CHCH2CH=CH 

pg.263: Discussion on evaluated data

Entry 47,20(a): Recommended rate coefficient is estimated from the addition of alkyl

radicals to C2H2.  Author notes that this could be used as an upper limit for
cyclopentadiene formation.
MRH 31-Aug-2009
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 48
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.505E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.90,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (363,577),
  rank = 5,
  old_id = "327",
  short_comment = "Kerr et al. [163] literature review.",
  long_comment = 
"""
[163] Kerr et al. literature review. Pressure 0.07-0.13 atm. C.D.W divided rate expression by 2, to get rate of addition per site.
C2H2 + Iso-C3H7 --> (CH3)2CHCH=CH
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 49
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.505E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.31,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (373,493),
  rank = 4,
  old_id = "328",
  short_comment = "Dominguez et al. [164] Data derived from fitting a complex mechanism.",
  long_comment = 
"""
[164] Dominguez et al. Data derived from fitting to a complex mechanism. Pressure 0.01-0.32 atm. Excitation : direct photolysis, analysis : GC. 
C2H2 + Tert-C4H9 --> (CH3)3CCH=CH C.D.W divided rate expression by 2, to get rate of addition per site.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 50
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 Cd 1 {2,D}, {3,S}
2 Cd 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.255E+05,A_UNITS,"+-",0.0),
                 n=(1.90,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.11,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "329",
  short_comment = "Weissman et al. [121] Transition state theory.",
  long_comment = 
"""
[121] Weissman et al. Transition state theory. C.D.W divided rate expression by 2, to get rate of addition per site.	
C2H2 + C2H3 --> CH2=CHCH=CH.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 51
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 Cd 1 {2,D}, {3,S}
2 Cd 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.155E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (700,1300),
  rank = 3,
  old_id = "330",
  short_comment = "Duran et al. [165] Ab initio.",
  long_comment = 
"""
[165] Duran et al. Ab initio. C.D.W divided rate expression by 2, to get rate of addition per site.
C2H2 + C2H3 --> CH2=CHCH=CH. (Rxn. -5?)

Verified by Greg Magoon: note: NIST seems to have values (http://kinetics.nist.gov/kinetics/Detail?id=1988DUR/AMO636:5 , which agree with RMG\'s original values) that are slightly diferent than this paper\'s values (p. 637); I can\'t seem to figure out where the NIST values are coming from (maybe Table 3?); therefore, I have changed rateLibrary to use paper parameters of 10^8.8 (/2) and 4.9 kcal/mol (these values seem to actually be taken from other publications, however), which I am assuming to be high-pressure values; also note that values from other sources are available in the NIST Kinetics Database
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 52
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
Ct_rad
1 *3 Ct 1 {2,T}
2 Ct 0 {1,T}
""",
  kf = Arrhenius(A=(5.0E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (700,1300),
  rank = 3,
  old_id = "331",
  short_comment = "Duran et al. [165] Ab initio.",
  long_comment = 
"""
[165] Duran et al. Ab initio. C.D.W divided rate expression by 2, to get rate of addition per site.
C2H2 + CCH --> HC(tb)CCH=CH. (Rxn. 18?) 

NIST Record: http://kinetics.nist.gov/kinetics/Detail?id=1988DUR/AMO636:4
Verified by Greg Magoon: it looks like value is taken from Rxn 18 of Table 3 (1E10), and is apparently non-pressure dependent (and non-temp dependent); based on the table, it looks like Ref. 42 in this paper may be the ultimate source of the value?
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 53
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.05E+11,A_UNITS,"*/",10.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.46,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,1100),
  rank = 4,
  old_id = "332",
  short_comment = "Baulch et al. [95] literature review.",
  long_comment = 
"""
[95] Baulch et al. literature review. C.D.W divided rate expression by 2, to get rate of addition per site.
C2H2 + .OH --> HOCH=CH

pg.583-584: Discussion on evaluated data

OH+C2H2(+m) --> C2H2OH(+m): \"At temperatures below ~1100K and at atmospheric pressure,

the addition channel becomes important and shows a strong pressure dependence.
The following parameters give a reasonable representation of the high temperature data
for k and are also compatible with Atkinson\'s analysis at low temperature ...\"
RMG stores the recommended high-pressure limit rate coefficient, k_inf.

MRH 31-Aug-2009
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 54
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
O_pri_rad
1 *3 O 1 {2,S}
2 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.60E+07,A_UNITS,"+-",0.0),
                 n=(1.70,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.00,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (250,2500),
  rank = 3,
  old_id = "333",
  short_comment = "Miller et al. [166] Transition state theory.",
  long_comment = 
"""
[166] Miller et al. Transition State Theory. C.D.W divided rate expression by 2, to get rate of addition per site. 
Same reaction as #332, #333 ranked as more accurate in rate library than #332, but they are both from relatively old sources from the early \'90s.  

C2H2 + .OH --> HOCH=CH
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 55
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
O_sec_rad
1 *3 O 1 {2,S}
2 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(5.2E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 4,
  old_id = "334",
  short_comment = "Bozzelli et al. [144] based on CH3 addition to C2H2 (NIST)",
  long_comment = 
"""
[144] Bozzelli et al. Based upon CH3 addition to C2H2 (NIST)
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 56
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.69E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "335",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment. 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 57
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.94E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "336",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 58
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.78E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "337",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 59
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.92E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "338",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 60
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.58E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "339",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 61
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.10E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "340",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 62
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H/NonDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.73E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "341",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 63
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.41E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "342",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 64
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.05E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "343",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 65
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/Cs3
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(8.54E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "344",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 66
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cd
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  kf = Arrhenius(A=(7.15E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "345",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 67
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cd
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  kf = Arrhenius(A=(5.44E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "346",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 68
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H/OneDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.38E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "347",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 69
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H/TwoDe
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
""",
  kf = Arrhenius(A=(1.08E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(20.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "348",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 70
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/Cs2
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.83E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "349",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 71
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/Cs
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(3.78E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(19.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "350",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 72
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Cb
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cb 0 {1,S}
""",
  kf = Arrhenius(A=(6.64E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "351",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 73
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H2/Ct
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Ct 0 {1,S}
""",
  kf = Arrhenius(A=(1.25E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "352",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 74
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/H/OneDeC
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.36E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "353",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 75
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_rad/Cs2
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 {Cd,Ct,Cb,CO} 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(7.00E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "354",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 76
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
Cd_pri_rad
1 *3 Cd 1 {2,D}, {3,S}
2 Cd 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.09E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "355",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 77
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
Cd_rad/NonDe
1 *3 Cd 1 {2,D}, {3,S}
2 Cd 0 {1,D}
3 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(7.24E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "356",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 78
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
Cb_rad
1 *3 Cb 1 {2,B}, {3,B}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
""",
  kf = Arrhenius(A=(1.68E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "357",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 79
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.69E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "358",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 80
rate(
  group1 = 
"""
Cd/H2_Cd/H/Nd
1 *1 Cd     0 {2,D}, {3,S}, {4,S}
2 *2 Cd     0 {1,D}, {5,S}, {6,S}
3    H      0 {1,S}
4    H      0 {1,S}
5    H      0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.10E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "359",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 81
rate(
  group1 = 
"""
Cd/H2_Cd/H/Nd
1 *1 Cd     0 {2,D}, {3,S}, {4,S}
2 *2 Cd     0 {1,D}, {5,S}, {6,S}
3    H      0 {1,S}
4    H      0 {1,S}
5    H      0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.22E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "360",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 82
rate(
  group1 = 
"""
Cd/H2_Cd/Nd2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4       H   0 {1,S}
5    {Cs,O} 0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.99E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "361",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 83
rate(
  group1 = 
"""
Cd/H2_Cd/H/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       H       0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.30E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "362",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 84
rate(
  group1 = 
"""
Ca_Cd/H/Nd
1 *1 Cdd 0 {2,D}, {3,D}
2 *2 Cd 0 {1,D} {4,S} {5,S}
3 C 0 {1,D}
4 H 0 {2,S}
5 {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(7.59E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "363",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 85
rate(
  group1 = 
"""
Cd/H2_Cd/Nd/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       {Cs,O}  0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.12E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "364",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 86
rate(
  group1 = 
"""
Ca_Cd/Nd2
1 *1 Cdd 0 {2,D}, {3,D}
2 *2 Cd 0 {1,D} {4,S} {5,S}
3 C 0 {1,D}
4 {Cs,O} 0 {2,S}
5 {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.86E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(8.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "365",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 87
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.06E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "366",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 88
rate(
  group1 = 
"""
Ct/H_Ct/Nd
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 {Cs,O} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.32E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "367",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 89
rate(
  group1 = 
"""
Cd/H2_Ca
1 *1 Cd  0 {2,D}, {3,S}, {4,S}
2 *2 Cdd  0 {1,D}, {5,D}
3 H  0 {1,S}
4 H  0 {1,S}
5 C 0 {2,D}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.23E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "368",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 90
rate(
  group1 = 
"""
Cd/H2_Cd/H/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       H       0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.18E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "369",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 91
rate(
  group1 = 
"""
Cd/H2_Cd/De2
1 *1        Cd  0 {2,D}, {3,S}, {4,S}
2 *2        Cd  0 {1,D}, {5,S}, {6,S}
3           H   0 {1,S}
4           H   0 {1,S}
5 {Cd,Ct,Cb,CO} 0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.58E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "370",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 92
rate(
  group1 = 
"""
Cd/H2_Cd/H/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       H       0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.37E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "371",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 93
rate(
  group1 = 
"""
Ct/H_Ct/De
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.19E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "372",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 94
rate(
  group1 = 
"""
Cd/H2_Cd/Nd/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       {Cs,O}  0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.32E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "373",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment. 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 95
rate(
  group1 = 
"""
Cd/H/Nd_Cd/H2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4   {Cs,O}  0 {1,S}
5       H   0 {2,S}
6       H   0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.36E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "374",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 96
rate(
  group1 = 
"""
Cd/H/Nd_Cd/H2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4   {Cs,O}  0 {1,S}
5       H   0 {2,S}
6       H   0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.42E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "375",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 97
rate(
  group1 = 
"""
Cd/H/Nd_Cd/H2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4   {Cs,O}  0 {1,S}
5       H   0 {2,S}
6       H   0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.24E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "376",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 98
rate(
  group1 = 
"""
Cd/Nd2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.12E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "377",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 99
rate(
  group1 = 
"""
Cd/Nd2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.51E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "378",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 100
rate(
  group1 = 
"""
Cd/H/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.43E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "379",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 101
rate(
  group1 = 
"""
Cd/Nd/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 {Cs,O} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.92E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "380",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 102
rate(
  group1 = 
"""
Ct/Nd_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 {Cs,O} 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.99E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(12.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "381",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 103
rate(
  group1 = 
"""
Cd/H/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.24E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "382",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 104
rate(
  group1 = 
"""
Cd/H/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.36E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(9.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "383",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 105
rate(
  group1 = 
"""
Ct/De_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.49E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "384",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 106
rate(
  group1 = 
"""
Ct/De_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.71E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "385",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 107
rate(
  group1 = 
"""
Cd/Nd/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 {Cs,O} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
C_methyl
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.07E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "386",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 108
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(7.74E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "387",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 109
rate(
  group1 = 
"""
Cd/H2_Cd/H/Nd
1 *1 Cd     0 {2,D}, {3,S}, {4,S}
2 *2 Cd     0 {1,D}, {5,S}, {6,S}
3    H      0 {1,S}
4    H      0 {1,S}
5    H      0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.01E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "388",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 110
rate(
  group1 = 
"""
Cd/H2_Cd/Nd2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4       H   0 {1,S}
5    {Cs,O} 0 {2,S}
6    {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(4.94E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "389",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 111
rate(
  group1 = 
"""
Cd/H2_Cd/H/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       H       0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(3.71E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "390",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 112
rate(
  group1 = 
"""
Cd/H2_Cd/Nd/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       {Cs,O}  0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.13E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "391",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 113
rate(
  group1 = 
"""
Cd/H2_Cd/H/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       H       0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.97E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "392",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 114
rate(
  group1 = 
"""
Cd/H2_Cd/Nd/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       {Cs,O}  0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.43E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "393",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 115
rate(
  group1 = 
"""
Ct/H_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.16E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "394",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 116
rate(
  group1 = 
"""
Ct/H_Ct/Nd
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.12E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "395",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 117
rate(
  group1 = 
"""
Cd/H2_Ca
1 *1 Cd  0 {2,D}, {3,S}, {4,S}
2 *2 Cdd  0 {1,D}, {5,D}
3 H  0 {1,S}
4 H  0 {1,S}
5 C 0 {2,D}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(6.68E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "396",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 118
rate(
  group1 = 
"""
Cd/H2_Ca
1 *1 Cd  0 {2,D}, {3,S}, {4,S}
2 *2 Cdd  0 {1,D}, {5,D}
3 H  0 {1,S}
4 H  0 {1,S}
5 C 0 {2,D}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.73E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "397",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 119
rate(
  group1 = 
"""
Cd/H2_Cd/H/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       H       0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.82E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "398",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 120
rate(
  group1 = 
"""
Cd/H2_Cd/Nd/De
1 *1    Cd      0 {2,D}, {3,S}, {4,S}
2 *2    Cd      0 {1,D}, {5,S}, {6,S}
3       H       0 {1,S}
4       H       0 {1,S}
5       {Cs,O}  0 {2,S}
6 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.04E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "399",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 121
rate(
  group1 = 
"""
Ct/H_Ct/De
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.86E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "400",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 122
rate(
  group1 = 
"""
Ct/H_Ct/De
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.69E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "401",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 123
rate(
  group1 = 
"""
Ct/H_Ct/De
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.07E+14,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "402",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 124
rate(
  group1 = 
"""
Ca_Cd/H/Nd
1 *1 Cdd 0 {2,D}, {3,D}
2 *2 Cd 0 {1,D} {4,S} {5,S}
3 C 0 {1,D}
4 H 0 {2,S}
5 {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.28E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "403",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 125
rate(
  group1 = 
"""
Ca_Cd/Nd2
1 *1 Cdd 0 {2,D}, {3,D}
2 *2 Cd 0 {1,D} {4,S} {5,S}
3 C 0 {1,D}
4 {Cs,O} 0 {2,S}
5 {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.57E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "404",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 126
rate(
  group1 = 
"""
Cd/H/Nd_Cd/H2
1 *1    Cd  0 {2,D}, {3,S}, {4,S}
2 *2    Cd  0 {1,D}, {5,S}, {6,S}
3       H   0 {1,S}
4   {Cs,O}  0 {1,S}
5       H   0 {2,S}
6       H   0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.18E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "405",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 127
rate(
  group1 = 
"""
Cd/Nd2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.17E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "406",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 128
rate(
  group1 = 
"""
Cd/H/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.85E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "407",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 129
rate(
  group1 = 
"""
Cd/Nd/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 {Cs,O} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(9.30E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "408",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 130
rate(
  group1 = 
"""
Ct/Nd_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 {Cs,O} 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(5.77E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "409",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 131
rate(
  group1 = 
"""
Cd/H/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.67E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "410",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 132
rate(
  group1 = 
"""
Cd/Nd/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 {Cs,O} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(5.62E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "411",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 133
rate(
  group1 = 
"""
Cd/H/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(1.33E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "412",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 134
rate(
  group1 = 
"""
Cd/Nd/De_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3 {Cs,O} 0 {1,S}
4 {Cd,Ct,Cb,CO} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(8.20E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "413",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 135
rate(
  group1 = 
"""
Ct/De_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(2.25E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "414",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 136
rate(
  group1 = 
"""
Ct/De_Ct/H
1 *1 Ct 0 {2,T}, {3,S}
2 *2 Ct 0 {1,T}, {4,S}
3 {Cd,Ct,Cb,CO} 0 {1,S}
4 H 0 {2,S}
""",
  group2 = 
"""
H_rad
1 *3 H 1
""",
  kf = Arrhenius(A=(9.99E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 3,
  old_id = "415",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations,without hindered rotor treatment.
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 137
rate(
  group1 = 
"""
Cd/H2_Cd/H2
1 *1 Cd 0 {2,D}, {3,S}, {4,S}
2 *2 Cd 0 {1,D}, {5,S}, {6,S}
3    H  0 {1,S}
4    H  0 {1,S}
5    H  0 {2,S}
6    H  0 {2,S}
""",
  group2 = 
"""
Cd_pri_rad-Cd/Cd
1 *3 Cd 1 {2,D} {3,S}
2 Cd 0 {1,D} {4,D}
3 H 0 {1,S}
4 Cd 0 {2,D}
""",
  kf = Arrhenius(A=(149.0295,A_UNITS,"+-",0.0),
                 n=(3.0074,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10.1708,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1600),
  rank = 3,
  old_id = "416",
  short_comment = "Sandeep CBS-QB3 calculations",
  long_comment = 
"""
Sandeep CBS-QB3 calculations 
""",
   history = [("2010-06-22","Generated from current RMG library.","rwest@mit.edu")]
)


