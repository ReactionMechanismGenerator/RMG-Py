# encoding: utf-8
header = """
R_Addition_COm

COm  + Y_rad -> YC.=O


Reverse name: COM_Elimination_From_Carbonyl

(1) FORM_BOND		{*1,S,*2}
(2) LOSE_RADICAL 	{*1,1}
(3) LOSE_RADICAL 	{*2,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "R_Addition_COm"

# determines permitted units for rate expression:
reaction_order = 2

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f04: radical addition to COm

// jing, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// Catherina Wijaya thesis, pg 155

// f04_radical_addition_to_COm
// rate constants from rate_library_4.txt, Cath, 03/07/28

//No.		COm		Y_rad		Temp.		A			n		a		E0		DA		Dn		Da		DE0		Rank	Comments
//422.		COm		Cb_rad		295-500		8.51E+11	0		0		2.99	*1.5	0		0		0.22	3		Nam et al [104].

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
COm
1 *1 C 2 {2,D}
2 O 0 {1,D}
""",
  group2 = 
"""
Y_rad
1 *2 R 1
""",
  kf = Arrhenius(A=(1E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "416",
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
COm
1 *1 C 2 {2,D}
2 O 0 {1,D}
""",
  group2 = 
"""
H_rad
1 *2 H 1
""",
  kf = Arrhenius(A=(1.18E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(2.72,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (345,449),
  rank = 4,
  old_id = "417",
  short_comment = "Arai et al [102].",
  long_comment = 
"""
[102] Arai, H.; Nagai, S.; Hatada, M.; Radiat. Phys. Chem. 1981, 17, 211.
CO + H --> HCO. Data estimated

pg.215, Table 2: Estimated values of k2 and the reference value of k3 used for the estimation of k2

Raw data is (Temperature [=] degC, k2 [=] cm3/molecule/s):

(72, 3.9x10^-15), (120, 5.6x10^-15), (176, 9.8x10^-15)
Plotting ln(k) vs. 1000/T[=K] and performing a \"Linear\" regression in Microsoft Excel results

in \"y = -1.3657x - 29.258\" with an R^2 value of 0.9762.  The A and Ea values calculated
by MRH are thus: A=1.18x10^11 cm3/mol/s, Ea=2.71 kcal/mol, in agreement w/database.
The authors performed an electron beam irradiation of a CH4 gas stream, containing small

amounts of CO, in a flow system at 1atm.  Authors observe a large decrease in H2 with 
the addition of small amounts of CO.  They assume that this observation must be due to 
H+CO-->HCO.  They propose the following mechanism:
(1) CH4 = H + CH3		k1
(2) H + CO = HCO		k2
(3) H + CH4 = H2 + CH3	k3
The rate of H2 formation:

d[H2]/dt = RM(H2) + k3[H][CH4]
where RM(H2) is the production of H2 through reactions NOT involving H atoms.

Using PSSA on [H]:

d[H]/dt = k1*I*[CH4] - k2[H][CO] - k3[H][CH4] = 0
where I is the dosage rate.

Solving for [H] and substituting into the rate of [H2] formation:

d[H2]/dt = RM(H2) + k1*k3*I*[CH4]^2 / (k2[CO] + k3[CH4])
Subtracting RM(H2) from both sides, taking the inverse of the expression, and rearrangement yields:

{d[H2]/dt - RM(H2)}^-1 = {1 + (k2/k3)*([CO]/[CH4])} / {k1*I*[CH4]}
The authors then introduce this \"G-value\" (proportional to how they detect H2 and CH4???):

{G(H2) - GM(H2)}^-1 = {1 + (k2/k3)*([CO]/[CH4])} / G(H)
The authors present a plot of {G(H2) - GM(H2)}^-1 vs. [CO]/[CH4] to show it is linear.

*** NOTE: The authors assume a value of GM(H2) of 4.63, according to Okazaki et al. ***

From the plot, they extract a (k2/k3) ratio for each temperature tested.  Using the k3 values

reported by Sepehrad et al., they estimate a value of k2.
*** NOTE: Value of k3 used: (72C, 1.52x10^-18 cm3/molecule/s), (120C, 1.51x10^-17 cm3/molecule/s),

(176C, 1.23x10^-16 cm3/molecule/s). ***
MRH 1-Sept-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
COm
1 *1 C 2 {2,D}
2 O 0 {1,D}
""",
  group2 = 
"""
H_rad
1 *2 H 1
""",
  kf = Arrhenius(A=(1.87E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.53,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (305,375),
  rank = 3,
  old_id = "418",
  short_comment = "Gordon et al [103].",
  long_comment = 
"""
[103] Gordon, E.B.; Ivanov, B.I; Perminov, A.P; Balalaev, V.E. Chem. Phys. 1978, 35, 79.
CO + H --> HCO.

Absolute value measured directly with flash photolysis technique. Rate constant is an upper limit.

pg.86, Table 1: H atom reactions with CO and SO2.  Experimentally determined are line shifts

dv and line broadening deltav/2; calculated are rate constants k and complex lifetimes tau_C.
Raw data is (Temperature [=] K, k [=] cm3/molecule/s):

(305, >2.5x10^-14), (375, >4.0x10^-14)
Plotting ln(k) vs. 1000/T[=K] and performing a \"Linear\" regression in Microsoft Excel results

in \"y = -0.768x - 28.802\" with an R^2 value of 1.  The A and Ea values calculated
by MRH are thus: A=1.87x10^11 cm3/mol/s, Ea=1.53 kcal/mol, in agreement w/database.
*** NOTE: MRH interprets table and \"H + CO --> HCO\" discussion to mean that the rate

coefficients reported are LOWER LIMITS.  The discussion appears to suggest that 
the authors suspect oxygen contamination; they further note that the reaction between
H-atom and O2 is 10^4 times faster than the H+CO-->HCO rxn. ***
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
COm
1 *1 C 2 {2,D}
2 O 0 {1,D}
""",
  group2 = 
"""
C_methyl
1 *2 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.06E+11,A_UNITS,"*/",3.16),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(6.88,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,500),
  rank = 4,
  old_id = "419",
  short_comment = "Baulch et al. [94]",
  long_comment = 
"""
[94] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Frank, P.; Hayman, G,; Just, T.; Kerr, J.A.; Murrells, T.; Pilling, M.J.; 
Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1994, 23, 847.

CO + CH3 --> CH3C0. Extensive literature review.

pg 871 Evaluated Kinetic Data for Combustion Modelling Supplement 1, Table 3. Combination reactions.

RMG data matches reference data for k(infinity).

Verified by Karma James

pg.973-974: Discussion on evaluated data

CH3+CO(+m) --> CH3CO(+m): RMG stores k_inf rate coefficient.  The recommended rate

coefficient comes from the preferred (from this reference) rxn rate and the equilibrium
constant (from Bencsura et al.)
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
COm
1 *1 C 2 {2,D}
2 O 0 {1,D}
""",
  group2 = 
"""
C_rad/H2/Cs
1 *2 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.51E+11,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.81,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "420",
  short_comment = "Tsang et al [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J.Phys. Chem. Ref. Data 1986, 15, 1087.
CO + C2H5 --> C2H5CO.

pg 1096, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 17,14.

Verified by Karma James

NOTE: Reported rate coefficients are for k_inf (MRH 11Aug2009)

pg. 1178-1179: Discussion on evaluated data

Recommended data (in the form of k_inf) comes from expression given by Watkins and Thompson

Fall-off corrections and collision efficiencies are also available
(although we do not store them in RMG_database)
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
COm
1 *1 C 2 {2,D}
2 O 0 {1,D}
""",
  group2 = 
"""
Cd_pri_rad
1 *2 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.51E+11,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(4.81,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "421",
  short_comment = "Tsang et al [89] literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F. J.Phys. Chem. Ref. Data 1986, 15, 1087.
CO + C2H3 --> CH2=CHCO.

pg 1099, Chemical Kinetic Database For Combustion Chemistry, 2. Index of Reactions and Summary of Recommended Rate Expressions. No. 19,14.

Verified by Karma James

NOTE: Reported rate coefficients are for k_inf (MRH 11Aug2009)

pg. 1198-1199: Discussion of evaluated data

Recommended data (in the form of k_inf) is assumed to be equal to the rate expression

for CO+C2H5-->H3C-CH2-C=O.  Authors note the rxn is in the fall-off region
under all conditions.
Fall-off corrections and collision efficiencies are also available
(although we do not store them in RMG_database).
MRH 28-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
COm
1 *1 C 2 {2,D}
2 O 0 {1,D}
""",
  group2 = 
"""
Cb_rad
1 *2 Cb 1 {2,B}, {3,B}
2 {Cb,Cbf} 0 {1,B}
3 {Cb,Cbf} 0 {1,B}
""",
  kf = Arrhenius(A=(1.48E+12,A_UNITS,"*/",1.5),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.33,E_UNITS,"+-",0.3)
                 ),
  temperature_range = (295,500),
  rank = 3,
  old_id = "422",
  short_comment = "Nam et al [104].",
  long_comment = 
"""
[104] Nam, G.-J.; Xia, W.; Park, J.; Lin, M. Phys. Chem. A 2000, 104, 1233.	
Phenyl + CO --> Benzoyl. Original deltaA = 2.8E+11

Absolute value measrued directly. Rate constant is high pressure limit. 

Pressure 0.02-0.16 atm. Excitation: flash photolysis, analysis: Vis-UV absorption.

Authors use a Beer-Lambert law type expression:

1/tc = 1/tc_0 + (c*l*epsilon / n*L) * [A](t)
where tc and tc_0 are the decay times of the injected probing photons in the presence

and absence of absorbing species, c is the speed of light, l is the length of the
absorbing medium, epsilon is the extinction coefficient, n is the refractive index
of the medium, L is the length of the cavity, and [A](t) is the concentration of
the absorbing species at time t.
Assuming a simple association rxn, A decays exponentially: [A](t) = [A](0)*exp(-k\'*t).

Combining this with the previous expression yields:
ln(1/tc - 1/tc_0) = B - k\'*t		eq. (*)
However, the authors assume the reverse rxn will be significant (C6H5 + CO <--> C6H5CO).

Thus, they propose the following rate equation:
dx/dt = kf([A](0) - x)[CO] - kr*x
where x is defined as [A](0) - [A](t), [A](t) is the concentration of
the C6H5CO radical at time t, kf is the rate coefficient for C6H5+CO-->C6H5CO,
and kr is the rate coefficient for C6H5CO-->C6H5+CO.
Integrating the above differential equation, assuming constant [CO], yields:

x = (a/b) * (1-exp(-b*t))
where a = kf*[CO]*[A](0) and b = kf*[CO] + kr
Recalling that x = [A](0) - [A](t):

[A](t) = [A](0) - x = [A](0) * {kr + kf*[CO]*exp(-b*t)} / b
Substituting this into the Beer-Lambert law expression:

1/tc - 1/tc_0 = [A](0) * {kr + kf*[CO]*exp(-b*t)} / b		eq. (**)
C6H5 radical was generated from C6H5NO.  The rate coefficient for the C6H5+CO reaction

was measured in the temperature range 295-500K at 12-120 torr, with Ar as the
carrier gas.  The authors note that plots of ln(1/tc - 1/tc_0) vs. t exhibited
linear behavior (for a given Temperature and [CO] concentration).  The slope of
the plot, computed using a \"standard weighted least-squares analysis\", yielded k\',
the pseudo first-order rate coefficient {eq. (*)}.  The authors also note that above 400K,
the plots became nonlinear with time, which the authors attribute to C6H5
re-generation from the reverse rxn C6H5CO --> C6H5 + CO.  This data was analyzed
using eq. (**), to yield b.  The pseudo first-order rate coefficients (either k\' or b)
were plotted against [CO] to yield the second-order rate coefficient for C6H5+CO.
The authors note that the evaluated kf calculated above and below 400K differ greatly.
The authors performed a \"weighted least-squares analysis\" on all data to arrive at
the reported bimolecular rate coefficient:
k1 = 10^11.93+/-0.14 * exp[(-1507+/-109)/T] cm3/mole/s
valid from 295-500K at 40 torr Ar pressure.
The authors also investigated the pressure dependence of the rxn at 347K, from 12-120 torr.

At 347K, the authors do not observe any significant difference.  However, at higher
temperatures, pressure effects become significant.  The authors performed RRKM
calculations to account for falloff effects, and report the adjusted second-order
rate coefficient as:
k1_inf = 10^12.17+/-0.18 * exp[(-1676+/-149)/T] cm3/mole/s
*** NOTE: RMG database was storing reported k1 value.  MRH has changed this so that RMG

now stores the k1_inf value. ***
MRH 1-Sept-2009		
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
COm
1 *1 C 2 {2,D}
2 O 0 {1,D}
""",
  group2 = 
"""
O_rad/NonDe
1 *2 O 1 {2,S}
2 {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(3.41E+07,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(3.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (250,2500),
  rank = 5,
  old_id = "423",
  short_comment = "Wang et al. [105].",
  long_comment = 
"""
[105] Wang, B.; Hou, H.; Gu, Y. Phys. Chem. A 1999, 103, 8021.
RRK(M) extrapolation. CH3O + CO --> CH3OCO, 250K and 2500K

Data stored in RMG appears to be linear fit of the following data, presented on pg.8028

in the right-hand column under the section heading \"3.Implications for Atmospheric
and Combustion Chemistry.\": (250K, 5torr, 1.39x10^-19 cm3/molecule/s) and 
(2500K, 760torr, 3.10x10^-17 cm3/molecule/s).
Plotting ln(k) vs. 1000/T[=K] and performing a \"Linear\" regression in Microsoft Excel results

in \"y = -1.502x - 37.412\" with an R^2 value of 1.  The A and Ea values calculated
by MRH are thus: A=3.40x10^7 cm3/mol/s, Ea=2.98 kcal/mol, in agreement w/database.
MRH 1-Sept-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


