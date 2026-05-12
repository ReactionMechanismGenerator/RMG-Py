#!/usr/bin/env python
# encoding: utf-8

name = "Disproportionation/rules"
shortDesc = ""
longDesc = """

"""
entry(
    index = 1,
    label = "Root",
    kinetics = ArrheniusBM(A=(211.065,'m^3/(mol*s)'), n=1.40533, w0=(576945,'J/mol'), E0=(16272,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.01126358049348453, var=2.2286817499296356, Tref=1000.0, N=137, data_mean=0.0, correlation='Root',), comment="""BM rule fitted to 137 training reactions at node Root
    Total Standard Deviation in ln(k): 3.021123350230915"""),
    rank = 11,
    shortDesc = """BM rule fitted to 137 training reactions at node Root
Total Standard Deviation in ln(k): 3.021123350230915""",
    longDesc = 
"""
BM rule fitted to 137 training reactions at node Root
Total Standard Deviation in ln(k): 3.021123350230915
""",
)

entry(
    index = 2,
    label = "Root_Ext-1R!H-R",
    kinetics = ArrheniusBM(A=(9445.08,'m^3/(mol*s)'), n=0.508694, w0=(543500,'J/mol'), E0=(-990.596,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-1.8021051701158732, var=30.01095873466917, Tref=1000.0, N=45, data_mean=0.0, correlation='Root_Ext-1R!H-R',), comment="""BM rule fitted to 45 training reactions at node Root_Ext-1R!H-R
    Total Standard Deviation in ln(k): 15.510294010456205"""),
    rank = 11,
    shortDesc = """BM rule fitted to 45 training reactions at node Root_Ext-1R!H-R
Total Standard Deviation in ln(k): 15.510294010456205""",
    longDesc = 
"""
BM rule fitted to 45 training reactions at node Root_Ext-1R!H-R
Total Standard Deviation in ln(k): 15.510294010456205
""",
)

entry(
    index = 3,
    label = "Root_Ext-2R!H-R",
    kinetics = ArrheniusBM(A=(219.804,'m^3/(mol*s)'), n=1.42205, w0=(552000,'J/mol'), E0=(50333.1,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.07130535463745602, var=3.640669304482096, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_Ext-2R!H-R',), comment="""BM rule fitted to 6 training reactions at node Root_Ext-2R!H-R
    Total Standard Deviation in ln(k): 4.004301565335125"""),
    rank = 11,
    shortDesc = """BM rule fitted to 6 training reactions at node Root_Ext-2R!H-R
Total Standard Deviation in ln(k): 4.004301565335125""",
    longDesc = 
"""
BM rule fitted to 6 training reactions at node Root_Ext-2R!H-R
Total Standard Deviation in ln(k): 4.004301565335125
""",
)

entry(
    index = 4,
    label = "Root_4R->H",
    kinetics = ArrheniusBM(A=(17085.5,'m^3/(mol*s)'), n=1.00205, w0=(591750,'J/mol'), E0=(67061.8,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.007188943348771603, var=0.2709148414024463, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_4R->H',), comment="""BM rule fitted to 12 training reactions at node Root_4R->H
    Total Standard Deviation in ln(k): 1.0615168637031318"""),
    rank = 11,
    shortDesc = """BM rule fitted to 12 training reactions at node Root_4R->H
Total Standard Deviation in ln(k): 1.0615168637031318""",
    longDesc = 
"""
BM rule fitted to 12 training reactions at node Root_4R->H
Total Standard Deviation in ln(k): 1.0615168637031318
""",
)

entry(
    index = 5,
    label = "Root_N-4R->H",
    kinetics = ArrheniusBM(A=(93.1776,'m^3/(mol*s)'), n=1.4801, w0=(596905,'J/mol'), E0=(23.239,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.006180387938361382, var=2.098032515364293, Tref=1000.0, N=74, data_mean=0.0, correlation='Root_N-4R->H',), comment="""BM rule fitted to 74 training reactions at node Root_N-4R->H
    Total Standard Deviation in ln(k): 2.919304514507989"""),
    rank = 11,
    shortDesc = """BM rule fitted to 74 training reactions at node Root_N-4R->H
Total Standard Deviation in ln(k): 2.919304514507989""",
    longDesc = 
"""
BM rule fitted to 74 training reactions at node Root_N-4R->H
Total Standard Deviation in ln(k): 2.919304514507989
""",
)

entry(
    index = 6,
    label = "Root_Ext-1R!H-R_4R->O",
    kinetics = ArrheniusBM(A=(2.75017e+18,'m^3/(mol*s)'), n=-3.9301, w0=(563000,'J/mol'), E0=(80386.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.09276657605815747, var=3.011300276946607, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O',), comment="""BM rule fitted to 12 training reactions at node Root_Ext-1R!H-R_4R->O
    Total Standard Deviation in ln(k): 3.711918376666456"""),
    rank = 11,
    shortDesc = """BM rule fitted to 12 training reactions at node Root_Ext-1R!H-R_4R->O
Total Standard Deviation in ln(k): 3.711918376666456""",
    longDesc = 
"""
BM rule fitted to 12 training reactions at node Root_Ext-1R!H-R_4R->O
Total Standard Deviation in ln(k): 3.711918376666456
""",
)

entry(
    index = 7,
    label = "Root_Ext-1R!H-R_N-4R->O",
    kinetics = ArrheniusBM(A=(13564.2,'m^3/(mol*s)'), n=0.470009, w0=(536409,'J/mol'), E0=(-16264.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-2.1185860084328203, var=34.56885654617875, Tref=1000.0, N=33, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O',), comment="""BM rule fitted to 33 training reactions at node Root_Ext-1R!H-R_N-4R->O
    Total Standard Deviation in ln(k): 17.10997764410757"""),
    rank = 11,
    shortDesc = """BM rule fitted to 33 training reactions at node Root_Ext-1R!H-R_N-4R->O
Total Standard Deviation in ln(k): 17.10997764410757""",
    longDesc = 
"""
BM rule fitted to 33 training reactions at node Root_Ext-1R!H-R_N-4R->O
Total Standard Deviation in ln(k): 17.10997764410757
""",
)

entry(
    index = 8,
    label = "Root_Ext-2R!H-R_2R!H->C",
    kinetics = ArrheniusBM(A=(4.08618e+20,'m^3/(mol*s)'), n=-4.37582, w0=(551000,'J/mol'), E0=(103613,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.2202160600853793, var=10.39330756070133, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-2R!H-R_2R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-2R!H-R_2R!H->C
    Total Standard Deviation in ln(k): 9.528865376815661"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-2R!H-R_2R!H->C
Total Standard Deviation in ln(k): 9.528865376815661""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-2R!H-R_2R!H->C
Total Standard Deviation in ln(k): 9.528865376815661
""",
)

entry(
    index = 9,
    label = "Root_Ext-2R!H-R_N-2R!H->C",
    kinetics = ArrheniusBM(A=(1273.64,'m^3/(mol*s)'), n=1.17018, w0=(552500,'J/mol'), E0=(39987.4,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.05673151786327586, var=3.6789860426599743, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_Ext-2R!H-R_N-2R!H->C',), comment="""BM rule fitted to 4 training reactions at node Root_Ext-2R!H-R_N-2R!H->C
    Total Standard Deviation in ln(k): 3.9877603244998157"""),
    rank = 11,
    shortDesc = """BM rule fitted to 4 training reactions at node Root_Ext-2R!H-R_N-2R!H->C
Total Standard Deviation in ln(k): 3.9877603244998157""",
    longDesc = 
"""
BM rule fitted to 4 training reactions at node Root_Ext-2R!H-R_N-2R!H->C
Total Standard Deviation in ln(k): 3.9877603244998157
""",
)

entry(
    index = 10,
    label = "Root_4R->H_Sp-2R!H-1R!H",
    kinetics = ArrheniusBM(A=(134417,'m^3/(mol*s)'), n=0.760068, w0=(591944,'J/mol'), E0=(73328,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.017301751919963533, var=0.2727982590527378, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H',), comment="""BM rule fitted to 9 training reactions at node Root_4R->H_Sp-2R!H-1R!H
    Total Standard Deviation in ln(k): 1.090546729162876"""),
    rank = 11,
    shortDesc = """BM rule fitted to 9 training reactions at node Root_4R->H_Sp-2R!H-1R!H
Total Standard Deviation in ln(k): 1.090546729162876""",
    longDesc = 
"""
BM rule fitted to 9 training reactions at node Root_4R->H_Sp-2R!H-1R!H
Total Standard Deviation in ln(k): 1.090546729162876
""",
)

entry(
    index = 11,
    label = "Root_4R->H_N-Sp-2R!H-1R!H",
    kinetics = ArrheniusBM(A=(16039.5,'m^3/(mol*s)'), n=0.960818, w0=(591167,'J/mol'), E0=(59116.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.02977097333687675, var=0.0, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_4R->H_N-Sp-2R!H-1R!H',), comment="""BM rule fitted to 3 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H
    Total Standard Deviation in ln(k): 0.07480144054491646"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H
Total Standard Deviation in ln(k): 0.07480144054491646""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H
Total Standard Deviation in ln(k): 0.07480144054491646
""",
)

entry(
    index = 12,
    label = "Root_N-4R->H_4CNOS-u1",
    kinetics = ArrheniusBM(A=(33.2765,'m^3/(mol*s)'), n=1.57455, w0=(595677,'J/mol'), E0=(2495.51,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.005896437090739888, var=1.8220385479226, Tref=1000.0, N=62, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1',), comment="""BM rule fitted to 62 training reactions at node Root_N-4R->H_4CNOS-u1
    Total Standard Deviation in ln(k): 2.720864875743989"""),
    rank = 11,
    shortDesc = """BM rule fitted to 62 training reactions at node Root_N-4R->H_4CNOS-u1
Total Standard Deviation in ln(k): 2.720864875743989""",
    longDesc = 
"""
BM rule fitted to 62 training reactions at node Root_N-4R->H_4CNOS-u1
Total Standard Deviation in ln(k): 2.720864875743989
""",
)

entry(
    index = 13,
    label = "Root_N-4R->H_N-4CNOS-u1",
    kinetics = ArrheniusBM(A=(9958.36,'m^3/(mol*s)'), n=1.02956, w0=(603250,'J/mol'), E0=(54279,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.03669059532945314, var=0.23702390055722886, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1',), comment="""BM rule fitted to 12 training reactions at node Root_N-4R->H_N-4CNOS-u1
    Total Standard Deviation in ln(k): 1.068194711584249"""),
    rank = 11,
    shortDesc = """BM rule fitted to 12 training reactions at node Root_N-4R->H_N-4CNOS-u1
Total Standard Deviation in ln(k): 1.068194711584249""",
    longDesc = 
"""
BM rule fitted to 12 training reactions at node Root_N-4R->H_N-4CNOS-u1
Total Standard Deviation in ln(k): 1.068194711584249
""",
)

entry(
    index = 14,
    label = "Root_Ext-1R!H-R_4R->O_Ext-4O-R",
    kinetics = ArrheniusBM(A=(0.631851,'m^3/(mol*s)'), n=1.38334, w0=(563000,'J/mol'), E0=(4656.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.29630937279034325, var=1.3223468183437184, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O_Ext-4O-R',), comment="""BM rule fitted to 9 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R
    Total Standard Deviation in ln(k): 3.0498077298705226"""),
    rank = 11,
    shortDesc = """BM rule fitted to 9 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R
Total Standard Deviation in ln(k): 3.0498077298705226""",
    longDesc = 
"""
BM rule fitted to 9 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R
Total Standard Deviation in ln(k): 3.0498077298705226
""",
)

entry(
    index = 15,
    label = "Root_Ext-1R!H-R_4R->O_Sp-5R!H-1R!H",
    kinetics = ArrheniusBM(A=(1.70765e+07,'m^3/(mol*s)'), n=4.66546e-07, w0=(563000,'J/mol'), E0=(56300,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=0.9494596051172368, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O_Sp-5R!H-1R!H',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_4R->O_Sp-5R!H-1R!H
    Total Standard Deviation in ln(k): 1.9534182263699047"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_4R->O_Sp-5R!H-1R!H
Total Standard Deviation in ln(k): 1.9534182263699047""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_4R->O_Sp-5R!H-1R!H
Total Standard Deviation in ln(k): 1.9534182263699047
""",
)

entry(
    index = 16,
    label = "Root_Ext-1R!H-R_4R->O_N-Sp-5R!H-1R!H",
    kinetics = ArrheniusBM(A=(6.03e+06,'m^3/(mol*s)'), n=0, w0=(563000,'J/mol'), E0=(56300,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O_N-Sp-5R!H-1R!H',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_N-Sp-5R!H-1R!H
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_N-Sp-5R!H-1R!H
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_N-Sp-5R!H-1R!H
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 17,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R",
    kinetics = ArrheniusBM(A=(3.79584e+06,'m^3/(mol*s)'), n=-0.196833, w0=(535591,'J/mol'), E0=(53559.1,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0037403777898211976, var=0.7799706060483732, Tref=1000.0, N=11, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R',), comment="""BM rule fitted to 11 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R
    Total Standard Deviation in ln(k): 1.779898653326379"""),
    rank = 11,
    shortDesc = """BM rule fitted to 11 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R
Total Standard Deviation in ln(k): 1.779898653326379""",
    longDesc = 
"""
BM rule fitted to 11 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R
Total Standard Deviation in ln(k): 1.779898653326379
""",
)

entry(
    index = 18,
    label = "Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H",
    kinetics = ArrheniusBM(A=(0.433929,'m^3/(mol*s)'), n=1.96758, w0=(539000,'J/mol'), E0=(70245.9,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.1692521391030103, var=6.079427974794933, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H',), comment="""BM rule fitted to 6 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H
    Total Standard Deviation in ln(k): 5.368230882682922"""),
    rank = 11,
    shortDesc = """BM rule fitted to 6 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H
Total Standard Deviation in ln(k): 5.368230882682922""",
    longDesc = 
"""
BM rule fitted to 6 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H
Total Standard Deviation in ln(k): 5.368230882682922
""",
)

entry(
    index = 19,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H",
    kinetics = ArrheniusBM(A=(5958.46,'m^3/(mol*s)'), n=0.568208, w0=(536000,'J/mol'), E0=(-32910.1,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-3.2867690498643043, var=55.8098478782717, Tref=1000.0, N=16, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H',), comment="""BM rule fitted to 16 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H
    Total Standard Deviation in ln(k): 23.23478535087143"""),
    rank = 11,
    shortDesc = """BM rule fitted to 16 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H
Total Standard Deviation in ln(k): 23.23478535087143""",
    longDesc = 
"""
BM rule fitted to 16 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H
Total Standard Deviation in ln(k): 23.23478535087143
""",
)

entry(
    index = 20,
    label = "Root_Ext-2R!H-R_2R!H->C_4R->C",
    kinetics = ArrheniusBM(A=(50000,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(26972.8,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-2R!H-R_2R!H->C_4R->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_2R!H->C_4R->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_2R!H->C_4R->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_2R!H->C_4R->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 21,
    label = "Root_Ext-2R!H-R_2R!H->C_N-4R->C",
    kinetics = ArrheniusBM(A=(7.23e+06,'m^3/(mol*s)'), n=0, w0=(563000,'J/mol'), E0=(97648.2,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-2R!H-R_2R!H->C_N-4R->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_2R!H->C_N-4R->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_2R!H->C_N-4R->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_2R!H->C_N-4R->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 22,
    label = "Root_Ext-2R!H-R_N-2R!H->C_4R->H",
    kinetics = ArrheniusBM(A=(480,'m^3/(mol*s)'), n=1.5, w0=(557500,'J/mol'), E0=(46201.9,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-2R!H-R_N-2R!H->C_4R->H',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_4R->H
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_4R->H
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_4R->H
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 23,
    label = "Root_Ext-2R!H-R_N-2R!H->C_N-4R->H",
    kinetics = ArrheniusBM(A=(3.55396e+06,'m^3/(mol*s)'), n=0.0868444, w0=(550833,'J/mol'), E0=(77421.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.000652761627062249, var=0.3548234455961088, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_Ext-2R!H-R_N-2R!H->C_N-4R->H',), comment="""BM rule fitted to 3 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H
    Total Standard Deviation in ln(k): 1.1958018205112724"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H
Total Standard Deviation in ln(k): 1.1958018205112724""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H
Total Standard Deviation in ln(k): 1.1958018205112724
""",
)

entry(
    index = 24,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1",
    kinetics = ArrheniusBM(A=(240554,'m^3/(mol*s)'), n=0.681806, w0=(599125,'J/mol'), E0=(74070.9,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.11637774989113545, var=0.5068220816987241, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1',), comment="""BM rule fitted to 8 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1
    Total Standard Deviation in ln(k): 1.7196061325740737"""),
    rank = 11,
    shortDesc = """BM rule fitted to 8 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1
Total Standard Deviation in ln(k): 1.7196061325740737""",
    longDesc = 
"""
BM rule fitted to 8 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1
Total Standard Deviation in ln(k): 1.7196061325740737
""",
)

entry(
    index = 25,
    label = "Root_4R->H_Sp-2R!H-1R!H_N-2R!H-u1",
    kinetics = ArrheniusBM(A=(480,'m^3/(mol*s)'), n=1.5, w0=(534500,'J/mol'), E0=(53450,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_N-2R!H-u1',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_N-2R!H-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_N-2R!H-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_N-2R!H-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 26,
    label = "Root_4R->H_N-Sp-2R!H-1R!H_1R!H->C",
    kinetics = ArrheniusBM(A=(240,'m^3/(mol*s)'), n=1.5, w0=(557500,'J/mol'), E0=(55750,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_N-Sp-2R!H-1R!H_1R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_1R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_1R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_1R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 27,
    label = "Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C",
    kinetics = ArrheniusBM(A=(16039.5,'m^3/(mol*s)'), n=0.960818, w0=(608000,'J/mol'), E0=(60800,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.4478648794535599, var=0.0, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C
    Total Standard Deviation in ln(k): 1.1252886418431154"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C
Total Standard Deviation in ln(k): 1.1252886418431154""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C
Total Standard Deviation in ln(k): 1.1252886418431154
""",
)

entry(
    index = 28,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O",
    kinetics = ArrheniusBM(A=(4.86137,'m^3/(mol*s)'), n=1.88298, w0=(647000,'J/mol'), E0=(-2980.15,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.2591287082868978, var=3.746723415601589, Tref=1000.0, N=22, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O',), comment="""BM rule fitted to 22 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O
    Total Standard Deviation in ln(k): 4.531533543179482"""),
    rank = 11,
    shortDesc = """BM rule fitted to 22 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O
Total Standard Deviation in ln(k): 4.531533543179482""",
    longDesc = 
"""
BM rule fitted to 22 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O
Total Standard Deviation in ln(k): 4.531533543179482
""",
)

entry(
    index = 29,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O",
    kinetics = ArrheniusBM(A=(63.3186,'m^3/(mol*s)'), n=1.4714, w0=(567450,'J/mol'), E0=(-1407.03,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.004295895700507421, var=1.692299905031995, Tref=1000.0, N=40, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O',), comment="""BM rule fitted to 40 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O
    Total Standard Deviation in ln(k): 2.6187220517465106"""),
    rank = 11,
    shortDesc = """BM rule fitted to 40 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O
Total Standard Deviation in ln(k): 2.6187220517465106""",
    longDesc = 
"""
BM rule fitted to 40 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O
Total Standard Deviation in ln(k): 2.6187220517465106
""",
)

entry(
    index = 30,
    label = "Root_N-4R->H_N-4CNOS-u1_1R!H->O",
    kinetics = ArrheniusBM(A=(6.06985e+07,'m^3/(mol*s)'), n=-0.0727699, w0=(665667,'J/mol'), E0=(73198.9,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.0730570062729807, var=4.7084908975045465, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_1R!H->O',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O
    Total Standard Deviation in ln(k): 7.046209272340443"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O
Total Standard Deviation in ln(k): 7.046209272340443""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O
Total Standard Deviation in ln(k): 7.046209272340443
""",
)

entry(
    index = 31,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O",
    kinetics = ArrheniusBM(A=(10895.1,'m^3/(mol*s)'), n=1.01432, w0=(582444,'J/mol'), E0=(55106.1,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.009342386299330369, var=0.25836087721149187, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O',), comment="""BM rule fitted to 9 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O
    Total Standard Deviation in ln(k): 1.042464370909332"""),
    rank = 11,
    shortDesc = """BM rule fitted to 9 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O
Total Standard Deviation in ln(k): 1.042464370909332""",
    longDesc = 
"""
BM rule fitted to 9 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O
Total Standard Deviation in ln(k): 1.042464370909332
""",
)

entry(
    index = 32,
    label = "Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H",
    kinetics = ArrheniusBM(A=(1.02976e+08,'m^3/(mol*s)'), n=-1.08436, w0=(563000,'J/mol'), E0=(46128.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.09151076886860776, var=0.2905474551235391, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H',), comment="""BM rule fitted to 8 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H
    Total Standard Deviation in ln(k): 1.3105279586067948"""),
    rank = 11,
    shortDesc = """BM rule fitted to 8 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H
Total Standard Deviation in ln(k): 1.3105279586067948""",
    longDesc = 
"""
BM rule fitted to 8 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H
Total Standard Deviation in ln(k): 1.3105279586067948
""",
)

entry(
    index = 33,
    label = "Root_Ext-1R!H-R_4R->O_Ext-4O-R_N-Sp-5R!H-1R!H",
    kinetics = ArrheniusBM(A=(602200,'m^3/(mol*s)'), n=0, w0=(563000,'J/mol'), E0=(40463,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O_Ext-4O-R_N-Sp-5R!H-1R!H',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_N-Sp-5R!H-1R!H
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_N-Sp-5R!H-1R!H
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_N-Sp-5R!H-1R!H
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 34,
    label = "Root_Ext-1R!H-R_4R->O_Sp-5R!H-1R!H_Ext-1R!H-R",
    kinetics = ArrheniusBM(A=(1.21e+07,'m^3/(mol*s)'), n=0, w0=(563000,'J/mol'), E0=(56300,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O_Sp-5R!H-1R!H_Ext-1R!H-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_Sp-5R!H-1R!H_Ext-1R!H-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_Sp-5R!H-1R!H_Ext-1R!H-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_Sp-5R!H-1R!H_Ext-1R!H-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 35,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C",
    kinetics = ArrheniusBM(A=(3.54362e+06,'m^3/(mol*s)'), n=-0.187345, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0005435713900316853, var=1.2012760832834504, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C',), comment="""BM rule fitted to 8 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C
    Total Standard Deviation in ln(k): 2.1986103517020235"""),
    rank = 11,
    shortDesc = """BM rule fitted to 8 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C
Total Standard Deviation in ln(k): 2.1986103517020235""",
    longDesc = 
"""
BM rule fitted to 8 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C
Total Standard Deviation in ln(k): 2.1986103517020235
""",
)

entry(
    index = 36,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C",
    kinetics = ArrheniusBM(A=(4.55971e+06,'m^3/(mol*s)'), n=-0.222135, w0=(526500,'J/mol'), E0=(52650,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.01226519497806337, var=0.009633055594024382, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C',), comment="""BM rule fitted to 3 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C
    Total Standard Deviation in ln(k): 0.22757807355262039"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C
Total Standard Deviation in ln(k): 0.22757807355262039""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C
Total Standard Deviation in ln(k): 0.22757807355262039
""",
)

entry(
    index = 37,
    label = "Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R",
    kinetics = ArrheniusBM(A=(1.7265,'m^3/(mol*s)'), n=1.78155, w0=(539000,'J/mol'), E0=(71042.8,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.2864715707319417, var=8.102447089509363, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R',), comment="""BM rule fitted to 5 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R
    Total Standard Deviation in ln(k): 6.426215660910921"""),
    rank = 11,
    shortDesc = """BM rule fitted to 5 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R
Total Standard Deviation in ln(k): 6.426215660910921""",
    longDesc = 
"""
BM rule fitted to 5 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R
Total Standard Deviation in ln(k): 6.426215660910921
""",
)

entry(
    index = 38,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R",
    kinetics = ArrheniusBM(A=(17549.8,'m^3/(mol*s)'), n=0.428311, w0=(534500,'J/mol'), E0=(-14543.2,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-3.582895962257546, var=62.60541265629812, Tref=1000.0, N=13, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R',), comment="""BM rule fitted to 13 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R
    Total Standard Deviation in ln(k): 24.8644332369221"""),
    rank = 11,
    shortDesc = """BM rule fitted to 13 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R
Total Standard Deviation in ln(k): 24.8644332369221""",
    longDesc = 
"""
BM rule fitted to 13 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R
Total Standard Deviation in ln(k): 24.8644332369221
""",
)

entry(
    index = 39,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C",
    kinetics = ArrheniusBM(A=(4.56235e+06,'m^3/(mol*s)'), n=-0.16, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-4.305460834090209e-17, var=0.26130883078297484, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C
    Total Standard Deviation in ln(k): 1.0247880034798968"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C
Total Standard Deviation in ln(k): 1.0247880034798968""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C
Total Standard Deviation in ln(k): 1.0247880034798968
""",
)

entry(
    index = 40,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_N-4CHNS->C",
    kinetics = ArrheniusBM(A=(1.81e+06,'m^3/(mol*s)'), n=0, w0=(549500,'J/mol'), E0=(54950,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_N-4CHNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_N-4CHNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_N-4CHNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_N-4CHNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 41,
    label = "Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_4CNO->O",
    kinetics = ArrheniusBM(A=(2.4,'m^3/(mol*s)'), n=2, w0=(571000,'J/mol'), E0=(57100,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_4CNO->O',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_4CNO->O
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_4CNO->O
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_4CNO->O
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 42,
    label = "Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O",
    kinetics = ArrheniusBM(A=(23777.5,'m^3/(mol*s)'), n=0.682531, w0=(540750,'J/mol'), E0=(64614.1,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.6075970314782776, var=0.6553537765114736, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O
    Total Standard Deviation in ln(k): 3.149537412505626"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O
Total Standard Deviation in ln(k): 3.149537412505626""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O
Total Standard Deviation in ln(k): 3.149537412505626
""",
)

entry(
    index = 43,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O",
    kinetics = ArrheniusBM(A=(1151.49,'m^3/(mol*s)'), n=1.37766, w0=(657250,'J/mol'), E0=(60077.6,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.03633496085430713, var=0.028750657654443026, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O
    Total Standard Deviation in ln(k): 0.43121712987894956"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O
Total Standard Deviation in ln(k): 0.43121712987894956""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O
Total Standard Deviation in ln(k): 0.43121712987894956
""",
)

entry(
    index = 44,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O",
    kinetics = ArrheniusBM(A=(358904,'m^3/(mol*s)'), n=0.629663, w0=(579750,'J/mol'), E0=(76357.6,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.18777365918498876, var=0.5057209361704464, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O',), comment="""BM rule fitted to 6 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O
    Total Standard Deviation in ln(k): 1.897441595632223"""),
    rank = 11,
    shortDesc = """BM rule fitted to 6 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O
Total Standard Deviation in ln(k): 1.897441595632223""",
    longDesc = 
"""
BM rule fitted to 6 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O
Total Standard Deviation in ln(k): 1.897441595632223
""",
)

entry(
    index = 45,
    label = "Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C_2R!H->C",
    kinetics = ArrheniusBM(A=(240,'m^3/(mol*s)'), n=1.5, w0=(545000,'J/mol'), E0=(54500,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 46,
    label = "Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(240,'m^3/(mol*s)'), n=1.5, w0=(671000,'J/mol'), E0=(67100,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C_N-2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C_N-2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_N-Sp-2R!H-1R!H_N-1R!H->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 47,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C",
    kinetics = ArrheniusBM(A=(7.43489e+07,'m^3/(mol*s)'), n=-0.0998224, w0=(662885,'J/mol'), E0=(44766.6,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.7298065405965491, var=4.858257525622271, Tref=1000.0, N=13, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C',), comment="""BM rule fitted to 13 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C
    Total Standard Deviation in ln(k): 6.252412638771265"""),
    rank = 11,
    shortDesc = """BM rule fitted to 13 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C
Total Standard Deviation in ln(k): 6.252412638771265""",
    longDesc = 
"""
BM rule fitted to 13 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C
Total Standard Deviation in ln(k): 6.252412638771265
""",
)

entry(
    index = 48,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C",
    kinetics = ArrheniusBM(A=(0.847891,'m^3/(mol*s)'), n=2.07148, w0=(624056,'J/mol'), E0=(-9108.61,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.3066934194206129, var=3.646168388500196, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C',), comment="""BM rule fitted to 9 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C
    Total Standard Deviation in ln(k): 4.598616635309227"""),
    rank = 11,
    shortDesc = """BM rule fitted to 9 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C
Total Standard Deviation in ln(k): 4.598616635309227""",
    longDesc = 
"""
BM rule fitted to 9 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C
Total Standard Deviation in ln(k): 4.598616635309227
""",
)

entry(
    index = 49,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O",
    kinetics = ArrheniusBM(A=(336.66,'m^3/(mol*s)'), n=1.32779, w0=(590706,'J/mol'), E0=(-8150.35,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00558337269774861, var=1.4985896104748897, Tref=1000.0, N=17, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O',), comment="""BM rule fitted to 17 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O
    Total Standard Deviation in ln(k): 2.4681630030595403"""),
    rank = 11,
    shortDesc = """BM rule fitted to 17 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O
Total Standard Deviation in ln(k): 2.4681630030595403""",
    longDesc = 
"""
BM rule fitted to 17 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O
Total Standard Deviation in ln(k): 2.4681630030595403
""",
)

entry(
    index = 50,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O",
    kinetics = ArrheniusBM(A=(32.606,'m^3/(mol*s)'), n=1.47916, w0=(550261,'J/mol'), E0=(62814.9,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.024499240070432513, var=1.3407175726784197, Tref=1000.0, N=23, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O',), comment="""BM rule fitted to 23 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O
    Total Standard Deviation in ln(k): 2.38282578132325"""),
    rank = 11,
    shortDesc = """BM rule fitted to 23 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O
Total Standard Deviation in ln(k): 2.38282578132325""",
    longDesc = 
"""
BM rule fitted to 23 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O
Total Standard Deviation in ln(k): 2.38282578132325
""",
)

entry(
    index = 51,
    label = "Root_N-4R->H_N-4CNOS-u1_1R!H->O_4CNOS->C",
    kinetics = ArrheniusBM(A=(1.21e+06,'m^3/(mol*s)'), n=0, w0=(655500,'J/mol'), E0=(65550,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_1R!H->O_4CNOS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_4CNOS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_4CNOS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_4CNOS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 52,
    label = "Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C",
    kinetics = ArrheniusBM(A=(1.99016e+08,'m^3/(mol*s)'), n=-0.225563, w0=(670750,'J/mol'), E0=(75059.9,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.2522028424953557, var=6.0581551465221715, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C
    Total Standard Deviation in ln(k): 8.080556867668387"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C
Total Standard Deviation in ln(k): 8.080556867668387""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C
Total Standard Deviation in ln(k): 8.080556867668387
""",
)

entry(
    index = 53,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS",
    kinetics = ArrheniusBM(A=(15938.6,'m^3/(mol*s)'), n=0.994035, w0=(571333,'J/mol'), E0=(58835.8,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.032911533750999235, var=0.2853175032213902, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS',), comment="""BM rule fitted to 6 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS
    Total Standard Deviation in ln(k): 1.153523940801922"""),
    rank = 11,
    shortDesc = """BM rule fitted to 6 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 1.153523940801922""",
    longDesc = 
"""
BM rule fitted to 6 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 1.153523940801922
""",
)

entry(
    index = 54,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS",
    kinetics = ArrheniusBM(A=(11361.3,'m^3/(mol*s)'), n=0.960818, w0=(604667,'J/mol'), E0=(60466.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.02977097491065517, var=1.805559322863034e-35, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS
    Total Standard Deviation in ln(k): 0.0748014444991336"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 0.0748014444991336""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 0.0748014444991336
""",
)

entry(
    index = 55,
    label = "Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R",
    kinetics = ArrheniusBM(A=(10000,'m^3/(mol*s)'), n=-2.70943e-08, w0=(563000,'J/mol'), E0=(19737.6,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=6.435288686813596e-10, var=5.162895344862459e-19, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R',), comment="""BM rule fitted to 4 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R
    Total Standard Deviation in ln(k): 3.0573748226362496e-09"""),
    rank = 11,
    shortDesc = """BM rule fitted to 4 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R
Total Standard Deviation in ln(k): 3.0573748226362496e-09""",
    longDesc = 
"""
BM rule fitted to 4 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R
Total Standard Deviation in ln(k): 3.0573748226362496e-09
""",
)

entry(
    index = 56,
    label = "Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-1R!H-R",
    kinetics = ArrheniusBM(A=(10974.5,'m^3/(mol*s)'), n=2.73044e-07, w0=(563000,'J/mol'), E0=(19809.2,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-3.2290956255676566e-17, var=0.06917824979652494, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-1R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-1R!H-R
    Total Standard Deviation in ln(k): 0.5272805777722495"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-1R!H-R
Total Standard Deviation in ln(k): 0.5272805777722495""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-1R!H-R
Total Standard Deviation in ln(k): 0.5272805777722495
""",
)

entry(
    index = 57,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R",
    kinetics = ArrheniusBM(A=(3.28526e+06,'m^3/(mol*s)'), n=-0.168394, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0006212264085765338, var=1.4218690013356239, Tref=1000.0, N=7, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R',), comment="""BM rule fitted to 7 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R
    Total Standard Deviation in ln(k): 2.3920500512879865"""),
    rank = 11,
    shortDesc = """BM rule fitted to 7 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R
Total Standard Deviation in ln(k): 2.3920500512879865""",
    longDesc = 
"""
BM rule fitted to 7 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R
Total Standard Deviation in ln(k): 2.3920500512879865
""",
)

entry(
    index = 58,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_4HS->H",
    kinetics = ArrheniusBM(A=(904000,'m^3/(mol*s)'), n=0, w0=(549500,'J/mol'), E0=(54950,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_4HS->H',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_4HS->H
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_4HS->H
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_4HS->H
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 59,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_N-4HS->H",
    kinetics = ArrheniusBM(A=(1.02405e+07,'m^3/(mol*s)'), n=-0.333202, w0=(515000,'J/mol'), E0=(51500,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.27677043112298655, var=0.0, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_N-4HS->H',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_N-4HS->H
    Total Standard Deviation in ln(k): 0.6954030932738355"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_N-4HS->H
Total Standard Deviation in ln(k): 0.6954030932738355""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_N-4HS->H
Total Standard Deviation in ln(k): 0.6954030932738355
""",
)

entry(
    index = 60,
    label = "Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-6R!H-R",
    kinetics = ArrheniusBM(A=(84300,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(76690.3,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-6R!H-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-6R!H-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-6R!H-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-6R!H-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 61,
    label = "Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-4CHNS-R",
    kinetics = ArrheniusBM(A=(5.7233e-06,'m^3/(mol*s)'), n=3.63493, w0=(539000,'J/mol'), E0=(10543.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.004220313209984355, var=8.94173911161422, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-4CHNS-R',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-4CHNS-R
    Total Standard Deviation in ln(k): 6.005311154003196"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-4CHNS-R
Total Standard Deviation in ln(k): 6.005311154003196""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-4CHNS-R
Total Standard Deviation in ln(k): 6.005311154003196
""",
)

entry(
    index = 62,
    label = "Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Sp-6R!H-4CHNS",
    kinetics = ArrheniusBM(A=(964000,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(92497.8,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Sp-6R!H-4CHNS',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Sp-6R!H-4CHNS
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Sp-6R!H-4CHNS
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Sp-6R!H-4CHNS
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 63,
    label = "Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_N-Sp-6R!H-4CHNS",
    kinetics = ArrheniusBM(A=(2.41e+06,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_N-Sp-6R!H-4CHNS',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_N-Sp-6R!H-4CHNS
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_N-Sp-6R!H-4CHNS
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_N-Sp-6R!H-4CHNS
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 64,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_6R!H->S",
    kinetics = ArrheniusBM(A=(1.27667e-13,'m^3/(mol*s)'), n=4.8323, w0=(527000,'J/mol'), E0=(25354.4,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=10.238614422136449, var=279.3313223037818, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_6R!H->S',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_6R!H->S
    Total Standard Deviation in ln(k): 59.23071623718318"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_6R!H->S
Total Standard Deviation in ln(k): 59.23071623718318""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_6R!H->S
Total Standard Deviation in ln(k): 59.23071623718318
""",
)

entry(
    index = 65,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S",
    kinetics = ArrheniusBM(A=(365253,'m^3/(mol*s)'), n=0.0895424, w0=(535864,'J/mol'), E0=(4633.75,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-1.8281757095360616, var=15.207890755227291, Tref=1000.0, N=11, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S',), comment="""BM rule fitted to 11 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S
    Total Standard Deviation in ln(k): 12.411330975974792"""),
    rank = 11,
    shortDesc = """BM rule fitted to 11 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S
Total Standard Deviation in ln(k): 12.411330975974792""",
    longDesc = 
"""
BM rule fitted to 11 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S
Total Standard Deviation in ln(k): 12.411330975974792
""",
)

entry(
    index = 66,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C_4C-u1",
    kinetics = ArrheniusBM(A=(1.15e+07,'m^3/(mol*s)'), n=-0.32, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C_4C-u1',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C_4C-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C_4C-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C_4C-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 67,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C_N-4C-u1",
    kinetics = ArrheniusBM(A=(1.81e+06,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C_N-4C-u1',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C_N-4C-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C_N-4C-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_4CHNS->C_N-4C-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 68,
    label = "Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O_4CN->C",
    kinetics = ArrheniusBM(A=(1.6,'m^3/(mol*s)'), n=1.87, w0=(547000,'J/mol'), E0=(47358.9,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O_4CN->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O_4CN->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O_4CN->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O_4CN->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 69,
    label = "Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O_N-4CN->C",
    kinetics = ArrheniusBM(A=(1.8,'m^3/(mol*s)'), n=1.94, w0=(534500,'J/mol'), E0=(53450,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O_N-4CN->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O_N-4CN->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O_N-4CN->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-2R!H-R_N-2R!H->C_N-4R->H_N-4CNO->O_N-4CN->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 70,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O_2R!H->C",
    kinetics = ArrheniusBM(A=(2e+07,'m^3/(mol*s)'), n=0, w0=(666000,'J/mol'), E0=(66600,'J/mol'), Tmin=(295,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 71,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O_N-2R!H->C",
    kinetics = ArrheniusBM(A=(480,'m^3/(mol*s)'), n=1.5, w0=(648500,'J/mol'), E0=(59319.9,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O_N-2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O_N-2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_1R!H->O_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 72,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C",
    kinetics = ArrheniusBM(A=(48750,'m^3/(mol*s)'), n=0.958373, w0=(589333,'J/mol'), E0=(58933.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-1.47580866012507, var=4.807609562897205, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C',), comment="""BM rule fitted to 3 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C
    Total Standard Deviation in ln(k): 8.103696573709435"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C
Total Standard Deviation in ln(k): 8.103696573709435""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C
Total Standard Deviation in ln(k): 8.103696573709435
""",
)

entry(
    index = 73,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C",
    kinetics = ArrheniusBM(A=(12490.1,'m^3/(mol*s)'), n=1.0397, w0=(570167,'J/mol'), E0=(68188.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.023704145862347308, var=0.7277245335318464, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C',), comment="""BM rule fitted to 3 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C
    Total Standard Deviation in ln(k): 1.7697329354950213"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C
Total Standard Deviation in ln(k): 1.7697329354950213""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C
Total Standard Deviation in ln(k): 1.7697329354950213
""",
)

entry(
    index = 74,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R",
    kinetics = ArrheniusBM(A=(2.64485e+07,'m^3/(mol*s)'), n=-0.132719, w0=(662045,'J/mol'), E0=(32116.4,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0510078349133705, var=2.263935109688879, Tref=1000.0, N=11, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R',), comment="""BM rule fitted to 11 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R
    Total Standard Deviation in ln(k): 3.144560699232883"""),
    rank = 11,
    shortDesc = """BM rule fitted to 11 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R
Total Standard Deviation in ln(k): 3.144560699232883""",
    longDesc = 
"""
BM rule fitted to 11 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R
Total Standard Deviation in ln(k): 3.144560699232883
""",
)

entry(
    index = 75,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_4CNOS->C",
    kinetics = ArrheniusBM(A=(8.49e+07,'m^3/(mol*s)'), n=0, w0=(655500,'J/mol'), E0=(65550,'J/mol'), Tmin=(298,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_4CNOS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_4CNOS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_4CNOS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_4CNOS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 76,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_N-4CNOS->C",
    kinetics = ArrheniusBM(A=(2.41e+07,'m^3/(mol*s)'), n=0, w0=(679500,'J/mol'), E0=(67950,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_N-4CNOS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_N-4CNOS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_N-4CNOS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_N-4CNOS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 77,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O",
    kinetics = ArrheniusBM(A=(127.138,'m^3/(mol*s)'), n=1.57348, w0=(653000,'J/mol'), E0=(4452.52,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-3.9155328033681194, var=28.66315200497933, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O
    Total Standard Deviation in ln(k): 20.570968575960745"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O
Total Standard Deviation in ln(k): 20.570968575960745""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O
Total Standard Deviation in ln(k): 20.570968575960745
""",
)

entry(
    index = 78,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O",
    kinetics = ArrheniusBM(A=(0.864362,'m^3/(mol*s)'), n=2.04909, w0=(609583,'J/mol'), E0=(-311.902,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.3524242618291431, var=5.122711251797606, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O',), comment="""BM rule fitted to 6 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O
    Total Standard Deviation in ln(k): 5.422886644907272"""),
    rank = 11,
    shortDesc = """BM rule fitted to 6 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O
Total Standard Deviation in ln(k): 5.422886644907272""",
    longDesc = 
"""
BM rule fitted to 6 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O
Total Standard Deviation in ln(k): 5.422886644907272
""",
)

entry(
    index = 79,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R",
    kinetics = ArrheniusBM(A=(202.633,'m^3/(mol*s)'), n=1.393, w0=(597000,'J/mol'), E0=(-5139.23,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.06911874917151527, var=4.397702883000043, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R',), comment="""BM rule fitted to 8 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R
    Total Standard Deviation in ln(k): 4.377735130161602"""),
    rank = 11,
    shortDesc = """BM rule fitted to 8 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R
Total Standard Deviation in ln(k): 4.377735130161602""",
    longDesc = 
"""
BM rule fitted to 8 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R
Total Standard Deviation in ln(k): 4.377735130161602
""",
)

entry(
    index = 80,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS",
    kinetics = ArrheniusBM(A=(625.943,'m^3/(mol*s)'), n=1.27875, w0=(575333,'J/mol'), E0=(57533.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.03990583248034886, var=0.24629950025387834, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS',), comment="""BM rule fitted to 6 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS
    Total Standard Deviation in ln(k): 1.0951872704960075"""),
    rank = 11,
    shortDesc = """BM rule fitted to 6 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 1.0951872704960075""",
    longDesc = 
"""
BM rule fitted to 6 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 1.0951872704960075
""",
)

entry(
    index = 81,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS",
    kinetics = ArrheniusBM(A=(330.615,'m^3/(mol*s)'), n=1.27907, w0=(604667,'J/mol'), E0=(60466.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.03980613352150107, var=0.0, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS
    Total Standard Deviation in ln(k): 0.10001541085804289"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 0.10001541085804289""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 0.10001541085804289
""",
)

entry(
    index = 82,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R",
    kinetics = ArrheniusBM(A=(373.459,'m^3/(mol*s)'), n=1.26122, w0=(547200,'J/mol'), E0=(46931,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.09318681774672263, var=0.1421306232878116, Tref=1000.0, N=10, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R',), comment="""BM rule fitted to 10 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R
    Total Standard Deviation in ln(k): 0.9899271731946009"""),
    rank = 11,
    shortDesc = """BM rule fitted to 10 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R
Total Standard Deviation in ln(k): 0.9899271731946009""",
    longDesc = 
"""
BM rule fitted to 10 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R
Total Standard Deviation in ln(k): 0.9899271731946009
""",
)

entry(
    index = 83,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C",
    kinetics = ArrheniusBM(A=(30.9211,'m^3/(mol*s)'), n=1.45757, w0=(548688,'J/mol'), E0=(78196.5,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.2556974699678075, var=0.6728411382279629, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C',), comment="""BM rule fitted to 8 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C
    Total Standard Deviation in ln(k): 2.2868778768615523"""),
    rank = 11,
    shortDesc = """BM rule fitted to 8 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C
Total Standard Deviation in ln(k): 2.2868778768615523""",
    longDesc = 
"""
BM rule fitted to 8 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C
Total Standard Deviation in ln(k): 2.2868778768615523
""",
)

entry(
    index = 84,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C",
    kinetics = ArrheniusBM(A=(238.944,'m^3/(mol*s)'), n=1.2433, w0=(558900,'J/mol'), E0=(34984.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0012581813448375246, var=0.5817728263131192, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C',), comment="""BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C
    Total Standard Deviation in ln(k): 1.5322535742673742"""),
    rank = 11,
    shortDesc = """BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C
Total Standard Deviation in ln(k): 1.5322535742673742""",
    longDesc = 
"""
BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C
Total Standard Deviation in ln(k): 1.5322535742673742
""",
)

entry(
    index = 85,
    label = "Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C_2R!H->C",
    kinetics = ArrheniusBM(A=(9.04e+07,'m^3/(mol*s)'), n=0, w0=(679500,'J/mol'), E0=(67950,'J/mol'), Tmin=(298,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 86,
    label = "Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(330,'m^3/(mol*s)'), n=1.5, w0=(662000,'J/mol'), E0=(47976.3,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C_N-2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C_N-2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_1R!H->O_N-4CNOS->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 87,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C",
    kinetics = ArrheniusBM(A=(33706.5,'m^3/(mol*s)'), n=0.959594, w0=(564500,'J/mol'), E0=(56450,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.1948601205724478, var=0.35256865674022686, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C
    Total Standard Deviation in ln(k): 1.6799597049858794"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 1.6799597049858794""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 1.6799597049858794
""",
)

entry(
    index = 88,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C",
    kinetics = ArrheniusBM(A=(9237.74,'m^3/(mol*s)'), n=1.04888, w0=(574750,'J/mol'), E0=(56238.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.010334634855692922, var=0.2652013577163164, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C',), comment="""BM rule fitted to 4 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C
    Total Standard Deviation in ln(k): 1.058358967033491"""),
    rank = 11,
    shortDesc = """BM rule fitted to 4 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 1.058358967033491""",
    longDesc = 
"""
BM rule fitted to 4 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 1.058358967033491
""",
)

entry(
    index = 89,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_1CNS->C",
    kinetics = ArrheniusBM(A=(170,'m^3/(mol*s)'), n=1.5, w0=(571000,'J/mol'), E0=(57100,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 90,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C",
    kinetics = ArrheniusBM(A=(11361.3,'m^3/(mol*s)'), n=0.960818, w0=(621500,'J/mol'), E0=(62150,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.4478648794535599, var=0.0, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C
    Total Standard Deviation in ln(k): 1.1252886418431154"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 1.1252886418431154""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 1.1252886418431154
""",
)

entry(
    index = 91,
    label = "Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R_Ext-1R!H-R",
    kinetics = ArrheniusBM(A=(10000,'m^3/(mol*s)'), n=-2.10203e-08, w0=(563000,'J/mol'), E0=(21006.8,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=0.0, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R_Ext-1R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R_Ext-1R!H-R
    Total Standard Deviation in ln(k): 0.0"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R_Ext-1R!H-R
Total Standard Deviation in ln(k): 0.0""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R_Ext-1R!H-R
Total Standard Deviation in ln(k): 0.0
""",
)

entry(
    index = 92,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_Sp-7R!H#4C",
    kinetics = ArrheniusBM(A=(6.03e+06,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_Sp-7R!H#4C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_Sp-7R!H#4C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_Sp-7R!H#4C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_Sp-7R!H#4C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 93,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C",
    kinetics = ArrheniusBM(A=(2.96901e+06,'m^3/(mol*s)'), n=-0.196459, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0007247616677251713, var=0.7512272601497324, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C',), comment="""BM rule fitted to 6 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C
    Total Standard Deviation in ln(k): 1.7393924065119997"""),
    rank = 11,
    shortDesc = """BM rule fitted to 6 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C
Total Standard Deviation in ln(k): 1.7393924065119997""",
    longDesc = 
"""
BM rule fitted to 6 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C
Total Standard Deviation in ln(k): 1.7393924065119997
""",
)

entry(
    index = 94,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_N-4HS->H_Ext-4S-R_Ext-7R!H-R",
    kinetics = ArrheniusBM(A=(763000,'m^3/(mol*s)'), n=0, w0=(515000,'J/mol'), E0=(51500,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_N-4HS->H_Ext-4S-R_Ext-7R!H-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_N-4HS->H_Ext-4S-R_Ext-7R!H-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_N-4HS->H_Ext-4S-R_Ext-7R!H-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_N-4CHNS->C_N-4HS->H_Ext-4S-R_Ext-7R!H-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 95,
    label = "Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-4CHNS-R_Ext-4CHNS-R",
    kinetics = ArrheniusBM(A=(2.89e+07,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(85952.7,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-4CHNS-R_Ext-4CHNS-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-4CHNS-R_Ext-4CHNS-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-4CHNS-R_Ext-4CHNS-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Sp-5R!H=1R!H_Ext-4CHNS-R_Ext-4CHNS-R_Ext-4CHNS-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 96,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_6R!H->S_Ext-2R!H-R",
    kinetics = ArrheniusBM(A=(644,'m^3/(mol*s)'), n=1.19, w0=(515000,'J/mol'), E0=(39525.8,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_6R!H->S_Ext-2R!H-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_6R!H->S_Ext-2R!H-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_6R!H->S_Ext-2R!H-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_6R!H->S_Ext-2R!H-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 97,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C",
    kinetics = ArrheniusBM(A=(779957,'m^3/(mol*s)'), n=0.0248232, w0=(537950,'J/mol'), E0=(22331.6,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.6732134459096422, var=1.1082236828190586, Tref=1000.0, N=10, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C',), comment="""BM rule fitted to 10 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C
    Total Standard Deviation in ln(k): 3.8019198602869193"""),
    rank = 11,
    shortDesc = """BM rule fitted to 10 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C
Total Standard Deviation in ln(k): 3.8019198602869193""",
    longDesc = 
"""
BM rule fitted to 10 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C
Total Standard Deviation in ln(k): 3.8019198602869193
""",
)

entry(
    index = 98,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_N-4CHNS->C",
    kinetics = ArrheniusBM(A=(7.19e-07,'m^3/(mol*s)'), n=3.13, w0=(515000,'J/mol'), E0=(51500,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_N-4CHNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_N-4CHNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_N-4CHNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_N-4CHNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 99,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_2R!H->C",
    kinetics = ArrheniusBM(A=(3.61e+06,'m^3/(mol*s)'), n=0, w0=(549500,'J/mol'), E0=(54950,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 100,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(48483.4,'m^3/(mol*s)'), n=0.959594, w0=(609250,'J/mol'), E0=(60925,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.6761713695632989, var=0.41113237262952895, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C
    Total Standard Deviation in ln(k): 2.984351249045606"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C
Total Standard Deviation in ln(k): 2.984351249045606""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C
Total Standard Deviation in ln(k): 2.984351249045606
""",
)

entry(
    index = 101,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_2R!H->C",
    kinetics = ArrheniusBM(A=(400,'m^3/(mol*s)'), n=1.5, w0=(564000,'J/mol'), E0=(56400,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 102,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(26.4856,'m^3/(mol*s)'), n=1.82402, w0=(573250,'J/mol'), E0=(52116.6,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.20147811603791788, var=0.27148248712845674, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C
    Total Standard Deviation in ln(k): 1.5507732128102802"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C
Total Standard Deviation in ln(k): 1.5507732128102802""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C
Total Standard Deviation in ln(k): 1.5507732128102802
""",
)

entry(
    index = 103,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS",
    kinetics = ArrheniusBM(A=(7.38112e+07,'m^3/(mol*s)'), n=-3.24554e-09, w0=(655500,'J/mol'), E0=(7107.41,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.06280608228446091, var=6.8952486502763435, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS
    Total Standard Deviation in ln(k): 5.421999069670699"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS
Total Standard Deviation in ln(k): 5.421999069670699""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS
Total Standard Deviation in ln(k): 5.421999069670699
""",
)

entry(
    index = 104,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS",
    kinetics = ArrheniusBM(A=(1.33617e+07,'m^3/(mol*s)'), n=-0.103414, w0=(663500,'J/mol'), E0=(25127,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.04090083298548687, var=1.1654817577770016, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS',), comment="""BM rule fitted to 9 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS
    Total Standard Deviation in ln(k): 2.2670273905941976"""),
    rank = 11,
    shortDesc = """BM rule fitted to 9 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS
Total Standard Deviation in ln(k): 2.2670273905941976""",
    longDesc = 
"""
BM rule fitted to 9 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS
Total Standard Deviation in ln(k): 2.2670273905941976
""",
)

entry(
    index = 105,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R",
    kinetics = ArrheniusBM(A=(8.34831e-40,'m^3/(mol*s)'), n=13.7834, w0=(648500,'J/mol'), E0=(-48242.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-2.1311443344629586, var=3.4442370952847847, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R
    Total Standard Deviation in ln(k): 9.075152857245785"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R
Total Standard Deviation in ln(k): 9.075152857245785""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R
Total Standard Deviation in ln(k): 9.075152857245785
""",
)

entry(
    index = 106,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_4CN->C",
    kinetics = ArrheniusBM(A=(1.6,'m^3/(mol*s)'), n=1.87, w0=(638000,'J/mol'), E0=(77381.5,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_4CN->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_4CN->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_4CN->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_4CN->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 107,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C",
    kinetics = ArrheniusBM(A=(44.6418,'m^3/(mol*s)'), n=1.56985, w0=(603900,'J/mol'), E0=(22091.6,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.1610562744436888, var=4.093886631501787, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C',), comment="""BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C
    Total Standard Deviation in ln(k): 4.46091569890246"""),
    rank = 11,
    shortDesc = """BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C
Total Standard Deviation in ln(k): 4.46091569890246""",
    longDesc = 
"""
BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C
Total Standard Deviation in ln(k): 4.46091569890246
""",
)

entry(
    index = 108,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0",
    kinetics = ArrheniusBM(A=(41.769,'m^3/(mol*s)'), n=1.71947, w0=(595400,'J/mol'), E0=(23393.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.059443513509909215, var=0.4132725572855109, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0',), comment="""BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0
    Total Standard Deviation in ln(k): 1.4381251318640216"""),
    rank = 11,
    shortDesc = """BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0
Total Standard Deviation in ln(k): 1.4381251318640216""",
    longDesc = 
"""
BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0
Total Standard Deviation in ln(k): 1.4381251318640216
""",
)

entry(
    index = 109,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0",
    kinetics = ArrheniusBM(A=(564726,'m^3/(mol*s)'), n=-0.242464, w0=(599667,'J/mol'), E0=(42555.6,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=2.455757428679977, var=11.56486943123057, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0
    Total Standard Deviation in ln(k): 12.987779484083884"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0
Total Standard Deviation in ln(k): 12.987779484083884""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0
Total Standard Deviation in ln(k): 12.987779484083884
""",
)

entry(
    index = 110,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C",
    kinetics = ArrheniusBM(A=(1004.69,'m^3/(mol*s)'), n=1.27744, w0=(576500,'J/mol'), E0=(57650,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.2983662469973947, var=5.131928286299004, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C
    Total Standard Deviation in ln(k): 7.803705422168636"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 7.803705422168636""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 7.803705422168636
""",
)

entry(
    index = 111,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C",
    kinetics = ArrheniusBM(A=(556.026,'m^3/(mol*s)'), n=1.27907, w0=(574750,'J/mol'), E0=(57475,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.03980613384668196, var=0.21353467318894948, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C',), comment="""BM rule fitted to 4 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C
    Total Standard Deviation in ln(k): 1.0263997235150666"""),
    rank = 11,
    shortDesc = """BM rule fitted to 4 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 1.0263997235150666""",
    longDesc = 
"""
BM rule fitted to 4 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 1.0263997235150666
""",
)

entry(
    index = 112,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_1CNS->C",
    kinetics = ArrheniusBM(A=(1.2,'m^3/(mol*s)'), n=2, w0=(571000,'J/mol'), E0=(57100,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 113,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C",
    kinetics = ArrheniusBM(A=(330.615,'m^3/(mol*s)'), n=1.27907, w0=(621500,'J/mol'), E0=(62150,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.5988305691570073, var=0.0, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C
    Total Standard Deviation in ln(k): 1.5045994199924806"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 1.5045994199924806""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 1.5045994199924806
""",
)

entry(
    index = 114,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_Sp-5R!H#4CCCNNNSSS",
    kinetics = ArrheniusBM(A=(3.61e+06,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_Sp-5R!H#4CCCNNNSSS',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_Sp-5R!H#4CCCNNNSSS
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_Sp-5R!H#4CCCNNNSSS
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_Sp-5R!H#4CCCNNNSSS
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 115,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS",
    kinetics = ArrheniusBM(A=(369.15,'m^3/(mol*s)'), n=1.26281, w0=(548111,'J/mol'), E0=(46919.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.19390939481463373, var=0.18238905447925002, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS',), comment="""BM rule fitted to 9 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS
    Total Standard Deviation in ln(k): 1.3433723769367325"""),
    rank = 11,
    shortDesc = """BM rule fitted to 9 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS
Total Standard Deviation in ln(k): 1.3433723769367325""",
    longDesc = 
"""
BM rule fitted to 9 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS
Total Standard Deviation in ln(k): 1.3433723769367325
""",
)

entry(
    index = 116,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C",
    kinetics = ArrheniusBM(A=(265.595,'m^3/(mol*s)'), n=1.19634, w0=(550667,'J/mol'), E0=(55066.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.03752710350769617, var=2.3527370323157824, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C
    Total Standard Deviation in ln(k): 3.1692790336585617"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C
Total Standard Deviation in ln(k): 3.1692790336585617""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C
Total Standard Deviation in ln(k): 3.1692790336585617
""",
)

entry(
    index = 117,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C",
    kinetics = ArrheniusBM(A=(12.3143,'m^3/(mol*s)'), n=1.5703, w0=(547500,'J/mol'), E0=(78127.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.3006228376348046, var=0.597653593422656, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C',), comment="""BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C
    Total Standard Deviation in ln(k): 2.3051555325691324"""),
    rank = 11,
    shortDesc = """BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C
Total Standard Deviation in ln(k): 2.3051555325691324""",
    longDesc = 
"""
BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C
Total Standard Deviation in ln(k): 2.3051555325691324
""",
)

entry(
    index = 118,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS",
    kinetics = ArrheniusBM(A=(260.548,'m^3/(mol*s)'), n=1.2433, w0=(537333,'J/mol'), E0=(31366,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.010636409468894204, var=1.4552981397233726, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS
    Total Standard Deviation in ln(k): 2.445151611974016"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 2.445151611974016""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 2.445151611974016
""",
)

entry(
    index = 119,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS",
    kinetics = ArrheniusBM(A=(209.849,'m^3/(mol*s)'), n=1.2433, w0=(591250,'J/mol'), E0=(59125,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.5787018105298811, var=0.0, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS
    Total Standard Deviation in ln(k): 1.4540246495725655"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 1.4540246495725655""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 1.4540246495725655
""",
)

entry(
    index = 120,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C_2R!H->C",
    kinetics = ArrheniusBM(A=(3.01e+07,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 121,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(500,'m^3/(mol*s)'), n=1.5, w0=(590000,'J/mol'), E0=(59000,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C_N-2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C_N-2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 122,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->C",
    kinetics = ArrheniusBM(A=(330,'m^3/(mol*s)'), n=1.5, w0=(577500,'J/mol'), E0=(57750,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 123,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(4747.78,'m^3/(mol*s)'), n=1.1268, w0=(573833,'J/mol'), E0=(53979.1,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.015461369399676534, var=0.5264999263267074, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C
    Total Standard Deviation in ln(k): 1.4934897420245925"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 1.4934897420245925""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 1.4934897420245925
""",
)

entry(
    index = 124,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C_2R!H->C",
    kinetics = ArrheniusBM(A=(170,'m^3/(mol*s)'), n=1.5, w0=(558500,'J/mol'), E0=(55850,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 125,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(170,'m^3/(mol*s)'), n=1.5, w0=(684500,'J/mol'), E0=(68450,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_N-Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 126,
    label = "Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R_Ext-1R!H-R_Ext-8R!H-R",
    kinetics = ArrheniusBM(A=(10000,'m^3/(mol*s)'), n=0, w0=(563000,'J/mol'), E0=(20706.2,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R_Ext-1R!H-R_Ext-8R!H-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R_Ext-1R!H-R_Ext-8R!H-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R_Ext-1R!H-R_Ext-8R!H-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_4R->O_Ext-4O-R_Sp-5R!H-1R!H_Ext-5R!H-R_Ext-1R!H-R_Ext-8R!H-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 127,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C",
    kinetics = ArrheniusBM(A=(4.90603e+06,'m^3/(mol*s)'), n=-0.235751, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0008697138710928922, var=0.39521319869048593, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C',), comment="""BM rule fitted to 5 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C
    Total Standard Deviation in ln(k): 1.2624816505233833"""),
    rank = 11,
    shortDesc = """BM rule fitted to 5 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C
Total Standard Deviation in ln(k): 1.2624816505233833""",
    longDesc = 
"""
BM rule fitted to 5 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C
Total Standard Deviation in ln(k): 1.2624816505233833
""",
)

entry(
    index = 128,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_N-7R!H->C",
    kinetics = ArrheniusBM(A=(241000,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_N-7R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_N-7R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_N-7R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_N-7R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 129,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C",
    kinetics = ArrheniusBM(A=(6.03e+06,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 130,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C",
    kinetics = ArrheniusBM(A=(772178,'m^3/(mol*s)'), n=0.0249446, w0=(537833,'J/mol'), E0=(23950.2,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.5224495907336478, var=0.6251865196623992, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C',), comment="""BM rule fitted to 9 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C
    Total Standard Deviation in ln(k): 2.8978061230720504"""),
    rank = 11,
    shortDesc = """BM rule fitted to 9 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C
Total Standard Deviation in ln(k): 2.8978061230720504""",
    longDesc = 
"""
BM rule fitted to 9 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C
Total Standard Deviation in ln(k): 2.8978061230720504
""",
)

entry(
    index = 131,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C_2NO->N",
    kinetics = ArrheniusBM(A=(720,'m^3/(mol*s)'), n=1.5, w0=(576500,'J/mol'), E0=(57650,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C_2NO->N',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C_2NO->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C_2NO->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C_2NO->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 132,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C_N-2NO->N",
    kinetics = ArrheniusBM(A=(1.81e+07,'m^3/(mol*s)'), n=0, w0=(642000,'J/mol'), E0=(64200,'J/mol'), Tmin=(300,'K'), Tmax=(1000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C_N-2NO->N',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C_N-2NO->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C_N-2NO->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_1CN->C_N-2R!H->C_N-2NO->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 133,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C_2NO->N",
    kinetics = ArrheniusBM(A=(240,'m^3/(mol*s)'), n=1.5, w0=(534500,'J/mol'), E0=(53450,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C_2NO->N',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C_2NO->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C_2NO->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C_2NO->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 134,
    label = "Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C_N-2NO->N",
    kinetics = ArrheniusBM(A=(480,'m^3/(mol*s)'), n=1.5, w0=(612000,'J/mol'), E0=(63724.7,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C_N-2NO->N',), comment="""BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C_N-2NO->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C_N-2NO->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_4R->H_Sp-2R!H-1R!H_2R!H-u1_N-1R!H->O_N-1CN->C_N-2R!H->C_N-2NO->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 135,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS_5R!H->C",
    kinetics = ArrheniusBM(A=(3.01e+07,'m^3/(mol*s)'), n=0, w0=(655500,'J/mol'), E0=(65550,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS_5R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS_5R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS_5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS_5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 136,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS_N-5R!H->C",
    kinetics = ArrheniusBM(A=(1.81e+08,'m^3/(mol*s)'), n=0, w0=(655500,'J/mol'), E0=(56853,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS_N-5R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS_N-5R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS_N-5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_Sp-5R!H=4CCNNOOSS_N-5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 137,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-4CNOS-R",
    kinetics = ArrheniusBM(A=(2.85561e+07,'m^3/(mol*s)'), n=-0.375, w0=(655500,'J/mol'), E0=(65550,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-4.305460834090209e-17, var=0.06912283083491383, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-4CNOS-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-4CNOS-R
    Total Standard Deviation in ln(k): 0.5270693322083513"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-4CNOS-R
Total Standard Deviation in ln(k): 0.5270693322083513""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-4CNOS-R
Total Standard Deviation in ln(k): 0.5270693322083513
""",
)

entry(
    index = 138,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-5R!H-R",
    kinetics = ArrheniusBM(A=(1.81e+07,'m^3/(mol*s)'), n=0, w0=(655500,'J/mol'), E0=(56999.2,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-5R!H-R',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-5R!H-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-5R!H-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-5R!H-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 139,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C",
    kinetics = ArrheniusBM(A=(5.1898e+06,'m^3/(mol*s)'), n=-3.43792e-07, w0=(655500,'J/mol'), E0=(65550,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-1.4852771291652849e-07, var=1.4637540075814492, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C
    Total Standard Deviation in ln(k): 2.4254431787653163"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C
Total Standard Deviation in ln(k): 2.4254431787653163""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C
Total Standard Deviation in ln(k): 2.4254431787653163
""",
)

entry(
    index = 140,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C",
    kinetics = ArrheniusBM(A=(2.63176e+08,'m^3/(mol*s)'), n=-0.401293, w0=(679500,'J/mol'), E0=(37186.8,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-1.1741100101612896, var=4.186068487337028, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C
    Total Standard Deviation in ln(k): 7.051689842241445"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C
Total Standard Deviation in ln(k): 7.051689842241445""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C
Total Standard Deviation in ln(k): 7.051689842241445
""",
)

entry(
    index = 141,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R_2NOS->N",
    kinetics = ArrheniusBM(A=(0.0294,'m^3/(mol*s)'), n=2.69, w0=(662000,'J/mol'), E0=(33464.5,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R_2NOS->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R_2NOS->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R_2NOS->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R_2NOS->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 142,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R_N-2NOS->N",
    kinetics = ArrheniusBM(A=(0.029,'m^3/(mol*s)'), n=2.69, w0=(635000,'J/mol'), E0=(6270.39,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R_N-2NOS->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R_N-2NOS->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R_N-2NOS->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_4CNOS->O_Ext-4O-R_N-2NOS->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 143,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R",
    kinetics = ArrheniusBM(A=(1.76809e+09,'m^3/(mol*s)'), n=-0.613668, w0=(598500,'J/mol'), E0=(43328.7,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.47302926670052275, var=7.002677954598716, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R
    Total Standard Deviation in ln(k): 6.493560675869729"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R
Total Standard Deviation in ln(k): 6.493560675869729""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R
Total Standard Deviation in ln(k): 6.493560675869729
""",
)

entry(
    index = 144,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_2NOS->N",
    kinetics = ArrheniusBM(A=(1.8,'m^3/(mol*s)'), n=1.94, w0=(625500,'J/mol'), E0=(51859.9,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_2NOS->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_2NOS->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_2NOS->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_2NOS->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 145,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_N-2NOS->N",
    kinetics = ArrheniusBM(A=(0.92,'m^3/(mol*s)'), n=1.94, w0=(598500,'J/mol'), E0=(35268.2,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_N-2NOS->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_N-2NOS->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_N-2NOS->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_N-2NOS->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 146,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS",
    kinetics = ArrheniusBM(A=(55.3682,'m^3/(mol*s)'), n=1.72068, w0=(573833,'J/mol'), E0=(26330,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.16555948575625073, var=0.14398858117357702, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS
    Total Standard Deviation in ln(k): 1.1766919183137672"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 1.1766919183137672""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 1.1766919183137672
""",
)

entry(
    index = 147,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS",
    kinetics = ArrheniusBM(A=(27.368,'m^3/(mol*s)'), n=1.71766, w0=(627750,'J/mol'), E0=(62775,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.8809875394276374, var=0.011502560091510204, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS
    Total Standard Deviation in ln(k): 2.4285443457657907"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 2.4285443457657907""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 2.4285443457657907
""",
)

entry(
    index = 148,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_Sp-2R!H-1CNS",
    kinetics = ArrheniusBM(A=(7.23e+06,'m^3/(mol*s)'), n=0, w0=(563000,'J/mol'), E0=(92345.7,'J/mol'), Tmin=(700,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_Sp-2R!H-1CNS',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_Sp-2R!H-1CNS
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 149,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS",
    kinetics = ArrheniusBM(A=(587946,'m^3/(mol*s)'), n=-0.248363, w0=(618000,'J/mol'), E0=(21389.1,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.394125347182226, var=15.830042173788614, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS
    Total Standard Deviation in ln(k): 11.479064056591858"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 11.479064056591858""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS
Total Standard Deviation in ln(k): 11.479064056591858
""",
)

entry(
    index = 150,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C_2R!H->C",
    kinetics = ArrheniusBM(A=(2.41e+07,'m^3/(mol*s)'), n=0, w0=(563000,'J/mol'), E0=(56300,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 151,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(3.6,'m^3/(mol*s)'), n=2, w0=(590000,'J/mol'), E0=(59000,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C_N-2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C_N-2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 152,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N",
    kinetics = ArrheniusBM(A=(467.561,'m^3/(mol*s)'), n=1.27907, w0=(548000,'J/mol'), E0=(54800,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.5988305691570074, var=0.9609060278364029, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N
    Total Standard Deviation in ln(k): 3.469757305105129"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N
Total Standard Deviation in ln(k): 3.469757305105129""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N
Total Standard Deviation in ln(k): 3.469757305105129
""",
)

entry(
    index = 153,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N",
    kinetics = ArrheniusBM(A=(661.23,'m^3/(mol*s)'), n=1.27907, w0=(601500,'J/mol'), E0=(60150,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.5988305691570073, var=0.0, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N
    Total Standard Deviation in ln(k): 1.5045994199924806"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N
Total Standard Deviation in ln(k): 1.5045994199924806""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N
Total Standard Deviation in ln(k): 1.5045994199924806
""",
)

entry(
    index = 154,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C_2R!H->C",
    kinetics = ArrheniusBM(A=(1.2,'m^3/(mol*s)'), n=2, w0=(558500,'J/mol'), E0=(55850,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 155,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(1.2,'m^3/(mol*s)'), n=2, w0=(684500,'J/mol'), E0=(68450,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_N-Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 156,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_2R!H->S",
    kinetics = ArrheniusBM(A=(979000,'m^3/(mol*s)'), n=0, w0=(537500,'J/mol'), E0=(38210.6,'J/mol'), Tmin=(298,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_2R!H->S',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_2R!H->S
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_2R!H->S
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_2R!H->S
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 157,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S",
    kinetics = ArrheniusBM(A=(363.33,'m^3/(mol*s)'), n=1.26517, w0=(549438,'J/mol'), E0=(41379.5,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.1533575328403071, var=0.10389045665815419, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S',), comment="""BM rule fitted to 8 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S
    Total Standard Deviation in ln(k): 1.031487497319548"""),
    rank = 11,
    shortDesc = """BM rule fitted to 8 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S
Total Standard Deviation in ln(k): 1.031487497319548""",
    longDesc = 
"""
BM rule fitted to 8 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S
Total Standard Deviation in ln(k): 1.031487497319548
""",
)

entry(
    index = 158,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C",
    kinetics = ArrheniusBM(A=(460.84,'m^3/(mol*s)'), n=1.19515, w0=(552500,'J/mol'), E0=(55250,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.1546561267678964, var=0.37853885495848577, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C
    Total Standard Deviation in ln(k): 1.6220067411055652"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C
Total Standard Deviation in ln(k): 1.6220067411055652""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C
Total Standard Deviation in ln(k): 1.6220067411055652
""",
)

entry(
    index = 159,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_N-Sp-2R!H-1C",
    kinetics = ArrheniusBM(A=(0.81,'m^3/(mol*s)'), n=1.87, w0=(547000,'J/mol'), E0=(54700,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_N-Sp-2R!H-1C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_N-Sp-2R!H-1C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_N-Sp-2R!H-1C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_N-Sp-2R!H-1C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 160,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C",
    kinetics = ArrheniusBM(A=(69.6524,'m^3/(mol*s)'), n=1.34293, w0=(544000,'J/mol'), E0=(54400,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.43780050013999694, var=0.5999111817527254, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C
    Total Standard Deviation in ln(k): 2.6527474307051118"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 2.6527474307051118""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C
Total Standard Deviation in ln(k): 2.6527474307051118
""",
)

entry(
    index = 161,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(0.725143,'m^3/(mol*s)'), n=1.93933, w0=(549833,'J/mol'), E0=(76446.4,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.008417366290780803, var=1.809976781745819, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C',), comment="""BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C
    Total Standard Deviation in ln(k): 2.718227067151954"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 2.718227067151954""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C
Total Standard Deviation in ln(k): 2.718227067151954
""",
)

entry(
    index = 162,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N",
    kinetics = ArrheniusBM(A=(207.556,'m^3/(mol*s)'), n=1.2433, w0=(511500,'J/mol'), E0=(51150,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.5787018105298811, var=3.7227133182354435, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N
    Total Standard Deviation in ln(k): 5.322027504076752"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N
Total Standard Deviation in ln(k): 5.322027504076752""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N
Total Standard Deviation in ln(k): 5.322027504076752
""",
)

entry(
    index = 163,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_N-2R!H->N",
    kinetics = ArrheniusBM(A=(1.8,'m^3/(mol*s)'), n=1.94, w0=(589000,'J/mol'), E0=(45404.8,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_N-2R!H->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_N-2R!H->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_N-2R!H->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_N-2R!H->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 164,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS_1CNS->C",
    kinetics = ArrheniusBM(A=(0.92,'m^3/(mol*s)'), n=1.94, w0=(534500,'J/mol'), E0=(53450,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS_1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS_1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 165,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS_N-1CNS->C",
    kinetics = ArrheniusBM(A=(0.92,'m^3/(mol*s)'), n=1.94, w0=(648000,'J/mol'), E0=(64800,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS_N-1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS_N-1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_N-Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 166,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1",
    kinetics = ArrheniusBM(A=(400.418,'m^3/(mol*s)'), n=1.43269, w0=(586750,'J/mol'), E0=(45173.2,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.23284945485790864, var=0.6476388655306268, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1
    Total Standard Deviation in ln(k): 2.1983797414238784"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1
Total Standard Deviation in ln(k): 2.1983797414238784""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1
Total Standard Deviation in ln(k): 2.1983797414238784
""",
)

entry(
    index = 167,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_N-2NO-u1",
    kinetics = ArrheniusBM(A=(330,'m^3/(mol*s)'), n=1.5, w0=(548000,'J/mol'), E0=(54800,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_N-2NO-u1',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_N-2NO-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_N-2NO-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_N-2NO-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 168,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-4C-R",
    kinetics = ArrheniusBM(A=(5.25814e+07,'m^3/(mol*s)'), n=-0.55, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=3.5036392791388526, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-4C-R',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-4C-R
    Total Standard Deviation in ln(k): 3.7524652808647927"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-4C-R
Total Standard Deviation in ln(k): 3.7524652808647927""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-4C-R
Total Standard Deviation in ln(k): 3.7524652808647927
""",
)

entry(
    index = 169,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-7C-R",
    kinetics = ArrheniusBM(A=(783000,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-7C-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-7C-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-7C-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-7C-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 170,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Sp-7C-4C",
    kinetics = ArrheniusBM(A=(843000,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Sp-7C-4C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Sp-7C-4C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Sp-7C-4C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Sp-7C-4C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 171,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_N-Sp-7C-4C",
    kinetics = ArrheniusBM(A=(843000,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_N-Sp-7C-4C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_N-Sp-7C-4C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_N-Sp-7C-4C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_N-Sp-7C-4C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 172,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C",
    kinetics = ArrheniusBM(A=(773968,'m^3/(mol*s)'), n=0.0250683, w0=(537688,'J/mol'), E0=(22302.1,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.6712802518990716, var=0.8880216259988452, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C',), comment="""BM rule fitted to 8 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C
    Total Standard Deviation in ln(k): 3.5757938816348322"""),
    rank = 11,
    shortDesc = """BM rule fitted to 8 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C
Total Standard Deviation in ln(k): 3.5757938816348322""",
    longDesc = 
"""
BM rule fitted to 8 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C
Total Standard Deviation in ln(k): 3.5757938816348322
""",
)

entry(
    index = 173,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_N-6BrCClFINOPSi->C",
    kinetics = ArrheniusBM(A=(482000,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_N-6BrCClFINOPSi->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_N-6BrCClFINOPSi->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_N-6BrCClFINOPSi->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_N-6BrCClFINOPSi->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 174,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-4CNOS-R_Ext-4CNOS-R",
    kinetics = ArrheniusBM(A=(3.47e+08,'m^3/(mol*s)'), n=-0.75, w0=(655500,'J/mol'), E0=(65550,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-4CNOS-R_Ext-4CNOS-R',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-4CNOS-R_Ext-4CNOS-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-4CNOS-R_Ext-4CNOS-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_Ext-4CNOS-R_Ext-4CNOS-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 175,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C",
    kinetics = ArrheniusBM(A=(3.40826e+06,'m^3/(mol*s)'), n=-1.95659e-07, w0=(655500,'J/mol'), E0=(65550,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=0.9609060278364027, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C
    Total Standard Deviation in ln(k): 1.9651578851126479"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C
Total Standard Deviation in ln(k): 1.9651578851126479""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C
Total Standard Deviation in ln(k): 1.9651578851126479
""",
)

entry(
    index = 176,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_N-Sp-5R!H-4C",
    kinetics = ArrheniusBM(A=(1.20333e+07,'m^3/(mol*s)'), n=0, w0=(655500,'J/mol'), E0=(65550,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_N-Sp-5R!H-4C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_N-Sp-5R!H-4C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_N-Sp-5R!H-4C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_N-Sp-5R!H-4C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 177,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0",
    kinetics = ArrheniusBM(A=(2.78282e+08,'m^3/(mol*s)'), n=-0.359057, w0=(679500,'J/mol'), E0=(75512.8,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.12167865343175763, var=0.6437215165799257, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0
    Total Standard Deviation in ln(k): 1.914169472146607"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0
Total Standard Deviation in ln(k): 1.914169472146607""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0
Total Standard Deviation in ln(k): 1.914169472146607
""",
)

entry(
    index = 178,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_N-5R!H-u0",
    kinetics = ArrheniusBM(A=(5.7209e+06,'m^3/(mol*s)'), n=0, w0=(679500,'J/mol'), E0=(20296,'J/mol'), Tmin=(298,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_N-5R!H-u0',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_N-5R!H-u0
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_N-5R!H-u0
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_N-5R!H-u0
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 179,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_5R!H->N",
    kinetics = ArrheniusBM(A=(0.92,'m^3/(mol*s)'), n=1.94, w0=(598500,'J/mol'), E0=(24711.5,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_5R!H->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_5R!H->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_5R!H->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_5R!H->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 180,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N",
    kinetics = ArrheniusBM(A=(39.3893,'m^3/(mol*s)'), n=1.71766, w0=(598500,'J/mol'), E0=(7731.75,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.8611054555760191, var=1.2141850405896415, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N
    Total Standard Deviation in ln(k): 4.372600429820919"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N
Total Standard Deviation in ln(k): 4.372600429820919""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N
Total Standard Deviation in ln(k): 4.372600429820919
""",
)

entry(
    index = 181,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N",
    kinetics = ArrheniusBM(A=(55.3682,'m^3/(mol*s)'), n=1.72068, w0=(548000,'J/mol'), E0=(41967.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-1.0426036023134155, var=0.11276809873671895, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N
    Total Standard Deviation in ln(k): 3.2928163591177713"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N
Total Standard Deviation in ln(k): 3.2928163591177713""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N
Total Standard Deviation in ln(k): 3.2928163591177713
""",
)

entry(
    index = 182,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_N-2R!H->N",
    kinetics = ArrheniusBM(A=(0.029,'m^3/(mol*s)'), n=2.69, w0=(625500,'J/mol'), E0=(26705.8,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_N-2R!H->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_N-2R!H->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_N-2R!H->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_N-2R!H->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 183,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS_1CNS->C",
    kinetics = ArrheniusBM(A=(0.014,'m^3/(mol*s)'), n=2.69, w0=(571000,'J/mol'), E0=(57100,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS_1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS_1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 184,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS_N-1CNS->C",
    kinetics = ArrheniusBM(A=(0.014,'m^3/(mol*s)'), n=2.69, w0=(684500,'J/mol'), E0=(68450,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS_N-1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS_N-1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_N-Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 185,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS_1CNS->C",
    kinetics = ArrheniusBM(A=(2.6e+09,'m^3/(mol*s)'), n=-1.26, w0=(551500,'J/mol'), E0=(34913.7,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS_1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS_1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 186,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS_N-1CNS->C",
    kinetics = ArrheniusBM(A=(1.2e+06,'m^3/(mol*s)'), n=-0.34, w0=(684500,'J/mol'), E0=(61174.3,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS_N-1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS_N-1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_N-5R!H-u0_N-Sp-2R!H-1CNS_N-1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 187,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N_2N-u1",
    kinetics = ArrheniusBM(A=(1.2,'m^3/(mol*s)'), n=2, w0=(548000,'J/mol'), E0=(54800,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N_2N-u1',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N_2N-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N_2N-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N_2N-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 188,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N_N-2N-u1",
    kinetics = ArrheniusBM(A=(2.4,'m^3/(mol*s)'), n=2, w0=(548000,'J/mol'), E0=(54800,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N_N-2N-u1',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N_N-2N-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N_N-2N-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_2R!H->N_N-2N-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 189,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N_2CO->C",
    kinetics = ArrheniusBM(A=(2.4,'m^3/(mol*s)'), n=2, w0=(577500,'J/mol'), E0=(57750,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N_2CO->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N_2CO->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N_2CO->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N_2CO->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 190,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N_N-2CO->C",
    kinetics = ArrheniusBM(A=(2.4,'m^3/(mol*s)'), n=2, w0=(625500,'J/mol'), E0=(55947.8,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N_N-2CO->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N_N-2CO->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N_N-2CO->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->N_N-2CO->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 191,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O",
    kinetics = ArrheniusBM(A=(334.46,'m^3/(mol*s)'), n=1.27744, w0=(593500,'J/mol'), E0=(32251.1,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.2786278929612721, var=0.671686074395935, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O
    Total Standard Deviation in ln(k): 2.3430799122511208"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O
Total Standard Deviation in ln(k): 2.3430799122511208""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O
Total Standard Deviation in ln(k): 2.3430799122511208
""",
)

entry(
    index = 192,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O",
    kinetics = ArrheniusBM(A=(1.47985e+07,'m^3/(mol*s)'), n=-0.311932, w0=(534750,'J/mol'), E0=(37491,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.14761419944759951, var=0.2702865095049491, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O',), comment="""BM rule fitted to 6 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O
    Total Standard Deviation in ln(k): 1.4131333979764096"""),
    rank = 11,
    shortDesc = """BM rule fitted to 6 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O
Total Standard Deviation in ln(k): 1.4131333979764096""",
    longDesc = 
"""
BM rule fitted to 6 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O
Total Standard Deviation in ln(k): 1.4131333979764096
""",
)

entry(
    index = 193,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C_2R!H->C",
    kinetics = ArrheniusBM(A=(2.19e+08,'m^3/(mol*s)'), n=-0.68, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C_2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C_2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C_2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 194,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C_N-2R!H->C",
    kinetics = ArrheniusBM(A=(2.4,'m^3/(mol*s)'), n=1.87, w0=(566000,'J/mol'), E0=(56600,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C_N-2R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C_N-2R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_1CNS->C_Sp-2R!H-1C_N-2R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 195,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C_Sp-2C-1N",
    kinetics = ArrheniusBM(A=(1.6,'m^3/(mol*s)'), n=1.87, w0=(553500,'J/mol'), E0=(55350,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C_Sp-2C-1N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C_Sp-2C-1N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C_Sp-2C-1N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C_Sp-2C-1N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 196,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C_N-Sp-2C-1N",
    kinetics = ArrheniusBM(A=(0.82,'m^3/(mol*s)'), n=1.87, w0=(534500,'J/mol'), E0=(53450,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C_N-Sp-2C-1N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C_N-Sp-2C-1N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C_N-Sp-2C-1N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_2R!H->C_N-Sp-2C-1N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 197,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1",
    kinetics = ArrheniusBM(A=(1.9471e-05,'m^3/(mol*s)'), n=3.27783, w0=(562750,'J/mol'), E0=(48062.1,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.8568478925825855, var=0.083327900484, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1
    Total Standard Deviation in ln(k): 2.73158245570468"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1
Total Standard Deviation in ln(k): 2.73158245570468""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1
Total Standard Deviation in ln(k): 2.73158245570468
""",
)

entry(
    index = 198,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_N-2NO-u1",
    kinetics = ArrheniusBM(A=(1.6,'m^3/(mol*s)'), n=1.87, w0=(524000,'J/mol'), E0=(52400,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_N-2NO-u1',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_N-2NO-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_N-2NO-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_N-2NO-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 199,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N_2N-u1",
    kinetics = ArrheniusBM(A=(0.46,'m^3/(mol*s)'), n=1.94, w0=(511500,'J/mol'), E0=(51150,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N_2N-u1',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N_2N-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N_2N-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N_2N-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 200,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N_N-2N-u1",
    kinetics = ArrheniusBM(A=(1.8,'m^3/(mol*s)'), n=1.94, w0=(511500,'J/mol'), E0=(51150,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N_N-2N-u1',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N_N-2N-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N_N-2N-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_N-4CNS->C_Sp-2R!H-1CNS_2R!H->N_N-2N-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 201,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1_2NO->N",
    kinetics = ArrheniusBM(A=(170,'m^3/(mol*s)'), n=1.5, w0=(548000,'J/mol'), E0=(54800,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1_2NO->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1_2NO->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1_2NO->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1_2NO->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 202,
    label = "Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1_N-2NO->N",
    kinetics = ArrheniusBM(A=(330,'m^3/(mol*s)'), n=1.5, w0=(625500,'J/mol'), E0=(52128.2,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1_N-2NO->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1_N-2NO->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1_N-2NO->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_N-4CNOS-u1_N-1R!H->O_Sp-2R!H-1CNS_N-1CNS->C_N-2R!H->C_2NO-u1_N-2NO->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 203,
    label = "Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-4C-R_Ext-4C-R",
    kinetics = ArrheniusBM(A=(1.08e+08,'m^3/(mol*s)'), n=-0.75, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-4C-R_Ext-4C-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-4C-R_Ext-4C-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-4C-R_Ext-4C-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_Ext-1R!H-R_4CHNS->C_Ext-4C-R_N-Sp-7R!H#4C_7R!H->C_Ext-4C-R_Ext-4C-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 204,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_1R!H-inRing",
    kinetics = ArrheniusBM(A=(0.000675,'m^3/(mol*s)'), n=2.7, w0=(483500,'J/mol'), E0=(35785.7,'J/mol'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_1R!H-inRing',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_1R!H-inRing
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_1R!H-inRing
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_1R!H-inRing
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 205,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing",
    kinetics = ArrheniusBM(A=(6.31109e+06,'m^3/(mol*s)'), n=-0.199393, w0=(545429,'J/mol'), E0=(32200.5,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.1700346401571909, var=0.5832966950308113, Tref=1000.0, N=7, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing',), comment="""BM rule fitted to 7 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing
    Total Standard Deviation in ln(k): 1.9583163355851787"""),
    rank = 11,
    shortDesc = """BM rule fitted to 7 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing
Total Standard Deviation in ln(k): 1.9583163355851787""",
    longDesc = 
"""
BM rule fitted to 7 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing
Total Standard Deviation in ln(k): 1.9583163355851787
""",
)

entry(
    index = 206,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C_5R!H->C",
    kinetics = ArrheniusBM(A=(2.41e+06,'m^3/(mol*s)'), n=0, w0=(655500,'J/mol'), E0=(65550,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C_5R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C_5R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C_5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C_5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 207,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C_N-5R!H->C",
    kinetics = ArrheniusBM(A=(4.82e+06,'m^3/(mol*s)'), n=0, w0=(655500,'J/mol'), E0=(65550,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C_N-5R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C_N-5R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C_N-5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_4CNOS->C_Sp-5R!H-4C_N-5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 208,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0_5R!H->C",
    kinetics = ArrheniusBM(A=(2.41e+07,'m^3/(mol*s)'), n=0, w0=(679500,'J/mol'), E0=(67950,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0_5R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0_5R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0_5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0_5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 209,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0_N-5R!H->C",
    kinetics = ArrheniusBM(A=(1.21e+07,'m^3/(mol*s)'), n=0, w0=(679500,'J/mol'), E0=(55695.5,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0_N-5R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0_N-5R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0_N-5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_2R!H->C_Ext-4CNOS-R_N-Sp-5R!H=4CCNNOOSS_N-4CNOS->C_5R!H-u0_N-5R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 210,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N_5BrCClFIOPSSi->C",
    kinetics = ArrheniusBM(A=(0.014,'m^3/(mol*s)'), n=2.69, w0=(598500,'J/mol'), E0=(17330.2,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N_5BrCClFIOPSSi->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N_5BrCClFIOPSSi->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N_5BrCClFIOPSSi->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N_5BrCClFIOPSSi->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 211,
    label = "Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N_N-5BrCClFIOPSSi->C",
    kinetics = ArrheniusBM(A=(0.029,'m^3/(mol*s)'), n=2.69, w0=(598500,'J/mol'), E0=(13374.2,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N_N-5BrCClFIOPSSi->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N_N-5BrCClFIOPSSi->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N_N-5BrCClFIOPSSi->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_1R!H->O_N-2R!H->C_N-4CNOS->O_N-4CN->C_Ext-4N-R_N-5R!H->N_N-5BrCClFIOPSSi->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 212,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N_2N-u1",
    kinetics = ArrheniusBM(A=(0.029,'m^3/(mol*s)'), n=2.69, w0=(548000,'J/mol'), E0=(39565.9,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N_2N-u1',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N_2N-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N_2N-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N_2N-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 213,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N_N-2N-u1",
    kinetics = ArrheniusBM(A=(0.029,'m^3/(mol*s)'), n=2.69, w0=(548000,'J/mol'), E0=(54800,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N_N-2N-u1',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N_N-2N-u1
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N_N-2N-u1
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_4CNOS->O_Ext-4O-R_5R!H-u0_Sp-2R!H-1CNS_2R!H->N_N-2N-u1
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 214,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O_1CNS->C",
    kinetics = ArrheniusBM(A=(2.89e+06,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O_1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O_1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O_1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 215,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O_N-1CNS->C",
    kinetics = ArrheniusBM(A=(1.2,'m^3/(mol*s)'), n=2, w0=(648000,'J/mol'), E0=(54170.7,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O_N-1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O_N-1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O_N-1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_5R!H->O_N-1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 216,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS",
    kinetics = ArrheniusBM(A=(1.05319e+07,'m^3/(mol*s)'), n=-0.250519, w0=(533900,'J/mol'), E0=(39536.6,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.1636050951277374, var=0.25119831907169365, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS',), comment="""BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS
    Total Standard Deviation in ln(k): 1.4158350573264697"""),
    rank = 11,
    shortDesc = """BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS
Total Standard Deviation in ln(k): 1.4158350573264697""",
    longDesc = 
"""
BM rule fitted to 5 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS
Total Standard Deviation in ln(k): 1.4158350573264697
""",
)

entry(
    index = 217,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_N-Sp-5CS-4CCNSS",
    kinetics = ArrheniusBM(A=(1.52e+08,'m^3/(mol*s)'), n=-0.7, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_N-Sp-5CS-4CCNSS',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_N-Sp-5CS-4CCNSS
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_N-Sp-5CS-4CCNSS
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_N-Sp-5CS-4CCNSS
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 218,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1_2NO->N",
    kinetics = ArrheniusBM(A=(0.82,'m^3/(mol*s)'), n=1.87, w0=(524000,'J/mol'), E0=(52400,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1_2NO->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1_2NO->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1_2NO->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1_2NO->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 219,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1_N-2NO->N",
    kinetics = ArrheniusBM(A=(1.6,'m^3/(mol*s)'), n=1.87, w0=(601500,'J/mol'), E0=(75336.7,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1_N-2NO->N',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1_N-2NO->N
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1_N-2NO->N
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_4CNS->C_N-1CNS->C_N-2R!H->C_2NO-u1_N-2NO->N
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 220,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R",
    kinetics = ArrheniusBM(A=(2.36058e+07,'m^3/(mol*s)'), n=-0.372184, w0=(550250,'J/mol'), E0=(35580.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.25457934710988045, var=1.6100541888330673, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R',), comment="""BM rule fitted to 4 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R
    Total Standard Deviation in ln(k): 3.1834130560690546"""),
    rank = 11,
    shortDesc = """BM rule fitted to 4 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R
Total Standard Deviation in ln(k): 3.1834130560690546""",
    longDesc = 
"""
BM rule fitted to 4 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R
Total Standard Deviation in ln(k): 3.1834130560690546
""",
)

entry(
    index = 221,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Sp-6C-4C",
    kinetics = ArrheniusBM(A=(1.97085e+06,'m^3/(mol*s)'), n=-0.0393785, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0327092327690802, var=0.00213978781668374, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Sp-6C-4C',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Sp-6C-4C
    Total Standard Deviation in ln(k): 0.1749187175813773"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Sp-6C-4C
Total Standard Deviation in ln(k): 0.1749187175813773""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Sp-6C-4C
Total Standard Deviation in ln(k): 0.1749187175813773
""",
)

entry(
    index = 222,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_N-Sp-6C-4C",
    kinetics = ArrheniusBM(A=(1.21e+06,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_N-Sp-6C-4C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_N-Sp-6C-4C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_N-Sp-6C-4C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_N-Sp-6C-4C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 223,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R",
    kinetics = ArrheniusBM(A=(7.76827e+08,'m^3/(mol*s)'), n=-0.9, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-1.0763652085225523e-17, var=0.04891149884417046, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R
    Total Standard Deviation in ln(k): 0.4433660913390881"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R
Total Standard Deviation in ln(k): 0.4433660913390881""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R
Total Standard Deviation in ln(k): 0.4433660913390881
""",
)

entry(
    index = 224,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_1CNS->C",
    kinetics = ArrheniusBM(A=(3.11937e+07,'m^3/(mol*s)'), n=-0.389378, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.032709232769080235, var=0.0016076635746038646, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_1CNS->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_1CNS->C
    Total Standard Deviation in ln(k): 0.16256521857885725"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_1CNS->C
Total Standard Deviation in ln(k): 0.16256521857885725""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_1CNS->C
Total Standard Deviation in ln(k): 0.16256521857885725
""",
)

entry(
    index = 225,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_N-1CNS->C",
    kinetics = ArrheniusBM(A=(644,'m^3/(mol*s)'), n=1.19, w0=(513500,'J/mol'), E0=(42507.8,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_N-1CNS->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_N-1CNS->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_N-1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_N-1CNS->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 226,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C",
    kinetics = ArrheniusBM(A=(1.05265e+08,'m^3/(mol*s)'), n=-0.55, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=3.513977146582821, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C
    Total Standard Deviation in ln(k): 3.75799723098171"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C
Total Standard Deviation in ln(k): 3.75799723098171""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C
Total Standard Deviation in ln(k): 3.75799723098171
""",
)

entry(
    index = 227,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_N-2R!H->C",
    kinetics = ArrheniusBM(A=(1.9053e+06,'m^3/(mol*s)'), n=-0.0575531, w0=(561500,'J/mol'), E0=(16527.3,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.2144940523654643, var=1.5169644222011907, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_N-2R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_N-2R!H->C
    Total Standard Deviation in ln(k): 3.008063935012343"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_N-2R!H->C
Total Standard Deviation in ln(k): 3.008063935012343""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_N-2R!H->C
Total Standard Deviation in ln(k): 3.008063935012343
""",
)

entry(
    index = 228,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Sp-6C-4C_Ext-6C-R",
    kinetics = ArrheniusBM(A=(1.45e+06,'m^3/(mol*s)'), n=0, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Sp-6C-4C_Ext-6C-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Sp-6C-4C_Ext-6C-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Sp-6C-4C_Ext-6C-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Sp-6C-4C_Ext-6C-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 229,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R_Ext-4CNS-R",
    kinetics = ArrheniusBM(A=(2.86e+09,'m^3/(mol*s)'), n=-1.1, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R_Ext-4CNS-R',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R_Ext-4CNS-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R_Ext-4CNS-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R_Ext-4CNS-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 230,
    label = "Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_1CNS->C_Ext-5CS-R",
    kinetics = ArrheniusBM(A=(2.29e+07,'m^3/(mol*s)'), n=-0.35, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_1CNS->C_Ext-5CS-R',), comment="""BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_1CNS->C_Ext-5CS-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_1CNS->C_Ext-5CS-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_1CNS->C_Ext-5CS-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 231,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C_Ext-4C-R",
    kinetics = ArrheniusBM(A=(2.16e+08,'m^3/(mol*s)'), n=-0.75, w0=(539000,'J/mol'), E0=(53900,'J/mol'), Tmin=(300,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C_Ext-4C-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C_Ext-4C-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C_Ext-4C-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C_Ext-4C-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 232,
    label = "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_N-2R!H->C_Ext-7R!H-R_Ext-6C-R",
    kinetics = ArrheniusBM(A=(1.94e+06,'m^3/(mol*s)'), n=0, w0=(561500,'J/mol'), E0=(36677.6,'J/mol'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_N-2R!H->C_Ext-7R!H-R_Ext-6C-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_N-2R!H->C_Ext-7R!H-R_Ext-6C-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_N-2R!H->C_Ext-7R!H-R_Ext-6C-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_N-2R!H->C_Ext-7R!H-R_Ext-6C-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

