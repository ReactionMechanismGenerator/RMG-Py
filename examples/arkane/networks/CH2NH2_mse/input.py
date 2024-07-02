#!/usr/bin/env python3
# encoding: utf-8

title = 'CH2NH2 PES'
description = """
Calculations done using ARC at the following levels of theory:
opt: wb97xd/6-311++g(d,p)
freq: wb97xd/6-311++g(d,p)
sp: ccsd(t)-f12/cc-pvqz-f12
1D rotor scans: wb97xd/6-311++g(d,p)
"""

modelChemistry = "CCSD(T)-F12/cc-pvqz-f12"

useHinderedRotors = True
useBondCorrections = False


species('CH2NH2', 'yaml_files/CH2NH2.yml',
        collisionModel = TransportData(sigma=(3.626,'angstrom'), epsilon=(481.8,'J/mol')),
        energyTransferModel = SingleExponentialDown(alpha0=(133,'cm^-1'), T0=(300,'K'), n=0.85),  # C3H4/N2
)

species('CH3NH', 'yaml_files/CH3NH.yml',
        collisionModel = TransportData(sigma=(3.626,'angstrom'), epsilon=(481.8,'J/mol')),
        energyTransferModel = SingleExponentialDown(alpha0=(133,'cm^-1'), T0=(300,'K'), n=0.85),  # C3H4/N2
)

species('CH2NH', 'yaml_files/CH2NH.yml',
        collisionModel = TransportData(sigma=(3.690,'angstrom'), epsilon=(417.0,'J/mol')),
        energyTransferModel = SingleExponentialDown(alpha0=(133,'cm^-1'), T0=(300,'K'), n=0.85),  # C3H4/N2
)

species('H', 'yaml_files/H.yml',
        collisionModel = TransportData(sigma=(2.050,'angstrom'), epsilon=(145.0,'J/mol')),
        energyTransferModel = SingleExponentialDown(alpha0=(133,'cm^-1'), T0=(300,'K'), n=0.85),  # C3H4/N2
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(133,'cm^-1'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,2.64613e-14,-3.72536e-17,1.7192e-20,-2.44483e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3802.52,'K')), NASAPolynomial(coeffs=[2.5,1.04239e-10,-3.81845e-14,6.18592e-18,-3.73869e-22,-745,4.3663], Tmin=(3802.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState('CH2NH_to_CH2NH2','yaml_files/CH2NH_to_CH2NH2.yml')

transitionState('CH2NH_to_CH3NH','yaml_files/CH2NH_to_CH3NH.yml')

transitionState('CH3NH_to_CH2NH2','yaml_files/CH3NH_to_CH2NH2.yml')

reaction(
    label = 'CH3NH = CH2NH2',
    reactants = ['CH3NH'],
    products = ['CH2NH2'],
    transitionState = 'CH3NH_to_CH2NH2',
    tunneling = 'Eckart',
)

reaction(
    label = 'CH2NH + H = CH2NH2',
    reactants = ['CH2NH', 'H'],
    products = ['CH2NH2'],
    transitionState = 'CH2NH_to_CH2NH2',
    tunneling = 'Eckart',
)

reaction(
    label = 'CH2NH + H = CH3NH',
    reactants = ['CH2NH', 'H'],
    products = ['CH3NH'],
    transitionState = 'CH2NH_to_CH3NH',
    tunneling = 'Eckart',
)

network(
    label = 'CH2NH2',
    isomers = ['CH2NH2', 'CH3NH'],
    reactants = [('CH2NH', 'H')],
    bathGas = {'Ar': 1}
)

pressureDependence(
    label='CH2NH2',
    Tmin=(500.0,'K'), Tmax=(2500.0,'K'), Tcount=25,
    Pmin=(0.01,'bar'), Pmax=(100.0,'bar'), Pcount=15,
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'chemically-significant eigenvalues',  # use the CSE method which is more expensive, less robust, yet more accurate, see: http://reactionmechanismgenerator.github.io/RMG-Py/theory/pdep/methods.html#the-chemically-signficant-eigenvalues-method
    interpolationModel = ('chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
)
