#!/usr/bin/env python
# encoding: utf-8

# Define the level of theory ("model chemistry"):
modelChemistry = LevelOfTheory(method="CCSD(T)-F12", basis="cc-pVTZ-F12", software="molpro")
useHinderedRotors = True
useBondCorrections = True

# FIXME: While switching to pytest, the commented-out examples stopped working. I suspect they have not
# been run in some time and fell vitim to tech. debt.

# Define the species:
species("methoxy", "data/methoxy.py", structure=SMILES("C[O]"))
# species("1,2-butadiene", "data/1,2-butadiene.py", structure=SMILES("C=C=CC"))
species("aziridine", "data/aziridine.py", structure=SMILES("C1NC1"))
species("hydrazino", "data/hydrazino.py", structure=SMILES("N[NH]"))
species("1-propene-12-diol", "data/1-propene-12-diol.py", structure=SMILES("CC(O)=CO"))
species("hydroxyiminomethyl", "data/hydroxyiminomethyl.py", structure=SMILES("[CH]=NO"))
species("2-methyl-2-propanamine", "data/2-methyl-2-propanamine.py", structure=SMILES("CC(C)(C)N"))
# species("nitrosodioxaziridine", "data/nitrosodioxaziridine.py", structure=SMILES("O=NN1OO1"))
species("ethynol", "data/ethynol.py", structure=SMILES("C#CO"))
species("2-iminoethyl", "data/2-iminoethyl.py", structure=SMILES("[CH2]C=N"))

# Request thermodynamic property calculation with a NASA polynomial output:
thermo("methoxy", "NASA")
# thermo("1,2-butadiene", "NASA")
thermo("aziridine", "NASA")
thermo("hydrazino", "NASA")
thermo("1-propene-12-diol", "NASA")
thermo("hydroxyiminomethyl", "NASA")
thermo("2-methyl-2-propanamine", "NASA")
# thermo("nitrosodioxaziridine", "NASA")
thermo("ethynol", "NASA")
thermo("2-iminoethyl", "NASA")
