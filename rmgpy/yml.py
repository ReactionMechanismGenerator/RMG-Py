#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################
"""
This file defines functions for outputting the RMG generated mechanism to a .rms yaml file, which can be read by the
Reaction Mechanism Simulator (RMS)
"""

import os
import yaml

from rmgpy.chemkin import loadChemkinFile
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.thermo.nasa import NASAPolynomial, NASA
from rmgpy.thermo.wilhoit import Wilhoit
from rmgpy.kinetics.arrhenius import Arrhenius, PDepArrhenius, MultiArrhenius, MultiPDepArrhenius
from rmgpy.kinetics.falloff import Troe, ThirdBody, Lindemann
from rmgpy.kinetics.chebyshev import Chebyshev
from rmgpy.data.solvation import SolventData
from rmgpy.kinetics.surface import StickingCoefficient
from rmgpy.util import makeOutputSubdirectory


def convertChemkin2yml(chemkinPath, spcDictPath=None, output="chem.rms"):
    if spcDictPath:
        spcs,rxns = loadChemkinFile(chemkinPath, dictionaryPath=spcDictPath)
    else:
        spcs,rxns = loadChemkinFile(chemkinPath)
    writeyml(spcs, rxns, path=output)


def writeyml(spcs, rxns, solvent=None, solventData=None, path="chem.yml"):
    D = getMechDict(spcs, rxns, solvent=solvent, solventData=solventData)
    with open(path, 'w') as f:
        yaml.dump(D, stream=f)


def getMechDict(spcs, rxns, solvent='solvent', solventData=None):
    D = dict()
    D["Units"] = dict()
    D["Phases"] = [dict()]
    D["Phases"][0]["name"] = "phase"
    D["Phases"][0]["Species"] = [obj2dict(x, spcs) for x in spcs]
    D["Reactions"] = [obj2dict(x, spcs) for x in rxns]
    if solventData:
        D["Solvents"] = [obj2dict(solventData, spcs, label=solvent)]
    return D


def getRadicals(spc):
    if spc.molecule[0].toSMILES() == "[O][O]":  # treat oxygen as stable to improve radical analysis
        return 0
    else:
        return spc.molecule[0].multiplicity-1


def obj2dict(obj, spcs, label="solvent"):
    D = dict()
    if isinstance(obj, Species):
        D["name"] = obj.label
        D["type"] = "Species"
        D["smiles"] = obj.molecule[0].toSMILES()
        D["thermo"] = obj2dict(obj.thermo, spcs)
        D["radicalelectrons"] = getRadicals(obj)
    elif isinstance(obj, NASA):
        D["polys"] = [obj2dict(k, spcs) for k in obj.polynomials]
        D["type"] = "NASA"
    elif isinstance(obj, NASAPolynomial):
        D["type"] = "NASApolynomial"
        D["coefs"] = obj.coeffs.tolist()
        D["Tmax"] = obj.Tmax.value_si
        D["Tmin"] = obj.Tmin.value_si
    elif isinstance(obj, Reaction):
        D["reactants"] = [x.label for x in obj.reactants]
        D["products"] = [x.label for x in obj.products]
        D["kinetics"] = obj2dict(obj.kinetics, spcs)
        D["type"] = "ElementaryReaction"
        D["radicalchange"] = sum([getRadicals(x) for x in obj.products]) - sum([getRadicals(x) for x in obj.reactants])
    elif isinstance(obj, Arrhenius):
        D["type"] = "Arrhenius"
        D["A"] = obj.A.value_si
        D["Ea"] = obj.Ea.value_si
        D["n"] = obj.n.value_si
    elif isinstance(obj, StickingCoefficient):
        D["type"] = "StickingCoefficient"
        D["A"] = obj.A.value_si
        D["Ea"] = obj.Ea.value_si
        D["n"] = obj.n.value_si
        D["T0"] = obj.T0.value_si
    elif isinstance(obj, PDepArrhenius):
        D["type"] = "PdepArrhenius"
        D["Ps"] = obj.pressures.value_si.tolist()
        D["arrs"] = [obj2dict(x, spcs) for x in obj.arrhenius]
    elif isinstance(obj, MultiArrhenius):
        D["type"] = "MultiArrhenius"
        D["arrs"] = [obj2dict(x, spcs) for x in obj.arrhenius]
    elif isinstance(obj, MultiPDepArrhenius):
        D["type"] = "MultiPdepArrhenius"
        D["parrs"] = [obj2dict(x, spcs) for x in obj.arrhenius]
    elif isinstance(obj, ThirdBody):
        D["type"] = "ThirdBody"
        D["arr"] = obj2dict(obj.arrheniusLow, spcs)
        D["efficiencies"] = {spcs[i].label: float(val) for i, val in enumerate(obj.getEffectiveColliderEfficiencies(spcs)) if val != 1}
    elif isinstance(obj, Lindemann):
        D["type"] = "Lindemann"
        D["arrhigh"] = obj2dict(obj.arrheniusHigh, spcs)
        D["arrlow"] = obj2dict(obj.arrheniusLow, spcs)
        D["efficiencies"] = {spcs[i].label: float(val) for i, val in enumerate(obj.getEffectiveColliderEfficiencies(spcs)) if val != 1}
    elif isinstance(obj, Troe):
        D["type"] = "Troe"
        D["arrhigh"] = obj2dict(obj.arrheniusHigh, spcs)
        D["arrlow"] = obj2dict(obj.arrheniusLow, spcs)
        D["efficiencies"] = {spcs[i].label: float(val) for i, val in enumerate(obj.getEffectiveColliderEfficiencies(spcs)) if val != 1}
        D["a"] = obj.alpha
        D["T1"] = obj.T1.value_si
        if obj.T2:
            D["T2"] = obj.T2.value_si
        else:
            D["T2"] = 0.0
        D["T3"] = obj.T3.value_si
    elif isinstance(obj, Chebyshev):
        D["type"] = "Chebyshev"
        D["coefs"] = obj.coeffs.value_si.tolist()
        D["Tmin"] = obj.Tmin.value_si
        D["Tmax"] = obj.Tmax.value_si
        D["Pmin"] = obj.Pmin.value_si
        D["Pmax"] = obj.Pmax.value_si
    elif isinstance(obj, Wilhoit):
        D["type"] = "Wilhoit"
        D["coefs"] = [obj.a0, obj.a1, obj.a2, obj.a3]
        D["Cp0"] = obj.Cp0.value_si
        D["Cpinf"] = obj.CpInf.value_si
        D["H0"] = obj.H0.value_si
        D["S0"] = obj.S0.value_si
        D["B"] = obj.B.value_si
    elif isinstance(obj, SolventData):
        D["type"] = "Solvent"
        D["name"] = label
        dsub = dict()
        dsub["type"] = "RiedelViscosity"
        dsub["A"] = float(obj.A)
        dsub["B"] = float(obj.B)
        dsub["C"] = float(obj.C)
        dsub["D"] = float(obj.D)
        dsub["E"] = float(obj.E)
        D["mu"] = dsub
    elif obj is None:
        return None
    else:
        raise ValueError("Object of type {} does not have a defined conversion to ReactionMechanismSimulator format".format(type(obj)))
    return D


class RMSWriter(object):
    """
    This class listens to a RMG subject
    and writes an rms file with the current state of the RMG model,
    to a rms subfolder.


    A new instance of the class can be appended to a subject as follows:

    rmg = ...
    listener = RMSWriter(outputDirectory)
    rmg.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    rmg.detach(listener)

    """
    def __init__(self, outputDirectory=''):
        super(RMSWriter, self).__init__()
        self.outputDirectory = outputDirectory
        makeOutputSubdirectory(outputDirectory, 'rms')

    def update(self, rmg):
        solventData = None
        if rmg.solvent:
            solventData = rmg.database.solvation.get_solvent_data(rmg.solvent)
        writeyml(rmg.reactionModel.core.species, rmg.reactionModel.core.reactions, solvent=rmg.solvent, solventData=solventData,
                 path=os.path.join(self.outputDirectory, 'rms', 'chem{}.rms').format(len(rmg.reactionModel.core.species)))
