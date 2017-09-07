#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import logging
from numpy import isclose
from rmgpy.molecule.element import getElement
from rmgpy.species import Species

def failsSpeciesConstraints(species):
    """
    Pass in either a `Species` or `Molecule` object and checks whether it passes 
    the speciesConstraints set by the user.  If not, returns `True` for failing speciesConstraints.
    """
    
    from rmgpy.rmg.input import getInput

    try:
        speciesConstraints = getInput('speciesConstraints')
    except Exception:
        logging.debug('Species constraints could not be found.')
        speciesConstraints = {}
    
    if isinstance(species, Species):
        struct = species.molecule[0]
    else:
        # expects a molecule here
        struct = species

    explicitlyAllowedMolecules = speciesConstraints.get('explicitlyAllowedMolecules', [])
    for molecule in explicitlyAllowedMolecules:
        if struct.isIsomorphic(molecule):
            return False  
    
    maxCarbonAtoms = speciesConstraints.get('maximumCarbonAtoms', -1)          
    if maxCarbonAtoms != -1:
        if struct.getNumAtoms('C') > maxCarbonAtoms:
            return True

    maxOxygenAtoms = speciesConstraints.get('maximumOxygenAtoms', -1)
    if maxOxygenAtoms != -1:
        if struct.getNumAtoms('O') > maxOxygenAtoms:
            return True

    maxNitrogenAtoms = speciesConstraints.get('maximumNitrogenAtoms', -1)
    if maxNitrogenAtoms != -1:
        if struct.getNumAtoms('N') > maxNitrogenAtoms:
            return True

    maxSiliconAtoms = speciesConstraints.get('maximumSiliconAtoms', -1)
    if maxSiliconAtoms != -1:
        if struct.getNumAtoms('Si') > maxSiliconAtoms:
            return True

    maxSulfurAtoms = speciesConstraints.get('maximumSulfurAtoms', -1)
    if maxSulfurAtoms != -1:
        if struct.getNumAtoms('S') > maxSulfurAtoms:
            return True

    maxHeavyAtoms = speciesConstraints.get('maximumHeavyAtoms', -1)
    if maxHeavyAtoms != -1:
        if struct.getNumAtoms() - struct.getNumAtoms('H') > maxHeavyAtoms:
            return True

    maxRadicals = speciesConstraints.get('maximumRadicalElectrons', -1)
    if maxRadicals != -1:
        if (struct.getRadicalCount() > maxRadicals):
            return True

    maxCarbenes = speciesConstraints.get('maximumSingletCarbenes', 1)
    if maxRadicals != -1:
        if struct.getSingletCarbeneCount() > maxCarbenes:
            return True

    maxCarbeneRadicals = speciesConstraints.get('maximumCarbeneRadicals', 0)
    if maxCarbeneRadicals != -1:
        if struct.getSingletCarbeneCount() > 0 and struct.getRadicalCount() > maxCarbeneRadicals:
            return True

    maxIsotopes = speciesConstraints.get('maximumIsotopicAtoms', -1)
    if maxIsotopes != -1:
        counter = 0
        for atom in struct.atoms:
            if not isclose(atom.mass, getElement(atom.symbol).mass, atol=1e-04):
                counter += 1
            if counter > maxIsotopes: return True

    return False
