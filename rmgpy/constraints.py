#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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

from rmgpy.species import Species

def failsSpeciesConstraints(species, verbose=False):
    """
    Pass in either a `Species` or `Molecule` object and checks whether it passes 
    the speciesConstraints set by the user.  If not, returns `True` for failing speciesConstraints.


    If `verbose` is set to ``True``, debugging information about why
    the molecule is forbidden will be logged.
    """
    
    from rmgpy.rmg.input import getInput

    try:
        speciesConstraints = getInput('speciesConstraints')
    except Exception, e:
        logging.debug('Species constraints could not be found.')
        speciesConstraints = {}
    

    explicitlyAllowedMolecules = speciesConstraints.get('explicitlyAllowedMolecules', [])
    maxCarbonAtoms = speciesConstraints.get('maximumCarbonAtoms', 1000000)
    maxHydrogenAtoms = speciesConstraints.get('maximumHydrogenAtoms', 1000000)
    maxOxygenAtoms = speciesConstraints.get('maximumOxygenAtoms', 1000000)
    maxNitrogenAtoms = speciesConstraints.get('maximumNitrogenAtoms', 1000000)
    maxSiliconAtoms = speciesConstraints.get('maximumSiliconAtoms', 1000000)
    maxSulfurAtoms = speciesConstraints.get('maximumSulfurAtoms', 1000000)
    maxHeavyAtoms = speciesConstraints.get('maximumHeavyAtoms', 1000000)
    maxRadicals = speciesConstraints.get('maximumRadicalElectrons', 1000000)
    
    if isinstance(species, Species):
        struct = species.molecule[0]
    else:
        # expects a molecule here
        struct = species
    for molecule in explicitlyAllowedMolecules:
        if struct.isIsomorphic(molecule):
            return False        
    H = struct.getNumAtoms('H')
    if struct.getNumAtoms('C') > maxCarbonAtoms:
        if verbose: 
            logging.error("Species {0} exceeded 'maximumCarbonAtoms = {1}' species constraint.".format(struct.toSMILES(), maxCarbonAtoms))
        return True
    if H > maxHydrogenAtoms:
        if verbose: 
            logging.error("Species {0} exceeded 'maximumHydrogenAtoms = {1}' species constraint.".format(struct.toSMILES(), maxHydrogenAtoms))
        return True
    if struct.getNumAtoms('O') > maxOxygenAtoms:
        if verbose: 
            logging.error("Species {0} exceeded 'maximumOxygenAtoms = {1}' species constraint.".format(struct.toSMILES(), maxOxygenAtoms))
        return True
    if struct.getNumAtoms('N') > maxNitrogenAtoms:
        if verbose:
            logging.error("Species {0} exceeded 'maximumNitrogenAtoms = {1}' species constraint.".format(struct.toSMILES(), maxNitrogenAtoms))
        return True
    if struct.getNumAtoms('Si') > maxSiliconAtoms:
        if verbose:
            logging.error("Species {0} exceeded 'maximumSiliconAtoms = {1}' species constraint.".format(struct.toSMILES(), maxSiliconAtoms))
        return True
    if struct.getNumAtoms('S') > maxSulfurAtoms:
        if verbose:
            logging.error("Species {0} exceeded 'maximumSulfurAtoms = {1}' species constraint.".format(struct.toSMILES(), maxSulfurAtoms))
        return True
    if len(struct.atoms) - H > maxHeavyAtoms:
        if verbose: 
            logging.error("Species {0} exceeded 'maximumHeavyAtoms = {1}' species constraint.".format(struct.toSMILES(), maxHeavyAtoms))
        return True
    if (struct.getNumberOfRadicalElectrons() > maxRadicals):
        if verbose: 
            logging.error("Species {0} exceeded 'maximumRadicalElectrons = {1}' species constraint.".format(struct, maxRadicals))
        return True
    return False