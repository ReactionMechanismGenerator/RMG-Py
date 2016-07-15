
import re
import logging
from rmgpy.quantity import Energy, Mass, Length, Frequency

class QMData(object):
    """
    General class for data extracted from a QM calculation
    """
    def __init__(self,
                 groundStateDegeneracy = -1,
                 numberOfAtoms = None,
                 stericEnergy = None,
                 molecularMass = None,
                 energy = 0,
                 atomicNumbers = None,
                 rotationalConstants = None,
                 atomCoords = None,
                 frequencies = None,
                 source = None,
                 ):
        #: Electronic ground state degeneracy in RMG taken as number of radicals +1
        self.groundStateDegeneracy = groundStateDegeneracy
        self.numberOfAtoms = numberOfAtoms #: Number of atoms.
        self.stericEnergy = Energy(stericEnergy)
        self.molecularMass = Mass(molecularMass)
        self.energy = Energy(energy)
        self.atomicNumbers = atomicNumbers
        self.rotationalConstants = Frequency(rotationalConstants)
        self.atomCoords = Length(atomCoords)
        self.frequencies = Frequency(frequencies)
        self.source = source
        
        self.testValid()

    def testValid(self):
        assert self.groundStateDegeneracy > 0
    
    def __repr__(self):
        things=[]
        for attribute in ['groundStateDegeneracy',
                          'numberOfAtoms',
                          'stericEnergy',
                          'molecularMass',
                          'energy',
                          'atomicNumbers',
                          'rotationalConstants',
                          'atomCoords',
                          'frequencies',
                          'source',
                          ]:
            things.append("{0!s}={1!r}".format(attribute, getattr(self,attribute)))
        string = ', '.join(things)
        string = re.sub('\s+',' ',string)
        return 'QMData({0!s})'.format(string)
        
def parseCCLibData(cclibData, groundStateDegeneracy):
    """
    Parses a CCLib data object and returns QMData object.
    """
    try:
        numberOfAtoms = cclibData.natom
        molecularMass = (cclibData.molmass, 'amu')
        energy = (cclibData.scfenergies[-1], 'eV/molecule')
        atomicNumbers = cclibData.atomnos
        rotationalConstants = ([i * 1e9 for i in cclibData.rotcons[-1]],'hertz')
        atomCoords = (cclibData.atomcoords[-1], 'angstrom')
        frequencies = (cclibData.vibfreqs, 'cm^-1')

    except AttributeError, e:
        logging.error("The passed in cclibData has these attributes: {0!r}".format(cclibData._attrlist))
        raise e

    if hasattr(cclibData, 'stericenergy'):
        stericEnergy = (cclibData.stericenergy, 'eV/molecule')
    else:
        stericEnergy = None

    return QMData(
        groundStateDegeneracy=groundStateDegeneracy,
        numberOfAtoms=numberOfAtoms,
        stericEnergy=stericEnergy,
        molecularMass=molecularMass,
        energy=energy,
        atomicNumbers=atomicNumbers,
        rotationalConstants=rotationalConstants,
        atomCoords=atomCoords,
        frequencies=frequencies
    )
