
import re
import logging
from rmgpy.quantity import Energy, Mass, Length, Frequency

class QMData:
    """
    General class for data extracted from a QM calculation
    """
    def __init__(self,
                 groundStateDegeneracy = None,
                 numberOfAtoms = None,
                 stericEnergy = None,
                 molecularMass = None,
                 energy = 0,
                 atomicNumbers = [],
                 rotationalConstants = [],
                 atomCoords = [],
                 frequencies = [],
                 source = None,
                 ):
        #: Electronic ground state degeneracy in RMG taken as number of radicals +1
        self.groundStateDegeneracy = groundStateDegeneracy
        self.numberOfAtoms = numberOfAtoms #: Number of atoms.
        self.stericEnergy = Energy(stericEnergy)
        """
        Steric energy.
        """
        self.molecularMass = Mass(molecularMass)
        self.energy = Energy(energy)
        self.atomicNumbers = atomicNumbers
        #: Rotational constants, in Hz.
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
        
class CCLibData(QMData):
    """
    QM Data extracted from a cclib data object
    
    This data objects collects information that CCLib was able to
    retrieve from a quantum chemistry output file.
    """

    def __init__(self,
                 cclib_data,
                 groundStateDegeneracy
                 ):
        
        #: data object returned by a parsing tool like CCLib.parse()
        self.cclib_data = cclib_data
        try:
            numberOfAtoms = cclib_data.natom
            molecularMass = (cclib_data.molmass,'amu')
            energy = (cclib_data.scfenergies[-1],'eV/molecule') # final optimized PM3 energy (cclib gives this in eV)
            atomicNumbers = cclib_data.atomnos
            rotationalConstants = ([i * 1e9 for i in cclib_data.rotcons[-1]],'hertz') # Hz. #print the final rotational constants (note that ideally we would use next to last value ([-2]) as this has more significant digits and is for the same geometry, but there is a complication for linear molecules (labeled as "Rotational constant" rather than "...constants"...there might be some ways around this like looking for "Rotational constant" string instead, but it is probably not a big deal to just use rounded values
            atomCoords = (cclib_data.atomcoords[-1],"angstrom") # I assume Angstrom (not Bohr?)
            frequencies = (cclib_data.vibfreqs, "cm^-1") # 1/cm

        except AttributeError, e:
            logging.error("The passed in cclib_data has these attributes: {0!r}".format(cclib_data._attrlist))
            raise e
        
        QMData.__init__( self,
            groundStateDegeneracy = groundStateDegeneracy,
            numberOfAtoms = numberOfAtoms,
            molecularMass = molecularMass,
            energy = energy,
            atomicNumbers = atomicNumbers,
            rotationalConstants = rotationalConstants,
            atomCoords = atomCoords,
            frequencies = frequencies,
            )
        if hasattr(cclib_data, 'stericenergy'):
            # steric energy
            self.stericEnergy = (cclib_data.stericenergy, 'eV/molecule') # this is eV.  /27.2113845 for Hartrees
