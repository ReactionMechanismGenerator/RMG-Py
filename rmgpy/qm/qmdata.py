class CCLibData:
    """
    This data objects collects information that CCLib was able to
    retrieve from a quantum chemistry output file.
    @author nmvdewie
    """

    def __init__(self, cclib_data, groundStateDegeneracy = -1, numberOfAtoms = -1, stericEnergy = -1, molecularMass = -1, energy = 0, atomicNumbers = [], rotationalConstants = [], atomCoords = [], frequencies = []):
        
        """
        * Electronic ground state degeneracy
        * 
        * in RMG taken as number of radicals +1
        """
        self.groundStateDegeneracy = groundStateDegeneracy
       
        self.cclib_data = cclib_data#data object returned by a parsing tool like CCLib.parse()
        
        self.numberOfAtoms = cclib_data.natom
        if hasattr(cclib_data, 'stericenergy'):
            self.stericEnergy = cclib_data.stericenergy/27.2113845#steric energy (in Hartree)

        
        self.molecularMass = cclib_data.molmass #print the molecular mass (in amu)
        
        self.energy = cclib_data.scfenergies[-1]/27.2113845 #print the final optimized PM3 energy (cclib gives this in eV, but I have converted here to Hartree); conversion value taken from cclib code; compare CODATA 2006 value 27.21138386(68) eV/Hartree
        
        self.atomicNumbers = cclib_data.atomnos
        
        self.rotationalConstants = [i * 1000000000 for i in cclib_data.rotcons[-1]] #print the final rotational constants (note that ideally we would use next to last value ([-2]) as this has more significant digits and is for the same geometry, but there is a complication for linear molecules (labeled as "Rotational constant" rather than "...constants"...there might be some ways around this like looking for "Rotational constant" string instead, but it is probably not a big deal to just use rounded values
        
        self.atomCoords = cclib_data.atomcoords[-1]
        
        self.frequencies = cclib_data.vibfreqs
