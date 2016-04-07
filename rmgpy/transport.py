'''
Created on Feb 26, 2013

@author: Jake
'''
from rmgpy import quantity
from rmgpy.quantity import DipoleMoment, Energy, Length, Volume
import rmgpy.constants as constants
import numpy

class TransportData:
    """
    A set of transport properties.
    
    The attributes are:
    =================  ============================================================
    Attribute          Description
    =================  ============================================================
    `shapeIndex`        0 for monoatomic, 1 for linear, 2 for nonlinear
    `epsilon`           Lennard-Jones well depth
    `sigma`             Lennard-Jones collision diameter
    `dipoleMoment`      Dipole Moment
    `polarizability`    Polarizability Volume
    `rotrelaxcollnum`   Rotational relaxation number at 298 K, saved as a double.
    =================  ============================================================
    
    """

    def __init__(self, shapeIndex=None, epsilon=None, sigma=None, dipoleMoment=None, polarizability=None, rotrelaxcollnum=None, comment = ''):
        self.shapeIndex = shapeIndex
        try:
            self.epsilon = Energy(epsilon)
        except quantity.QuantityError:
                self.epsilon = quantity.Temperature(epsilon)
                self.epsilon.value_si *= constants.R
                self.epsilon.units = 'kJ/mol'
        self.sigma = Length(sigma)
        self.dipoleMoment = DipoleMoment(dipoleMoment)
        self.polarizability = Volume(polarizability)
        self.rotrelaxcollnum = rotrelaxcollnum
        self.comment = comment

    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        TransportData object.
        """
        attributes=[]
        if self.shapeIndex is not None:
            attributes.append('shapeIndex={0!r}'.format(self.shapeIndex))
        if self.epsilon is not None:
            attributes.append('epsilon={0!r}'.format(self.epsilon))
        if self.sigma is not None:
            attributes.append('sigma={0!r}'.format(self.sigma))
        if self.dipoleMoment is not None:
            attributes.append('dipoleMoment={0!r}'.format(self.dipoleMoment))
        if self.polarizability is not None:
            attributes.append('polarizability={0!r}'.format(self.polarizability))
        if self.rotrelaxcollnum is not None:
            attributes.append('rotrelaxcollnum={0!r}'.format(self.rotrelaxcollnum))
        if self.comment:
            attributes.append('comment="""{0!s}"""'.format(self.comment))
        string = 'TransportData({0!s})'.format(', '.join(attributes))
        return string

    def __reduce__(self):
        """
        A helper function used when picking a TransportData object.
        """
        return (TransportData, (self.shapeIndex, self.epsilon, self.sigma, self.dipoleMoment, self.polarizability, self.rotrelaxcollnum, self.comment))

    def getCollisionFrequency(self, T, M, mu):
        """
        Return the value of the Lennard-Jones collision frequency in Hz at the
        given temperature `T` in K for colliders with the given concentration
        `M` in mol/m^3 and reduced mass `mu` in amu.
        
        This seems to also exist in rmgpy.pdep.configuration.calculateCollisionFrequency
        Why the redundancy?
        """
        sigma = self.sigma.value_si
        epsilon = self.epsilon.value_si
        M *= constants.Na       # mol/m^3 -> molecules/m^3
        Tred = constants.R*T / epsilon
        omega22 = 1.16145 * Tred**(-0.14874) + 0.52487 * numpy.exp(-0.77320 * Tred) + 2.16178 * numpy.exp(-2.43787 * Tred)
        mu *= constants.amu
        return omega22 * numpy.sqrt(8 * constants.kB * T / constants.pi / mu) * constants.pi * sigma * sigma * M
    
    def toCantera(self):
        """
        Returns a Cantera GasTransportData object.
    
        The Cantera usage is as follows:
        
        GasTransportData().set_customary_units(self, geometry, diameter, well_depth, dipole=0.0, polarizability=0.0, rotational_relaxation=0.0, acentric_factor=0.0)
        Set the parameters using customary units: diameter in Angstroms, well depth in Kelvin, dipole in Debye, rotational relaxiation at 298 K, and polarizability in Angstroms^3. 
        These are the units used in in CK-style input files.
        """
        import cantera as ct
        
        ctTransport = ct.GasTransportData()
        
        
        if self.shapeIndex == 0:
            geometry = 'atom'
        elif self.shapeIndex == 1:
            geometry = 'linear'
        elif self.shapeIndex == 2:
            geometry = 'nonlinear' 
        
        # collision diameter in angstroms
        diameter = self.sigma.value_si * 1e10
        # Well depth in Kelvins
        well_depth = self.epsilon.value_si / constants.R   
        # Dipole in debye
        dipole = self.dipoleMoment.value_si * constants.c * 1e21 if self.dipoleMoment else 0.0
        # polarizability in cubic angstroms
        polarizability = self.polarizability.value_si * 1e30 if self.polarizability else 0.0
        rotational_relaxation = self.rotrelaxcollnum if self.rotrelaxcollnum else 0.0
        acentric_factor = 0.0
        
        ctTransport.set_customary_units(geometry, diameter, well_depth, dipole, polarizability, rotational_relaxation, acentric_factor)

        return ctTransport
