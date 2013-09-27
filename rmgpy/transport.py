'''
Created on Feb 26, 2013

@author: Jake
'''

from rmgpy.quantity import DipoleMoment, Energy, Length, Polarizability
import rmgpy.constants as constants
import numpy

class TransportData:
    """
    A heat capacity model based on a set of discrete heat capacity data points.
    The attributes are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `shapeIndex`         
    `epsilon`        
    `sigma`          
    `dipoleMoment`          
    `polarizability`           
    `rotrelaxcollnum`         
    =============== ============================================================
    
    """

    def __init__(self, shapeIndex=None, epsilon=None, sigma=None, dipoleMoment=None, polarizability=None, rotrelaxcollnum=None, comment = ''):
        self.shapeIndex = shapeIndex
        self.epsilon = Energy(epsilon)
        self.sigma = Length(sigma)
        self.dipoleMoment = DipoleMoment(dipoleMoment)
        self.polarizability = Polarizability(polarizability)
        self.rotrelaxcollnum = rotrelaxcollnum
        self.comment = comment
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        TransportData object.
        """
        string = 'TransportData(shapeIndex={0!r}, epsilon={1!r}, sigma={2!r}, dipoleMoment={3!r}, polarizability={4!r}, rotrelaxcollnum={5!r}'.format(self.shapeIndex, self.epsilon, self.sigma, self.dipoleMoment, self.polarizability, self.rotrelaxcollnum)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string
    
    def getLennardJones(self):
        """
        Return a string representation that can be used for collision Frequencies
        """
        string = 'sigma={0!r}, epsilon={1!r}'.format(self.sigma, self.epsilon)

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
        """
        sigma = self.sigma.value_si
        epsilon = self.epsilon.value_si
        M *= constants.Na       # mol/m^3 -> molecules/m^3
        Tred = constants.R*T / epsilon
        omega22 = 1.16145 * Tred**(-0.14874) + 0.52487 * numpy.exp(-0.77320 * Tred) + 2.16178 * numpy.exp(-2.43787 * Tred)
        mu *= constants.amu
        return omega22 * numpy.sqrt(8 * constants.kB * T / constants.pi / mu) * constants.pi * sigma * sigma * M