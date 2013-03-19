'''
Created on Feb 26, 2013

@author: Jake
'''
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

    def __init__(self, shapeIndex=None, epsilon=None, sigma=None, dipoleMoment=None, polarizability=None, rotrelaxcollnum=None):
        self.shapeIndex = shapeIndex
        self.epsilon = epsilon
        self.sigma = sigma
        self.dipoleMoment = dipoleMoment
        self.polarizability = polarizability
        self.rotrelaxcollnum = rotrelaxcollnum
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        TransportData object.
        """
        string = 'TransportData(shapeIndex={0!r}, epsilon={1!r}, sigma={2!r}, dipoleMoment={3!r}, polarizability={4!r}, rotrelaxcollnum={5!r}'.format(self.shapeIndex, self.epsilon, self.sigma, self.dipoleMoment, self.polarizability, self.rotrelaxcollnum)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a TransportData object.
        """
        return (TransportData, (self.shapeIndex, self.epsilon, self.sigma, self.dipoleMoment, self.polarizability, self.rotrelaxcollnum))

    property shapeIndex:
        """."""
        def __get__(self):
            return self._shapeIndex
        def _set_(self,value):
            self.shapeIndex = value
        
    property epsilon:
        """."""
        def __get__(self):
            return self._epsilon
        def _set_(self,value):
            self.epsilon = value 
            
    property sigma:
        """."""
        def __get__(self):
            return self._sigma
        def _set_(self,value):
            self.sigma = value
             
    property dipoleMoment:
        """."""
        def __get__(self):
            return self._dipoleMoment
        def _set_(self,value):
            self.dipoleMoment = value
             
    property polarizability:
        """."""
        def __get__(self):
            return self._polarizability
        def _set_(self,value):
            self.polarizability = value
             
    property rotrelaxcollnum:
        """."""
        def __get__(self):
            return self._rotrelaxcollnum
        def _set_(self,value):
            self.rotrelaxcollnum = value
