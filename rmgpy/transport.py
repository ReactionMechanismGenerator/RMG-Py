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

    def __init__(self, shapeIndex=None, epsilon=None, sigma=None, dipoleMoment=None, polarizability=None, rotrelaxcollnum=None, comment = ''):
        self.shapeIndex = shapeIndex
        self.epsilon = epsilon
        self.sigma = sigma
        self.dipoleMoment = dipoleMoment
        self.polarizability = polarizability
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

    def __reduce__(self):
        """
        A helper function used when picking a TransportData object.
        """
        return (TransportData, (self.shapeIndex, self.epsilon, self.sigma, self.dipoleMoment, self.polarizability, self.rotrelaxcollnum, self.comment))
