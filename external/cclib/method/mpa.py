"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 733 $"

import random

import numpy

from population import Population


class MPA(Population):
    """The Mulliken population analysis."""
    
    def __init__(self, *args):

        # Call the __init__ method of the superclass.
        super(MPA, self).__init__(logname="MPA", *args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "MPA of" % (self.data)

    def __repr__(self):
        """Return a representation of the object."""
        return 'MPA("%s")' % (self.data)
    
    def calculate(self, indices=None, fupdate=0.05):
        """Perform a Mulliken population analysis."""
    
        # Do we have the needed attributes in the data object?
        if not hasattr(self.data, "mocoeffs"):
            self.logger.error("Missing mocoeffs")
            return False
        if not (hasattr(self.data, "aooverlaps") \
             or hasattr(self.data, "fooverlaps") ):
            self.logger.error("Missing overlap matrix")
            return False
        if not hasattr(self.data, "nbasis"):
            self.logger.error("Missing nbasis")
            return False
        if not hasattr(self.data, "homos"):
            self.logger.error("Missing homos")
            return False


        # Determine number of steps, and whether process involves beta orbitals.
        self.logger.info("Creating attribute aoresults: [array[2]]")
        nbasis = self.data.nbasis
        alpha = len(self.data.mocoeffs[0])
        self.aoresults = [ numpy.zeros([alpha, nbasis], "d") ]
        nstep = alpha
        unrestricted = (len(self.data.mocoeffs) == 2)
        if unrestricted:
            beta = len(self.data.mocoeffs[1])
            self.aoresults.append(numpy.zeros([beta, nbasis], "d"))
            nstep += beta

        # Intialize progress if available.
        if self.progress:
            self.progress.initialize(nstep)

        step = 0
        for spin in range(len(self.data.mocoeffs)):

            for i in range(len(self.data.mocoeffs[spin])):

                if self.progress and random.random() < fupdate:
                    self.progress.update(step, "Mulliken Population Analysis")

                #X_{ai} = \sum_b c_{ai} c_{bi} S_{ab}
                #       = c_{ai} \sum_b c_{bi} S_{ab}
                #       = c_{ai} C(i) \cdot S(a)
                # X = C(i) * [C(i) \cdot S]
                # C(i) is 1xn and S is nxn, result of matrix mult is 1xn

                ci = self.data.mocoeffs[spin][i]
                if hasattr(self.data, "aooverlaps"):
                    temp = numpy.dot(ci, self.data.aooverlaps)
                elif hasattr(self.data, "fooverlaps"):
                    temp = numpy.dot(ci, self.data.fooverlaps)

                self.aoresults[spin][i] = numpy.multiply(ci, temp).astype("d")

                step += 1

        if self.progress:
            self.progress.update(nstep, "Done")

        retval = super(MPA, self).partition(indices)

        if not retval:
            self.logger.error("Error in partitioning results")
            return False

        # Create array for mulliken charges.
        self.logger.info("Creating fragcharges: array[1]")
        size = len(self.fragresults[0][0])
        self.fragcharges = numpy.zeros([size], "d")
        
        for spin in range(len(self.fragresults)):

            for i in range(self.data.homos[spin] + 1):

                temp = numpy.reshape(self.fragresults[spin][i], (size,))
                self.fragcharges = numpy.add(self.fragcharges, temp)
        
        if not unrestricted:
            self.fragcharges = numpy.multiply(self.fragcharges, 2)

        return True

if __name__ == "__main__":
    import doctest, mpa
    doctest.testmod(mpa, verbose=False)
