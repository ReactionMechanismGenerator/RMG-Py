"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 733 $"

import random # For sometimes running the progress updater

import numpy

from population import Population


class CSPA(Population):
    """The C-squared population analysis."""
    
    def __init__(self, *args):

        # Call the __init__ method of the superclass.
        super(CSPA, self).__init__(logname="CSPA", *args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "CSPA of" % (self.data)

    def __repr__(self):
        """Return a representation of the object."""
        return 'CSPA("%s")' % (self.data)
    
    def calculate(self, indices=None, fupdate=0.05):
        """Perform the C squared population analysis.
        
        Inputs:
           indices - list of lists containing atomic orbital indices of fragments
        """
    
        # Do we have the needed info in the parser?
        if not hasattr(self.data, "mocoeffs"):
            self.logger.error("Missing mocoeffs")
            return False
        if not hasattr(self.data, "nbasis"):
            self.logger.error("Missing nbasis")
            return False
        if not hasattr(self.data, "homos"):
            self.logger.error("Missing homos")
            return False

        self.logger.info("Creating attribute aoresults: array[3]")

        # Determine number of steps, and whether process involves beta orbitals.
        unrestricted = (len(self.data.mocoeffs)==2)
        nbasis = self.data.nbasis
        self.aoresults = []
        alpha = len(self.data.mocoeffs[0])
        self.aoresults.append(numpy.zeros([alpha, nbasis], "d"))
        nstep = alpha
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
                    self.progress.update(step, "C^2 Population Analysis")

                submocoeffs = self.data.mocoeffs[spin][i]
                scale = numpy.inner(submocoeffs, submocoeffs)
                tempcoeffs = numpy.multiply(submocoeffs, submocoeffs)
                tempvec = tempcoeffs/scale
                self.aoresults[spin][i] = numpy.divide(tempcoeffs, scale).astype("d")

                step += 1

        if self.progress:
            self.progress.update(nstep, "Done")

        retval = super(CSPA, self).partition(indices)

        if not retval:
            self.logger.error("Error in partitioning results")
            return False

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
    import doctest, cspa
    doctest.testmod(cspa, verbose=False)
