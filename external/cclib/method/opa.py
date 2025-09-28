"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 739 $"

import random

import numpy

from calculationmethod import Method


def func(x):
    if x==1:
        return 1
    else:
        return x+func(x-1)


class OPA(Method):
    """The overlap population analysis."""
    
    def __init__(self, *args):

        # Call the __init__ method of the superclass.
        super(OPA, self).__init__(logname="OPA", *args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "OPA of" % (self.data)

    def __repr__(self):
        """Return a representation of the object."""
        return 'OPA("%s")' % (self.data)
    
    def calculate(self, indices=None, fupdate=0.05):
        """Perform an overlap population analysis given the results of a parser"""
    
        # Do we have the needed info in the ccData object?
        if not hasattr(self.data, "mocoeffs") \
          and not ( hasattr(self.data, "aooverlaps") \
                    or hasattr(self.data, "fooverlaps") ) \
          and not hasattr(self.data, "nbasis"):
            self.logger.error("Missing mocoeffs, aooverlaps/fooverlaps or nbasis")
            return False #let the caller of function know we didn't finish

        if not indices:

            # Build list of groups of orbitals in each atom for atomresults.
            if hasattr(self.data, "aonames"):
                names = self.data.aonames
            elif hasattr(self.data, "foonames"):
                names = self.data.fonames

            atoms = []
            indices = []

            name = names[0].split('_')[0]
            atoms.append(name)
            indices.append([0])

            for i in range(1, len(names)):
                name = names[i].split('_')[0]
                try:
                    index = atoms.index(name)
                except ValueError: #not found in atom list
                    atoms.append(name)
                    indices.append([i])
                else:
                    indices[index].append(i)

        # Determine number of steps, and whether process involves beta orbitals.
        nfrag = len(indices) #nfrag
        nstep = func(nfrag - 1)
        unrestricted = (len(self.data.mocoeffs) == 2)
        alpha = len(self.data.mocoeffs[0])
        nbasis = self.data.nbasis

        self.logger.info("Creating attribute results: array[4]")
        results= [ numpy.zeros([nfrag, nfrag, alpha], "d") ]
        if unrestricted:
            beta = len(self.data.mocoeffs[1])
            results.append(numpy.zeros([nfrag, nfrag, beta], "d"))
            nstep *= 2
            
        if hasattr(self.data, "aooverlaps"):
            overlap = self.data.aooverlaps
        elif hasattr(self.data,"fooverlaps"):
            overlap = self.data.fooverlaps

        #intialize progress if available
        if self.progress:
            self.progress.initialize(nstep)

        size = len(self.data.mocoeffs[0])
        step = 0

        preresults = []
        for spin in range(len(self.data.mocoeffs)):
            two = numpy.array([2.0]*len(self.data.mocoeffs[spin]),"d")


            # OP_{AB,i} = \sum_{a in A} \sum_{b in B} 2 c_{ai} c_{bi} S_{ab}

            for A in range(len(indices)-1):

                for B in range(A+1, len(indices)):

                    if self.progress: #usually only a handful of updates, so remove random part
                        self.progress.update(step, "Overlap Population Analysis")

                    for a in indices[A]:

                        ca = self.data.mocoeffs[spin][:,a]

                        for b in indices[B]:
                            
                            cb = self.data.mocoeffs[spin][:,b]
                            temp = ca * cb * two *overlap[a,b]
                            results[spin][A,B] = numpy.add(results[spin][A,B],temp)
                            results[spin][B,A] = numpy.add(results[spin][B,A],temp)

                    step += 1

        temparray2 = numpy.swapaxes(results[0],1,2)
        self.results = [ numpy.swapaxes(temparray2,0,1) ]
        if unrestricted:
            temparray2 = numpy.swapaxes(results[1],1,2)
            self.results.append(numpy.swapaxes(temparray2, 0, 1))

        if self.progress:
            self.progress.update(nstep, "Done")

        return True


if __name__ == "__main__":
    import doctest, opa
    doctest.testmod(opa, verbose=False)
