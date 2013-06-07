"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 733 $"

import logging

import numpy

from calculationmethod import Method


class Population(Method):
    """A base class for all population-type methods."""
    
    def __init__(self, data, progress=None, \
                 loglevel=logging.INFO, logname="Log"):

        # Call the __init__ method of the superclass.
        super(Population, self).__init__(data, progress, loglevel, logname)
        self.fragresults = None
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Population"

    def __repr__(self):
        """Return a representation of the object."""
        return "Population"
    
    def partition(self, indices=None):

        if not hasattr(self, "aoresults"):
            self.calculate()

        if not indices:

            # Build list of groups of orbitals in each atom for atomresults.
            if hasattr(self.data, "aonames"):
                names = self.data.aonames
            elif hasattr(self.data, "fonames"):
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

        natoms = len(indices)
        nmocoeffs = len(self.aoresults[0])
        
        # Build results numpy array[3].
        alpha = len(self.aoresults[0])
        results = []
        results.append(numpy.zeros([alpha, natoms], "d"))

        if len(self.aoresults) == 2:
            beta = len(self.aoresults[1])
            results.append(numpy.zeros([beta, natoms], "d"))
        
        # For each spin, splice numpy array at ao index,
        #   and add to correct result row.
        for spin in range(len(results)):

            for i in range(natoms): # Number of groups.

                for j in range(len(indices[i])): # For each group.
                
                    temp = self.aoresults[spin][:, indices[i][j]]
                    results[spin][:, i] = numpy.add(results[spin][:, i], temp)

        self.logger.info("Saving partitioned results in fragresults: [array[2]]")
        self.fragresults = results

        return True


if __name__ == "__main__":
    import doctest, population
    doctest.testmod(population, verbose=False)
