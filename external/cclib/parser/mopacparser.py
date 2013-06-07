"""
gmagoon 07/06/09: new class for MOPAC parsing, based on gaussianparser.py from cclib, described below:
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 814 $"


#import re

import numpy

import utils
import logfileparser


def symbol2int(symbol):
    t = utils.PeriodicTable()
    return t.number[symbol]
    #if symbol == 'C': return 6
    #elif symbol == 'H': return 1
    #elif symbol == 'O': return 8
    #elif symbol == 'Si': return 14
    #else: return -1

class Mopac(logfileparser.Logfile):
    """A MOPAC2009 output file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Mopac, self).__init__(logname="Mopac", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Mopac log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Mopac("%s")' % (self.filename)

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""
        
        # Number of atoms. (I think this section of code may be redundant and not needed)
        # Example:            Empirical Formula: C H2 O  =     4 atoms
        if line.find("Empirical Formula:") > -1:

            self.updateprogress(inputfile, "Attributes", self.fupdate)
            #locate the component that beg        
            natom = int(line.split()[-2]) #second to last component should be number of atoms (last element is "atoms" (or possibly "atom"?))
            if hasattr(self, "natom"):
                assert self.natom == natom
            else:
                self.natom = natom
        
        # Extract the atomic numbers and coordinates from the optimized geometry
        # note that cartesian coordinates section occurs multiple times in the file, and we want to end up using the last instance
        # also, note that the section labeled cartesian coordinates doesn't have as many decimal places as the one used here
        # Example 1 (not used):
#          CARTESIAN COORDINATES 
#
#    NO.       ATOM               X         Y         Z
#
#     1         O                  4.7928   -0.8461    0.3641
#     2         O                  5.8977   -0.3171    0.0092
#     3         C                  3.8616    0.0654    0.8629
#     4         O                  2.9135    0.0549   -0.0719
#     5        Si                 -0.6125   -0.0271    0.0487
#     6         O                  0.9200    0.2818   -0.6180
#     7         O                 -1.3453   -1.2462   -0.8684
#     8         O                 -1.4046    1.4708    0.0167
#     9         O                 -0.5716   -0.5263    1.6651
#    10         C                  1.8529    1.0175    0.0716
#    11         C                 -1.5193   -1.0359   -2.2416
#    12         C                 -2.7764    1.5044    0.2897
#    13         C                 -0.0136   -1.7640    2.0001
#    14         C                  2.1985    2.3297   -0.6413
#    15         C                 -2.2972   -2.2169   -2.8050
#    16         C                 -3.2205    2.9603    0.3151
#    17         C                  1.2114   -1.5689    2.8841
#    18         H                  4.1028    0.8832    1.5483
# ...
         # Example 2 (used):
#   ATOM   CHEMICAL          X               Y               Z
#  NUMBER    SYMBOL      (ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)
# 
#     1       O          4.79280259  *  -0.84610232  *   0.36409474  *
#     2       O          5.89768035  *  -0.31706418  *   0.00917035  *
#     3       C          3.86164836  *   0.06535206  *   0.86290800  *
#     4       O          2.91352871  *   0.05485130  *  -0.07194851  *
#     5      Si         -0.61245484  *  -0.02707117  *   0.04871188  *
#     6       O          0.91999240  *   0.28181302  *  -0.61800545  *
#     7       O         -1.34526429  *  -1.24617340  *  -0.86844046  *
#     8       O         -1.40457125  *   1.47080489  *   0.01671181  *
#     9       O         -0.57162101  *  -0.52628027  *   1.66508989  *
#    10       C          1.85290140  *   1.01752620  *   0.07159039  *
#    11       C         -1.51932072  *  -1.03592573  *  -2.24160046  *
#    12       C         -2.77644395  *   1.50443941  *   0.28973441  *
#    13       C         -0.01360776  *  -1.76397803  *   2.00010724  *
#    14       C          2.19854080  *   2.32966388  *  -0.64131311  *
#    15       C         -2.29721668  *  -2.21688022  *  -2.80495545  *
#    16       C         -3.22047132  *   2.96028967  *   0.31511890  *
#    17       C          1.21142471  *  -1.56886315  *   2.88414255  *
#    18       H          4.10284938  *   0.88318846  *   1.54829483  *
#    19       H          1.60266809  *   1.19314394  *   1.14931859  *
#    20       H         -2.06992519  *  -0.08909329  *  -2.41564011  *
#    21       H         -0.53396028  *  -0.94280520  *  -2.73816125  *
#    22       H         -2.99280631  *   1.01386560  *   1.25905636  *
#    23       H         -3.32412961  *   0.94305635  *  -0.49427315  *
#    24       H         -0.81149878  *  -2.30331548  *   2.54543351  *
#    25       H          0.24486568  *  -2.37041735  *   1.10943219  *
#    26       H          2.46163770  *   2.17667287  *  -1.69615441  *
#    27       H          1.34364456  *   3.01690600  *  -0.61108044  *
#    28       H          3.04795301  *   2.82487051  *  -0.15380555  *
#    29       H         -1.76804185  *  -3.16646015  *  -2.65234745  *
#    30       H         -3.28543199  *  -2.31880074  *  -2.33789659  *
#    31       H         -2.45109195  *  -2.09228197  *  -3.88420787  *
#    32       H         -3.02567427  *   3.46605770  *  -0.63952294  *
#    33       H         -4.29770055  *   3.02763638  *   0.51281387  *
#    34       H         -2.70317481  *   3.53302115  *   1.09570604  *
#    35       H          2.01935375  *  -1.03805729  *   2.35810565  *
#    36       H          1.60901654  *  -2.53904354  *   3.20705714  *
#    37       H          0.97814118  *  -0.98964976  *   3.78695207  *
        if (line.find("NUMBER    SYMBOL      (ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)") > -1 or line.find("NUMBER   SYMBOL      (ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)") > -1):



            self.updateprogress(inputfile, "Attributes", self.cupdate)
                    
            self.inputcoords = []
            self.inputatoms = []
            
            blankline = inputfile.next()
            
            atomcoords = []
            line = inputfile.next()
           # while line != blankline:
            while len(line.split()) > 0:
                broken = line.split()
                self.inputatoms.append(symbol2int(broken[1]))
                xc = float(broken[2])
                yc = float(broken[4])
                zc = float(broken[6])
                atomcoords.append([xc,yc,zc])
                line = inputfile.next()

            self.inputcoords.append(atomcoords)

            if not hasattr(self, "natom"):
                self.atomnos = numpy.array(self.inputatoms, 'i')
                self.natom = len(self.atomnos)

#read energy (in kcal/mol, converted to eV)
#       Example:           FINAL HEAT OF FORMATION =       -333.88606 KCAL =   -1396.97927 KJ
        if line[0:35] == '          FINAL HEAT OF FORMATION =':
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            self.scfenergies.append(utils.convertor(self.float(line.split()[5])/627.5095, "hartree", "eV")) #note conversion from kcal/mol to hartree

        #molecular mass parsing (units will be amu)
        #Example:          MOLECULAR WEIGHT        =
        if line[0:35] == '          MOLECULAR WEIGHT        =':
            self.molmass = self.float(line.split()[3])
        
	  #rotational constants (converted to GHZ)
        #Example:

#          ROTATIONAL CONSTANTS IN CM(-1)
#
#          A =    0.01757641   B =    0.00739763   C =    0.00712013
        #could also read in moment of inertia, but this should just differ by a constant: rot cons= h/(8*Pi^2*I)
        #note that the last occurence of this in the thermochemistry section has reduced precision, so we will want to use the 2nd to last instance
        if line[0:40] == '          ROTATIONAL CONSTANTS IN CM(-1)':
	    blankline = inputfile.next();
            rotinfo=inputfile.next();
            if not hasattr(self, "rotcons"):
                self.rotcons = []
            broken = rotinfo.split()
            sol = 29.9792458 #speed of light in vacuum in 10^9 cm/s, cf. http://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=universal_in!
            a = float(broken[2])*sol 
            b = float(broken[5])*sol
            c = float(broken[8])*sol
            self.rotcons.append([a, b, c]) 

        # Start of the IR/Raman frequency section.
#Example:
# VIBRATION    1    1A       ATOM PAIR        ENERGY CONTRIBUTION    RADIAL
# FREQ.        15.08        C 12 --  C 16           +7.9% (999.0%)     0.0%
# T-DIPOLE    0.2028        C 16 --  H 34           +5.8% (999.0%)    28.0%
# TRAVEL      0.0240        C 16 --  H 32           +5.6% (999.0%)    35.0%
# RED. MASS   1.7712        O  1 --  O  4           +5.2% (999.0%)     0.4%
# EFF. MASS7752.8338
#
# VIBRATION    2    2A       ATOM PAIR        ENERGY CONTRIBUTION    RADIAL
# FREQ.        42.22        C 11 --  C 15           +9.0% (985.8%)     0.0%
# T-DIPOLE    0.1675        C 15 --  H 31           +6.6% (843.6%)     3.3%
# TRAVEL      0.0359        C 15 --  H 29           +6.0% (802.8%)    24.5%
# RED. MASS   1.7417        C 13 --  C 17           +5.8% (792.7%)     0.0%
# EFF. MASS1242.2114
        if line[1:10] == 'VIBRATION':
	    line = inputfile.next()
            self.updateprogress(inputfile, "Frequency Information", self.fupdate)
      
            if not hasattr(self, 'vibfreqs'):
                self.vibfreqs = []
            freq = self.float(line.split()[1])
            #self.vibfreqs.extend(freqs)
            self.vibfreqs.append(freq)


if __name__ == "__main__":
    import doctest, mopacparser
    doctest.testmod(mopacparser, verbose=False)
