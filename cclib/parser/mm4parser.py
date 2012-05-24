"""
gmagoon 05/03/10: new class for MM4 parsing, based on mopacparser.py, which, in turn, is based on gaussianparser.py from cclib, described below:
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 814 $"


#import re

import numpy
import math
import utils
import logfileparser


def symbol2int(symbol):
    t = utils.PeriodicTable()
    return t.number[symbol]

class MM4(logfileparser.Logfile):
    """An MM4 output file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(MM4, self).__init__(logname="MM4", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "MM4 log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'MM4("%s")' % (self.filename)

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""
        
        # Number of atoms.
        # Example:          THE COORDINATES OF    20 ATOMS ARE READ IN.
        if line[0:28] == '          THE COORDINATES OF':

            self.updateprogress(inputfile, "Attributes", self.fupdate)       
            natom = int(line.split()[-5]) #fifth to last component should be number of atoms
            if hasattr(self, "natom"):
                assert self.natom == natom
            else:
                self.natom = natom
        
        # Extract the atomic numbers and coordinates from the optimized (final) geometry
        
        # Example:
#	      FINAL ATOMIC COORDINATE
#           ATOM          X           Y           Z      TYPE
#         C(    1)    -3.21470    -0.22058     0.00000   (  1)
#         H(    2)    -3.30991    -0.87175     0.89724   (  5)
#         H(    3)    -3.30991    -0.87174    -0.89724   (  5)
#         H(    4)    -4.08456     0.47380     0.00000   (  5)
#         C(    5)    -1.88672     0.54893     0.00000   (  1)
#         H(    6)    -1.84759     1.21197    -0.89488   (  5)
#         H(    7)    -1.84759     1.21197     0.89488   (  5)
#         C(    8)    -0.66560    -0.38447     0.00000   (  1)
#         H(    9)    -0.70910    -1.04707    -0.89471   (  5)
#         H(   10)    -0.70910    -1.04707     0.89471   (  5)
#         C(   11)     0.66560     0.38447     0.00000   (  1)
#         H(   12)     0.70910     1.04707     0.89471   (  5)
#         H(   13)     0.70910     1.04707    -0.89471   (  5)
#         C(   14)     1.88672    -0.54893     0.00000   (  1)
#         H(   15)     1.84759    -1.21197    -0.89488   (  5)
#         H(   16)     1.84759    -1.21197     0.89488   (  5)
#         C(   17)     3.21470     0.22058     0.00000   (  1)
#         H(   18)     3.30991     0.87174     0.89724   (  5)
#         H(   19)     4.08456    -0.47380     0.00000   (  5)
#         H(   20)     3.30991     0.87175    -0.89724   (  5)

        if line[0:29] == '      FINAL ATOMIC COORDINATE':


            self.updateprogress(inputfile, "Attributes", self.cupdate)
                    
            self.inputcoords = []
            self.inputatoms = []
            
            headerline = inputfile.next()
            
            atomcoords = []
            line = inputfile.next()
            while len(line.split()) > 0:
                broken = line.split()
                self.inputatoms.append(symbol2int(line[0:10].strip()))
                xc = float(line[17:29])
                yc = float(line[29:41])
                zc = float(line[41:53])
                atomcoords.append([xc,yc,zc])
                line = inputfile.next()

            self.inputcoords.append(atomcoords)

	    if not hasattr(self, "atomnos"):
		self.atomnos = numpy.array(self.inputatoms, 'i')
            if not hasattr(self, "natom"):
                self.natom = len(self.atomnos)


#read energy (in kcal/mol, converted to eV)
#       Example:     HEAT OF FORMATION (HFN) AT  298.2 K       =       -42.51 KCAL/MOLE
        if line[0:31] == '     HEAT OF FORMATION (HFN) AT':
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            self.scfenergies.append(utils.convertor(self.float(line.split()[-2])/627.5095, "hartree", "eV")) #note conversion from kcal/mol to hartree

        #molecular mass parsing (units will be amu); note that this can occur multiple times in the file, but all values should be the same
        #Example:               FORMULA WEIGHT   :     86.112
        if line[0:33] == '               FORMULA WEIGHT   :':
            self.updateprogress(inputfile, "Attributes", self.fupdate)
	    molmass = self.float(line.split()[-1])
	    if hasattr(self, "molmass"):
                assert self.molmass == molmass #check that subsequent occurences match the original value
            else:
                self.molmass = molmass
        
	  #rotational constants (converted to GHZ)
        #Example:
#         THE MOMENTS OF INERTIA CALCULATED FROM R(g), R(z) VALUES
#                  (also from R(e), R(alpha), R(s) VALUES)
#
#         Note: (1) All calculations are based on principle isotopes.
#               (2) R(z) values include harmonic vibration (Coriolis)
#                   contribution indicated in parentheses.
#
#
#   (1)  UNIT = 10**(-39) GM*CM**2
#
#                    IX                   IY                   IZ
#
#   R(e)         5.7724              73.4297              76.0735
#   R(z)         5.7221(-0.0518)     74.0311(-0.0285)     76.7102(-0.0064)
#
#   (2)  UNIT = AU A**2
#
#                    IX                   IY                   IZ
#
#   R(e)        34.7661             442.2527             458.1757
#   R(z)        34.4633(-0.3117)    445.8746(-0.1714)    462.0104(-0.0385)
        #moments of inertia converted into rotational constants via rot cons= h/(8*Pi^2*I)
	#we will use the equilibrium values (R(e)) in units of 10**-39 GM*CM**2 (these units are less precise (fewer digits) than AU A**2 units but it is simpler as it doesn't require use of Avogadro's number
        #***even R(e) may include temperature dependent effects, though, and maybe the one I actually want is r(mm4) (not reported)
	if line[0:33] == '   (1)  UNIT = 10**(-39) GM*CM**2':
	    dummyline = inputfile.next();
	    dummyline = inputfile.next();
	    dummyline = inputfile.next();
            rotinfo=inputfile.next();
            if not hasattr(self, "rotcons"):
                self.rotcons = []
            broken = rotinfo.split()
	    h = 6.62606896E3 #Planck constant in 10^-37 J-s = 10^-37 kg m^2/s cf. http://physics.nist.gov/cgi-bin/cuu/Value?h#mid
            a = h/(8*math.pi*math.pi*float(broken[1]))
            b = h/(8*math.pi*math.pi*float(broken[2]))
            c = h/(8*math.pi*math.pi*float(broken[3]))
            self.rotcons.append([a, b, c]) 

        # Start of the IR/Raman frequency section.
#Example:
#0       FUNDAMENTAL NORMAL VIBRATIONAL FREQUENCIES
#                ( THEORETICALLY  54 VALUES )
#
#             Frequency :  in 1/cm
#             A(i)      :  IR intensity (vs,s,m,w,vw,-- or in 10**6 cm/mole)
#             A(i) = -- :  IR inactive
#
#
#             no       Frequency   Symmetry      A(i)
#
#             1.          2969.6     (Bu  )        s
#             2.          2969.6     (Bu  )        w
#             3.          2967.6     (Bu  )        w
#             4.          2967.6     (Bu  )        s
#             5.          2931.2     (Au  )       vs
#             6.          2927.8     (Bg  )       --
#             7.          2924.9     (Au  )        m
#             8.          2923.6     (Bg  )       --
#             9.          2885.8     (Ag  )       --
#            10.          2883.9     (Bu  )        w
#            11.          2879.8     (Ag  )       --
#            12.          2874.6     (Bu  )        w
#            13.          2869.6     (Ag  )       --
#            14.          2869.2     (Bu  )        s
#            15.          1554.4     (Ag  )       --
#            16.          1494.3     (Bu  )        w
#            17.          1449.7     (Bg  )       --
#            18.          1449.5     (Au  )        w
#            19.          1444.8     (Ag  )       --
#            20.          1438.5     (Bu  )        w
#            21.          1421.5     (Ag  )       --
#            22.          1419.3     (Ag  )       --
#            23.          1416.5     (Bu  )        w
#            24.          1398.8     (Bu  )        w
#            25.          1383.9     (Ag  )       --
#            26.          1363.7     (Bu  )        m
#            27.          1346.3     (Ag  )       --
#            28.          1300.2     (Au  )       vw
#            29.          1298.7     (Bg  )       --
#            30.          1283.4     (Bu  )        m
#            31.          1267.4     (Bg  )       --
#            32.          1209.6     (Au  )        w
#            33.          1132.2     (Bg  )       --
#            34.          1094.4     (Ag  )       --
#            35.          1063.4     (Bu  )        w
#            36.          1017.8     (Bu  )        w
#            37.          1011.6     (Ag  )       --
#            38.          1004.2     (Au  )        w
#            39.           990.2     (Ag  )       --
#            40.           901.8     (Ag  )       --
#            41.           898.4     (Bg  )       --
#            42.           875.9     (Bu  )        w
#	     43.           795.4     (Au  )        w
#            44.           725.0     (Bg  )       --
#            45.           699.6     (Au  )        w
#            46.           453.4     (Bu  )        w
#            47.           352.1     (Ag  )       --
#            48.           291.1     (Ag  )       --
#            49.           235.9     (Au  )       vw
#            50.           225.2     (Bg  )       --
#            51.           151.6     (Bg  )       --
#            52.           147.7     (Bu  )        w
#            53.           108.0     (Au  )       vw
#            54.            77.1     (Au  )       vw
#            55.     (       0.0)    (t/r )
#            56.     (       0.0)    (t/r )
#            57.     (       0.0)    (t/r )
#            58.     (       0.0)    (t/r )
#            59.     (       0.0)    (t/r )
#            60.     (       0.0)    (t/r )

        if line[0:52] == '             no       Frequency   Symmetry      A(i)':
	    blankline = inputfile.next()
            self.updateprogress(inputfile, "Frequency Information", self.fupdate)
      
            if not hasattr(self, 'vibfreqs'):
                self.vibfreqs = []
	    line = inputfile.next()
	    while(line[15:31].find('(') < 0):#terminate once we reach zero frequencies (which include parentheses)
		    freq = self.float(line[15:31])
		    self.vibfreqs.append(freq)
		    line = inputfile.next()
	#parsing of final steric energy in eV (for purposes of providing a baseline for possible subsequent hindered rotor calculations)
	#example line:"    FINAL STERIC ENERGY IS                0.8063 KCAL/MOL."
	if line[6:28] == 'FINAL STERIC ENERGY IS':
            stericenergy = utils.convertor(self.float(line.split()[4])/627.5095, "hartree", "eV") #note conversion from kcal/mol to hartree
	    if hasattr(self, "stericenergy"):
                assert self.stericenergy == stericenergy #check that subsequent occurences match the original value
            else:
                self.stericenergy = stericenergy


if __name__ == "__main__":
    import doctest, mm4parser
    doctest.testmod(mm4parser, verbose=False)
