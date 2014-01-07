"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 861 $"


import re

import numpy

import logfileparser
import utils


class GAMESSUK(logfileparser.Logfile):
    """A GAMESS UK log file"""
    SCFRMS, SCFMAX, SCFENERGY = range(3) # Used to index self.scftargets[]
    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(GAMESSUK, self).__init__(logname="GAMESSUK", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "GAMESS UK log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'GAMESSUK("%s")' % (self.filename)
    
    def normalisesym(self, label):
        """Use standard symmetry labels instead of GAMESS UK labels.

        >>> t = GAMESSUK("dummyfile.txt")
        >>> labels = ['a', 'a1', 'ag', "a'", 'a"', "a''", "a1''", 'a1"']
        >>> labels.extend(["e1+", "e1-"])
        >>> answer = [t.normalisesym(x) for x in labels]
        >>> answer
        ['A', 'A1', 'Ag', "A'", 'A"', 'A"', 'A1"', 'A1"', 'E1', 'E1']
        """
        label = label.replace("''", '"').replace("+", "").replace("-", "")
        ans = label[0].upper() + label[1:]
        
        return ans

    def before_parsing(self):

        # This will be used to detect the first set of "nuclear coordinates" in
        # a geometry-optimization
        self.firstnuccoords = True

        # used for determining whether to add a second mosyms, etc.
        self.betamosyms = self.betamoenergies = self.betamocoeffs = False

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        if line[1:22] == "total number of atoms":
            if not hasattr(self, "natom"):
                self.natom = int(line.split()[-1])

        if line[3:44] == "convergence threshold in optimization run":
            # Assuming that this is only found in the case of OPTXYZ
            # (i.e. an optimization in Cartesian coordinates)
            self.geotargets = [float(line.split()[-2])]

        if line[32:61] == "largest component of gradient":
            # This is the geotarget in the case of OPTXYZ
            if not hasattr(self, "geovalues"):
                self.geovalues = []
            self.geovalues.append([float(line.split()[4])])

        if line[37:49] == "convergence?":
            # Get the geovalues and geotargets for OPTIMIZE
            if not hasattr(self, "geovalues"):
                self.geovalues = []
                self.geotargets = []
            geotargets = []
            geovalues = []
            for i in range(4):
                temp = line.split()
                geovalues.append(float(temp[2]))
                if not self.geotargets:
                    geotargets.append(float(temp[-2]))
                line = inputfile.next()
            self.geovalues.append(geovalues)
            if not self.geotargets:
                self.geotargets = geotargets
        
        if line[40:58] == "molecular geometry":
            # Only one set of atomcoords is taken from this section
            # For geo-opts, more coordinates are taken from the "nuclear coordinates"
            if not hasattr(self, "atomcoords"):
                self.atomcoords = []
            self.atomnos = []
            
            stop = " "*9 + "*"*79
            line = inputfile.next()
            while not line.startswith(stop):
                line = inputfile.next()
            line = inputfile.next()
            while not line.startswith(stop):
                line = inputfile.next()
            empty = inputfile.next()

            atomcoords = []
            empty = inputfile.next()
            while not empty.startswith(stop):
                line = inputfile.next().split() # the coordinate data
                atomcoords.append(map(float,line[3:6]))
                self.atomnos.append(int(round(float(line[2]))))
                while line!=empty:
                    line = inputfile.next()
                # at this point, line is an empty line, right after
                # 1 or more lines containing basis set information
                empty = inputfile.next()
                # empty is either a row of asterisks or the empty line
                # before the row of coordinate data
            
            self.atomcoords.append(atomcoords)
            self.atomnos = numpy.array(self.atomnos, "i")

        if line[40:59] == "nuclear coordinates":
            # We need not remember the first geometry in the geo-opt as this will
            # be recorded already, in the "molecular geometry" section
            # (note: single-point calculations have no "nuclear coordinates" only
            # "molecular geometry")
            if self.firstnuccoords:
                self.firstnuccoords = False
                return
                # This was continue (in loop) before parser refactoring.
                # continue
            if not hasattr(self, "atomcoords"):
                self.atomcoords = []
                self.atomnos = []
                
            asterisk = inputfile.next()
            blank = inputfile.next()
            colmname = inputfile.next()
            equals = inputfile.next()

            atomcoords = []
            atomnos = []
            line = inputfile.next()
            while line != equals:
                temp = line.strip().split()
                atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom") for x in temp[0:3]])
                if not hasattr(self, "atomnos") or len(self.atomnos) == 0:
                    atomnos.append(int(float(temp[3])))
                    
                line = inputfile.next()

            self.atomcoords.append(atomcoords)
            if not hasattr(self, "atomnos") or len(self.atomnos) == 0:
                self.atomnos = atomnos

        if line[1:32] == "total number of basis functions":
            self.nbasis = int(line.split()[-1])
            while line.find("charge of molecule")<0:
                line = inputfile.next()
            self.charge = int(line.split()[-1])
            self.mult = int(inputfile.next().split()[-1])

            alpha = int(inputfile.next().split()[-1])-1
            beta = int(inputfile.next().split()[-1])-1
            if self.mult==1:
                self.homos = numpy.array([alpha], "i")
            else:
                self.homos = numpy.array([alpha,beta], "i")

        if line[37:69] == "s-matrix over gaussian basis set":
            self.aooverlaps = numpy.zeros((self.nbasis, self.nbasis), "d")

            minus = inputfile.next()
            blank = inputfile.next()
            i = 0
            while i < self.nbasis:
                blank = inputfile.next()
                blank = inputfile.next()
                header = inputfile.next()
                blank = inputfile.next()
                blank = inputfile.next()

                for j in range(self.nbasis):
                    temp = map(float, inputfile.next().split()[1:])
                    self.aooverlaps[j,(0+i):(len(temp)+i)] = temp
                    
                i += len(temp)

        if line[18:43] == 'EFFECTIVE CORE POTENTIALS':
            self.coreelectrons = numpy.zeros(self.natom, 'i')
            asterisk = inputfile.next()
            line = inputfile.next()
            while line[15:46]!="*"*31:
                if line.find("for atoms ...")>=0:
                    atomindex = []
                    line = inputfile.next()
                    while line.find("core charge")<0:
                        broken = line.split()
                        atomindex.extend([int(x.split("-")[0]) for x in broken])
                        line = inputfile.next()
                    charge = float(line.split()[4])
                    for idx in atomindex:
                        self.coreelectrons[idx-1] = self.atomnos[idx-1] - charge
                line = inputfile.next()
                            
        if line[3:27] == "Wavefunction convergence":
            self.scftarget = float(line.split()[-2])
            self.scftargets = []

        if line[11:22] == "normal mode":
            if not hasattr(self, "vibfreqs"):
                self.vibfreqs = []
                self.vibirs = []
            
            units = inputfile.next()
            xyz = inputfile.next()
            equals = inputfile.next()
            line = inputfile.next()
            while line!=equals:
                temp = line.split()
                self.vibfreqs.append(float(temp[1]))
                self.vibirs.append(float(temp[-2]))
                line = inputfile.next()
            # Use the length of the vibdisps to figure out
            # how many rotations and translations to remove
            self.vibfreqs = self.vibfreqs[-len(self.vibdisps):]
            self.vibirs = self.vibirs[-len(self.vibdisps):]

        if line[44:73] == "normalised normal coordinates":
            self.vibdisps = []
            equals = inputfile.next()
            blank = inputfile.next()
            blank = inputfile.next()
            freqnum = inputfile.next()
            while freqnum.find("=")<0:
                blank = inputfile.next()
                equals = inputfile.next()
                freqs = inputfile.next()
                equals = inputfile.next()
                blank = inputfile.next()
                header = inputfile.next()
                equals = inputfile.next()
                p = [ [] for x in range(9) ]
                for i in range(len(self.atomnos)):
                    brokenx = map(float, inputfile.next()[25:].split())
                    brokeny = map(float, inputfile.next()[25:].split())            
                    brokenz = map(float, inputfile.next()[25:].split())
                    for j,x in enumerate(zip(brokenx, brokeny, brokenz)):
                        p[j].append(x)
                self.vibdisps.extend(p)
        
                blank = inputfile.next()
                blank = inputfile.next()
                freqnum = inputfile.next()                    

        if line[26:36] == "raman data":
            self.vibramans = []

            stars = inputfile.next()
            blank = inputfile.next()
            header = inputfile.next()

            blank = inputfile.next()
            line = inputfile.next()
            while line[1]!="*":
                self.vibramans.append(float(line.split()[3]))
                blank = inputfile.next()
                line = inputfile.next()
            # Use the length of the vibdisps to figure out
            # how many rotations and translations to remove
            self.vibramans = self.vibramans[-len(self.vibdisps):]
                        
        if line[3:11] == "SCF TYPE":
            self.scftype = line.split()[-2]
            assert self.scftype in ['rhf', 'uhf', 'gvb'], "%s not one of 'rhf', 'uhf' or 'gvb'" % self.scftype

        if line[15:31] == "convergence data":
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []
            self.scftargets.append([self.scftarget]) # Assuming it does not change over time
            while line[1:10] != "="*9:
                line = inputfile.next()
            line = inputfile.next()
            tester = line.find("tester") # Can be in a different place depending
            assert tester>=0
            while line[1:10] != "="*9: # May be two or three lines (unres)
                line = inputfile.next()
            
            scfvalues = []
            line = inputfile.next()
            while line.strip():
                if line[2:6]!="****":
            # e.g. **** recalulation of fock matrix on iteration  4 (examples/chap12/pyridine.out)
                    scfvalues.append([float(line[tester-5:tester+6])])
                line = inputfile.next()
            self.scfvalues.append(scfvalues)   

        if line[10:22] == "total energy" and len(line.split()) == 3:
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            scfenergy = utils.convertor(float(line.split()[-1]), "hartree", "eV")
            self.scfenergies.append(scfenergy)
        
        # Total energies after Moller-Plesset corrections
        # Second order correction is always first, so its first occurance
        #   triggers creation of mpenergies (list of lists of energies)
        # Further corrections are appended as found
        # Note: GAMESS-UK sometimes prints only the corrections,
        #   so they must be added to the last value of scfenergies
        if line[10:32] == "mp2 correlation energy" or \
           line[10:42] == "second order perturbation energy":
            if not hasattr(self, "mpenergies"):
                self.mpenergies = []
            self.mpenergies.append([])
            self.mp2correction = self.float(line.split()[-1])
            self.mp2energy = self.scfenergies[-1] + self.mp2correction
            self.mpenergies[-1].append(utils.convertor(self.mp2energy, "hartree", "eV"))
        if line[10:41] == "third order perturbation energy":
            self.mp3correction = self.float(line.split()[-1])
            self.mp3energy = self.mp2energy + self.mp3correction
            self.mpenergies[-1].append(utils.convertor(self.mp3energy, "hartree", "eV"))

        if line[40:59] == "molecular basis set":
            self.gbasis = []
            line = inputfile.next()
            while line.find("contraction coefficients")<0:
                line = inputfile.next()
            equals = inputfile.next()
            blank = inputfile.next()
            atomname = inputfile.next()
            basisregexp = re.compile("\d*(\D+)") # Get everything after any digits
            shellcounter = 1
            while line!=equals:
                gbasis = [] # Stores basis sets on one atom
                blank = inputfile.next()
                blank = inputfile.next()
                line = inputfile.next()
                shellno = int(line.split()[0])
                shellgap = shellno - shellcounter
                shellsize = 0
                while len(line.split())!=1 and line!=equals:
                    if line.split():
                        shellsize += 1
                    coeff = {}
                    # coefficients and symmetries for a block of rows
                    while line.strip() and line!=equals:
                        temp = line.strip().split()
                    # temp[1] may be either like (a) "1s" and "1sp", or (b) "s" and "sp"
                    # See GAMESS-UK 7.0 distribution/examples/chap12/pyridine2_21m10r.out
                    # for an example of the latter
                        sym = basisregexp.match(temp[1]).groups()[0]
                        assert sym in ['s', 'p', 'd', 'f', 'sp'], "'%s' not a recognized symmetry" % sym
                        if sym == "sp":
                            coeff.setdefault("S", []).append( (float(temp[3]), float(temp[6])) )
                            coeff.setdefault("P", []).append( (float(temp[3]), float(temp[10])) )
                        else:
                            coeff.setdefault(sym.upper(), []).append( (float(temp[3]), float(temp[6])) )
                        line = inputfile.next()
                    # either a blank or a continuation of the block
                    if coeff:
                        if sym == "sp":
                            gbasis.append( ('S', coeff['S']))
                            gbasis.append( ('P', coeff['P']))
                        else:
                            gbasis.append( (sym.upper(), coeff[sym.upper()]))
                    if line==equals:
                        continue
                    line = inputfile.next()
                    # either the start of the next block or the start of a new atom or
                    # the end of the basis function section (signified by a line of equals)
                numtoadd = 1 + (shellgap / shellsize)
                shellcounter = shellno + shellsize
                for x in range(numtoadd):
                    self.gbasis.append(gbasis)

        if line[50:70] == "----- beta set -----":
            self.betamosyms = True
            self.betamoenergies = True
            self.betamocoeffs = True
            # betamosyms will be turned off in the next
            # SYMMETRY ASSIGNMENT section
                
        if line[31:50] == "SYMMETRY ASSIGNMENT":
            if not hasattr(self, "mosyms"):
                self.mosyms = []

            multiple = {'a':1, 'b':1, 'e':2, 't':3, 'g':4, 'h':5}
            
            equals = inputfile.next()
            line = inputfile.next()
            while line != equals: # There may be one or two lines of title (compare mg10.out and duhf_1.out)
                line = inputfile.next()

            mosyms = []
            line = inputfile.next()
            while line != equals:
                temp = line[25:30].strip()
                if temp[-1]=='?':
                    # e.g. e? or t? or g? (see example/chap12/na7mg_uhf.out)
                    # for two As, an A and an E, and two Es of the same energy respectively.
                    t = line[91:].strip().split()
                    for i in range(1,len(t),2):
                        for j in range(multiple[t[i][0]]): # add twice for 'e', etc.
                            mosyms.append(self.normalisesym(t[i]))
                else:
                    for j in range(multiple[temp[0]]):
                        mosyms.append(self.normalisesym(temp)) # add twice for 'e', etc.
                line = inputfile.next()
            assert len(mosyms) == self.nmo, "mosyms: %d but nmo: %d" % (len(mosyms), self.nmo)
            if self.betamosyms:
                # Only append if beta (otherwise with IPRINT SCF
                # it will add mosyms for every step of a geo opt)
                self.mosyms.append(mosyms)
                self.betamosyms = False
            elif self.scftype=='gvb':
                # gvb has alpha and beta orbitals but they are identical
                self.mosysms = [mosyms, mosyms]
            else:
                self.mosyms = [mosyms]

        if line[50:62] == "eigenvectors":
        # Mocoeffs...can get evalues from here too
        # (only if using FORMAT HIGH though will they all be present)                
            if not hasattr(self, "mocoeffs"):
                self.aonames = []
                aonames = []
            minus = inputfile.next()

            mocoeffs = numpy.zeros( (self.nmo, self.nbasis), "d")
            readatombasis = False
            if not hasattr(self, "atombasis"):
                self.atombasis = []
                for i in range(self.natom):
                    self.atombasis.append([])
                readatombasis = True

            blank = inputfile.next()
            blank = inputfile.next()
            evalues = inputfile.next()

            p = re.compile(r"\d+\s+(\d+)\s*(\w+) (\w+)")
            oldatomname = "DUMMY VALUE"

            mo = 0
            while mo < self.nmo:
                blank = inputfile.next()
                blank = inputfile.next()
                nums = inputfile.next()
                blank = inputfile.next()
                blank = inputfile.next()
                for basis in range(self.nbasis):
                    line = inputfile.next()
                    # Fill atombasis only first time around.
                    if readatombasis:
                        orbno = int(line[1:5])-1
                        atomno = int(line[6:9])-1
                        self.atombasis[atomno].append(orbno)
                    if not self.aonames:
                        pg = p.match(line[:18].strip()).groups()
                        atomname = "%s%s%s" % (pg[1][0].upper(), pg[1][1:], pg[0])
                        if atomname!=oldatomname:
                            aonum = 1
                        oldatomname = atomname
                        name = "%s_%d%s" % (atomname, aonum, pg[2].upper())
                        if name in aonames:
                            aonum += 1
                        name = "%s_%d%s" % (atomname, aonum, pg[2].upper())
                        aonames.append(name) 
                    temp = map(float, line[19:].split())
                    mocoeffs[mo:(mo+len(temp)), basis] = temp
                # Fill atombasis only first time around.
                readatombasis = False
                if not self.aonames:
                    self.aonames = aonames

                line = inputfile.next() # blank line
                while line==blank:
                    line = inputfile.next()
                evalues = line
                if evalues[:17].strip(): # i.e. if these aren't evalues
                    break # Not all the MOs are present
                mo += len(temp)
            mocoeffs = mocoeffs[0:(mo+len(temp)), :] # In case some aren't present
            if self.betamocoeffs:
                self.mocoeffs.append(mocoeffs)
            else:
                self.mocoeffs = [mocoeffs]

        if line[7:12] == "irrep":
            ########## eigenvalues ###########
            # This section appears once at the start of a geo-opt and once at the end
            # unless IPRINT SCF is used (when it appears at every step in addition)
            if not hasattr(self, "moenergies"):
                self.moenergies = []

            equals = inputfile.next()
            while equals[1:5] != "====": # May be one or two lines of title (compare duhf_1.out and mg10.out)
                equals = inputfile.next()

            moenergies = []
            line = inputfile.next()
            if not line.strip(): # May be a blank line here (compare duhf_1.out and mg10.out)
                line = inputfile.next()

            while line.strip() and line != equals: # May end with a blank or equals
                temp = line.strip().split()
                moenergies.append(utils.convertor(float(temp[2]), "hartree", "eV"))
                line = inputfile.next()
            self.nmo = len(moenergies)
            if self.betamoenergies:
                self.moenergies.append(moenergies)
                self.betamoenergies = False
            elif self.scftype=='gvb':
                self.moenergies = [moenergies, moenergies]
            else:
                self.moenergies = [moenergies]
                
             
if __name__ == "__main__":
    import doctest
    doctest.testmod()
