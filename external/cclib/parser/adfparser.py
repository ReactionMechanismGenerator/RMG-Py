"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 861 $"


import numpy

import logfileparser
import utils


class ADF(logfileparser.Logfile):
    """An ADF log file"""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(ADF, self).__init__(logname="ADF", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "ADF log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'ADF("%s")' % (self.filename)

    def normalisesym(self, label):
        """Use standard symmetry labels instead of ADF labels.

        To normalise:
        (1) any periods are removed (except in the case of greek letters)
        (2) XXX is replaced by X, and a " added.
        (3) XX is replaced by X, and a ' added.
        (4) The greek letters Sigma, Pi, Delta and Phi are replaced by
            their lowercase equivalent.

        >>> sym = ADF("dummyfile").normalisesym
        >>> labels = ['A','s','A1','A1.g','Sigma','Pi','Delta','Phi','Sigma.g','A.g','AA','AAA','EE1','EEE1']
        >>> map(sym,labels)
        ['A', 's', 'A1', 'A1g', 'sigma', 'pi', 'delta', 'phi', 'sigma.g', 'Ag', "A'", 'A"', "E1'", 'E1"']
        """
        greeks = ['Sigma', 'Pi', 'Delta', 'Phi']
        for greek in greeks:
            if label.startswith(greek):
                return label.lower()

        ans = label.replace(".", "")
        if ans[1:3] == "''":
            temp = ans[0] + '"'
            ans = temp

        l = len(ans)
        if l > 1 and ans[0] == ans[1]: # Python only tests the second condition if the first is true
            if l > 2 and ans[1] == ans[2]:
                ans = ans.replace(ans[0]*3, ans[0]) + '"'
            else:
                ans = ans.replace(ans[0]*2, ans[0]) + "'"
        return ans

    def normalisedegenerates(self, label, num, ndict=None):
        """Generate a string used for matching degenerate orbital labels

        To normalise:
        (1) if label is E or T, return label:num
        (2) if label is P or D, look up in dict, and return answer
        """

        if not ndict:
          ndict = { 'P': {0:"P:x", 1:"P:y", 2:"P:z"},\
                    'D': {0:"D:z2", 1:"D:x2-y2", 2:"D:xy", 3:"D:xz", 4:"D:yz"}}

        if ndict.has_key(label):
            if ndict[label].has_key(num):
                return ndict[label][num]
            else:
                return "%s:%i"%(label,num+1)
        else:
            return "%s:%i"%(label,num+1)

    def before_parsing(self):

        # Used to avoid extracting the final geometry twice in a GeoOpt
        self.NOTFOUND, self.GETLAST, self.NOMORE = range(3)
        self.finalgeometry = self.NOTFOUND 

        # Used for calculating the scftarget (variables names taken from the ADF manual)
        self.accint = self.SCFconv = self.sconv2 = None

        # keep track of nosym and unrestricted case to parse Energies since it doens't have an all Irreps section
        self.nosymflag = False
        self.unrestrictedflag = False

        SCFCNV, SCFCNV2 = range(2) #used to index self.scftargets[]
        maxelem, norm = range(2) # used to index scf.values

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        if line.find("INPUT FILE") >= 0:
        #check to make sure we aren't parsing Create jobs
            while line:

                self.updateprogress(inputfile, "Unsupported Information", self.fupdate)

                if line.find("INPUT FILE") >=0 and hasattr(self,"scftargets"):
                #does this file contain multiple calculations?
                #if so, print a warning and skip to end of file
                    self.logger.warning("Skipping remaining calculations")
                    inputfile.seek(0,2)
                    break

                if line.find("INPUT FILE") >= 0:
                    line2 = inputfile.next()
                else:
                    line2 = None

                if line2 and len(line2) <= 2:
                #make sure that it's not blank like in the NiCO4 regression
                    line2 = inputfile.next()

                if line2 and (line2.find("Create") < 0 and line2.find("create") < 0):
                    break

                line = inputfile.next()

        if line[1:10] == "Symmetry:":
            info = line.split()
            if info[1] == "NOSYM":
                self.nosymflag = True

        # Use this to read the subspecies of irreducible representations.
        # It will be a list, with each element representing one irrep.
        if line.strip() == "Irreducible Representations, including subspecies":
            dashes = inputfile.next()
            self.irreps = []
            line = inputfile.next()
            while line.strip() != "":
                self.irreps.append(line.split())
                line = inputfile.next()

        if line[4:13] == 'Molecule:':
            info = line.split()
            if info[1] == 'UNrestricted':
                self.unrestrictedflag = True

        if line[1:6] == "ATOMS":
        # Find the number of atoms and their atomic numbers
        # Also extract the starting coordinates (for a GeoOpt anyway)
            self.updateprogress(inputfile, "Attributes", self.cupdate)

            self.atomnos = []
            self.atomcoords = []
            self.coreelectrons = []

            underline = inputfile.next()  #clear pointless lines
            label1 = inputfile.next()     # 
            label2 = inputfile.next()     #
            line = inputfile.next()
            atomcoords = []
            while len(line)>2: #ensure that we are reading no blank lines
                info = line.split()
                element = info[1].split('.')[0]
                self.atomnos.append(self.table.number[element])
                atomcoords.append(map(float, info[2:5]))
                self.coreelectrons.append(int(float(info[5]) - float(info[6])))
                line = inputfile.next()
            self.atomcoords.append(atomcoords)

            self.natom = len(self.atomnos)
            self.atomnos = numpy.array(self.atomnos, "i")

        if line[1:10] == "FRAGMENTS":
            header = inputfile.next()

            self.frags = []
            self.fragnames = []

            line = inputfile.next()
            while len(line) > 2: #ensure that we are reading no blank lines
                info = line.split()

                if len(info) == 7: #fragment name is listed here
                    self.fragnames.append("%s_%s"%(info[1],info[0]))
                    self.frags.append([])
                    self.frags[-1].append(int(info[2]) - 1)

                elif len(info) == 5: #add atoms into last fragment
                    self.frags[-1].append(int(info[0]) - 1)

                line = inputfile.next()

        # Extract charge
        if line[1:11] == "Net Charge":
            self.charge = int(line.split()[2])
            line = inputfile.next()
            if len(line.strip()):
                #  Spin polar: 1 (Spin_A minus Spin_B electrons)
                self.mult = int(line.split()[2]) + 1
                 # (Not sure about this for higher multiplicities)
            else:
                self.mult = 1

        if line[1:22] == "S C F   U P D A T E S":
        # find targets for SCF convergence

            if not hasattr(self,"scftargets"):
                self.scftargets = []

            #underline, blank, nr
            for i in range(3):
                inputfile.next()

            line = inputfile.next()
            self.SCFconv = float(line.split()[-1])
            line = inputfile.next()
            self.sconv2 = float(line.split()[-1])

        if line[1:11] == "CYCLE    1":

            self.updateprogress(inputfile, "QM convergence", self.fupdate)

            newlist = []
            line = inputfile.next()

            if not hasattr(self,"geovalues"):
                # This is the first SCF cycle
                self.scftargets.append([self.sconv2*10, self.sconv2])
            elif self.finalgeometry in [self.GETLAST, self.NOMORE]:
                # This is the final SCF cycle
                self.scftargets.append([self.SCFconv*10, self.SCFconv])
            else:
                # This is an intermediate SCF cycle
                oldscftst = self.scftargets[-1][1]
                grdmax = self.geovalues[-1][1]
                scftst = max(self.SCFconv, min(oldscftst, grdmax/30, 10**(-self.accint)))
                self.scftargets.append([scftst*10, scftst])

            while line.find("SCF CONVERGED") == -1 and line.find("SCF not fully converged, result acceptable") == -1 and line.find("SCF NOT CONVERGED") == -1:
                if line[4:12] == "SCF test":
                    if not hasattr(self, "scfvalues"):
                        self.scfvalues = []

                    info = line.split()
                    newlist.append([float(info[4]), abs(float(info[6]))])
                try:
                    line = inputfile.next()
                except StopIteration: #EOF reached?
                    self.logger.warning("SCF did not converge, so attributes may be missing")
                    break            

            if line.find("SCF not fully converged, result acceptable") > 0:
                self.logger.warning("SCF not fully converged, results acceptable")

            if line.find("SCF NOT CONVERGED") > 0:
                self.logger.warning("SCF did not converge! moenergies and mocoeffs are unreliable")

            if hasattr(self, "scfvalues"):
                self.scfvalues.append(newlist)

        # Parse SCF energy for SP calcs from bonding energy decomposition section.
        # It seems ADF does not print it earlier for SP calcualtions.
        # If it does (does it?), parse that instead.
        # Check that scfenergies does not exist, becuase gopt runs also print this,
        #   repeating the values in the last "Geometry Convergence Tests" section.
        if "Total Bonding Energy:" in line:
            if not hasattr(self, "scfenergies"):
                energy = utils.convertor(float(line.split()[3]), "hartree", "eV")
                self.scfenergies = [energy]            

        if line[51:65] == "Final Geometry":
            self.finalgeometry = self.GETLAST

        if line[1:24] == "Coordinates (Cartesian)" and self.finalgeometry in [self.NOTFOUND, self.GETLAST]:
            # Get the coordinates from each step of the GeoOpt
            if not hasattr(self, "atomcoords"):
                self.atomcoords = []
            equals = inputfile.next()
            blank = inputfile.next()
            title = inputfile.next()
            title = inputfile.next()
            hyphens = inputfile.next()

            atomcoords = []
            line = inputfile.next()
            while line != hyphens:
                atomcoords.append(map(float, line.split()[5:8]))
                line = inputfile.next()
            self.atomcoords.append(atomcoords)
            if self.finalgeometry == self.GETLAST: # Don't get any more coordinates
                self.finalgeometry = self.NOMORE

        if line[1:27] == 'Geometry Convergence Tests':
        # Extract Geometry convergence information
            if not hasattr(self, "geotargets"):
                self.geovalues = []
                self.geotargets = numpy.array([0.0, 0.0, 0.0, 0.0, 0.0], "d")
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            equals = inputfile.next()
            blank = inputfile.next()
            line = inputfile.next()
            temp = inputfile.next().strip().split()
            self.scfenergies.append(utils.convertor(float(temp[-1]), "hartree", "eV"))
            for i in range(6):
                line = inputfile.next()
            values = []
            for i in range(5):
                temp = inputfile.next().split()
                self.geotargets[i] = float(temp[-3])
                values.append(float(temp[-4]))
            self.geovalues.append(values)

        if line[1:27] == 'General Accuracy Parameter':
            # Need to know the accuracy of the integration grid to
            # calculate the scftarget...note that it changes with time
            self.accint = float(line.split()[-1])

        if line.find('Orbital Energies, per Irrep and Spin') > 0 and not hasattr(self, "mosyms") and self.nosymflag and not self.unrestrictedflag:
        #Extracting orbital symmetries and energies, homos for nosym case
        #Should only be for restricted case because there is a better text block for unrestricted and nosym

            self.mosyms = [[]]

            self.moenergies = [[]]

            underline = inputfile.next()
            header = inputfile.next()
            underline = inputfile.next()
            label = inputfile.next()
            line = inputfile.next()

            info = line.split()

            if not info[0] == '1':
                self.logger.warning("MO info up to #%s is missing" % info[0])

            #handle case where MO information up to a certain orbital are missing
            while int(info[0]) - 1 != len(self.moenergies[0]):
                self.moenergies[0].append(99999)
                self.mosyms[0].append('A')

            homoA = None

            while len(line) > 10:
                info = line.split()
                self.mosyms[0].append('A')
                self.moenergies[0].append(utils.convertor(float(info[2]), 'hartree', 'eV'))
                if info[1] == '0.000' and not hasattr(self, 'homos'):
                    self.homos = [len(self.moenergies[0]) - 2]
                line = inputfile.next()

            self.moenergies = [numpy.array(self.moenergies[0], "d")]
            self.homos = numpy.array(self.homos, "i")

        if line[1:29] == 'Orbital Energies, both Spins' and not hasattr(self, "mosyms") and self.nosymflag and self.unrestrictedflag:
        #Extracting orbital symmetries and energies, homos for nosym case
        #should only be here if unrestricted and nosym

            self.mosyms = [[], []]

            moenergies = [[], []]

            underline = inputfile.next()
            blank = inputfile.next()
            header = inputfile.next()
            underline = inputfile.next()
            line = inputfile.next()

            homoa = 0
            homob = None

            while len(line) > 5:
                info = line.split()
                if info[2] == 'A': 
                    self.mosyms[0].append('A')
                    moenergies[0].append(utils.convertor(float(info[4]), 'hartree', 'eV'))
                    if info[3] != '0.00':
                        homoa = len(moenergies[0]) - 1
                elif info[2] == 'B':
                    self.mosyms[1].append('A')
                    moenergies[1].append(utils.convertor(float(info[4]), 'hartree', 'eV'))
                    if info[3] != '0.00':
                        homob = len(moenergies[1]) - 1
                else:
                    print "Error reading line: %s" % line

                line = inputfile.next()

            self.moenergies = [numpy.array(x, "d") for x in moenergies]
            self.homos = numpy.array([homoa, homob], "i")


        if line[1:29] == 'Orbital Energies, all Irreps' and not hasattr(self, "mosyms"):
        #Extracting orbital symmetries and energies, homos
            self.mosyms = [[]]
            self.symlist = {}

            self.moenergies = [[]]

            underline = inputfile.next()
            blank = inputfile.next()
            header = inputfile.next()
            underline2 = inputfile.next()
            line = inputfile.next()

            homoa = None
            homob = None

            #multiple = {'E':2, 'T':3, 'P':3, 'D':5}
            # The above is set if there are no special irreps
            names = [irrep[0].split(':')[0] for irrep in self.irreps]
            counts = [len(irrep) for irrep in self.irreps]
            multiple = dict(zip(names, counts))
            irrepspecies = {}
            for n in range(len(names)):
                indices = range(counts[n])
                subspecies = self.irreps[n]
                irrepspecies[names[n]] = dict(zip(indices, subspecies))

            while line.strip():
                info = line.split()
                if len(info) == 5: #this is restricted
                    #count = multiple.get(info[0][0],1)
                    count = multiple.get(info[0],1)
                    for repeat in range(count): # i.e. add E's twice, T's thrice
                        self.mosyms[0].append(self.normalisesym(info[0]))
                        self.moenergies[0].append(utils.convertor(float(info[3]), 'hartree', 'eV'))

                        sym = info[0]
                        if count > 1: # add additional sym label
                            sym = self.normalisedegenerates(info[0],repeat,ndict=irrepspecies)

                        try:
                            self.symlist[sym][0].append(len(self.moenergies[0])-1)
                        except KeyError:
                            self.symlist[sym]=[[]]
                            self.symlist[sym][0].append(len(self.moenergies[0])-1)

                    if info[2] == '0.00' and not hasattr(self, 'homos'):
                        self.homos = [len(self.moenergies[0]) - (count + 1)] #count, because need to handle degenerate cases
                    line = inputfile.next()
                elif len(info) == 6: #this is unrestricted
                    if len(self.moenergies) < 2: #if we don't have space, create it
                        self.moenergies.append([])
                        self.mosyms.append([])
#                    count = multiple.get(info[0][0], 1)
                    count = multiple.get(info[0], 1)
                    if info[2] == 'A':
                        for repeat in range(count): # i.e. add E's twice, T's thrice
                            self.mosyms[0].append(self.normalisesym(info[0]))
                            self.moenergies[0].append(utils.convertor(float(info[4]), 'hartree', 'eV'))

                            sym = info[0]
                            if count > 1: #add additional sym label
                                sym = self.normalisedegenerates(info[0],repeat)

                            try:
                                self.symlist[sym][0].append(len(self.moenergies[0])-1)
                            except KeyError:
                                self.symlist[sym]=[[],[]]
                                self.symlist[sym][0].append(len(self.moenergies[0])-1)

                        if info[3] == '0.00' and homoa == None:
                            homoa = len(self.moenergies[0]) - (count + 1) #count because degenerate cases need to be handled

                    if info[2] == 'B':
                        for repeat in range(count): # i.e. add E's twice, T's thrice
                            self.mosyms[1].append(self.normalisesym(info[0]))
                            self.moenergies[1].append(utils.convertor(float(info[4]), 'hartree', 'eV'))

                            sym = info[0]
                            if count > 1: #add additional sym label
                                sym = self.normalisedegenerates(info[0],repeat)

                            try:
                                self.symlist[sym][1].append(len(self.moenergies[1])-1)
                            except KeyError:
                                self.symlist[sym]=[[],[]]
                                self.symlist[sym][1].append(len(self.moenergies[1])-1)

                        if info[3] == '0.00' and homob == None:
                            homob = len(self.moenergies[1]) - (count + 1)

                    line = inputfile.next()

                else: #different number of lines
                    print "Error", info

            if len(info) == 6: #still unrestricted, despite being out of loop
                self.homos = [homoa, homob]

            self.moenergies = [numpy.array(x, "d") for x in self.moenergies]
            self.homos = numpy.array(self.homos, "i")

        if line[1:28] == "Vibrations and Normal Modes":
            # Section on extracting vibdisps
            # Also contains vibfreqs, but these are extracted in the
            # following section (see below)
            self.vibdisps = []
            equals = inputfile.next()
            blank = inputfile.next()
            header = inputfile.next()
            header = inputfile.next()
            blank = inputfile.next()
            blank = inputfile.next()

            freqs = inputfile.next()
            while freqs.strip()!="":
                minus = inputfile.next()
                p = [ [], [], [] ]
                for i in range(len(self.atomnos)):
                    broken = map(float, inputfile.next().split()[1:])
                    for j in range(0, len(broken), 3):
                        p[j/3].append(broken[j:j+3])
                self.vibdisps.extend(p[:(len(broken)/3)])
                blank = inputfile.next()
                blank = inputfile.next()
                freqs = inputfile.next()
            self.vibdisps = numpy.array(self.vibdisps, "d")

        if line[1:24] == "List of All Frequencies":
        # Start of the IR/Raman frequency section
            self.updateprogress(inputfile, "Frequency information", self.fupdate)

        #                 self.vibsyms = [] # Need to look into this a bit more
            self.vibirs = []
            self.vibfreqs = []
            for i in range(8):
                line = inputfile.next()
            line = inputfile.next().strip()
            while line:
                temp = line.split()
                self.vibfreqs.append(float(temp[0]))                    
                self.vibirs.append(float(temp[2])) # or is it temp[1]?
                line = inputfile.next().strip()
            self.vibfreqs = numpy.array(self.vibfreqs, "d")
            self.vibirs = numpy.array(self.vibirs, "d")
            if hasattr(self, "vibramans"):
                self.vibramans = numpy.array(self.vibramans, "d")


        #******************************************************************************************************************8
        #delete this after new implementation using smat, eigvec print,eprint?
        if line[1:49] == "Total nr. of (C)SFOs (summation over all irreps)":
        # Extract the number of basis sets
            self.nbasis = int(line.split(":")[1].split()[0])

        # now that we're here, let's extract aonames

            self.fonames = []
            self.start_indeces = {}

            blank = inputfile.next()
            note = inputfile.next()
            symoffset = 0

            blank = inputfile.next() 
            blank = inputfile.next()
            if len(blank) > 2: #fix for ADF2006.01 as it has another note
                blank = inputfile.next()
                blank = inputfile.next()
            blank = inputfile.next()

            self.nosymreps = []
            while len(self.fonames) < self.nbasis:

                symline = inputfile.next()
                sym = symline.split()[1]
                line = inputfile.next()
                num = int(line.split(':')[1].split()[0])
                self.nosymreps.append(num)

                #read until line "--------..." is found
                while line.find('-----') < 0:
                    line = inputfile.next()

                line = inputfile.next() # the start of the first SFO

                while len(self.fonames) < symoffset + num:
                    info = line.split()

                    #index0 index1 occ2 energy3/4 fragname5 coeff6 orbnum7 orbname8 fragname9
                    if not sym in self.start_indeces.keys():
                    #have we already set the start index for this symmetry?
                        self.start_indeces[sym] = int(info[1])

                    orbname = info[8]
                    orbital = info[7] + orbname.replace(":", "")

                    fragname = info[5]
                    frag = fragname + info[9]

                    coeff = float(info[6])

                    line = inputfile.next()
                    while line.strip() and not line[:7].strip(): # while it's the same SFO
                        # i.e. while not completely blank, but blank at the start
                        info = line[43:].split()
                        if len(info)>0: # len(info)==0 for the second line of dvb_ir.adfout
                            frag += "+" + fragname + info[-1]
                            coeff = float(info[-4])
                            if coeff < 0:
                                orbital += '-' + info[-3] + info[-2].replace(":", "")
                            else:
                                orbital += '+' + info[-3] + info[-2].replace(":", "")
                        line = inputfile.next()
                    # At this point, we are either at the start of the next SFO or at
                    # a blank line...the end

                    self.fonames.append("%s_%s" % (frag, orbital))
                symoffset += num

                # blankline blankline
                inputfile.next(); inputfile.next()

        if line[1:32] == "S F O   P O P U L A T I O N S ,":
        #Extract overlap matrix

            self.fooverlaps = numpy.zeros((self.nbasis, self.nbasis), "d")

            symoffset = 0

            for nosymrep in self.nosymreps:

                line = inputfile.next()
                while line.find('===') < 10: #look for the symmetry labels
                    line = inputfile.next()
                #blank blank text blank col row
                for i in range(6):
                    inputfile.next()

                base = 0
                while base < nosymrep: #have we read all the columns?

                    for i in range(nosymrep - base):

                        self.updateprogress(inputfile, "Overlap", self.fupdate)
                        line = inputfile.next()
                        parts = line.split()[1:]
                        for j in range(len(parts)):
                            k = float(parts[j])
                            self.fooverlaps[base + symoffset + j, base + symoffset +i] = k
                            self.fooverlaps[base + symoffset + i, base + symoffset + j] = k

                    #blank, blank, column
                    for i in range(3):
                        inputfile.next()

                    base += 4

                symoffset += nosymrep
                base = 0

# The commented code below makes the atombasis attribute based on the BAS function in ADF,
#   but this is probably not so useful, since SFOs are used to build MOs in ADF.
#        if line[1:54] == "BAS: List of all Elementary Cartesian Basis Functions":
#
#            self.atombasis = []
#
#            # There will be some text, followed by a line:
#            #       (power of) X  Y  Z  R     Alpha  on Atom
#            while not line[1:11] == "(power of)":
#                line = inputfile.next()
#            dashes = inputfile.next()
#            blank = inputfile.next()
#            line = inputfile.next()
#            # There will be two blank lines when there are no more atom types.
#            while line.strip() != "":
#                atoms = [int(i)-1 for i in line.split()[1:]]
#                for n in range(len(atoms)):
#                    self.atombasis.append([])
#                dashes = inputfile.next()
#                line = inputfile.next()
#                while line.strip() != "":
#                    indices = [int(i)-1 for i in line.split()[5:]]
#                    for i in range(len(indices)):
#                        self.atombasis[atoms[i]].append(indices[i])
#                    line = inputfile.next()
#                line = inputfile.next()

        if line[48:67] == "SFO MO coefficients":

            self.mocoeffs = [numpy.zeros((self.nbasis, self.nbasis), "d")]
            spin = 0
            symoffset = 0
            lastrow = 0

            # Section ends with "1" at beggining of a line.
            while line[0] != "1":
                line = inputfile.next()

                # If spin is specified, then there will be two coefficient matrices. 
                if line.strip() == "***** SPIN 1 *****":
                    self.mocoeffs = [numpy.zeros((self.nbasis, self.nbasis), "d"),
                                     numpy.zeros((self.nbasis, self.nbasis), "d")]

                # Bump up the spin.
                if line.strip() == "***** SPIN 2 *****":
                    spin = 1
                    symoffset = 0
                    lastrow = 0

                # Next symmetry.
                if line.strip()[:4] == "=== ":
                    sym = line.split()[1]
                    if self.nosymflag:
                        aolist = range(self.nbasis)
                    else:
                        aolist = self.symlist[sym][spin]
                    # Add to the symmetry offset of AO ordering.
                    symoffset += lastrow

                # Blocks with coefficient always start with "MOs :".
                if line[1:6] == "MOs :":
                    # Next line has the MO index contributed to.
                    monumbers = [int(n) for n in line[6:].split()]
                    occup = inputfile.next()
                    label = inputfile.next()
                    line = inputfile.next()
                    # The table can end with a blank line or "1".
                    row = 0
                    while not line.strip() in ["", "1"]:
                        info = line.split()

                        if int(info[0]) < self.start_indeces[sym]:
                        #check to make sure we aren't parsing CFs
                            line = inputfile.next()
                            continue

                        self.updateprogress(inputfile, "Coefficients", self.fupdate)
                        row += 1
                        coeffs = [float(x) for x in info[1:]]
                        moindices = [aolist[n-1] for n in monumbers]
                        # The AO index is 1 less than the row.
                        aoindex = symoffset + row - 1
                        for i in range(len(monumbers)):
                            self.mocoeffs[spin][moindices[i],aoindex] = coeffs[i]
                        line = inputfile.next()
                    lastrow = row

        if line[4:53] == "Final excitation energies from Davidson algorithm":

            # move forward in file past some various algorthm info

            # *   Final excitation energies from Davidson algorithm                    *
            # *                                                                        *
            # **************************************************************************

            #     Number of loops in Davidson routine     =   20                    
            #     Number of matrix-vector multiplications =   24                    
            #     Type of excitations = SINGLET-SINGLET 

            inputfile.next(); inputfile.next(); inputfile.next()
            inputfile.next(); inputfile.next(); inputfile.next()
            inputfile.next(); inputfile.next()

            symm = self.normalisesym(inputfile.next().split()[1])

            # move forward in file past some more txt and header info

            # Excitation energies E in a.u. and eV, dE wrt prev. cycle,
            # oscillator strengths f in a.u.

            # no.  E/a.u.        E/eV      f           dE/a.u.
            # -----------------------------------------------------

            inputfile.next(); inputfile.next(); inputfile.next()
            inputfile.next(); inputfile.next(); inputfile.next()

            # now start parsing etenergies and etoscs

            etenergies = []
            etoscs = []
            etsyms = []

            line = inputfile.next()
            while len(line) > 2:
                info = line.split()
                etenergies.append(utils.convertor(float(info[2]), "eV", "cm-1"))
                etoscs.append(float(info[3]))
                etsyms.append(symm)
                line = inputfile.next()

            # move past next section
            while line[1:53] != "Major MO -> MO transitions for the above excitations":
                line = inputfile.next()

            # move past headers

            #  Excitation  Occupied to virtual  Contribution                         
            #   Nr.          orbitals           weight        contribibutions to      
            #                                   (sum=1) transition dipole moment   
            #                                             x       y       z       

            inputfile.next(), inputfile.next(), inputfile.next()
            inputfile.next(), inputfile.next(), inputfile.next()

            # before we start handeling transitions, we need
            # to create mosyms with indices
            # only restricted calcs are possible in ADF

            counts = {}
            syms = []
            for mosym in self.mosyms[0]:
                if counts.keys().count(mosym) == 0:
                    counts[mosym] = 1
                else:
                    counts[mosym] += 1

                syms.append(str(counts[mosym]) + mosym)

            import re
            etsecs = []
            printed_warning = False 

            for i in range(len(etenergies)):
                etsec = []
                line = inputfile.next()
                info = line.split()
                while len(info) > 0:

                    match = re.search('[^0-9]', info[1])
                    index1 = int(info[1][:match.start(0)])
                    text = info[1][match.start(0):]
                    symtext = text[0].upper() + text[1:]
                    sym1 = str(index1) + self.normalisesym(symtext)

                    match = re.search('[^0-9]', info[3])
                    index2 = int(info[3][:match.start(0)])
                    text = info[3][match.start(0):]
                    symtext = text[0].upper() + text[1:]
                    sym2 = str(index2) + self.normalisesym(symtext)

                    try:
                        index1 = syms.index(sym1)
                    except ValueError:
                        if not printed_warning:
                            self.logger.warning("Etsecs are not accurate!")
                            printed_warning = True

                    try:
                        index2 = syms.index(sym2)
                    except ValueError:
                        if not printed_warning:
                            self.logger.warning("Etsecs are not accurate!")
                            printed_warning = True

                    etsec.append([(index1, 0), (index2, 0), float(info[4])])

                    line = inputfile.next()
                    info = line.split()

                etsecs.append(etsec)


            if not hasattr(self, "etenergies"):
                self.etenergies = etenergies
            else:
                self.etenergies += etenergies

            if not hasattr(self, "etoscs"):
                self.etoscs = etoscs
            else:
                self.etoscs += etoscs

            if not hasattr(self, "etsyms"):
                self.etsyms = etsyms
            else:
                self.etsyms += etsyms

            if not hasattr(self, "etsecs"):
                self.etsecs = etsecs
            else:
                self.etsecs += etsecs

if __name__ == "__main__":
    import doctest, adfparser
    doctest.testmod(adfparser, verbose=False)
