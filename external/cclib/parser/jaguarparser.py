"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 861 $"


import re

import numpy

import logfileparser
import utils


class Jaguar(logfileparser.Logfile):
    """A Jaguar output file"""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Jaguar, self).__init__(logname="Jaguar", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Jaguar output file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Jaguar("%s")' % (self.filename)

    def normalisesym(self, label):
        """Normalise the symmetries used by Jaguar.

        To normalise, three rules need to be applied:
        (1) To handle orbitals of E symmetry, retain everything before the /
        (2) Replace two p's by "
        (2) Replace any remaining single p's by '

        >>> t = Jaguar("dummyfile").normalisesym
        >>> labels = ['A', 'A1', 'Ag', 'Ap', 'App', "A1p", "A1pp", "E1pp/Ap"]
        >>> answers = map(t, labels)
        >>> print answers
        ['A', 'A1', 'Ag', "A'", 'A"', "A1'", 'A1"', 'E1"']
        """
        ans = label.split("/")[0].replace("pp", '"').replace("p", "'")
        return ans

    def before_parsing(self):

        self.geoopt = False # Is this a GeoOpt? Needed for SCF targets/values.

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""
            
        if line[0:4] == "etot":
        # Get SCF convergence information
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []
                self.scftargets = [[5E-5, 5E-6]]
            values = []
            while line[0:4] == "etot":
        # Jaguar 4.2
        # etot   1  N  N  0  N  -382.08751886450           2.3E-03  1.4E-01
        # etot   2  Y  Y  0  N  -382.27486023153  1.9E-01  1.4E-03  5.7E-02
        # Jaguar 6.5
        # etot   1  N  N  0  N    -382.08751881733           2.3E-03  1.4E-01
        # etot   2  Y  Y  0  N    -382.27486018708  1.9E-01  1.4E-03  5.7E-02
                temp = line.split()[7:]
                if len(temp)==3:
                    denergy = float(temp[0])
                else:
                    denergy = 0 # Should really be greater than target value
                                # or should we just ignore the values in this line
                ddensity = float(temp[-2])
                maxdiiserr = float(temp[-1])
                if not self.geoopt:
                    values.append([denergy, ddensity])
                else:
                    values.append([ddensity])
                line = inputfile.next()
            self.scfvalues.append(values)

        # Hartree-Fock energy after SCF
        if line[1:18] == "SCFE: SCF energy:":
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            temp = line.strip().split()
            scfenergy = float(temp[temp.index("hartrees") - 1])
            scfenergy = utils.convertor(scfenergy, "hartree", "eV")
            self.scfenergies.append(scfenergy)

        # Energy after LMP2 correction
        if line[1:18] == "Total LMP2 Energy":
            if not hasattr(self, "mpenergies"):
                self.mpenergies = [[]]
            lmp2energy = float(line.split()[-1])
            lmp2energy = utils.convertor(lmp2energy, "hartree", "eV")
            self.mpenergies[-1].append(lmp2energy)

        if line[2:14] == "new geometry" or line[1:21] == "Symmetrized geometry" or line.find("Input geometry") > 0:
        # Get the atom coordinates
            if not hasattr(self, "atomcoords") or line[1:21] == "Symmetrized geometry":
                # Wipe the "Input geometry" if "Symmetrized geometry" present
                self.atomcoords = []
            p = re.compile("(\D+)\d+") # One/more letters followed by a number
            atomcoords = []
            atomnos = []
            angstrom = inputfile.next()
            title = inputfile.next()
            line = inputfile.next()
            while line.strip():
                temp = line.split()
                element = p.findall(temp[0])[0]
                atomnos.append(self.table.number[element])
                atomcoords.append(map(float, temp[1:]))
                line = inputfile.next()
            self.atomcoords.append(atomcoords)
            self.atomnos = numpy.array(atomnos, "i")
            self.natom = len(atomcoords)

        # Extract charge and multiplicity
        if line[2:22] == "net molecular charge":
            self.charge = int(line.split()[-1])
            self.mult = int(inputfile.next().split()[-1])

        if line[2:24] == "start of program geopt":
            if not self.geoopt:
                # Need to keep only the RMS density change info
                # if this is a geoopt
                self.scftargets = [[self.scftargets[0][0]]]
                if hasattr(self, "scfvalues"):
                    self.scfvalues[0] = [[x[0]] for x in self.scfvalues[0]]
                self.geoopt = True
            else:
                self.scftargets.append([5E-5])

        if line[2:28] == "geometry optimization step":
        # Get Geometry Opt convergence information
            if not hasattr(self, "geovalues"):
                self.geovalues = []
                self.geotargets = numpy.zeros(5, "d")
            gopt_step = int(line.split()[-1])
            energy = inputfile.next()
            # quick hack for messages of the sort:
            #   ** restarting optimization from step    2 **
            # as found in regression file ptnh3_2_H2O_2_2plus.out
            if inputfile.next().strip():
                blank = inputfile.next()
            line = inputfile.next()
            values = []
            target_index = 0                
            if gopt_step == 1:
                # The first optimization step does not produce an energy change
                values.append(0.0)
                target_index = 1
            while line.strip():
                if len(line) > 40 and line[41] == "(":
                    # A new geo convergence value
                    values.append(float(line[26:37]))
                    self.geotargets[target_index] = float(line[43:54])
                    target_index += 1
                line = inputfile.next()
            self.geovalues.append(values)

        if line.find("number of occupied orbitals") > 0:
        # Get number of MOs
            occs = int(line.split()[-1])
            line = inputfile.next()
            virts = int(line.split()[-1])
            self.nmo = occs + virts
            self.homos = numpy.array([occs-1], "i")

            self.unrestrictedflag = False

        if line.find("number of alpha occupied orb") > 0:
        # Get number of MOs for an unrestricted calc

            aoccs = int(line.split()[-1])
            line = inputfile.next()
            avirts = int(line.split()[-1])
            line = inputfile.next()
            boccs = int(line.split()[-1])
            line = inputfile.next()
            bvirt = int(line.split()[-1])

            self.nmo = aoccs + avirts
            self.homos = numpy.array([aoccs-1,boccs-1], "i")
            self.unrestrictedflag = True

        # MO energies and symmetries.
        # Jaguar 7.0: provides energies and symmetries for both
        #   restricted and unrestricted calculations, like this:
        #     Alpha Orbital energies/symmetry label: 
        #     -10.25358 Bu  -10.25353 Ag  -10.21931 Bu  -10.21927 Ag     
        #     -10.21792 Bu  -10.21782 Ag  -10.21773 Bu  -10.21772 Ag     
        #     ...
        # Jaguar 6.5: prints both only for restricted calculations,
        #   so for unrestricted calculations the output it looks like this:
        #     Alpha Orbital energies: 
        #     -10.25358  -10.25353  -10.21931  -10.21927  -10.21792  -10.21782
        #     -10.21773  -10.21772  -10.21537  -10.21537   -1.02078   -0.96193
        #     ...
        # Presence of 'Orbital energies' is enough to catch all versions.
        if "Orbital energies" in line:

            # Parsing results is identical for restricted/unrestricted
            #   calculations, just assert later that alpha/beta order is OK.
            spin = int(line[2:6] == "Beta")

            # Check if symmetries are printed also.
            issyms = "symmetry label" in line

            if not hasattr(self, "moenergies"):
                self.moenergies = []
            if issyms and not hasattr(self, "mosyms"):
                    self.mosyms = []
            
            # Grow moeneriges/mosyms and make sure they are empty when
            #   parsed multiple times - currently cclib returns only
            #   the final output (ex. in a geomtry optimization).
            if len(self.moenergies) < spin+1:
                self.moenergies.append([])
            self.moenergies[spin] = []
            if issyms:
                if len(self.mosyms) < spin+1:
                    self.mosyms.append([])
                self.mosyms[spin] = []
            
            line = inputfile.next().split()
            while len(line) > 0:
                if issyms:
                    energies = [float(line[2*i]) for i in range(len(line)/2)]
                    syms = [line[2*i+1] for i in range(len(line)/2)]
                else:
                    energies = [float(e) for e in line]
                energies = [utils.convertor(e, "hartree", "eV") for e in energies]
                self.moenergies[spin].extend(energies)
                if issyms:
                    syms = [self.normalisesym(s) for s in syms]
                    self.mosyms[spin].extend(syms)
                line = inputfile.next().split()
            
            # There should always be an extra blank line after all this.
            line = inputfile.next()

        if line.find("Occupied + virtual Orbitals- final wvfn") > 0:
            
            blank = inputfile.next()
            stars = inputfile.next()
            blank = inputfile.next()
            blank = inputfile.next()
            
            if not hasattr(self,"mocoeffs"):
                if self.unrestrictedflag:
                    spin = 2
                else:
                    spin = 1

                self.mocoeffs = []
                
            
            aonames = []
            lastatom = "X"
            
            readatombasis = False
            if not hasattr(self, "atombasis"):
                self.atombasis = []
                for i in range(self.natom):
                    self.atombasis.append([])
                readatombasis = True

            offset = 0

            for s in range(spin):
                mocoeffs = numpy.zeros((len(self.moenergies[s]), self.nbasis), "d")

                if s == 1: #beta case
                    stars = inputfile.next()
                    blank = inputfile.next()
                    title = inputfile.next()
                    blank = inputfile.next()
                    stars = inputfile.next()
                    blank = inputfile.next()
                    blank = inputfile.next()

                for k in range(0,len(self.moenergies[s]),5):

                    numbers = inputfile.next()
                    eigens = inputfile.next()
                    line = inputfile.next()

                    for i in range(self.nbasis):

                        info = line.split()
                        
                        # Fill atombasis only first time around.
                        if readatombasis and k == 0:
                            orbno = int(info[0])
                            atom = info[1]
                            if atom[1].isalpha():
                                atomno = int(atom[2:])
                            else:
                                atomno = int(atom[1:])
                            self.atombasis[atomno-1].append(orbno-1)

                        if not hasattr(self,"aonames"):
                            if lastatom != info[1]:
                                scount = 1
                                pcount = 3
                                dcount = 6 #six d orbitals in Jaguar

                            if info[2] == 'S':
                                aonames.append("%s_%i%s"%(info[1], scount, info[2]))
                                scount += 1
                        
                            if info[2] == 'X' or info[2] == 'Y' or info[2] == 'Z':
                                aonames.append("%s_%iP%s"%(info[1], pcount / 3, info[2]))
                                pcount += 1
                        
                            if info[2] == 'XX' or info[2] == 'YY' or info[2] == 'ZZ' or \
                               info[2] == 'XY' or info[2] == 'XZ' or info[2] == 'YZ':

                                aonames.append("%s_%iD%s"%(info[1], dcount / 6, info[2]))
                                dcount += 1

                            lastatom = info[1]

                        for j in range(len(info[3:])):
                            mocoeffs[j+k,i] = float(info[3+j])

                        line = inputfile.next()

                    if not hasattr(self,"aonames"):
                        self.aonames = aonames

                    offset += 5
                self.mocoeffs.append(mocoeffs)
                        
                        
        if line[2:6] == "olap":
            if line[6]=="-":
                return
                # This was continue (in loop) before parser refactoring.
                # continue # avoid "olap-dev"
            self.aooverlaps = numpy.zeros((self.nbasis, self.nbasis), "d")

            for i in range(0, self.nbasis, 5):
                blank = inputfile.next()
                header = inputfile.next()
                for j in range(i, self.nbasis):
                    temp = map(float, inputfile.next().split()[1:])
                    self.aooverlaps[j, i:(i+len(temp))] = temp
                    self.aooverlaps[i:(i+len(temp)), j] = temp
            
        if line[1:28] == "number of occupied orbitals":
            self.homos = numpy.array([float(line.strip().split()[-1])-1], "i")

        if line[2:27] == "number of basis functions":
            self.nbasis = int(line.strip().split()[-1])

        # IR output looks like this:
        #   frequencies        72.45   113.25   176.88   183.76   267.60   312.06
        #   symmetries       Au       Bg       Au       Bu       Ag       Bg      
        #   intensities         0.07     0.00     0.28     0.52     0.00     0.00
        #   reduc. mass         1.90     0.74     1.06     1.42     1.19     0.85
        #   force const         0.01     0.01     0.02     0.03     0.05     0.05
        #   C1       X     0.00000  0.00000  0.00000 -0.05707 -0.06716  0.00000
        #   C1       Y     0.00000  0.00000  0.00000  0.00909 -0.02529  0.00000
        #   C1       Z     0.04792 -0.06032 -0.01192  0.00000  0.00000  0.11613
        #   C2       X     0.00000  0.00000  0.00000 -0.06094 -0.04635  0.00000
        #   ... etc. ...
        # This is a complete ouput, some files will not have intensities,
        #   and older Jaguar versions sometimes skip the symmetries.
        if line[2:23] == "start of program freq":

            self.vibfreqs = []
            self.vibdisps = []
            forceconstants = False
            intensities = False
            blank = inputfile.next()
            line = inputfile.next()
            while line.strip():
                if "force const" in line:
                    forceconstants = True
                if "intensities" in line:
                    intensities = True
                line = inputfile.next()
            freqs = inputfile.next()
            
            # The last block has an extra blank line after it - catch it.
            while freqs.strip():

                # Number of modes (columns printed in this block).
                nmodes = len(freqs.split())-1

                # Append the frequencies.
                self.vibfreqs.extend(map(float, freqs.split()[1:]))
                line = inputfile.next().split()
                
                # May skip symmetries (older Jaguar versions).
                if line[0] == "symmetries":
                    if not hasattr(self, "vibsyms"):
                        self.vibsyms = []
                    self.vibsyms.extend(map(self.normalisesym, line[1:]))
                    line = inputfile.next().split()                                
                if intensities:
                    if not hasattr(self, "vibirs"):
                        self.vibirs = []
                    self.vibirs.extend(map(float, line[1:]))
                    line = inputfile.next().split()                                
                if forceconstants:
                    line = inputfile.next()

                # Start parsing the displacements.
                # Variable 'q' holds up to 7 lists of triplets.
                q = [ [] for i in range(7) ]
                for n in range(self.natom):
                    # Variable 'p' holds up to 7 triplets.
                    p = [ [] for i in range(7) ]
                    for i in range(3):
                        line = inputfile.next()
                        disps = [float(disp) for disp in line.split()[2:]]
                        for j in range(nmodes):
                            p[j].append(disps[j])
                    for i in range(nmodes):
                        q[i].append(p[i])

                self.vibdisps.extend(q[:nmodes])
                blank = inputfile.next()
                freqs = inputfile.next()

            # Convert new data to arrays.
            self.vibfreqs = numpy.array(self.vibfreqs, "d")
            self.vibdisps = numpy.array(self.vibdisps, "d")
            if hasattr(self, "vibirs"):
                self.vibirs = numpy.array(self.vibirs, "d")
                
        # Parse excited state output (for CIS calculations).
        # Jaguar calculates only singlet states.
        if line[2:15] == "Excited State":
            if not hasattr(self, "etenergies"):
                self.etenergies = []
            if not hasattr(self, "etoscs"):
                self.etoscs = []
            if not hasattr(self, "etsecs"):
                self.etsecs = []
                self.etsyms = []
            etenergy = float(line.split()[3])
            etenergy = utils.convertor(etenergy, "eV", "cm-1")
            self.etenergies.append(etenergy)
            # Skip 4 lines
            for i in range(5):
                line = inputfile.next()
            self.etsecs.append([])
            # Jaguar calculates only singlet states.
            self.etsyms.append('Singlet-A')
            while line.strip() != "":
                fromMO = int(line.split()[0])-1
                toMO = int(line.split()[2])-1
                coeff = float(line.split()[-1])
                self.etsecs[-1].append([(fromMO,0),(toMO,0),coeff])
                line = inputfile.next()
            # Skip 3 lines
            for i in range(4):
                line = inputfile.next()
            strength = float(line.split()[-1])
            self.etoscs.append(strength)

        
if __name__ == "__main__":
    import doctest, jaguarparser
    doctest.testmod(jaguarparser, verbose=False)
